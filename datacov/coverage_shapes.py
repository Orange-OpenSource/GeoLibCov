from shapely import MultiPoint, Polygon
from shapely import voronoi_polygons
from shapely.ops import linemerge, unary_union, polygonize
import geopandas as gpd
import pandas as pd
import numpy as np
from datacov import TopoGen
from shapely import concave_hull
from shapely.affinity import scale

class Geom:
    def __init__(self, topo=None) -> None:
        self.topo = None
        self.sites = None
        self.cells = None
        
        if topo is not None:
            self._init_topo_attr(topo)
        
        self.site_shapes = None
        self.cell_shapes = None
        self.radius_scale = None 

    def _init_topo_attr(self, topo):
        self.topo = topo
        self.sites = self.topo.sites.copy()
        self.cells = self.topo.cells.copy()      
        
    def plot(self, band=None, **kwargs):
        if band is None:
            self.cell_shapes.plot(**kwargs)
        else:
            self.cell_shapes[self.cell_shapes.band == band].plot(**kwargs)
            
    def set_shape_radius_scale(self):
        """sets an approximation of the radius of the bounding circle of the polygons. Approximation is the largest side of the bounding box / 2
        
        Returns:
            _type_: _description_
        """
        b = self.cell_shapes.geometry.bounds
        bx_len = b.maxx - b.minx
        by_len = b.maxy - b.miny
        approx_radius = np.max(np.array([bx_len, by_len]).T, axis=1) / 2
        radius_scale = pd.DataFrame({'cell_id': self.cell_shapes.cell_id, 'model_radius': approx_radius})
        self.radius_scale = radius_scale
        return radius_scale
        
        
    def create_shapes(self, topo=None):
        """Initialization method.

        Args:
            topo (TopoGen, optional): generated topography object. Defaults to None.

        Raises:
            RuntimeError: if topo has not been specified at constructor call, it must be specified when calling this method.

        Returns:
            self: returns itself
        """
        if (self.topo is None) and (topo is None):
            raise RuntimeError(f'topo must be specified at instanciation or in this call.')
        
        if self.topo is None:
            self._init_topo_attr(topo)
            
        self.site_shapes = self.create_site_shapes()
        self.cell_shapes = self.create_cell_shapes()
        
        self.set_shape_radius_scale()
        return self
        
    def create_site_shapes(self):
        raise NotImplementedError('Method has not been overriden.')
        
    def scale_cell_shapes(self, coverage, sfact=1.1):
        """Initialization method

        Args:
            coverage (Geom): An initialized Geom object on which the scaling of cell_shapes is based.
            sfact (float, optional): scaling factor. Defaults to 1.1.

        Raises:
            RuntimeError: if topo has not been specified at constructor call, it must be specified when calling this method.

        Returns:
            self: returns itself
        """
        if (self.topo is None) and (coverage.topo is None):
            raise RuntimeError(f'Cannot set topo because topo of coverage is None.')
        
        if self.topo is None:
            self._init_topo_attr(coverage.topo)
            
        if isinstance(sfact, list):
            if len(sfact) != coverage.cell_shapes.shape[0]:
                raise AssertionError(f'List of factors with length {len(sfact)} does not match length {coverage.cell_shapes.shape[0]} of coverage.cell_shapes.')
        else:
            sfact = [sfact] * coverage.cell_shapes.shape[0]
        
        s = []
        cell_shapes = coverage.cell_shapes.copy()
        for i, r in enumerate(cell_shapes.iterrows()):
            values = r[1]
            s.append(scale(values.geometry, xfact=sfact[i], yfact=sfact[i], origin=(values.x, values.y)))
        cell_shapes['geometry'] = s
        self.cell_shapes = cell_shapes
        
        self.set_shape_radius_scale()
        return self
        
    def create_cell_shapes(self):       
        az_start = self.topo.azimuts.reshape((self.topo.n_sites, self.topo.n_azimuts))
        az_end = np.hstack([az_start[:, 1:], (az_start[:, 0] + 360).reshape((-1, 1))])
        sector_azimuts = ((az_start + az_end) / 2)%360

        cut_length = (self.topo.maxlim - self.topo.minlim) * 2
        sector_geometries = self.topo.create_az_lines(sector_azimuts, {700: cut_length, 800: cut_length, 1800: cut_length, 2100: cut_length, 2600: cut_length})
        sector_geometries = sector_geometries.merge(self.site_shapes[['bs_id', 'band']], on=['bs_id', 'band'])[['bs_id', 'band', 'geometry']]
        sector_geometries = sector_geometries.dissolve(by=['bs_id', 'band']).reset_index()
        sector_geometries = sector_geometries.sort_values(by=['bs_id', 'band'])
        
        polygons = []
        bands = []
        for r1, r2 in zip(sector_geometries.iterrows(), self.site_shapes.iterrows()):
            band = r1[1].band
            g1 = list(r1[1].geometry.geoms)
            g2 = r2[1].geometry
            p = list(self.split_polygon(g1, g2))
            polygons += p
            bands = bands + [band] * self.topo.n_azimuts
        vorocells = gpd.GeoDataFrame({'band': bands}, geometry=polygons)
        
        cell_shapes = []
        sect_len = 0.001
        section = self.topo.create_az_lines(ue_dist_loc={700: sect_len, 800: sect_len, 1800: sect_len, 2100: sect_len, 2600: sect_len}, use_cid=True)
        for b in self.topo.bands:
            # TODO: may not need a loop here
            vor = vorocells[vorocells.band==b].drop('band', axis=1)
            sec = section[section.band == b]
            cell_shapes.append(vor.sjoin(sec, predicate='contains'))
            
        cell_shapes = pd.concat(cell_shapes)
        cell_shapes = cell_shapes.merge(self.topo.sites[['bs_id', 'x', 'y']], on='bs_id')
        assert cell_shapes.shape[0] == self.cells.shape[0], f"cell_shapes.shape[0]={cell_shapes.shape[0]} and self.cells.shape[0]={self.cells.shape[0]}"
        return cell_shapes[['cell_id', 'bs_id', 'x', 'y', 'band', 'az_id', 'azimut', 'geometry']].sort_values(by=['bs_id', 'band', 'az_id'])
    
    def split_polygon(self, lines, polygon):
        """split polygon with lines

        Args:
            lines (list): list of linestrings
            polygon (Polygon): shapely polygon
        """
        merged = linemerge(lines + [polygon.boundary])
        borders = unary_union(merged)
        polygons = polygonize(borders)
        return polygons
    

class VoroGeom(Geom):
    def __init__(self, topo=None) -> None:
        super().__init__(topo=topo)
    
    def create_site_shapes(self):
        geometries = []
        for b in self.topo.bands:
            site_band = self.sites.merge(self.cells[self.cells.band == b][['bs_id', 'band']], on='bs_id')
            site_band = site_band[['bs_id', 'x', 'y', 'height', 'band', 'geometry']].drop_duplicates()
            
            if site_band.shape[0] == 1:
                gs_shape = [Polygon(((0., 0.), (0., 1.), (1., 1.), (1., 0.), (0., 0.)))]
            else:
                g = voronoi_polygons(MultiPoint(site_band.geometry.to_list()))
                gs_shape = gpd.GeoSeries(g).explode(index_parts=False)
            gdf_shape = gpd.GeoDataFrame(geometry=gs_shape)

            gdf_shape = gpd.sjoin(gdf_shape, site_band, predicate='contains')
            gdf_shape = gdf_shape.drop(['index_right'], axis=1)[['bs_id', 'x', 'y', 'height', 'band', 'geometry']]
            geometries.append(gdf_shape)
        geometries = pd.concat(geometries).sort_values(by=['bs_id', 'band']).reset_index(drop=True)
        assert self.cells[['bs_id', 'band']].drop_duplicates().shape[0] == geometries.shape[0]
        return geometries
    

class VoroSiteGeom(Geom):
    def __init__(self, topo=None) -> None:
        super().__init__(topo)
        
        self.site_shapes = None
        self.cell_shapes = None
        
    def create_shapes(self, coverage):
        if (self.topo is None) and (coverage.topo is None):
            raise RuntimeError(f'topo must be specified at instanciation or in this call.')
        
        if self.topo is None:
            self._init_topo_attr(coverage.topo)
            
        self.site_shapes = self.create_site_shapes(coverage)
        self.cell_shapes = self.create_cell_shapes(coverage)
        
        self.set_shape_radius_scale()
        return self
    
    def create_cell_shapes(self, coverage):
        cell_shapes = coverage.site_shapes.merge(coverage.cell_shapes.drop('geometry', axis=1), on=['bs_id', 'x', 'y', 'band'])
        assert cell_shapes.shape[0] == self.cells.shape[0]
        return cell_shapes
    
    def create_site_shapes(self, coverage):
        return coverage.site_shapes
    
    
class CircleGeom(Geom):
    def __init__(self, topo=None) -> None:
        super().__init__(topo=topo)
        self.shape_radius = 0.01
    
    def create_site_shapes(self):
        geometries = []
        for b in self.topo.bands:
            site_band = self.sites.merge(self.cells[self.cells.band == b][['bs_id', 'band']], on='bs_id')
            site_band = site_band[['bs_id', 'x', 'y', 'height', 'band', 'geometry']].drop_duplicates()
            
            gs_shape = site_band.geometry.buffer(self.shape_radius)
            gdf_shape = gpd.GeoDataFrame(geometry=gs_shape)

            gdf_shape = gpd.sjoin(gdf_shape, site_band, predicate='contains')
            gdf_shape = gdf_shape.drop(['index_right'], axis=1)[['bs_id', 'x', 'y', 'height', 'band', 'geometry']]
            geometries.append(gdf_shape)
        geometries = pd.concat(geometries).sort_values(by=['bs_id', 'band']).reset_index(drop=True)
        assert self.cells[['bs_id', 'band']].drop_duplicates().shape[0] == geometries.shape[0]
        return geometries
    

class EmpiricalGeom(Geom):
    def __init__(self, topo=None) -> None:
        self.ues = None
        super().__init__(topo=topo)
        
    def _init_topo_attr(self, topo):
        super()._init_topo_attr(topo)
        self.ues = topo.ues.copy()
        
    def create_shapes(self, topo=None, hull='convex'):
        if (self.topo is None) and (topo is None):
            raise RuntimeError(f'topo must be specified at instanciation or in this call.')
        
        if self.topo is None:
            self._init_topo_attr(topo)
            
        cell_shapes = self.ues.merge(self.cells.drop(['geometry'], axis=1), on='cell_id')
        if hull == 'convex':
            cell_shapes['geometry'] = self.ues.geometry.convex_hull
        if hull == 'concave':
            cell_shapes['geometry'] = concave_hull(self.ues.geometry, ratio=0)
        self.cell_shapes = cell_shapes
        return self
    

class DistGen:
    def __init__(self, topo) -> None:
        self.ue_dist_loc = topo.ue_dist_loc
        self.ue_dist_scale = topo.ue_dist_scale
        self.cells = topo.cells.copy()
        
        self.distances = None
    
    def generate_distances(self):
        c = self.cells.groupby('band').count()['cell_id']
        dists = []
        for band, size in zip(c.index, c):
            loc = self.ue_dist_loc[band] + 2 * self.ue_dist_scale[band]
            scale = self.ue_dist_scale[band] * 0.1
            dists.append(np.abs(np.random.normal(loc, scale, size)))
        dists = np.hstack(dists)
        self.distances = pd.DataFrame({'cell_id': self.cells.cell_id, 'band': self.cells.band, 'proxy_radius': dists})
        return self
    
    
if __name__ == '__main__':
    topo_test = TopoGen()
    gdf = topo_test.generate_topo()
    
    dg = DistGen(topo_test).generate_distances()
    
    v = VoroGeom(topo_test)
    v.create_shapes()
    print(v.radius_scale)
    c = CircleGeom(topo_test)
    c.create_shapes()
    
    eg = EmpiricalGeom().create_shapes(topo_test)
    eg2 = EmpiricalGeom(topo_test)
    eg2.create_shapes()

    v2 = VoroGeom(topo_test).scale_cell_shapes(v)
    v3 = VoroGeom().create_shapes(topo_test)
    
    v4 = VoroGeom().scale_cell_shapes(v)
    v4.plot()
    
    v5 = VoroSiteGeom().create_shapes(v)