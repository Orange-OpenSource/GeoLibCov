# Software Name: DataCov
# Version: 0.1.0
# SPDX-FileCopyrightText: Copyright (c) 2023 Orange
# SPDX-License-Identifier: BSD-3-Clause
#
# This software is distributed under the BSD 3-Clause "New" or "Revised" License,
# see the "LICENSE.txt" file for more details.
#
# Author: Danny Qiu <danny.qiu@orange.com>

from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import singledispatchmethod
from typing import overload
from typing_extensions import override
from shapely import MultiPoint, Polygon
from shapely import voronoi_polygons, minimum_bounding_radius
from shapely.ops import linemerge, unary_union, polygonize
import geopandas as gpd
import pandas as pd
import numpy as np
from datacov import TopoGen
from shapely import concave_hull
from shapely.affinity import scale

@dataclass
class Geom(ABC):
    topo: TopoGen
    site_shapes: pd.DataFrame
    cell_shapes: pd.DataFrame
    radius_scale: pd.DataFrame

    @classmethod
    def from_topo(cls, topo: TopoGen):
        site_shapes = cls._create_site_shapes(topo)
        cell_shapes = cls._create_cell_shapes(topo, site_shapes)
        radius_scale = cls._shape_radius_scale(site_shapes, cell_shapes)
        
        return cls(
            topo, 
            site_shapes,
            cell_shapes, 
            radius_scale,
        )

    def copy(self):
        return self.__class__(
            self.topo,
            self.site_shapes.copy(),
            self.cell_shapes.copy(),
            self.radius_scale.copy(),
        )
    
    @property
    def sites(self):
        return self.topo.sites
    
    @property
    def cells(self):
        return self.topo.cells
        
    @staticmethod
    def compute_scale(geom1, geom2):
        """get scale factor of geom2.model_radius / geom1.model_radius

        Args:
            geom1 (Geom): _description_
            geom2 (Geom): _description_

        Returns:
            Series: scale factor for each cell
        """
        compare = geom1.radius_scale.merge(geom2.radius_scale[['cell_id', 'model_radius']], on='cell_id', suffixes=['_geom1', '_geom2'])
        scale_factors = compare['model_radius_geom2'] / compare['model_radius_geom1']
        return scale_factors
        
    def plot(self, band=None, **kwargs):
        if band is None:
            self.cell_shapes.plot(**kwargs)
        else:
            self.cell_shapes[self.cell_shapes.band == band].plot(**kwargs)
            
    @staticmethod
    def _shape_radius_scale(site_shapes, cell_shapes):
        """sets an approximation of the radius of the bounding circle of site_shapes. Approximation is the largest side of the bounding box / 2
        
        Returns:
            _type_: _description_
        """
        radius = site_shapes.geometry.apply(minimum_bounding_radius)
        #bx_len = b.maxx - b.minx
        #by_len = b.maxy - b.miny
        #approx_radius = np.max(np.array([bx_len, by_len]).T, axis=1) / 2
        radius_scale = pd.DataFrame({'bs_id': site_shapes.bs_id, 'band': site_shapes.band, 'model_radius': radius})
        radius_scale = radius_scale.merge(cell_shapes, on=['bs_id', 'band'])
        assert cell_shapes.shape[0] == radius_scale.shape[0]
        return radius_scale
    
    @classmethod
    @abstractmethod
    def _create_site_shapes(cls, topo):
        pass
    
    @staticmethod
    def retrieve_scale_factors(r_l_model, radius_scale):
        df_scale = r_l_model.data.merge(radius_scale, on='cell_id')
        df_scale['scale_factor'] = df_scale['r_l'] / df_scale['model_radius']
        return list(df_scale['scale_factor'])
    
    @classmethod
    def scale_cell_shapes(cls, coverage, sfact=1.1, radius_model=None):
        """Initialization method

        Args:
            coverage (Geom): An initialized Geom object on which the scaling of cell_shapes is based.
            sfact (float, optional): scaling factor. Defaults to 1.1.

        Raises:
            RuntimeError: if topo has not been specified at constructor call, it must be specified when calling this method.

        Returns:
            self: returns itself
        """
        if radius_model is not None:
            sfact = cls.retrieve_scale_factors(radius_model, coverage.radius_scale)
        elif isinstance(sfact, list):
            if len(sfact) != coverage.cell_shapes.shape[0]:
                raise AssertionError(f'List of factors with length {len(sfact)} does not match length {coverage.cell_shapes.shape[0]} of coverage.cell_shapes.')
        else:
            sfact = [sfact] * coverage.cell_shapes.shape[0]

        self = coverage.copy()
        
        scaled_geom = []
        cell_shapes = coverage.cell_shapes.copy()
        for i, r in enumerate(cell_shapes.iterrows()):
            values = r[1]
            scaled_geom.append(scale(values.geometry, xfact=sfact[i], yfact=sfact[i], origin=(values.x, values.y)))
        cell_shapes['geometry'] = scaled_geom
        cell_shapes['scale'] = sfact
        self.cell_shapes = cell_shapes.copy().drop('scale', axis=1)
        
        sfact_site = cell_shapes.groupby(['bs_id', 'band']).scale.mean().reset_index()
        site_shapes = coverage.site_shapes.copy()
        site_shapes = site_shapes.merge(sfact_site, on=['bs_id', 'band'])
        
        scaled_site_geom = []
        for r in site_shapes.iterrows():
            values = r[1]
            scaled_site_geom.append(scale(values.geometry, xfact=values.scale, yfact=values.scale, origin=(values.x, values.y)))
        
        assert site_shapes.shape[0] == coverage.site_shapes.shape[0]
        self.site_shapes = site_shapes.drop('scale', axis=1)
        self.site_shapes.loc[:, 'geometry'] = gpd.GeoSeries(scaled_site_geom)
        self.radius_scale = cls._shape_radius_scale(self.site_shapes, self.cell_shapes)
        return self
        
    @staticmethod
    def _create_cell_shapes(topo, site_shapes):       
        az_start = topo.azimuts.reshape((topo.n_sites, topo.n_azimuts))
        az_end = np.hstack([az_start[:, 1:], (az_start[:, 0] + 360).reshape((-1, 1))])
        sector_azimuts = ((az_start + az_end) / 2)%360

        cut_length = (topo.maxlim - topo.minlim) * 2
        sector_geometries = topo.create_az_lines(sector_azimuts, {700: cut_length, 800: cut_length, 1800: cut_length, 2100: cut_length, 2600: cut_length})
        sector_geometries = sector_geometries.merge(site_shapes[['bs_id', 'band']], on=['bs_id', 'band'])[['bs_id', 'band', 'geometry']]
        sector_geometries = sector_geometries.dissolve(by=['bs_id', 'band']).reset_index()
        sector_geometries = sector_geometries.sort_values(by=['bs_id', 'band'])
        
        polygons = []
        bands = []
        for r1, r2 in zip(sector_geometries.iterrows(), site_shapes.iterrows()):
            band = r1[1].band
            g1 = list(r1[1].geometry.geoms)
            g2 = r2[1].geometry
            p = list(Geom.split_polygon(g1, g2))
            polygons += p
            bands = bands + [band] * topo.n_azimuts
        vorocells = gpd.GeoDataFrame({'band': bands}, geometry=polygons)
        
        cell_shapes = []
        sect_len = 0.001
        section = topo.create_az_lines(ue_dist_loc={700: sect_len, 800: sect_len, 1800: sect_len, 2100: sect_len, 2600: sect_len}, use_cid=True)
        for b in topo.bands:
            # TODO: may not need a loop here
            vor = vorocells[vorocells.band==b].drop('band', axis=1)
            sec = section[section.band == b]
            cell_shapes.append(vor.sjoin(sec, predicate='contains'))
            
        cell_shapes = pd.concat(cell_shapes)
        cell_shapes = cell_shapes.merge(topo.sites[['bs_id', 'x', 'y']], on='bs_id')
        assert cell_shapes.shape[0] == topo.cells.shape[0], f"cell_shapes.shape[0]={cell_shapes.shape[0]} and self.cells.shape[0]={topo.cells.shape[0]}"
        return cell_shapes[['cell_id', 'bs_id', 'x', 'y', 'band', 'az_id', 'azimut', 'geometry']].sort_values(by=['bs_id', 'band', 'az_id'])
    
    @staticmethod
    def split_polygon(lines, polygon):
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
    @override
    @classmethod
    def _create_site_shapes(cls, topo):
        geometries = []
        for b in topo.bands:
            site_band = topo.sites.merge(topo.cells[topo.cells.band == b][['bs_id', 'band']], on='bs_id')
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
        assert topo.cells[['bs_id', 'band']].drop_duplicates().shape[0] == geometries.shape[0]
        return geometries
    

class VoroSiteGeom(Geom):
    @classmethod
    @override
    def from_topo(cls, topo):
        voro_geom = VoroGeom.from_topo(topo)
        return VoroSiteGeom.from_voro_geom(voro_geom)

    @classmethod
    def from_voro_geom(cls, voro_geom: VoroGeom):
        self = voro_geom.copy()
        self.cell_shapes = self.site_shapes.merge(self.cell_shapes.drop('geometry', axis=1), on=['bs_id', 'x', 'y', 'band'])
        return self


class CircleGeom(Geom):
    SHAPE_RADIUS = 0.01
    
    @override
    @classmethod
    def _create_site_shapes(cls, topo):
        geometries = []
        for b in topo.bands:
            site_band = topo.sites.merge(topo.cells[topo.cells.band == b][['bs_id', 'band']], on='bs_id')
            site_band = site_band[['bs_id', 'x', 'y', 'height', 'band', 'geometry']].drop_duplicates()
            
            gs_shape = site_band.geometry.buffer(cls.SHAPE_RADIUS)
            gdf_shape = gpd.GeoDataFrame(geometry=gs_shape)

            gdf_shape = gpd.sjoin(gdf_shape, site_band, predicate='contains')
            gdf_shape = gdf_shape.drop(['index_right'], axis=1)[['bs_id', 'x', 'y', 'height', 'band', 'geometry']]
            geometries.append(gdf_shape)
        geometries = pd.concat(geometries).sort_values(by=['bs_id', 'band']).reset_index(drop=True)
        assert topo.cells[['bs_id', 'band']].drop_duplicates().shape[0] == geometries.shape[0]
        return geometries
    

class EmpiricalGeom(Geom):
    @property
    def ues(self):
        return self.topo.ues
    
    @override
    @classmethod
    def _create_site_shapes(cls, topo):
        return None
    
    @override
    @classmethod
    def _shape_radius_scale(cls, site_shapes, cell_shapes):
        return None


class ConvexGeom(EmpiricalGeom):
    @override
    @classmethod
    def _create_cell_shapes(cls, topo, site_shapes):
        cell_shapes = topo.ues.merge(topo.cells.drop(['geometry'], axis=1), on='cell_id')
        cell_shapes['geometry'] = topo.ues.geometry.convex_hull
            
        return cell_shapes


class ConcaveGeom(EmpiricalGeom):
    @override
    @classmethod
    def _create_cell_shapes(cls, topo, site_shapes):
        cell_shapes = topo.ues.merge(topo.cells.drop(['geometry'], axis=1), on='cell_id')
        cell_shapes['geometry'] = concave_hull(topo.ues.geometry, ratio=0)
            
        return cell_shapes
    

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