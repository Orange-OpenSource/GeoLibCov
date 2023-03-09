import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
import geopandas as gpd
import pandas as pd
from shapely import MultiPoint
import matplotlib.patheffects as pe

class TopoGen:
    def __init__(self, n_sites=5, minlim=0.2, maxlim=0.8, 
                  bs_dist_loc=0.5, bs_dist_scale=0.15, bs_height=[10, 100],
                  ue_dist_loc=[0.1, 0.05], 
                  azimut_scale=10, n_azimuts=3,
                  n_ue=15,
                  ue_dist_scale=[0.05, 0.025],
                  bands=[800, 2600],
                  band_probs=[1., 0.5],
                ) -> None:
        
        """This is a function.

        Args:
            n_sites (int, optional): number of BS to generate . Defaults to 5.
            minlim (float, optional): minimum (X, Y)=(minlim, minlim) coordinates. Defaults to 0.2.
            maxlim (float, optional): maximum (X, Y)=(maxlim, maxlim) coordinates. Defaults to 0.8.
            bs_dist_loc (float, optional): mean of normal distribution of the BS locations. Defaults to 0.5.
            bs_dist_scale (float, optional): standard deviation of normal distribution of the BS locations. Defaults to 0.15.
            ue_avg_dist (float, optional): avg distance of UE from BS. Used in the generation of UE position. Defaults to 0.1.
            azimut_scale (int, optional): standard deviation of the azimut distribution. Defaults to 10.
        """
        self.n_sites = n_sites
        self.minlim = minlim
        self.maxlim = maxlim
        self.bs_dist_loc = bs_dist_loc
        self.bs_dist_scale = bs_dist_scale
        self.bs_height = bs_height
        self.azimut_scale = azimut_scale
        self.n_azimuts = n_azimuts
        self.n_ue = n_ue
        self.bands = bands
        self.ue_dist_scale = self.to_band_dict(ue_dist_scale)
        self.ue_dist_loc = self.to_band_dict(ue_dist_loc)
        self.band_probs = self.to_band_dict(band_probs)
        
        self.sites = None
        self.azimuts = None
        self.cells = None
        self.ues = None
        
    def to_band_dict(self, values):
        d = {}
        for b, v in zip(self.bands, values):
            d[b] = v
        return d
        
    def generate_topo(self):
        self.sites = self.generate_bs()
        self.azimuts = self.generate_azimuts(self.n_sites, self.azimut_scale, self.n_azimuts)
        df_az_geometry = self.create_az_lines()
        self.cells = self.generate_cells(df_az_geometry)
        gdf_ue = self.generate_ue_positions(df_az_geometry)
        self.ues = gdf_ue
        return self.cells, gdf_ue

    def generate_bs(self):
        x_bs = np.random.normal(loc=self.bs_dist_loc, scale=self.bs_dist_scale, size=self.n_sites)
        y_bs = np.random.normal(loc=self.bs_dist_loc, scale=self.bs_dist_scale, size=self.n_sites)
        height_bs = np.random.randint(low=self.bs_height[0], high=self.bs_height[1], size=self.n_sites)
        
        # rescale distribution to (minlim, minlim), (maxlim, maxlim) bounds
        x_std = (x_bs - x_bs.min()) / (x_bs.max() - x_bs.min())
        y_std = (y_bs - y_bs.min()) / (y_bs.max() - y_bs.min())
        x_bs = x_std * (self.maxlim-self.minlim) + self.minlim
        y_bs = y_std * (self.maxlim-self.minlim) + self.minlim

        bs_ids = [f'bs{i}' for i in range(self.n_sites)]
        gdf_bs = gpd.GeoDataFrame({'bs_id': bs_ids, 'x': x_bs, 'y': y_bs, 'height': height_bs}, geometry=gpd.points_from_xy(x_bs, y_bs))
        
        return gdf_bs

    def generate_cells(self, azimuts_geometry):
        bs_ids = self.sites.bs_id
        x_bs = self.sites.x
        y_bs = self.sites.y

        df_cells = []
        for b in self.bands:
            p = self.band_probs[b]
            df_cells.append(pd.DataFrame({'bs_id': bs_ids, 
                                'x': x_bs, 'y': y_bs, 
                                'band': b * np.random.binomial(1, p, self.n_sites)}))
        df_cells = pd.concat(df_cells)
        df_cells = df_cells[df_cells.band != 0]
        
        df_cells = azimuts_geometry.merge(df_cells, on=['bs_id', 'band'])
        df_cells['cell_id'] = df_cells.bs_id + '_' + df_cells.az_id + '_' + df_cells.band.astype(str)
        
        return df_cells[['cell_id', 'bs_id', 'az_id', 'band', 'azimut', 'x', 'y', 'geometry']]

    def plot_topography(self, frequency=None):
        x_bs = self.sites.x
        y_bs = self.sites.y

        colormaps = [
            'Greys', 
            'Purples', 
            'Blues', 
            'Greens', 
            'Oranges', 
            'Reds'
        ]
        
        band_marker = {}
        legend_marker = []
        markers = ['o', '^', '*', 'X', 's']
        for i, b in enumerate(self.bands):
            band_marker[b] = markers[i]
            legend_marker.append(mlines.Line2D([], [], color='grey', marker=markers[i], markersize=15, label=f'LTE{b}', linestyle='None'))
        
        site_color = {}
        for i, bs_id in enumerate(self.sites.bs_id.unique()):
            site_color[bs_id] = mpl.colormaps[colormaps[i % 6]].resampled(360)
        
        fig, ax = plt.subplots(figsize=(10, 10))
        max_loc = np.array(list(self.ue_dist_loc.values())).max()
        ax.set_xlim([self.minlim - max_loc, self.maxlim + max_loc])
        ax.set_ylim([self.minlim - max_loc, self.maxlim + max_loc])
        
        if frequency is not None:
            frequencies_to_plot = [frequency]
        else:
            frequencies_to_plot = self.bands
            
        site_cells = self.cells.copy()
        az_min = 0
        az_max = 359
        az_std = (site_cells['azimut'] - az_min) / (az_max - az_min)
        az_scaled = az_std * (0.9 - 0.1) + 0.1
        site_cells['azimut_color'] = az_scaled
            
        for i, bs_id in enumerate(self.sites.bs_id):
            c = site_color[bs_id]
            ax.scatter(x_bs[i], y_bs[i], marker='1', color=c(0.5))
            
            for b in frequencies_to_plot:
                cells_to_plot = site_cells[(site_cells.bs_id == bs_id) & (site_cells.band == b)].sort_values(by='cell_id')
                site_ues = self.ues.merge(cells_to_plot[['cell_id']], on='cell_id').sort_values(by='cell_id')
                azimuts = cells_to_plot.azimut_color
                az_geom = cells_to_plot.geometry
                ue_geom = site_ues.geometry
                for az, az_g, ue_g in zip(azimuts, az_geom, ue_geom):
                    x_az, y_az = az_g.xy
                    ax.plot(x_az, y_az, color=c(az), 
                            path_effects=[pe.Stroke(linewidth=2, foreground='grey'), pe.Normal()])
                    
                    x_ue = []
                    y_ue = []
                    for p in list(ue_g.geoms):
                        x_ue.append(p.x)
                        y_ue.append(p.y)
                    ax.scatter(x_ue, y_ue, color=c(az), marker=band_marker[b],
                               path_effects=[pe.Stroke(linewidth=2, foreground='grey'), pe.Normal()])
                    
            plt.legend(handles=legend_marker)
        plt.show()

    def generate_ue_positions(self, azimuts_geometry):
        df_cells = azimuts_geometry[['bs_id', 'band', 'azimut', 'x_ue_loc', 'y_ue_loc', 'ue_scale']].merge(self.cells, on=['bs_id', 'band', 'azimut'])
        cid = df_cells.cell_id
        x_ue_loc = df_cells.x_ue_loc
        y_ue_loc = df_cells.y_ue_loc
        ue_scale = df_cells.ue_scale
        
        gdf_ue = []
        ue_geometries = []
        for c, x_loc, y_loc, ue_s in zip(cid, x_ue_loc, y_ue_loc, ue_scale):
            x_ue = np.random.normal(x_loc, scale=ue_s, size=self.n_ue)
            y_ue = np.random.normal(y_loc, scale=ue_s, size=self.n_ue)
            gdf_ue.append({
                'cell_id': c,
            })
            pts = MultiPoint(np.vstack([x_ue, y_ue]).T)
            ue_geometries.append(pts)
        gdf_ue = pd.DataFrame(gdf_ue)
        gdf_ue = gpd.GeoDataFrame(gdf_ue, geometry=ue_geometries)
        assert gdf_ue.shape[0] == self.cells.shape[0]
        return gdf_ue
        
    def create_az_lines(self, azimuts=None, ue_dist_loc=None, use_cid=False):
        """creates geometry lines in the direction of azimuts.

        Args:
            azimuts (array, optional): (n_sites, n_azimuts). Defaults to None.
            ue_dist_loc (dict, optional): keys are frequency band, values are the length of the lines. Defaults to None.

        Returns:
            GeoDataFrame: geodataframe with geometry of azimut
        """
        bs_ids = self.sites.bs_id
        x_bs = self.sites.x
        y_bs = self.sites.y
        
        if azimuts is None:
            azimuts = self.azimuts
        
        if ue_dist_loc is None:
            ue_dist_loc = self.ue_dist_loc
        
        rad_azimuts = np.deg2rad(azimuts)

        df_azimuts = []
        geometries = []
        for b in self.bands:
            ue_loc = ue_dist_loc[b]
            ue_s = self.ue_dist_scale[b]
            for bs, x, y, rad_az_site, az_site in zip(bs_ids, x_bs, y_bs, rad_azimuts, azimuts):
                for az_id, (rad_az, az) in enumerate(zip(rad_az_site, az_site)):
                    x_az_end = x + ue_loc * np.cos(rad_az)
                    y_az_end = y + ue_loc * np.sin(rad_az)
                    wkt_line = f'LINESTRING({x} {y}, {x_az_end} {y_az_end})'
                    
                    df_azimuts.append({
                        'bs_id': bs, 
                        'band': b,
                        'az_id': f'az{az_id}',
                        'azimut': az,
                        'x_ue_loc': x_az_end,
                        'y_ue_loc': y_az_end,
                        'ue_scale': ue_s,})
                    geometries.append(wkt_line)
        df_azimuts = pd.DataFrame(df_azimuts)
        df_azimuts = gpd.GeoDataFrame(df_azimuts, geometry=gpd.GeoSeries.from_wkt(geometries))
        
        if self.cells is not None:
            if use_cid:
                df_azimuts = df_azimuts.merge(self.cells[['cell_id', 'bs_id', 'band', 'azimut']].drop_duplicates(), on=['bs_id', 'band', 'azimut'])
            else:
                df_azimuts = df_azimuts.merge(self.cells[['bs_id', 'band']].drop_duplicates(), on=['bs_id', 'band'])
        return df_azimuts

    def generate_azimuts(self, n_sites, azimut_scale, n_azimuts):
        azimuts = np.zeros((n_azimuts, n_sites))
        for az_idx in range(n_azimuts):
            if az_idx == 0:
                az = np.random.randint(low=0, high=360/(2*n_azimuts), size=n_sites)
            else:
                az_loc = azimuts[az_idx-1] + (360 / n_azimuts)
                az = np.round(np.random.normal(loc = az_loc, scale=azimut_scale, size=n_sites))
            azimuts[az_idx] = az
        azimuts = azimuts.T
        return azimuts

if __name__ == '__main__':
    gen = TopoGen()
    gdf = gen.generate_topo()
    gen.plot_topography()

