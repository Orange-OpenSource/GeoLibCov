import geopandas as gpd
import pandas as pd
import numpy as np
from datacov import Geom
from tqdm import tqdm

def precision_recall(coverage):
    ue_pts = coverage.topo.ues.merge(coverage.topo.cells.drop('geometry', axis=1), on='cell_id')
    ue_pts = ue_pts.explode(index_parts=False)
    total = ue_pts.groupby('cell_id').bs_id.count()

    overlaps = gpd.sjoin(coverage.cell_shapes, ue_pts, predicate='intersects')
    positives = overlaps[(overlaps.band_left == overlaps.band_right)].groupby('cell_id_left').bs_id_left.count()
    tp = overlaps[(overlaps.cell_id_left == overlaps.cell_id_right)].groupby('cell_id_left').bs_id_left.count()
    pr_stats = pd.DataFrame({'true_positives': tp, 'positives': positives, 'n_measurement': total})
    pr_stats['true_positives'] = pr_stats.true_positives.fillna(0)
    pr_stats['precision'] = pr_stats['true_positives'] / pr_stats['positives']
    pr_stats['recall'] = pr_stats['true_positives'] / pr_stats['n_measurement']
    return pr_stats

def get_pr_band(v_scaled):
    stats = precision_recall(v_scaled)
    band = stats.index.str.split('_').str[2].astype(int)
    
    stats['band'] = band
    p = stats.groupby('band').precision.mean().to_numpy()
    r = stats.groupby('band').recall.mean().to_numpy()
    return p, r


def precision_recall_site(coverage):
    vorosite_shape = coverage.cell_shapes[['bs_id', 'band', 'geometry']].drop_duplicates()
    ues = coverage.topo.ues.copy()
    ues = ues.explode(index_parts=False)
    ue_pts = ues.merge(coverage.cells.drop('geometry', axis=1), on='cell_id')[['bs_id', 'band', 'geometry']]
    total = ue_pts.groupby(['bs_id', 'band']).geometry.count().reset_index()

    overlaps = gpd.sjoin(vorosite_shape, ue_pts, predicate='intersects')
    positives = overlaps[(overlaps.band_left == overlaps.band_right)].groupby(['bs_id_left', 'band_left']).geometry.count().reset_index(drop=False)

    tp = overlaps[(overlaps.bs_id_left == overlaps.bs_id_right) & (overlaps.band_left == overlaps.band_right)].groupby(['bs_id_left', 'band_left']).geometry.count().reset_index(drop=False)
    pr_stats = pd.DataFrame({'bs_id': tp.bs_id_left, 'band':tp.band_left, 'true_positives': tp.geometry, 'positives': positives.geometry, 'n_measurement': total.geometry})
    pr_stats['true_positives'] = pr_stats.true_positives.fillna(0)
    pr_stats['precision'] = pr_stats['true_positives'] / pr_stats['positives']
    pr_stats['recall'] = pr_stats['true_positives'] / pr_stats['n_measurement']
    return pr_stats

def get_pr_band_site(coverage):
    stats = precision_recall_site(coverage)
    
    p = stats.groupby('band').precision.mean().to_numpy()
    r = stats.groupby('band').recall.mean().to_numpy()
    return p, r

def evaluate(coverage, scales=np.arange(0.3, 2.1, 0.1), pr_func=get_pr_band):
    p = np.zeros((len(scales), len(coverage.topo.bands)))
    r = np.zeros((len(scales), len(coverage.topo.bands)))

    for j, s in enumerate(tqdm(scales)):
        v_scaled = Geom().scale_cell_shapes(coverage, sfact=s)
    
        p[j], r[j] = pr_func(v_scaled)
    return p.T, r.T

def retrieve_scale_factors(r_l_model, coverage):
    df_scale = r_l_model.data.merge(coverage.radius_scale, on='cell_id')
    df_scale['scale_factor'] = df_scale['r_l'] / df_scale['model_radius']
    return list(df_scale['scale_factor'])

def evaluate_mapl_models(n_prb_range, mapl_model, topo_mapl, mapl_context, rev_model, coverage_shape):
    p = np.zeros((len(n_prb_range), len(coverage_shape.topo.bands)))
    r = np.zeros((len(n_prb_range), len(coverage_shape.topo.bands)))
    
    for i, n_prb in enumerate(tqdm(n_prb_range)):
        mapl = mapl_model(topo_mapl, n_prb=n_prb, context_array=mapl_context)
        radius_model = rev_model(mapl)
        rescale_factors = retrieve_scale_factors(radius_model, coverage_shape)
        scaled_shapes = Geom().scale_cell_shapes(coverage_shape, rescale_factors)
        p[i], r[i] = get_pr_band(scaled_shapes)
    return p.T, r.T

def tune_scale(dist_gen, coverage, scales=np.arange(0.1, 2.1, 0.1)):
    errors = []
    bands = coverage.cell_shapes.band.unique()
    for s in scales:
        scaled_cov = Geom().scale_cell_shapes(coverage, s)
        cmp = dist_gen.distances.merge(scaled_cov.radius_scale, on='cell_id')
        cmp['ae'] = np.abs(cmp['proxy_radius'] - cmp['model_radius'])
        mae = dict(cmp.groupby('band').ae.mean())
        mae['scale'] = s
        errors.append(mae)
    errors = pd.DataFrame(errors)
    
    opt_scale = []
    for b in bands:
        opt_scale.append({'band': b, 'scale': scales[np.argmin(errors[b])], 'scale_idx': np.argmin(errors[b]), 'mae': np.min(errors[b])})
    return pd.DataFrame(opt_scale)