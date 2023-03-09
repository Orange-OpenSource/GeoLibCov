from datacov import GenMAPL, TopoGen

band_params = {
    'bands': [800, 2600],
    'band_probs': [1., 0.6],
    'ue_dist_scale': [0.08, 0.05],
    'ue_dist_loc': [0.1, 0.05]
}
topo = TopoGen(**band_params)
_, _ = topo.generate_topo()

topo_config = topo.cells.drop('geometry', axis=1).merge(topo.sites.drop('geometry', axis=1), on=['bs_id', 'x', 'y'])
mapl = GenMAPL(topo_config, n_prb=1)