import numpy
import pandas
import connalysis

import tqdm

def get_simplex_dataframes(M, dimensions):
    simplices = connalysis.network.topology.list_simplices_by_dimension(M.matrix)

    dfs = [pandas.DataFrame(M.gids[_smpl])
           for _smpl in simplices.loc[dimensions]]
    for _df in dfs: _df.index.name = "Simplex"
    ugids = [numpy.unique(_df.values) for _df in dfs]
    ugid_lsts = [list(_ugids) for _ugids in ugids]
    df_idxs = [
        _df.map(lambda _x: _ugid_lst.index(_x))
        for _df, _ugid_lst in tqdm.tqdm(list(zip(dfs, ugid_lsts)))
    ]
    return simplices, dfs, df_idxs, ugids


def get_divergences(simplices):
    return simplices.apply(lambda _smpl: [len(numpy.unique(_x)) / len(_x) for _x in _smpl.transpose()])


def exc_inh_sparce_matrices(M, ugids_smpl, ugids_nrn):
    in_nrn = numpy.in1d(M.gids, ugids_nrn)
    in_smpl = numpy.in1d(M.gids, ugids_smpl)

    assert (M.gids[in_nrn] == ugids_nrn).all()
    assert (M.gids[in_smpl] == ugids_smpl).all()

    m = M.matrix.tocsr()
    m_nrn_smpl = m[numpy.ix_(in_nrn,
                             in_smpl)]
    m_smpl_nrn = m[numpy.ix_(in_smpl,
                             in_nrn)]
    return m_nrn_smpl, m_smpl_nrn


def get_simplex_neuron_path_df(m_nrn_smpl, m_smpl_nrn, df_idxs, str_neuron_population):
    def to_adjacency_series_1(subm):
        subm = subm.tocoo()
        idx = pandas.MultiIndex.from_frame(pandas.DataFrame({"Position": subm.row, str_neuron_population: subm.col}))
        return pandas.Series(subm.data, index=idx)

    def to_adjacency_series_2(subm):
        subm = subm.tocoo()
        idx = pandas.MultiIndex.from_frame(pandas.DataFrame({"Position": subm.col, str_neuron_population: subm.row}))
        return pandas.Series(subm.data, index=idx)

    smplx_nrn_paths = df_idxs.apply(lambda _row: to_adjacency_series_1(m_smpl_nrn[_row]), axis=1)
    nrn_smplx_paths = df_idxs.apply(lambda _row: to_adjacency_series_2(m_nrn_smpl[:, _row]), axis=1)

    return smplx_nrn_paths.stack(["Position", str_neuron_population]), nrn_smplx_paths.stack(["Position", str_neuron_population])


def simplex_specific_inout_degrees(s_n_paths, n_s_paths, str_neuron_population):
    return pandas.concat([s_n_paths.groupby(str_neuron_population).count(),
                          n_s_paths.groupby(str_neuron_population).count()],
                          keys=["In", "Out"], axis=1).fillna(0)


def get_disynaptic_path_sum(s_n_path, n_s_path, tgt_dim, use_weight=False):
    mulA = s_n_path.unstack("Position", fill_value=0).reindex(columns=numpy.arange(tgt_dim + 1))
    mulB = n_s_path.unstack("Position", fill_value=0).reindex(columns=numpy.arange(tgt_dim + 1))
    reidx = mulB.index.union(mulA.index)
    mulA = mulA.reindex(reidx, fill_value=0).values
    mulB = mulB.reindex(reidx, fill_value=0).values

    if not use_weight:
        mulA = (mulA > 0).astype(int)
        mulB = (mulB > 0).astype(int)

    return numpy.dot(mulA.transpose(),
                     mulB)

def sum_disynaptic_path_sum(s_n_paths, n_s_paths, tgt_dim, str_neuron_population,
                            min_degree={}, use_weight=False):
    import tqdm

    m_path_count = numpy.zeros((tgt_dim + 1, tgt_dim + 1))
    smpl_deg = simplex_specific_inout_degrees(s_n_paths, n_s_paths, str_neuron_population)
    for k, v in min_degree.items():
        smpl_deg = smpl_deg.loc[smpl_deg[k] >= v]

    _s_n = s_n_paths.reorder_levels([2, 1, 0])
    _n_s = n_s_paths.reorder_levels([2, 1, 0])

    for i in tqdm.tqdm(smpl_deg.index.values):
        m_path_count += get_disynaptic_path_sum(_s_n[i], _n_s[i], tgt_dim, use_weight=use_weight)

    return m_path_count / len(smpl_deg.index)

def to_coo_matrix(c_df, tgt_shape, str_neuron_population, shuffle_cols=[]):
    from scipy import sparse

    for col in shuffle_cols:
        c_df[col] = numpy.random.permutation(c_df[col].values)

    return sparse.coo_matrix((c_df[0].values,
                              (c_df["Simplex"].values, c_df[str_neuron_population].values)),
                             shape=tgt_shape)


def normalize_coo_matrix(mtrx, axis):
    t_sums = numpy.array(mtrx.mean(axis=axis)).flatten()
    attrs = ["col", "row"]
    mtrx.data = mtrx.data / t_sums[getattr(mtrx, attrs[axis])]


def get_disynaptic_con_mat(s_n_paths, n_s_paths, str_neuron_population, tgt_shape,
                           use_weight=False, normalize=False, shuffle_cols=[]):
    if use_weight:
        m_smpl_i = to_coo_matrix(s_n_paths.groupby(["Simplex", str_neuron_population]).sum().reset_index(),
                                 tgt_shape, str_neuron_population, shuffle_cols=shuffle_cols)
        m_i_smpl = to_coo_matrix(n_s_paths.groupby(["Simplex", str_neuron_population]).sum().reset_index(),
                                 tgt_shape, str_neuron_population, shuffle_cols=shuffle_cols).transpose()
    else:
        m_smpl_i = to_coo_matrix(s_n_paths.groupby(["Simplex", str_neuron_population]).count().reset_index(),
                                 tgt_shape, str_neuron_population, shuffle_cols=shuffle_cols)
        m_i_smpl = to_coo_matrix(n_s_paths.groupby(["Simplex", str_neuron_population]).count().reset_index(),
                                 tgt_shape, str_neuron_population, shuffle_cols=shuffle_cols).transpose()

    if normalize:
        normalize_coo_matrix(m_smpl_i, 1)
        normalize_coo_matrix(m_i_smpl, 0)

    m_smpl_smpl = numpy.array((m_smpl_i * m_i_smpl).todense())

    if normalize: m_smpl_smpl = numpy.sqrt(m_smpl_smpl)
    return m_smpl_smpl


def simplex_overlap_in_positions(data_df, pos_idx, ugid):
    from scipy import sparse
    use_idx = data_df[pos_idx]
    shape = (len(data_df), len(ugid))

    mbr_mat = sparse.csr_matrix((numpy.ones(numpy.prod(use_idx.shape), dtype=int),
                                use_idx.values.flatten(),
                                numpy.arange(len(use_idx) + 1, dtype=int) * use_idx.shape[1]),
                                shape=shape)
    assert numpy.all(numpy.array(mbr_mat.sum(axis=1)).flatten() == use_idx.shape[1])
    return numpy.array((mbr_mat * mbr_mat.transpose()).todense())


def overlaps_vs_disyn(ol_src, ol_tgt, disyn_nrn_mat):
    def offdiagonal(mat):
        return mat[numpy.eye(mat.shape[0]) == 0]
    
    master_df = pandas.DataFrame({"src": offdiagonal(ol_src),
                                  "tgt": offdiagonal(ol_tgt),
                                  "inh": offdiagonal(disyn_nrn_mat)})
    ol_disyn_mat = master_df.groupby(["src", "tgt"])["inh"].mean().unstack("tgt")
    ol_disyn_cc = numpy.corrcoef(master_df.values, rowvar=False)
    return ol_disyn_mat, ol_disyn_cc, master_df.columns

def merge_small_clusters(grp_df, frac_merge_clst):
    grp_df = grp_df.copy()
    for col in grp_df.columns:
        lbls_counts = grp_df[col].value_counts()
        lbls_counts = lbls_counts < (lbls_counts[0] * frac_merge_clst)
        lbls_merge = lbls_counts[lbls_counts].index.values

        grp_df[col] = numpy.where(numpy.isin(grp_df[col].values, lbls_merge),
                                        numpy.max(lbls_counts.index) + 1,
                                        grp_df[col].values)
    return grp_df

def simplex_clustering(ol_src, ol_tgt, clst_param):
    from sknetwork.clustering import Louvain

    lbls_src = Louvain(clst_param).fit_predict(ol_src)
    print(len(numpy.unique(lbls_src)))
    lbls_tgt = Louvain(clst_param).fit_predict(ol_tgt)
    print(len(numpy.unique(lbls_tgt)))
    return pandas.DataFrame({"src_grp": lbls_src, "tgt_grp": lbls_tgt})


def simplex_cluster_connectivity(grp_df, normalization):
    values = grp_df.value_counts().sort_index().unstack("tgt_grp", fill_value=0)

    if normalization == "None":
        pass
    elif normalization == "target":
        values = values / grp_df["tgt_grp"].value_counts().sort_index()
    elif normalization == "source":
        values = (values.transpose() / grp_df["src_grp"].value_counts().sort_index()).transpose()
    elif normalization == "pairs":
        pairs = grp_df["src_grp"].value_counts().sort_index().values.reshape((-1, 1)) * grp_df["tgt_grp"].value_counts().sort_index().values.reshape((1, -1))
        values = values / pairs
    return values

def simplex_cluster_disyn(disyn_mat, grp_df, normalization, cutoff=None):
    lbls_src = grp_df["src_grp"]
    lbls_tgt = grp_df["tgt_grp"]
    
    if cutoff is not None:
        disyn_mat = disyn_mat >= cutoff

    values_inh = numpy.array([[(disyn_mat[numpy.ix_(lbls_src == i, lbls_tgt == j)]).sum()
                                for j in numpy.unique(lbls_tgt)]
                                for i in numpy.unique(lbls_src)])
    if normalization == "None":
        pass
    elif normalization == "target":
        values_inh = values_inh / grp_df["tgt_grp"].value_counts().sort_index().values
    elif normalization == "source":
        values_inh = values_inh / grp_df["src_grp"].value_counts().sort_index().values.reshape((-1, 1))
    elif normalization == "pairs":
        pairs = grp_df["src_grp"].value_counts().sort_index().values.reshape((-1, 1)) * grp_df["tgt_grp"].value_counts().sort_index().values.reshape((1, -1))
        values_inh = values_inh / pairs
    return values_inh



