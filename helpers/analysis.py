import numpy
import pandas
import connalysis

import tqdm

def compress_id_ranges(dfs):
    """
    Reindexes values in a DataFrame, compressing them to have [0 ... n] as values, where
    n is the number of unique values in the DataFrame.

    Args:
      dfs: list of DataFrames with integer values

    Returns: 
      df_idxs: list of DataFrames with values adjusted as described
      ugids: list of arrays, each holding the unique values in the corresponding DataFrame.
      That is, dfs[i][col][j] == ugids[i][df_idxs[i][col][j]]
    """
    ugids = [numpy.unique(_df.values) for _df in dfs]
    los = [pandas.Series(_gids, name="gid").reset_index().set_index("gid")["index"]
       for _gids in ugids]
    df_idxs = [
        _df.apply(lambda _col: _lo[_col].values, axis=0)
        for _df, _lo in tqdm.tqdm(list(zip(dfs, los)))
    ]
    return df_idxs, ugids

def get_simplex_dataframes(M, dimensions):
    """
    Get the list of simplices in a ConnectivityMatrix for given dimensions.

    Args:
      M: a conntility.ConnectivityMatrix object
      dimensions: list of dimensions to consider (list of ints)
    
    Returns:
      simplices: The output of connalysis.network.topology.list_simplices_by_dimension
      dfs: list of DataFrames. One entry per dimension. Each entry is a n x dim + 1
        DataFrame, where dim is the corresponding dimension and n the number of simplices in
        that dimension. Each entry specifies the "gid" of the node at the corresponding position
        of the corresponding simplex..
      df_idxs: Like dfs, but instead of returning the node "gid"s it returns the indices into
      the sorted list of unique "gid"s.
      ugids: Lists of sorted "gid"s.

    See also: compress_id_range  
    """
    simplices = connalysis.network.topology.list_simplices_by_dimension(M.matrix,threads=20)

    dfs = [pandas.DataFrame(M.gids[_smpl])
           for _smpl in simplices.loc[dimensions]]
    for _df in dfs: _df.index.name = "Simplex"
    df_idxs, ugids = compress_id_ranges(dfs)
    return simplices, dfs, df_idxs, ugids


def get_divergences(simplices):
    """
    From a simplex list calculate the simplex divergence. That is, the relative number of unique node ids
    in each simplex position.

    Args:
      simplices: n x dim + 1 DataFrame of node ids of simplices, where n is the number of simplices and dim
        the considered dimension.
    
    Returns:
      A pandas.Series of length dim + 1 with the relative number of unique nodes in each simplex position.
    """
    return simplices.apply(lambda _smpl: [len(numpy.unique(_x)) / len(_x) for _x in _smpl.transpose()])


def exc_inh_sparce_matrices(M, ugids_smpl, ugids_nrn):
    """
    Returns the connectivity matrices from the "simplex" population to the "neuron" population and the other
    way around.

    Args:
      M: A conntility.ConnectivityMatrix.
      ugids_smpl: The "gid"s of nodes in the "simplex" population.
      ugids_nrn: The "gid"s of nodes in the "neuron" population.

    Returns:
      m_nrn_smpl: A sparse matrix that is the adjacency matrix of connectivity from the "neuron" to the "simplex"
        population.
      m_smpl_nrn: A sparse matrix that is the adjacency matrix of connectivity from the "simplex" to the "neuron"
        population.
    """
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
    """
    Get the DataFrames of disynaptic paths from a simplex via an (inhibitory) neuron to another simplex.
    This is essentially a 3-dimensional sparse connectivity matrix between simplices and (inhibitory) neurons.
    That is, a list of triplets, (i, j, k) such that entry is 1 if there is a connections from the neuron in
    position j of simplex i to (inhibitory) neuron k; and another for the other direction (from neurons to simplices).

    Args:
      m_nrn_smpl: sparse connectivity matrix from the "neuron" population to the "simplex" population. 
      m_smpl_nrn: sparse connectivity matrix from the "simplex" population to the "neuron" population.
        Outputs of exc_inh_sparce_matrices
      df_idxs: DataFrame of nodes in simplices (of a given dimension). Output of  get_simplex_dataframes
      str_neuron_population: A string denoting the name of the column corresponding to the "neuron" population
        in the output.

    Returns:
      Two DataFrames with MultiIndices. MultiIndex levels are called "Simplex", "Position" and str_neuron_population.
      Entries exist if there is a connection from the node at the indicated position of
      the indicated simplex to the indicated node; in m_smpl_nrn (first output) or m_nrn_smpl (second output).
      Values of the DataFrame are the values in the sparse matrices for the corresponding connection.
    """
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
    """
    Counts the number of simplex neurons innervated by or innervating a given (inhibitory) neuron. 

    Args:
      s_n_paths:
      n_s_paths: Outputs of get_simplex_neuron_path_df. See that function.
      str_neuron_population: Input used for get_simplex_neuron_path_df. See that function.

    Returns:
      pandas.DataFrame; one row per (inhibitory) neuron. Two columns "In" and "Out" denoting the simplex
      in- and out-degree respectively.
    """
    return pandas.concat([s_n_paths.groupby(str_neuron_population).count(),
                          n_s_paths.groupby(str_neuron_population).count()],
                          keys=["In", "Out"], axis=1).fillna(0)


def get_disynaptic_path_sum(s_n_path, n_s_path, tgt_dim, use_weight=False):
    """
    Count number of disynaptic paths between nodes in each position of a given simplex via an (inhibitory) 
    neuron. 

    Args:
      s_n_path: pandas.DataFrame specifying connections from simplices to (inhibitory) neurons. MultiIndex
        specifies simplex position and (inhibitory) neuron identifier. Entries exist if a connection exists
        from simplex node to (inhibitory) neuron.
      n_s_path: As s_n_path, but from neurons to simplices. See also: get_simplex_neuron_path_df.
      tgt_dim: The dimension of the simplices considered.
      use_weight (optional; default: False): If False, counts number of paths; if True, sum over weights
      of existing paths.

    Returns:
      numpy.matrix, of shape n x n, where n is the number of simplices.
    """
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
    """
    As get_disynaptic_path_sum, but counting paths for all simplices, not just a single one.
    Output is normalized by the number of (inhibitory) neurons.
    """
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
    """
    Convert a DataFrame of the connections from simplices to neurons to a sparse.coo_matrix.

    Args:
      c_df: Simplex to neuron connectivity DataFrame. See get_simplex_neuron_path_df
      tgt_shape: The shape of the output matrix (see below).
      str_neuron_population: Name of the neuron population that was used as the input of 
        get_simplex_neuron_path_df.
      shuffle_cols (optional): Names of the columns to randomly shuffle to create a control.
        Can be any combination of: "Simplex" and / or str_neuron_population. Default: No shuffle.

    Returns:
      n x m sparse matrix, where n is the number of simplices and m the number of neurons in the
      specified (inhibitory) population.
    """
    from scipy import sparse

    for col in shuffle_cols:
        c_df[col] = numpy.random.permutation(c_df[col].values)

    return sparse.coo_matrix((c_df[0].values,
                              (c_df["Simplex"].values, c_df[str_neuron_population].values)),
                             shape=tgt_shape)


def normalize_coo_matrix(mtrx, axis):
    """
    Normalize values of a coo_matrix object, such that the sum over a specified axis is 1

    Args:
      mtrx: a coo_matrix object.
      axis: Either 0 or 1, specifying the axis over which to normalize.

    Returns: None, the coo_matrix is normalized in-place.
    """
    t_sums = numpy.array(mtrx.mean(axis=axis)).flatten()
    attrs = ["col", "row"]
    mtrx.data = mtrx.data / t_sums[getattr(mtrx, attrs[axis])]


def get_disynaptic_con_mat(s_n_paths, n_s_paths, str_neuron_population, tgt_shape,
                           use_weight=False, normalize=False, shuffle_cols=[]):
    """
    Get matrix of the number of disynaptic paths from and to each simplex.

    Args:
      s_n_paths:
      n_s_paths: Outputs of get_simplex_neuron_path_df. See that function.
      str_neuron_population: Input used for get_simplex_neuron_path_df. See that function.
      tgt_shape: Shape of the output matrix (see below).
      use_weight (optional, default: False): If False, count number of paths; if True, sum over
        weights of paths.
      normalize: (optional, default: False): If True, normalize with respect to the total number
        of paths from / to each simplex
      shuffle_cols (optional, defalt: no shuffle): See to_coo_matrix.

    Returns:
      A n x n array, where n is the number of simplices. 
    """
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


def histogram_data_and_ctrl(m_data, m_ctrl, bins):
    """
    Create histograms of values in a "data" and a "control" distribution.

    Args:
      m_data: array with "data" values.
      m_ctrl: array with "control" values.
      bins: The binning to use

    Returns:
      A DataFrame with two columns. One with histogram counts for the "Data" and for the "Control"
      distribution
    """
    if not hasattr(bins, "__iter__"):
        bins = numpy.linspace(0, m_data.max(), bins)
    H_data = numpy.histogram(m_data.flatten(), bins=bins)[0]
    H_ctrl = numpy.histogram(m_ctrl.flatten(), bins=bins)[0]
    return pandas.DataFrame({"Data": H_data, "Control": H_ctrl}, index=bins[:-1])


def simplex_overlap_in_positions(data_df, pos_idx, ugid):
    """
    Calculates the number of overlapping nodes for each pair of simplices in the specified positions.

    Args:
      data_df: DataFrame of node "gid"s of simplices
      pos_idx: Positions to consider the overlap in. I.e., [0, 1, 2] would calculate the overlap in the
        first three simplex positions.
      ugid: List of unique node "gid"s in simplices.

    Returns:
      n x n array specifying the number of overlapping nodes for each pair of simplices. 
    """
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
    """
    Calculates correlation between simplex overlaps at the source and target side and the strength
    of disynaptic inhibition.

    Args:
      ol_src: n x n array with integer values. (Meant to specify overlap at the source side for each
        pair of n simplices)
      ol_tgt: n x n array with integer values. (Meant to specify overlap at the target side for each
        pair of n simplices)
      disyn_nrn_mat: n x n array. (Meant to specify strength of disynaptic inhibition between each
        pair of n simplices)

    Returns:
      ol_disyn_mat: a x b array, where a is the number of unique values in ol_src and b is the number
        of unique values in ol_tgt. For each combination of unique values the mean value of 
        disyn_nrn_mat where those values are obtained.
      ol_disyn_cc: Correlation between the values of the three inputs. Shape 3 x 3. Output of numpy.corrcoef
      columns: List of string specifying the order of rows / columns in ol_disyn_cc.
    """
    def offdiagonal(mat):
        return mat[numpy.eye(mat.shape[0]) == 0]
    
    master_df = pandas.DataFrame({"src": offdiagonal(ol_src),
                                  "tgt": offdiagonal(ol_tgt),
                                  "inh": offdiagonal(disyn_nrn_mat)})
    ol_disyn_mat = master_df.groupby(["src", "tgt"])["inh"].mean().unstack("tgt")
    ol_disyn_cc = numpy.corrcoef(master_df.values, rowvar=False)
    return ol_disyn_mat, ol_disyn_cc, master_df.columns

def merge_small_clusters(grp_df, frac_merge_clst):
    """
    Postprocessing of the output of a clustering run. All clusters smaller than a cutoff
    are merged into a single "dump" cluster. The id of that "dump" cluster will be larger
    than that of all proper clusters, placing it at the end.

    Args:
      grp_df: Raw clustering results. DataFrame with one clustering output in each column.
        All columns will be independently post-processed.
      frac_merge_clst: Determines the cutoff for small clusters. Clusters smaller than this
        fraction of the largest cluster will be merged.
    """
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
    """
    Generates clusters of simplices from their overlaps.

    Args:
      ol_src: array of overlaps on the source side. Shape n x n, where n is the number of simplices.
      ol_tgt: array of overlaps on the target side.
      clst_param: Cluster parameter given to the Louvain clusterer.

    Returns:
      DataFrame with cluster assignments. Two columns, one based on source side overlaps, the other 
      based on target side overlap.
    """
    from sknetwork.clustering import Louvain

    lbls_src = Louvain(clst_param).fit_predict(ol_src)
    print(len(numpy.unique(lbls_src)))
    lbls_tgt = Louvain(clst_param).fit_predict(ol_tgt)
    print(len(numpy.unique(lbls_tgt)))
    return pandas.DataFrame({"src_grp": lbls_src, "tgt_grp": lbls_tgt})


def simplex_cluster_connectivity(data_df, grp_df, ugid, normalization):
    """
    Calculates the average total overlap between pairs of simplices in simplex groups.
    Since nodes in a simplex are highly connected, this is a proxy of excitatory connection
    strengths between simplex groups. 

    Args:
      data_df: DataFrame of node "gid"s of simplices
      grp_df: DataFrame of simplex clustering results. See simplex_clustering.
      ugid: list of unique "gid"s of nodes in simplices. 
      normalization: Specify how the results are normalized. One of: "None", "target", "source",
        or "pairs".

    Returns: 
      Array of total overlap strength between simplex clusters. Shape is u x v, where u is
      the number of source clusters, and v is the number of target clusters.
    """
    full_overlap = simplex_overlap_in_positions(data_df, numpy.arange(data_df.shape[1]), ugid)
    return simplex_cluster_disyn(full_overlap, grp_df, normalization, cutoff=None)

def simplex_cluster_disyn(disyn_mat, grp_df, normalization, cutoff=None):
    """
    Calculates total strength of connectivity between simplex groups.

    Args:
      disyn_mat: Array of connections strengths between simplices. E.g., the output of get_disynaptic_con_mat.
      grp_df: DataFrame of simplex clustering results. See simplex_clustering.
      normalization: Specify how the results are normalized. One of: "None", "target", "source",
        or "pairs".
      cutoff (optional, default: None): If specified, values of disyn_mat below cutoff are ignored.

    Returns: 
      Array of connection strength (according to disyn_mat) between simplex clusters. Shape is u x v, 
      where u is the number of source clusters, and v is the number of target clusters.
    """
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


def edge_participation_df(M, max_dim=None):
    """
    Calculates edge partipation of each edge in a ConnectivityMatrix.

    Args:
      M: conntility.ConnectivityMatrix object.
      max_dim: highest simplex dimension to consider.

    Returns:
      DataFrame with one column per simplex dimension. Index is a MultiIndex with two levels:
      "row" and "col" that specify the index of the source and target node of the corresponding
      edge. The order of rows (i.e. edges) matches their order in M.
    """
    ep = connalysis.network.topology.edge_participation(M.matrix, threads=20, max_dim=max_dim)
    idx = pandas.DataFrame(numpy.vstack(ep.index), columns=["row", "col"])
    ep.index = pandas.MultiIndex.from_frame(idx)
    ep = ep.loc[pandas.MultiIndex.from_frame(M._edge_indices)]
    return ep

