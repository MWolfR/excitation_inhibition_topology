import numpy

from matplotlib import pyplot as plt

"""
Various plots used in the notebooks.
"""

def plot_simplex_divergences(divergences):
    """
    Fig. 2B
    """
    fig = plt.figure(figsize=(7, 2.5))

    ax = fig.add_subplot(1, 2, 1)
    for dim, vals in divergences.items():
        if dim > 0:
            ax.plot(numpy.linspace(0, 1, len(vals)), vals, marker='o', label="Dim: {0}".format(dim))
    ax.set_xticks([0, 1]); plt.gca().set_xticklabels(["Source", "Sink"])
    ax.set_ylabel("Frac. unique nodes")
    plt.legend()

    ax = fig.add_subplot(1, 2, 2)
    for dim, vals in divergences.items():
        if dim > 0:
            ax.plot(numpy.linspace(0, 1, len(vals)), numpy.array(vals) / vals[0], marker='o', label="Dim: {0}".format(dim))
    ax.set_xticks([0, 1]); plt.gca().set_xticklabels(["Source", "Sink"])
    ax.set_ylabel("Rel. Frac. unique nodes")
    plt.legend()
    return fig


def indegree_outdegree_figure(indeg, outdeg, nbins):
    """
    Not used in final version.
    """
    from scipy.interpolate import RegularGridInterpolator
    fig = plt.figure(figsize=(4, 4))
    
    xbins = numpy.linspace(0, indeg.max() + 1, nbins)
    ybins = numpy.linspace(0, outdeg.max() + 1, nbins)
    H = numpy.histogram2d(indeg, outdeg, bins=(xbins, ybins))[0]
    
    ip = RegularGridInterpolator((xbins[:-1], ybins[:-1]), H, bounds_error=False, fill_value=1)
    c = ip((indeg, outdeg))
    
    ax = fig.add_axes((0.1, 0.1, 0.65, 0.65))
    ax.scatter(indeg, outdeg, c=numpy.sqrt(c), s=2)
    ax.set_frame_on(False)
    ax.set_xlim([xbins[0], xbins[-1]])
    ax.set_ylim([ybins[0], ybins[-1]])
    ax.set_xlabel("Indegree")
    ax.set_ylabel("Outdegree")
    
    ax = fig.add_axes((0.1, 0.8, 0.65, 0.15))
    ax.plot(xbins[:-1], H.sum(axis=1) + 1, color="black")
    ax.set_xlim([xbins[0], xbins[-1]])
    ax.set_frame_on(False)
    ax.set_xticks([])
    ax.set_yscale("log")
    
    ax = fig.add_axes((0.8, 0.1, 0.15, 0.65))
    ax.plot(H.sum(axis=0) + 1, ybins[:-1], color="black")
    ax.set_ylim([ybins[0], ybins[-1]])
    ax.set_frame_on(False)
    ax.set_yticks([])
    ax.set_xscale("log")

    return fig

def plot_disynaptic_path_sum(mat_disyn):
    """
    Not used in final version.
    """
    ticks = [0, (mat_disyn.shape[0] - 1) / 2, mat_disyn.shape[0] - 1]
    fig = plt.figure(figsize=(3.5, 3.5))
    ax = fig.gca()

    img = ax.imshow(mat_disyn) # Number of E-I-E paths where both E are in the same simpl. in indicated positions.
    plt.colorbar(img)
    ax.set_xlabel("I -> E position")
    ax.set_ylabel("E -> I position")
    ax.set_xticks(ticks); ax.set_xticklabels(["Source", "->", "Sink"])
    ax.set_yticks(ticks); ax.set_yticklabels(["Source", "->", "Sink"])
    return fig

def plot_position_degrees(paths_df):
    """
    Fig. 3C
    """
    fig = plt.figure(figsize=(3.5, 3.5))
    ax = fig.gca()
    I = paths_df.groupby(["Simplex", "Position"]).count().unstack("Position", fill_value=0).values
    idxx = numpy.argsort(I.sum(axis=1)); I = I[idxx]
    img = ax.imshow(I, aspect="auto", interpolation="nearest")
    plt.colorbar(img)
    tcks = [0, (I.shape[1] - 1) / 2, I.shape[1] - 1]
    ax.set_xticks(tcks); ax.set_xticklabels(["Source", "->", "Target"])
    ax.set_ylabel("Simplex")
    ax.set_xlabel("Position")
    return fig

def summary_position_degrees(all_classic, all_smpl, labels=("", ""), add_title=False):
    """
    Related to Fig. 3D. But plotting only results for a single subnetwork.
    """
    fig = plt.figure(figsize=(2.5, 3.5))
    ax = fig.add_subplot(2, 1, 1)
    cols = ["red", "blue"]
    offset = [-0.01, 0.01]

    for data, col, o, lbl in zip([all_classic["In"], all_smpl["In"]], cols, offset, labels):
        mn = data.mean(axis=0)
        sd = data.std(axis=0)
        x = numpy.linspace(0, 1, len(mn)) + o
        ax.errorbar(x, mn, yerr=sd, color=col, marker='o', label=lbl)
    ax.set_frame_on(False)
    ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["Source", "->", "Sink"])
    if add_title:
        ax.set_title("In-degrees")

    ax = fig.add_subplot(2, 1, 2)

    for data, col, o, lbl in zip([all_classic["Out"], all_smpl["Out"]], cols, offset, labels):
        mn = data.mean(axis=0)
        sd = data.std(axis=0)
        x = numpy.linspace(0, 1, len(mn)) + o
        ax.errorbar(x, mn, yerr=sd, color=col, marker='o', label=lbl)
    ax.set_frame_on(False)
    ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["Source", "->", "Sink"])
    if add_title:
        ax.set_title("Out-degrees")
    return fig
    
def plot_position_mean(paths_df):
    """
    Suppl. Fig. S1
    """
    fig = plt.figure(figsize=(4.5, 3.5))
    ax = fig.gca()
    plt.colorbar(ax.imshow(paths_df.groupby(["Simplex", "Position"]).mean().unstack("Position", fill_value=0),
                           aspect="auto", interpolation="nearest"))
    ax.set_ylabel("Simplex")
    ax.set_xlabel("Position")
    return fig

def compare_disyn_inhibition(hists):
    """
    Fig. 4C
    """
    fig = plt.figure(figsize=(3.5, 3))
    ax = fig.gca()

    for col in hists.columns:
        ax.plot(hists[col], label=col)
    
    ax.set_xlabel("Disyn. connection strength")
    ax.set_ylabel("Simplex pairs")
    plt.legend()
    return fig
    
def plot_overlaps_vs_disyn(ol_disyn_mat):
    """
    Fig. 4F
    """
    fig = plt.figure(figsize=(3, 3))
    ax_main = fig.add_axes([0.1, 0.1, 0.6, 0.6])
    ax_tgt = fig.add_axes([0.1, 0.75, 0.6, 0.2])
    ax_src = fig.add_axes([0.75, 0.1, 0.2, 0.6])

    ax_main.imshow(ol_disyn_mat)
    ax_main.set_xticks([0, 1, 2, 3]); ax_main.set_yticks([0, 1, 2, 3])
    ax_main.set_xlabel("Target side overlap"); ax_main.set_ylabel("Source side overlap")

    d_src = ol_disyn_mat.mean(axis=1) 
    d_tgt = ol_disyn_mat.mean(axis=0) 
    ax_src.plot(d_src.values, d_src.index)
    ax_tgt.plot(d_tgt.index, d_tgt.values)

    lims = numpy.vstack([ax_src.get_xlim(), ax_tgt.get_ylim()])
    lims = [lims[:, 0].min(), lims[:, 1].max()]
    ax_src.set_xlim(lims); ax_tgt.set_ylim(lims)
    ax_src.set_ylim(ax_main.get_ylim()); ax_tgt.set_xlim(ax_main.get_xlim())
    ax_src.set_yticks([]); ax_tgt.set_xticks([])
    ax_src.set_frame_on(False); ax_tgt.set_frame_on(False)

    return fig

def plot_overlap_vs_disyn_corr(ol_disyn_cc, column_names):
    """
    Fig. 4H
    """
    fig = plt.figure(figsize=(.75, 1.75))
    ax = fig.gca()

    baseline = 0.2
    pairs = [("src", "inh"), ("tgt", "inh"), ("src", "tgt")]
    cols = ["red", "blue", "grey"]
    lbls = []
    for i, pair in enumerate(pairs):
        idx = list(column_names).index(pair[0])
        idy = list(column_names).index(pair[1])
        ax.bar(i, ol_disyn_cc[idx, idy] - baseline, bottom=baseline, color=cols[i],
            label="{0} - {1}".format(*pair))
        lbls.append("{0} - {1}".format(*pair))
    ax.set_xticks(range(len(lbls))); ax.set_xticklabels(lbls, rotation="vertical")
    ax.set_ylabel("Pearson correlation")
    ax.set_frame_on(False)
    return fig

def plot_overlap_vs_disyn_corr_l(ol_disyn_cc, column_names):
    """
    Alternate version.
    """
    fig = plt.figure(figsize=(3.5, 3.5))
    ax = fig.gca()
    plt.colorbar(ax.imshow(ol_disyn_cc, clim=[0, 1.0]))
    ax.set_xticks(range(len(column_names)))
    ax.set_xticklabels(column_names)
    ax.set_yticks(range(len(column_names)))
    ax.set_yticklabels(column_names)
    return fig

def _position_normalize_by_percentile(pos):
    pos = pos - pos.mean(axis=0)
    nrm = numpy.linalg.norm(pos, axis=1)
    percs = numpy.percentile(nrm, numpy.arange(101))
    tgt_nrm = numpy.interp(nrm, percs, numpy.arange(101))
    tgt_pos = tgt_nrm.reshape((-1, 1)) * pos / nrm.reshape((-1, 1))
    return tgt_pos.values

def _position_horizontalize(pos, mean_loc, M, node_is_exc):
    pos = pos - pos.mean(axis=0)
    new_y = numpy.arctan2(pos["x"], pos["y"])
    new_x = numpy.zeros(len(new_y))
    new_x[node_is_exc] = mean_loc.values

    tmp_mat = (M.array[numpy.ix_(node_is_exc, ~node_is_exc)] > 0).astype(int).transpose()

    new_x[~node_is_exc] = numpy.dot(tmp_mat, mean_loc.values.reshape((-1, 1)))[:, 0] / tmp_mat.sum(axis=1)
    new_y = new_y * (new_x + new_x.max())
    return numpy.vstack([new_x, new_y]).transpose()

def _position_normalize_by_mean_pos(pos, mean_loc, M, node_is_exc):
    pos = _position_horizontalize(pos, mean_loc, M, node_is_exc)
    x = pos[:, 0]; y = pos[:, 1]
    circ_pos = numpy.vstack([numpy.cos(y) * x,
                         numpy.sin(y) * x])
    return circ_pos.transpose()

def plot_simplex_group_as_network(smplx_tgt_ids, M, Msmpl, Mnrn, simplices, s_n_paths, n_s_paths, name_nrn,
                                  tgt_dim, cfg):
    """
    Fig. 6A
    """
    import networkx
    import pandas
    import connalysis

    nrn_min_indeg = cfg["nrn_minimum_indegree"]
    nrn_min_outdeg = cfg["nrn_minimum_outdegree"]
    nrn_pick_n = cfg["nrn_num_picked"]
    position_strategies = cfg.get("positioning", [None])

    # All gids associated with the simplices
    exc_gids = Msmpl.gids[numpy.unique(simplices[smplx_tgt_ids.values])]

    # Find exemplary connected neurons from population nrn
    i_degs = (s_n_paths.loc[smplx_tgt_ids].groupby(name_nrn).count(),
              n_s_paths.loc[smplx_tgt_ids].groupby(name_nrn).count())
    is_strongly_connected = (i_degs[0] > nrn_min_indeg) & (i_degs[1] > nrn_min_outdeg)
    inh_gids = Mnrn.gids[numpy.random.choice(is_strongly_connected.index[is_strongly_connected],
                                             nrn_pick_n, replace=False)]

    # Relevant subnetworks
    Msub = M.subpopulation(exc_gids)
    Mfl = M.subpopulation(numpy.hstack([inh_gids, exc_gids]))
    node_is_exc = numpy.in1d(Mfl.gids, exc_gids)

    # For the simplex subnetwork: For each node its mean simplex position. Determines color
    loc_df = connalysis.network.topology.list_simplices_by_dimension(Msub.matrix)[tgt_dim]
    loc_df = pandas.DataFrame(loc_df, columns=numpy.arange(tgt_dim + 1))
    mean_loc = pandas.concat([loc_df[_col] for _col in loc_df.columns], axis=0,
                             keys=loc_df.columns).reset_index().groupby(0)["level_0"].mean()
    cols = plt.cm.coolwarm(mean_loc.values / tgt_dim)

    # All colors: The non-simplex neurons in green
    cols_full = numpy.zeros((len(Mfl), cols.shape[1]))
    cols_full[:, 1] = 1.0
    cols_full[:, 3] = 1.0
    cols_full[node_is_exc, :] = cols

    # Create DiGraph
    grph = networkx.DiGraph(Mfl.array)

    # simplex edges in black, non-simplex edges in green
    e = numpy.vstack(list(grph.edges))
    edge_is_exc = node_is_exc[e[:, 0]]
    edge_cols = numpy.zeros((len(edge_is_exc), 3))
    edge_cols[~edge_is_exc, 1] = 1

    pos = pandas.DataFrame(networkx.spring_layout(grph, iterations=50), index=["x", "y"]).transpose()
    figs = []
    for positioning in position_strategies:
        if positioning is not None:
            if positioning == "percentile":
                use_pos = _position_normalize_by_percentile(pos)
            elif positioning == "horizontal":
                use_pos = _position_horizontalize(pos, mean_loc, Mfl, node_is_exc)
            elif positioning == "mean_pos":
                use_pos = _position_normalize_by_mean_pos(pos, mean_loc, Mfl, node_is_exc)
            elif positioning == "raw":
                use_pos = pos.values
            else:
                raise ValueError(positioning)
        else:
            use_pos = pos.values

        figs.append(plt.figure(figsize=(7, 7)))
        networkx.draw_networkx(grph, node_size=50, node_color=cols_full, pos=use_pos,
                            width=0.25, with_labels=False, edge_color=edge_cols)
    return figs

def save_plot(figs, root, name):
    """
    Saves a figure to specified locations.
    """
    import os
    if not os.path.exists(root):
        os.makedirs(root)
    if not isinstance(figs, list):
        figs = [figs]
    for i, fig in enumerate(figs):
        fn = os.path.join(root, name + str(i) + ".pdf")
        fig.savefig(fn)
