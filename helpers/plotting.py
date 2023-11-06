import numpy

from matplotlib import pyplot as plt

def plot_simplex_divergences(divergences):
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

def plot_inout_degrees(m_nrn_smpl, m_smpl_nrn):
    outdeg = numpy.array(m_nrn_smpl.sum(axis=1)).flatten()
    indeg = numpy.array(m_smpl_nrn.sum(axis=0)).flatten()

    fig = plt.figure(figsize=(3.5, 3))
    ax = fig.gca()
    ax.plot(outdeg, indeg, '.', ms=2)
    # plt.axis("equal")
    plt.gca().set_frame_on(False)

    plt.gca().set_xlabel("Outdegree (to simplices)")
    plt.gca().set_ylabel("Indegree (from simplices)")
    print(numpy.corrcoef(indeg, outdeg))
    return fig, indeg, outdeg

def plot_disynaptic_path_sum(mat_disyn):
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
    fig = plt.figure(figsize=(3.5, 3.5))
    ax = fig.gca()
    ax.imshow(paths_df.groupby(["Simplex", "Position"]).count().unstack("Position", fill_value=0),
              aspect="auto", interpolation="nearest")
    ax.set_ylabel("Simplex")
    ax.set_xlabel("Position")
    return fig
    
def plot_position_mean(paths_df):
    fig = plt.figure(figsize=(4.5, 3.5))
    ax = fig.gca()
    plt.colorbar(ax.imshow(paths_df.groupby(["Simplex", "Position"]).mean().unstack("Position", fill_value=0),
                           aspect="auto", interpolation="nearest"))
    ax.set_ylabel("Simplex")
    ax.set_xlabel("Position")
    return fig

def compare_disyn_inhibition(m_data, m_ctrl):
    bins = numpy.linspace(0, m_data.max(), 101)
    H_data = numpy.histogram(m_data.flatten(), bins=bins)[0]
    H_ctrl = numpy.histogram(m_ctrl.flatten(), bins=bins)[0]

    plt.plot(bins[:-1], H_data, label="Data")
    plt.plot(bins[:-1], H_ctrl, label="Shuffled")

    plt.gca().set_xlabel("Disyn. connection strength")
    plt.gca().set_ylabel("Simplex pairs")
    plt.legend()
    return plt.gcf()
    
def plot_overlaps_vs_disyn(ol_disyn_mat):
    fig = plt.figure(figsize=(6, 6))
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
    fig = plt.figure(figsize=(3.5, 3.5))
    ax = fig.gca()
    plt.colorbar(ax.imshow(ol_disyn_cc, clim=[0, 1.0]))
    ax.set_xticks(range(len(column_names)))
    ax.set_xticklabels(column_names)
    ax.set_yticks(range(len(column_names)))
    ax.set_yticklabels(column_names)
    return fig

