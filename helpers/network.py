import conntility
import pandas
import numpy
import bluepysnap as snap
import json

"""
Helper functions to load the connectivity of a network from multiple sources.
To filter them into a subnetwork. And generate random control networks.
"""

def load_network(cfg):
    """
    Load the connectivity of a network from sources specified in cfg.

    Returns it as a conntility.ConnectivityMatrix object.
    """
    if "conntility" in cfg.keys():
        fn = cfg["conntility"]
        return conntility.ConnectivityMatrix.from_h5(fn, *cfg.get("args", []), **cfg.get("kwargs", {}))
    elif "snap" in cfg.keys():
        load_cfg = {
            "loading":{    
                "properties": ["x", "y", "z", "mtype", "layer", "synapse_class"],
                "atlas": [
                    {
                        "data": cfg["snap"]["flatmap"],
                        "properties": ["ss_flat_x", "depth", "ss_flat_y"]
                    }
                ]
            }
        }
        connectome = cfg["snap"]["connectome"]
        fn_circ = cfg["snap"]["circuit"]
        circ = snap.Circuit(fn_circ)
        cmat = conntility.ConnectivityMatrix.from_bluepy(circ, load_config=load_cfg,
                                                        connectome=connectome,
                                                        load_full=True)
        return cmat
    elif "npz" in cfg.keys():
        from scipy import sparse
        mat = sparse.load_npz(cfg["npz"])
        vertices = None
        if "vertices" in cfg.keys():
            vertices = pandas.read_pickle(cfg["vertices"])
        return conntility.ConnectivityMatrix(mat, vertex_properties=vertices)

def filter_network(M, cfg):
    """
    Filter a ConnectivityMatrix according the specs in cfg.

    Returns a filtered ConnectivityMatrix.
    """
    for fltr in cfg:
        idx = M.index(fltr["column"])
        M = getattr(idx, fltr["function"])(*fltr.get("args", []), **fltr.get("kwargs", {}))
    return M

def _pathway_idxx(M, cfg):
    if isinstance(cfg, dict):
        Msub = filter_network(M, cfg["pre"])
        idxx_pre = numpy.nonzero(numpy.in1d(M.gids, Msub.gids))[0]
        Msub = filter_network(M, cfg["post"])
        idxx_post = numpy.nonzero(numpy.in1d(M.gids, Msub.gids))[0]
    else:
        Msub = filter_network(M, cfg)
        idxx_pre = numpy.nonzero(numpy.in1d(M.gids, Msub.gids))[0]
        idxx_post = idxx_pre
    return idxx_pre, idxx_post

def load_override_connectome(cfg):
    """
    Loads a connectivity matrix from a pickle file that can be used to override parts
    of a connectome.

    Returns a sparse matrix.
    """
    import pickle
    from scipy import sparse
    assert "pkl" in cfg.keys()
    with open(cfg["pkl"], 'rb') as f: 
        mat=pickle.load(f) 
    if "key" in cfg.keys():
        lst_keys = cfg["key"]
        if not isinstance(lst_keys, list): lst_keys = [lst_keys]
        for _k in lst_keys:
            mat = mat[_k]
    assert isinstance(mat, sparse.spmatrix)
    return mat

def override_connectivity(M, cfg):
    """
    Overrides parts of a connectome with a custom (rewired) connectome as specified in cfg.

    Returns None, but changes the connectivity of M in place.
    """
    for override in cfg:
        print("Executing override...")
        idxx_pre, idxx_post = _pathway_idxx(M, override["pathway"])
        kwargs={}
        if "node_subselection" in override.keys():
            kwargs['restrict']=override['node_subselection']['restrict_nodes']
            kwargs['per']= override['node_subselection']['top_percentile_to_choose'] 
        if "matrix" in override.keys():
            repl_mat = load_override_connectome(override["matrix"])
        elif "rewire" in override.keys():
            from .rewiring import rewire_step

            repl_mat = filter_network(M, override["pathway"])
            for rw in override["rewire"]:
                print("Executing rewiring...")
                rewire_step(repl_mat, rw["dims_add"],
                            rw["dims_remove"], rw["n"],
                            rw["positions"], 
                           **kwargs)
            repl_mat = repl_mat.matrix
        assert len(idxx_pre) == repl_mat.shape[0]
        assert len(idxx_post) == repl_mat.shape[1]

        is_pw = M._edge_indices["row"].isin(idxx_pre) & M._edge_indices["col"].isin(idxx_post)
        repl_mat = repl_mat.tocoo()
        repl_df = pandas.DataFrame({"row": idxx_pre[repl_mat.row], "col": idxx_post[repl_mat.col]})
        M._edge_indices = pandas.concat([M._edge_indices.loc[~is_pw], repl_df], axis=0).reset_index(drop=True)

        assert len(M.edge_properties) == 1
        M._edges = pandas.concat([M.edges.loc[~is_pw], pandas.DataFrame({M.edge_properties[0]: repl_mat.data})],
                                axis=0).reset_index(drop=True)


