import conntility
import pandas
import numpy
import bluepysnap as snap
import json

def load_network(cfg):
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

def filter_network(M, cfg):
    for fltr in cfg:
        idx = M.index(fltr["column"])
        M = getattr(idx, fltr["function"])(*fltr.get("args", []), **fltr.get("kwargs", {}))
    return M

