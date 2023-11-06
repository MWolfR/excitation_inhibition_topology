import conntility
import pandas
import numpy
import json

def load_network(cfg):
    if "conntility" in cfg.keys():
        fn = cfg["conntility"]
        return conntility.ConnectivityMatrix.from_h5(fn, *cfg.get("args", []), **cfg.get("kwargs", {}))
    elif "snap" in cfg.keys():
        raise NotImplementedError()

def filter_network(M, cfg):
    for fltr in cfg:
        idx = M.index(fltr["column"])
        M = getattr(idx, fltr["function"])(*fltr.get("args", []), **fltr.get("kwargs", {}))
    return M

