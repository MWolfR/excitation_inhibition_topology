import json

from helpers import network


def main():
    import sys
    with open(sys.argv[1], "r") as fid:
        cfg = json.load(fid)

    M = network.load_network(cfg["connectome"]["loading"])
    network.override_connectivity(M, cfg["connectome"]["override"])
    M.to_h5(cfg["connectome"]["save"])


if __name__ == "__main__":
    main()
