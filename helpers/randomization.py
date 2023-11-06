import pandas
import numpy
import connalysis
import conntility

from .network import filter_network

def create_dd_control_networks(M, Msmpl, Mnrn, cols_coords, model_name, cfg):
    is_smpl = numpy.in1d(M.gids, Msmpl.gids)
    is_nrn = numpy.in1d(M.gids, Mnrn.gids)

    mdl_smpl_smpl = connalysis.modelling.conn_prob_2nd_order_model(Msmpl.matrix.tocsc(),
                                               Msmpl.vertices,
                                               bin_size_um=25000, coord_names=cols_coords)
    mdl_nrn_nrn = connalysis.modelling.conn_prob_2nd_order_model(Mnrn.matrix.tocsc(),
                                                Mnrn.vertices,
                                                bin_size_um=25000, coord_names=cols_coords)
    m = M.matrix.tocsc()
    mdl_nrn_smpl = connalysis.modelling.conn_prob_2nd_order_pathway_model(m[numpy.ix_(is_nrn, is_smpl)],
                                                                        Mnrn.vertices,
                                                                        Msmpl.vertices,
                                                                        bin_size_um=25000,
                                                                        coord_names=cols_coords)
    mdl_smpl_nrn = connalysis.modelling.conn_prob_2nd_order_pathway_model(m[numpy.ix_(is_smpl, is_nrn)],
                                                                        Msmpl.vertices,
                                                                        Mnrn.vertices,
                                                                        bin_size_um=25000,
                                                                        coord_names=cols_coords)
    
    n = numpy.sum(is_smpl) + numpy.sum(is_nrn)
    params = numpy.array([[mdl_smpl_smpl.values[0], mdl_smpl_nrn.values[0]],
            [mdl_nrn_smpl.values[0], mdl_nrn_nrn.values[0]]])
    blocks = numpy.hstack([numpy.zeros(numpy.sum(is_smpl)), numpy.ones(numpy.sum(is_smpl))]).astype(int)
    xyz = pandas.concat([M.vertices[cols_coords].loc[is_smpl], M.vertices[cols_coords].loc[is_nrn]], axis=0).values.astype(float)

    edges = connalysis.randomization.run_DD2_block(
        n,
        params,
        blocks,
        xyz,
        8,
    )
    edges = pandas.DataFrame(edges)
    _vp = pandas.concat([M._vertex_properties.loc[is_smpl],
                        M._vertex_properties.loc[is_nrn]], axis=0).reset_index(drop=True)
    _ep = pandas.DataFrame({"value": numpy.ones(len(edges))})
    Mctrl = conntility.ConnectivityMatrix(edges, vertex_properties=_vp, edge_properties=_ep)

    Msmpl = filter_network(Mctrl, cfg["connectome"]["subnetwork"]["simplices"])
    Mnrn = filter_network(Mctrl, cfg["connectome"]["subnetwork"]["neurons"]["filters"])
    return Msmpl, Mnrn
