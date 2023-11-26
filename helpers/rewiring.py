import connalysis
import pandas
import numpy

from .analysis import edge_participation_df

def find_low_participation(ep_df, dim, n):
    max_dim = (ep_df > 0).sum(axis=1) - 1
    cands = numpy.nonzero(max_dim.isin(dim))[0]
    n = numpy.minimum(n, len(cands))
    return numpy.random.choice(cands, n, replace=False)

def find_high_participation(M, smpl, dims, n, positions, nmax=50000):
    count_df = pandas.DataFrame(numpy.empty((0, numpy.max(dims) + 1)))
    for dim in dims:
        tmp = pandas.DataFrame(smpl[dim])
        if len(tmp) > nmax:
            _smpl = numpy.random.choice(len(tmp), nmax, replace=False)
            tmp = tmp.iloc[_smpl].reset_index(drop=True)
        tmp = tmp.apply(lambda _col: _col.value_counts(), axis=0).fillna(0)
        count_df = count_df.add(tmp, fill_value=0)
    
    count_df = count_df / count_df.sum(axis=0)
    
    combos = numpy.vstack([(positions[i], positions[j])
                           for i in range(len(positions))
                           for j in range(i + 1, len(positions))])
    n_per_combo = int(n / len(combos))
    picked = []
    
    for combo in combos:
        a = count_df[combo[0]].dropna()
        b = count_df[combo[1]].dropna()
        i = numpy.random.choice(a.index.values, n_per_combo, p=a.values)
        j = numpy.random.choice(b.index.values, n_per_combo, p=b.values)
        v = (i != j)
        i = i[v]; j = j[v]
        
        cands = numpy.vstack([i, j]).transpose()
        cands = cands[~numpy.array(M.matrix.tocsc()[cands[:, 0], cands[:, 1]]).flatten()]
        picked.append(cands)
    
    return numpy.vstack(picked)

def add_and_remove(M, to_add, to_remove):
    print("Removing {0} edges...".format(len(to_remove)))
    idx_keep = numpy.setdiff1d(range(len(M.edges)), to_remove)
    M._edges = M._edges.iloc[idx_keep]
    M._edge_indices = M._edge_indices.iloc[idx_keep]

    assert(len(M.edge_properties) == 1)
    _e_col = M.edge_properties[0]
    e_dtype = M.edges["data"].dtype
    print("Adding {0} edges...".format(len(to_add)))
    M._edges = pandas.concat([
        M._edges,
        pandas.DataFrame({_e_col: numpy.ones(len(to_add), dtype=e_dtype)})
    ], axis=0).reset_index(drop=True)
    M._edge_indices = pandas.concat([
        M._edge_indices,
        pandas.DataFrame({"row": to_add[:, 0], "col": to_add[:, 1]})
    ], axis=0).reset_index(drop=True)

def rewire_step(M, detect_at_dims, remove_at_dim, n, positions, **kwargs):
    smpl = connalysis.network.topology.list_simplices_by_dimension(M.matrix)
    edge_part = edge_participation_df(M)
    print("""
        Before: 
{0}
    """.format(str(((edge_part > 0).sum(axis=1)).value_counts().sort_index())))

    edges_to_remove = find_low_participation(edge_part, remove_at_dim, n)
    edges_to_add = find_high_participation(M, smpl, detect_at_dims, len(edges_to_remove), positions, **kwargs)
    edges_to_remove = edges_to_remove[:len(edges_to_add)]
    add_and_remove(M, edges_to_add, edges_to_remove)

