import igraph as ig
import louvain
import numpy as np
from random import randint
import sys


def isclose(a, b):
    return abs(a - b) <= 1e-10


def calculate_coefficient(com_vec, adj_matrix):
    '''
    For a given connection matrix and set of community labels, calculate the coefficient
    for plane/line associated with that connectivity matrix
    :param com_vec: list or vector with community membership for each element of network
    ordered the same as the rows/col of adj_matrix
    :param adj_matrix: adjacency matrix for connections to calculate coefficients for
    (i.e. A_ij, P_ij, C_ij, etc..) ordered the same as com_vec
    :return:
    '''

    assert com_vec.shape[0] == adj_matrix.shape[0]
    idx_sort = np.argsort(com_vec)
    _, idx_start = np.unique(com_vec[idx_sort], return_index=True)
    c_inds = np.split(idx_sort, idx_start[1:])

    # adjacency matrix sum over the community blocks
    return sum(np.sum(adj_matrix[np.ix_(c_ind, c_ind)]) for c_ind in c_inds)


def test_diff_move():
    intraslice = ig.Graph.Read_Ncol("multilayer_SBM_interslice_edges.csv", directed=False)
    n = intraslice.vcount()
    layer_vec = [0] * n
    membership = list(range(n))

    part_rbc = louvain.RBConfigurationVertexPartition(intraslice, resolution_parameter=1.0,
                                                      initial_membership=membership)
    part_weighted_layers = louvain.RBConfigurationVertexPartitionWeightedLayers(intraslice, resolution_parameter=1.0,
                                                                               layer_vec=layer_vec,
                                                                               initial_membership=membership)

    # check diff_move() - quality() consistency across 100 random moves
    for repeat in range(100):
        v = randint(0, n - 1)
        c = randint(0, n - 1)
        old_quality = part_weighted_layers.quality()
        wl_diff = part_weighted_layers.diff_move(v, c)
        part_weighted_layers.move_node(v, c)
        true_diff = part_weighted_layers.quality() - old_quality

        rbc_diff = part_rbc.diff_move(v, c)
        part_rbc.move_node(v, c)

        assert isclose(wl_diff, true_diff), "WeightedLayers diff_move() inconsistent with quality()"
        assert isclose(wl_diff, rbc_diff), "WeightedLayers diff_move() inconsistent with single-layer"
        assert isclose(part_weighted_layers.quality(),
                       part_rbc.quality()), "WeightedLayers quality() inconsistent with single-layer"

    # check rng consistency between RBConfigurationVertexPartition and its WeightedLayers variant
    louvain.set_rng_seed(0)
    part_weighted_layers = louvain.RBConfigurationVertexPartitionWeightedLayers(intraslice, resolution_parameter=1.0,
                                                                               layer_vec=layer_vec)
    opt = louvain.Optimiser()
    opt.optimise_partition(partition=part_weighted_layers)

    louvain.set_rng_seed(0)
    part_rbc = louvain.RBConfigurationVertexPartition(intraslice, resolution_parameter=1.0)
    opt = louvain.Optimiser()
    opt.optimise_partition(partition=part_rbc)

    assert isclose(part_weighted_layers.quality(),
                   part_rbc.quality()), "WeightedLayers optimisation inconsistent with single-layer"


def test_multilayer_louvain():
    intraslice = ig.Graph.Read_Ncol("multilayer_SBM_intraslice_edges.csv", directed=False)
    interslice = ig.Graph.Read_Ncol("multilayer_SBM_interslice_edges.csv", directed=False)
    n_layers = 4
    n = intraslice.vcount() // n_layers
    layer_vec = np.array([i // n for i in range(n * n_layers)])

    intralayer_part = louvain.RBConfigurationVertexPartitionWeightedLayers(intraslice, resolution_parameter=1.0,
                                                                           layer_vec=layer_vec.tolist())
    interlayer_part = louvain.RBConfigurationVertexPartition(interslice, resolution_parameter=0.0)

    opt = louvain.Optimiser()
    opt.optimise_partition_multiplex(partitions=[intralayer_part, interlayer_part])

    louvain_mod = intralayer_part.quality(resolution_parameter=1.0) + interlayer_part.quality()
    A = np.array(intraslice.get_adjacency()._get_data())
    C = np.array(interslice.get_adjacency()._get_data())
    P = np.zeros((n_layers * n, n_layers * n))
    for i in range(n_layers):
        c_degrees = np.array(intraslice.degree(list(range(n * i, n * i + n))))
        c_inds = np.where(layer_vec == i)[0]
        P[np.ix_(c_inds, c_inds)] = np.outer(c_degrees, c_degrees.T) / (1.0 * np.sum(c_degrees))

    membership = np.array(intralayer_part.membership)
    true_mod = sum(calculate_coefficient(membership, X) for X in (A, -P, C))

    assert isclose(louvain_mod, true_mod), "WeightedLayers quality() inconsistent with alternate calculation"


def main():
    test_multilayer_louvain()
    test_diff_move()


if __name__ == '__main__':
    sys.exit(main())
