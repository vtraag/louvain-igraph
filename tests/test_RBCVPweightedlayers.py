from context import louvain
import champ
import igraph as ig
import numpy as np
import scipy.io as scio
import pandas as pd
import sys, os
import numpy as np
np.set_printoptions(precision=4, suppress=True)
import matplotlib.pyplot as plt
from time import time
import seaborn as sbn
import sklearn.metrics as skm
import champ
import modbp
import scipy.io as scio


def test_diff_move():
	n = 100
	q = 4
	nlayers = 4
	nblocks = q
	c = 4
	ep = .1
	eta = .1
	pin = (n * c / (2.0)) / ((n / float(q)) * (n / float(q) - 1) / 2.0 * float(q) + ep * (q * (q - 1) / 2.0) * (
		np.power(n / (q * 1.0), 2.0)))

	pout = ep * pin

	prob_mat = np.identity(nblocks) * pin + (np.ones((nblocks, nblocks)) - np.identity(nblocks)) * pout

	sbm = modbp.RandomSBMGraph(n, comm_prob_mat=prob_mat, use_gcc=False)
	layer_vec = [0 for _ in range(sbm.graph.vcount())]  # add extra node to use
	memvec=sbm.block
	intraslice, interslice = champ.create_multilayer_igraph_from_edgelist(intralayer_edges=sbm.get_edgelist(),
																		  interlayer_edges=[], layer_vec=layer_vec)

	# partobj = louvain.RBConfigurationVertexPartitionWeightedLayers(intraslice, resolution_parameter=1.0,
	# 															   layer_vec=layer_vec,initial_membership=memvec)
	# partobj2 = louvain.RBConfigurationVertexPartition(intraslice, resolution_parameter=1.0,
	# 															   initial_membership=memvec)
	# # partobj.diff_move(0,1)
	# # partobj2.diff_move(0,1)
	# #
	louvain.set_rng_seed(0)
	partobj = louvain.RBConfigurationVertexPartitionWeightedLayers(intraslice, resolution_parameter=1.0,
																   layer_vec=layer_vec)
	opt = louvain.Optimiser()
	opt.optimise_partition(partition=partobj)
	print(partobj.quality(resolution_parameter=1.0))

	louvain.set_rng_seed(0)
	partobj2 = louvain.RBConfigurationVertexPartition(intraslice, resolution_parameter=1.0)
	opt = louvain.Optimiser()
	opt.optimise_partition(partition=partobj2)
	print(partobj2.quality(resolution_parameter=1.0))


def test_multilayer_louvain():


	n = 100
	q = 4
	nlayers = 4
	nblocks = q
	c = 4
	ep = .1
	eta = .1
	pin = (n * c / (2.0)) / ((n / float(q)) * (n / float(q) - 1) / 2.0 * float(q) + ep * (q * (q - 1) / 2.0) * (
		np.power(n / (q * 1.0), 2.0)))

	pout = ep * pin


	prob_mat = np.identity(nblocks) * pin + (np.ones((nblocks, nblocks)) - np.identity(nblocks)) * pout


	ml_sbm = modbp.MultilayerSBM(n, comm_prob_mat=prob_mat, layers=nlayers, transition_prob=eta)



	mgraph = modbp.MultilayerGraph(intralayer_edges=ml_sbm.intraedges, interlayer_edges=ml_sbm.interedges,
								   layer_vec=ml_sbm.layer_vec,
								   comm_vec=ml_sbm.get_all_layers_block())
	intraslice, interslice = champ.create_multilayer_igraph_from_edgelist(intralayer_edges=mgraph.intralayer_edges,										   interlayer_edges=mgraph.interlayer_edges,layer_vec=mgraph.layer_vec)




	alldegrees=intraslice.degree()
	layers=np.unique(mgraph.layer_vec)
	degree_by_layer=np.zeros((len(alldegrees),len(layers)))

	for i,d in enumerate(alldegrees): #split out by layers
		degree_by_layer[i][mgraph.layer_vec[i]]=d

	print(np.sum(degree_by_layer, axis=0))



	RBCpartobj=louvain.RBConfigurationVertexPartitionWeightedLayers(intraslice, resolution_parameter=1.0,
																	layer_vec=mgraph.layer_vec.tolist())
	InterlayerPartobj=louvain.RBConfigurationVertexPartition(interslice,resolution_parameter=0.0)

	opt=louvain.Optimiser()
	opt.optimise_partition_multiplex(partitions=[RBCpartobj,InterlayerPartobj])


	print(RBCpartobj.quality(resolution_parameter=1.0))
	print()
	# gamma=1.0
	# omega=1
	# args = (layers, interslice, gamma, omega)
	# part = champ.louvain_ext._parallel_run_louvain_multimodularity(args)

	# print(skm.adjusted_mutual_info_score(part[0]['partition'],mgraph.comm_vec))

	# plt.close()
	# a=plt.subplot2grid((1,2),(0,0))
	# a.set_title('Ground')
	# champ.plot_multilayer_community(mgraph.comm_vec,mgraph.layer_vec,ax=a)
	# a=plt.subplot2grid((1,2),(0,1))
	# a.set_title("Multilayer gamma={:.3f},omega={:.3f}".format(gamma,omega))
	# champ.plot_multilayer_community(part[0]['partition'],mgraph.layer_vec,ax=a)
	# plt.gcf().set_size_inches((10,5))
	# plt.show()

def main():
	test_multilayer_louvain()

if __name__=='__main__':
	sys.exit(main())