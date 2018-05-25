#include "RBConfigurationVertexPartitionWeightedLayers.h"

RBConfigurationVertexPartitionWeightedLayers::RBConfigurationVertexPartitionWeightedLayers(Graph* graph,
      vector<size_t> const& _layer_vec, vector<vector<double> > const& _degree_by_layers,
      vector<size_t> const& membership, double resolution_parameter) :
        LinearResolutionParameterVertexPartition(graph,
        membership, resolution_parameter), _layer_vec(_layer_vec),_degree_by_layers(_degree_by_layers)
{
  this->setLayerVec(_layer_vec);
  this->setDegreeByLayers(_degree_by_layers);
}

RBConfigurationVertexPartitionWeightedLayers::RBConfigurationVertexPartitionWeightedLayers(Graph* graph,
      vector<size_t> const& _layer_vec, vector<vector<double> > const& degree_by_layers,
      vector<size_t> const& membership) :
        LinearResolutionParameterVertexPartition(graph,
        membership), _layer_vec(_layer_vec),_degree_by_layers(_degree_by_layers)
{
  this->setLayerVec(_layer_vec);
  this->setDegreeByLayers(_degree_by_layers);
}

RBConfigurationVertexPartitionWeightedLayers::RBConfigurationVertexPartitionWeightedLayers(Graph* graph,
      vector<size_t> const& _layer_vec, vector<vector<double> > const& _degree_by_layers,
      double resolution_parameter) :
        LinearResolutionParameterVertexPartition(graph, resolution_parameter),_layer_vec(_layer_vec),_degree_by_layers(_degree_by_layers)
{
  this->setLayerVec(_layer_vec);
  this->setDegreeByLayers(_degree_by_layers);
}

RBConfigurationVertexPartitionWeightedLayers::RBConfigurationVertexPartitionWeightedLayers(Graph* graph,
    vector<size_t> const& _layer_vec, vector<vector<double> > const& _degree_by_layers) :
        LinearResolutionParameterVertexPartition(graph), _layer_vec(_layer_vec),_degree_by_layers(_degree_by_layers)
{
  this->setLayerVec(_layer_vec);
  this->setDegreeByLayers(_degree_by_layers);
}

RBConfigurationVertexPartitionWeightedLayers::~RBConfigurationVertexPartitionWeightedLayers()
{ }

//CREATE METHODS
RBConfigurationVertexPartitionWeightedLayers* RBConfigurationVertexPartitionWeightedLayers::create(Graph* graph)
{
  return new RBConfigurationVertexPartitionWeightedLayers(graph, this->getLayerVec(), this->getDegreeByLayers(),
                                                          this->resolution_parameter);
}

RBConfigurationVertexPartitionWeightedLayers* RBConfigurationVertexPartitionWeightedLayers::create(Graph* graph, vector<size_t> const& membership)
{
  vector<vector<double> > new_degree_by_layer = this->_condense_degree_by_layer(membership);
  return new RBConfigurationVertexPartitionWeightedLayers(graph, this->getLayerVec(), new_degree_by_layer, membership, this->resolution_parameter);
}

//flatten the degree by layer matrix by the common
vector<vector<double> > RBConfigurationVertexPartitionWeightedLayers::_condense_degree_by_layer(vector<size_t> const& membership) {

  map<int, int> group_sizes;
  vector<vector<size_t> > ind_to_add;

  size_t cgroup;
  for (size_t v = 0; v < membership.size(); v++) {
    cgroup = membership[v];
    ind_to_add[cgroup].push_back(v);
  }
  size_t N = this->getDegreeByLayers().size();    //number of nodes
  size_t M = this->getDegreeByLayers()[0].size(); //number of layers
  vector<vector<double> > condensed_degrees(N, vector<double>(M, 0.0));
  // condensed_degrees.resize(ind_to_add.size())
  vector<size_t> cinds;
  for (size_t cgroup = 0; cgroup < ind_to_add.size(); cgroup++)
  {
    cinds = ind_to_add[cgroup]; //current rows of the degree by layers to add
    for (size_t i = 0; i < cinds.size(); i++)
    {
      for (size_t j = 0; j < M; j++)
      {
        condensed_degrees[cgroup][j] += this->getDegreeByLayers()[cinds[i]][j];
      }
    }
  }
  return condensed_degrees;
}

/*****************************************************************************
  Returns the difference in modularity if we move a node to a new community
*****************************************************************************/
double RBConfigurationVertexPartitionWeightedLayers::diff_move(size_t v, size_t new_comm)
{
  #ifdef DEBUG
    cerr << "double RBConfigurationVertexPartitionWeightedLayers::diff_move(" << v << ", " << new_comm << ")" << endl;
  #endif
  size_t old_comm = this->_membership[v];
  double diff = 0.0;
  double total_weight = this->graph->total_weight()*(2.0 - this->graph->is_directed());
  if (total_weight == 0.0)
    return 0.0;
  if (new_comm != old_comm)
  {
    #ifdef DEBUG
      cerr << "\t" << "old_comm: " << old_comm << endl;
    #endif
    double w_to_old = this->weight_to_comm(v, old_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_to_old: " << w_to_old << endl;
    #endif
    double w_from_old = this->weight_from_comm(v, old_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_from_old: " << w_from_old << endl;
    #endif
    double w_to_new = this->weight_to_comm(v, new_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_to_new: " << w_to_new << endl;
    #endif
    double w_from_new = this->weight_from_comm(v, new_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_from_new: " << w_from_new << endl;
    #endif
    double k_out = this->graph->strength(v, IGRAPH_OUT);
    #ifdef DEBUG
      cerr << "\t" << "k_out: " << k_out << endl;
    #endif
    double k_in = this->graph->strength(v, IGRAPH_IN);
    #ifdef DEBUG
      cerr << "\t" << "k_in: " << k_in << endl;
    #endif
    double self_weight = this->graph->node_self_weight(v);
    #ifdef DEBUG
      cerr << "\t" << "self_weight: " << self_weight << endl;
    #endif
    double K_out_old = this->total_weight_from_comm(old_comm);
    #ifdef DEBUG
      cerr << "\t" << "K_out_old: " << K_out_old << endl;
    #endif
    double K_in_old = this->total_weight_to_comm(old_comm);
    #ifdef DEBUG
      cerr << "\t" << "K_in_old: " << K_in_old << endl;
    #endif
    double K_out_new = this->total_weight_from_comm(new_comm) + k_out;
    #ifdef DEBUG
      cerr << "\t" << "K_out_new: " << K_out_new << endl;
    #endif
    double K_in_new = this->total_weight_to_comm(new_comm) + k_in;
    #ifdef DEBUG
      cerr << "\t" << "K_in_new: " << K_in_new << endl;
      cerr << "\t" << "total_weight: " << total_weight << endl;
    #endif
    double diff_old = (w_to_old - this->resolution_parameter*k_out*K_in_old/total_weight) + \
               (w_from_old - this->resolution_parameter*k_in*K_out_old/total_weight);
    #ifdef DEBUG
      cerr << "\t" << "diff_old: " << diff_old << endl;
    #endif
    double diff_new = (w_to_new + self_weight - this->resolution_parameter*k_out*K_in_new/total_weight) + \
               (w_from_new + self_weight - this->resolution_parameter*k_in*K_out_new/total_weight);
    #ifdef DEBUG
      cerr << "\t" << "diff_new: " << diff_new << endl;
    #endif
    diff = diff_new - diff_old;
    #ifdef DEBUG
      cerr << "\t" << "diff: " << diff << endl;
    #endif
  }
  #ifdef DEBUG
    cerr << "exit RBConfigurationVertexPartitionWeightedLayers::diff_move(" << v << ", " << new_comm << ")" << endl;
    cerr << "return " << diff << endl << endl;
  #endif
  return diff;
}

/*****************************************************************************
  Give the modularity of the partition.

  We here use the unscaled version of modularity, in other words, we don"t
  normalise by the number of edges.
******************************************************************************/
double RBConfigurationVertexPartitionWeightedLayers::quality(double resolution_parameter)
{
  #ifdef DEBUG
    cerr << "double ModularityVertexPartition::quality()" << endl;
  #endif
  double mod = 0.0;

  double m;
  if (this->graph->is_directed())
    m = this->graph->total_weight();
  else
    m = 2*this->graph->total_weight();

  if (m == 0)
    return 0.0;

  for (size_t c = 0; c < this->nb_communities(); c++)
  {
    double w = this->total_weight_in_comm(c);
    double w_out = this->total_weight_from_comm(c);
    double w_in = this->total_weight_to_comm(c);
    #ifdef DEBUG
      size_t csize = this->csize(c);
      cerr << "\t" << "Comm: " << c << ", size=" << csize << ", w=" << w << ", w_out=" << w_out << ", w_in=" << w_in << "." << endl;
    #endif
    mod += w - resolution_parameter*w_out*w_in/((this->graph->is_directed() ? 1.0 : 4.0)*this->graph->total_weight());
  }
  double q = (2.0 - this->graph->is_directed())*mod;
  #ifdef DEBUG
    cerr << "exit double RBConfigurationVertexPartitionWeightedLayers::quality()" << endl;
    cerr << "return " << q << endl << endl;
  #endif
  return q;
}