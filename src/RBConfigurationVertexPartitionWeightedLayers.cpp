#include "RBConfigurationVertexPartitionWeightedLayers.h"

RBConfigurationVertexPartitionWeightedLayers::RBConfigurationVertexPartitionWeightedLayers(
    Graph *graph, vector<size_t> const &_layer_vec, vector<vector<double> > const &_degree_by_layers,
    vector<size_t> const &membership, double resolution_parameter)
    : LinearResolutionParameterVertexPartition(graph, membership, resolution_parameter),
      _layer_vec(_layer_vec), _degree_by_layers(_degree_by_layers),
      _total_layer_weights(_compute_total_layer_weights(_degree_by_layers))
{ }

RBConfigurationVertexPartitionWeightedLayers::RBConfigurationVertexPartitionWeightedLayers(
    Graph *graph, vector<size_t> const &_layer_vec, vector<vector<double> > const &_degree_by_layers,
    vector<size_t> const &membership)
    : LinearResolutionParameterVertexPartition(graph, membership), _layer_vec(_layer_vec),
      _degree_by_layers(_degree_by_layers),
      _total_layer_weights(_compute_total_layer_weights(_degree_by_layers))
{ }

RBConfigurationVertexPartitionWeightedLayers::RBConfigurationVertexPartitionWeightedLayers(
    Graph *graph, vector<size_t> const &_layer_vec, vector<vector<double> > const &_degree_by_layers,
    double resolution_parameter)
    : LinearResolutionParameterVertexPartition(graph, resolution_parameter), _layer_vec(_layer_vec),
      _degree_by_layers(_degree_by_layers),
      _total_layer_weights(_compute_total_layer_weights(_degree_by_layers))
{ }

RBConfigurationVertexPartitionWeightedLayers::RBConfigurationVertexPartitionWeightedLayers(
    Graph *graph, vector<size_t> const &_layer_vec, vector<vector<double> > const &_degree_by_layers)
    : LinearResolutionParameterVertexPartition(graph), _layer_vec(_layer_vec),
      _degree_by_layers(_degree_by_layers),
      _total_layer_weights(_compute_total_layer_weights(_degree_by_layers))
{ }

RBConfigurationVertexPartitionWeightedLayers::~RBConfigurationVertexPartitionWeightedLayers()
{ }

//CREATE METHODS
RBConfigurationVertexPartitionWeightedLayers* RBConfigurationVertexPartitionWeightedLayers::create(Graph* graph)
{
  return new RBConfigurationVertexPartitionWeightedLayers(graph, this->_layer_vec, this->_degree_by_layers, this->resolution_parameter);
}

RBConfigurationVertexPartitionWeightedLayers* RBConfigurationVertexPartitionWeightedLayers::create(Graph* graph, vector<size_t> const& membership)
{
  vector<vector<double> > new_degree_by_layer = this->_condense_degree_by_layer(membership);
  return new RBConfigurationVertexPartitionWeightedLayers(graph, this->_layer_vec, new_degree_by_layer, membership, this->resolution_parameter);
}

vector<double> RBConfigurationVertexPartitionWeightedLayers::_compute_total_layer_weights(vector<vector<double> > const& degree_by_layers)
{
  vector<double> total_layer_weights;
  total_layer_weights.resize(degree_by_layers.size());

  for (size_t l = 0; l < degree_by_layers.size(); ++l)
    for (size_t v = 0; v < degree_by_layers[l].size(); ++v)
      total_layer_weights[l] += degree_by_layers[l][v];

  return total_layer_weights;
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
  size_t N = this->_degree_by_layers.size();    //number of nodes
  size_t M = this->_degree_by_layers[0].size(); //number of layers
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
        condensed_degrees[cgroup][j] += this->_degree_by_layers[cinds[i]][j];
      }
    }
  }
  return condensed_degrees;
}


void RBConfigurationVertexPartitionWeightedLayers::_clear_resize(vector<vector<double > > &input_vec, size_t N, size_t M) {

    input_vec.clear();
    input_vec.resize(N);
    for (size_t i =0; i<N; i++){
        input_vec[i].resize(M);
    }

}
//helper method for totalling weight vectors
vector <double> RBConfigurationVertexPartitionWeightedLayers::add_vectors(vector<double> &v1, vector<double> &v2){

    if (v1.size() != v2.size() ) {
            PyErr_SetString(PyExc_ValueError, "size of vectors to add must be equal.");
            return NULL;
        }
    vector<double> outvec(v1.size(),0);
    for (size_t i=0;i<outvec.size;i++){

        outvec[i]=v1[i]+v2[i];
    }

    return outvec;
}

/*****************************************************************************
  Overriden methods from Mutable Vertex Partition
*****************************************************************************/

void RBConfigurationVertexPartitionWeightedLayers::init_admin()
{
  #ifdef DEBUG
    cerr << "void RBConfigurationVertexPartitionWeightedLayers::init_admin()" << endl;
  #endif
  size_t n = this->graph->vcount();

  // First determine number of communities (assuming they are consecutively numbered
  size_t nb_comms = 0;
  for (size_t i = 0; i < n; i++)
  {
    if (this->_membership[i] + 1 > nb_comms)
      nb_comms = this->_membership[i] + 1;
  }

  size_t nb_layers=0;
  for (size_t i = 0; i < this->_layer_vec ; i++)
  {
    if (this->_layer_vec[i] + 1 > nb_layers)
      nb_layers = this->_layer_vec[i] + 1;
  }

  // Reset administration
  this->community.clear();
  for (size_t i = 0; i < nb_comms; i++)
    this->community.push_back(new set<size_t>());

  this->_clear_resize(_total_weight_in_comm_by_layer,nb_comms,nb_layers)
  this->_clear_resize(_total_weight_to_comm_by_layer,nb_comms,nb_layers)
  this->_clear_resize(_total_weight_from_comm_by_layer,nb_comms,nb_layers)

  this->_csize.clear();
  this->_csize.resize(nb_comms);

  this->_current_node_cache_community_from = n + 1; this->_clear_resize(_cached_weight_from_community,n,nb_layers);
  this->_current_node_cache_community_to = n + 1;   this->_clear_resize(_cached_weight_to_community,n,nb_layers);
  this->_current_node_cache_community_all = n + 1;  this->_clear_resize(_cached_weight_all_community,n,nb_layers);

  this->_total_weight_in_all_comms = 0.0;
  for (size_t v = 0; v < n; v++)
  {
    size_t v_comm = this->_membership[v];
    // Add this node to the community sets
    this->community[v_comm]->insert(v);
    // Update the community size
    this->_csize[v_comm] += this->graph->node_size(v);
  }

  size_t m = graph->ecount();
  for (size_t e = 0; e < m; e++)
  {
    pair<size_t, size_t> endpoints = this->graph->get_endpoints(e);
    size_t v = endpoints.first;
    size_t u = endpoints.second;

    size_t v_comm = this->_membership[v];
    size_t u_comm = this->_membership[u];

    // Get the weights of the edge
//    double w = this->graph->edge_weight(e);
    vector< double>  w_layers = this->graph->edge_weight_layers(e);


    // Add weight to the outgoing weight of community of v
    this->_total_weight_from_comm[v_comm]=this->add_vectors(_total_weight_from_comm[v_comm],w_layers);
    #ifdef DEBUG
      cerr << "\t" << "Add (" << v << ", " << u << ") weight " << w << " to from_comm " << v_comm <<  "." << endl;
    #endif
    // Add weight to the incoming weight of community of u
    this->_total_weight_to_comm[u_comm] = this->add_vectors(_total_weight_to_comm[u_comm],w_layers);
    #ifdef DEBUG
      cerr << "\t" << "Add (" << v << ", " << u << ") weight " << w << " to to_comm " << u_comm << "." << endl;
    #endif
    if (!this->graph->is_directed())
    {
      #ifdef DEBUG
        cerr << "\t" << "Add (" << u << ", " << v << ") weight " << w << " to from_comm " << u_comm <<  "." << endl;
      #endif
      this->_total_weight_from_comm[u_comm] = this->add_vectors(_total_weight_from_comm[u_comm],w_layers);
      #ifdef DEBUG
        cerr << "\t" << "Add (" << u << ", " << v << ") weight " << w << " to to_comm " << v_comm << "." << endl;
      #endif
      this->_total_weight_to_comm[v_comm] = this->add_vectors(_total_weight_to_comm[v_comm],w_layers);
    }
    // If it is an edge within a community
    if (v_comm == u_comm)
    {
      this->_total_weight_in_comm[v_comm] = this->add_vectors(_total_weight_in_comm[v_comm],w_layers);
      this->_total_weight_in_all_comms = this->add_vectors(_total_weight_in_all_comms,w_layers);
      #ifdef DEBUG
        cerr << "\t" << "Add (" << v << ", " << u << ") weight " << w << " to in_comm " << v_comm << "." << endl;
      #endif
    }
  }

  this->_total_possible_edges_in_all_comms = 0;
  for (size_t c = 0; c < nb_comms; c++)
  {
    size_t n_c = this->csize(c);
    size_t possible_edges = this->graph->possible_edges(n_c);

    #ifdef DEBUG
      cerr << "\t" << "c=" << c << ", n_c=" << n_c << ", possible_edges=" << possible_edges << endl;
    #endif

    this->_total_possible_edges_in_all_comms += possible_edges;

    // It is possible that some community have a zero size (if the order
    // is for example not consecutive. We add those communities to the empty
    // communities vector for consistency.
    if (this->community[c]->size() == 0)
      this->_empty_communities.push_back(c);
  }

  #ifdef DEBUG
    cerr << "exit MutableVertexPartition::init_admin()" << endl << endl;
  #endif

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