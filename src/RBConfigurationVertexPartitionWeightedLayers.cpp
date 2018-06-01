#include "RBConfigurationVertexPartitionWeightedLayers.h"

RBConfigurationVertexPartitionWeightedLayers::RBConfigurationVertexPartitionWeightedLayers(Graph *graph,
    vector<size_t> const &membership, double resolution_parameter)
    : LinearResolutionParameterVertexPartition(graph, membership, resolution_parameter),
      _layer_vec(graph->layer_memberships()), _degree_by_layers(graph->degrees_by_layers()),
      _total_layer_weights(_compute_total_layer_weights(_degree_by_layers))
{
    init_admin();
}

RBConfigurationVertexPartitionWeightedLayers::RBConfigurationVertexPartitionWeightedLayers(Graph *graph,
    vector<size_t> const &membership)
    : LinearResolutionParameterVertexPartition(graph, membership), _layer_vec(graph->layer_memberships()),
      _degree_by_layers(graph->degrees_by_layers()),
      _total_layer_weights(_compute_total_layer_weights(_degree_by_layers))
{
   init_admin();
}

RBConfigurationVertexPartitionWeightedLayers::RBConfigurationVertexPartitionWeightedLayers(Graph *graph,
    double resolution_parameter)
    : LinearResolutionParameterVertexPartition(graph, resolution_parameter), _layer_vec(graph->layer_memberships()),
      _degree_by_layers(graph->degrees_by_layers()),
      _total_layer_weights(_compute_total_layer_weights(_degree_by_layers))
{
 init_admin();
}


RBConfigurationVertexPartitionWeightedLayers::RBConfigurationVertexPartitionWeightedLayers(Graph *graph)
    : LinearResolutionParameterVertexPartition(graph), _layer_vec(graph->layer_memberships()),
      _degree_by_layers(graph->degrees_by_layers()),
      _total_layer_weights(_compute_total_layer_weights(_degree_by_layers))
{
init_admin();
}

RBConfigurationVertexPartitionWeightedLayers::~RBConfigurationVertexPartitionWeightedLayers()
{
 init_admin();
}

//CREATE METHODS
RBConfigurationVertexPartitionWeightedLayers* RBConfigurationVertexPartitionWeightedLayers::create(Graph* graph)
{
  return new RBConfigurationVertexPartitionWeightedLayers(graph, this->resolution_parameter);
}

RBConfigurationVertexPartitionWeightedLayers* RBConfigurationVertexPartitionWeightedLayers::create(Graph* graph, vector<size_t> const& membership)
{
  // TODO: verify we can remove this; degree_by_layer condensing should occur as a result of Graph::collapse_graph()
  // vector<vector<double> > new_degree_by_layer = this->_condense_degree_by_layer(membership);

  return new RBConfigurationVertexPartitionWeightedLayers(graph, membership, this->resolution_parameter);
}

//computes the 2 times the sum of all the edges within each layer
//used to compute the multilayer modularity
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

void RBConfigurationVertexPartitionWeightedLayers::_clear_resize(vector<vector<double> > &input_vec, size_t N, size_t M) {
  input_vec.clear();
  input_vec.resize(N);
  for (size_t i = 0; i < N; i++)
    input_vec[i].resize(M);
}

// helper method for totalling weight vectors
vector<double> RBConfigurationVertexPartitionWeightedLayers::add_vectors(vector<double> &v1, vector<double> &v2) {
  vector<double> outvec(v1.size(), 0);
  for (size_t i = 0; i < outvec.size(); i++)
    outvec[i] = v1[i] + v2[i];
  return outvec;
}

vector<double> RBConfigurationVertexPartitionWeightedLayers::multiply_vectors_elementwise(vector<double> &v1, vector<double> &v2) {
  vector<double> outvec(v1.size(), 0);
  for (size_t i = 0; i < outvec.size(); i++)
    outvec[i] = v1[i] * v2[i];
  return outvec;
}

vector<double> RBConfigurationVertexPartitionWeightedLayers::divide_vectors_elementwise(vector<double> &v1, vector<double> &v2) {
  vector<double> outvec(v1.size(), 0);
  for (size_t i = 0; i < outvec.size(); i++)
    outvec[i] = v1[i] / v2[i];
  return outvec;
}

double RBConfigurationVertexPartitionWeightedLayers::dot_product(vector<double> &v1, vector<double> &v2) {
  double outval = 0;
  for (size_t i = 0; i < v1.size(); i++)
    outval += v1[i] * v2[i];
  return outval;
}

// helper method for totalling weight vectors
vector<double> RBConfigurationVertexPartitionWeightedLayers::subtract_vectors(vector<double> &v1, vector<double> &v2) {
  vector<double> outvec(v1.size(), 0);
  for (size_t i = 0; i < outvec.size(); i++)
    outvec[i] = v1[i] - v2[i];
  return outvec;
}

// for multiply vector by scalar
vector<double> RBConfigurationVertexPartitionWeightedLayers::scalar_multiply(double scalar, vector<double> &v1) {
  vector<double> outvec(v1.size(), 0);
  for (size_t i = 0; i < v1.size(); i++)
    outvec[i] = scalar * v1[i];
  return outvec;
}

// sum all of the elements within a vector
double RBConfigurationVertexPartitionWeightedLayers::sum_over_vector(vector<double> &v1) {
  double output = 0;
  for (size_t i = 0; i < v1.size(); i++)
    output += v1[i];
  return output;
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
  for (size_t i = 0; i < this->_layer_vec.size() ; i++)
  {
    if (this->_layer_vec[i] + 1 > nb_layers)
      nb_layers = this->_layer_vec[i] + 1;
  }
  this->_nb_layers=nb_layers;
  // Reset administration
  this->community.clear();
  for (size_t i = 0; i < nb_comms; i++)
    this->community.push_back(new set<size_t>());

  this->_clear_resize(_total_weight_in_comm_by_layer,nb_comms,nb_layers);
  this->_clear_resize(_total_weight_to_comm_by_layer,nb_comms,nb_layers);
  this->_clear_resize(_total_weight_from_comm_by_layer,nb_comms,nb_layers);

  this->_csize.clear();
  this->_csize.resize(nb_comms);

  this->_current_node_cache_community_from = n + 1; this->_clear_resize(_cached_weight_from_community,n,nb_layers);
  this->_current_node_cache_community_to = n + 1;   this->_clear_resize(_cached_weight_to_community,n,nb_layers);
  this->_current_node_cache_community_all = n + 1;  this->_clear_resize(_cached_weight_all_community,n,nb_layers);

  this->_total_weight_in_all_comms_by_layer.resize(this->nb_layers());

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
    vector< double>  w_layers = this->graph->edge_layer_weights(e);


    // Add weight to the outgoing weight of community of v
    this->_total_weight_from_comm_by_layer[v_comm]=this->add_vectors(_total_weight_from_comm_by_layer[v_comm],w_layers);
    #ifdef DEBUG
      cerr << "\t" << "Add (" << v << ", " << u << ") weight " << w << " to from_comm " << v_comm <<  "." << endl;
    #endif
    // Add weight to the incoming weight of community of u
    this->_total_weight_from_comm_by_layer[u_comm] = this->add_vectors(_total_weight_from_comm_by_layer[u_comm],w_layers);
    #ifdef DEBUG
      cerr << "\t" << "Add (" << v << ", " << u << ") weight " << w << " to to_comm " << u_comm << "." << endl;
    #endif
    if (!this->graph->is_directed())
    {
      #ifdef DEBUG
        cerr << "\t" << "Add (" << u << ", " << v << ") weight " << w << " to from_comm " << u_comm <<  "." << endl;
      #endif
      this->_total_weight_from_comm_by_layer[u_comm] = this->add_vectors(_total_weight_from_comm_by_layer[u_comm],w_layers);
      #ifdef DEBUG
        cerr << "\t" << "Add (" << u << ", " << v << ") weight " << w << " to to_comm " << v_comm << "." << endl;
      #endif
      this->_total_weight_from_comm_by_layer[v_comm] = this->add_vectors(_total_weight_from_comm_by_layer[v_comm],w_layers);
    }
    // If it is an edge within a community
    if (v_comm == u_comm)
    {
      this->_total_weight_in_comm_by_layer[v_comm] = this->add_vectors(_total_weight_in_comm_by_layer[v_comm],w_layers);
      this->_total_weight_in_all_comms_by_layer = this->add_vectors(_total_weight_in_all_comms_by_layer,w_layers);
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

void RBConfigurationVertexPartitionWeightedLayers::move_node(size_t v,size_t new_comm)
{
  #ifdef DEBUG
    cerr << "void MutableVertexPartition::move_node(" << v << ", " << new_comm << ")" << endl;
    if (new_comm >= this->nb_communities())
      cerr << "ERROR: New community (" << new_comm << ") larger than total number of communities (" << this->nb_communities() << ")." << endl;
  #endif
  // Move node and update internal administration
  if (new_comm >= this->nb_communities())
  {
    if (new_comm < this->graph->vcount())
    {
      while (new_comm >= this->nb_communities())
        this->add_empty_community();
    }
    else
    {
      throw Exception("Cannot add new communities beyond the number of nodes.");
    }
  }

  // Keep track of all possible edges in all communities;
  size_t node_size = this->graph->node_size(v);
  size_t old_comm = this->_membership[v];
  #ifdef DEBUG
    cerr << "Node size: " << node_size << ", old comm: " << old_comm << ", new comm: " << new_comm << endl;
  #endif
  // Incidentally, this is independent of whether we take into account self-loops or not
  // (i.e. whether we count as n_c^2 or as n_c(n_c - 1). Be careful to do this before the
  // adaptation of the community sizes, otherwise the calculations are incorrect.
  if (new_comm != old_comm)
  {
    double delta_possible_edges_in_comms = 2.0*node_size*(ptrdiff_t)(this->_csize[new_comm] - this->_csize[old_comm] + node_size)/(2.0 - this->graph->is_directed());
    _total_possible_edges_in_all_comms += delta_possible_edges_in_comms;
    #ifdef DEBUG
      cerr << "Change in possible edges in all comms: " << delta_possible_edges_in_comms << endl;
    #endif
  }

  // Remove from old community
  #ifdef DEBUG
    cerr << "Removing from old community " << old_comm << ", community size: " << this->_csize[old_comm] << endl;
  #endif
  this->community[old_comm]->erase(v);
  this->_csize[old_comm] -= node_size;
  #ifdef DEBUG
    cerr << "Removed from old community." << endl;
  #endif

  // We have to use the size of the set of nodes rather than the csize
  // to account for nodes that have a zero size (i.e. community may not be empty, but
  // may have zero size).
  if (this->community[old_comm]->size() == 0)
  {
    #ifdef DEBUG
      cerr << "Adding community " << old_comm << " to empty communities." << endl;
    #endif
    this->_empty_communities.push_back(old_comm);
    #ifdef DEBUG
      cerr << "Added community " << old_comm << " to empty communities." << endl;
    #endif
  }

  if (this->community[new_comm]->size() == 0)
  {
    #ifdef DEBUG
      cerr << "Removing from empty communities (number of empty communities is " << this->_empty_communities.size() << ")." << endl;
    #endif
    vector<size_t>::reverse_iterator it_comm = this->_empty_communities.rbegin();
    while (it_comm != this->_empty_communities.rend() && *it_comm != new_comm)
    {
      #ifdef DEBUG
        cerr << "Empty community " << *it_comm << " != new community " << new_comm << endl;
      #endif
      it_comm++;
    }
    #ifdef DEBUG
      cerr << "Erasing empty community " << *it_comm << endl;
      if (it_comm == this->_empty_communities.rend())
        cerr << "ERROR: empty community does not exist." << endl;
    #endif
    if (it_comm != this->_empty_communities.rend())
      this->_empty_communities.erase( (++it_comm).base() );
  }

  #ifdef DEBUG
    cerr << "Adding to new community " << new_comm << ", community size: " << this->_csize[new_comm] << endl;
  #endif
  // Add to new community
  this->community[new_comm]->insert(v);
  this->_csize[new_comm] += this->graph->node_size(v);

  // Switch outgoing links
  #ifdef DEBUG
    cerr << "Added to new community." << endl;
  #endif

  // Use set for incident edges, because self loop appears twice
  igraph_neimode_t modes[2] = {IGRAPH_OUT, IGRAPH_IN};
  for (size_t mode_i = 0; mode_i < 2; mode_i++)
  {
    igraph_neimode_t mode = modes[mode_i];

    // Loop over all incident edges
    vector<size_t> const& neighbours = this->graph->get_neighbours(v, mode);
    vector<size_t> const& neighbour_edges = this->graph->get_neighbour_edges(v, mode);

    size_t degree = neighbours.size();

    #ifdef DEBUG
      if (mode == IGRAPH_OUT)
        cerr << "\t" << "Looping over outgoing links." << endl;
      else if (mode == IGRAPH_IN)
        cerr << "\t" << "Looping over incoming links." << endl;
      else
        cerr << "\t" << "Looping over unknown mode." << endl;
    #endif

    for (size_t idx = 0; idx < degree; idx++)
    {
      size_t u = neighbours[idx];
      size_t e = neighbour_edges[idx];

      size_t u_comm = this->_membership[u];
      // Get the weight of the edge
//      double w=this->graph->edge_weight(e)
      vector<double> w_layer = this->graph->edge_layer_weights(e);
      if (mode == IGRAPH_OUT)
      {
        // Remove the weight from the outgoing weights of the old community
        this->_total_weight_from_comm_by_layer[old_comm] = this->subtract_vectors(this->_total_weight_from_comm_by_layer[old_comm],w_layer);
        // Add the weight to the outgoing weights of the new community
        this->_total_weight_from_comm_by_layer[new_comm] = this->add_vectors(this->_total_weight_from_comm_by_layer[new_comm],w_layer);
        #ifdef DEBUG
          cerr << "\t" << "Moving link (" << v << "-" << u << ") "
               << "outgoing weight " << w
               << " from " << old_comm << " to " << new_comm
               << "." << endl;
        #endif
      }
      else if (mode == IGRAPH_IN)
      {
        // Remove the weight from the outgoing weights of the old community
        this->_total_weight_to_comm_by_layer[old_comm] = this->subtract_vectors(this->_total_weight_to_comm_by_layer[old_comm],w_layer);
        // Add the weight to the outgoing weights of the new community
        this->_total_weight_to_comm_by_layer[new_comm] = this->add_vectors(this->_total_weight_to_comm_by_layer[new_comm],w_layer);
        #ifdef DEBUG
          cerr << "\t" << "Moving link (" << v << "-" << u << ") "
               << "incoming weight " << w
               << " from " << old_comm << " to " << new_comm
               << "." << endl;
        #endif
      }
      else
        throw Exception("Incorrect mode for updating the admin.");
      // Get internal weight (if it is an internal edge)
      vector<double> int_weight = this->scalar_multiply(1.0/((this->graph->is_directed() ? 1.0 : 2.0)/( u == v ? 2.0 : 1.0)),w_layer);
      // If it is an internal edge in the old community
      if (old_comm == u_comm)
      {
        // Remove the internal weight
        this->_total_weight_in_comm_by_layer[old_comm] = this->subtract_vectors( this->_total_weight_in_comm_by_layer[old_comm],int_weight);
        this->_total_weight_in_all_comms_by_layer = this->subtract_vectors( _total_weight_in_all_comms_by_layer, int_weight);
        #ifdef DEBUG
          cerr << "\t" << "From link (" << v << "-" << u << ") "
               << "remove internal weight " << int_weight
               << " from " << old_comm << "." << endl;
        #endif
      }
      // If it is an internal edge in the new community
      // i.e. if u is in the new community, or if it is a self loop
      if ((new_comm == u_comm) || (u == v))
      {
        // Add the internal weight
        this->_total_weight_in_comm_by_layer[new_comm] = this->add_vectors( _total_weight_in_comm_by_layer[new_comm], int_weight);
        this->_total_weight_in_all_comms_by_layer = this->add_vectors(_total_weight_in_all_comms_by_layer, int_weight);
        #ifdef DEBUG
          cerr << "\t" << "From link (" << v << "-" << u << ") "
               << "add internal weight " << int_weight
               << " to " << new_comm << "." << endl;
        #endif
      }
    }
  }
  #ifdef DEBUG
    // Check this->_total_weight_in_all_comms
    vector <double> check_total_weight_in_all_comms_by_layer(self->nb_layers(),0);
    for (size_t c = 0; c < this->nb_communities(); c++)
      check_total_weight_in_all_comms_by_layer = this->add_vectors(check_total_weight_in_all_comms_by_layer, this->total_weight_in_comm_by_layer(c);
    cerr << "Internal _total_weight_in_all_comms=" << this->sum_over_vector(this->_total_weight_in_all_comms_by_layer_
         << ", calculated check_total_weight_in_all_comms=" << this->sum_over_vector(check_total_weight_in_all_comms_by_layer) << endl;
  #endif
  // Update the membership vector
  this->_membership[v] = new_comm;
  #ifdef DEBUG
    cerr << "exit MutableVertexPartition::move_node(" << v << ", " << new_comm << ")" << endl << endl;
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

  //if graph is directed
  vector <double> total_weight = this->scalar_multiply(2.0 - this->graph->is_directed(),this->_total_layer_weights);
//  vector <double> total_weight = this->_total_layer_weights;

  if (this->sum_over_vector(total_weight) == 0.0)
    return 0.0;
  if (new_comm != old_comm)
  {

    vector<double> w_to_old = this->weight_to_comm_by_layer(v, old_comm);
    vector<double> w_from_old = this->weight_from_comm_by_layer(v, old_comm);
    vector<double> w_to_new = this->weight_to_comm_by_layer(v, new_comm);
    vector<double> w_from_new = this->weight_from_comm_by_layer(v, new_comm);

    vector<double> k_out = this->graph->layer_strength_out[v];
    vector<double> k_in = this->graph->layer_strength_in[v];

    double self_weight = this->graph->node_self_weight(v);

    vector<double> K_out_old = this->total_weight_from_comm_by_layer(old_comm);

    vector<double> K_in_old = this->total_weight_to_comm_by_layer(old_comm);

    vector<double> K_out_new = this->add_vectors(this->total_weight_from_comm_by_layer(new_comm),k_out);

    vector<double> K_in_new = this->add_vectors(this->total_weight_to_comm_by_layer(new_comm), k_in);


    //vectorized versions
    double diff_old=this->sum_over_vector(w_to_old);
    diff_old-=this->resolution_parameter*this->sum_over_vector(this->divide_vectors_elementwise(this->multiply_vectors_elementwise(k_out,K_in_old),total_weight));
    diff_old+=(this->sum_over_vector(w_from_old));
    diff_old-=( this->resolution_parameter*this->sum_over_vector(this->divide_vectors_elementwise(this->multiply_vectors_elementwise(k_in,K_out_old),total_weight)) )

    double diff_new = this->sum_over_vector(w_to_new);
    diff_new -= this->resolution_parameter*this->sum_over_vector(this->divide_vectors_elementwise(this->multiply_vectors_elementwise(k_out,K_in_new),total_weight));
    diff_new += (this->sum_over_vector(w_from_new)+self_weight);//self_weight should be scalar
    diff_new -= ( this->resolution_parameter*this->sum_over_vector(this->divide_vectors_elementwise(this->multiply_vectors_elementwise(k_in,K_out_new),total_weight)) )


//    double diff_old = (w_to_old - this->resolution_parameter*k_out*K_in_old/total_weight) + \
//               (w_from_old - this->resolution_parameter*k_in*K_out_old/total_weight);

//    double diff_new = (w_to_new + self_weight - this->resolution_parameter*k_out*K_in_new/total_weight) + \
//               (w_from_new + self_weight - this->resolution_parameter*k_in*K_out_new/total_weight);

    diff = diff_new - diff_old;

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

  vector<double> m;
  if (this->graph->is_directed())
    m = this->graph->total_weight();
  else
    m = this->scalar_multiply(2,this->graph->total_weight());

  if (this->sum_over_vector(m) == 0)
    return 0.0;

  double mod=0.0;
  for (size_t c = 0; c < this->nb_communities(); c++)
  {
        vector <double> w = this->total_weight_in_comm_by_layer(c);
        vector <double> w_out = this->total_weight_from_comm_by_layer(c);
        vector <double> w_in = this->total_weight_to_comm_by_layer(c);
        #ifdef DEBUG
          size_t csize = this->csize(c);
          cerr << "\t" << "Comm: " << c << ", size=" << csize << ", w=" << w << ", w_out=" << w_out << ", w_in=" << w_in << "." << endl;
        #endif

        mod += this->sum_over_vector(this->subtract_vectors(w,this->divide_vectors_elementwise(this->multiply_vectors_elementwise(w_out,w_in),this->scalar_multiply(this->graph->is_directed() ? 1.0 : 4.0,m))));
    //    mod += w - resolution_parameter*w_out*w_in/((this->graph->is_directed() ? 1.0 : 4.0)*this->graph->total_weight());
  }
//  double q = mod; //factor of 2 should be account for
  double q = (2.0 - this->graph->is_directed())*mod;

  #ifdef DEBUG
    cerr << "exit double RBConfigurationVertexPartitionWeightedLayers::quality()" << endl;
    cerr << "return " << q << endl << endl;
  #endif
  return q;
}