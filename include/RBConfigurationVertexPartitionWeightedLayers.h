#ifndef RBConfigurationVertexPartitionWeightedLayers_H
#define RBConfigurationVertexPartitionWeightedLayers_H

#include "LinearResolutionParameterVertexPartition.h"

class RBConfigurationVertexPartitionWeightedLayers : public LinearResolutionParameterVertexPartition
{
  public:
    RBConfigurationVertexPartitionWeightedLayers(Graph* graph,
      vector<size_t> const& _layer_vec, vector<vector <double> > const& _degree_by_layers,
      vector<size_t> const& membership, double resolution_parameter);
    RBConfigurationVertexPartitionWeightedLayers(Graph* graph,
      vector<size_t> const& _layer_vec, vector<vector <double> > const& _degree_by_layers,
      vector<size_t> const& membership);
    RBConfigurationVertexPartitionWeightedLayers(Graph* graph,
      vector<size_t> const& _layer_vec, vector<vector <double> > const& _degree_by_layers,
      double resolution_parameter);
    RBConfigurationVertexPartitionWeightedLayers(Graph* graph,
      vector<size_t> const& _layer_vec, vector<vector <double> > const& _degree_by_layers);

    virtual ~RBConfigurationVertexPartitionWeightedLayers();
    virtual RBConfigurationVertexPartitionWeightedLayers* create(Graph* graph);
    virtual RBConfigurationVertexPartitionWeightedLayers* create(Graph* graph, vector<size_t> const& membership);

    virtual double diff_move(size_t v, size_t new_comm);
    virtual double quality(double resolution_parameter);
    virtual void move_node(size_t v,size_t new_comm); //override this to change how calculations are stored

    vector <double> weight_to_comm_by_layer(size_t v, size_t comm);
    vector <double> weight_from_comm_by_layer(size_t v, size_t comm);

    inline vector<double> total_weight_in_comm_by_layer(size_t comm) { return this->_total_weight_in_comm_by_layer[comm]; };
    inline vector<double> total_weight_from_comm_by_layer(size_t comm) { return this->_total_weight_from_comm_by_layer[comm]; };
    inline vector<double> total_weight_to_comm_by_layer(size_t comm) { return this->_total_weight_to_comm_by_layer[comm]; };
    inline vector<double> total_weight_in_all_comms_by_layer() { return this->_total_weight_in_all_comms_by_layer; };
    inline size_t nb_layers() {return this->_nb_layers;};

  protected:

    virtual void init_admin();

  private:
    //methods
    void _clear_resize(vector<vector<double > > &input_vec, size_t N, size_t M);
    vector <double> add_vectors(vector<double> &v1, vector<double> &v2);
    vector <double> subtract_vectors(vector<double> &v1, vector<double> &v2);
    vector <double> multiply_vectors_elementwise(vector<double> &v1, vector<double> &v2);
    vector <double> divide_vectors_elementwise(vector<double> &v1, vector<double> &v2);
    vector <double> scalar_multiply(double scalar, vector<double> &v1);
    double sum_over_vector(vector<double> &v1);
    double dot_product(vector<double> &v1, vector<double> &v2);
    vector<vector<double> > _condense_degree_by_layer(vector<size_t> const& membership);
    vector<double> _compute_total_layer_weights(vector<vector<double> > const& degree_by_layers);
    void cache_neigh_communities_by_layer(size_t v, igraph_neimode_t mode);

    //variables
    size_t _nb_layers;
    vector<size_t> const _layer_vec;
    vector<vector<double> > const _degree_by_layers;
    vector<double>  _total_layer_weights;




    // Keep track of the internal weight of each community for each layer
    vector<vector<double> > _total_weight_in_comm_by_layer;
    // Keep track of the total weight to a community
    vector<vector<double> > _total_weight_to_comm_by_layer;
    // Keep track of the total weight from a community
    vector<vector<double> > _total_weight_from_comm_by_layer;
    // Keep track of the total internal weight
    vector <double> _total_weight_in_all_comms_by_layer;
    size_t _total_possible_edges_in_all_comms;

    vector<size_t> _empty_communities;


    //these are the same as the inherited class
    size_t _current_node_cache_community_from; vector<vector<double> > _cached_weight_from_community; vector<size_t> _cached_neigh_comms_from;
    size_t _current_node_cache_community_to;   vector<vector<double> > _cached_weight_to_community;   vector<size_t> _cached_neigh_comms_to;
    size_t _current_node_cache_community_all;  vector<vector<double> > _cached_weight_all_community;  vector<size_t> _cached_neigh_comms_all;


};

#endif // RBConfigurationVertexPartitionWeightedLayers_H
