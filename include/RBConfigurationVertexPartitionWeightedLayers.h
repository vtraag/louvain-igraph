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


  protected:

    virtual void init_admin();

  private:
    void _clear_resize(vector<vector<double > > &input_vec, size_t N, size_t M);
    vector<double> add_vectors(vector<double> &v1, vector<double> &v2);
    vector<size_t> const _layer_vec;
    vector<vector<double> > const _degree_by_layers;
    vector<double> const _total_layer_weights;
    vector<vector<double> > _condense_degree_by_layer(vector<size_t> const& membership);
    vector<double> _compute_total_layer_weights(vector<vector<double> > const& degree_by_layers);

    // Keep track of the internal weight of each community for each layer
    vector<vector<double> > _total_weight_in_comm_by_layer;
    // Keep track of the total weight to a community
    vector<vector<double> > _total_weight_to_comm_by_layer;
    // Keep track of the total weight from a community
    vector<vector<double> > _total_weight_from_comm_by_layer;
    // Keep track of the total internal weight
    vector <double> _total_weight_in_all_comms_by_layer;

    vector<size_t> _empty_communities;

    void cache_neigh_communities_by_layer(size_t v, igraph_neimode_t mode);

    //these are being overloaded
    size_t _current_node_cache_community_from; vector<vector<double> > _cached_weight_from_community; vector<size_t> _cached_neigh_comms_from;
    size_t _current_node_cache_community_to;   vector<vector<double> > _cached_weight_to_community;   vector<size_t> _cached_neigh_comms_to;
    size_t _current_node_cache_community_all;  vector<vector<double> > _cached_weight_all_community;  vector<size_t> _cached_neigh_comms_all;


};

#endif // RBConfigurationVertexPartitionWeightedLayers_H
