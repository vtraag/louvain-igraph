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

  protected:
  private:
    vector<size_t> const _layer_vec;
    vector<vector<double> > const _degree_by_layers;
    vector<double> const _total_layer_weights;
    vector<vector<double> > _condense_degree_by_layer(vector<size_t> const& membership);
    vector<double> _compute_total_layer_weights(vector<vector<double> > const& degree_by_layers);
};

#endif // RBConfigurationVertexPartitionWeightedLayers_H
