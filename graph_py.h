#ifndef UTIL_GRAPH_PY_H_
#define UTIL_GRAPH_PY_H_

#include <boost/python.hpp>
#include "graph.h"


class PyIntGraph : public WeightedGraph<int>
{
  public:
    void PyAddEdges(boost::python::list& py_edges){
      int num_edges = len(py_edges);
      for (int i = 0; i < num_edges; ++i)
      {
        Vertex v_s = boost::python::extract<Vertex>(py_edges[i][0]);
        Vertex v_t = boost::python::extract<Vertex>(py_edges[i][1]);
        int w = 1;
        this->AddEdge(v_s, v_t, w);
      } 
    }
};


/* export as a python module */

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
//  ReadEdgeListDefault, PyIntGraph::ReadWeightedEdgeList, 1, 2);

BOOST_PYTHON_MODULE(cpp_graph)
{
  using namespace boost::python;

  class_<PyIntGraph>("int_graph")
    .def("add_edges", &PyIntGraph::PyAddEdges)
//    .def("read_weighted_edge_list", &PyIntGraph::ReadWeightedEdgeList, ReadEdgeListDefault(args("filename")))
    .def("read_weighted_edge_list", &PyIntGraph::ReadWeightedEdgeList)
    .def("read_edge_list", &PyIntGraph::ReadEdgeList)
    .def("invert_weight_sign", &PyIntGraph::InvertWeightSign)
    .def("make_adjacent_vertex_vector", &PyIntGraph::MakeAdjacentVertexVector)
    .def("make_adjacent_vertex_vector_directed", &WeightedGraph<int>::MakeAdjacentVertexVectorDirected)
    .def("min_vertex", &PyIntGraph::min_vertex)
    .def("max_vertex", &PyIntGraph::max_vertex)
    .def("num_vertices", &PyIntGraph::num_vertices)
    .def("num_edges", &PyIntGraph::num_edges)
    .def("sum_weights", &PyIntGraph::sum_weights)
    .def("degree", &PyIntGraph::degree)
    .def("ave_degree", &PyIntGraph::ave_degree);
}

#endif