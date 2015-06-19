#ifndef UTIL_GRAPH_H_
#define UTIL_GRAPH_H_

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <limits>

#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/range/irange.hpp>

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

typedef size_t Vertex;

struct Edge
{
  Vertex source;
  Vertex target;
  
  Edge(){}

  Edge(Vertex s, Vertex t) :
    source(s), target(t)
  {}

  Edge(const Edge &e) :
    source(e.source), target(e.target)
  {}
};

template<typename Weight>
struct WeightedVertex
{
  Vertex vertex;
  Weight weight;
  
  WeightedVertex(){}

  WeightedVertex(Vertex v, Weight w) :
    vertex(v), weight(w)
  {}

  WeightedVertex(const WeightedVertex &w_v) :
    vertex(w_v.vertex), weight(w_v.weight)
  {};
};

template<typename Weight>
struct WeightedEdge
{
  Edge edge;
  Weight weight;

  WeightedEdge(){}

  WeightedEdge(Edge e, Weight w) :
    edge(e), weight(w)
  {}

  WeightedEdge(const WeightedEdge &w_e) :
    edge(w_e.edge), weight(w_e.weight)
  {}
};

BOOST_FUSION_ADAPT_STRUCT( 
  Edge,
  (Vertex, source)
  (Vertex, target)
)

BOOST_FUSION_ADAPT_TPL_STRUCT(
  (Weight),
  (WeightedEdge) (Weight),
  (Edge, edge)
  (Weight, weight)
)


template<typename Weight>
using WeightedVertexVector = std::vector< WeightedVertex<Weight> >;
template<typename Weight>
using WeightedEdgeList = std::list< WeightedEdge<Weight> >;

std::ostream& operator<<(std::ostream &os, const Edge &e);


template<typename Weight>
std::ostream& operator<<(std::ostream &os, const WeightedEdge<Weight> &w_e);

template<typename Weight>
std::ostream& operator<<(std::ostream &os, const WeightedVertexVector<Weight> &w_edge_vec);

template<typename Weight>
std::ostream& operator<<(std::ostream &os, const WeightedVertex<Weight> &w_v);

template<typename Weight>
std::ostream& operator<<(std::ostream &os, const WeightedEdgeList<Weight> &w_vertex_list);


template<typename Weight>
class WeightedGraph
{
  typedef typename WeightedVertexVector<Weight>::const_iterator const_w_vertex_iterator;
  typedef typename WeightedEdgeList<Weight>::const_iterator const_w_edge_iterator;

  public:
    WeightedGraph() :
      min_vertex_(std::numeric_limits<Vertex>::max()),
      max_vertex_(std::numeric_limits<Vertex>::min())
    {}
    ~WeightedGraph(){}
    void AddEdge(Vertex vertex_a, Vertex vertex_b, Weight weight);
    void AddEdge(Edge edge, Weight weight);
    void AddEdge(WeightedEdge<Weight> w_edge);
    // template<typename Iter>
    // void AddEdges(Iter beg, Iter end){
    //   for (Iter i = beg; i != end; ++i)
    //   {
    //     this->AddEdge(*i);
    //   }
    // }
    // template<typename Iter>  //  for no weight case
    // void AddEdges(Iter beg, Iter end, Weight w){
    //   for (Iter i = beg; i != end; ++i)
    //   {
    //     this->AddEdge(*i, w);
    //   }
    // }
    void ReadEdgeList(const std::string& filename);
    void ReadWeightedEdgeList(const std::string& filename
      //, bool adjust_vertex_index_p = true
      );
    void MakeAdjacentVertexVector();
    void MakeAdjacentVertexVectorDirected();
    void InvertWeightSign();

    boost::integer_range<Vertex> vertices() const
    {
      return boost::irange(min_vertex_, max_vertex_ + 1);
    }
    std::pair<const_w_edge_iterator, const_w_edge_iterator> edges() const
    {
      return make_pair(weighted_edge_list_.begin(), weighted_edge_list_.end());
    }
    std::pair<const_w_vertex_iterator, const_w_vertex_iterator> adjacent_vertices(const Vertex v) const
    {
      return make_pair(adjacent_weighted_vertices_vector_[v].begin(),
        adjacent_weighted_vertices_vector_[v].end());
    };

//    void OutputAdjacencyMatrix(Eigen::Matrix<Weight, Eigen::Dynamic, Eigen::Dynamic> *out_matrix);
    Vertex min_vertex() const {return min_vertex_;}
    Vertex max_vertex() const {return max_vertex_;}
    size_t num_vertices() const {return max_vertex_ - min_vertex_ + 1;}
    size_t num_edges() const {return weighted_edge_list_.size();}

    Weight sum_weights() const;

    size_t degree(const Vertex i) const {return adjacent_weighted_vertices_vector_[i].size();}
    double ave_degree() const {
      return 2.0*num_edges()/num_vertices();
    };

  private:
    void AdjustVertexIndex();
    void UpdateMinMaxVertex(Vertex vertex);

    // member variables
    WeightedEdgeList<Weight> weighted_edge_list_;
    Vertex min_vertex_;
    Vertex max_vertex_;
    std::vector< WeightedVertexVector<Weight> > adjacent_weighted_vertices_vector_;
};


// template functions

template<typename Weight>
std::ostream& operator<<(std::ostream &os, const WeightedEdge<Weight> &w_e)
{
    os << w_e.edge << " " << w_e.weight;
    return os;
}

template<typename Weight>
std::ostream& operator<<(std::ostream &os, const WeightedEdgeList<Weight> &w_edge_list)
{
  for (typename WeightedEdgeList<Weight>::const_iterator i = w_edge_list.begin();
       i != w_edge_list.end(); ++i)
  {
    os << (*i) << "\n";
  }
  return os;
}

template<typename Weight>
std::ostream& operator<<(std::ostream &os, const WeightedVertex<Weight> &w_v)
{
  os << w_v.vertex << " " << w_v.weight;
  return os;
}

template<typename Weight>
std::ostream& operator<<(std::ostream &os, const WeightedVertexVector<Weight> &w_vertex_list)
{
  for (typename WeightedVertexVector<Weight>::const_iterator i = w_vertex_list.begin();
       i != w_vertex_list.end(); ++i)
  {
    os << (*i) << "\n";
  }
  return os;
}


template<typename Weight>
void WeightedGraph<Weight>::AddEdge(Vertex vertex_a, Vertex vertex_b, Weight weight)
{
  Edge edge = Edge(vertex_a, vertex_b);
  WeightedEdge<Weight> w_edge = WeightedEdge<Weight>(edge, weight);
  weighted_edge_list_.emplace_back(w_edge);

  this->UpdateMinMaxVertex(vertex_a);
  this->UpdateMinMaxVertex(vertex_b);
}

template<typename Weight>
void WeightedGraph<Weight>::AddEdge(Edge edge, Weight weight)
{
  WeightedEdge<Weight> w_edge = WeightedEdge<Weight>(edge, weight);
  weighted_edge_list_.emplace_back(w_edge);

  this->UpdateMinMaxVertex(edge.source);
  this->UpdateMinMaxVertex(edge.target);
}

template<typename Weight>
void WeightedGraph<Weight>::AddEdge(WeightedEdge<Weight> w_edge)
{
  weighted_edge_list_.emplace_back(w_edge);

  this->UpdateMinMaxVertex(w_edge.edge.source);
  this->UpdateMinMaxVertex(w_edge.edge.target);
}

template <typename Iterator, typename Weight>
struct weighted_edge_grammer : qi::grammar<Iterator, WeightedEdge<Weight>(), ascii::space_type>
{
  weighted_edge_grammer() : weighted_edge_grammer::base_type(w_edge)
  {
    using qi::uint_;
    using qi::auto_;
//    using qi::lexeme;
    using qi::omit;
    using ascii::space;

//    edge %= uint_ >> omit[+space] >> uint_;
//    w_edge %= lexeme[edge >> omit[+space] >> auto_];
    edge %= uint_ >> uint_;
    w_edge %= edge >> auto_;
  }

  qi::rule<Iterator, WeightedEdge<Weight>(), ascii::space_type> w_edge;
  // qi::rule<Iterator, Edge()> edge;
  qi::rule<Iterator, Edge(), ascii::space_type> edge;
};


template<typename Weight>
void WeightedGraph<Weight>::ReadWeightedEdgeList(const std::string& filename
  //, bool adjust_vertex_index_p
  )
{
  std::ifstream ifs;
  ifs.open(filename);
  if ( !ifs.is_open() )
  {
    std::cerr << "Cannot open the file: " << filename << std::endl;
    std::exit(1);
  }

  std::string buf;
  weighted_edge_grammer<std::string::iterator, Weight> grammer;

  while (getline(ifs, buf))
  {
    std::string::iterator itr = buf.begin();
    std::string::iterator end = buf.end();

    WeightedEdge<Weight> w_edge;
    bool rst = phrase_parse(itr, end, grammer, ascii::space, w_edge);
    
    if (rst && itr == end) {
      this->AddEdge(w_edge);
    }
  }

  ifs.close();

  // if (min_vertex_ != 0 && adjust_vertex_index_p)
  // {
  AdjustVertexIndex();
  // }
}

template <typename Iterator>
struct edge_grammer : qi::grammar<Iterator, Edge(), ascii::space_type>
{
  edge_grammer() : edge_grammer::base_type(edge)
  {
    using qi::uint_;
    edge %= uint_ >> uint_;
  }

  qi::rule<Iterator, Edge(), ascii::space_type> edge;
};

template<typename Weight>
void WeightedGraph<Weight>::ReadEdgeList(const std::string& filename)
{
  std::ifstream ifs;
  ifs.open(filename);
  if ( !ifs.is_open() )
  {
    std::cerr << "Cannot open the file: " << filename << std::endl;
    std::exit(1);
  }

  std::string buf;
  edge_grammer<std::string::iterator> grammer;

  while (getline(ifs, buf))
  {
    std::string::iterator itr = buf.begin();
    std::string::iterator end = buf.end();

    Edge edge;
    bool rst = phrase_parse(itr, end, grammer, ascii::space, edge);
    WeightedEdge<Weight> w_edge(edge, 1);

    if (rst && itr == end) {
      this->AddEdge(w_edge);
    }
  }

  ifs.close();

  AdjustVertexIndex();
}

template<typename Weight>
void WeightedGraph<Weight>::MakeAdjacentVertexVector()
{
  std::vector<int> degree_vec(max_vertex_ + 1, 0);
  for (typename WeightedEdgeList<Weight>::iterator itr = weighted_edge_list_.begin();
    itr != weighted_edge_list_.end(); ++itr)
  {
    Vertex vertex_a = (*itr).edge.source;
    Vertex vertex_b = (*itr).edge.target;
    ++degree_vec[vertex_a];
    ++degree_vec[vertex_b];
  }

  adjacent_weighted_vertices_vector_.resize(max_vertex_ + 1, WeightedVertexVector<Weight>());
  for (Vertex v = min_vertex_; v <= max_vertex_; ++v){
    adjacent_weighted_vertices_vector_[v].reserve(degree_vec[v]);
  }

  for (typename WeightedEdgeList<Weight>::iterator itr = weighted_edge_list_.begin();
    itr != weighted_edge_list_.end(); ++itr)
  {
    Vertex vertex_a = (*itr).edge.source;
    Vertex vertex_b = (*itr).edge.target;
    Weight weight = (*itr).weight;

    adjacent_weighted_vertices_vector_[vertex_a].push_back(WeightedVertex<Weight>(vertex_b, weight));
    adjacent_weighted_vertices_vector_[vertex_b].push_back(WeightedVertex<Weight>(vertex_a, weight));
  }
}

template<typename Weight>
void WeightedGraph<Weight>::MakeAdjacentVertexVectorDirected()
{
  std::vector<int> degree_vec(max_vertex_ + 1, 0);
  for (typename WeightedEdgeList<Weight>::iterator itr = weighted_edge_list_.begin();
    itr != weighted_edge_list_.end(); ++itr)
  {
    Vertex vertex_a = (*itr).edge.source;
    ++degree_vec[vertex_a];
  }

  adjacent_weighted_vertices_vector_.resize(max_vertex_ + 1, WeightedVertexVector<Weight>());
  for (Vertex v = min_vertex_; v <= max_vertex_; ++v){
    adjacent_weighted_vertices_vector_[v].reserve(degree_vec[v]);
  }

  for (typename WeightedEdgeList<Weight>::iterator itr = weighted_edge_list_.begin();
    itr != weighted_edge_list_.end(); ++itr)
  {
    Vertex vertex_a = (*itr).edge.source;
    Vertex vertex_b = (*itr).edge.target;
    Weight weight = (*itr).weight;

    adjacent_weighted_vertices_vector_[vertex_a].push_back(WeightedVertex<Weight>(vertex_b, weight));
  }
}


template<typename Weight>
void WeightedGraph<Weight>::InvertWeightSign()
{
    for (typename WeightedEdgeList<Weight>::iterator itr = weighted_edge_list_.begin();
        itr != weighted_edge_list_.end(); ++itr)
    {
        Weight &w = (*itr).weight;
        w = -w;
    }
}

template<typename Weight>
Weight WeightedGraph<Weight>::sum_weights() const
{
  Weight weight_sum = 0;
  for (typename WeightedEdgeList<Weight>::const_iterator itr = weighted_edge_list_.begin();
    itr != weighted_edge_list_.end(); ++itr)
  {
    weight_sum += (*itr).weight;
  }
  return weight_sum;
}

template<typename Weight>
void WeightedGraph<Weight>::AdjustVertexIndex()
{
  for (typename WeightedEdgeList<Weight>::iterator itr = weighted_edge_list_.begin();
        itr != weighted_edge_list_.end(); ++itr)
  {
      Vertex &v_a = (*itr).edge.source;
      Vertex &v_b = (*itr).edge.target;
      v_a -= min_vertex_;
      v_b -= min_vertex_;
  }
  max_vertex_ -= min_vertex_;
  min_vertex_ = 0;
}

template<typename Weight>
void WeightedGraph<Weight>::UpdateMinMaxVertex(Vertex vertex)
{
  if (vertex < min_vertex_)
  {
    min_vertex_ = vertex;
  }
  else if (vertex > max_vertex_)
  {
    max_vertex_ = vertex;
  }
}



#endif