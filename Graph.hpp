
#include <iostream>
#include <vector>
#include <list>

namespace CSC586C {

// The first edition of data structure
// Baseline_graph with AoS and use std::vector<Vertex*> to store the edges and their relationships
namespace baseline_graph{
// data structure to store graph edges
struct Edge
{
	unsigned src, dest;
};

class Vertex
{
public:
   Vertex(unsigned index_) : index(index_) 
   {
      page_rank = 0.0;
      pre_pagerank = 0.0;
      num_inward_edges = 0;
      num_outward_edges = 0;
      neighbors.resize(0);
   }
   void addNeighbor(Vertex* vertex) 
   {
      neighbors.push_back(vertex);
      ++num_inward_edges;
   }
   double pageRank() const{
      return page_rank;
   }
   double prePageRank() const{
      return pre_pagerank;
   }
   void setPageRank(double page_rank_){
      page_rank = page_rank_;
   }
   void setPrePageRank(double page_rank_){
      pre_pagerank = page_rank_;
   }
   unsigned numInwardEdges() const{
      return num_inward_edges;
   }
   unsigned numOutwardEdges() const{
      return num_outward_edges;
   }
   unsigned getIndex() const{
      return index;
   }
   void outCountIncrease() {
      ++num_outward_edges;
   }
   // For debugging only
   void printVertexInfo() const
   { 
      std::cout << "Vertex index: " << index << " outward count " << num_outward_edges
                << "page rank " << page_rank << '\n'
                << "previous page rank " << pre_pagerank <<'\n';
      std::cout << "Neighbors: ";
      for(std::list<Vertex*>::const_iterator it=neighbors.begin(); it!= neighbors.end(); ++it) 
      {
         std::cout << (*it)->index << " ";
      }
      std::cout << '\n';
   }
   std::list<Vertex*> neighbors; //list of inward neighbors (sources)
private:
   double page_rank;
   double pre_pagerank;
   unsigned num_inward_edges;
   unsigned num_outward_edges;
   unsigned index;
};

class Graph
{
public:
   Graph(unsigned num_vertices_, std::vector<Edge> edges) : num_vertices(num_vertices_){
      vertices.resize(0);
      //Initialize all the vertices on this graph
      for(int i=0; i<num_vertices; ++i){
         Vertex* vertex = new Vertex(i);
         vertices.push_back(vertex);
      }
      //Add all the edges info
      for(int i=0; i<edges.size(); ++i){
         vertices[edges[i].dest]->addNeighbor(vertices[edges[i].src]);
         vertices[edges[i].src]->outCountIncrease();
      }
   }
   unsigned VertexesNum() const{
      return num_vertices;
   }
   void printGraphInfo() const 
   { //For debugging only
      std::cout << "==================PRINT GRAPH=====================\n";
      std::cout << "Number of vertices: "<<num_vertices<<'\n';
      double pagerank_sum = 0.0;
      for(int i=0; i<num_vertices; ++i) {
         vertices[i]->printVertexInfo();
         pagerank_sum += vertices[i]->pageRank();
      }
      std::cout << "Total pagerank sum is " << pagerank_sum << '\n';
      std::cout << "=========================END======================\n";
   }
   ~Graph() {
      for(int i=0; i<num_vertices; ++i){
         Vertex* vertex  = vertices[i];
         delete vertex;
      }
      vertices.resize(0);
      num_vertices = 0;
   }
   std::vector<Vertex*> vertices; //All vertices in this graph
private:
   unsigned num_vertices;
};
}// End of namespace baseline_graph





// The second edition of data structure
// Used cold/hot data and std::vector<std::vector<int> > to store the edges and their relationships
namespace optimize_graph {
// data structure to store cold data: edges (pairs of src and dest)
struct ColdEdge 
{
	unsigned int src;
	unsigned int dest;
};
// data structure to store hot data: number of outgoing links each node and its pagerank values 
struct HotData
{
	double pagerank;
	double pre_pagerank;
	int outgoing_edges_num;
};
// class to represent a graph object
class Optimized_Graph
{
public:
	// construct a vector of vectors to represent each edge index and its inward edges (hot data)
	std::vector<std::vector<int> > adjEdges;
	// construct a vector to store hot data structure
	std::vector<HotData> nodes;

	// Graph Constructor
	Optimized_Graph(unsigned Num, std::vector<optimize_graph::ColdEdge> &input): vertex_num(Num), edges(input)
	{
		adjEdges.resize(Num);
		nodes.resize(Num);
		// add edges to the directed graph
		for (auto &edge: edges)
		{
			adjEdges[edge.dest].push_back(edge.src);
			nodes[edge.src].outgoing_edges_num++;
		}
	}
	unsigned VertexesNum() {return vertex_num;}

   //For debugging only
   void printGraph() const 
   {
      std::cout << "Number of vertices is: " << vertex_num <<'\n';
      std::cout << "Hot data: \n";
      std::cout << "Pre_pagerank: ";
      for(int i = 0; i < vertex_num; ++i){
         std::cout << nodes[i].pre_pagerank << ' ';
      }
      std::cout << '\n' << "Pagerank: ";
      for(int i = 0; i < vertex_num; ++i){
         std::cout << nodes[i].pagerank << ' ';
      }
      std::cout << '\n' << "OutgoingEdgeCount: ";
      for(int i = 0; i < vertex_num; ++i){
         std::cout << nodes[i].outgoing_edges_num << ' ';
      }
      std::cout << '\n' << "The adj matrix:\n";
      for(int i = 0; i < vertex_num; ++i){
         for(int j = 0; j < adjEdges[i].size(); ++j){
            std::cout << adjEdges[i][j] << "--->" << i << " ";
         }
         std::cout<<'\n';
      }
   }
   void printHotData() const 
   {
      std::cout << "Hot data: \n";
      std::cout << "Pre_pagerank: ";
      for(int i = 0; i < vertex_num; ++i){
         std::cout << nodes[i].pre_pagerank <<' ';
      }
      std::cout << '\n' << "Pagerank: ";
      for(int i = 0; i < vertex_num; ++i){
         std::cout << nodes[i].pagerank << ' ';
      }
   }

private:
	unsigned int vertex_num;	
	// construct a vector to store cold data structure
	std::vector<optimize_graph::ColdEdge> const &edges;
};
} // End of namespace optimize_graph





// The third edition of data structure
// Used SoA and std::vector<std::vector<int> > to store the edges and their relationships
namespace soa_graph {

// data structure to store cold data: edges (pairs of src and dest)
struct ColdEdge 
{
	std::vector<unsigned> src;
   std::vector<unsigned> dest;
};
// data structure to store hot data: number of outgoing links each node and its pagerank values 
struct HotData
{
	std::vector<double> pagerank;
	std::vector<double> pre_pagerank;
   std::vector<int> outgoing_edges_num; 
};
// class to represent a graph object
class SoA_Graph
{
public:
	// construct a vector of vectors to represent each edge index and its inward edges (hot data)
	std::vector<std::vector<int> > adjEdges;
	// construct a data structure to store hot data arrays
	HotData hotData;

	// Graph Constructor
	SoA_Graph(unsigned Num, const ColdEdge &input): vertex_num(Num), edges(input)
	{
		adjEdges.resize(Num);
      hotData.pagerank.resize(Num);
      hotData.pre_pagerank.resize(Num);
      hotData.outgoing_edges_num.resize(Num);
		// add edges to the directed graph
      for (auto i = 0; i < edges.src.size(); i++)
      {
         adjEdges[edges.dest[i]].push_back(edges.src[i]);
         hotData.outgoing_edges_num[edges.src[i]]++;
      }
	}
	unsigned VertexesNum() {return vertex_num;}
   //For debugging only
   void printGraph() const 
   {
      std::cout << "Number of vertices is: " << vertex_num <<'\n';
      std::cout << "Hot data: \n";
      std::cout << "Pre_pagerank: ";
      for(int i = 0; i < vertex_num; ++i){
         std::cout << hotData.pre_pagerank[i] << ' ';
      }
      std::cout << '\n' << "Pagerank: ";
      for(int i = 0; i < vertex_num; ++i){
         std::cout << hotData.pagerank[i] << ' ';
      }
      std::cout << '\n' << "OutgoingEdgeCount: ";
      for(int i = 0; i < vertex_num; ++i){
         std::cout << hotData.outgoing_edges_num[i] << ' ';
      }
      std::cout << '\n' << "The adj matrix:\n";
      for(int i = 0; i < vertex_num; ++i){
         for(int j = 0; j < adjEdges[i].size(); ++j){
            std::cout << adjEdges[i][j] << "--->" << i << " ";
         }
         std::cout<<'\n';
      }
   }
   void printHotData() const 
   {
      std::cout << "Hot data: \n";
      std::cout << "Pre_pagerank: ";
      for(int i = 0; i < vertex_num; ++i){
         std::cout << hotData.pre_pagerank[i] <<' ';
      }
      std::cout << '\n' << "Pagerank: ";
      for(int i = 0; i < vertex_num; ++i){
         std::cout << hotData.pagerank[i] << ' ';
      }
   }

private:
	unsigned int vertex_num;	
	// construct data structure to store cold data arrays
	soa_graph::ColdEdge edges;
};
} // End of namespace soa_graph

} // End of namespace CSC586C
