#include <iostream>
#include <fstream> // std::ifstream
#include <sstream> // std::stringstream
#include <string> // std::string, std::stoi
#include <limits.h> // UINT_MAX
#include <cstring> // std::strcmp
#include <cmath>
#include <vector>
#include <chrono>
#include <ctime>
#include "Graph.hpp"

using namespace CSC586C::optimize_algorithm;

extern const double damping_factor = 0.85;
extern const unsigned max_iterations = 100;
extern const double tolerance = 1e-8;

// Read Input from file with format:
// N (#vertex)
// src_index dest_index ... src_index dest_index (pairs of source and destination links)
std::vector<ColdEdge> ReadInputFromTextFile(const char* input_file, unsigned& num_vertices)
{
    std::vector<ColdEdge> input;
    std::string input_line, str;
    std::ifstream myfile (input_file);
    if (myfile.is_open())
    {
        while (getline(myfile,input_line))
        {
           if (input_line.size() < std::to_string(UINT_MAX).size())
           {
                // Record the number of vertexes
               num_vertices = std::stoi(input_line);
               // std::cout << input_line << std::endl;
           } 
           else
           {
                std::stringstream temp_line(input_line);
                unsigned source = 0;
                unsigned destination = 0;
                unsigned count = 0;
                while (getline(temp_line, str, ' '))
                {
                    ++count;
                    unsigned temp_value = std::stoi(str);
                    if (count == 1)
                    {
                        source = temp_value;
                    }
                    else if (count == 2)
                    {
                        count = 0;
                        destination = temp_value;
                        // Writing out the edges into the graph
                        // Example: std::vector<Edge> input{ Edge{1,0}, Edge{4,2}, Edge{2,5}, Edge{1,3}};
                        if (source == destination) break; 
                        input.push_back(ColdEdge{source, destination});  
                        // std::cout << source << std::endl;
                        // std::cout << destination << std::endl;
                    }
                }
           }
        }
        myfile.close();
    }
    return input;
}

bool ToleranceCheck(const unsigned& num_v, std::vector<HotData>& nodes)
{
    // Sum up the pagerank
    double pr_sum = 0.0;
    for (unsigned i = 0; i < num_v; i++) 
    {
        pr_sum += nodes[i].pagerank;
    }
    // Calculate the cur_toleranceor
    pr_sum = 1.0 / pr_sum;
    double cur_tolerance = 0.0;
    for (unsigned i = 0; i < num_v; i++)
    {
        nodes[i].pagerank *= pr_sum;
        // norm 1
        cur_tolerance += std::fabs(nodes[i].pagerank - nodes[i].pre_pagerank);
    }

    if (cur_tolerance < tolerance)
    {
        std::cout << "Current toleranceor: " << cur_tolerance << std::endl;
        return true;
    }
    return false;
}

void PageRank(Algorithm_Graph *graph)
{
    const unsigned num_v = graph->VertexesNum();
    double init_rank = double(1.0 / num_v);
    double pr_random = (1.0 - damping_factor) / num_v;
    
    for (unsigned i = 0; i < num_v; i++)
    {
        graph->nodes[i].pagerank = init_rank;
        graph->nodes[i].pre_pagerank = 0.0;
    }

    unsigned iter = 0;
    while (iter++ < max_iterations)
    {
        // Update the pagerank values in every iteration
        for (unsigned i = 0; i < num_v; i++)
        {
            graph->nodes[i].pre_pagerank = graph->nodes[i].pagerank;
            graph->nodes[i].pagerank = 0.0;
        }

        // Distribute the pr_sum of all dangling nodes(no outer edges) to all nodes.
        double dangling_pr_sum = 0.0;
        for (unsigned i = 0; i < num_v; i++)
        {
           // Remove branch prediction
            dangling_pr_sum += graph->nodes[i].pre_pagerank * (graph->adjEdges[i].size() == 0);
        }
        double pr_dangling = damping_factor * dangling_pr_sum / num_v;

        // Iterater all the vertexes and calculate its adjacency function l(pi,pj) of all inward links
        // Update its pagerank value by adding pr_eigenvector from its inward links separately
        for (unsigned i = 0; i < num_v; i++)
        {
            unsigned outgoing_edges_num = graph->adjEdges[i].size();
            double pr_eigenvector = damping_factor * graph->nodes[i].pre_pagerank / outgoing_edges_num;
            for (unsigned edge_index = 0; edge_index < outgoing_edges_num; edge_index++)
            {
                unsigned update_edge = graph->adjEdges[i].at(edge_index);
                graph->nodes[update_edge].pagerank += pr_eigenvector;
            }
            graph->nodes[i].pagerank += pr_random + pr_dangling;
        }
       

        // finish when cur_toleranceor is smaller than tolerance we set
        if(ToleranceCheck(num_v, graph->nodes)) 
        {
            std::cout << "Iteration time: " << iter << std::endl;
            break;
        }
    }
}

void PrintBenchmark(std::chrono::time_point<std::chrono::steady_clock> start_t, std::chrono::time_point<std::chrono::steady_clock> const end_t, const int loop_t)
{
    auto const avg_time = std::chrono::duration_cast<std::chrono::microseconds>( end_t - start_t ).count() / double(loop_t);
    std::cout << "Average total running time  = " << avg_time << " us" << std::endl;
}

int main(int argc, char *argv[])
{
    if(argc == 3)
    {
        int loop_times = 10;
        unsigned num_vertices = 0;
        const char* test_mode = argv[2];

        std::vector<ColdEdge> input = ReadInputFromTextFile(argv[1], num_vertices);

        if(std::strcmp(test_mode, "total") == 0)
        {
            auto const start_time = std::chrono::steady_clock::now();

            for (int i = 0; i < loop_times; i++)
            {
                Algorithm_Graph graph(num_vertices, input);

                PageRank(&graph);
            }  
            auto const end_time = std::chrono::steady_clock::now(); 
            PrintBenchmark(start_time, end_time, loop_times);
        }
        else if(std::strcmp(test_mode, "graph") == 0 )
        {
            auto const start_time = std::chrono::steady_clock::now();
            Algorithm_Graph graph(num_vertices, input);
            auto const end_time = std::chrono::steady_clock::now(); 

            PageRank(&graph);  
            PrintBenchmark(start_time, end_time, 1);          
        }
        else if(std::strcmp(test_mode, "pagerank") == 0)
        {
            Algorithm_Graph graph(num_vertices, input);
            auto const start_time = std::chrono::steady_clock::now();
            for (int i = 0; i < loop_times; i++)
            {
                PageRank(&graph);
            }
            auto const end_time = std::chrono::steady_clock::now(); 
            PrintBenchmark(start_time, end_time, loop_times);
        }
        else
        {
            std::cout << "Invalid Input!" << std::endl;
        }
    }
    else
    {
        std::cout << "Invalid Input: " << std::endl;
        std::cout << "Please input the input text file name wanted in argc[1]" << std::endl;
        std::cout << "Please input the time mode(total/graph/pangerank) to be record in argc[2]" << std::endl;
    }

    return 0;
}
