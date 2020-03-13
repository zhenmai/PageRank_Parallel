#include <iostream>
#include <fstream> // std::ifstream
#include <sstream> // std::stringstream
#include <string> // std::string, std::stoi
#include <limits.h> // UINT_MAX
#include <cstring> // std::strcmp
#include <cmath>  // std::fabs
#include <numeric>  // std::accumulate()
#include <vector>
#include <chrono>  // timing libraries
#include <omp.h>
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

bool ToleranceCheck(std::vector<HotData>& nodes)
{
    const unsigned& num_v = nodes.size();

    // Normalize every pagerank value to make the probabilities sum to 1
    double pr_sum = 0.0;
    #pragma omp parallel for reduction(+:pr_sum)
    for (unsigned i = 0; i < num_v; i++) 
    {
        pr_sum += nodes[i].pagerank;
    }
    pr_sum = 1.0 / pr_sum;

    #pragma omp parallel for
    for (unsigned i = 0; i < num_v; i++)
    {
        nodes[i].pagerank *= pr_sum;
    }

    // Calculate the cur_toleranceor 
    auto const num_threads = omp_get_max_threads();
    std::vector<double> cur_tolerance(num_threads, 0.0);
    #pragma omp parallel for
    for (unsigned i = 0; i < num_v; i++)
    {
        // norm 1
        auto const t = omp_get_thread_num();
        cur_tolerance[t] += std::fabs(nodes[i].pagerank - nodes[i].pre_pagerank);
    }
    double tolerance_sum = std::accumulate(std::begin(cur_tolerance), std::end(cur_tolerance), 0.0 );
    
    // If we meet our tolerance then we break
    if (tolerance_sum < tolerance)
    {
        std::cout << "Current toleranceor: " << tolerance_sum << std::endl;
        return true;
    }
    return false;
}

void PageRank(Algorithm_Graph *graph)
{
    const unsigned num_v = graph->VertexesNum();
    double init_rank = double(1.0 / num_v);
    double pr_random = (1.0 - damping_factor) / num_v;
    
    #pragma omp parallel for
    for (unsigned i = 0; i < num_v; i++)
    {
        graph->nodes[i].pagerank = init_rank;
        graph->nodes[i].pre_pagerank = 0.0;
    }

    unsigned iter = 0;
    while (iter++ < max_iterations)
    {
        #pragma omp parallel for
        // Update the pagerank values in every iteration
        for (unsigned i = 0; i < num_v; i++)
        {
            graph->nodes[i].pre_pagerank = graph->nodes[i].pagerank;
            graph->nodes[i].pagerank = 0.0;
        }

        // Distribute the pr_sum of all dangling nodes(no outer edges) to all nodes.
        double dangling_pr_sum = 0.0;
        #pragma omp parallel for reduction(+:dangling_pr_sum)
        for (unsigned i = 0; i < num_v; i++)
        {
           // Remove branch prediction
            dangling_pr_sum += graph->nodes[i].pre_pagerank * (graph->adjEdges[i].size() == 0);
        }
        double pr_dangling = damping_factor * dangling_pr_sum / num_v;

        // Iterater all the vertexes and calculate its adjacency function l(pi,pj) of all inward links
        // Update its pagerank value by adding pr_eigenvector from its inward links separately
        // In this step we have done pr[p]← pr[p]
        #pragma omp parallel for
        for (unsigned i = 0; i < num_v; i++)
        {
            unsigned outgoing_edges_num = graph->adjEdges[i].size();
            for (unsigned edge_index = 0; edge_index < outgoing_edges_num; edge_index++)
            {
                double pr_eigenvector = damping_factor * graph->nodes[i].pre_pagerank / outgoing_edges_num;
                unsigned update_edge = graph->adjEdges[i].at(edge_index);
                #pragma omp atomic
                graph->nodes[update_edge].pagerank += pr_eigenvector;
            }
        }

        // Calculate the total pagerank value by adding the pr[p]← 1-d/N + pr_dangling/N
        #pragma omp parallel for
        for (unsigned i = 0; i < num_v; ++i)
        {
            graph->nodes[i].pagerank += pr_random + pr_dangling;
        }

        // finish when cur_toleranceor is smaller than tolerance we set
        if(ToleranceCheck(graph->nodes)) 
        {
            std::cout << "Total iteration time: " << iter << std::endl;
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
    int loop_times = 10;
    unsigned num_vertices = 0;

    if(argc >= 4)
    {
        std::vector<ColdEdge> input = ReadInputFromTextFile(argv[ 1 ], num_vertices);
        const char* test_mode = argv[ 2 ];
        omp_set_num_threads( atoi( argv[ 3 ] ) );
        std::cout << "The number of threads used: " << omp_get_max_threads() << std::endl;

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
    else if (argc >= 2 && argc < 4)
    {
        std::vector<ColdEdge> input = ReadInputFromTextFile(argv[ 1 ], num_vertices);
        std::cout << "The number of threads used: " << omp_get_max_threads() << std::endl;
        auto const start_time = std::chrono::steady_clock::now();
        for (int i = 0; i < loop_times; i++)
        {
            Algorithm_Graph graph(num_vertices, input);

            PageRank(&graph);
        }  
        auto const end_time = std::chrono::steady_clock::now(); 
        PrintBenchmark(start_time, end_time, loop_times);
    }
    else
    {
        std::cout << "Invalid Input: " << std::endl;
        std::cout << "Please input the input text file name wanted in argc[1]" << std::endl;
        std::cout << "Please input the time mode(total/graph/pangerank) to be record in argc[2]" << std::endl;
        std::cout << "Please input the number of threads wanted to use in argc[3]" << std::endl;
    }

    return 0;
}
