#include <iostream>
#include <fstream>  // std::ofstream
#include <random> 	// std::rand, std::srand, std::default_random_engine
#include <algorithm> // std::generate

// from benchmarking.cpp in class
template < typename RandomGenerator >
	auto build_rand_vec( RandomGenerator gen, size_t const size )
	{
		// see timing.hpp for a description of decltype().
		// Initialises a vector with the appropriate size, in which we will create the random data.
		std::vector< decltype((gen)()) > random_data( size );

		// std::generate() = STL. This could be done with a loop from 0,...,size, calling gen() each
		// iteration. See https://www.fluentcpp.com/2017/01/05/the-importance-of-knowing-stl-algorithms/
		// for a compelling discussion of why using STL algorithms like this is better than explicit loops.
		// Spoiler: it better self-documents (in code) intent and is less prone to off-by-one-errors.
		// For an explanation of syntax, see the use of std::for_each() in timing.hpp
		std::generate( random_data.begin(), random_data.end(), gen );

		return random_data;
	}

int main (int argc, char *argv[])
{
	if(argc == 4)
	{
		unsigned num_vertexes = std::stoi(argv[1]);
		unsigned test_size = std::stoi(argv[2]) * 2;
		std::ofstream myfile (argv[3]);

		if(num_vertexes == 0 || test_size == 0 || !myfile.is_open())
		{
			std::cout << "Invalid Input! " << std::endl;
			return 0;
		}

		if (myfile.is_open())
	 	{
	  	 	myfile << num_vertexes << "\n";

	  	 	auto output = build_rand_vec([](){return static_cast< unsigned >( std::rand());}, test_size);
															
	  	 	for(unsigned i = 0; i < output.size(); i++)
	  	 	{
	  	 		// output the size of random generation data between [0, num_vertexes)
	  	 		myfile << output[i] % num_vertexes << ' ';
	  	 	}
	  	  	myfile.close();
		}
  	}
  	else
  	{
  		std::cout << "Invalid Input: " << std::endl;
  		std::cout << "Please input the number of vertexed wanted in argc[1]" << std::endl;
  		std::cout << "Please input the size of graph links wanted in argc[2]" << std::endl;
  		std::cout << "Please input the output text file name wanted in argc[3]" << std::endl;
  	}
  
  	return 0;
}