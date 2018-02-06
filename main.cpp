// main.cpp
//
// ICS 46 Spring 2017
// Project #4: Rock and Roll Stops the Traffic
//
// This is the program's main() function, which is the entry point for your
// console user interface.

#include "RoadMap.hpp"
#include "RoadMapReader.hpp"
#include "RoadMapWriter.hpp"
#include "InputReader.hpp"
#include "Digraph.hpp"
#include "TripReader.hpp"
#include "Trip.hpp"
#include <iostream>


int main()
{

RoadMap newMap;
RoadMapReader reader;

InputReader ireader(std::cin);


RoadMapWriter rwriter;
rwriter.writeRoadMap(std::cout, reader.readRoadMap(ireader));
	
    return 0;
}

