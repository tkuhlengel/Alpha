/*
 * rundelaunay.h
 *
 *  Created on: Aug 23, 2012
 *      Author: winnen
 */

#ifndef RUNDELAUNAY_H_
#define RUNDELAUNAY_H_
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <stdio.h>
#include "vertex.h"
#include "Geom.h"
#include "mathfunc.h"
#include "set.h"
	


using namespace std;

/** Executes a command line call given in the parameter
 \param cmd The full command that needs to be called
 \return The output from the command to stdout

 Code courtesy of the users of stackoverflow (http://stackoverflow.com/a/478960)
 */
std::string exec(char * cmd);

/** Constructs the command line call to qHull from the names provided
 \param fileinput Name of the input file to be read from.
 \param params Options to be passed to qDelaunay.
 \return The string output of qDelaunay.

 */
std::string buildCommand(std::string fileinput, std::string params = "Qt o");
std::string buildAndPrint(int argc, char** arg);

Alpha::Vertex* parsePoints(std::string pyOutput);


set_head_t* parseLinks(std::string qhullOutput);

std::string getNext(std::string buffer, int& pos, bool eol=false);



/** Parses the double-type coordinate data from a string line and returns C++ containers
 *
 * \return A double array of length dim (if radii is false), or dim+1 (if radii is true).
 */
double* parsePoints(std::string pyline /**<The output from the python script pdb2coord using -s option*/
		,int dim /**<Number of dimensions to the data (not including radius of the sphere)*/
		,bool radii /**<Indicates if there is a radius column after the coordinates*/
		);
/** Function to break the text stream down into data structures
 *MAKE SURE THIS FAILS IF A RADIUS IS NOT FOUND
 */
void parseData(std::string qhullOutput, std::string pyOutput);






#endif /* RUNDELAUNAY_H_ */
