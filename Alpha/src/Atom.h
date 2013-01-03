/*
 * Atom.h
 *
 *  Created on: Jun 29, 2012
 *      Author: Trevor Kuhlengel
 */

#ifndef ATOM_H_
#define ATOM_H_

#include <iostream>
#include <string>
#include <map>
#include <cstdlib>

#include "vertex.h"
//#include "parser.cpp"
//global static string RADIUS[]={"H", "C", "N", "O", "F", "P", "SD"};

namespace Alpha{
//char* cleanspaces(char* str);

class Atom {
	friend std::ostream& operator <<(std::ostream&, const Atom&);
private:
	void setRadius(); //Sets the radius according to the element's Van der Waals radius
protected:
	Vertex *location;
	std::string rname; //record name 0-5
	int serial; //Atom serial number 6-10

	std::string atom; //Atom name, 4 chars 12-15
	//Empty char 16
	std::string residue; //residue name (3 chars) 17-19
	char chain; //Chain identifier 21
	int rsequence; //Residue sequence number 22-25
	//Skipping insertion code 26
	//30-37 = X coordinate in Angstroms 8.3f
	//38-45 = Y Coordinate in Angstroms 8.3f
	//46-53 = Z Coordinate in Angstroms 8.3f
	double occupancy;  //54-59
	double tfactor; //Temperature factor, defaults to 0.0 Real 6.2 60-65
	std::string segIdent; //Segment Identifier 72-75

	std::string element; //Element symbol right just 76-77
	std::string charge; //Charge on the atom 78-79
	double radius; //Van der Waals radius


public:
	Atom();
	Atom(double x, double y, double z, std::string element);
	Atom(std::string buffer);
	virtual ~Atom();
	double Radius() const;//Returns the Van der Waals radius of the atom
	double operator [](int index) const; //Returns the x,y, or z coordinate.
		//index must be a value from 0-2 0=X, 1=Y, 2=Z
	Vertex* getCenter() const;

};

}
#endif /* ATOM_H_ */
