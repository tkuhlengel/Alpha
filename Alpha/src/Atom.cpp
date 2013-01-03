/*
 * Atom.cpp
 *
 *  Created on: Jun 29, 2012
 *      Author: winnen
 */

#include "Atom.h"
namespace Alpha {
Alpha::Atom::Atom() {
	// TODO Auto-generated constructor stub
	location = new Vertex();
	setRadius();
}

Alpha::Atom::Atom(double x, double y, double z, std::string element) {
	this->location = new Vertex(x, y, z);
	this->element = element;
	setRadius();
}
Alpha::Atom::Atom(std::string buffer) { //Comments are Field name and index in the string
	this->rname = buffer.substr(0, 6); //Record name 0-5
	this->serial = atoi(buffer.substr(6, 5).data()); //Atom serial number 6-10
	this->atom = buffer.substr(12, 4); //Atom name, 4 chars 12-15
	//Empty char 16
	this->residue = buffer.substr(17, 3); //residue name (3 chars) 17-19
	this->chain = buffer[21]; //Chain identifier 21
	this->rsequence = atoi(buffer.substr(22, 4).data()); //Residue sequence number INT 22-25
	//Skipping insertion code 26
	//Using Vertex class to store coordinates
	this->location = new Vertex(atof(buffer.substr(30, 8).data()), //30-37 = X coordinate in Angstroms 8.3f
	atof(buffer.substr(38, 8).data()), //38-45 = Y Coordinate in Angstroms 8.3f
	atof(buffer.substr(46, 8).data()) //46-53 = Z Coordinate in Angstroms 8.3f
			);
	this->occupancy = atof(buffer.substr(54, 6).data()); //54-59
	this->tfactor = atof(buffer.substr(60, 6).data()); //Temperature factor, defaults to 0.0 Real 6.2 60-65
	this->segIdent = buffer.substr(72, 4); //Segment Identifier 72-75
	this->element = buffer.substr(76, 2); //Element symbol right justified 76-77
	this->charge = buffer.substr(78, 6); //Charge on the atom 78-79

	/*TODO: Make a lookup table for the Van der Waals radius using the above information
	 *  Will require use of element, charge, and occupancy to determine proper radius.
	 *  May have to make use of Gromacs source code to do this*/
	if (element.length() > 0) {
		setRadius();
	}

}
Alpha::Atom::~Atom() {
	// TODO Auto-generated destructor stub
	delete location;
}
void Alpha::Atom::setRadius() {
	//Based on the Element entry of the PDB, sets the radius to match.
	std::string elements = "CHONPS";
	int pos = elements.find(this->element, 0);
	switch (pos) {
	case 0: //Carbon
		this->radius = 1.7;
		break;
	case 1: //Hydrogen
		this->radius = 1.2;
		break;
	case 2: //Oxygen
		this->radius = 1.52;
		break;
	case 3: //Nitrogen
		this->radius = 1.55;
		break;
	case 4: //Phosphorus
		//this->radius=1.8;
	case 5: //Sulfur
		this->radius = 1.8;
		break;

		//TODO update this, because it needs more values.
	default: //Set it to the median value if it is not known.
		this->radius = 1.8;
		break;
	}
}
double Alpha::Atom::Radius() const {
	return this->radius;
}

double Alpha::Atom::operator [](int index) const {
	return (*(this->location))[index];
}
Alpha::Vertex* Alpha::Atom::getCenter() const {
	return this->location;
}
/*
char* cleanspaces(char *strn) {
//Removes whitespace from the beginning and end of the string
	std::string result;
	std::string str=std::string(strn);
	bool first = false;
	int start = 0, end = str.length(), len = str.length();
	;
	for (int i = 0; i < len; i++) {
		if (' ' == str.data()[i] && !first) {
			start = i + 1;
		} else if (!first) {
			first = true;
		}
		if (first && ' ' != str.data()[i]) {
			end = i;
		}
	}
	//str=
	return str.substr(start, (end - start)).data();
}*/
}
