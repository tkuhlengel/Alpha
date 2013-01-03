/*
 * parser.cpp
 *
 *  Created on: Jun 29, 2012
 *      Author: winnen
 */

#include "parser.h"
#include "Atom.h"
/*
Node::Node() {
	head = NULL;
	tail = NULL;
	data = NULL;
}
template <class T>
Node::Node(Node *head, Node *tail, T *obj) {
	this->head=head;
	this->tail=tail;
	this->data=obj;
}
template<class T>
T* Node::getItem(int counter, int desired) const {
	if (counter == desired) {
		return this->data;
	} else if (counter < desired) {
		if (this->tail != NULL) {
			return this->tail->getItem(counter + 1, desired);
		} else
			return NULL;
	}
}
template<class T>
void Node::append(T *item) {
	if (this->tail != NULL) {
		this->tail->append(item);
	} else {
		this->tail = new Node;
		this->tail->data = item;
	}
}
};
*/


/*
 ifstream& openFile(string filename){
 ifstream pdb;
 if (filename.length()>0){

 if (pdb.is_open() && pdb.good()){

 }
 }
 else{
 throw "InputError: Please specify a filename to read from.";
 }
 return pdb;
 }
 */
/*
void readPDB(const string name) {
	//Loads the PDB into memory and splits out the ATOM data for use in the components of the program
	//output is an Atom vector which stores all the atoms in the structure
	//Atom *first=new Atom();
	char buffer[100];
	//TODO figure this out
	//vector<*Atom> numbers(1000, first, Atom);
	string atm = "ATOM  ";
	string data;
	fstream input;
	input.open(name.data(), ifstream::in);
	while (!input.eof() && input.good()) {
		input.getline(buffer, 100);
		if (atm.compare(0, 5, buffer) != 0) { //If its not an ATOM line, i probably don't care
			continue;
		} else {
			data = string(buffer);

		}
	}

}

string clean(const string &str) {
//Removes whitespace from the beginning and end of the string
string result;
bool first = false;
int start = 0, end = str.length(), len = str.length();
const char *dat = str.data();
for (int i = 0; i < str.length(); i++) {
	if (' ' == dat[i] && !first) {
		start = i + 1;
	} else if (!first) {
		first = true;
	}
	if (first && ' ' != dat[i]) {
		end = i;
	}
}
return str.substr(start, (end - start));
}
*/
