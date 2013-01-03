/*
 * vertex.cpp
 *
 *  Created on: Jul 19, 2012
 *      Author: winnen
 */
#include "vertex.h"
namespace Alpha {

bool doubleEqual(double a, double b){
	/*Compares two doubles to see if they are within tolerances of each other
	 *
	 */
	return fabs(b-a)<=SMALLNUM;
}
void Vertex::allocArray() {
	if (this->arry == NULL) {
		//this->arry = new double[3];
		std::cout<<"We've got big problems here"<<std::endl;
	}
}
Vertex::Vertex() {
	//allocArray();
	for (int i = 0; i < 3; i++) {
		arry[i] = 0.0;
	}
}

Vertex::Vertex(double x, double y, double z) {
	//allocArray();
	arry[0] = x;
	arry[1] = y;
	arry[2] = z;
}
Vertex::Vertex(double * arr) {
	//allocArray();
	arry[0] = arr[0];
	arry[1] = arr[1];
	arry[2] = arr[2];
}
Vertex::Vertex(const Vertex * orig) {
	//allocArray();
	for (int i = 0; i < 3; i++) {
		arry[i] = orig->arry[i];
	}

}

Vertex::~Vertex() {
	//delete this.array;
	//delete[] arry;
}

double* Vertex::getPos() {
	//Returns a heap copy of the location array as a three member
	// double array for use with other code.
	//double result[3] = {arry[0],arry[1],arry[2]};
	//TODO: Figure out how to keep this clean.
	return arry;
}

double Vertex::fastDistance(const Vertex * a, const Vertex * b) {
	// Returns the square of the coordinate difference between
	// vector a and vector b, (x1-x2)^2 + ...
	double sum = 0.0;
	if (b->arry != NULL) {
		double x = 0.0;
		for (int i = 0; i < 3; i++) {
			x = a->arry[i] - b->arry[i];
			sum += (x * x);
		}
	}
	return sum;
}
double Vertex::qDist(const Vertex *b) const {
	//Quick distance, or the square of the difference between two positions
	return Vertex::fastDistance(this, b);
}
double Vertex::getX() const {
	//Get the x coordinate of the vector
	return arry[0];
}

double Vertex::getY() const {
	//Get the Y coordinate of the vector
	return arry[1];
}

double Vertex::getZ() const {
	//Get the Z coordinate of the vector
	return arry[2];
}

double Vertex::distance(const Vertex *b) const {
	//Returns the absolute distance to vector B.
	if (b->arry != NULL && arry != NULL) {
		return sqrt((double) this->qDist(b));
	} else {
		return -1;
	}
}

bool Vertex::gt(const Vertex* other, int dimension) const {
	/*Function primarily for sorting on a given dimension
	 Returns true if this vector is greater than the vector
	 passed in the parameters.
	 "dimension" is the number of the XYZ dimension that the
	 comparison is being tested on, with X = 0
	 */
	return (this->arry[dimension] > other->arry[dimension]);
}

void Vertex::operator =(const Vertex& other){
	allocArray();
	for (int i=0;i<3;i++){
		arry[i]=other.arry[i];
	}
}
bool Vertex::equals(Vertex * other){
	bool isEqual=true;
	for (int i=0; i<3; i++){
		if (!doubleEqual(arry[i], other->arry[i])){
			isEqual=false;
		}
	}
	return isEqual;
}
bool Vertex::operator ==(Vertex * other) {
	return this->equals(other);
}

double Vertex::operator [](int index) const {
	//Returns the axis coordinate of the respective axis, X,Y, or Z
	//Index must be between 0 and 2;  Error code 50 thrown otherwise;
	if (!((index >= 0) && (index <= 2))) {
		/*std::error(1,51,
		"Overflow Error: parameter index to Vertex::operator [] must be "
		"between 0 and 2.\n Other values such as %d will cause an overflow",
		index);*/
		throw 51;
	}
	return arry[index];
	
}
Vertex* Vertex::minus(const Vertex * b) const{
	Vertex *result=new Vertex(arry[0] - b->arry[0],
			arry[1] - b->arry[1],
			arry[2] - b->arry[2]);
	//std::cout<<"A =  "<<this<<"\nB =  "<<b<<"\nB-A= "<<result<<std::endl;
	return result;
}
Vertex* Vertex::plus(const Vertex * b) const{
	Vertex *result=new Vertex(arry[0] + b->arry[0],
			arry[1] + b->arry[1],
			arry[2] + b->arry[2]);
	//std::cout<<"A =  "<<this<<"\nB =  "<<b<<"\nB-A= "<<result<<std::endl;
	return result;
}
std::ostream& operator <<(std::ostream& out, Alpha::Vertex* that) {
	//char x[100];
	//sprintf(x,"{%.2f, %.2f, %.2f}",that->getX(),that->getY(), that->getZ());
	//out<<x;
	out << that->to_string();
	return out;
}

std::string Vertex::to_string() {
	//%[flags][width][.precision][length]
	char intermediate[50];
	sprintf(intermediate, "%6.3f, %6.3f, %6.3f", this->arry[0], this->arry[1], this->arry[2]);
	return intermediate;
}
std::ostream& operator <<(std::ostream& out, Alpha::Vertex &that) {
	out << that.to_string();
	return out;
}
}
