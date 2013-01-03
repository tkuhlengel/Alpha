/**
 *
 * Author: Trevor Kuhlengel
 *  Date: June 17, 2012
 *  Purpose: A simple 3 point vector class for use with the alpha shape
 *  generating program
 */
#ifndef VERTEX_H_
#define VERTEX_H_
#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
//#include <error.h>
namespace Alpha{
//#define X_AXIS=0
//#define Y_AXIS=1
//#define Z_AXIS=2
#define SMALLNUM 0.0000001
//Tolerance equality
bool doubleEqual(double a, double b);

class Vertex {
	friend std::ostream& operator <<(std::ostream& out,	Vertex* that);
	//friend std::ostream& operator <<(std::ostream& out,	Vertex* that);
private:
	void allocArray();
protected:
	

public:
	double arry[3];
	Vertex();
	Vertex(double x, double y, double z);
	Vertex(double* arr);
	Vertex(const Vertex* orig);
	~Vertex();

	double* getPos();
	//Returns a the location array as a three member

	double static fastDistance(const Vertex* a, const Vertex* b);
	// Returns the square of the coordinate difference between
	// vector a and vector b, (x1-x2)^2 + ...

	double qDist(const Vertex * b) const;
	//Quick distance, or the square of the difference between two positions
	//

	double getX() const;
	//Get the x coordinate of the vector

	double getY() const;
	//Get the Y coordinate of the vector

	double getZ() const;
	//Get the Z coordinate of the vector

	double distance(const Vertex* b) const;
	//Returns the absolute distance to vector B.

	bool gt(const Vertex* other, int dimension = 0) const;
	/*Function primarily for sorting on a given dimension
	 Returns true if this vector is greater than the vector
	 passed in the parameters.
	 "dimension" is the number of the XYZ dimension that the
	 comparison is being tested on, with X = 0
	 */
	void operator =(const Vertex& other); //Copies other onto this vertex
	bool equals(Vertex* other);
	bool operator ==(Vertex* other);

	/**Subtraction method
	 * \return The difference between two vertices, as a Vertex storing a vector.
	 */
	Vertex* minus(const Vertex* b) const;
	Vertex* plus(const Vertex* b) const;

	double operator [](int index) const;
	//Returns the axis coordinate of the respective axis, X,Y, or Z
	//Index must be between 0 and 2;  Error code 50 thrown otherwise;	std::ostream& operator <<(std::ostream&, Alpha::Vertex&);
	std::string to_string();

};

}
#endif
