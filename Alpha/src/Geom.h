/** \file Geom.h
 * \brief Class container for Geometrical objects in cartesian coordinates.
 * 
 * Geometry classes for abtraction from the other classes of geometric shapes used
 * in the Alpha shapes program.
 *
 * \author Trevor Kuhlengel
 *  Institution: Lehigh University
 */

#ifndef GEOM_H_
#define GEOM_H_

#include <iostream>
#include <cmath>
#include "vertex.h"
#include "Atom.h"
#include <cstring>
#include "mathfunc.h"
//#include "triIntersect.h"

/** \namespace Alpha 
 * 	\brief Namespace container to separate classes and functions used in Alpha
 * 		shapes calculations from others, due to a number of name collisions.
 */
namespace Alpha {
#define PROBE_SIZE = 1.0


/** \class Geom
 * \brief An abstract base class for the rest of the Geometry objects in this project.
 * 
 * 
 */
class Geom {
public:
	virtual bool ptInside(const Vertex*) =0;
	virtual Vertex* segtest(Vertex * A, Vertex * B){
		std::cout<<"Error: Unexpected use of Geom::segTest"<<std::endl;
		return NULL;
	};
	/**< Virtual function for handling the unexpected bass class segtest call.
	 *
	 * @param A[in] Vertex pointer denoting the start of a line.
	 * @param B[in] Second Vertex pointer B which denote a line from A->B.
	 * @return Vertex indicating the intersection of the Geom object
	 * 		closest to B farthest from A that intersects as a double vector, or
	 *  will return NULL if no intersections are found.  Tests to see if there is
	 *  an intersection of the line denoted  by a->b.  Will return the vector
	 *  location of the farthest point from A that intersects the object.
	 *  */
	virtual ~Geom(){};


};
/** @class Sphere
 * @brief Class to represent a spherical structure.  
 * 
 * Represents a sphere used in Alpha Shape generator, which can determine
 * if material is inside or outside the sphere rapidly.
 * 
 * @see Geom
 */
class Sphere: public Geom {
	

protected:
	Vertex *center;
	double radius;
public:
	/**@brief Empty Constructor
	 * @return Sets all values to 0, and pointer to center to NULL.
	 */
	Sphere();
	/**@brief Constructor
	 * @param center[in] A Vertex object. Is copied for storage. 
	 * @return Sphere object copies the center passed in and stores the copy, and pointer to center to NULL.
	 */
	Sphere(const Vertex* center, double radius);
	
	/** @deprecated  Due to the unfinished implementation of Atom and the other available PDB loaders, 
	 * 		this constructor should not be used.
	 */
	Sphere(const Alpha::Atom * a);
	
	/**Copy constructor
	 */
	Sphere(const Sphere*);
	
	/** Destructor */
	~Sphere();
	
	/** Set the radius of the sphere
	 * @param radius[in] The value of the new radius.
	 */
	void setRadius(double radius);
	
	/** Set or change the center point of the Sphere
	 * 
	 * Memory safe replacement of the center of the sphere. 
	 * @param center[in] The desired new center point to be copied into the Sphere. 
	 */ 
	void setCenter(const Vertex* center);
	
	/** Radius accession method.
	 * 
	 * @return The radius of the sphere
	 */
	double getRadius() const;
	
	/** Center accession method.
	 * 
	 * @return A mutable pointer to the center of the sphere.
	 */
	Vertex* getCenter() const;
	
	/** @brief Enclosing space tester. 
	 * 
	 * Tests whether or not the input point is within the confines of the 
	 * 	structure. Points ON the boundary of the surface are considered to be
	 * 	inside the surface.
	 * 
	 * Inputs: Two vectors A and B which denote a line from A->B.
	 * 
	 * @param A[in] Vertex object pointer denoting the start point of a line.
	 * 		preferred if point A is inside the object, but not required.
	 * @param B[in] Vertex object pointer denoting an endpoint along the 
	 * line A->B
	 * @return A boolean denoting if the point is inside the Sphere
	 */
	bool ptInside(const Vertex*);
	
	/** @brief Line-sphere intersection calculator.
	 * 
	 * Calculates the intersection point of a line A->B with the surface of the
	 * 3-dimensional structure.  If more than one point is present, it returns
	 * the point that is closest to point B (second Vertex parameter).
	 * 
	 * @param A[in] A Vertex denoting the start of line A->B. 
	 * @param B[in] A Vertex denoting the end of the line A->B.
	 * @return If the line A->B intersects the surface of the structure, 
	 * a heap pointer to a Vertex is returned.
	 * 
	 * @note It is preferred that one of the points is inside the structure and returns true when 
	 * ptInside() is called. However, if this is not the case, the function the point is known
	 */
	Vertex* segTest(Vertex* A, Vertex* B);

	/** Counts the number of intersections of the line A->B with the Sphere.
	 * 
	 * Returns the number of times the segment intersects the sphere
	 * retrieving multiple intersection points
	 * 
	 * @param A[in] A Vertex denoting the start of line A->B. 
	 * @param B[in] A Vertex denoting the end of the line A->B.
	 * 
	 * @note Used in Sphere::segTest and Cup::segTest.
	 */
	int intersectionCount(const Vertex *A, const Vertex *B);
	
	/** @brief Quadratic Root Solver 
	 *
	 * Calculates the quadratic roots in the form
	 * \f
	 * a*x^2 + b*x + c = 0
	 * \f
	 * using the quadratic formula: 
	 * \f
	 * b\pm\sqrt{(b^2-4*a*c)}
	 * /(2*a)
	 * \f
	 * @note Used by Sphere::segTest and Cup::segTest.
	 * @note T
	 */
	double* getQuadraticRoots(double* abc);

	/** @brief Returns the two intersections of a line with the sphere.
	 * 
	 * @param A[in] A Vertex denoting the start of line A->B. 
	 * @param B[in] A Vertex denoting the end of the line A->B. 
	 * 
	 * @return A pointer to a heap array of two Vertex objects which indicate
	 * the point of intersection with the line.
	 * 
	 * @note \ref intersectionCount must return a value of two or this function 
	 * will throw an error. 
	 */
	Vertex* getBothIntersections(const Vertex *A, const Vertex *B);
	
	/**@brief Returns the two intersections of a line with the sphere.
	 * 
	 * @param A[in] A Vertex denoting the start of line A->B. 
	 * @param B[in] A Vertex denoting the end of the line A->B. 
	 * 
	 * @return An array of length 3 containing the quadratic values of a,b,c
	 * which can be put into getQuadraticRoots to get the coordinates
	 * 		
	 * 		determine the intersection points
	 */
	double* getUValues(const Vertex *A, const Vertex *B);
	
	
	Vertex getMidPoint(const Vertex* A, const Vertex* B, double u) const;
	/* Used in Sphere::segTest
	 * Returns the point on line A->B that corresponds to the multiplication factor u,
	 * When used correctly, denotes the closest point on a line  to the shell of the sphere.
	 * See http://paulbourke.net/geometry/sphereline/ for more details.
	 */
	
	/** @brief Equality tester between two spheres
	 * 
	 * Inputs: Two vectors A and B which denote a line from A->B.
	 * 
	 * @param A[in] Vertex object pointer denoting the start point of a line.
	 * 		preferred if point A is inside the object, but not required.
	 * @param B[in] Vertex object pointer denoting an endpoint along the line
	 * @return The point farthest from A that intersects as a double vector, or
	 * 		will return NULL if no intersections are found.  
	 * 	@note Tests to see if there is
	 * 		an intersection of the line denoted  by a->b.  Will return the vector
	 * 		location of the farthest point from A that intersects the object.
	 */
	bool operator ==(const Sphere*);

	
};
class Triangle {
	/* Serves as the basis for a triangle in the system.  Consists of three
	 * vectors, and a number of basic computational functions to utilize
	 * the triangles effectively.  Does not inherit from Geom.
	 */
	//friend bool operator ==(const Triangle&); 
	bool operator ==(const Triangle&);
protected:
	//protected variables
	Vertex vertices[3]; //Stores the vertices of the triangle
public:
	Triangle();
	Triangle(const Vertex* corner, int size);
	Triangle(const Vertex* a, const Vertex* b, const Vertex* c);
	//Allows creation of triangles from 3 points. Maintains
	//the reference to the memory location, rather than creating
	//a new one.  This is a conservation choice, since Vertex
	//is currently immutable.
	~Triangle(); //Destructor
	bool sameSide(const Vertex* origin, const Vertex* test);

	const Vertex& getVertex(int index);
	//Returns a copy of a single vertex
	const Vertex** getVertices();
	//Returns the list of 3 vertices.

	bool sameSide(Vertex& testPt);
	//Returns whether two points are on the same side of the plane.
	Vertex* segTest(Vertex* a, Vertex* b);

	//Calculates the intersection points of the
};

class Tetrahedron: public Geom {
protected:
	Vertex *vert[4]; //4 Vertices of the tetrahedron that define it
	//Degeneracy is not checked until the determinant function is called
	// for the first time
	static double* buildMatrix(Vertex** pts, int count, bool filler = true);
	/**
	 * Builds a matrix for use with the determinant from the points of the
	 * tetrahedron.  
	 * Precondition:
	 * @param[in] pts An array of Vertex pointers with the number of items
	 * 	 contained in it equal to the value of count (ideally count=4).
	 * @param[in] filler a boolean trigger, which signals to add a column of 0's
	 * 	 in a 4th column to make a square matrix.
	 * \return The array of pts is unchanged.  A heap array of length
	 * 		count*(3+filler) (nominally 16 items) is returned.
	 */
	double determinant(const Vertex* matrix, int pos);
	bool checkTriangleIntersect(
		double* a, double* b, ///Points defining the line segment
		int* vertices, ///Len=3, indicies of vertices to test a triangle of.
		Vertex** results, /// Stores the intersecting points intersects
		int &resCount);
	/**
	 * Determines IF and WHERE the line segment A->B (where B is assumed
	 * to be an outside point) intersects the tetrahedron.  This is a wrapper
	 * for the triIntersect code provided by Brian Chen.
	 * Postcondition: If an intersection is found, resCount is incremented
	 * 	and the coordinates are stored in results, and true is returned.
	 */
public:
	Vertex* segTest(Vertex* a, Vertex* b);
	/*Returns the intersection of the tetrahedron with the line A->B
	 closest to point B with the tetrahedron.

	 */
	bool segTestIntersections(Vertex* A, Vertex* B,
			Vertex** results, //Array of pointers *res[2] to store results,
								//must start empty and uninitialized
			int& resCount /*Integer counter to indicate the number of intersection
						   contained in results**, must start = 0*/
			);
	/**
	 * Returns all intersections of the line A->B with the
	 *  tetrahedron if more than 1 exists.  For use with other classes that
	 *  require access to the intersections for their own sorting purposes
	 *  (namely Cup).
	 * Precondition:  A and B are Vertexes defining a line A->B
	 * 	"results" is an uninitialized array of pointers of size 2
	 * 	resCount is an integer counter for the number
	 * 	@return A boolean indicating whether any intersections were present.
	 *
	 */
	bool isInitialized();
	/*Returns true if 4 Vertex points exist within the tetrahedron, and all
	 * are unique.
	 */
	bool ptInside(const Vertex* point);
	//Returns true if the point is inside the tetrahedron
	Tetrahedron(); // Allows for empty set.
	Tetrahedron(Vertex** vertexi, int vcount);
	/* Constructor, Vertexi is the vertices of the tetrahedron, and the 
 	 * 
 	 */
	Tetrahedron(const Tetrahedron*); //Copy constructor
	Tetrahedron(Vertex*,Vertex*,Vertex*,Vertex*);
	~Tetrahedron();
	void addVertex(Vertex*);
	
};

/** Object Which has both a Tetrahedron "base" and a spherical probe which
 * defines the negative space of the cup.  The tetrahedron connects the 3
 * adjacent Sphere centers with that of a Sphere probe, which is positioned to be in contact
 * with all 3 Spheres and no other objects.  
 */
class Cup: public Geom{

protected:
	Sphere *probe;
	Tetrahedron *base;
public:
	Cup();
	Cup(Vertex** tetra_pts, int count, int sphere_base);
	Cup(Tetrahedron * tetra , Sphere * probe);
	/**
	 * Constructor transforming a Tetrahedron and a Sphere object into a
	 *  cup structure.
	 *  @param tetra A Tetrahedron object sharing at least one Vertex with
	 *  	the midpoint of probe
	 *  @param probe A Sphere that shares at least
	 */
	Cup(const Cup*);

	bool segTestIntersections(
			Vertex* a, Vertex* b, Vertex** results,
			char[],//resultSource[4] stores a letter indicating source
			int& baseCount, int& sphereCount,
			int& resultCount) const;
	Vertex* segTest(Vertex* a, Vertex* b);
	bool ptInside(const Vertex* pt);
	~Cup();


	
};


class Spindle{//public Geom {
	/*Two atoms + probe
	-Major triangle
	-Axis between two atoms
	-Probe sits tangental to atoms surfaces, as close to Axis as possible
	      o
	     /|
	    / |
	 (pO )|
	    \ |
	     \|
	      O


	Requires 3 points and 3 radii
	Plane of rotation not necessarily*/
protected:
	Sphere *probe;
	Sphere *ends[2];//Stores the endpoints of the spindle
public:
	Spindle();
	

	Spindle(Sphere* probe, Sphere* a, Sphere* b);

	/** Constructor for 3 spheres.
	 * \param centers The centers for the three Spheres defining the spindle
	 * \param radii The radii of the spheres
	 * \param probeIndex Indicates which of the indices represents the probe 
	 */
	Spindle(Vertex** centers, double* radii, int probeIndex);
	
	///Copy constructor
	Spindle(const Spindle* orig);
	//Copy Constructor
	~Spindle();

	/** An interpolation function for radius along the spindle
	 * \param A Center of a sphere endpoint
	 * \param B Center of a different sphere endpoint than A
	 * \param ptOnLine Vertex of a point you wish to find the radius at.
	 * 		Must be colinear with line AB.
	 * 	\return double value indicating the radius at ptOnLine.
	 */
	double interpolateRadius(Vertex* ptOnLine);
	bool ptInside(Vertex* pt);
	/* 
	 */
	
	
};

void check(bool isTrue, std::string errMsg, int exitcode=1, bool continu=false);

/** A function to find the closest point along a line with respect to the
 * input points.
 * \param	A[in]	Unordered endpoint of a line segment with B to test against
 * \param 	B[in]	Unordered endpoint of a line segment with A to test against
 * \param	pt[in]	Vertex denoting the point to find the closest point on the
 * 		segment.
 * \param	result[out]	Pointer to the result of the calculation. Should be
 * 		equal to NULL upon call, or the memory will be deleted during function.
 * \return	Boolean denoting whether the point has a normal to the line segment
 * 		If true is returned, the perpendicular point on the line is on or
 * 		between the endpoints of the line.
 * 		If false, the point lies beyond the endpoints, and one of the endpoints
 * 		is returned in result
 */
bool getClosestPointFromLine(Vertex* A, Vertex* B, Vertex* pt, Vertex*& result);

}//End Alpha namespace
#endif
