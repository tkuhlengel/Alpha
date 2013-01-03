/*
 * Delaunay.h
 *
 *  Created on: Aug 20, 2012
 *      Author: winnen
 */

#ifndef DELAUNAY_H_
#define DELAUNAY_H_
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>


namespace Alpha {


class ConnectGraph
{
public:
	int numPoints;
	int numTriangles;

	///mandatory memory allocation
	double * surfacePoints;
	double * surfaceNormals;		///these averaged normals are used for rendering
	int * triangles;
	double * triangleNormals;	///these normals are used for flat shading and interior checking
	double * centroid;

	///optional memory allocation
	int * highlights;
	int * edges; 				///edge highlights
	double * colors;			///optional colors, initialized only through separate request.  One vector for each pt.


	///////////////////////////////////////////////////////
	///Constructors / Destructors
	SurfaceObject();
	SurfaceObject(int numPts, double * ptsAndNorms, int numTris, int * inputTris, bool elimIntCavities);
	SurfaceObject(int numPts, double * pts, double * norms, int numTris, int * tris, double * triangleNorms);
	SurfaceObject(int* surfSet );
	virtual ~SurfaceObject();
	///////////////////////////////////////////////////////


	///////////////////////////////////////////////////////
	///returns the coords of the requested triangle
	///returns null if out of range.
	///notice that the triangle must be formatted in this way in the array (3 then 3 then 3)
	double * getTriangle(int num);
	///////////////////////////////////////////////////////



	///Computes the centroid based on the point positions
	double * getCentroid();




	/**Adds the geometry of another SurfaceObject to this SurfaceObject.
	ASSUMES THAT THIS SURFACEOBJECT DOES NOT COLLIDE WITH THE OTHER SURFACE OBJECT.
	DOES NOT CHECK FOR COLLISION, SO MUST BE DONE BEFORE THIS FUNCTION IS CALLED.
	*/
	void addObject(SurfaceObject * obj);

	/**Copies this surfaceObject into another surfaceObject
	If indices is non-NULL, returns a subset of the surface
	using only the points that are the indices of this surface.*/
	SurfaceObject * copy(set_t indices);

	///helper function tells you how int elements are identical between two sets.
	int countNumCommon(int* set1, int* set2);


	/**This function flips the normals of the surfaceObject backwards
	so that we can treat it as a negative volume
	*/
	void flipNormals();



	///////////////////////////////////////////////////////
	///add colors
	void addColors(double * c);
	///////////////////////////////////////////////////////


	///print out details
	void toString();
	void printSummary();
};


} /* namespace Alpha */
#endif /* DELAUNAY_H_ */
