/* Author: Trevor Kuhlengel
 * Geom.cpp
 * Purpose: This program was created to serve as the computational geometry
 * basis of the alpha shapes generator.  Its functions are integral in the 
 * translation of Alpha Shapes into meshes for use with other components of
 * this program.
 */

#include "Geom.h"

namespace Alpha {
bool Sphere::operator ==(const Sphere* that) {
	//Determines if the two spheres have equal centers and radii.
	if (this->center != NULL && that->center != NULL)
		return (*(this->center) == that->center)
				&& (this->radius == that->radius);
	else
		throw 49;
}

Sphere::Sphere() {
	this->radius = 0.0;
	this->center = new Vertex();

}

Sphere::Sphere(const Vertex *center, double radius) {
	this->radius = radius;
	this->center = new Vertex(center);
}
/*
 Sphere::Sphere(const Alpha::Atom *atm){
 this->center=new Vertex(atm->getCenter());
 this->radius=atm->Radius();
 }*/
Sphere::Sphere(const Sphere * that) {
	this->radius = that->radius;
	this->center = new Vertex(that->center);
}

Sphere::~Sphere() {
	delete this->center;
}

void Sphere::setRadius(double radius) {
	this->radius = radius;
}

void Sphere::setCenter(const Vertex *newcenter) {
	if (center != NULL) {
		delete center;
	}
	if (newcenter != NULL){
		center = new Vertex(*newcenter);
	}
}

double Sphere::getRadius() const {
	return radius;
}

Vertex* Sphere::getCenter() const {

	return center;
}


bool Sphere::ptInside(const Vertex *a) {
	//Uses difference squared, instead of using square root.
	double dist, rsq;
	dist = this->center->qDist(a);
	rsq = this->radius * this->radius;
	//std::cout<<"Distance = "<<dist<<"\nRadius Squared = "<<rsq<<std::endl;
	return ((dist < rsq) || doubleEqual(dist,rsq));
}

Vertex* Sphere::segTest(Vertex *A, Vertex *B) {
	Vertex *result;
	double a, b, c, sqr, u, u1;
	double *abc;

	bool a_inside, b_inside;
	//First make sure that both points are not inside the Sphere.  If they are
	// then there can be no intersection.
	a_inside = ptInside(A);
	b_inside = ptInside(B);
	Vertex *C;
	C = this->center;
	if (a_inside && b_inside)
		return NULL;
	//Source: http://paulbourke.net/geometry/sphereline/
	//Line Segment defined
	//Finding the intersection closest to B/P2
	//A=P1, B=P2
	//Find the closest point to the center of the sphere, defined by:
	//a*X^2+b*X+

	/*=================MOVED================
	 double a,b,c, *result=new double[3];
	 a = pow((B->getX() - A->getX()), 2) + pow((B->getY() - A->getY()), 2)
	 + pow((B->getZ() - A->getZ()), 2);
	 b = 2
	 * ((B->getX() - A->getX()) * (A->getX() - C->getX())
	 + (B->getY() - A->getY()) * (A->getY() - C->getY())
	 + (B->getZ() - A->getZ()) * (A->getZ() - C->getZ()));
	 c = pow(C->getX(), 2) + pow(C->getY(), 2) + pow(C->getZ(), 2)
	 + pow(A->getX(), 2) + pow(A->getY(), 2) + pow(A->getZ(), 2)
	 - 2
	 * (C->getX() * A->getX() + C->getY() * A->getY()
	 + C->getZ() * A->getZ())
	 - this->radius * this->radius;

	 return result;

	 ===================MOVED==================*/
	abc = getUValues(A, B);
	a = abc[0];
	b = abc[1];
	c = abc[2];

//The sqrt component from the quadratic solution of the form a*U^2 + b*U + c;
	sqr = b * b - 4 * a * c;

	if (0 > sqr) { //Imaginary intersection - Then the line doesn't interesect
				   //the sphere at all
		result = NULL;
	} else if (0 == sqr) { // Then the line is tangent to the sphere
		//And intersects the sphere at u=-b/2a
		u = -1 * b / 2 * a;
		result = new Vertex(getMidPoint(A, B, u));
	} else if (0 < sqr) { //Then the line intersects in 2 places
		// This could be replaced with getBothIntersections(A,B)
		u = (-b + sqrt(sqr)) / (2 * a);
		u1 = (-b - sqrt(sqr)) / (2 * a);

		//Calculate the vertices of intersection
		Vertex t, t1;
		t = getMidPoint(A, B, u);
		t1 = getMidPoint(A, B, u1);

		//Distances from the points to B, since only want to return one of them
		double t_d = t.qDist(B);
		double t1_d = t1.qDist(B);

		//Find the closer one to B
		if (doubleEqual(t_d,t1_d) || t_d < t1_d) {
			result = new Vertex(&t);
		}
		else{
			result = new Vertex(&t1);
		}
	}
	return result;
}

int Sphere::intersectionCount(const Vertex *A, const Vertex *B) {

	double sqr, *u, *d = getUValues(A, B);
	int number=3;

	//sqr = b * b - 4 * a * c;
	sqr = d[1] * d[1] - 4 * d[0] * d[2];

	if (0 > sqr) { //Imaginary numbers, and is a miss
		return 0;
	}
	else if (0 == sqr){ //Exactly 1 point tangent to sphere
		return 1;
	}
	else { // 2 points along the infinite line. This should always be the case during a single intersection
		u=getQuadraticRoots(d);
		//Line segment is tangential to the sphere, in which case both values
		//of u will be the same and between 0 and 1
		if (doubleEqual(u[0],u[1]) && (u[0]>=0 && u[0]<=1)){
			return 1;
		}
		//Line segment doesn't intersect and on outside of sphere
		else if((u[0]<0 && u[1]<0)||(u[0]>1 && u[1]>1)){
			return 0;
		}
		//Line segment doesn't intersect and is inside sphere,
		//in which case one value of u will be negative and the other greater than 1.
		else if((u[0]<0 && u[1]>1)||(u[0]>1 && u[1]<0)){
			return 0;
		}
		//Line segment intersects at two points, in which case both values of u will be between 0 and 1.
		else if((u[0]>0 && u[0]<1)&&(u[1]>0 && u[1]<1)){
			return 2;
		}
		//Line segment intersects at one point, in which case one value of u will be between 0 and 1 and the other not.
		else if((u[0]>0 && u[0]<1) xor (u[1]>0 && u[0]<1)){
			return 1;
		}
		delete u;
	}
	return number;
}
double* Sphere::getQuadraticRoots(double *abc){
	//Must be deleted
	double *g=abc, *result = new double[2];
	//-b +/- sqrt(b^2 - 4ac)/2a
	result[0]=(-g[1] + sqrt(g[1]*g[1]-4*g[0]*g[2]))/(2*g[0]);
	result[1]=(-g[1] - sqrt(g[1]*g[1]-4*g[0]*g[2]))/(2*g[0]);
	return result;
}
Vertex* Sphere::getBothIntersections(const Vertex *A, const Vertex *B) {
	//Returns a list of 2 vertexes.
	if (2 != intersectionCount(A, B)) {
		return NULL;
	}
	double sqr, a, b, c, //Quadratic formula components
			u, u1, //u values for point determination
			*d = getUValues(A, B);
	a = d[0];
	b = d[1];
	c = d[2]; //mapping them for readability
	sqr = b * b - 4 * a * c;

	//Get U values for the two points
	u = (-b + sqrt(sqr)) / (2 * a);
	u1 = (-b - sqrt(sqr)) / (2 * a);

	//Calculate vertices from U values, gives intersection points
	Vertex t, t1;
	t = getMidPoint(A, B, u);
	t1 = getMidPoint(A, B, u1);
	Vertex *result;
	result = new Alpha::Vertex[2];
	result[0] = t;
	result[1] = t1;
	return result;

}
double* Sphere::getUValues(const Vertex *A, const Vertex *B) { //C defines the center of the Sphere
	/*Finding the intersection closest to B/P2
	 *A=P1, B=P2
	 *Find the closest point to the center of the sphere, point C, defined by:
	 * a*X^2+b*X+c=0
	 * Returns the a, b and c values for the calculations in Quadratic Square roots
	 *
	 */
	Vertex *C = this->center;
	double a, b, c, *result = new double[3];
	a = pow((B->getX() - A->getX()), 2) + pow((B->getY() - A->getY()), 2)
			+ pow((B->getZ() - A->getZ()), 2);
	b = 2* ((B->getX() - A->getX()) * (A->getX() - C->getX())
					+ (B->getY() - A->getY()) * (A->getY() - C->getY())
					+ (B->getZ() - A->getZ()) * (A->getZ() - C->getZ()));
	c = pow(C->getX(), 2) + pow(C->getY(), 2) + pow(C->getZ(), 2)
			+ pow(A->getX(), 2) + pow(A->getY(), 2) + pow(A->getZ(), 2)
			- 2 * (C->getX() * A->getX() + C->getY() * A->getY()
							+ C->getZ() * A->getZ())
			- this->radius * this->radius;
	result[0] = a;
	result[1] = b;
	result[2] = c;
	return result;
}

Vertex Sphere::getMidPoint(const Vertex * A, const Vertex * B, double u) const {
	double x, y, z;

	x = A->getX() + u * (B->getX() - A->getX());
	y = A->getY() + u * (B->getY() - A->getY());
	z = A->getZ() + u * (B->getZ() - A->getZ());
	Vertex result = Vertex(x, y, z);
	return result;

}

Tetrahedron::Tetrahedron() {
	for (int i = 0; i < 4; i++) {
		vert[i] = new Vertex();
	}
}
Tetrahedron::Tetrahedron(Vertex** vertexi, int vcount) {
	if (vcount > 4)
		throw 47;
	else {
		for (int i = 0; i < vcount; i++) {
			vert[i] = new Vertex(vertexi[i]);
		}
	}
}
Tetrahedron::Tetrahedron(const Tetrahedron * t) {
	for (int i = 0; i < 4; i++) {
		vert[i] = new Vertex(t->vert[i]);
	}
}
Tetrahedron::Tetrahedron(Vertex * A, Vertex * B, Vertex * C, Vertex * D) {
	vert[0] = new Vertex(A);
	vert[1] = new Vertex(B);
	vert[2] = new Vertex(C);
	vert[3] = new Vertex(D);
}
Tetrahedron::~Tetrahedron() {
	for (int i = 0; i < 4; i++) {

		delete vert[i];
	}
}
bool Tetrahedron::isInitialized() {
	for (int i = 0; i < 4; i++) {
		if (vert[i] == NULL) {
			return false;
		}
		for (int i = 0; i < 4; i++) {
			for (int j = i + 1; j < 4; j++) {
				if (vert[i] == vert[j]) {
					return false;
				}
			}
		}
	}
	return true;

}
bool Tetrahedron::ptInside(const Vertex *t) {
	/*Source: http://steve.hollasch.net/cgindex/geometry/ptintet.html
	 * Determines if a Vertex is inside or outside the indicated tetrahedron.
	 * Precondition: Tetrahedron must be instantiated with 4 Vertex objects.
	 */
	double result[5];
	Vertex *test = new Vertex(t); //Something not const to pass around

	bool x,  //Test for whether the determinants have the same sign.
			b_result = true; //Boolean Output from the function
	for (int i = 0; i < 5; i++) {
		result[i] = this->determinant(test, i); //Check the 5 determinants to see if they have the same sign
		if (0 == i) {
			if (result[i] == 0.0) {
				throw "Error: Tetrahedron is degenerate and all points are coplanar";
			}
			x = (result[i] < 0);

		} else if ((result[i] < 0.0) != x) {
			if (result[i] != 0.0) { //Coplanarity check. If coplanar, want to consider it to be inside.
				b_result = false;
			}
		}
	}
	delete test;
	return b_result;
}

double Tetrahedron::determinant(const Vertex *t, int pos) {
	Vertex *v0 = vert[0], *v1 = vert[1], *v2 = vert[2], *v3 = vert[3], *test =
			new Vertex(t);
	switch (pos) {
	case 0:  //Degeneracy check as well as initial check
		break;
	case 1:
		v0 = test;
		break;
	case 2:
		v1 = test;
		break;
	case 3:
		v2 = test;
		break;
	case 4:
		v3 = test;
		break;
	default:
		break;
	}
	Vertex *v[4] = { v0, v1, v2, v3 };
	double *matrix = Tetrahedron::buildMatrix(v, 4, true);
	double result = determinant4x4(matrix);
	delete matrix;
	delete test;
	return result;
	/*
	 return (*v0)[0]*(*v1)[1]*(*v2)[2] - (*v0)[0]*(*v1)[1]*(*v3)[2] -
	 (*v0)[0]*(*v2)[1]*(*v1)[2] + (*v0)[0]*(*v2)[1]*(*v3)[2] +
	 (*v0)[0]*(*v3)[1]*(*v1)[2] - (*v0)[0]*(*v3)[1]*(*v2)[2] -
	 (*v1)[0]*(*v0)[1]*(*v2)[2] + (*v1)[0]*(*v0)[1]*(*v3)[2] +
	 (*v1)[0]*(*v2)[1]*(*v0)[2]  - (*v1)[0]*(*v2)[1]*(*v3)[2] -
	 (*v1)[0]*(*v3)[1]*(*v0)[2] + (*v1)[0]*(*v3)[1]*(*v2)[2] +
	 (*v2)[0]*(*v0)[1]*(*v1)[2] - (*v2)[0]*(*v0)[1]*(*v3)[2] -
	 (*v2)[0]*(*v1)[1]*(*v0)[2] + (*v2)[0]*(*v1)[1]*(*v3)[2] +
	 (*v2)[0]*(*v3)[1]*(*v0)[2] - (*v2)[0]*(*v3)[1]*(*v1)[2] -
	 (*v3)[0]*(*v0)[1]*(*v1)[2] + (*v3)[0]*(*v0)[1]*(*v2)[2] +
	 (*v3)[0]*(*v1)[1]*(*v0)[2] - (*v3)[0]*(*v1)[1]*(*v2)[2] -
	 (*v3)[0]*(*v2)[1]*(*v0)[2] + (*v3)[0]*(*v2)[1]*(*v1)[2];
	 */
}
double* Tetrahedron::buildMatrix(Vertex **pts, int count, bool filler) {
	int col = 3;

	if (filler)
		col++;
	double * matrix = new double[col * count];
	for (int i = 0; i < count; i++) { //Iterating through the points
		for (int j = 0; j < col; j++) { //Iterating through the axes, or columns of the matrix
			if (filler && 3 == j) {
				matrix[(col * i) + j] = 1; //Makes the 4th column of ones.
			} else {
				matrix[(col * i) + j] = (*(pts[i]))[j];
			}
		}
	}
	return matrix;
}
bool Tetrahedron::segTestIntersections(Vertex * a, Vertex * b, Vertex **results, //Array of pointers *res[4] to store results,
																				 //must start empty and uninitialized
		int &resCount /*Integer counter to indicate the number of intersection
		 contained in results**, must start = 0*/
		) {
	/* Tests for the intersection of a line segment A->B with the Tetrahedron.
	 * If more than one intersection is found, it returns the intersection closest
	 * to point B.
	 *
	 */
	/* Triangles of the Tetrahedron
	 0,1,2
	 1,2,3
	 0,1,3
	 0,2,3
	 */
	int pts[4][3] = { //Triangles of tetrahedron by index
			{ 0, 1, 2 }, { 1, 2, 3 }, { 0, 1, 3 }, { 0, 2, 3 } };
	double A[3], B[3];
	A[0] = a->getX();
	A[1] = a->getY();
	A[2] = a->getZ();
	B[0] = b->getX();
	B[1] = b->getY();
	B[2] = b->getZ();
	resCount = 0;
	bool intersects = false;

	for (int i = 0; i < 4; i++) {
		if (this->checkTriangleIntersect(A, B, pts[i], results, resCount)) {
			intersects = true;
		}
	}
	return intersects;

}
Vertex* Tetrahedron::segTest(Vertex *a, Vertex *b) {
	Vertex *results[2], *return_ptr;
	int resCount = 0;
	if (segTestIntersections(a, b, results, resCount)) {
		if (1 < resCount) { //More than 1 intersection
			double dist, minVal = results[0]->qDist(b); //Arbitrary large value
			int minIndex = 0;
			for (int i = 0; i < resCount; i++) { //Find the one closest to b
				dist = results[i]->qDist(b);
				if (dist < minVal) {
					minVal = dist;
					minIndex = i;
				}
			}
			return_ptr = results[minIndex];
		} else {
			return_ptr = results[0];
		}
	} else {
		return_ptr = NULL;
	}
	return return_ptr;

}
bool Tetrahedron::checkTriangleIntersect(double *a, double *b, //Points defining the line segment
		int *v, //List of integers, indexing the vertices[] to test (e.g. this->vert[v[0]])
		Vertex **results, // Stores the intersections on the heap.
		int &resCount) //Number of intersections with the tetrahedron
		{
	double *m;
	/*NOTE: This is where the main requirement of Brian Chen's code is,
	 and the mathfunc.h include*/
	m = intersect_Seg_Triangle(a, b, this->vert[v[0]]->getPos(),
			this->vert[v[1]]->getPos(), this->vert[v[2]]->getPos());
	if (m != NULL) {
		results[resCount++] = new Vertex(m);
		return true;
	} else
		return false;
}

Cup::Cup() {
	this->probe = new Sphere();
	this->base = new Tetrahedron();
}

Cup::Cup(Vertex **tetra_pts, int count, int sphere_base) {
	if (sphere_base > (count - 1)) {
		throw 30;  //the sphere must be centered in the
	}
	this->probe = new Sphere(tetra_pts[sphere_base], 1.0); //PROBE_SIZE);
	this->base = new Tetrahedron(tetra_pts, count);
}
Cup::Cup(Tetrahedron * tetra, Sphere * probe) {
	this->probe = new Sphere(probe);
	this->base = new Tetrahedron(tetra);
}
Cup::Cup(const Cup* orig) {
	this->probe = new Alpha::Sphere(orig->probe);
	this->base = new Alpha::Tetrahedron(orig->base);
}
Cup::~Cup() {
	delete this->probe;
	delete this->base;
}

bool Cup::ptInside(const Vertex * pt) {
	return ((!this->probe->ptInside(pt)) && this->base->ptInside(pt));
}

bool Cup::segTestIntersections(Vertex * a, Vertex * b,
		Vertex **results, //Must be length 4 *results[4], stores output intersections
		char resultSource[], //Stores 'b' for tetrahedron, 's' for sphere
					//Aligns to Vertex in results
		int &baseCount, int &sphereCount, int &resultCount
		) const {
	Vertex *baseResults[2], *sphereResults; ///*results[4];
	baseCount = 0; sphereCount=0; resultCount = 0;

	/*Indicator to track where vertexes came from for case evaluations
	 * It needs to start blank*/
	for (int i=0; i<4; i++){
		resultSource[i] = ' ';
	}

	//Tetrahedron intersections
	if (this->base->segTestIntersections(a, b, baseResults, baseCount)) {

		for (int i = 0; i < baseCount; i++) {
			resultSource[resultCount] = 'b';
			results[resultCount++] = baseResults[i];
		}
	}
	/*If it never intersects the base, there is no intersection with the cup
	 structure*/

	//Sphere intersections
	sphereCount = this->probe->intersectionCount(a, b);
	if (sphereCount > 1) { //Two intersection case
		sphereResults = this->probe->getBothIntersections(a, b);
		for (int i = 0; i < sphereCount; i++) {
			resultSource[resultCount] = 's';
			results[resultCount++] = &sphereResults[i];
		}
	//One intersection case
	} else if (sphereCount == 1) {
		resultSource[resultCount] = 's';
		results[resultCount++] = this->probe->segTest(a, b);

	}
	if (resultCount>0)
		return true;
	else
		return false;
}
Vertex* Cup::segTest(Vertex *a, Vertex *b){
	/* If there is an error in logic, it happens in this function
	 * if there is an error in points, it happens in segTestIntersections()
	 */

	Vertex *results[4];
	char resultSource[4];
	int baseCount, sphereCount, resultCount;

	//Confirm that the inputs are valid.
	check((a!=NULL && b!=NULL),
		"Error, one of the inputs to Cup::segTest is null", 40, false);
	/*Get the intersection points in results, and all the other data
	 * Necessary to evaluate the
	 */

	bool intersects=segTestIntersections(a,b,
			results,
			resultSource,
			baseCount,
			sphereCount,
			resultCount);
	if (!intersects){
		return NULL;
	}
	double *distance=new double[resultCount];
	for (int i=0; i<resultCount; i++){
		distance[i]=b->qDist(results[i]);
	}

	//One intersection cases
	if (resultCount==1){

		bool insideTet; //Is the point inside the tetrahedron
		insideTet=this->base->ptInside(results[0]);

		//Sphere intersection: Intersection pt is inside the Tetrahedron
		if (insideTet && resultSource[0]=='s'){
			return results[0];
		}
		//Point on the tetrahedron, not contained in the sphere
		else if(insideTet && !this->probe->ptInside(results[0])){
			return results[0];
		}
		//TODO Finish
	}
	//Further evaluations need to be sorted
	//Bubble Sort the arrays sorting data in place
	else if (resultCount>1){
		bool swapped = true;
		int j = 0;
		char tempc;
		Vertex *tempv;
		double tempd;
		while (swapped) {
			swapped = false;
			j++;
			for (int i = 0; i < resultCount - j; i++) {
				if (distance[i] < distance[i + 1]) {
					//Temps
					tempd = distance[i];
					tempv = results[i];
					tempc = resultSource[i];

					//First Swaps
					distance[i] = distance[i + 1];
					results[i]=results[i+1];
					resultSource[i]=resultSource[i+1];

					//Final Swap
					distance[i + 1] = tempd;
					results[i+1] = tempv;
					resultSource[i+1]=tempc;
					swapped = true;
				}
			}
		}
	}


	//For simplicity, I'm renaming the result array for easy programming
	char *rS=resultSource;

	//Case evaluation working from B to A
	/*
	if (this->probe->ptInside(b)){
		if (rS[0]=='s'){
			if (rS[1]=='b'){

			}
			//Intersection is with free sphere
			else if (rS[1]==' '){
				return NULL;

			}

		}

		else if (rS[0]=='b'){

		}
		else if (rS[0]==' '){

		}
	}*/

	//TODO Tangental hits sbb bsb bbs
	//Tetrahedron intersection is closest to B
	//Tree levels are denoted by comments, easier to follow
	//1
	if(rS[0]=='b'){
		//2
		if (rS[1]=='s'){
			//3 b s Tetrahedron Outside boundary case
			if (rS[2]=='s'){
				return results[0];
				}

			//3 b s *
			else if (rS[2]=='b'){

				//4 Tetrahedron boundary case
				//b s b *
				if (rS[3]=='s'){
					return results[0];
				}
				//4 b s b *
				else if (rS[3]==' '){
				//Tangental Hit?
					//This is only necessary if ptInside()
					//is not strictly necessary for execution.
		//TODO A or B must end in Probe sphere, or there is a Tangent
					if (this->probe->ptInside(a)){
						return NULL;
					}
					else{
						check(this->ptInside(a)||this->ptInside(b),
					"Error: At least one point must be inside the cup surface",
							42, false);
					}

				}
			}
			//3 b s *
			else if (rS[2]==' '){
				if (probe->ptInside(b)){
					return results[1];
				}

			}
		}
		//2: b *
		else if (rS[1]=='b'){
			//3 b b *
			if (rS[2]=='s'){
				//4 b b s *
				if (rS[3]=='s'){
					return results[0];
				}
				//4 b b s *
				else if (rS[3]==' '){
					//Make sure that b isn't inside the sphere
					if (probe->ptInside(b)){
						return NULL;
					}
					else{
						return results[0];
					}
				}

			}
			//3 b b *
			else if (rS[2]==' '){
				return results[0];
				/*This case should be impossible
				 * Due to the requirement that pieces have an endpoint inside the
				 * structure, and the grid size optimization
				 *
				Check to make sure that the endpoints are not inside probe
				 *
				if (probe->ptInside(a) || probe->ptInside(b)){
					//Remove COUT after testing
					std::cout<<"Unexpected case occurred in Cup::segTest "
							<<"'b b _' :"
							<<"\n\tat least one endpoint was inside the probe"
							<<" sphere"<<std::endl;
					return NULL;
				}
				else{
					return results[0];
				}
				*/
			}
		}
		//2 b
		else if (rS[1]==' '){
			std::cout<<"Error: Case b _ should never happen. Single intersection after sorting."
				<<std::endl;
			if (base->ptInside(results[0])){
				return results[0];
			}
			else{
				return NULL;
			}
		}

		//return results[0];
	}
	//1 Sphere intersection closest to B cases
	else if(rS[0]=='s'){
		//2 s *
		if (rS[1]=='s'){
			//3 s s *
			if (rS[2]=='b'){
				//4 s s b *
				if (rS[3]=='b'){
					return results[2];
				}
				//4 s s b *
				else if (rS[3]==' '){
					//Make sure that a or b is inside the tetrahedron
					if(base->ptInside(a)){
						return results[2];
					}
					else if (base->ptInside(b)){
						return results[0];
					}
					else{
					check(false, "Error in Cup.segTest: At least one point\
	(a or b) must be in the tetrahedral base\n\
	 in case ssb.", 41, false );
					}
				}
			}
			//3 s s *
			else if (rS[2]==' '){
				//2 internal points on grid
				if (base->ptInside(a)){
					if (base->ptInside(b)){
						return results[0];
					}
				}
				//default else due to return
				//Miss case
				return NULL;
			}
		}
		//2 s *
		else if (rS[1]=='b'){
			//3 s b *
			if (rS[2]=='s'){
				//4 s b s *
				if (rS[3]=='b'){
					//Sphere intersection inside tetrahedron case
					return results[2];
				}
				//4 s b s *
				else if (rS[3]==' '){
					if(base->ptInside(b)){
						return results[0];
					}
					else if(base->ptInside(a)){
						return results[2];
					}

				}
			}
			//3 s b *
			else if (rS[2]=='b'){
				//4 s b b *
				if (rS[3]=='s'){
					//Miss Case
					return NULL;
				}
				//4 s b b *
				else if (rS[3]==' '){
					if(probe->ptInside(a)){
						return NULL;
					}
					else{
						return results[1];
					}

				}
			}
			//3 s b *
			else if (rS[2]==' '){
				//TODO Need to solve this.
				//This is a confusing case, because there are only 2 points
				//to work with in order to figure out. So I need to use inside
				//checks in order to figure out where everything is.
				bool bInTet=base->ptInside(b),//is pt B in the Tetrahedron base?
					 bInSph=probe->ptInside(b),

					 aInTet=base->ptInside(a),
					 aInSph=probe->ptInside(a);//Is pt A in the probe sphere?
					 
				if(bInTet){ //Both subcases involving B are identical
					return results[0];

				}
				else if(aInTet){
					if (bInSph){
						return results[0];
					}
					else if (aInSph){

					}

				}

			}
		}
		/*
		//2 s *
		else if (rS[1]==' '){

		}
		*/

	}

/*
if (rS[1]=='s'){

}

else if (rS[1]=='b'){

}
else if (rS[1]==' '){

}
 */

	throw 59;
}
Spindle::Spindle(){
	
}




Spindle::Spindle(Sphere *probe, Sphere *a, Sphere *b){
	check(probe!=NULL && a!=NULL && b!=NULL, 
		"Spindle->Spindle(Sphere*, Sphere*, Sphere*) Pointers must not be null",
		61, false);
	this->probe=new Sphere(probe);
	this->ends[0]=new Sphere(a);
	this->ends[1]=new Sphere(b);
}

Spindle::Spindle(Vertex **centers, double *radii, int probeIndex){
	bool first=true;
	for (int i=0; i<3; i++){
		if (i==probeIndex){
			this->probe=new Sphere(centers[i], radii[i]);
		}
		else if(first){
			this->ends[0]=new Sphere(centers[i], radii[i]);
			first = false;
		}
		else{
			this->ends[1]=new Sphere(centers[i], radii[i]);
		}
		
	}
}

Spindle::~Spindle(){
	delete probe;
	delete ends[0];
	delete ends[1];
}
double Spindle::interpolateRadius(Vertex * ptOnLine){
	Vertex * result=ptOnLine;
	Vertex *A=ends[0]->getCenter(), *B=ends[1]->getCenter();
	double Ar=ends[0]->getRadius(), Br=ends[1]->getRadius();
	double apDist, abDist, ap_ab_ratio, radiusRatio, radiusAtPt;

	//Get distances and ratio
	apDist= A->distance(result);
	abDist=A->distance(B);
	ap_ab_ratio=apDist/abDist;

	check(ap_ab_ratio<1.0 || doubleEqual(ap_ab_ratio, 1.000),
		"Unexpected Error in Spindle::ptInside: Distance from A to pt must be less than or equal to distance AB",
		64, false);
	//Interpolate the radius of the truncated cylinder at the result point
	radiusRatio=(Br-Ar)/abDist;
	radiusAtPt=Ar+radiusRatio*ap_ab_ratio;
	return radiusAtPt;
}

bool Spindle::ptInside(Vertex *pt){
	//TODO Spindle::ptInside
	/*=======================================================================
	  Step 1: Figure out if the point is inside the Cylinder/Truncated Cone
	  	  1a: Using the truncated cone model, find the nearest point on the
	  	  	  line AB to the point p (the normal to AB intersecting with P).
	  	  1b: Interpolate the radius of the cone at that point by using the
	  	  	  z difference between A and B, and the radii Ar and Br,
	  	  	  lets call it Zr.
	  	  1c: If the distance d(p,AB) is less than Zr, then it is inside the
	  	      cone and it can be used.
	 ========================================================================*/

	//Making lazy assumption that the circle centers are the endpoint of the axis
	Vertex *A=ends[0]->getCenter(), *B=ends[1]->getCenter(), *ptMidline=NULL;
	double  radiusAtPt,	linePtDist, probeRadius=probe->getRadius();
	
	bool insideEnds=getClosestPointFromLine(A,B,pt,ptMidline);
	if(!insideEnds){
		return false;
	}
	//Implicit else - Find the interpolation of the spheres at the result
	radiusAtPt=this->interpolateRadius(ptMidline);
	linePtDist=ptMidline->distance(pt);

	//Eliminate ones that are greater than the cone radius
	if(linePtDist>radiusAtPt && !doubleEqual(linePtDist, radiusAtPt)){
		return false;
	}
	/*=======================================================================
	 * Step 2: If it is inside the cylinder, we have a few things to do before
		we know if it is in fact inside the spindle
		2a: Define a Plane F which is coplanar with both endpoints and the
			test point pt.
		2b: Set a sphere H coplanar with the Torus and Plane F at the correct
			position with respect to atoms A and B. Sphere H represents an
			instance of the probe sphere along the plane of the intended line.
		2c: Test if the point is inside H. If it is, then we know it is not
			inside the spindle.
	=========================================================================*/
	Vertex *C=probe->getCenter(),
			*probeMid=NULL, //The point on AB axis closest to probe
	*directionToPt=pt->minus(ptMidline);//Vector denoting direction from
		//result to test point

	bool probeInsideEnds=getClosestPointFromLine(A,B,C,probeMid);
	check(probeInsideEnds,
		"Unexpected Error in Spindle::ptInside: Probe is expected to have "
		"a closest point within the boundaries of the endpoints"
		,65,true);
	double x,y,z,
	scaler,//Ratio of distance to probe over distance to test point
	probeDist=probeMid->distance(C);
	scaler=probeDist/linePtDist;

	//Scale the points in the direction of test point, starting from the probe
	// test point, so that they share a plane. Purpose is to create an instance
	// circle (sphere) of the and reduce this to a pseudo-2D problem
	double *dDirectionToPt=directionToPt->getPos(),
			*dProbeMid=probeMid->getPos();
	x=dProbeMid[0]+scaler*dDirectionToPt[0];
	y=dProbeMid[1]+scaler*dDirectionToPt[1];
	z=dProbeMid[2]+scaler*dDirectionToPt[2];
	Vertex probeInstanceCenter //Center point of the probe sharing the plane
		=Vertex(x,y,z);
	Sphere probeInstance=Sphere(&probeInstanceCenter, probeRadius);
	bool ptInsideTorus=probeInstance.ptInside(pt);
	if (ptInsideTorus){
		delete directionToPt;
		return false;
	}

	delete directionToPt;
	return true;
}
Vertex* pointAlongLine(Vertex *A, Vertex *B, double scaler){
	/*Returns a point along a line which scales the relative location of B,
	With A remaining fixed in space. As an example, if the scaler is 2.0,
	the function will return the location of C shown below.
	Old:   A--------->B
	Return A---------B---------->C
	\param scaler A ratio of the distance from A to the new point divided
			by the original distance to B
	*/
	Vertex *directionToB=B->minus(A);

	//Double copies. For some reason, C++ doesn't like using dereferencing
	// AND array calls
	double *dDirectionToB=directionToB->getPos(),
			*dA=A->getPos()//Vector location of A
			,x,y,z;
	x=dA[0]+scaler*dDirectionToB[0];
	y=dA[1]+scaler*dDirectionToB[1];
	z=dA[2]+scaler*dDirectionToB[2];
	return new Vertex(x,y,z);
}
bool getClosestPointFromLine(Vertex *A, Vertex *B, Vertex *pt, Vertex *&result){
	//First, find the equation of the line in Mx+B form.
	Vertex *M=B->minus(A)
		  ,*N=pt->minus(A); //Gets a heap vector denoting the vector from A->B

	
	check(!(A->equals(B)), "Error in getClosestPointFromLine: A and B must not be coincident",
			62, false);


	//std::cout<<"Point tested: "<<pt<<std::endl;


	double u, x,y,z;
	u=DOT((*N),(*M))/DOT((*M),(*M));
	delete N; delete M;
	double *a=A->getPos(), *b=B->getPos();
	x=a[0]+u*(b[0]-a[0]);
	y=a[1]+u*(b[1]-a[1]);
	z=a[2]+u*(b[2]-a[2]);

	result=new Vertex(x, y, z);


	//std::cout<<"Translated Point: "<<result<<std::endl;


	double dist=A->distance(B);
	double apDist=A->distance(result), bpDist=B->distance(result);
	//A + B = C test.  If this fits, then it lies on or between the two points.
	if (doubleEqual((apDist+bpDist),dist)){
		return true;
	}
	else{//It is past one of the endpoints
		//If closer to A, then it is past A, and we return A
		if (apDist < bpDist){
			delete result;
			result= new Vertex(A);
			return false;
		}
		//If it is closer to B, then it is past B, and we return B
		else if(bpDist < apDist){
			delete result;
			result=new Vertex(B);
			return false;
		}
	}
	
	return true;
}
void check(bool isTrue, std::string errMsg, int exitcode, bool continu){
	if (!isTrue){
		std::cout<<errMsg<<std::endl;
		if (!continu){
			throw exitcode;
		}
	}
}

} //Alpha Namespace
