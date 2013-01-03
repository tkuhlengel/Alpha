#include <iostream>
#include <cstdlib>
#include "rundelaunay.h"
#include "Geom.h"
#include "vertex.h"
using namespace std;

bool vertexTest(){
	Alpha::Vertex *A=new Alpha::Vertex(1.0f,0.0f,0.0f)
		,*B=new Alpha::Vertex(0.5f, 0.5f, 0.5f)
		,*C=new Alpha::Vertex(1.0f, 2.5f, -1.5f)
		,*D//=new Alpha::Vertex(0.5f, 0.5f, 0.5f)
		,*pt=new Alpha::Vertex(0.5f, 0.5f, 0.5f)
	;
	cout<<"A = "<<A<<"\nB = "<<B<<endl
		<<"Distance from A to B = "
		<<A->distance(B)<<endl
		<<"Vector from A to B is D. \nD = ";
	D= B->minus(A);
	cout<<D<<endl;
	cout<<"A plus D = "<<A->plus(D)<<endl;

	return true;
	delete A; delete B;delete C; delete D; delete pt;
}
Alpha::Tetrahedron* tetTest() {
	/*Creates a simple test to see if functions are working as intended
	 * for the tetrahedron.
	 */
	int nsteps=20;
	Alpha::Vertex *A, *B, *C, *D, *pt, *ptA, *ptB, *pts;

	Alpha::Tetrahedron *test;
	A = new Alpha::Vertex(0.0, 0.0, 0.0);
	B = new Alpha::Vertex(1.0, 0.0, 0.0);
	cout << A << endl;
	cout << "Distance from {0,0,0} to {1,0,0} = " << A->distance(B) << endl;
	C = new Alpha::Vertex(0.0, 1.0, 0.0);
	cout << "Distance from {1,0,0} to {0,1,0} = " << B->distance(C) << endl;
	D = new Alpha::Vertex(0.0, 0.0, 1.0);
	test = new Alpha::Tetrahedron(A, B, C, D);

	//Boundary Case
	pt = new Alpha::Vertex(1.0, 0.0, 0.0);

	//Tests for segtest
	ptA = new Alpha::Vertex(0.1, 0.1, 0.1);
	ptB = new Alpha::Vertex(1.0, 1.0, 1.0);


	// *pts=new Alpha::Vertex[10];
	pts=new Alpha::Vertex[20];
	double x=2.0/nsteps,a=1.0, b=-1.0;
	for (int i = 0; i < nsteps; i++) {
		pts[i] = Alpha::Vertex(b + a* x * i,
								b + a * x * i,
								b + a * x * i);
	}
	for (int i = 0; i < nsteps; i++) {
			cout << "Testing if the point pt=" << &pts[i]
					<< "is inside the z=1-x-y tetrahedron." << endl
					<< test->ptInside(&pts[i]) << endl;
			//delete pts[i];
		}



	//delete test;
	delete A;
	delete B;
	delete C;
	delete D;
	delete pt;
	delete ptA;
	delete ptB;
	delete[] pts;

	return test;
}
Alpha::Sphere* sphereTest() {
	cout<<"Testing Sphere Functions... "<<endl;
	Alpha::Sphere *sph;
	Alpha::Vertex *B, *ptA, *ptB, *ptC, *ptD, *delMePlz;
	B = new Alpha::Vertex(1.0, 0.0, 0.0);
	sph=new Alpha::Sphere(B, 0.9);
	ptA = new Alpha::Vertex(0.1,0.0, 0.0);
	ptB = new Alpha::Vertex(1.0, 1.0, 1.0);
	ptC = new Alpha::Vertex(0.0,0.0,0.0);
	ptD = new Alpha::Vertex(0.5,0.5,0.5);

	cout<<"\tTesting known inside point... "<<sph->ptInside(B)
		<<"\n\tTesting known boundary point... "<<sph->ptInside(ptA)
		<<"\n\tTesting known outside point... "<<sph->ptInside(ptC)
		<<"\n\tTesting unknown point... "<<sph->ptInside(ptD)<<endl;



	cout<<"\tTesting known inside-outside line... ";
	delMePlz=sph->segTest(B, ptC);
	cout<<"\n\t\tExpected: "<<ptA<<"\n\t\tResult:   "<<delMePlz<<endl;
	cout<<"\t\t Number of Intersections: "<<sph->intersectionCount(B, ptC)<<endl;
	cout<<"\n\t\t Coordinates of intersection: "<<delMePlz<<endl;
	delete delMePlz;

	cout<<"\tTesting known outside-outside line passing through sphere... ";
	int ct=sph->intersectionCount(ptB,ptC);
	cout<<(ct==2)<<endl;

	if (ct==2){
		delMePlz=sph->getBothIntersections(ptB,ptC);
		cout<<"\t Intersection points are: \n\t\tFirst:  "
			<<&delMePlz[0]<<endl<<"\t\tSecond: "
			<<&delMePlz[1]<<endl;
		delete[] delMePlz;
	}



	delete ptA; delete ptB; delete ptC;delete B;
	return sph;
}
void cupTest(Alpha::Cup * cup){
	cout<<"Testing Cup Functions... "<<endl;
	Alpha::Vertex *ptA= new Alpha::Vertex(0.1f,0.1f,0.1f), //Point inside the tetrahedron and not circle
		*ptB = new Alpha::Vertex(0.222f,0.222f,0.222f),//Point inside Tetrahedron and circl
		*ptC=new Alpha::Vertex(1.1f,0.0f,0.1f), //Point inside circle but not tetrahedron
		*ptD=new Alpha::Vertex(0.9,1.0,0.0), //Point outside Circle and Tet
		*ptE = new Alpha::Vertex(0.9,-1.0,0.0),
		*ptF = new Alpha::Vertex(-0.1,0,0), //Point -x of both, closer to Tet
		*ptG = new Alpha::Vertex(0,0.1, 0); //Boundary point on the tetrahedron
	cout<<"\tTesting known inside point... "<<cup->ptInside(ptA)
			<<"\n\tTesting known boundary point on tetrahedron... "<<cup->ptInside(ptG)
			<<"\n\tTesting known outside point to tetrahedron inside to circle... "
				<<cup->ptInside(ptC)
			<<"\n\tTesting unknown point... "<<cup->ptInside(ptD)
			<<"\n\tTesting known pt. Tet: inoutside point to tetrahedron inside to circle... "
			<<cup->ptInside(ptC)
			<<endl;
}
void spindleTest(){
	cout<<"Testing Spindle class..."<<endl;
	Alpha::Vertex ac=Alpha::Vertex(0,0,0)
			,bc=Alpha::Vertex(1,0,0)
			,pc=Alpha::Vertex(0.5,0.4,0)
			,test1=Alpha::Vertex(-1,0,0) //Known outside past A
	,test2=Alpha::Vertex(0.5,.1,0) //Known inside spindle
	,test3=Alpha::Vertex(0.5,0.3,0)  //Known inside torus and cylinder
	;

	Alpha::Sphere a=Alpha::Sphere(&ac, 0.4);
	Alpha::Sphere b=Alpha::Sphere(&bc, 0.4);
	Alpha::Sphere probe=Alpha::Sphere(&pc, 0.240312);
	Alpha::Spindle *spin=new Alpha::Spindle(&probe, &a, &b);
	cout<<"\tTesting known outside past A\n\t\tExpected = 0, Result = "
		<<spin->ptInside(&test1)<<endl;
	cout<<"\tTesting known inside spindle\n\t\tExpected = 1, Result = "
		<<spin->ptInside(&test2)<<endl;
	cout<<"\n\tTesting known inside torus and spindle"
		<<"\n\t\tExpected = 0, Result = "
		<<spin->ptInside(&test3)<<endl;

}
int main(int argc, char** arg) {
	std::string qhullOutput=buildAndPrint(argc, arg);
	parseData(qhullOutput);
	return 0;
	/*Alpha::Sphere *sph;
	Alpha::Tetrahedron *tet;
	Alpha::Cup *cp;
	*/

	//vertexTest();
	//spindleTest();
	//sph=sphereTest();
	//tet=tetTest();

	//cp=new Alpha::Cup(tet,sph);
	//cupTest(cp);


	/*
	//double probe_size=1.0;
	//bool spherepassthroughs=true;
	//string dat=" This is some text ";
	//string fname="/home/winnen/workspace/1TE6.pdb";

	//cout << clean(dat)<<endl;
	Alpha::Vertex *A, *B, *C, *D, *pt;
	Alpha::Vertex *pts[20];
	// *pts=new Alpha::Vertex[10];
	for (int i = 0; i < 10; i++) {
		pts[i] = new Alpha::Vertex(-1 + 0.2 * i, -1 + 0.2 * i, -1 + 0.2 * i);
	}

	Alpha::Tetrahedron *test;
	A = new Alpha::Vertex(0.0, 0.0, 0.0);
	B = new Alpha::Vertex(1.0, 0.0, 0.0);
	cout << A << endl;
	cout << "Distance from {0,0,0} to {1,0,0} = " << A->distance(B) << endl;
	C = new Alpha::Vertex(0.0, 1.0, 0.0);
	cout << "Distance from {1,0,0} to {0,1,0} = " << B->distance(C) << endl;
	D = new Alpha::Vertex(0.0, 0.0, 1.0);
	test = new Alpha::Tetrahedron(A, B, C, D);
	pt = new Alpha::Vertex(1.0, 0.0, 0.0);
	cout << "Testing if the point pt=" << pt
			<< "is inside the z=1-x-y tetrahedron. " << endl
			<< test->ptInside(pt) << endl;
	Alpha::Sphere *sph = new Alpha::Sphere(A, 1.0);

	for (int i = 0; i < 10; i++) {
		cout << "Testing if the point pt=" << pts[i]
				<< "is inside the z=1-x-y tetrahedron." << endl
				<< test->ptInside(pts[i]) << endl;
		delete pts[i];
	}
	//pt->arry[0]=0.0; pt->arry[1]=1.0; pt
	cout << "Testing if the point pt=" << pt << "is inside the z=1-r*r sphere. "
			<< endl << sph->ptInside(pt) << endl;
	//delete[] pts;
	//delete[](A,B,C,D,pt,test);
	delete test;
	delete A;
	delete B;
	delete C;
	delete D;
	delete pt;
	*/
	return 0;


}
