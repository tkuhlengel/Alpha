/*
 * mathfunc.cpp
 *
 *  Created on: Jul 24, 2012
 *      Author: winnen
 */
#include "mathfunc.h"

double vectorSize(double x, double y, double z) {
	return sqrt((x * x) + (y * y) + (z * z));
}
double vectorSize(double * A, double * B) {
	double x = B[0] - A[0];
	double y = B[1] - A[1];
	double z = B[2] - A[2];

	return sqrt((x * x) + (y * y) + (z * z));
}

////simple determinant for 3x3 matrices
double determinant3x3(double * matrix) {

	double det = 0;

	////calculate determinant
	det = matrix[0] * ((matrix[4] * matrix[8]) - (matrix[7] * matrix[5]));
	det = det - matrix[3] * ((matrix[1] * matrix[8]) - (matrix[7] * matrix[2]));
	det = det + matrix[6] * ((matrix[1] * matrix[5]) - (matrix[4] * matrix[2]));

	return det;

}

////simple determinant for 4x4 matrices
double determinant4x4(double * m) {
	double det = 0;

	// 0  4  8  12
	// 1  5  9  13
	// 2  6  10 14
	// 3  7  11 15

	double * d0 = new double[9];
	d0[0] = m[5];
	d0[1] = m[6];
	d0[2] = m[7];
	d0[3] = m[9];
	d0[4] = m[10];
	d0[5] = m[11];
	d0[6] = m[13];
	d0[7] = m[14];
	d0[8] = m[15];
	double * d1 = new double[9];
	d1[0] = m[1];
	d1[1] = m[2];
	d1[2] = m[3];
	d1[3] = m[9];
	d1[4] = m[10];
	d1[5] = m[11];
	d1[6] = m[13];
	d1[7] = m[14];
	d1[8] = m[15];
	double * d2 = new double[9];
	d2[0] = m[1];
	d2[1] = m[2];
	d2[2] = m[3];
	d2[3] = m[5];
	d2[4] = m[6];
	d2[5] = m[7];
	d2[6] = m[13];
	d2[7] = m[14];
	d2[8] = m[15];
	double * d3 = new double[9];
	d3[0] = m[1];
	d3[1] = m[2];
	d3[2] = m[3];
	d3[3] = m[5];
	d3[4] = m[6];
	d3[5] = m[7];
	d3[6] = m[9];
	d3[7] = m[10];
	d3[8] = m[11];

	////calculate determinant
	if (m[0] != 0) {
		det += m[0] * determinant3x3(d0);
	}
	if (m[4] != 0) {
		det -= m[4] * determinant3x3(d1);
	}
	if (m[8] != 0) {
		det += m[8] * determinant3x3(d2);
	}
	if (m[12] != 0) {
		det -= m[12] * determinant3x3(d3);
	}

	delete[] (d0);
	delete[] (d1);
	delete[] (d2);
	delete[] (d3);

	return det;
}

////simple determinant for 5x5 matrices
double determinant5x5(double * m) {
	double det = 0;

	// 0  5  10 15 20
	// 1  6  11 16 21
	// 2  7  12 17 22
	// 3  8  13 18 23
	// 4  9  14 19 24

	double * d0 = new double[16];
	d0[0] = m[6];
	d0[1] = m[7];
	d0[2] = m[8];
	d0[3] = m[9];
	d0[4] = m[11];
	d0[5] = m[12];
	d0[6] = m[13];
	d0[7] = m[14];
	d0[8] = m[16];
	d0[9] = m[17];
	d0[10] = m[18];
	d0[11] = m[19];
	d0[12] = m[21];
	d0[13] = m[22];
	d0[14] = m[23];
	d0[15] = m[24];
	double * d1 = new double[16];
	d1[0] = m[1];
	d1[1] = m[2];
	d1[2] = m[3];
	d1[3] = m[4];
	d1[4] = m[11];
	d1[5] = m[12];
	d1[6] = m[13];
	d1[7] = m[14];
	d1[8] = m[16];
	d1[9] = m[17];
	d1[10] = m[18];
	d1[11] = m[19];
	d1[12] = m[21];
	d1[13] = m[22];
	d1[14] = m[23];
	d1[15] = m[24];
	double * d2 = new double[16];
	d2[0] = m[1];
	d2[1] = m[2];
	d2[2] = m[3];
	d2[3] = m[4];
	d2[4] = m[6];
	d2[5] = m[7];
	d2[6] = m[8];
	d2[7] = m[9];
	d2[8] = m[16];
	d2[9] = m[17];
	d2[10] = m[18];
	d2[11] = m[19];
	d2[12] = m[21];
	d2[13] = m[22];
	d2[14] = m[23];
	d2[15] = m[24];
	double * d3 = new double[16];
	d3[0] = m[1];
	d3[1] = m[2];
	d3[2] = m[3];
	d3[3] = m[4];
	d3[4] = m[6];
	d3[5] = m[7];
	d3[6] = m[8];
	d3[7] = m[9];
	d3[8] = m[11];
	d3[9] = m[12];
	d3[10] = m[13];
	d3[11] = m[14];
	d3[12] = m[21];
	d3[13] = m[22];
	d3[14] = m[23];
	d3[15] = m[24];
	double * d4 = new double[16];
	d4[0] = m[1];
	d4[1] = m[2];
	d4[2] = m[3];
	d4[3] = m[4];
	d4[4] = m[6];
	d4[5] = m[7];
	d4[6] = m[8];
	d4[7] = m[9];
	d4[8] = m[11];
	d4[9] = m[12];
	d4[10] = m[13];
	d4[11] = m[14];
	d4[12] = m[16];
	d4[13] = m[17];
	d4[14] = m[18];
	d4[15] = m[19];

	////calculate determinant
	if (m[0] != 0) {
		det += m[0] * determinant4x4(d0);
	}
	if (m[5] != 0) {
		det -= m[5] * determinant4x4(d1);
	}
	if (m[10] != 0) {
		det += m[10] * determinant4x4(d2);
	}
	if (m[15] != 0) {
		det -= m[15] * determinant4x4(d3);
	}
	if (m[20] != 0) {
		det += m[20] * determinant4x4(d4);
	}

	delete[] (d0);
	delete[] (d1);
	delete[] (d2);
	delete[] (d3);
	delete[] (d4);

	return det;
}
double * intersect_Seg_Triangle(double * r0, double * r1, double * t0,
		double * t1, double * t2) {
	double * result = NULL;

	double ba0 = t1[0] - t0[0];
	double ba1 = t1[1] - t0[1];
	double ba2 = t1[2] - t0[2];
	double ca0 = t2[0] - t0[0];
	double ca1 = t2[1] - t0[1];
	double ca2 = t2[2] - t0[2];
	double norm0 = (ba1 * ca2) - (ba2 * ca1);
	double norm1 = -((ba0 * ca2) - (ba2 * ca0));
	double norm2 = (ba0 * ca1) - (ba1 * ca0);

	/////////////////////////////////////////////////////////////////////////
	///DEGENERACY CHECK FIRST:///////////////////////////////////////////////
	///deal with degenerate triangles: triangles that are essentially lines or points
	double check_area = vectorSize(norm0, norm1, norm2);
	double check_d1 = vectorSize(t0, t1);
	double check_d2 = vectorSize(t1, t2);
	double check_d3 = vectorSize(t2, t0);

//	printf("area: %f d1: %f d2: %f d3: %f\n", check_area, check_d1, check_d2, check_d3 );
//	printf("checkarea: %i  check_d1 %i  check_d2 %i  check_d3 %i\n", check_area<SMALL_NUM, check_d1<SMALL_NUM, check_d2<SMALL_NUM, check_d3<SMALL_NUM );
//	printf("checkarea: %f  check_d1 %f  check_d2 %f  check_d3 %f\n", check_area, check_d1, check_d2, check_d3 );

	if (check_area < SMALL_NUM) {
		//if this is a point-triangle, take any point of the triangle, see if it intersects.
		if (check_d1 < SMALL_NUM && check_d2 < SMALL_NUM && check_d3 < SMALL_NUM) {
			bool check_pt = pointSegmentOverlapTest(t0, r0, r1);
			if (check_pt) {
				double * check_result = new double[3];
				check_result[0] = t0[0];
				check_result[1] = t0[1];
				check_result[2] = t0[2];
				return check_result;
			} else {
				return result;
			}
		}
		//if this is a line-triangle, then
		else {
			result = intersect_Triangle_Segments(t0, t1, t2, r0, r1);
			return result;
		}
	}
	///DEGENERACY CHECK FIRST:///////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	/*
	 //	t1
	 //   | \
	//	|  \
	//   u   \
	//	|    \
	//	t0-v-t2
	 */
	///////////////////////////////////////////////////////////////////////
	///first check if the segment even intersects the plane of the triangle.
	double mag;
	mag = sqrt((ba0 * ba0) + (ba1 * ba1) + (ba2 * ba2));
	double u0 = ba0 / mag;
	double u1 = ba1 / mag;
	double u2 = ba2 / mag;
	mag = sqrt((ca0 * ca0) + (ca1 * ca1) + (ca2 * ca2));
	double v0 = ca0 / mag;
	double v1 = ca1 / mag;
	double v2 = ca2 / mag;

	double uvCross0 = (u1 * v2) - (u2 * v1);
	double uvCross1 = -((u0 * v2) - (u2 * v0));
	double uvCross2 = (u0 * v1) - (u1 * v0);

	mag = sqrt(
			(uvCross0 * uvCross0) + (uvCross1 * uvCross1)
					+ (uvCross2 * uvCross2));
	double tnorm0 = uvCross0 / mag;
	double tnorm1 = uvCross1 / mag;
	double tnorm2 = uvCross2 / mag;

	mag = sqrt(
			((r0[0] - t0[0]) * (r0[0] - t0[0]))
					+ ((r0[1] - t0[1]) * (r0[1] - t0[1]))
					+ ((r0[2] - t0[2]) * (r0[2] - t0[2])));
	if (mag == 0) {
		double * result = new double[3];
		result[0] = t0[0];
		result[1] = t0[1];
		result[0] = t0[2];
		return result;
	}
	double r0t0 = (r0[0] - t0[0]) / mag; ///vect0
	double r0t1 = (r0[1] - t0[1]) / mag;
	double r0t2 = (r0[2] - t0[2]) / mag;
	mag = sqrt(
			((r1[0] - t0[0]) * (r1[0] - t0[0]))
					+ ((r1[1] - t0[1]) * (r1[1] - t0[1]))
					+ ((r1[2] - t0[2]) * (r1[2] - t0[2])));
	if (mag == 0) {
		double * result = new double[3];
		result[0] = t0[0];
		result[1] = t0[1];
		result[0] = t0[2];
		return result;
	}
	double r1t0 = (r1[0] - t0[0]) / mag; ///vect1
	double r1t1 = (r1[1] - t0[1]) / mag;
	double r1t2 = (r1[2] - t0[2]) / mag;

//	printf("ASD tnorm: %f %f %f   vect0: %f %f %f (normed from: %f)  vect1: %f %f %f (normed from: %f)\n",
//		tnorm[0], tnorm[1], tnorm[2],
//		vect0[0], vect0[1], vect0[2], vectorSize( r0[0]-t0[0], r0[1]-t0[1], r0[2]-t0[2] ),
//		vect1[0], vect1[1], vect1[2], vectorSize( r1[0]-t0[0], r1[1]-t0[1], r1[2]-t0[2] )  );

	///////Use the dot product to dtermine if the vector is with, against, or perpendicular to the normal
	double test0 = (tnorm0 * r0t0) + (tnorm1 * r0t1) + (tnorm2 * r0t2);
	double test1 = (tnorm0 * r1t0) + (tnorm1 * r1t1) + (tnorm2 * r1t2);

	///eliminate issues with small numbers
	if (fabs(test0) < SMALL_NUM) {
		test0 = 0.0;
	}
	if (fabs(test1) < SMALL_NUM) {
		test1 = 0.0;
	}

//	printf("test0: %f test1: %f\n", test0, test1);

	///////if the dot product is 0, then it is perpendicular, and the point is on the plane.
	///////if the dot product is positive, then it is on the same side as the normal
	///////if the dot product is negative, then it is on the opposite side of the triangle from the normal.
	int triangleCase = -1;
	if ((test0 < 0) && (test1 < 0)) {
		triangleCase = 0;
	}
	if ((test0 < 0) && (test1 == 0)) {
		triangleCase = 1;
	}
	if ((test0 < 0) && (test1 > 0)) {
		triangleCase = 2;
	}
	if ((test0 == 0) && (test1 < 0)) {
		triangleCase = 3;
	}
	if ((test0 == 0) && (test1 == 0)) {
		triangleCase = 4;
	}
	if ((test0 == 0) && (test1 > 0)) {
		triangleCase = 5;
	}
	if ((test0 > 0) && (test1 < 0)) {
		triangleCase = 6;
	}
	if ((test0 > 0) && (test1 == 0)) {
		triangleCase = 7;
	}
	if ((test0 > 0) && (test1 > 0)) {
		triangleCase = 8;
	}

	if (triangleCase == -1) {
		printf("WTF!! triangleCase == -1 TEST0: [%f] TEST1: [%f]\n", test0,
				test1);
		printf("ba0 %f  ba1 %f  ba2 %f  ca0 %f  ca1 %f  ca2 %f\n", ba0, ba1,
				ba2, ca0, ca1, ca2);
		printf(
				"TEST0== tnorm0: %f r0t0: %f tnorm1: %f r0t1: %f tnorm2: %f r0t2: %f   mag: %f\n",
				tnorm0, r0t0, tnorm1, r0t1, tnorm2, r0t2,
				sqrt(
						((r0[0] - t0[0]) * (r0[0] - t0[0]))
								+ ((r0[1] - t0[1]) * (r0[1] - t0[1]))
								+ ((r0[2] - t0[2]) * (r0[2] - t0[2]))));
		printf(
				"TEST1== tnorm0: %f r1t0: %f tnorm1: %f r1t1: %f tnorm2: %f r1t2: %f   mag: %f\n",
				tnorm0, r1t0, tnorm1, r1t1, tnorm2, r1t2,
				sqrt(
						((r1[0] - t0[0]) * (r1[0] - t0[0]))
								+ ((r1[1] - t0[1]) * (r1[1] - t0[1]))
								+ ((r1[2] - t0[2]) * (r1[2] - t0[2]))));
		exit(1);
	}
	//printf("TRIANGLE CASE: %i  test0: %f test1: %f\n", triangleCase, test0, test1);	///debug instrumentation - DO NOT DELETE

	////////return NULL for any cases where the segment is clearly off plane (i.e. no intersection)
	if (triangleCase == 0 || triangleCase == 8) {
		return NULL;
	}

	///Deal with the case where the segment is totally coplanar
	if (triangleCase == 4) {
		result = intersect_Triangle_Segments(t0, t1, t2, r0, r1);
		return result;
	}

	///////////////////////////////////////////////////////////////////////
	///second, identify the point of intersection on the plane.
	double planePt0, planePt1, planePt2;

	///in the cases where there is exactly one point in the plane, this is well defined.
	if (triangleCase == 3 || triangleCase == 5) { ///in this case, r0 is in the plane
		planePt0 = r0[0];
		planePt1 = r0[1];
		planePt2 = r0[2];
	}

	if (triangleCase == 1 || triangleCase == 7) { ///in this case, r1 is in the plane
		planePt0 = r1[0];
		planePt1 = r1[1];
		planePt2 = r1[2];
	}

	///The most common case is when the points are on either side of the triangle.
	if (triangleCase == 2 || triangleCase == 6) {
		///
		///treat the segment as a ray.  We will identify l, the paramaterization variable along this ray,
		///where this ray hits the plane.  We will say the ray hits the plane at point p.
		///
		///This computed by solving two simultaneous equations:
		///Equation 1:
		///   tnorm dot (p-t0) = 0     "the dot product of the normal and an in-plane vector is zero."
		///
		///Equation 2:
		///   p = r0 + l(r1-r0)        "the point at intersection is along the line of parametrization"
		///
		///Substituting the right side of equation 2 for p in equation 1, gives
		///   tnorm dot ( (r0 + l(r1-r0)) - t0 ) = 0
		///
		///Or, more simply:
		///   tnorm dot( r0 + l(r1-r0) ) = tnorm dot t0
		///
		///Solving for l, this gives
		///
		///          tnorm dot (t0-r0)
		///   l = -----------------------
		///          tnorm dot (r1-r0)
		///

		double topValue = (tnorm0 * (t0[0] - r0[0]))
				+ (tnorm1 * (t0[1] - r0[1])) + (tnorm2 * (t0[2] - r0[2]));
		double botValue = (tnorm0 * (r1[0] - r0[0]))
				+ (tnorm1 * (r1[1] - r0[1])) + (tnorm2 * (r1[2] - r0[2]));
		double l = (topValue) / (botValue);

		//	if( fabs(topValue) < SMALL_NUM ){ topValue = 0; printf("topValue is TINY [%f] (triangle area: %f)\n", topValue, check_area); }
		//	if( fabs(botValue) < SMALL_NUM ){ botValue = 0; printf("botValue is TINY [%f] (triangle area: %f)\n", botValue, check_area); }
		//	if( fabs(l) < SMALL_NUM ){ l = 0; printf("l = [%f] top: [%f] bot: [%f] (triangle area: %f\n", l, topValue, botValue, check_area); }

//		///eliminate issues with small numbers
		if (fabs(l) < SMALL_NUM) {
			test0 = 0.0;
		}

		///if l<0 or l>1, then the parametrization is outside the segment, and they do not touch.
		if ((l >= 0) && (l <= 1)) {
			///we solve the parametrization, if it is within [0,1]
			///simply plug in l, now that we have it.
			planePt0 = r0[0] + l * (r1[0] - r0[0]);
			planePt1 = r0[1] + l * (r1[1] - r0[1]);
			planePt2 = r0[2] + l * (r1[2] - r0[2]);
		} else {
			return NULL;
		}
	}

//	printf("planePt: %f %f %f\n", planePt[0], planePt[1], planePt[2] );

	///////////////////////////////////////////////////////////////////////
	///third identify if the point is inside the triangle.
	///
	///We will identify if the point is inside the triangle by checking half-planes.
	///If the point is on the correct side of each of the triangle edges, then it is inside.
	///
	///tnorm is (t1-t0) CROSS (t2-t0).  Thus, to determine if the point is
	///inside the triangle,
	///
	///i0 = tnorm CROSS (t0-t1) gives an in-plane vector pointing inwards of the triangle.
	/// (p-t0) DOT (VECTOR) is positive if the point is on the correct side, negative otherwise.
	///
	///Repeat this process for all edges.
	///
	double hp0, hp1, hp2;
	double tmp0, tmp1, tmp2;

	///crossproduct with first half plane
	tmp0 = ((t0[1] - t1[1]) * tnorm2) - ((t0[2] - t1[2]) * tnorm1);
	tmp1 = -(((t0[0] - t1[0]) * tnorm2) - ((t0[2] - t1[2]) * tnorm0));
	tmp2 = ((t0[0] - t1[0]) * tnorm1) - ((t0[1] - t1[1]) * tnorm0);
	hp0 = (tmp0 * (planePt0 - t0[0])) + (tmp1 * (planePt1 - t0[1]))
			+ (tmp2 * (planePt2 - t0[2]));

	///crossproduct with second half plane
	tmp0 = ((t1[1] - t2[1]) * tnorm2) - ((t1[2] - t2[2]) * tnorm1);
	tmp1 = -(((t1[0] - t2[0]) * tnorm2) - ((t1[2] - t2[2]) * tnorm0));
	tmp2 = ((t1[0] - t2[0]) * tnorm1) - ((t1[1] - t2[1]) * tnorm0);
	hp1 = (tmp0 * (planePt0 - t1[0])) + (tmp1 * (planePt1 - t1[1]))
			+ (tmp2 * (planePt2 - t1[2]));

	///crossproduct with third half plane
	tmp0 = ((t2[1] - t0[1]) * tnorm2) - ((t2[2] - t0[2]) * tnorm1);
	tmp1 = -(((t2[0] - t0[0]) * tnorm2) - ((t2[2] - t0[2]) * tnorm0));
	tmp2 = ((t2[0] - t0[0]) * tnorm1) - ((t2[1] - t0[1]) * tnorm0);
	hp2 = (tmp0 * (planePt0 - t2[0])) + (tmp1 * (planePt1 - t2[1]))
			+ (tmp2 * (planePt2 - t2[2]));

	///fix small-value precision issues
	if (fabs(hp0) < SMALL_NUM) {
		hp0 = 0.0;
	}
	if (fabs(hp1) < SMALL_NUM) {
		hp1 = 0.0;
	}
	if (fabs(hp2) < SMALL_NUM) {
		hp2 = 0.0;
	}

	////this might be a numerical issue with hitting the edge (hp0=0, for example)
	if ((hp0 >= 0) && (hp1 >= 0) && (hp2 >= 0)) {
		result = new double[3];
		result[0] = planePt0;
		result[1] = planePt1;
		result[2] = planePt2;
	} else {
		///if not all cases are true, then we are not inside the triangle.  clear all and return NULL.
		result = NULL;
	}

	return result;
}

bool pointSegmentOverlapTest(Alpha::Vertex * pt, Alpha::Vertex * s1, Alpha::Vertex * s2) {
	Alpha::Vertex * t1;
	t1 = pt->minus(s1);


	Alpha::Vertex * t2;
	t2 = s2->minus(s1);

	bool result = false;
	double dprod = DOT(t1->arry, t2->arry);

	delete t1;
	delete t2;
	if (dprod > 0 && dprod < 1) {
		result = true;
	}

	return result;
}
double * intersect_Triangle_Segments(Alpha::Vertex * t1, Alpha::Vertex * t2, Alpha::Vertex * t3,
		Alpha::Vertex * s0, Alpha::Vertex * s1) {

	//Copy a few Vertices to make them triangle relevant
	Alpha::Vertex * s1p1, *s1p2, *s2p1, *s2p2, *s3p1, *s3p2;
	s1p1=Alpha::Vertex(t1);
	s1p2=Alpha::Vertex(t2);

	s2p1=Alpha::Vertex(t2);
	s2p2=Alpha::Vertex(t3);

	s3p1=Alpha::Vertex(t3);
	s3p2=Alpha::Vertex(t1);


	Alpha::Vertex * closestPt = new Alpha::Vertex();
	Alpha::Vertex * closePt = new Alpha::Vertex();

	double dist;
	double minDist = HUGE_VAL;

	if (s1p1->distance(s1p2) > 0) { ///check for positive distance, to avoid degeneracy.
		dist = dist3D_Segment_to_Segment(s1p1, s1p2, s0, s1, closePt);
		if (dist < minDist) {
			closestPt = closePt;
			//losestPt[1] = closePt[1];
			//closestPt[2] = closePt[2];
			minDist = dist;
		}
	}
	if (vectorSize(s2p1, s2p2) > 0) { ///check for positive distance, to avoid degeneracy.
		dist = dist3D_Segment_to_Segment(s2p1, s2p2, s0, s1, closePt);
		if (dist < minDist) {
			closestPt[0] = closePt[0];
			closestPt[1] = closePt[1];
			closestPt[2] = closePt[2];
			minDist = dist;
		}
	}
	if (vectorSize(s3p1, s3p2) > 0) { ///check for positive distance, to avoid degeneracy.
		dist = dist3D_Segment_to_Segment(s3p1->arry, s3p2->arry, s0->arry, s1->arry, closePt->arry);
		if (dist < minDist) {
			closestPt[0] = closePt[0];
			closestPt[1] = closePt[1];
			closestPt[2] = closePt[2];
			minDist = dist;
		}
	}

	delete[] (s1p1);
	delete[] (s1p2);
	delete[] (s2p1);
	delete[] (s2p2);
	delete[] (s3p1);
	delete[] (s3p2);
	delete[] (closePt);

	if (minDist < SMALL_NUM) {
		return closestPt;
	}
	delete[] (closestPt);
	return NULL;

}

double dist3D_Segment_to_Segment(double * s1p0, double * s1p1, double * s2p0,
		double * s2p1, double * closest) {
	//Vector   u = S1.P1 - S1.P0;
	double * u = new double[3];
	u[0] = s1p1[0] - s1p0[0];
	u[1] = s1p1[1] - s1p0[1];
	u[2] = s1p1[2] - s1p0[2];

	//Vector   v = S2.P1 - S2.P0;
	double * v = new double[3];
	v[0] = s2p1[0] - s2p0[0];
	v[1] = s2p1[1] - s2p0[1];
	v[2] = s2p1[2] - s2p0[2];

	//Vector   w = S1.P0 - S2.P0;
	double * w = new double[3];
	w[0] = s1p0[0] - s2p0[0];
	w[1] = s1p0[1] - s2p0[1];
	w[2] = s1p0[2] - s2p0[2];

	double a = DOT(u,u);        // always >= 0
	double b = DOT(u,v);
	double c = DOT(v,v);        // always >= 0
	double d = DOT(u,w);
	double e = DOT(v,w);
	double D = a * c - b * b;       // always >= 0
	double sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
	double tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0

	// compute the line parameters of the two closest points
	if (D < SMALL_NUM) {  // the lines are almost parallel
		sN = 0.0;        // force using point P0 on segment S1
		sD = 1.0;        // to prevent possible division by 0.0 later
		tN = e;
		tD = c;
	}
	// get the closest points on the infinite lines
	else {
		sN = (b * e - c * d);
		tN = (a * e - b * d);
		if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
			sN = 0.0;
			tN = e;
			tD = c;
		} else {
			if (sN > sD) {  // sc > 1 => the s=1 edge is visible
				sN = sD;
				tN = e + b;
				tD = c;
			}
		}
	}

	if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
		tN = 0.0;
		// recompute sc for this edge
		if (-d < 0.0) {
			sN = 0.0;
		} else {
			if (-d > a) {
				sN = sD;
			} else {
				sN = -d;
				sD = a;
			}
		}
	} else {
		// tc > 1 => the t=1 edge is visible
		if (tN > tD) {
			tN = tD;
			// recompute sc for this edge
			if ((-d + b) < 0.0) {
				sN = 0;
			} else {
				if ((-d + b) > a) {
					sN = sD;
				} else {
					sN = (-d + b);
					sD = a;
				}
			}
		}
	}

	// finally do the division to get sc and tc
	if (fabs(sN) < SMALL_NUM) {
		sc = 0.0;
	} else {
		sc = (sN / sD);
	}

	if (fabs(tN) < SMALL_NUM) {
		tc = 0.0;
	} else {
		tc = (tN / tD);
	}

	// get the difference of the two closest points
	//Vector   dP = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)
	double * dP = new double[3];
	dP[0] = w[0] + sc * u[0] - tc * v[0];
	dP[1] = w[1] + sc * u[1] - tc * v[1];
	dP[2] = w[2] + sc * u[2] - tc * v[2];

	// Now we assign the closest point
	if (closest != NULL) {
		closest[0] = tc * v[0];
		closest[1] = tc * v[1];
		closest[2] = tc * v[2];
	}

	// return the closest distance
	double result = vectorSize(dP[0], dP[1], dP[2]);

	delete[] (u);
	delete[] (v);
	delete[] (w);
	delete[] (dP);

	return result;
}
//Meno: Prof. Pasko's handout pp.8 in "Skelton Model"
double distanceToSegment(double * A, double * B, double * P) {
	double result;
	result = distanceToSegmentFast(A, B, P);
	return sqrt(result);

}

//Meno: Prof. Pasko's handout pp.8 in "Skelton Model"
double distanceToSegmentFast(double * A, double * B, double * P) {
	double * AC = new double[3];
	double result;

	double AB2 = (B[0] - A[0]) * (B[0] - A[0])
			+ (B[1] - A[1]) * (B[1] - A[1])
			+ (B[2] - A[2]) * (B[2] - A[2]);
	double dot = (P[0] - A[0]) * (B[0] - A[0])
			+ (P[1] - A[1]) * (B[1] - A[1])
			+ (P[2] - A[2]) * (B[2] - A[2]);

	AC[0] = (double) (dot * (B[0] - A[0]) / AB2);
	AC[1] = (double) (dot * (B[1] - A[1]) / AB2);
	AC[2] = (double) (dot * (B[2] - A[2]) / AB2);

	double dot2 = AC[0] * (B[0] - A[0])
			+ AC[1] * (B[1] - A[1])
			+ AC[2] * (B[2] - A[2]);

	if (dot2 < 0) {
		result = (A[0] - P[0]) * (A[0] - P[0])
				+ (A[1] - P[1]) * (A[1] - P[1])
				+ (A[2] - P[2]) * (A[2] - P[2]);
	} else {
		if (AB2 < AC[0] * AC[0] + AC[1] * AC[1] + AC[2] * AC[2]) {
			result = (B[0] - P[0]) * (B[0] - P[0])
					+ (B[1] - P[1]) * (B[1] - P[1])
					+ (B[2] - P[2]) * (B[2] - P[2]);
		} else {
			result = (AC[0] + A[0] - P[0]) * (AC[0] + A[0] - P[0])
					+ (AC[1] + A[1] - P[1]) * (AC[1] + A[1] - P[1])
					+ (AC[2] + A[2] - P[2]) * (AC[2] + A[2] - P[2]);
		}
	}

	delete[] (AC);
	return result;
}

