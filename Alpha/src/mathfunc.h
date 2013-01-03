/*
 * mathfunc.h
 *
 *  Created on: Jul 24, 2012
 *      Author: winnen
 *
 *      Functions taken from Brian Chen's code (byc210@lehigh.edu)
 *      to make compilation easier
 */

#ifndef MATHFUNC_H_
#define MATHFUNC_H_
#include <cmath>
#include <cstdio>
#include <cstdlib>
#ifndef SMALL_NUM
#define SMALL_NUM  0.000001
#endif

#ifndef NULL
#define NULL   ((void *) 0)
#endif
typedef unsigned int word_t;

typedef unsigned char byte_t;
typedef char * string_t;
typedef byte_t * chunk_t;
typedef int boolean_t_new;

double determinant3x3(double * matrix);
double determinant4x4(double * m);
double * intersect_Seg_Triangle(double * r0, double * r1, double * t0,
		double * t1, double * t2);
double determinant5x5(double * m);
double vectorSize(double * A, double * B);
double vectorSize(double x, double y, double z);
bool pointSegmentOverlapTest(double * pt, double * s1, double * s2);
double * intersect_Triangle_Segments(double * t1, double * t2, double * t3, double * s0, double * s1);
double dist3D_Segment_to_Segment( double * s1p0, double * s1p1, double * s2p0, double * s2p1, double * closest);
/**
 * \return The distance of the closest points between the two lines
 */
/** \brief A distance square calculator
 * \param A[in] Array of length 3 denoting a line
 * \param
 * \return double representing the square of the distance between line AB and P.
 */
double distanceToSegmentFast(double * A, double * B, double * P);

double distanceToSegment(double * A, double * B, double * P);
/**
 * Returns the absolute shortest distance between a segment and
 */
typedef unsigned char byte_t;
#define get_set_head(s) (set_head_t *) (((byte_t *) s) - sizeof(set_head_t))
typedef struct
{
  int capacity;
  int size;

  int table_size;

  char* properties;

  void * header;
} set_head_t;
void free_set(int* s){
  set_head_t *sh = get_set_head(s);

  free(sh);
}
#ifndef CROSS
#define CROSS(dest,v1,v2)                       \
               dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
               dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
               dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
#endif

///Dot Product
#ifndef DOT
#define DOT(v1,v2) ((v1[0]*v2[0]) + (v1[1]*v2[1]) + (v1[2]*v2[2]))
#endif
#endif /* MATHFUNC_H_ */
