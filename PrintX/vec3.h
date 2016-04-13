/*
* Original work Copyright (c) 2011, Jernej Barbic, Yili Zhao, University of Southern California
* Modified work Copyright (c) 2015, Miaojun Yao
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of University of Southern California, nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY JERNEJ BARBIC, YILI ZHAO AND UNIVERSITY OF SOUTHERN CALIFORNIA
* ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL JERNEJ BARBIC, YILI ZHAO OR UNIVERSITY OF SOUTHERN CALIFORNIA BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

3d vector operations (and a few 3d triangle operations)
*/

#ifndef __PRINTX_VEC3_H__
#define __PRINTX_VEC3_H__

#include <math.h>
#include <algorithm>
#include <iostream>

namespace PRINTX_LIB
{
	// a = b
#define VECTOR_SET3(a,b)\
	{(a)[0] = (b)[0];\
	(a)[1] = (b)[1];\
	(a)[2] = (b)[2];}

	template <class TYPE>
	inline void VectorSet3(TYPE *a, TYPE *b)
	{
		a[0] = b[0];
		a[1] = b[1];
		a[2] = b[2];
	}
#define VectorSet3d VectorSet3<double>
#define VectorSet3f VectorSet3<float>
#define VectorSet3i VectorSet3<int>


	// v[0] = a; v[1] = b; v[2] = c;
#define VECTOR_SET3S(v,a,b,c)\
	{(v)[0] = (a);\
	(v)[1] = (b);\
	(v)[2] = (c);}

	template <class TYPE>
	inline void VectorSet3S(TYPE *v, TYPE a, TYPE b, TYPE c)
	{
		v[0] = a;
		v[1] = b;
		v[2] = c;
	}
#define VectorSet3Sd VectorSet3S<double>
#define VectorSet3Sf VectorSet3S<float>
#define VectorSet3Si VectorSet3S<int>


#define VECTOR_SET_CONST3(a,b)\
	{(a)[0] = (b);\
	(a)[1] = (b);\
	(a)[2] = (b);}

	template <class TYPE>
	inline void VectorSetConst3(TYPE *a, TYPE b)
	{
		a[0] = b;
		a[1] = b;
		a[2] = b;
	}
#define VectorSetConst3d VectorSetConst3<double>
#define VectorSetConst3f VectorSetConst3<float>
#define VectorSetConst3i VectorSetConst3<int>


#define VECTOR_MIN3(a,b,c)\
	{(c)[0] = std::min((a)[0], (b)[0]);\
	(c)[1] = std::min((a)[1], (b)[1]);\
	(c)[2] = std::min((a)[2], (b)[2]);}

	template <class TYPE>
	inline void VectorMin3(TYPE *a, TYPE *b, TYPE *c)
	{
		(c)[0] = std::min((a)[0], (b)[0]);
		(c)[1] = std::min((a)[1], (b)[1]);
		(c)[2] = std::min((a)[2], (b)[2]);
	}
#define VectorMin3d VectorMin3<double>
#define VectorMin3f VectorMin3<float>
#define VectorMin3i VectorMin3<int>



#define VECTOR_MAX3(a,b,c)\
	{(c)[0] = std::max((a)[0], (b)[0]);\
	(c)[1] = std::max((a)[1], (b)[1]);\
	(c)[2] = std::max((a)[2], (b)[2]);}

	template <class TYPE>
	inline void VectorMax3(TYPE *a, TYPE *b, TYPE *c)
	{
		(c)[0] = std::max((a)[0], (b)[0]);
		(c)[1] = std::max((a)[1], (b)[1]);
		(c)[2] = std::max((a)[2], (b)[2]);
	}
#define VectorMax3d VectorMax3<double>
#define VectorMax3f VectorMax3<float>
#define VectorMax3i VectorMax3<int>



	// c = a + b
#define VECTOR_ADD3(a,b,c)\
	{(c)[0] = (a)[0] + (b)[0];\
	(c)[1] = (a)[1] + (b)[1];\
	(c)[2] = (a)[2] + (b)[2];}

	template <class TYPE>
	inline void VectorAdd3(TYPE *a, TYPE *b, TYPE *c)
	{
		(c)[0] = (a)[0] + (b)[0];
		(c)[1] = (a)[1] + (b)[1];
		(c)[2] = (a)[2] + (b)[2];
	}
#define VectorAdd3d VectorAdd3<double>
#define VectorAdd3f VectorAdd3<float>
#define VectorAdd3i VectorAdd3<int>

	// a += b
#define VECTOR_ADD_EQUAL3(a,b)\
			{(a)[0] += (b)[0];\
	(a)[1] += (b)[1];\
	(a)[2] += (b)[2];}

	template <class TYPE>
	inline void VectorAddEqual3(TYPE *a, TYPE *b)
	{
		(a)[0] += (b)[0];
		(a)[1] += (b)[1];
		(a)[2] += (b)[2];
	}
#define VectorAddEqual3d VectorAddEqual3<double>
#define VectorAddEqual3f VectorAddEqual3<float>
#define VectorAddEqual3i VectorAddEqual3<int>

	// c = a + b
#define VECTOR_ADD_CONST3(a,b,c)\
	{(c)[0] = (a)[0] + (b);\
	(c)[1] = (a)[1] + (b);\
	(c)[2] = (a)[2] + (b);}

	template <class TYPE>
	inline void VectorAddConst3(TYPE *a, TYPE b, TYPE *c)
	{
		(c)[0] = (a)[0] + (b);
		(c)[1] = (a)[1] + (b);
		(c)[2] = (a)[2] + (b);
	}
#define VectorAddConst3d VectorAddConst3<double>
#define VectorAddConst3f VectorAddConst3<float>
#define VectorAddConst3i VectorAddConst3<int>


	// a += b
#define VECTOR_ADD_CONST_EQUAL3(a,b)\
					{(a)[0] += (b);\
	(a)[1] += (b);\
	(a)[2] += (b);}

	template <class TYPE>
	inline void VectorAddConstEqual3(TYPE *a, TYPE *b)
	{
		(a)[0] += (b);
		(a)[1] += (b);
		(a)[2] += (b);
	}
#define VectorAddConstEqual3d VectorAddConstEqual3<double>
#define VectorAddConstEqual3f VectorAddConstEqual3<float>
#define VectorAddConstEqual3i VectorAddConstEqual3<int>


	// c = a + scalar * b
#define VECTOR_SCALE_ADD3(a,scalar,b,c)\
	{(c)[0] = (a)[0] + (scalar) * (b)[0];\
	(c)[1] = (a)[1] + (scalar) * (b)[1];\
	(c)[2] = (a)[2] + (scalar) * (b)[2];}

	template <class TYPE>
	inline void VectorScaleAdd3(TYPE *a, TYPE scalar, TYPE *b, TYPE *c)
	{
		(c)[0] = (a)[0] + (scalar)* (b)[0];
		(c)[1] = (a)[1] + (scalar)* (b)[1];
		(c)[2] = (a)[2] + (scalar)* (b)[2];
	}
#define VectorScaleAdd3d VectorScaleAdd3<double>
#define VectorScaleAdd3f VectorScaleAdd3<float>
#define VectorScaleAdd3i VectorScaleAdd3<int>


	// c = a + scalar * b
#define VECTOR_SCALE_ADD_EQUAL3(a,scalar,b)\
			{(a)[0] += (scalar) * (b)[0];\
	(a)[1] += (scalar) * (b)[1];\
	(a)[2] += (scalar) * (b)[2];}

	template <class TYPE>
	inline void VectorScaleAddEqual3(TYPE *a, TYPE scalar, TYPE *b)
	{
		(a)[0] += (scalar)* (b)[0];
		(a)[1] += (scalar)* (b)[1];
		(a)[2] += (scalar)* (b)[2];
	}
#define VectorScaleAddEqual3d VectorScaleAddEqual3<double>
#define VectorScaleAddEqual3f VectorScaleAddEqual3<float>
#define VectorScaleAddEqual3i VectorScaleAddEqual3<int>


	// c = scalar1 * a + scalar2 * b
#define VECTOR_SCALE2_ADD3(scalar1, a,scalar2,b,c)\
	{(c)[0] = (scalar1) * (a)[0] + (scalar2) * (b)[0];\
	(c)[1] = (scalar1) * (a)[1] + (scalar2) * (b)[1];\
	(c)[2] = (scalar1) * (a)[2] + (scalar2) * (b)[2];}

	template <class TYPE>
	inline void VectorScale2Add3(TYPE scalar1, TYPE *a, TYPE scalar2, TYPE *b, TYPE *c)
	{
		(c)[0] = (scalar1) * (a)[0] + (scalar2)* (b)[0];
		(c)[1] = (scalar1) * (a)[1] + (scalar2)* (b)[1];
		(c)[2] = (scalar1) * (a)[2] + (scalar2)* (b)[2];
	}
#define VectorScale2Add3d VectorScale2Add3<double>
#define VectorScale2Add3f VectorScale2Add3<float>
#define VectorScale2Add3i VectorScale2Add3<int>


	// c = a - b
#define VECTOR_SUBTRACT3(a,b,c)\
	{(c)[0] = (a)[0] - (b)[0];\
	(c)[1] = (a)[1] - (b)[1];\
	(c)[2] = (a)[2] - (b)[2];}

	template <class TYPE>
	inline void VectorSubtract3(TYPE *a, TYPE *b, TYPE *c)
	{
		(c)[0] = (a)[0] - (b)[0];
		(c)[1] = (a)[1] - (b)[1];
		(c)[2] = (a)[2] - (b)[2];
	}
#define VectorSubtract3d VectorSubtract3<double>
#define VectorSubtract3f VectorSubtract3<float>
#define VectorSubtract3i VectorSubtract3<int>

	// a -= b
#define VECTOR_SUBTRACT_EQUAL3(a,b)\
			{(a)[0] -= (b)[0];\
	(a)[1] -= (b)[1];\
	(a)[2] -= (b)[2];}

	template <class TYPE>
	inline void VectorSubtractEqual3(TYPE *a, TYPE *b)
	{
		(a)[0] -= (b)[0];
		(a)[1] -= (b)[1];
		(a)[2] -= (b)[2];
	}
#define VectorSubtractEqual3d VectorSubtractEqual3<double>
#define VectorSubtractEqual3f VectorSubtractEqual3<float>
#define VectorSubtractEqual3i VectorSubtractEqual3<int>


	// c = a - b
#define VECTOR_SUBTRACT_CONST3(a,b,c)\
	{(c)[0] = (a)[0] - (b);\
	(c)[1] = (a)[1] - (b);\
	(c)[2] = (a)[2] - (b);}

	template <class TYPE>
	inline void VectorSubtractConst3(TYPE *a, TYPE b, TYPE *c)
	{
		(c)[0] = (a)[0] - (b);
		(c)[1] = (a)[1] - (b);
		(c)[2] = (a)[2] - (b);
	}
#define VectorSubtractConst3d VectorSubtractConst3<double>
#define VectorSubtractConst3f VectorSubtractConst3<float>
#define VectorSubtractConst3i VectorSubtractConst3<int>


	// a -= b
#define VECTOR_SUBTRACT_CONST_EQUAL3(a,b)\
	{(a)[0] -= (b);\
	(a)[1] -= (b);\
	(a)[2] -= (b);}

	template <class TYPE>
	inline void VectorSubtractConstEqual3(TYPE *a, TYPE *b)
	{
		(a)[0] -= (b);
		(a)[1] -= (b);
		(a)[2] -= (b);
	}
#define VectorSubtractConstEqual3d VectorSubtractConstEqual3<double>
#define VectorSubtractConstEqual3f VectorSubtractConstEqual3<float>
#define VectorSubtractConstEqual3i VectorSubtractConstEqual3<int>


	// square of length
#define VECTOR_LENGTH_SQR3(a) ((a)[0]*(a)[0]+(a)[1]*(a)[1]+(a)[2]*(a)[2])
#define VECTOR_LENGTH_SQR3S(a, b, c) ((a)*(a)+(b)*(b)+(c)*(c))
	template <class TYPE>
	inline TYPE VectorLengthSqr3(TYPE *a)
	{
		return ((a)[0] * (a)[0] + (a)[1] * (a)[1] + (a)[2] * (a)[2]);
	}
	template <class TYPE>
	inline TYPE VectorLengthSqr3S(TYPE a, TYPE b, TYPE c)
	{
		return ((a)*(a)+(b)*(b)+(c)*(c));
	}
#define VectorLengthSqr3d VectorLengthSqr3<double>
#define VectorLengthSqr3f VectorLengthSqr3<float>
#define VectorLengthSqr3i VectorLengthSqr3<int>
#define VectorLengthSqr3Sd VectorLengthSqr3S<double>
#define VectorLengthSqr3Sf VectorLengthSqr3S<float>
#define VectorLengthSqr3Si VectorLengthSqr3S<int>

	// length
#define VECTOR_LENGTH3(a)\
	sqrt((a)[0]*(a)[0]+(a)[1]*(a)[1]+(a)[2]*(a)[2])
#define VECTOR_LENGTH3S(a, b, c)\
	sqrt((a)*(a)+(b)*(b)+(c)*(c))

	template <class TYPE>
	inline TYPE VectorLength3(TYPE *a)
	{
		return sqrt((a)[0] * (a)[0] + (a)[1] * (a)[1] + (a)[2] * (a)[2]);
	}
	template <class TYPE>
	inline TYPE VectorLength3S(TYPE a, TYPE b, TYPE c)
	{
		return sqrt((a)*(a)+(b)*(b)+(c)*(c));
	}
#define VectorLength3d VectorLength3<double>
#define VectorLength3f VectorLength3<float>
#define VectorLength3Sd VectorLength3S<double>
#define VectorLength3Sf VectorLength3S<float>


	// a *= scalar
#define VECTOR_SCALE3(a,scalar)\
	{(a)[0] *= (scalar);\
	(a)[1] *= (scalar);\
	(a)[2] *= (scalar);}
	// a *= scalar
#define VECTOR_SCALE3S(a,scalar1,scalar2,scalar3)\
	{(a)[0] *= (scalar1);\
	(a)[1] *= (scalar2);\
	(a)[2] *= (scalar3);}

	template <class TYPE>
	inline void VectorScale3(TYPE *a, TYPE scalar)
	{
		(a)[0] *= (scalar);
		(a)[1] *= (scalar);
		(a)[2] *= (scalar);
	}
	template <class TYPE>
	inline void VectorScale3S(TYPE *a, TYPE scalar1, TYPE scalar2, TYPE scalar3)
	{
		(a)[0] *= (scalar1);
		(a)[1] *= (scalar2);
		(a)[2] *= (scalar3);
	}
#define VectorScale3d VectorScale3<double>
#define VectorScale3f VectorScale3<float>
#define VectorScale3i VectorScale3<int>
#define VectorScale3Sd VectorScale3S<double>
#define VectorScale3Sf VectorScale3S<float>
#define VectorScale3Si VectorScale3S<int>


	// c = a DOT b
#define VECTOR_DOT_PRODUCT3(a,b)	( (a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2] )
	template <class TYPE>
	inline TYPE VectorDotProduct3(TYPE *a, TYPE *b)
	{
		return ((a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2]);
	}
#define VectorDotProduct3d VectorDotProduct3<double>
#define VectorDotProduct3f VectorDotProduct3<float>
#define VectorDotProduct3i VectorDotProduct3<int>


	// c = a x b
#define VECTOR_CROSS_PRODUCT3(a,b,c)\
	{(c)[0] = ((a)[1]) * ((b)[2]) - ((b)[1]) * ((a)[2]);\
	(c)[1] = ((b)[0]) * ((a)[2]) - ((a)[0]) * ((b)[2]);\
	(c)[2] = ((a)[0]) * ((b)[1]) - ((b)[0]) * ((a)[1]);}

	template <class TYPE>
	inline void VectorCrossProduct3(TYPE * a, TYPE * b, TYPE * c)
	{
		c[0] = a[1] * b[2] - a[2] * b[1];
		c[1] = a[2] * b[0] - a[0] * b[2];
		c[2] = a[0] * b[1] - a[1] * b[0];
	}
#define VectorCrossProduct3d VectorCrossProduct3<double>
#define VectorCrossProduct3f VectorCrossProduct3<float>
#define VectorCrossProduct3i VectorCrossProduct3<int>


	//	rotate around x axis
#define VECTOR_ROTATE_X3(rad, A, a)\
	{(a)[0] = ((A)[0]);\
	(a)[1] = cos(rad)*((A)[1])-sin(rad)*((A)[2]);\
	(a)[2] = sin(rad)*((A)[1])+cos(rad)*((A)[2]);}

	template <class TYPE>
	inline void VectorRotateX3(TYPE rad, TYPE * A, TYPE * a)
	{
		(a)[0] = ((A)[0]);
		(a)[1] = cos(rad)*((A)[1]) - sin(rad)*((A)[2]);
		(a)[2] = sin(rad)*((A)[1]) + cos(rad)*((A)[2]);
	}
#define VectorRotateX3d VectorRotateX3<double>
#define VectorRotateX3f VectorRotateX3<float>


	//	rotate around y axis
#define VECTOR_ROTATE_Y3(rad, A, a)\
	{(a)[0] = cos(rad)*((A)[0])+sin(rad)*((A)[2]);\
	(a)[1] = ((A)[1]);\
	(a)[2] = -sin(rad)*((A)[0])+cos(rad)*((A)[2]);}

	template <class TYPE>
	inline void VectorRotateY3(TYPE rad, TYPE * A, TYPE * a)
	{
		(a)[0] = cos(rad)*((A)[0]) + sin(rad)*((A)[2]);
		(a)[1] = ((A)[1]);
		(a)[2] = -sin(rad)*((A)[0]) + cos(rad)*((A)[2]);
	}
#define VectorRotateY3d VectorRotateY3<double>
#define VectorRotateY3f VectorRotateY3<float>


	//	rotate around z axis
#define VECTOR_ROTATE_Z3(rad, A, a)\
	{(a)[0] = cos(rad)*((A)[0])-sin(rad)*((A)[1]);\
	(a)[1] = sin(rad)*((A)[0])+cos(rad)*((A)[1]);\
	(a)[2] = ((A)[2]);}

	template <class TYPE>
	inline void VectorRotateZ3(TYPE rad, TYPE * A, TYPE * a)
	{
		(a)[0] = cos(rad)*((A)[0]) - sin(rad)*((A)[1]);
		(a)[1] = sin(rad)*((A)[0]) + cos(rad)*((A)[1]);
		(a)[2] = ((A)[2]);
	}
#define VectorRotateZ3d VectorRotateZ3<double>
#define VectorRotateZ3f VectorRotateZ3<float>


	// normalize a vector and return length, must input non-zero vector
#define VECTOR_NORMALIZE3(v, u, len)\
	{len=(v)[0]*(v)[0]+(v)[1]*(v)[1]+(v)[2]*(v)[2];\
	len=sqrt(len);\
	(u)[0]=(v)[0]/len;(u)[1]=(v)[1]/len;(u)[2]=(v)[2]/len;}\

	template <class TYPE>
	inline TYPE VectorNormalize3(TYPE *v, TYPE *u)
	{
		TYPE len = (v)[0] * (v)[0] + (v)[1] * (v)[1] + (v)[2] * (v)[2];
		len = sqrt(len);
		(u)[0] = (v)[0] / len; (u)[1] = (v)[1] / len; (u)[2] = (v)[2] / len;
		return len;
	}
	template <class TYPE>
	inline TYPE VectorNormalize3(TYPE *v)
	{
		TYPE len = (v)[0] * (v)[0] + (v)[1] * (v)[1] + (v)[2] * (v)[2];
		len = sqrt(len);
		(v)[0] = (v)[0] / len; (v)[1] = (v)[1] / len; (v)[2] = (v)[2] / len;
		return len;
	}
#define VectorNormalize3d VectorNormalize3<double>
#define VectorNormalize3f VectorNormalize3<float>


#define VECTOR_PRINT3(v)	std::cout << '[' << (v)[0] << ', ' << (v)[1] << ', ' << (v)[2] << ']\n'
	template <class TYPE>
	inline void VectorPrint3(TYPE *v)
	{
		std::cout << '[' << (v)[0] << ', ' << (v)[1] << ', ' << (v)[2] << ']\n';
	}
#define VectorPrint3d VectorPrint3<double>
#define VectorPrint3f VectorPrint3<float>



#define VECTOR_DISTANCE3(p1, p2)	sqrt(((p1)[0]-(p2)[0])*((p1)[0]-(p2)[0])+((p1)[1]-(p2)[1])*((p1)[1]-(p2)[1])+((p1)[2]-(p2)[2])*((p1)[2]-(p2)[2]))
	template <class TYPE>
	inline TYPE VectorDistance3(TYPE *p1, TYPE *p2)
	{
		return sqrt(((p1)[0] - (p2)[0])*((p1)[0] - (p2)[0]) + ((p1)[1] - (p2)[1])*((p1)[1] - (p2)[1]) + ((p1)[2] - (p2)[2])*((p1)[2] - (p2)[2]));
	}
#define VectorDistance3d VectorDistance3<double>
#define VectorDistance3f VectorDistance3<float>


#define VECTOR_DISTANCE_SQR3(p1, p2)	(((p1)[0]-(p2)[0])*((p1)[0]-(p2)[0])+((p1)[1]-(p2)[1])*((p1)[1]-(p2)[1])+((p1)[2]-(p2)[2])*((p1)[2]-(p2)[2]))
	template <class TYPE>
	inline TYPE VectorDistanceSqr3(TYPE *p1, TYPE *p2)
	{
		return (((p1)[0] - (p2)[0])*((p1)[0] - (p2)[0]) + ((p1)[1] - (p2)[1])*((p1)[1] - (p2)[1]) + ((p1)[2] - (p2)[2])*((p1)[2] - (p2)[2]));
	}
#define VectorDistanceSqr3d VectorDistanceSqr3<double>
#define VectorDistanceSqr3f VectorDistanceSqr3<float>


#define VECTOR_MIN_INDEX3(v)	((v)[0]>(v)[1] ? ((v)[1]>(v)[2] ? 2 : 1) : ((v)[0]>(v)[2] ? 2 : 0))
	template <class TYPE>
	inline int VectorMinIndex3(TYPE *v)
	{
		return ((v)[0] > (v)[1] ? ((v)[1] > (v)[2] ? 2 : 1) : ((v)[0] > (v)[2] ? 2 : 0));
	}
#define VectorMinIndex3d VectorMinIndex3<double>
#define VectorMinIndex3f VectorMinIndex3<float>
#define VectorMinIndex3i VectorMinIndex3<int>


#define VECTOR_MAX_INDEX3(v)	((v)[0]<(v)[1] ? ((v)[1]<(v)[2] ? 2 : 1) : ((v)[0]<(v)[2] ? 2 : 0))
	template <class TYPE>
	inline int VectorMaxIndex3(TYPE *v)
	{
		return ((v)[0] < (v)[1] ? ((v)[1] < (v)[2] ? 2 : 1) : ((v)[0] < (v)[2] ? 2 : 0));
	}
#define VectorMaxIndex3d VectorMaxIndex3<double>
#define VectorMaxIndex3f VectorMaxIndex3<float>
#define VectorMaxIndex3i VectorMaxIndex3<int>



#define VECTOR_MAX_INDEX3S(a, b, c) ((a)<(b) ? ((b)<(c) ? 2 : 1) : ((a)<(c) ? 2 : 0))
#define VECTOR_MIN_INDEX3S(a, b, c) ((a)>(b) ? ((b)>(c) ? 2 : 1) : ((a)>(c) ? 2 : 0))

	template <class TYPE>
	inline void VectorFindOrthoNormal3(TYPE *v, TYPE *result)
	{
		// find smallest abs component of v
		int smallestIndex = 0;
		for (int dim = 1; dim<3; dim++)
			if (fabs(v[dim]) < fabs(v[smallestIndex]))
				smallestIndex = dim;

		TYPE axis[3] = { 0.0, 0.0, 0.0 };
		axis[smallestIndex] = 1.0;

		// this cross-product will be non-zero (as long as v is not zero)
		VECTOR_CROSS_PRODUCT3(v, axis, result);
		TYPE len;
		VECTOR_NORMALIZE3(result, result, len);
	}
#define VectorFindOrthoNormal3d VectorFindOrthoNormal3<double>
#define VectorFindOrthoNormal3f VectorFindOrthoNormal3<float>


template <class TYPE>
class Vec3
{
public:
	inline Vec3() { v[0] = 0; v[1] = 0; v[2] = 0; }
	inline Vec3(TYPE x, TYPE y, TYPE z) { v[0] = x; v[1] = y; v[2] = z; }
	inline Vec3(TYPE entry){ v[0] = entry; v[1] = entry; v[2] = entry; }
	inline Vec3(TYPE * vec){ v[0] = vec[0]; v[1] = vec[1]; v[2] = vec[2]; }
	inline Vec3(Vec3 & vec){ v[0] = vec.v[0]; v[1] = vec.v[1]; v[2] = vec.v[2]; }

	inline void set(TYPE x, TYPE y, TYPE z){ v[0] = x; v[1] = y; v[2] = z; }
	inline void set(TYPE entry){ v[0] = entry; v[1] = entry; v[2] = entry; }

	inline Vec3 & operator=(Vec3 & vec){ v[0] = vec.v[0]; v[1] = vec.v[1]; v[2] = vec.v[2]; }
	inline TYPE & operator[] (int index) // M[i] returns i-th row
	{
		return v[index];
	}

public:
	TYPE v[3];
};
typedef Vec3<int> Vec3i;
typedef Vec3<float> Vec3f;
typedef Vec3<double> Vec3d;


// get the normal of a 3d triangle, return area of the triangle
template <class TYPE>
inline TYPE TriangleNormal3(TYPE *v1, TYPE *v2, TYPE *v3, TYPE *norm)
{
	TYPE d12[3], d13[3];
	VECTOR_SUBTRACT3(v2, v1, d12);
	VECTOR_SUBTRACT3(v3, v1, d13);
	VECTOR_CROSS_PRODUCT3(d12, d13, norm);
	TYPE len;
	VECTOR_NORMALIZE3(norm, norm, len);
	return len;
}
#define TriangleNormal3d TriangleNormal3<double>
#define TriangleNormal3f TriangleNormal3<float>

// get the center of a 3d triangle,
template <class TYPE>
inline void TriangleCenter3(TYPE *v1, TYPE *v2, TYPE *v3, TYPE *center)
{
	center[0] = (v1[0] + v2[0] + v3[0]) / 3.0;
	center[1] = (v1[1] + v2[1] + v3[1]) / 3.0;
	center[2] = (v1[2] + v2[2] + v3[2]) / 3.0;
}
#define TriangleCenter3d TriangleCenter3<double>
#define TriangleCenter3f TriangleCenter3<float>

}

#endif