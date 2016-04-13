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

3D matrix operations
*/


#ifndef __PRINTX_MAT3_H__
#define __PRINTX_MAT3_H__

#include "vec3.h"

namespace PRINTX_LIB
{
#define MATRIX_SET_IDENTITY3X3(a)\
	{	(a)[0] = 1;(a)[1] = 0;(a)[2] = 0;\
	(a)[3] = 0;(a)[4] = 1;(a)[5] = 0;\
	(a)[6] = 0;(a)[7] = 0;(a)[8] = 1;}

	template <class TYPE>
	inline void MatrixSetIdentity3X3(TYPE *a)
	{
		(a)[0] = 1; (a)[1] = 0; (a)[2] = 0;
		(a)[3] = 0; (a)[4] = 1; (a)[5] = 0;
		(a)[6] = 0; (a)[7] = 0; (a)[8] = 1;
	}
#define MatrixSetIdentity3X3d MatrixSetIdentity3X3<double>
#define MatrixSetIdentity3X3f MatrixSetIdentity3X3<float>


#define MATRIX_SET_DIAG3X3(a, diag)\
	{	(a)[0] = diag;(a)[1] = 0;(a)[2] = 0;\
	(a)[3] = 0;(a)[4] = diag;(a)[5] = 0;\
	(a)[6] = 0;(a)[7] = 0;(a)[8] = diag;}

	template <class TYPE>
	inline void MatrixSetDiag3X3(TYPE *a, TYPE diag)
	{
		(a)[0] = diag; (a)[1] = 0; (a)[2] = 0;
		(a)[3] = 0; (a)[4] = diag; (a)[5] = 0;
		(a)[6] = 0; (a)[7] = 0; (a)[8] = diag;
	}
#define MatrixSetDiag3X3d MatrixSetDiag3X3<double>
#define MatrixSetDiag3X3f MatrixSetDiag3X3<float>



#define MATRIX_FROM_2D3X3(a, b)\
	{	(a)[0] = b[0][0];(a)[1] = b[0][1];(a)[2] = b[0][2];\
	(a)[3] = b[1][0];(a)[4] = b[1][1];(a)[5] = b[1][2];\
	(a)[6] = b[2][0];(a)[7] = b[2][1];(a)[8] = b[2][2];}



#define MATRIX_SET_CONST3X3(a, c)\
	{	(a)[0] = c;(a)[1] = c;(a)[2] = c;\
	(a)[3] = c;(a)[4] = c;(a)[5] = c;\
	(a)[6] = c;(a)[7] = c;(a)[8] = c;}

	template <class TYPE>
	inline void MatrixSetConst3X3(TYPE *a, TYPE c)
	{
		(a)[0] = c; (a)[1] = c; (a)[2] = c;
		(a)[3] = c; (a)[4] = c; (a)[5] = c;
		(a)[6] = c; (a)[7] = c; (a)[8] = c;
	}
#define MatrixSetConst3X3d MatrixSetConst3X3<double>
#define MatrixSetConst3X3f MatrixSetConst3X3<float>


	// a = b
#define MATRIX_SET3X3(a,b)\
	{(a)[0] = (b)[0]; (a)[1] = (b)[1]; (a)[2] = (b)[2];\
	(a)[3] = (b)[3]; (a)[4] = (b)[4]; (a)[5] = (b)[5];\
	(a)[6] = (b)[6]; (a)[7] = (b)[7]; (a)[8] = (b)[8];}

	template <class TYPE>
	inline void MatrixSet3X3(TYPE *a, TYPE *b)
	{
		(a)[0] = (b)[0]; (a)[1] = (b)[1]; (a)[2] = (b)[2];
		(a)[3] = (b)[3]; (a)[4] = (b)[4]; (a)[5] = (b)[5];
		(a)[6] = (b)[6]; (a)[7] = (b)[7]; (a)[8] = (b)[8];
	}
#define MatrixSet3X3d MatrixSet3X3<double>
#define MatrixSet3X3f MatrixSet3X3<float>


	// a *= scalar
#define MATRIX_SCALE3X3(a,scalar)\
	{  (a)[0] *= (scalar); (a)[1] *= (scalar); (a)[2] *= (scalar);\
	(a)[3] *= (scalar); (a)[4] *= (scalar); (a)[5] *= (scalar);\
	(a)[6] *= (scalar); (a)[7] *= (scalar); (a)[8] *= (scalar);}

	template <class TYPE>
	inline void MatrixScale3X3(TYPE *a, TYPE scalar)
	{
		(a)[0] *= (scalar); (a)[1] *= (scalar); (a)[2] *= (scalar);
		(a)[3] *= (scalar); (a)[4] *= (scalar); (a)[5] *= (scalar);
		(a)[6] *= (scalar); (a)[7] *= (scalar); (a)[8] *= (scalar);
	}
#define MatrixScale3X3d MatrixScale3X3<double>
#define MatrixScale3X3f MatrixScale3X3<float>


	// c = a + b
#define MATRIX_ADD3X3(a,b,c)\
	{(c)[0] = (a)[0] + (b)[0]; (c)[1] = (a)[1] + (b)[1]; (c)[2] = (a)[2] + (b)[2];\
	(c)[3] = (a)[3] + (b)[3]; (c)[4] = (a)[4] + (b)[4]; (c)[5] = (a)[5] + (b)[5];\
	(c)[6] = (a)[6] + (b)[6]; (c)[7] = (a)[7] + (b)[7]; (c)[8] = (a)[8] + (b)[8];}

	template <class TYPE>
	inline void MatrixAdd3X3(TYPE *a, TYPE *b, TYPE *c)
	{
		(c)[0] = (a)[0] + (b)[0]; (c)[1] = (a)[1] + (b)[1]; (c)[2] = (a)[2] + (b)[2];
		(c)[3] = (a)[3] + (b)[3]; (c)[4] = (a)[4] + (b)[4]; (c)[5] = (a)[5] + (b)[5];
		(c)[6] = (a)[6] + (b)[6]; (c)[7] = (a)[7] + (b)[7]; (c)[8] = (a)[8] + (b)[8];
	}
#define MatrixAdd3X3d MatrixAdd3X3<double>
#define MatrixAdd3X3f MatrixAdd3X3<float>

	// c = a * b
#define MATRIX_MULTIPLY3X3(a,b,c)\
	{(c)[0]=(a)[0]*(b)[0]+(a)[1]*(b)[3]+(a)[2]*(b)[6];\
	(c)[1]=(a)[0]*(b)[1]+(a)[1]*(b)[4]+(a)[2]*(b)[7];\
	(c)[2]=(a)[0]*(b)[2]+(a)[1]*(b)[5]+(a)[2]*(b)[8];\
	(c)[3]=(a)[3]*(b)[0]+(a)[4]*(b)[3]+(a)[5]*(b)[6];\
	(c)[4]=(a)[3]*(b)[1]+(a)[4]*(b)[4]+(a)[5]*(b)[7];\
	(c)[5]=(a)[3]*(b)[2]+(a)[4]*(b)[5]+(a)[5]*(b)[8];\
	(c)[6]=(a)[6]*(b)[0]+(a)[7]*(b)[3]+(a)[8]*(b)[6];\
	(c)[7]=(a)[6]*(b)[1]+(a)[7]*(b)[4]+(a)[8]*(b)[7];\
	(c)[8]=(a)[6]*(b)[2]+(a)[7]*(b)[5]+(a)[8]*(b)[8];}

	template <class TYPE>
	inline void MatrixMultiply3X3(TYPE *a, TYPE *b, TYPE *c)
	{
		(c)[0] = (a)[0] * (b)[0] + (a)[1] * (b)[3] + (a)[2] * (b)[6];
		(c)[1] = (a)[0] * (b)[1] + (a)[1] * (b)[4] + (a)[2] * (b)[7];
		(c)[2] = (a)[0] * (b)[2] + (a)[1] * (b)[5] + (a)[2] * (b)[8];
		(c)[3] = (a)[3] * (b)[0] + (a)[4] * (b)[3] + (a)[5] * (b)[6];
		(c)[4] = (a)[3] * (b)[1] + (a)[4] * (b)[4] + (a)[5] * (b)[7];
		(c)[5] = (a)[3] * (b)[2] + (a)[4] * (b)[5] + (a)[5] * (b)[8];
		(c)[6] = (a)[6] * (b)[0] + (a)[7] * (b)[3] + (a)[8] * (b)[6];
		(c)[7] = (a)[6] * (b)[1] + (a)[7] * (b)[4] + (a)[8] * (b)[7];
		(c)[8] = (a)[6] * (b)[2] + (a)[7] * (b)[5] + (a)[8] * (b)[8];
	}
#define MatrixMultiply3X3d MatrixMultiply3X3<double>
#define MatrixMultiply3X3f MatrixMultiply3X3<float>


	// c = a * diag(b), can be self scaled
#define MATRIX_MULTIPLY_DIAG_RIGHT3X3(a,b,c)\
	{(c)[0]=(a)[0]*(b)[0];(c)[3]=(a)[3]*(b)[0];(c)[6]=(a)[6]*(b)[0];\
	(c)[1]=(a)[1]*(b)[1];(c)[4]=(a)[4]*(b)[1];(c)[7]=(a)[7]*(b)[1];\
	(c)[2]=(a)[2]*(b)[2];(c)[5]=(a)[5]*(b)[2];(c)[8]=(a)[8]*(b)[2];}

	template <class TYPE>
	inline void MatrixMultiplyDiagRight3X3(TYPE *a, TYPE *b, TYPE *c)
	{
		(c)[0] = (a)[0] * (b)[0]; (c)[3] = (a)[3] * (b)[0]; (c)[6] = (a)[6] * (b)[0];
		(c)[1] = (a)[1] * (b)[1]; (c)[4] = (a)[4] * (b)[1]; (c)[7] = (a)[7] * (b)[1];
		(c)[2] = (a)[2] * (b)[2]; (c)[5] = (a)[5] * (b)[2]; (c)[8] = (a)[8] * (b)[2];
	}
#define MatrixMultiplyDiagRight3X3d MatrixMultiplyDiagRight3X3<double>
#define MatrixMultiplyDiagRight3X3f MatrixMultiplyDiagRight3X3<float>


	// c = diag(b) * a, can be self scaled
#define MATRIX_MULTIPLY_DIAG_LEFT3X3(a,b,c)\
	{(c)[0]=(a)[0]*(b)[0];(c)[1]=(a)[1]*(b)[0];(c)[2]=(a)[2]*(b)[0];\
	(c)[3]=(a)[3]*(b)[1];(c)[4]=(a)[4]*(b)[1];(c)[5]=(a)[5]*(b)[1];\
	(c)[6]=(a)[6]*(b)[2];(c)[7]=(a)[7]*(b)[2];(c)[8]=(a)[8]*(b)[2];}

	template <class TYPE>
	inline void MatrixMultiplyDiagLeft3X3(TYPE *a, TYPE *b, TYPE *c)
	{
		(c)[0] = (a)[0] * (b)[0]; (c)[1] = (a)[1] * (b)[0]; (c)[2] = (a)[2] * (b)[0];
		(c)[3] = (a)[3] * (b)[1]; (c)[4] = (a)[4] * (b)[1]; (c)[5] = (a)[5] * (b)[1];
		(c)[6] = (a)[6] * (b)[2]; (c)[7] = (a)[7] * (b)[2]; (c)[8] = (a)[8] * (b)[2];
	}
#define MatrixMultiplyDiagLeft3X3d MatrixMultiplyDiagLeft3X3<double>
#define MatrixMultiplyDiagLeft3X3f MatrixMultiplyDiagLeft3X3<float>


	//c = a * bT
#define MATRIX_MULTIPLY3X3ABT(a,b,c)\
	{(c)[0]=(a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2];\
	(c)[1]=(a)[0]*(b)[3]+(a)[1]*(b)[4]+(a)[2]*(b)[5];\
	(c)[2]=(a)[0]*(b)[6]+(a)[1]*(b)[7]+(a)[2]*(b)[8];\
	(c)[3]=(a)[3]*(b)[0]+(a)[4]*(b)[1]+(a)[5]*(b)[2];\
	(c)[4]=(a)[3]*(b)[3]+(a)[4]*(b)[4]+(a)[5]*(b)[5];\
	(c)[5]=(a)[3]*(b)[6]+(a)[4]*(b)[7]+(a)[5]*(b)[8];\
	(c)[6]=(a)[6]*(b)[0]+(a)[7]*(b)[1]+(a)[8]*(b)[2];\
	(c)[7]=(a)[6]*(b)[3]+(a)[7]*(b)[4]+(a)[8]*(b)[5];\
	(c)[8]=(a)[6]*(b)[6]+(a)[7]*(b)[7]+(a)[8]*(b)[8];}

	template <class TYPE>
	inline void MatrixMultiply3X3ABT(TYPE *a, TYPE *b, TYPE *c)
	{
		(c)[0] = (a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2];
		(c)[1] = (a)[0] * (b)[3] + (a)[1] * (b)[4] + (a)[2] * (b)[5];
		(c)[2] = (a)[0] * (b)[6] + (a)[1] * (b)[7] + (a)[2] * (b)[8];
		(c)[3] = (a)[3] * (b)[0] + (a)[4] * (b)[1] + (a)[5] * (b)[2];
		(c)[4] = (a)[3] * (b)[3] + (a)[4] * (b)[4] + (a)[5] * (b)[5];
		(c)[5] = (a)[3] * (b)[6] + (a)[4] * (b)[7] + (a)[5] * (b)[8];
		(c)[6] = (a)[6] * (b)[0] + (a)[7] * (b)[1] + (a)[8] * (b)[2];
		(c)[7] = (a)[6] * (b)[3] + (a)[7] * (b)[4] + (a)[8] * (b)[5];
		(c)[8] = (a)[6] * (b)[6] + (a)[7] * (b)[7] + (a)[8] * (b)[8];
	}
#define MatrixMultiply3X3ABTd MatrixMultiply3X3ABT<double>
#define MatrixMultiply3X3ABTf MatrixMultiply3X3ABT<float>


	//c = aT * b
#define MATRIX_MULTIPLY3X3ATB(a,b,c)\
	{(c)[0]=(a)[0]*(b)[0]+(a)[3]*(b)[3]+(a)[6]*(b)[6];\
	(c)[1]=(a)[0]*(b)[1]+(a)[3]*(b)[4]+(a)[6]*(b)[7];\
	(c)[2]=(a)[0]*(b)[2]+(a)[3]*(b)[5]+(a)[6]*(b)[8];\
	(c)[3]=(a)[1]*(b)[0]+(a)[4]*(b)[3]+(a)[7]*(b)[6];\
	(c)[4]=(a)[1]*(b)[1]+(a)[4]*(b)[4]+(a)[7]*(b)[7];\
	(c)[5]=(a)[1]*(b)[2]+(a)[4]*(b)[5]+(a)[7]*(b)[8];\
	(c)[6]=(a)[2]*(b)[0]+(a)[5]*(b)[3]+(a)[8]*(b)[6];\
	(c)[7]=(a)[2]*(b)[1]+(a)[5]*(b)[4]+(a)[8]*(b)[7];\
	(c)[8]=(a)[2]*(b)[2]+(a)[5]*(b)[5]+(a)[8]*(b)[8];}

	template <class TYPE>
	inline void MatrixMultiply3X3ATB(TYPE *a, TYPE *b, TYPE *c)
	{
		(c)[0] = (a)[0] * (b)[0] + (a)[3] * (b)[3] + (a)[6] * (b)[6];
		(c)[1] = (a)[0] * (b)[1] + (a)[3] * (b)[4] + (a)[6] * (b)[7];
		(c)[2] = (a)[0] * (b)[2] + (a)[3] * (b)[5] + (a)[6] * (b)[8];
		(c)[3] = (a)[1] * (b)[0] + (a)[4] * (b)[3] + (a)[7] * (b)[6];
		(c)[4] = (a)[1] * (b)[1] + (a)[4] * (b)[4] + (a)[7] * (b)[7];
		(c)[5] = (a)[1] * (b)[2] + (a)[4] * (b)[5] + (a)[7] * (b)[8];
		(c)[6] = (a)[2] * (b)[0] + (a)[5] * (b)[3] + (a)[8] * (b)[6];
		(c)[7] = (a)[2] * (b)[1] + (a)[5] * (b)[4] + (a)[8] * (b)[7];
		(c)[8] = (a)[2] * (b)[2] + (a)[5] * (b)[5] + (a)[8] * (b)[8];
	}
#define MatrixMultiply3X3ATBd MatrixMultiply3X3ATB<double>
#define MatrixMultiply3X3ATBf MatrixMultiply3X3ATB<float>


	// w = A * v
#define MATRIX_VECTOR_MULTIPLY3X3(A,v,w)\
	{(w)[0] = (A)[0]*(v)[0]+(A)[1]*(v)[1]+(A)[2]*(v)[2];\
	(w)[1] = (A)[3]*(v)[0]+(A)[4]*(v)[1]+(A)[5]*(v)[2];\
	(w)[2] = (A)[6]*(v)[0]+(A)[7]*(v)[1]+(A)[8]*(v)[2];}

	template <class TYPE>
	inline void MatrixVectorMultiply3X3(TYPE *A, TYPE *v, TYPE *w)
	{
		(w)[0] = (A)[0] * (v)[0] + (A)[1] * (v)[1] + (A)[2] * (v)[2];
		(w)[1] = (A)[3] * (v)[0] + (A)[4] * (v)[1] + (A)[5] * (v)[2];
		(w)[2] = (A)[6] * (v)[0] + (A)[7] * (v)[1] + (A)[8] * (v)[2];
	}
#define MatrixVectorMultiply3X3d MatrixVectorMultiply3X3<double>
#define MatrixVectorMultiply3X3f MatrixVectorMultiply3X3<float>



	// w = A^T * v
#define MATRIX_VECTOR_MULTIPLY3X3T(A,v,w)\
	{(w)[0] = (A)[0]*(v)[0]+(A)[3]*(v)[1]+(A)[6]*(v)[2];\
	(w)[1] = (A)[1]*(v)[0]+(A)[4]*(v)[1]+(A)[7]*(v)[2];\
	(w)[2] = (A)[2]*(v)[0]+(A)[5]*(v)[1]+(A)[8]*(v)[2];}

	template <class TYPE>
	inline void MatrixVectorMultiply3X3T(TYPE *A, TYPE *v, TYPE *w)
	{
		(w)[0] = (A)[0] * (v)[0] + (A)[3] * (v)[1] + (A)[6] * (v)[2];
		(w)[1] = (A)[1] * (v)[0] + (A)[4] * (v)[1] + (A)[7] * (v)[2];
		(w)[2] = (A)[2] * (v)[0] + (A)[5] * (v)[1] + (A)[8] * (v)[2];
	}
#define MatrixVectorMultiply3X3Td MatrixVectorMultiply3X3T<double>
#define MatrixVectorMultiply3X3Tf MatrixVectorMultiply3X3T<float>


// B = A^T
#define MATRIX_TRANSPOSE3X3(A,B)\
	{(B)[0] = (A)[0]; (B)[1] = (A)[3]; (B)[2] = (A)[6];\
	(B)[3] = (A)[1]; (B)[4] = (A)[4]; (B)[5] = (A)[7];\
	(B)[6] = (A)[2]; (B)[7] = (A)[5]; (B)[8] = (A)[8];}

	template <class TYPE>
	inline void MatrixTranspose3X3(TYPE *A, TYPE *B)
	{
		(B)[0] = (A)[0]; (B)[1] = (A)[3]; (B)[2] = (A)[6];
		(B)[3] = (A)[1]; (B)[4] = (A)[4]; (B)[5] = (A)[7];
		(B)[6] = (A)[2]; (B)[7] = (A)[5]; (B)[8] = (A)[8];
	}
#define MatrixTranspose3X3d MatrixTranspose3X3<double>
#define MatrixTranspose3X3f MatrixTranspose3X3<float>



// A = A^T
#define MATRIX_SELF_TRANSPOSE3X3(A) {std::swap((A)[1], (A)[3]); std::swap((A)[2], (A)[6]);std::swap((A)[5], (A)[7]);}
	template <class TYPE>
	inline void MatrixSelfTranspose3X3(TYPE *A)
	{
		std::swap((A)[1], (A)[3]); std::swap((A)[2], (A)[6]); std::swap((A)[5], (A)[7]);
	}
#define MatrixSelfTranspose3X3d MatrixSelfTranspose3X3<double>
#define MatrixSelfTranspose3X3f MatrixSelfTranspose3X3<float>



// A = u * v^T
#define VECTOR_TENSOR_PRODUCT3X3(u,v,A)\
	{(A)[0] = (u)[0] * (v)[0]; (A)[1] = (u)[0] * (v)[1]; (A)[2] = (u)[0] * (v)[2];\
	(A)[3] = (u)[1] * (v)[0]; (A)[4] = (u)[1] * (v)[1]; (A)[5] = (u)[1] * (v)[2];\
	(A)[6] = (u)[2] * (v)[0]; (A)[7] = (u)[2] * (v)[1]; (A)[8] = (u)[2] * (v)[2];}

	template <class TYPE>
	inline void MatrixTensorProduct3X3(TYPE *u, TYPE *v, TYPE *A)
	{
		(A)[0] = (u)[0] * (v)[0]; (A)[1] = (u)[0] * (v)[1]; (A)[2] = (u)[0] * (v)[2];
		(A)[3] = (u)[1] * (v)[0]; (A)[4] = (u)[1] * (v)[1]; (A)[5] = (u)[1] * (v)[2];
		(A)[6] = (u)[2] * (v)[0]; (A)[7] = (u)[2] * (v)[1]; (A)[8] = (u)[2] * (v)[2];
	}
#define MatrixTensorProduct3X3d MatrixTensorProduct3X3<double>
#define MatrixTensorProduct3X3f MatrixTensorProduct3X3<float>


// A += u * v^T
#define VECTOR_TENSOR_PRODUCT_ADD3X3(u,v,A)\
	{(A)[0] += (u)[0] * (v)[0]; (A)[1] += (u)[0] * (v)[1]; (A)[2] += (u)[0] * (v)[2];\
	(A)[3] += (u)[1] * (v)[0]; (A)[4] += (u)[1] * (v)[1]; (A)[5] += (u)[1] * (v)[2];\
	(A)[6] += (u)[2] * (v)[0]; (A)[7] += (u)[2] * (v)[1]; (A)[8] += (u)[2] * (v)[2];}

	template <class TYPE>
	inline void VectorTensorProductAdd3X3(TYPE *u, TYPE *v, TYPE *A)
	{
		(A)[0] += (u)[0] * (v)[0]; (A)[1] += (u)[0] * (v)[1]; (A)[2] += (u)[0] * (v)[2];
		(A)[3] += (u)[1] * (v)[0]; (A)[4] += (u)[1] * (v)[1]; (A)[5] += (u)[1] * (v)[2];
		(A)[6] += (u)[2] * (v)[0]; (A)[7] += (u)[2] * (v)[1]; (A)[8] += (u)[2] * (v)[2];
	}
#define VectorTensorProductAdd3X3d VectorTensorProductAdd3X3<double>
#define VectorTensorProductAdd3X3f VectorTensorProductAdd3X3<float>


// A dot B
#define MATRIX_DOT_PRODUCT3X3(A,B)\
	{((A)[0] * (B)[0] + (A)[1] * (B)[1] + (A)[2] * (B)[2] +\
	(A)[3] * (B)[3] + (A)[4] * (B)[4] + (A)[5] * (B)[5] +\
	(A)[6] * (B)[6] + (A)[7] * (B)[7] + (A)[8] * (B)[8] )}

	template <class TYPE>
	inline TYPE MatrixDotProduct3X3(TYPE *A, TYPE *B)
	{
		return ((A)[0] * (B)[0] + (A)[1] * (B)[1] + (A)[2] * (B)[2] +
			(A)[3] * (B)[3] + (A)[4] * (B)[4] + (A)[5] * (B)[5] +
			(A)[6] * (B)[6] + (A)[7] * (B)[7] + (A)[8] * (B)[8]);
	}
#define MatrixDotProduct3X3d MatrixDotProduct3X3<double>
#define MatrixDotProduct3X3f MatrixDotProduct3X3<float>


// A=skew(a) so that AT=-A
#define SKEW_MATRIX3X3(a, A)\
	{(A)[0] =       0;  (A)[1] = -(a)[2]; (A)[2] =  (a)[1];\
	(A)[3] =  (a)[2];  (A)[4] =       0; (A)[5] = -(a)[0];\
	(A)[6] = -(a)[1];  (A)[7] =  (a)[0]; (A)[8] =       0;}

	template <class TYPE>
	inline void SkewMatrix3X3(TYPE *a, TYPE *A)
	{
		(A)[0] = 0;  (A)[1] = -(a)[2]; (A)[2] = (a)[1];
		(A)[3] = (a)[2];  (A)[4] = 0; (A)[5] = -(a)[0];
		(A)[6] = -(a)[1];  (A)[7] = (a)[0]; (A)[8] = 0;
	}
#define SkewMatrix3X3d SkewMatrix3X3<double>
#define SkewMatrix3X3f SkewMatrix3X3<float>


	// A=sym(a) so that AT=A
#define SYM_MATRIX3X3(a, A)\
	{(A)[0] =  (a)[0];  (A)[1] =  (a)[1]; (A)[2] =  (a)[2];\
	(A)[3] =  (a)[1];  (A)[4] =  (a)[3]; (A)[5] =  (a)[4];\
	(A)[6] =  (a)[2];  (A)[7] =  (a)[4]; (A)[8] =  (a)[5];}

	template <class TYPE>
	inline void SymMatrix3X3(TYPE *a, TYPE *A)
	{
		(A)[0] = (a)[0];  (A)[1] = (a)[1]; (A)[2] = (a)[2];
		(A)[3] = (a)[1];  (A)[4] = (a)[3]; (A)[5] = (a)[4];
		(A)[6] = (a)[2];  (A)[7] = (a)[4]; (A)[8] = (a)[5];
	}
#define SymMatrix3X3d SymMatrix3X3<double>
#define SymMatrix3X3f SymMatrix3X3<float>




	// det(A)
#define MATRIX_DETERMINANT3X3(A)\
	( (A)[0] * ( (A)[4] * (A)[8] - (A)[5] * (A)[7] ) +\
	(A)[1] * ( (A)[5] * (A)[6] - (A)[3] * (A)[8] ) +\
	(A)[2] * ( (A)[3] * (A)[7] - (A)[4] * (A)[6] ) )

	template <class TYPE>
	inline TYPE MatrixDeterminant3X3(TYPE *A)
	{
		return 	((A)[0] * ((A)[4] * (A)[8] - (A)[5] * (A)[7]) +
			(A)[1] * ((A)[5] * (A)[6] - (A)[3] * (A)[8]) +
			(A)[2] * ((A)[3] * (A)[7] - (A)[4] * (A)[6]));
	}
#define MatrixDeterminant3X3d MatrixDeterminant3X3<double>
#define MatrixDeterminant3X3f MatrixDeterminant3X3<float>



// a = skew( 0.5 * (A - A^T) )
#define MATRIX_SKEW_PART3X3(A, a)\
	{(a)[0] = 0.5 * ((A)[7] - (A)[5]);\
	(a)[1] = 0.5 * ((A)[2] - (A)[6]);\
	(a)[2] = 0.5 * ((A)[3] - (A)[1]);}

	template <class TYPE>
	inline void MatrixSkewPart3X3(TYPE *A, TYPE *a)
	{
		(a)[0] = 0.5 * ((A)[7] - (A)[5]);
		(a)[1] = 0.5 * ((A)[2] - (A)[6]);
		(a)[2] = 0.5 * ((A)[3] - (A)[1]);
	}
#define MatrixSkewPart3X3d MatrixSkewPart3X3<double>
#define MatrixSkewPart3X3f MatrixSkewPart3X3<float>


// a = upper-triangle( 0.5 * (A + A^T) )
#define MATRIX_SYM_PART3X3(A, a)\
	{(a)[0] = ((A)[0]);\
	(a)[1] = 0.5 * ((A)[1] + (A)[3]);\
	(a)[2] = 0.5 * ((A)[2] + (A)[6]);\
	(a)[3] = ((A)[4]);\
	(a)[4] = 0.5 * ((A)[5] + (A)[7]);\
	(a)[5] = ((A)[8]);}

	template <class TYPE>
	inline void MatrixSymPart3X3(TYPE *A, TYPE *a)
	{
		(a)[0] = ((A)[0]);
		(a)[1] = 0.5 * ((A)[1] + (A)[3]);
		(a)[2] = 0.5 * ((A)[2] + (A)[6]);
		(a)[3] = ((A)[4]);
		(a)[4] = 0.5 * ((A)[5] + (A)[7]);
		(a)[5] = ((A)[8]);
	}
#define MatrixSymPart3X3d MatrixSymPart3X3<double>
#define MatrixSymPart3X3f MatrixSymPart3X3<float>


// inverse a 3x3 matrix, must input non-singular matrix
#define MATRIX_INV3X3(a, inv, det)\
	{\
	(inv)[0]=((a)[4]*(a)[8]-(a)[7]*(a)[5]);\
	(inv)[3]=-((a)[3]*(a)[8]-(a)[5]*(a)[6]);\
	(inv)[6]=((a)[3]*(a)[7]-(a)[6]*(a)[4]);\
	(inv)[1]=-((a)[1]*(a)[8]-(a)[2]*(a)[7]);\
	(inv)[4]=((a)[0]*(a)[8]-(a)[2]*(a)[6]);\
	(inv)[7]=-((a)[0]*(a)[7]-(a)[6]*(a)[1]);\
	(inv)[2]=((a)[1]*(a)[5]-(a)[2]*(a)[4]);\
	(inv)[5]=-((a)[0]*(a)[5]-(a)[3]*(a)[2]);\
	(inv)[8]=((a)[0]*(a)[4]-(a)[3]*(a)[1]);\
	det = ((a)[0] *(inv)[0] + (a)[1] * (inv)[3] + (a)[2] * (inv)[6]);\
	(inv)[0]=(inv)[0]/det;\
	(inv)[3]=(inv)[3]/det;\
	(inv)[6]=(inv)[6]/det;\
	(inv)[1]=(inv)[1]/det;\
	(inv)[4]=(inv)[4]/det;\
	(inv)[7]=(inv)[7]/det;\
	(inv)[2]=(inv)[2]/det;\
	(inv)[5]=(inv)[5]/det;\
	(inv)[8]=(inv)[8]/det;}

	template <class TYPE>
	inline TYPE MatrixInv3X3(TYPE *a, TYPE *inv)
	{
		(inv)[0] = ((a)[4] * (a)[8] - (a)[7] * (a)[5]);
		(inv)[3] = -((a)[3] * (a)[8] - (a)[5] * (a)[6]);
		(inv)[6] = ((a)[3] * (a)[7] - (a)[6] * (a)[4]);
		(inv)[1] = -((a)[1] * (a)[8] - (a)[2] * (a)[7]);
		(inv)[4] = ((a)[0] * (a)[8] - (a)[2] * (a)[6]);
		(inv)[7] = -((a)[0] * (a)[7] - (a)[6] * (a)[1]);
		(inv)[2] = ((a)[1] * (a)[5] - (a)[2] * (a)[4]);
		(inv)[5] = -((a)[0] * (a)[5] - (a)[3] * (a)[2]);
		(inv)[8] = ((a)[0] * (a)[4] - (a)[3] * (a)[1]);
		TYPE det = ((a)[0] * (inv)[0] + (a)[1] * (inv)[3] + (a)[2] * (inv)[6]);
		(inv)[0] = (inv)[0] / det;
		(inv)[3] = (inv)[3] / det;
		(inv)[6] = (inv)[6] / det;
		(inv)[1] = (inv)[1] / det;
		(inv)[4] = (inv)[4] / det;
		(inv)[7] = (inv)[7] / det;
		(inv)[2] = (inv)[2] / det;
		(inv)[5] = (inv)[5] / det;
		(inv)[8] = (inv)[8] / det;
		return det;
	}
#define MatrixInv3X3d MatrixInv3X3<double>
#define MatrixInv3X3f MatrixInv3X3<float>


#define MATRIX_PRINT3X3(a)\
			{std::cout<<"[" << (a)[0]<<", "<<(a)[1]<<", "<<(a)[2]<<";\n";\
	std::cout<<(a)[3]<<", "<<(a)[4]<<", "<<(a)[5]<<";\n";\
	std::cout<<(a)[6]<<", "<<(a)[7]<<", "<<(a)[8]<<"]";}

	template <class TYPE>
	inline void MatrixPrint3X3(TYPE *a)
	{
		std::cout << "[" << (a)[0] << ", " << (a)[1] << ", " << (a)[2] << "; ";
		std::cout << (a)[3] << ", " << (a)[4] << ", " << (a)[5] << "; ";
		std::cout << (a)[6] << ", " << (a)[7] << ", " << (a)[8] << "]\n";
	}
#define MatrixPrint3X3d MatrixPrint3X3<double>
#define MatrixPrint3X3f MatrixPrint3X3<float>



#define MATRIX_AT3X3(a, x, y)	((a)[3*(x)+y])

	template <class TYPE>
	inline TYPE MatrixAt3X3(TYPE *a, int x, int y)
	{
		return ((a)[3 * (x)+y]);
	}
#define MatrixAt3X3d MatrixAt3X3<double>
#define MatrixAt3X3f MatrixAt3X3<float>


template <class TYPE>
class Mat3 {
public:
	inline Mat3()
	{
		MATRIX_SET_IDENTITY3X3(v);
	}

	inline Mat3(TYPE x0, TYPE x1, TYPE x2,
		TYPE x3, TYPE x4, TYPE x5,
		TYPE x6, TYPE x7, TYPE x8)
	{
		v[0] = x0; v[1] = x1; v[2] = x2;
		v[3] = x3; v[4] = x4; v[5] = x5;
		v[6] = x6; v[7] = x7; v[8] = x8;
	}

	inline Mat3(TYPE *mat)
	{
		MATRIX_SET3X3(v, mat);
	}

	inline Mat3(Vec3<TYPE> row0, Vec3<TYPE> row1, Vec3<TYPE> row2)
	{
		v[0] = row0[0]; v[1] = row0[1]; v[2] = row0[2];
		v[3] = row1[0]; v[4] = row1[1]; v[5] = row1[2];
		v[6] = row2[0]; v[7] = row2[1]; v[8] = row2[2];
	}
	inline Mat3(TYPE diag) // create a diagonal matrix with all entries "diag" (can create zero matrix by passing 0.0)
	{
		MATRIX_SET_DIAG3X3(v, diag);
	}

	inline void set(TYPE x0, TYPE x1, TYPE x2,
		TYPE x3, TYPE x4, TYPE x5,
		TYPE x6, TYPE x7, TYPE x8)
	{
		v[0] = x0; v[1] = x1; v[2] = x2;
		v[3] = x3; v[4] = x4; v[5] = x5;
		v[6] = x6; v[7] = x7; v[8] = x8;
	}

	inline TYPE & operator[] (int index)
	{
		return v[index];
	}
	inline TYPE & operator() (int index1, int index2)
	{
		return v[3*index1+index2];
	}
public:
	TYPE v[9];
};
typedef Mat3<int> Mat3i;
typedef Mat3<float> Mat3f;
typedef Mat3<double> Mat3d;




/* Eigen-decomposition for symmetric 3x3 real matrices.
Public domain, copied from the public domain Java library JAMA. */
/* Symmetric matrix A => eigenvectors in columns of V, corresponding
eigenvalues in d. */
void eigenDecompose3X3(double *A, double V[3][3], double d[3]);
void testEigenDecomposition3X3();

/*
Standard SVD (modifiedSVD == 0):

Given a 3x3 matrix M, decomposes it using SVD so
that M = U Sigma V^T where U is a 3x3 rotation,
V is 3x3 orthonormal (V^T V = V V^T = I), and
Sigma is a diagonal matrix with non-negative entries,
in descending order, Sigma[0] >= Sigma[1] >= Sigma[2] >= 0.
Note that the matrix V may have a determinant equaling
to -1 (i.e., reflection).

singularValue_eps is a threshold to determine
when a singular value is deemed zero, and special handling is then invoked
to improve numerical robustness.

Modified SVD (modifiedSVD == 1):

The SVD is modified so that it gives the
following properties (useful in solid mechanics applications) :

1) Not just the determinant of U, but also the determinant of V is 1
(i.e., both U and V are rotations, not reflections).

2) The smallest diagonal value Sigma[2] may be negative. We have:
Sigma[0] >= Sigma[1] >= abs(Sigma[2]) >= 0 .
*/
void SVD3X3(double* M, double* U, double* Sigma, double* V, double singularValue_eps, int modifiedSVD);

double matrixOneNorm3X3(double *A);
double matrixInfNorm3X3(double *A);

// Computes the Polar Decomposition of a general 3x3 matrix M.
// M = Q * S
// M is 3x3 input matrix
// Q is 3x3 orthogonal output matrix, Q Q^T = Q^T Q = I 
// S is 3x3 symmetric positive-definite output matrix
// Note: det(Q)=sgn(det(M)); this sign can be 1 or -1, depending on M
// if forceRotation is 1, the matrix Q will be a rotation, S will be symmetric, but not necessarily positive-definite
// M is not modified
// All matrices are row-major
double polarDecompose3X3(double * M, double * Q, double * S, int forceRotation = 0, double tolerance = 1E-6);

}

#endif