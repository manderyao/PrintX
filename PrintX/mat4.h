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

4D matrix operations
*/


#ifndef __PRINTX_MAT4_H__
#define __PRINTX_MAT4_H__
#include <math.h>
#include <algorithm>
namespace PRINTX_LIB
{
#define MATRIX_SET_IDENTITY4X4(a)\
	{	(a)[0] = 1;(a)[1] = 0;(a)[2] = 0;(a)[3] = 0;\
	(a)[4] = 0;(a)[5] = 1;(a)[6] = 0;(a)[7] = 0;\
	(a)[8] = 0;(a)[9] = 0;(a)[10] = 1;(a)[11] = 0;\
	(a)[12] = 0;(a)[13] = 0;(a)[14] = 0;(a)[15] = 1;}

	template <class TYPE>
	inline void MatrixSetIdentity4X4(TYPE *a)
	{
		(a)[0] = 1; (a)[1] = 0; (a)[2] = 0; (a)[3] = 0;
		(a)[4] = 0; (a)[5] = 1; (a)[6] = 0; (a)[7] = 0;
		(a)[8] = 0; (a)[9] = 0; (a)[10] = 1; (a)[11] = 0;
		(a)[12] = 0; (a)[13] = 0; (a)[14] = 0; (a)[15] = 1;
	}
#define MatrixSetIdentity4X4d MatrixSetIdentity4X4<double>
#define MatrixSetIdentity4X4f MatrixSetIdentity4X4<float>


	// a = b
#define MATRIX_SET4X4(a,b)\
	{(a)[0] = (b)[0]; (a)[1] = (b)[1]; (a)[2] = (b)[2]; (a)[3] = (b)[3];\
	(a)[4] = (b)[4]; (a)[5] = (b)[5]; (a)[6] = (b)[6]; (a)[7] = (b)[7];\
	(a)[8] = (b)[8]; (a)[9] = (b)[9]; (a)[10] = (b)[10]; (a)[11] = (b)[11];\
	(a)[12] = (b)[12]; (a)[13] = (b)[13]; (a)[14] = (b)[14]; (a)[15] = (b)[15];}

	template <class TYPE>
	inline void MatrixSet4X4(TYPE *a, TYPE *b)
	{
		(a)[0] = (b)[0]; (a)[1] = (b)[1]; (a)[2] = (b)[2]; (a)[3] = (b)[3];
		(a)[4] = (b)[4]; (a)[5] = (b)[5]; (a)[6] = (b)[6]; (a)[7] = (b)[7];
		(a)[8] = (b)[8]; (a)[9] = (b)[9]; (a)[10] = (b)[10]; (a)[11] = (b)[11];
		(a)[12] = (b)[12]; (a)[13] = (b)[13]; (a)[14] = (b)[14]; (a)[15] = (b)[15];
	}
#define MatrixSet4X4d MatrixSet4X4<double>
#define MatrixSet4X4f MatrixSet4X4<float>


	// c = a * b
#define MATRIX_MULTIPLY4X4(a,b,c)\
	{(c)[0]=(a)[0]*(b)[0]+(a)[1]*(b)[4]+(a)[2]*(b)[8]+(a)[3]*(b)[12];\
	(c)[1]=(a)[0]*(b)[1]+(a)[1]*(b)[5]+(a)[2]*(b)[9]+(a)[3]*(b)[13];\
	(c)[2]=(a)[0]*(b)[2]+(a)[1]*(b)[6]+(a)[2]*(b)[10]+(a)[3]*(b)[14];\
	(c)[3]=(a)[0]*(b)[3]+(a)[1]*(b)[7]+(a)[2]*(b)[11]+(a)[3]*(b)[15];\
	(c)[4]=(a)[4]*(b)[0]+(a)[5]*(b)[4]+(a)[6]*(b)[8]+(a)[7]*(b)[12];\
	(c)[5]=(a)[4]*(b)[1]+(a)[5]*(b)[5]+(a)[6]*(b)[9]+(a)[7]*(b)[13];\
	(c)[6]=(a)[4]*(b)[2]+(a)[5]*(b)[6]+(a)[6]*(b)[10]+(a)[7]*(b)[14];\
	(c)[7]=(a)[4]*(b)[3]+(a)[5]*(b)[7]+(a)[6]*(b)[11]+(a)[7]*(b)[15];\
	(c)[8]=(a)[8]*(b)[0]+(a)[9]*(b)[4]+(a)[10]*(b)[8]+(a)[11]*(b)[12];\
	(c)[9]=(a)[8]*(b)[1]+(a)[9]*(b)[5]+(a)[10]*(b)[9]+(a)[11]*(b)[13];\
	(c)[10]=(a)[8]*(b)[2]+(a)[9]*(b)[6]+(a)[10]*(b)[10]+(a)[11]*(b)[14];\
	(c)[11]=(a)[8]*(b)[3]+(a)[9]*(b)[7]+(a)[10]*(b)[11]+(a)[11]*(b)[15];\
	(c)[12]=(a)[12]*(b)[0]+(a)[13]*(b)[4]+(a)[14]*(b)[8]+(a)[15]*(b)[12];\
	(c)[13]=(a)[12]*(b)[1]+(a)[13]*(b)[5]+(a)[14]*(b)[9]+(a)[15]*(b)[13];\
	(c)[14]=(a)[12]*(b)[2]+(a)[13]*(b)[6]+(a)[14]*(b)[10]+(a)[15]*(b)[14];\
	(c)[15]=(a)[12]*(b)[3]+(a)[13]*(b)[7]+(a)[14]*(b)[11]+(a)[15]*(b)[15];}

	template <class TYPE>
	inline void MatrixMultiply4X4(TYPE *a, TYPE *b, TYPE *c)
	{
		(c)[0] = (a)[0] * (b)[0] + (a)[1] * (b)[4] + (a)[2] * (b)[8] + (a)[3] * (b)[12];
		(c)[1] = (a)[0] * (b)[1] + (a)[1] * (b)[5] + (a)[2] * (b)[9] + (a)[3] * (b)[13];
		(c)[2] = (a)[0] * (b)[2] + (a)[1] * (b)[6] + (a)[2] * (b)[10] + (a)[3] * (b)[14];
		(c)[3] = (a)[0] * (b)[3] + (a)[1] * (b)[7] + (a)[2] * (b)[11] + (a)[3] * (b)[15];
		(c)[4] = (a)[4] * (b)[0] + (a)[5] * (b)[4] + (a)[6] * (b)[8] + (a)[7] * (b)[12];
		(c)[5] = (a)[4] * (b)[1] + (a)[5] * (b)[5] + (a)[6] * (b)[9] + (a)[7] * (b)[13];
		(c)[6] = (a)[4] * (b)[2] + (a)[5] * (b)[6] + (a)[6] * (b)[10] + (a)[7] * (b)[14];
		(c)[7] = (a)[4] * (b)[3] + (a)[5] * (b)[7] + (a)[6] * (b)[11] + (a)[7] * (b)[15];
		(c)[8] = (a)[8] * (b)[0] + (a)[9] * (b)[4] + (a)[10] * (b)[8] + (a)[11] * (b)[12];
		(c)[9] = (a)[8] * (b)[1] + (a)[9] * (b)[5] + (a)[10] * (b)[9] + (a)[11] * (b)[13];
		(c)[10] = (a)[8] * (b)[2] + (a)[9] * (b)[6] + (a)[10] * (b)[10] + (a)[11] * (b)[14];
		(c)[11] = (a)[8] * (b)[3] + (a)[9] * (b)[7] + (a)[10] * (b)[11] + (a)[11] * (b)[15];
		(c)[12] = (a)[12] * (b)[0] + (a)[13] * (b)[4] + (a)[14] * (b)[8] + (a)[15] * (b)[12];
		(c)[13] = (a)[12] * (b)[1] + (a)[13] * (b)[5] + (a)[14] * (b)[9] + (a)[15] * (b)[13];
		(c)[14] = (a)[12] * (b)[2] + (a)[13] * (b)[6] + (a)[14] * (b)[10] + (a)[15] * (b)[14];
		(c)[15] = (a)[12] * (b)[3] + (a)[13] * (b)[7] + (a)[14] * (b)[11] + (a)[15] * (b)[15];
	}
#define MatrixMultiply4X4d MatrixMultiply4X4<double>
#define MatrixMultiply4X4f MatrixMultiply4X4<float>



#define MATRIX_DETERMINANT4X4(a)\
	((a)[0]*((a)[5]*(a)[10]*(a)[15]-(a)[5]*(a)[11]*(a)[14]-(a)[9]*(a)[6]*(a)[15]+(a)[9]*(a)[7]*(a)[14]+(a)[13]*(a)[6]*(a)[11]-(a)[13]*(a)[7]*(a)[10])+\
	(a)[1]*(-(a)[4] * (a)[10] * (a)[15] +(a)[4] * (a)[11] * (a)[14] +(a)[8] * (a)[6] * (a)[15] -(a)[8] * (a)[7] * (a)[14] -(a)[12] * (a)[6] * (a)[11] +(a)[12] * (a)[7] * (a)[10])+\
	(a)[2]*((a)[4] * (a)[9] * (a)[15] -(a)[4] * (a)[11] * (a)[13] -(a)[8] * (a)[5] * (a)[15] +(a)[8] * (a)[7] * (a)[13] +(a)[12] * (a)[5] * (a)[11] -(a)[12] * (a)[7] * (a)[9])+\
	(a)[3]*(-(a)[4] * (a)[9] * (a)[14] +(a)[4] * (a)[10] * (a)[13] +(a)[8] * (a)[5] * (a)[14] -(a)[8] * (a)[6] * (a)[13] -(a)[12] * (a)[5] * (a)[10] +(a)[12] * (a)[6] * (a)[9]))

	template <class TYPE>
	inline TYPE MatrixDeterminant4X4(TYPE *a)
	{
		return ((a)[0] * ((a)[5] * (a)[10] * (a)[15] - (a)[5] * (a)[11] * (a)[14] - (a)[9] * (a)[6] * (a)[15] + (a)[9] * (a)[7] * (a)[14] + (a)[13] * (a)[6] * (a)[11] - (a)[13] * (a)[7] * (a)[10]) +
			(a)[1] * (-(a)[4] * (a)[10] * (a)[15] + (a)[4] * (a)[11] * (a)[14] + (a)[8] * (a)[6] * (a)[15] - (a)[8] * (a)[7] * (a)[14] - (a)[12] * (a)[6] * (a)[11] + (a)[12] * (a)[7] * (a)[10]) +
			(a)[2] * ((a)[4] * (a)[9] * (a)[15] - (a)[4] * (a)[11] * (a)[13] - (a)[8] * (a)[5] * (a)[15] + (a)[8] * (a)[7] * (a)[13] + (a)[12] * (a)[5] * (a)[11] - (a)[12] * (a)[7] * (a)[9]) +
			(a)[3] * (-(a)[4] * (a)[9] * (a)[14] + (a)[4] * (a)[10] * (a)[13] + (a)[8] * (a)[5] * (a)[14] - (a)[8] * (a)[6] * (a)[13] - (a)[12] * (a)[5] * (a)[10] + (a)[12] * (a)[6] * (a)[9]));

	}
#define MatrixDeterminant4X4d MatrixDeterminant4X4<double>
#define MatrixDeterminant4X4f MatrixDeterminant4X4<float>



	// B = A^T
#define MATRIX_TRANSPOSE4X4(A,B)\
	{(B)[0] = (A)[0]; (B)[1] = (A)[4]; (B)[2] = (A)[8];(B)[3] = (A)[12];\
	(B)[4] = (A)[1]; (B)[5] = (A)[5]; (B)[6] = (A)[9];(B)[7] = (A)[13];\
	(B)[8] = (A)[2]; (B)[9] = (A)[6]; (B)[10] = (A)[10];(B)[11] = (A)[14];\
	(B)[12] = (A)[3]; (B)[13] = (A)[7]; (B)[14] = (A)[11]; (B)[15] = (A)[15];}

	template <class TYPE>
	inline void MatrixTranspose4X4(TYPE *A, TYPE *B)
	{
		(B)[0] = (A)[0]; (B)[1] = (A)[4]; (B)[2] = (A)[8]; (B)[3] = (A)[12];
		(B)[4] = (A)[1]; (B)[5] = (A)[5]; (B)[6] = (A)[9]; (B)[7] = (A)[13];
		(B)[8] = (A)[2]; (B)[9] = (A)[6]; (B)[10] = (A)[10]; (B)[11] = (A)[14];
		(B)[12] = (A)[3]; (B)[13] = (A)[7]; (B)[14] = (A)[11]; (B)[15] = (A)[15];
	}
#define MatrixTranspose4X4d MatrixTranspose4X4<double>
#define MatrixTranspose4X4f MatrixTranspose4X4<float>



#define MATRIX_SELF_TRANSPOSE4X4(A)\
	{std::swap((A)[1], (A)[4]); std::swap((A)[2], (A)[8]); std::swap((A)[3], (A)[12]);\
	std::swap((A)[6], (A)[9]); std::swap((A)[7], (A)[13]); std::swap((A)[11], (A)[14]);}

	template <class TYPE>
	inline void MatrixSelfTranspose4X4(TYPE *A)
	{
		std::swap((A)[1], (A)[4]); std::swap((A)[2], (A)[8]); std::swap((A)[3], (A)[12]);
		std::swap((A)[6], (A)[9]); std::swap((A)[7], (A)[13]); std::swap((A)[11], (A)[14]);
	}
#define MatrixSelfTranspose4X4d MatrixSelfTranspose4X4<double>
#define MatrixSelfTranspose4X4f MatrixSelfTranspose4X4<float>



// inverse 4x4 matrix, must be non-singular
#define MATRIX_INV_4X4(a, inv, det)\
	{(inv)[0]=((a)[5]*(a)[10]*(a)[15]-(a)[5]*(a)[11]*(a)[14]-(a)[9]*(a)[6]*(a)[15]+(a)[9]*(a)[7]*(a)[14]+(a)[13]*(a)[6]*(a)[11]-(a)[13]*(a)[7]*(a)[10]);\
	(inv)[4]=(-(a)[4] * (a)[10] * (a)[15] +(a)[4] * (a)[11] * (a)[14] +(a)[8] * (a)[6] * (a)[15] -(a)[8] * (a)[7] * (a)[14] -(a)[12] * (a)[6] * (a)[11] +(a)[12] * (a)[7] * (a)[10]);\
	(inv)[8]=((a)[4] * (a)[9] * (a)[15] -(a)[4] * (a)[11] * (a)[13] -(a)[8] * (a)[5] * (a)[15] +(a)[8] * (a)[7] * (a)[13] +(a)[12] * (a)[5] * (a)[11] -(a)[12] * (a)[7] * (a)[9]);\
	(inv)[12]=(-(a)[4] * (a)[9] * (a)[14] +(a)[4] * (a)[10] * (a)[13] +(a)[8] * (a)[5] * (a)[14] -(a)[8] * (a)[6] * (a)[13] -(a)[12] * (a)[5] * (a)[10] +(a)[12] * (a)[6] * (a)[9]);\
	(inv)[1] = -(a)[1] * (a)[10] * (a)[15] +(a)[1] * (a)[11] * (a)[14] +(a)[9] * (a)[2] * (a)[15] -(a)[9] * (a)[3] * (a)[14] -(a)[13] * (a)[2] * (a)[11] +(a)[13] * (a)[3] * (a)[10];\
	(inv)[5] = (a)[0] * (a)[10] * (a)[15] -(a)[0] * (a)[11] * (a)[14] -(a)[8] * (a)[2] * (a)[15] +(a)[8] * (a)[3] * (a)[14] +(a)[12] * (a)[2] * (a)[11] -(a)[12] * (a)[3] * (a)[10];\
	(inv)[9] = -(a)[0] * (a)[9] * (a)[15] +(a)[0] * (a)[11] * (a)[13] +(a)[8] * (a)[1] * (a)[15] -(a)[8] * (a)[3] * (a)[13] -(a)[12] * (a)[1] * (a)[11] +(a)[12] * (a)[3] * (a)[9];\
	(inv)[13] = (a)[0] * (a)[9] * (a)[14] -(a)[0] * (a)[10] * (a)[13] -(a)[8] * (a)[1] * (a)[14] +(a)[8] * (a)[2] * (a)[13] +(a)[12] * (a)[1] * (a)[10] -(a)[12] * (a)[2] * (a)[9];\
	(inv)[2] = (a)[1] * (a)[6] * (a)[15] -(a)[1] * (a)[7] * (a)[14] -(a)[5] * (a)[2] * (a)[15] +(a)[5] * (a)[3] * (a)[14] +(a)[13] * (a)[2] * (a)[7] -(a)[13] * (a)[3] * (a)[6];\
	(inv)[6] = -(a)[0] * (a)[6] * (a)[15] +(a)[0] * (a)[7] * (a)[14] +(a)[4] * (a)[2] * (a)[15] -(a)[4] * (a)[3] * (a)[14] -(a)[12] * (a)[2] * (a)[7] +(a)[12] * (a)[3] * (a)[6];\
	(inv)[10] = (a)[0] * (a)[5] * (a)[15] -(a)[0] * (a)[7] * (a)[13] -(a)[4] * (a)[1] * (a)[15] +(a)[4] * (a)[3] * (a)[13] +(a)[12] * (a)[1] * (a)[7] -(a)[12] * (a)[3] * (a)[5];\
	(inv)[14] = -(a)[0] * (a)[5] * (a)[14] +(a)[0] * (a)[6] * (a)[13] +(a)[4] * (a)[1] * (a)[14] -(a)[4] * (a)[2] * (a)[13] -(a)[12] * (a)[1] * (a)[6] +(a)[12] * (a)[2] * (a)[5];\
	(inv)[3] = -(a)[1] * (a)[6] * (a)[11] +(a)[1] * (a)[7] * (a)[10] +(a)[5] * (a)[2] * (a)[11] -(a)[5] * (a)[3] * (a)[10] -(a)[9] * (a)[2] * (a)[7] +(a)[9] * (a)[3] * (a)[6];\
	(inv)[7] = (a)[0] * (a)[6] * (a)[11] -(a)[0] * (a)[7] * (a)[10] -(a)[4] * (a)[2] * (a)[11] +(a)[4] * (a)[3] * (a)[10] +(a)[8] * (a)[2] * (a)[7] -(a)[8] * (a)[3] * (a)[6];\
	(inv)[11] = -(a)[0] * (a)[5] * (a)[11] +(a)[0] * (a)[7] * (a)[9] +(a)[4] * (a)[1] * (a)[11] -(a)[4] * (a)[3] * (a)[9] -(a)[8] * (a)[1] * (a)[7] +(a)[8] * (a)[3] * (a)[5];\
	(inv)[15] = (a)[0] * (a)[5] * (a)[10] -(a)[0] * (a)[6] * (a)[9] -(a)[4] * (a)[1] * (a)[10] +(a)[4] * (a)[2] * (a)[9] +(a)[8] * (a)[1] * (a)[6] -(a)[8] * (a)[2] * (a)[5];\
	det = (a)[0] * (inv)[0] + (a)[1] * (inv)[4] + (a)[2] * (inv)[8] + (a)[3] * (inv)[12];\
	(inv)[0]=(inv)[0]/det;\
	(inv)[1]=(inv)[1]/det;\
	(inv)[2]=(inv)[2]/det;\
	(inv)[3]=(inv)[3]/det;\
	(inv)[4]=(inv)[4]/det;\
	(inv)[5]=(inv)[5]/det;\
	(inv)[6]=(inv)[6]/det;\
	(inv)[7]=(inv)[7]/det;\
	(inv)[8]=(inv)[8]/det;\
	(inv)[9] = (inv)[9] / det; \
	(inv)[10] = (inv)[10] / det; \
	(inv)[11] = (inv)[11] / det; \
	(inv)[12] = (inv)[12] / det; \
	(inv)[13] = (inv)[13] / det; \
	(inv)[14] = (inv)[14] / det; \
	(inv)[15] = (inv)[15] / det;}

	template <class TYPE>
	inline TYPE MatrixInv4X4(TYPE *a, TYPE *inv)
	{
		(inv)[0] = ((a)[5] * (a)[10] * (a)[15] - (a)[5] * (a)[11] * (a)[14] - (a)[9] * (a)[6] * (a)[15] + (a)[9] * (a)[7] * (a)[14] + (a)[13] * (a)[6] * (a)[11] - (a)[13] * (a)[7] * (a)[10]);
		(inv)[4] = (-(a)[4] * (a)[10] * (a)[15] + (a)[4] * (a)[11] * (a)[14] + (a)[8] * (a)[6] * (a)[15] - (a)[8] * (a)[7] * (a)[14] - (a)[12] * (a)[6] * (a)[11] + (a)[12] * (a)[7] * (a)[10]);
		(inv)[8] = ((a)[4] * (a)[9] * (a)[15] - (a)[4] * (a)[11] * (a)[13] - (a)[8] * (a)[5] * (a)[15] + (a)[8] * (a)[7] * (a)[13] + (a)[12] * (a)[5] * (a)[11] - (a)[12] * (a)[7] * (a)[9]);
		(inv)[12] = (-(a)[4] * (a)[9] * (a)[14] + (a)[4] * (a)[10] * (a)[13] + (a)[8] * (a)[5] * (a)[14] - (a)[8] * (a)[6] * (a)[13] - (a)[12] * (a)[5] * (a)[10] + (a)[12] * (a)[6] * (a)[9]);
		(inv)[1] = -(a)[1] * (a)[10] * (a)[15] + (a)[1] * (a)[11] * (a)[14] + (a)[9] * (a)[2] * (a)[15] - (a)[9] * (a)[3] * (a)[14] - (a)[13] * (a)[2] * (a)[11] + (a)[13] * (a)[3] * (a)[10];
		(inv)[5] = (a)[0] * (a)[10] * (a)[15] - (a)[0] * (a)[11] * (a)[14] - (a)[8] * (a)[2] * (a)[15] + (a)[8] * (a)[3] * (a)[14] + (a)[12] * (a)[2] * (a)[11] - (a)[12] * (a)[3] * (a)[10];
		(inv)[9] = -(a)[0] * (a)[9] * (a)[15] + (a)[0] * (a)[11] * (a)[13] + (a)[8] * (a)[1] * (a)[15] - (a)[8] * (a)[3] * (a)[13] - (a)[12] * (a)[1] * (a)[11] + (a)[12] * (a)[3] * (a)[9];
		(inv)[13] = (a)[0] * (a)[9] * (a)[14] - (a)[0] * (a)[10] * (a)[13] - (a)[8] * (a)[1] * (a)[14] + (a)[8] * (a)[2] * (a)[13] + (a)[12] * (a)[1] * (a)[10] - (a)[12] * (a)[2] * (a)[9];
		(inv)[2] = (a)[1] * (a)[6] * (a)[15] - (a)[1] * (a)[7] * (a)[14] - (a)[5] * (a)[2] * (a)[15] + (a)[5] * (a)[3] * (a)[14] + (a)[13] * (a)[2] * (a)[7] - (a)[13] * (a)[3] * (a)[6];
		(inv)[6] = -(a)[0] * (a)[6] * (a)[15] + (a)[0] * (a)[7] * (a)[14] + (a)[4] * (a)[2] * (a)[15] - (a)[4] * (a)[3] * (a)[14] - (a)[12] * (a)[2] * (a)[7] + (a)[12] * (a)[3] * (a)[6];
		(inv)[10] = (a)[0] * (a)[5] * (a)[15] - (a)[0] * (a)[7] * (a)[13] - (a)[4] * (a)[1] * (a)[15] + (a)[4] * (a)[3] * (a)[13] + (a)[12] * (a)[1] * (a)[7] - (a)[12] * (a)[3] * (a)[5];
		(inv)[14] = -(a)[0] * (a)[5] * (a)[14] + (a)[0] * (a)[6] * (a)[13] + (a)[4] * (a)[1] * (a)[14] - (a)[4] * (a)[2] * (a)[13] - (a)[12] * (a)[1] * (a)[6] + (a)[12] * (a)[2] * (a)[5];
		(inv)[3] = -(a)[1] * (a)[6] * (a)[11] + (a)[1] * (a)[7] * (a)[10] + (a)[5] * (a)[2] * (a)[11] - (a)[5] * (a)[3] * (a)[10] - (a)[9] * (a)[2] * (a)[7] + (a)[9] * (a)[3] * (a)[6];
		(inv)[7] = (a)[0] * (a)[6] * (a)[11] - (a)[0] * (a)[7] * (a)[10] - (a)[4] * (a)[2] * (a)[11] + (a)[4] * (a)[3] * (a)[10] + (a)[8] * (a)[2] * (a)[7] - (a)[8] * (a)[3] * (a)[6];
		(inv)[11] = -(a)[0] * (a)[5] * (a)[11] + (a)[0] * (a)[7] * (a)[9] + (a)[4] * (a)[1] * (a)[11] - (a)[4] * (a)[3] * (a)[9] - (a)[8] * (a)[1] * (a)[7] + (a)[8] * (a)[3] * (a)[5];
		(inv)[15] = (a)[0] * (a)[5] * (a)[10] - (a)[0] * (a)[6] * (a)[9] - (a)[4] * (a)[1] * (a)[10] + (a)[4] * (a)[2] * (a)[9] + (a)[8] * (a)[1] * (a)[6] - (a)[8] * (a)[2] * (a)[5];
		TYPE det = (a)[0] * (inv)[0] + (a)[1] * (inv)[4] + (a)[2] * (inv)[8] + (a)[3] * (inv)[12];
		(inv)[0] = (inv)[0] / det;
		(inv)[1] = (inv)[1] / det;
		(inv)[2] = (inv)[2] / det;
		(inv)[3] = (inv)[3] / det;
		(inv)[4] = (inv)[4] / det;
		(inv)[5] = (inv)[5] / det;
		(inv)[6] = (inv)[6] / det;
		(inv)[7] = (inv)[7] / det;
		(inv)[8] = (inv)[8] / det;
		(inv)[9] = (inv)[9] / det;
		(inv)[10] = (inv)[10] / det;
		(inv)[11] = (inv)[11] / det;
		(inv)[12] = (inv)[12] / det;
		(inv)[13] = (inv)[13] / det;
		(inv)[14] = (inv)[14] / det;
		(inv)[15] = (inv)[15] / det;
		return det;
	}
#define MatrixInv4X4d MatrixInv4X4<double>
#define MatrixInv4X4f MatrixInv4X4<float>


	// w = A * v
#define MATRIX_VECTOR_MULTIPLY4X3(A, v, w) \
	{(w)[0] = (A)[0]*(v)[0]+(A)[1]*(v)[1]+(A)[2]*(v)[2]+(A)[3];\
	(w)[1] = (A)[4]*(v)[0]+(A)[5]*(v)[1]+(A)[6]*(v)[2]+(A)[7];\
	(w)[2] = (A)[8]*(v)[0]+(A)[9]*(v)[1]+(A)[10]*(v)[2]+(A)[11];}

	template <class TYPE>
	inline void MatrixVectorMultiply4X3(TYPE *A, TYPE *v, TYPE *w)
	{
		(w)[0] = (A)[0] * (v)[0] + (A)[1] * (v)[1] + (A)[2] * (v)[2] + (A)[3];
		(w)[1] = (A)[4] * (v)[0] + (A)[5] * (v)[1] + (A)[6] * (v)[2] + (A)[7];
		(w)[2] = (A)[8] * (v)[0] + (A)[9] * (v)[1] + (A)[10] * (v)[2] + (A)[11];
	}
#define MatrixVectorMultiply4X3d MatrixVectorMultiply4X3<double>
#define MatrixVectorMultiply4X3f MatrixVectorMultiply4X3<float>


	// c = a * b
#define MATRIX_MULTIPLY4X3(a,b,c)\
	{(c)[0]=(a)[0]*(b)[0]+(a)[1]*(b)[3]+(a)[2]*(b)[6];\
	(c)[1]=(a)[0]*(b)[1]+(a)[1]*(b)[4]+(a)[2]*(b)[7];\
	(c)[2]=(a)[0]*(b)[2]+(a)[1]*(b)[5]+(a)[2]*(b)[8];\
	(c)[3]=(a)[3];\
	(c)[4]=(a)[4]*(b)[0]+(a)[5]*(b)[3]+(a)[6]*(b)[6];\
	(c)[5]=(a)[4]*(b)[1]+(a)[5]*(b)[4]+(a)[6]*(b)[7];\
	(c)[6]=(a)[4]*(b)[2]+(a)[5]*(b)[5]+(a)[6]*(b)[8];\
	(c)[7]=(a)[7];\
	(c)[8]=(a)[8]*(b)[0]+(a)[9]*(b)[3]+(a)[10]*(b)[6];\
	(c)[9]=(a)[8]*(b)[1]+(a)[9]*(b)[4]+(a)[10]*(b)[7];\
	(c)[10]=(a)[8]*(b)[2]+(a)[9]*(b)[5]+(a)[10]*(b)[8];\
	(c)[11]=(a)[11];\
	(c)[12]=(a)[12]*(b)[0]+(a)[13]*(b)[3]+(a)[14]*(b)[6];\
	(c)[13]=(a)[12]*(b)[1]+(a)[13]*(b)[4]+(a)[14]*(b)[7];\
	(c)[14]=(a)[12]*(b)[2]+(a)[13]*(b)[5]+(a)[14]*(b)[8];\
	(c)[15]=(a)[15];}

	template <class TYPE>
	inline void MatrixMultiply4X3(TYPE *a, TYPE *b, TYPE *c)
	{
		(c)[0] = (a)[0] * (b)[0] + (a)[1] * (b)[3] + (a)[2] * (b)[6];
		(c)[1] = (a)[0] * (b)[1] + (a)[1] * (b)[4] + (a)[2] * (b)[7];
		(c)[2] = (a)[0] * (b)[2] + (a)[1] * (b)[5] + (a)[2] * (b)[8];
		(c)[3] = (a)[3];
		(c)[4] = (a)[4] * (b)[0] + (a)[5] * (b)[3] + (a)[6] * (b)[6];
		(c)[5] = (a)[4] * (b)[1] + (a)[5] * (b)[4] + (a)[6] * (b)[7];
		(c)[6] = (a)[4] * (b)[2] + (a)[5] * (b)[5] + (a)[6] * (b)[8];
		(c)[7] = (a)[7];
		(c)[8] = (a)[8] * (b)[0] + (a)[9] * (b)[3] + (a)[10] * (b)[6];
		(c)[9] = (a)[8] * (b)[1] + (a)[9] * (b)[4] + (a)[10] * (b)[7];
		(c)[10] = (a)[8] * (b)[2] + (a)[9] * (b)[5] + (a)[10] * (b)[8];
		(c)[11] = (a)[11];
		(c)[12] = (a)[12] * (b)[0] + (a)[13] * (b)[3] + (a)[14] * (b)[6];
		(c)[13] = (a)[12] * (b)[1] + (a)[13] * (b)[4] + (a)[14] * (b)[7];
		(c)[14] = (a)[12] * (b)[2] + (a)[13] * (b)[5] + (a)[14] * (b)[8];
		(c)[15] = (a)[15];
	}
#define MatrixMultiply4X3d MatrixMultiply4X3<double>
#define MatrixMultiply4X3f MatrixMultiply4X3<float>


	// w = A * v
#define MATRIX_VECTOR_MULTIPLY4X4(A,v,w)\
	{(w)[0] = (A)[0]*(v)[0]+(A)[1]*(v)[1]+(A)[2]*(v)[2]+(A)[3]*(v)[3];\
	(w)[1] = (A)[4]*(v)[0]+(A)[5]*(v)[1]+(A)[6]*(v)[2]+(A)[7]*(v)[3];\
	(w)[2] = (A)[8]*(v)[0]+(A)[9]*(v)[1]+(A)[10]*(v)[2]+(A)[11]*(v)[3];\
	(w)[3] = (A)[12]*(v)[0]+(A)[13]*(v)[1]+(A)[14]*(v)[2]+(A)[15]*(v)[3];}

	template <class TYPE>
	inline void MatrixVectorMultiply4X4(TYPE *A, TYPE *v, TYPE *w)
	{
		(w)[0] = (A)[0] * (v)[0] + (A)[1] * (v)[1] + (A)[2] * (v)[2] + (A)[3] * (v)[3];
		(w)[1] = (A)[4] * (v)[0] + (A)[5] * (v)[1] + (A)[6] * (v)[2] + (A)[7] * (v)[3];
		(w)[2] = (A)[8] * (v)[0] + (A)[9] * (v)[1] + (A)[10] * (v)[2] + (A)[11] * (v)[3];
		(w)[3] = (A)[12] * (v)[0] + (A)[13] * (v)[1] + (A)[14] * (v)[2] + (A)[15] * (v)[3];
	}
#define MatrixVectorMultiply4X4d MatrixVectorMultiply4X4<double>
#define MatrixVectorMultiply4X4f MatrixVectorMultiply4X4<float>


	// w = A^T * v
#define MATRIX_VECTOR_MULTIPLY4X4T(A,v,w)\
	{(w)[0] = (A)[0]*(v)[0]+(A)[4]*(v)[1]+(A)[8]*(v)[2]+(A)[12]*(v)[3];\
	(w)[1] = (A)[1]*(v)[0]+(A)[5]*(v)[1]+(A)[9]*(v)[2]+(A)[13]*(v)[3];\
	(w)[2] = (A)[2]*(v)[0]+(A)[6]*(v)[1]+(A)[10]*(v)[2]+(A)[14]*(v)[3];\
	(w)[3] = (A)[3]*(v)[0]+(A)[7]*(v)[1]+(A)[11]*(v)[2]+(A)[15]*(v)[3];}

	template <class TYPE>
	inline void MatrixVectorMultiply4X4T(TYPE *A, TYPE *v, TYPE *w)
	{
		(w)[0] = (A)[0] * (v)[0] + (A)[4] * (v)[1] + (A)[8] * (v)[2] + (A)[12] * (v)[3];
		(w)[1] = (A)[1] * (v)[0] + (A)[5] * (v)[1] + (A)[9] * (v)[2] + (A)[13] * (v)[3];
		(w)[2] = (A)[2] * (v)[0] + (A)[6] * (v)[1] + (A)[10] * (v)[2] + (A)[14] * (v)[3];
		(w)[3] = (A)[3] * (v)[0] + (A)[7] * (v)[1] + (A)[11] * (v)[2] + (A)[15] * (v)[3];
	}
#define MatrixVectorMultiply4X4Td MatrixVectorMultiply4X4T<double>
#define MatrixVectorMultiply4X4Tf MatrixVectorMultiply4X4T<float>


// build translation matrix
#define MATRIX_FROM_TRANSLATION4X4(a, v)\
	{	(a)[0] = 1;(a)[1] = 0;(a)[2] = 0;(a)[3] = (v)[0];\
	(a)[4] = 0;(a)[5] = 1;(a)[6] = 0;(a)[7] = (v)[1];\
	(a)[8] = 0;(a)[9] = 0;(a)[10] = 1;(a)[11] = (v)[2];\
	(a)[12] = 0;(a)[13] = 0;(a)[14] = 0;(a)[15] = 1;}

	template <class TYPE>
	inline void MatrixFromTranslation4X4(TYPE *a, TYPE *v)
	{
		(a)[0] = 1; (a)[1] = 0; (a)[2] = 0; (a)[3] = (v)[0];
		(a)[4] = 0; (a)[5] = 1; (a)[6] = 0; (a)[7] = (v)[1];
		(a)[8] = 0; (a)[9] = 0; (a)[10] = 1; (a)[11] = (v)[2];
		(a)[12] = 0; (a)[13] = 0; (a)[14] = 0; (a)[15] = 1;
	}
#define MatrixFromTranslation4X4d MatrixFromTranslation4X4<double>
#define MatrixFromTranslation4X4f MatrixFromTranslation4X4<float>


#define MATRIX_FROM_TRANSLATION_INV4X4(a, v)\
			{	(a)[0] = 1;(a)[1] = 0;(a)[2] = 0;(a)[3] = -(v)[0];\
	(a)[4] = 0;(a)[5] = 1;(a)[6] = 0;(a)[7] = -(v)[1];\
	(a)[8] = 0;(a)[9] = 0;(a)[10] = 1;(a)[11] = -(v)[2];\
	(a)[12] = 0;(a)[13] = 0;(a)[14] = 0;(a)[15] = 1;}

	template <class TYPE>
	inline void MatrixFromTranslationInv4X4(TYPE *a, TYPE *v)
	{
		(a)[0] = 1; (a)[1] = 0; (a)[2] = 0; (a)[3] = -(v)[0];
			(a)[4] = 0; (a)[5] = 1; (a)[6] = 0; (a)[7] = -(v)[1];
			(a)[8] = 0; (a)[9] = 0; (a)[10] = 1; (a)[11] = -(v)[2];
			(a)[12] = 0; (a)[13] = 0; (a)[14] = 0; (a)[15] = 1;
	}
#define MatrixFromTranslationInv4X4d MatrixFromTranslationInv4X4<double>
#define MatrixFromTranslationInv4X4f MatrixFromTranslationInv4X4<float>



	// extend a 3x3 matrix to 4x4
#define MATRIX_FROM_3X3_MATRIX4X4(a, m)\
			{	(a)[0] = (m)[0];(a)[1] = (m)[1];(a)[2] = (m)[2];(a)[3] = 0;\
	(a)[4] = (m)[3];(a)[5] = (m)[4];(a)[6] = (m)[5];(a)[7] = 0;\
	(a)[8] = (m)[6];(a)[9] = (m)[7];(a)[10] = (m)[8];(a)[11] = 0;\
	(a)[12] = 0;(a)[13] = 0;(a)[14] = 0;(a)[15] = 1;}

	template <class TYPE>
	inline void MatrixFrom3X3Matrix4X4(TYPE *a, TYPE *m)
	{
		(a)[0] = (m)[0]; (a)[1] = (m)[1]; (a)[2] = (m)[2]; (a)[3] = 0;
		(a)[4] = (m)[3]; (a)[5] = (m)[4]; (a)[6] = (m)[5]; (a)[7] = 0;
		(a)[8] = (m)[6]; (a)[9] = (m)[7]; (a)[10] = (m)[8]; (a)[11] = 0;
		(a)[12] = 0; (a)[13] = 0; (a)[14] = 0; (a)[15] = 1;
	}
#define MatrixFrom3X3Matrix4X4d MatrixFrom3X3Matrix4X4<double>
#define MatrixFrom3X3Matrix4X4f MatrixFrom3X3Matrix4X4<float>


// a=translation(v)*a
#define MATRIX_TRANSLATE4X4(a, v)\
		{(a)[3]+=(a)[0]*(v)[0]+(a)[1]*(v)[1]+(a)[2]*(v)[2];\
	(a)[7]+=(a)[4]*(v)[0]+(a)[5]*(v)[1]+(a)[6]*(v)[2];\
	(a)[11]+=(a)[8]*(v)[0]+(a)[9]*(v)[1]+(a)[10]*(v)[2];\
	(a)[15]+=(a)[12]*(v)[0]+(a)[13]*(v)[1]+(a)[14]*(v)[2];}

	template <class TYPE>
	inline void MatrixTranslate4X4(TYPE *a, TYPE *v)
	{
		(a)[3] += (a)[0] * (v)[0] + (a)[1] * (v)[1] + (a)[2] * (v)[2];
		(a)[7] += (a)[4] * (v)[0] + (a)[5] * (v)[1] + (a)[6] * (v)[2];
		(a)[11] += (a)[8] * (v)[0] + (a)[9] * (v)[1] + (a)[10] * (v)[2];
		(a)[15] += (a)[12] * (v)[0] + (a)[13] * (v)[1] + (a)[14] * (v)[2];
	}
#define MatrixTranslate4X4d MatrixTranslate4X4<double>
#define MatrixTranslate4X4f MatrixTranslate4X4<float>



// a=translation(-v)*a
#define MATRIX_TRANSLATE_INV4X4(a, v)\
	{(a)[3]-=((a)[0]*(v)[0]+(a)[1]*(v)[1]+(a)[2]*(v)[2]);\
	(a)[7]-=((a)[4]*(v)[0]+(a)[5]*(v)[1]+(a)[6]*(v)[2]);\
	(a)[11]-=((a)[8]*(v)[0]+(a)[9]*(v)[1]+(a)[10]*(v)[2]);\
	(a)[15]-=((a)[12]*(v)[0]+(a)[13]*(v)[1]+(a)[14]*(v)[2]);}

	template <class TYPE>
	inline void MatrixTranslateInv4X4(TYPE *a, TYPE *v)
	{
		(a)[3] -= (a)[0] * (v)[0] + (a)[1] * (v)[1] + (a)[2] * (v)[2];
		(a)[7] -= (a)[4] * (v)[0] + (a)[5] * (v)[1] + (a)[6] * (v)[2];
		(a)[11] -= (a)[8] * (v)[0] + (a)[9] * (v)[1] + (a)[10] * (v)[2];
		(a)[15] -= (a)[12] * (v)[0] + (a)[13] * (v)[1] + (a)[14] * (v)[2];
	}
#define MatrixTranslateInv4X4d MatrixTranslateInv4X4<double>
#define MatrixTranslateInv4X4f MatrixTranslateInv4X4<float>
}

#endif