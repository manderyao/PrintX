/***********************************************************************
* Copyright (c) 2016, Miaojun Yao
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of Carnegie Mellon University, nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
***********************************************************************/
// MY_MATH library.
// Contents:
//		Min and Max functions, float2integer functions, random functions
//		Vector and matrix functions
//**************************************************************************************
// Quaternion library.

// Quaternion implementation, simple operations: convert from/to rotation, slerp...
//**************************************************************************************


#ifndef __PRINTX_QUATERNION_H__
#define __PRINTX_QUATERNION_H__

#include <math.h>

namespace PRINTX_LIB
{
class Quaternion
{
public:
	double w, x, y, z;

public:
	Quaternion(){ x = 0; y = 0; z = 0; w = 1; }

	// direct construction
	Quaternion(double w_, double x_, double y_, double z_);
	void init(double w_, double x_, double y_, double z_);

	// construct quaternion from 3 euler angles (in the order of Z->Y->X), in radians
	Quaternion(double theta_z, double theta_y, double theta_x);
	void init(double theta_z, double theta_y, double theta_x);

	// construct quaternion from rotation matrix
	void fromMatrix3(double *m);

	// convert to a 3D rotation matrix
	void toMatrix3(double* m);

	// construct quaternion from rotating between two vectors
	void fromRotate(double *v1, double *v2);		// from v1 to v2

	// construct quaternion from rotating around axis
	void fromAngleAxis(double *axis, double angle);

	// normalize
	void normalize();

	// length
	inline double length() {return sqrt(w*w + x*x + y*y + z*z);}
	inline double lengthSqr() { return (w*w + x*x + y*y + z*z); }

	// convert to a 4D matrix(transposed), translate part is 0 
	void toMatrix4T(double* m);

	//	rotate the quaternion by angular velocity W and time t
	void updateByRot(double *w, double t);

	//	rotate a vector
	void rotateVec(double *v, double *v_out);

	//	slerp between two quaternion
	static Quaternion slerp(Quaternion& q1, Quaternion& q2, double t);

	//	angle between two quaternion
	static double getAngle(Quaternion& q1, Quaternion& q2);

	inline Quaternion operator* (Quaternion & q2)
	{		
		return Quaternion(
			w * q2.w - x * q2.x - y    * q2.y - z * q2.z,
			w * q2.x + q2.w * x + y    * q2.z - q2.y * z,
			w * q2.y + q2.w * y + q2.x * z - x    * q2.z,
			w * q2.z + q2.w * z + x    * q2.y - q2.x * y);
	}
	inline Quaternion operator/ (Quaternion & q2)
	{
		// compute invQ2 = q2^{-1}
		Quaternion invQ2;
		double invNorm2 = 1.0 / q2.lengthSqr();
		invQ2.w = q2.w * invNorm2;
		invQ2.x = -q2.x * invNorm2;
		invQ2.y = -q2.y * invNorm2;
		invQ2.z = -q2.z * invNorm2;

		// result = *this * invQ2
		return (*this * invQ2);
	}
};


}

#endif