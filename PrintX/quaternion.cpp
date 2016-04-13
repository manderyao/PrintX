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

#include "quaternion.h"
#include "vec3.h"
using namespace PRINTX_LIB;

Quaternion::Quaternion(double w_, double x_, double y_, double z_)
	:w(w_), x(x_), y(y_), z(z_)
{
	normalize();
}

Quaternion::Quaternion(double theta_z, double theta_y, double theta_x)
{
	double cos_z_2 = cos(0.5*theta_z);
	double cos_y_2 = cos(0.5*theta_y);
	double cos_x_2 = cos(0.5*theta_x);

	double sin_z_2 = sin(0.5*theta_z);
	double sin_y_2 = sin(0.5*theta_y);
	double sin_x_2 = sin(0.5*theta_x);

	// and now compute quaternion
	w = cos_z_2*cos_y_2*cos_x_2 + sin_z_2*sin_y_2*sin_x_2;
	x = cos_z_2*cos_y_2*sin_x_2 - sin_z_2*sin_y_2*cos_x_2;
	y = cos_z_2*sin_y_2*cos_x_2 + sin_z_2*cos_y_2*sin_x_2;
	z = sin_z_2*cos_y_2*cos_x_2 - cos_z_2*sin_y_2*sin_x_2;

	normalize();
}

void Quaternion::init(double w_, double x_, double y_, double z_)
{
	w = w_;
	x = x_;
	y = y_;
	z = z_;
	normalize();
}
void Quaternion::init(double theta_z, double theta_y, double theta_x)
{
	double cos_z_2 = cos(0.5*theta_z);
	double cos_y_2 = cos(0.5*theta_y);
	double cos_x_2 = cos(0.5*theta_x);

	double sin_z_2 = sin(0.5*theta_z);
	double sin_y_2 = sin(0.5*theta_y);
	double sin_x_2 = sin(0.5*theta_x);

	// and now compute quaternion
	w = cos_z_2*cos_y_2*cos_x_2 + sin_z_2*sin_y_2*sin_x_2;
	x = cos_z_2*cos_y_2*sin_x_2 - sin_z_2*sin_y_2*cos_x_2;
	y = cos_z_2*sin_y_2*cos_x_2 + sin_z_2*cos_y_2*sin_x_2;
	z = sin_z_2*cos_y_2*cos_x_2 - cos_z_2*sin_y_2*sin_x_2;

	normalize();
}

void Quaternion::normalize()
{
	double len=length();
	w/=len;
	x/=len;
	y/=len;
	z/=len;
}

void Quaternion::toMatrix4T(double* m)
{
	m[0] = 1 - 2 * (y*y + z*z); m[4] = 2 * (x*y - w*z); m[8] = 2 * (x*z + w*y); m[12] = 0;

	m[1] = 2 * (x*y + w*z); m[5] = 1 - 2 * (x*x + z*z); m[9] = 2 * (y*z - w*x); m[13] = 0;

	m[2] = 2 * (x*z - w*y); m[6] = 2 * (y*z + w*x); m[10] = 1 - 2 * (x*x + y*y); m[14] = 0;

	m[3] = 0; m[7] = 0; m[11] = 0; m[15] = 1;
}
void Quaternion::toMatrix3(double* m)
{
	m[0]=1-2*(y*y+z*z);
	m[1]=2*(x*y-w*z);
	m[2]=2*(x*z+w*y);   
	m[3]=2*(x*y+w*z);
	m[4]=1-2*(x*x+z*z);
	m[5]=2*(y*z-w*x);
	m[6]=2*(x*z-w*y);
	m[7]=2*(y*z+w*x);
	m[8]=1-2*(x*x+y*y);
}

void Quaternion::fromAngleAxis(double *axis, double angle)
{
	double halfAngle=0.5*angle;
	double fSin = sin(halfAngle);
	w = cos(halfAngle);
	x = fSin*axis[0];
	y = fSin*axis[1];
	z = fSin*axis[2];
}

void Quaternion::fromRotate(double *v1, double *v2)		// from v1 to v2
{
	//normalize

	double d = VectorDotProduct3d(v1, v2);

	// If dot == 1, vectors are the same
	if (d >= 1.0f)
	{
		x = 0; y = 0; z = 0; w = 1;
	}
	if (d < (1e-6f - 1.0f))
	{
			// Generate an axis
		double axis[3];
		double ref[3] = { 1, 0, 0 };
		VectorCrossProduct3d(ref, v1, axis);
		if (VectorLength3d(axis) < 1e-12)
		{
			double ref2[3] = { 0, 1, 0 };
			VectorCrossProduct3d(ref2, v1, axis);
		}

		VectorNormalize3d(axis, axis);
#ifndef M_PI
#define M_PI 3.14159265358979
#endif
		fromAngleAxis(axis, M_PI);
	}
	else
	{
		double s = sqrt((1 + d) * 2);
		double invs = 1 / s;

		double c[3];
		VectorCrossProduct3d(v1, v2, c);

		x = c[0] * invs;
		y = c[1] * invs;
		z = c[2] * invs;
		w = s * 0.5f;
		normalize();
	}
}

void Quaternion::updateByRot(double *W, double t)
{
	//Update angular status.
	float temp_q[4];
	temp_q[0] = -W[0] * x - W[1] * y - W[2] * z;
	temp_q[1] = w * W[0] + W[1] * z - W[2] * y;
	temp_q[2] = w * W[1] + W[2] * x - W[0] * z;
	temp_q[3] = w * W[2] + W[0] * y - W[1] * x;
	w += 0.5*t*temp_q[0];
	x += 0.5*t*temp_q[1];
	y += 0.5*t*temp_q[2];
	z += 0.5*t*temp_q[3];
	//Normalize the quaternion
	normalize();
}

void Quaternion::rotateVec(double *v, double *v_out)
{
	// nVidia SDK implementation
	double uv[3], uuv[3];
	double qvec[3] = { x, y, z };
	VectorCrossProduct3d(qvec, v, uv);
	VectorCrossProduct3d(qvec, uv, uuv);
	VectorScale3d(uv, (2.0 * w));
	VectorScale3d(uuv, 2.0);

	VectorAdd3d(v, uv, v_out);
	VectorAdd3d(v_out, uuv, v_out);
}

void  Quaternion::fromMatrix3(double *a)
{
	double trace = a[0] + a[4] + a[8]; // I removed + 1.0f; see discussion with Ethan
	if (trace > 0) {// I changed M_EPSILON to 0
		double s = 0.5 / sqrt(trace + 1.0);
		w = 0.25f / s;
		x = (a[7] - a[5]) * s;
		y = (a[2] - a[6]) * s;
		z = (a[3] - a[1]) * s;
	}
	else {
		if (a[0] > a[4] && a[0] > a[8]) {
			double s = 2.0 * sqrt(1.0f + a[0] - a[4] - a[8]);
			w = (a[7] - a[5]) / s;
			x = 0.25 * s;
			y = (a[1] + a[3]) / s;
			z = (a[2] + a[6]) / s;
		}
		else if (a[4] > a[8]) {
			double s = 2.0 * sqrt(1.0f + a[4] - a[0] - a[8]);
			w = (a[2] - a[6]) / s;
			x = (a[1] + a[3]) / s;
			y = 0.25 * s;
			z = (a[5] + a[7]) / s;
		}
		else {
			double s = 2.0 * sqrt(1.0 + a[8] - a[0] - a[4]);
			w = (a[3] - a[1]) / s;
			x = (a[2] + a[6]) / s;
			y = (a[5] + a[7]) / s;
			z = 0.25 * s;
		}
	}

	normalize();
}

// Spherical linear interpolation between unit quaternions q1 and q2 with interpolation parameter t.
Quaternion Quaternion::slerp(Quaternion& q1, Quaternion& q2, double t)
{
	double w1, x1, y1, z1, w2, x2, y2, z2, w3, x3, y3, z3;
	Quaternion q2New;
	double theta, mult1, mult2;

	w1 = q1.w; x1 = q1.x; y1 = q1.y; z1 = q1.z;
	w2 = q2.w; x2 = q2.x; y2 = q2.y; z2 = q2.z;

	// Reverse the sign of q2 if q1.q2 < 0.
	if (w1*w2 + x1*x2 + y1*y2 + z1*z2 < 0)
	{
		w2 = -w2; x2 = -x2; y2 = -y2; z2 = -z2;
	}

	theta = acos(w1*w2 + x1*x2 + y1*y2 + z1*z2);

	if (theta > 0.000001)
	{
		mult1 = sin((1 - t)*theta) / sin(theta);
		mult2 = sin(t*theta) / sin(theta);
	}

	// To avoid division by 0 and by very small numbers the approximation of sin(angle)
	// by angle is used when theta is small (0.000001 is chosen arbitrarily).
	else
	{
		mult1 = 1 - t;
		mult2 = t;
	}

	w3 = mult1*w1 + mult2*w2;
	x3 = mult1*x1 + mult2*x2;
	y3 = mult1*y1 + mult2*y2;
	z3 = mult1*z1 + mult2*z2;

	Quaternion q(w3, x3, y3, z3);
	q.normalize();
	return q;
}

// Spherical linear interpolation between unit quaternions q1 and q2 with interpolation parameter t.
double Quaternion::getAngle(Quaternion& q1, Quaternion& q2)
{
	double w1, x1, y1, z1, w2, x2, y2, z2;
	Quaternion q2New;

	w1 = q1.w; x1 = q1.x; y1 = q1.y; z1 = q1.z;
	w2 = q2.w; x2 = q2.x; y2 = q2.y; z2 = q2.z;

	// Reverse the sign of q2 if q1.q2 < 0.
	if (w1*w2 + x1*x2 + y1*y2 + z1*z2 < 0)
	{
		w2 = -w2; x2 = -x2; y2 = -y2; z2 = -z2;
	}

	return acos(w1*w2 + x1*x2 + y1*y2 + z1*z2);
}