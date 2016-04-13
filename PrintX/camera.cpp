/*

* Original work Copyright (c) 2007, Carnegie Mellon University
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

Code author: Jernej Barbic, 2003-2007
Research: Jernej Barbic, Doug L. James
Version: 1.0

A "spherical" OpenGL camera class, with the ability to set the OpenGL modelview transformation matrix.

A "spherical" camera is a camera that is located at location (R, Phi, Theta), in spherical coordinates, away from some focus position. It is pointed towards the focus position. It is useful for orbiting a fixed location (focus position) in interactive applications.

cameraPosition = focusPosition +
[ R * cos(Phi) * cos (Theta), R * sin(Theta), -R * sin(Phi) * cos (Theta) ]

*/

#include "vec3.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "camera.h"

#ifndef M_PI
#define M_PI 3.14159265358979
#endif


SphericalCamera::SphericalCamera(float R, float Phi, float Theta, float * focusPosition_, float * up_, float movementSensitivity_, float camera2WorldScalingFactor_)
	: R0(R), Phi0(Phi), Theta0(Theta), movementSensitivity(movementSensitivity_), camera2WorldScalingFactor(camera2WorldScalingFactor_)
{
  focusPosition[0] = focusPosition_[0]; 
  focusPosition[1] = focusPosition_[1]; 
  focusPosition[2] = focusPosition_[2];
  up[0] = up_[0]; 
  up[1] = up_[1]; 
  up[2] = up_[2];
  origin[0] = 0.0; 
  origin[1] = 0.0; 
  origin[2] = 0.0;

	//default projection
  isPersp = false;
  ortho_left = -100;
  ortho_right = 100;
  ortho_bottom = -100;
  ortho_top = 100;
  Reset();
}

void SphericalCamera::setPersp(float fovy_, float aspect_, float zNear_, float zFar_)
{
	fovy = fovy_;
	aspect = aspect_;
	zNear = zNear_;
	zFar = zFar_;
	isPersp = true;
}

void SphericalCamera::setOrtho(float left, float right, float bottom, float top, float near, float far)
{
	ortho_left = left;
	ortho_right = right;
	ortho_bottom = bottom;
	ortho_top = top;
	ortho_near = near;
	ortho_far = far;
	isPersp = false;
}

void SphericalCamera::MoveRight(float amount)
{
  Phi += amount * movementSensitivity;
  ComputeCameraPosition();
  ComputeLocalCoordinateSystem();
}

void SphericalCamera::MoveUp(float amount)
{
  Theta += amount * movementSensitivity;

  if (Theta > 89.0 * M_PI / 180)
    Theta = 89.0 * M_PI / 180;

  if (Theta < -89.0 * M_PI / 180)
    Theta = -89.0 * M_PI / 180;

  ComputeCameraPosition();
  ComputeLocalCoordinateSystem();
}

void SphericalCamera::ZoomIn(float amount)
{
  R -= amount * movementSensitivity;
  
/*
  if (R < fabs(movementSensitivity))
    R = fabs(movementSensitivity);
*/

  ComputeCameraPosition();
  ComputeLocalCoordinateSystem();
}

void SphericalCamera::MoveIn(float amount)
{
  focusPosition[0] +=  amount * zAxis[0];
  focusPosition[1] +=  amount * zAxis[1];
  focusPosition[2] +=  amount * zAxis[2];
  ComputeCameraPosition();
  ComputeLocalCoordinateSystem();
}

void SphericalCamera::SetDefault()
{
  R0 = R;
  Phi0 = Phi;
  Theta0 = Theta;
}

void SphericalCamera::Reset() // resets position to defaults
{
  R = R0;
  Phi = Phi0;
  Theta = Theta0;
  ComputeCameraPosition();
  ComputeLocalCoordinateSystem();
}

void SphericalCamera::SetFocusPosition(float focusPosition_[3])
{ 
  focusPosition[0] = focusPosition_[0];
  focusPosition[1] = focusPosition_[1];
  focusPosition[2] = focusPosition_[2];
  ComputeCameraPosition();
  ComputeLocalCoordinateSystem();
}

void SphericalCamera::SetPosition(float r, float phi, float theta) 
{ 
  R = r; 
  Phi = phi; 
  Theta = theta; 
  ComputeCameraPosition();
  ComputeLocalCoordinateSystem();
}

void SphericalCamera::Frustum(float *matrix, float left, float right, float bottom, float top, float znear, float zfar)
{
	GLfloat temp, temp2, temp3, temp4;
    temp = 2.0 * znear;
    temp2 = right - left;
    temp3 = top - bottom;
    temp4 = zfar - znear;
    matrix[0] = temp / temp2;
    matrix[1] = 0.0;
    matrix[2] = 0.0;
    matrix[3] = 0.0;
    matrix[4] = 0.0;
    matrix[5] = temp / temp3;
    matrix[6] = 0.0;
    matrix[7] = 0.0;
    matrix[8] = (right + left) / temp2;
    matrix[9] = (top + bottom) / temp3;
    matrix[10] = (-zfar - znear) / temp4;
    matrix[11] = -1.0;
    matrix[12] = 0.0;
    matrix[13] = 0.0;
    matrix[14] = (-temp * zfar) / temp4;
    matrix[15] = 0.0;
}

void SphericalCamera::project(float * m)
{
	if (isPersp)
	{
		if (!m)
		{
			gluPerspective(fovy, aspect, zNear, zFar);
			return;
		}
	
		float ymax, xmax;
		ymax = zNear * tanf(fovy * M_PI / 360.0);
		xmax = ymax * aspect;
		Frustum(m, -xmax, xmax, -ymax, ymax, zNear, zFar);
	}
	else
	{
		if (!m)
		{
			glOrtho(ortho_left, ortho_right, ortho_bottom, ortho_top, ortho_near, ortho_far);
			return;
		}

		//TODO: complete the matrix computation
	}
}

void SphericalCamera::Look(float * m)
{
	if (!m)
	{
		gluLookAt(cameraPosition[0], cameraPosition[1], cameraPosition[2], focusPosition[0], focusPosition[1], focusPosition[2], up[0], up[1], up[2]);
		return;
	}

	float x[3], y[3], z[3];
	float mag;

	/* Make rotation matrix */

	/* Z vector */
	VECTOR_SUBTRACT3(cameraPosition, focusPosition, z);
	VECTOR_NORMALIZE3(z, z, mag);

	/* Y vector */
	y[0] = up[0];
	y[1] = up[1];
	y[2] = up[2];

	/* X vector = Y cross Z */
	x[0] = y[1] * z[2] - y[2] * z[1];
	x[1] = -y[0] * z[2] + y[2] * z[0];
	x[2] = y[0] * z[1] - y[1] * z[0];

	/* Recompute Y = Z cross X */
	y[0] = z[1] * x[2] - z[2] * x[1];
	y[1] = -z[0] * x[2] + z[2] * x[0];
	y[2] = z[0] * x[1] - z[1] * x[0];

	VECTOR_NORMALIZE3(x, x, mag);

	VECTOR_NORMALIZE3(y, y, mag);

#define M(row,col)  m[col*4+row]
	M(0, 0) = x[0];
	M(0, 1) = x[1];
	M(0, 2) = x[2];
	M(0, 3) = -(x[0]*cameraPosition[0]+x[1]*cameraPosition[1]+x[2]*cameraPosition[2]);
	M(1, 0) = y[0];
	M(1, 1) = y[1];
	M(1, 2) = y[2];
	M(1, 3) = -(y[0]*cameraPosition[0]+y[1]*cameraPosition[1]+y[2]*cameraPosition[2]);
	M(2, 0) = z[0];
	M(2, 1) = z[1];
	M(2, 2) = z[2];
	M(2, 3) = -(z[0]*cameraPosition[0]+z[1]*cameraPosition[1]+z[2]*cameraPosition[2]);
	M(3, 0) = 0.0;
	M(3, 1) = 0.0;
	M(3, 2) = 0.0;
	M(3, 3) = 1.0;
#undef M


}

void SphericalCamera::ComputeCameraPosition()
{
  //cameraPosition = focusPosition + 
  //  Vec3f(R * cos(Phi) * cos (Theta), R * sin(Theta), -R * sin(Phi) * cos (Theta));

  cameraPosition[0] = focusPosition[0] + R * cos(Phi) * cos (Theta);
  cameraPosition[1] = focusPosition[1] + R * sin(Theta);
  cameraPosition[2] = focusPosition[2] - R * sin(Phi) * cos (Theta);
}

void SphericalCamera::GetAbsWorldPosition(float & x, float & y, float & z)
{
  x = cameraPosition[0];
  y = cameraPosition[1];
  z = cameraPosition[2];
}

// get observer's ear positions, assuming given ear separation distance
// observer is located at the camera and faces the focus position
void SphericalCamera::GetStereoPosition(float earSeparation,
            float * leftEarX, float * leftEarY, float * leftEarZ,
            float * rightEarX, float * rightEarY, float * rightEarZ)
{
  // leftEar = cameraPosition - 0.5 * earSeparation * xAxis
  // rightEar = cameraPosition + 0.5 * earSeparation * xAxis

  *leftEarX = cameraPosition[0] - 0.5 * earSeparation * xAxis[0];
  *leftEarY = cameraPosition[1] - 0.5 * earSeparation * xAxis[1];
  *leftEarZ = cameraPosition[2] - 0.5 * earSeparation * xAxis[2];

  *rightEarX = cameraPosition[0] + 0.5 * earSeparation * xAxis[0];
  *rightEarY = cameraPosition[1] + 0.5 * earSeparation * xAxis[1];
  *rightEarZ = cameraPosition[2] + 0.5 * earSeparation * xAxis[2];
}

void SphericalCamera::ComputeLocalCoordinateSystem()
{
  xAxis[0] = -R*sin(Phi)*cos(Theta);
  xAxis[1] = 0;
  xAxis[2] = -R*cos(Phi)*cos(Theta);

  yAxis[0] = -R*cos(Phi)*sin(Theta);
  yAxis[1] = R*cos(Theta);
  yAxis[2] = R*sin(Phi)*sin(Theta);

  zAxis[0] = cameraPosition[0] - focusPosition[0];
  zAxis[1] = cameraPosition[1] - focusPosition[1];
  zAxis[2] = cameraPosition[2] - focusPosition[2];

  float length;

  length = sqrt(xAxis[0] * xAxis[0] + xAxis[1] * xAxis[1] + xAxis[2] * xAxis[2]);
  xAxis[0] /= length;
  xAxis[1] /= length;
  xAxis[2] /= length;

  length = sqrt(yAxis[0] * yAxis[0] + yAxis[1] * yAxis[1] + yAxis[2] * yAxis[2]);
  yAxis[0] /= length;
  yAxis[1] /= length;
  yAxis[2] /= length;

  length = sqrt(zAxis[0] * zAxis[0] + zAxis[1] * zAxis[1] + zAxis[2] * zAxis[2]);
  zAxis[0] /= length;
  zAxis[1] /= length;
  zAxis[2] /= length;

  xAxis2D[0] = -R*sin(Phi);
  xAxis2D[1] = -R*cos(Phi);
  xAxis2D[2] = 0;

  yAxis2D[0] = R*cos(Phi);
  yAxis2D[1] = -R*sin(Phi);
  yAxis2D[2] = 0;

  length = sqrt(xAxis2D[0] * xAxis2D[0] + xAxis2D[1] * xAxis2D[1] + xAxis2D[2] * xAxis2D[2]);
  xAxis2D[0] /= length;
  xAxis2D[1] /= length;
  xAxis2D[2] /= length;

  length = sqrt(yAxis2D[0] * yAxis2D[0] + yAxis2D[1] * yAxis2D[1] + yAxis2D[2] * yAxis2D[2]);
  yAxis2D[0] /= length;
  yAxis2D[1] /= length;
  yAxis2D[2] /= length;
}

void SphericalCamera::Get2DAxes(float xAxis2D[2], float yAxis2D[2])
{
  xAxis2D[0] = this->xAxis2D[0];
  xAxis2D[1] = this->xAxis2D[1];

  yAxis2D[0] = this->yAxis2D[0];
  yAxis2D[1] = this->yAxis2D[1];
}

void SphericalCamera::Get3DAxes(float xAxis3D[3], float yAxis3D[3], float zAxis3D[3])
{
  ComputeLocalCoordinateSystem();
  memcpy(xAxis3D, this->xAxis, sizeof(float) * 3);
  memcpy(yAxis3D, this->yAxis, sizeof(float) * 3);
  memcpy(zAxis3D, this->zAxis, sizeof(float) * 3);
}

void SphericalCamera::CameraRotation2WorldRotation2D(float * c, float * w)
{
  // c is in row-major format
  // 0  1  2
  // 3  4  5
  // 6  7  8 

  // w = cameraMatrix * c
  w[0] = c[0] * xAxis2D[0] + c[6] * yAxis2D[0];
  w[3] = c[3];
  w[6] = c[0] * xAxis2D[1] + c[6] * yAxis2D[1];

  w[1] = c[1] * xAxis2D[0] + c[7] * yAxis2D[0];
  w[4] = c[4];
  w[7] = c[1] * xAxis2D[1] + c[7] * yAxis2D[1];

  w[2] = c[2] * xAxis2D[0] + c[8] * yAxis2D[0];
  w[5] = c[5];
  w[8] = c[2] * xAxis2D[1] + c[8] * yAxis2D[1];
}

void SphericalCamera::CameraVector2WorldVector2D(float * c, float * w)
{ 
  w[0] = origin[0] + ( c[0] * xAxis2D[0] + c[2] * yAxis2D[0] ) * camera2WorldScalingFactor;
  w[1] = origin[1] + c[1] * camera2WorldScalingFactor;
  w[2] = origin[2] + ( c[0] * xAxis2D[1] + c[2] * yAxis2D[1] ) * camera2WorldScalingFactor;
}

void SphericalCamera::CameraVector2WorldVector2D(float c0, float c1, float c2, float * w)
{ 
  w[0] = origin[0] + ( c0 * xAxis2D[0] + c2 * yAxis2D[0] ) * camera2WorldScalingFactor;
  w[1] = origin[1] + c1 * camera2WorldScalingFactor;
  w[2] = origin[2] + ( c0 * xAxis2D[1] + c2 * yAxis2D[1] ) * camera2WorldScalingFactor;
}

void SphericalCamera::CameraVector2WorldVector_OrientationOnly2D(float * c, float * w)
{ 
  w[0] = ( c[0] * xAxis2D[0] + c[2] * yAxis2D[0] ) * camera2WorldScalingFactor;
  w[1] = c[1] * camera2WorldScalingFactor;
  w[2] = ( c[0] * xAxis2D[1] + c[2] * yAxis2D[1] ) * camera2WorldScalingFactor;
}

void SphericalCamera::CameraVector2WorldVector_OrientationOnly2D(float c0, float c1, float c2, float * w)
{ 
  w[0] = ( c0 * xAxis2D[0] + c2 * yAxis2D[0] ) * camera2WorldScalingFactor;
  w[1] = c1 * camera2WorldScalingFactor;
  w[2] = ( c0 * xAxis2D[1] + c2 * yAxis2D[1] ) * camera2WorldScalingFactor;
}

void SphericalCamera::CameraVector2WorldVector_OrientationOnly3D(float c0, float c1, float c2, float * w)
{ 
  // cameraPosition = focusPosition + 
  //   (R * cos(Phi) * cos (Theta), R * sin(Theta), -R * sin(Phi) * cos (Theta));

  // w = c0 * xAxis + c1 * yAxis + c2 * zAxis
  w[0] = (c0 * xAxis[0] + c1 * yAxis[0] + c2 * zAxis[0]) * camera2WorldScalingFactor;
  w[1] = (c0 * xAxis[1] + c1 * yAxis[1] + c2 * zAxis[1]) * camera2WorldScalingFactor;
  w[2] = (c0 * xAxis[2] + c1 * yAxis[2] + c2 * zAxis[2]) * camera2WorldScalingFactor;
}

// does not scale the axes
void SphericalCamera::WorldVector2CameraVector_OrientationOnly2D(float * w, float * c)
{
  c[0] = w[0] * xAxis2D[0] + w[2] * xAxis2D[1];
  c[1] = w[1];
  c[2] = w[0] * yAxis2D[0] + w[2] * yAxis2D[1];
}

void SphericalCamera::WorldVector2CameraVector_OrientationOnly2D(float w0, float w1, float w2, float * c)
{
  c[0] = w0 * xAxis2D[0] + w2 * xAxis2D[1];
  c[1] = w1;
  c[2] = w0 * yAxis2D[0] + w2 * yAxis2D[1];
}

void SphericalCamera::SavePosition(char * filename)
{
  FILE * fout;
  fout = fopen(filename, "w");
  fprintf(fout,"%G %G %G\n", R, Phi, Theta);
  fprintf(fout,"%G %G %G\n", focusPosition[0], focusPosition[1], focusPosition[2]);
  fclose(fout);
}

void SphericalCamera::LoadPosition(char * filename)
{
  FILE * fin;
  fin = fopen(filename,"r");
  if (fin)
  {
    if (fscanf(fin,"%lf %lf %lf\n",&R,&Phi,&Theta) < 3)
      printf("Warning: bad file syntax. Unable to read camera parameters.\n");
    if (fscanf(fin,"%lf %lf %lf\n", &focusPosition[0], &focusPosition[1], &focusPosition[2]) < 3)
      printf("Warning: bad file syntax. Unable to read camera parameters.\n");
    fclose(fin);
  }
  else
  {
    printf("Warning: file %s does not exist. Camera position not modified.\n",filename);
  }
  ComputeCameraPosition();
  ComputeLocalCoordinateSystem();
}

void SphericalCamera::PushPosition()
{
  savedState.R = R;
  savedState.Phi = Phi;
  savedState.Theta = Theta;

  savedState.focusPosition[0] = focusPosition[0];
  savedState.focusPosition[1] = focusPosition[1];
  savedState.focusPosition[2] = focusPosition[2];
}

void SphericalCamera::PopPosition()
{
  R = savedState.R;
  Phi = savedState.Phi;
  Theta = savedState.Theta;

  focusPosition[0] = savedState.focusPosition[0];
  focusPosition[1] = savedState.focusPosition[1];
  focusPosition[2] = savedState.focusPosition[2];

  ComputeCameraPosition();
  ComputeLocalCoordinateSystem();
}

void SphericalCamera::DetermineCameraParameters( float centerX, float centerY, float centerZ, float modelRadius, float * focusX, float * focusY, float * focusZ, float * cameraRadius, float * zNear, float * zFar)
{
  *focusX = centerX;
  *focusY = centerY;
  *focusZ = centerZ;

  *cameraRadius = modelRadius * 3;
  
  *zNear = *cameraRadius * 0.01;
  *zFar = *cameraRadius * 100;
}

// converts 4x4 row-major camera-space transform matrix to 4x4 row-major world-space transform matrix
void SphericalCamera::CameraTransform2WorldTransform2D(float * c, float * w)
{
 // c is in row-major format

  // 0  1   2   3
  // 4  5   6   7
  // 8  9  10  11
  // 12 13 14  15

  // w = cameraMatrix * c
  w[0] = c[0] * xAxis2D[0] + c[8] * yAxis2D[0];
  w[4] = c[4];
  w[8] = c[0] * xAxis2D[1] + c[8] * yAxis2D[1];
  w[12] = 0;

  w[1] = c[1] * xAxis2D[0] + c[9] * yAxis2D[0];
  w[5] = c[5];
  w[9] = c[1] * xAxis2D[1] + c[9] * yAxis2D[1];
  w[13] = 0;

  w[2] = c[2] * xAxis2D[0] + c[10] * yAxis2D[0];
  w[6] = c[6];
  w[10] = c[2] * xAxis2D[1] + c[10] * yAxis2D[1];
  w[14] = 0;

  w[3] = origin[0] + ( c[3] * xAxis2D[0] + c[11] * yAxis2D[0] ) * camera2WorldScalingFactor;
  w[7] = origin[1] + c[7] * camera2WorldScalingFactor;
  w[11] = origin[2] + ( c[3] * xAxis2D[1] + c[11] * yAxis2D[1] ) * camera2WorldScalingFactor;
  w[15] = 1;
}

// converts 4x4 column-major camera-space transform matrix to 4x4 column-major world-space transform matrix
void SphericalCamera::CameraTransform2WorldTransform2D_ColumnMajor(float * c, float * w)
{
 // c is in column-major format

  // 0  4   8  12
  // 1  5   9  13
  // 2  6  10  14
  // 3  7  11  15

  // w = cameraMatrix * c
  w[0] = c[0] * xAxis2D[0] + c[2] * yAxis2D[0];
  w[1] = c[1];
  w[2] = c[0] * xAxis2D[1] + c[2] * yAxis2D[1];
  w[3] = 0;

  w[4] = c[4] * xAxis2D[0] + c[6] * yAxis2D[0];
  w[5] = c[5];
  w[6] = c[4] * xAxis2D[1] + c[6] * yAxis2D[1];
  w[7] = 0;

  w[8] = c[8] * xAxis2D[0] + c[10] * yAxis2D[0];
  w[9] = c[9];
  w[10] = c[8] * xAxis2D[1] + c[10] * yAxis2D[1];
  w[11] = 0;

  w[12] = origin[0] + ( c[12] * xAxis2D[0] + c[14] * yAxis2D[0] ) * camera2WorldScalingFactor;
  w[13] = origin[1] + c[13] * camera2WorldScalingFactor;
  w[14] = origin[2] + ( c[12] * xAxis2D[1] + c[14] * yAxis2D[1] ) * camera2WorldScalingFactor;
  w[15] = 1;
}

// converts 4x4 column-major camera-space transform matrix to 4x4 column-major world-space transform matrix
// without scaling!!!
void SphericalCamera::CameraTransform2WorldTransform2D_NoScaling_ColumnMajor(float * c, float * w)
{
 // c is in column-major format

  // 0  4   8  12
  // 1  5   9  13
  // 2  6  10  14
  // 3  7  11  15

  // w = cameraMatrix * c
  w[0] = c[0] * xAxis2D[0] + c[2] * yAxis2D[0];
  w[1] = c[1];
  w[2] = c[0] * xAxis2D[1] + c[2] * yAxis2D[1];
  w[3] = 0;

  w[4] = c[4] * xAxis2D[0] + c[6] * yAxis2D[0];
  w[5] = c[5];
  w[6] = c[4] * xAxis2D[1] + c[6] * yAxis2D[1];
  w[7] = 0;

  w[8] = c[8] * xAxis2D[0] + c[10] * yAxis2D[0];
  w[9] = c[9];
  w[10] = c[8] * xAxis2D[1] + c[10] * yAxis2D[1];
  w[11] = 0;

  w[12] = origin[0] + ( c[12] * xAxis2D[0] + c[14] * yAxis2D[0] );
  w[13] = origin[1] + c[13];
  w[14] = origin[2] + ( c[12] * xAxis2D[1] + c[14] * yAxis2D[1] );
  w[15] = 1;
}

void SphericalCamera::GetFocusPosition(float focusPosition_[3])
{ 
  focusPosition_[0] = focusPosition[0];
  focusPosition_[1] = focusPosition[1];
  focusPosition_[2] = focusPosition[2];
}

void SphericalCamera::MoveFocusRight(float amount)
{
  float focusPos[3];
  GetFocusPosition(focusPos);
  focusPos[0] += movementSensitivity*amount * xAxis[0];
  focusPos[1] += movementSensitivity*amount * xAxis[1];
  focusPos[2] += movementSensitivity*amount * xAxis[2];
  SetFocusPosition(focusPos);
}

void SphericalCamera::MoveFocusUp(float amount)
{
  float focusPos[3];
  GetFocusPosition(focusPos);
  focusPos[0] += movementSensitivity*amount * yAxis[0];
  focusPos[1] += movementSensitivity*amount * yAxis[1];
  focusPos[2] += movementSensitivity*amount * yAxis[2];
  SetFocusPosition(focusPos);
}

void SphericalCamera::MoveFocusTo(float* focusPos)
{
	SetFocusPosition(focusPos);
}


