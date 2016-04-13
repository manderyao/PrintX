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

#include "mat3.h"


#define MAT3_N 3

static double hypot2(double x, double y) {
	return sqrt(x*x + y*y);
}

// Symmetric Householder reduction to tridiagonal form.

static void tred2(double V[MAT3_N][MAT3_N], double d[MAT3_N], double e[MAT3_N]) {

	//  This is derived from the Algol procedures tred2 by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	for (int j = 0; j < MAT3_N; j++) {
		d[j] = V[MAT3_N - 1][j];
	}

	// Householder reduction to tridiagonal form.

	for (int i = MAT3_N - 1; i > 0; i--) {

		// Scale to avoid under/overflow.

		double scale = 0.0;
		double h = 0.0;
		for (int k = 0; k < i; k++) {
			scale = scale + fabs(d[k]);
		}
		if (scale == 0.0) {
			e[i] = d[i - 1];
			for (int j = 0; j < i; j++) {
				d[j] = V[i - 1][j];
				V[i][j] = 0.0;
				V[j][i] = 0.0;
			}
		}
		else {

			// Generate Householder vector.

			for (int k = 0; k < i; k++) {
				d[k] /= scale;
				h += d[k] * d[k];
			}
			double f = d[i - 1];
			double g = sqrt(h);
			if (f > 0) {
				g = -g;
			}
			e[i] = scale * g;
			h = h - f * g;
			d[i - 1] = f - g;
			for (int j = 0; j < i; j++) {
				e[j] = 0.0;
			}

			// Apply similarity transformation to remaining columns.

			for (int j = 0; j < i; j++) {
				f = d[j];
				V[j][i] = f;
				g = e[j] + V[j][j] * f;
				for (int k = j + 1; k <= i - 1; k++) {
					g += V[k][j] * d[k];
					e[k] += V[k][j] * f;
				}
				e[j] = g;
			}
			f = 0.0;
			for (int j = 0; j < i; j++) {
				e[j] /= h;
				f += e[j] * d[j];
			}
			double hh = f / (h + h);
			for (int j = 0; j < i; j++) {
				e[j] -= hh * d[j];
			}
			for (int j = 0; j < i; j++) {
				f = d[j];
				g = e[j];
				for (int k = j; k <= i - 1; k++) {
					V[k][j] -= (f * e[k] + g * d[k]);
				}
				d[j] = V[i - 1][j];
				V[i][j] = 0.0;
			}
		}
		d[i] = h;
	}

	// Accumulate transformations.

	for (int i = 0; i < MAT3_N - 1; i++) {
		V[MAT3_N - 1][i] = V[i][i];
		V[i][i] = 1.0;
		double h = d[i + 1];
		if (h != 0.0) {
			for (int k = 0; k <= i; k++) {
				d[k] = V[k][i + 1] / h;
			}
			for (int j = 0; j <= i; j++) {
				double g = 0.0;
				for (int k = 0; k <= i; k++) {
					g += V[k][i + 1] * V[k][j];
				}
				for (int k = 0; k <= i; k++) {
					V[k][j] -= g * d[k];
				}
			}
		}
		for (int k = 0; k <= i; k++) {
			V[k][i + 1] = 0.0;
		}
	}
	for (int j = 0; j < MAT3_N; j++) {
		d[j] = V[MAT3_N - 1][j];
		V[MAT3_N - 1][j] = 0.0;
	}
	V[MAT3_N - 1][MAT3_N - 1] = 1.0;
	e[0] = 0.0;
}

// Symmetric tridiagonal QL algorithm.

static void tql2(double V[MAT3_N][MAT3_N], double d[MAT3_N], double e[MAT3_N]) {

	//  This is derived from the Algol procedures tql2, by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	for (int i = 1; i < MAT3_N; i++) {
		e[i - 1] = e[i];
	}
	e[MAT3_N - 1] = 0.0;

	double f = 0.0;
	double tst1 = 0.0;
	double eps = pow(2.0, -52.0);
	for (int l = 0; l < MAT3_N; l++) {

		// Find small subdiagonal element

		tst1 = std::max(tst1, fabs(d[l]) + fabs(e[l]));
		int m = l;
		while (m < MAT3_N) {
			if (fabs(e[m]) <= eps*tst1) {
				break;
			}
			m++;
		}

		// If m == l, d[l] is an eigenvalue,
		// otherwise, iterate.

		if (m > l) {
			int iter = 0;
			do {
				iter = iter + 1;  // (Could check iteration count here.)

				// Compute implicit shift

				double g = d[l];
				double p = (d[l + 1] - g) / (2.0 * e[l]);
				double r = hypot2(p, 1.0);
				if (p < 0) {
					r = -r;
				}
				d[l] = e[l] / (p + r);
				d[l + 1] = e[l] * (p + r);
				double dl1 = d[l + 1];
				double h = g - d[l];
				for (int i = l + 2; i < MAT3_N; i++) {
					d[i] -= h;
				}
				f = f + h;

				// Implicit QL transformation.

				p = d[m];
				double c = 1.0;
				double c2 = c;
				double c3 = c;
				double el1 = e[l + 1];
				double s = 0.0;
				double s2 = 0.0;
				for (int i = m - 1; i >= l; i--) {
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = hypot2(p, e[i]);
					e[i + 1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i + 1] = h + s * (c * g + s * d[i]);

					// Accumulate transformation.

					for (int k = 0; k < MAT3_N; k++) {
						h = V[k][i + 1];
						V[k][i + 1] = s * V[k][i] + c * h;
						V[k][i] = c * V[k][i] - s * h;
					}
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;

				// Check for convergence.

			} while (fabs(e[l]) > eps*tst1);
		}
		d[l] = d[l] + f;
		e[l] = 0.0;
	}

	// Sort eigenvalues and corresponding vectors.

	for (int i = 0; i < MAT3_N - 1; i++) {
		int k = i;
		double p = d[i];
		for (int j = i + 1; j < MAT3_N; j++) {
			if (d[j] < p) {
				k = j;
				p = d[j];
			}
		}
		if (k != i) {
			d[k] = d[i];
			d[i] = p;
			for (int j = 0; j < MAT3_N; j++) {
				p = V[j][i];
				V[j][i] = V[j][k];
				V[j][k] = p;
			}
		}
	}
}

/* Eigen-decomposition for symmetric 3x3 real matrices.
Public domain, copied from the public domain Java library JAMA. */
/* Symmetric matrix A => eigenvectors in columns of V, corresponding
eigenvalues in d. */
void PRINTX_LIB::eigenDecompose3X3(double *A, double V[MAT3_N][MAT3_N], double d[MAT3_N])
{
	double e[MAT3_N];
	for (int i = 0; i < MAT3_N; i++) {
		for (int j = 0; j < MAT3_N; j++) {
			V[i][j] = A[i*MAT3_N+j];
		}
	}
	tred2(V, d, e);
	tql2(V, d, e);
}

void PRINTX_LIB::testEigenDecomposition3X3()
{
	double A[9] = { 1, 2, 3, 2, 3, 4, 3, 4, 5 };
	double V[3][3];
	double d[3];
	eigenDecompose3X3(A, V, d);
	std::cout << d[0] << " " << d[1] << " " << d[2] << "\n";
}



double PRINTX_LIB::matrixOneNorm3X3(double *A)
{
	double norm = 0.0;
	for (int i = 0; i<3; i++)
	{
		double columnAbsSum = fabs(A[i + 0]) + fabs(A[i + 3]) + fabs(A[i + 6]);
		if (columnAbsSum > norm)
			norm = columnAbsSum;
	}
	return norm;
}

double PRINTX_LIB::matrixInfNorm3X3(double *A)
{
	double norm = 0.0;
	for (int i = 0; i<3; i++)
	{
		double rowSum = fabs(A[3 * i + 0]) + fabs(A[3 * i + 1]) + fabs(A[3 * i + 2]);
		if (rowSum > norm)
			norm = rowSum;
	}
	return norm;
}

void PRINTX_LIB::SVD3X3(double* F, double* U, double* Sigma, double* V, double singularValue_eps, int modifiedSVD)
{
	// form F^T F and do eigendecomposition
	double F_[9];
	MATRIX_TRANSPOSE3X3(F, F_);
	double normalEq[9];
	MATRIX_MULTIPLY3X3(F_, F, normalEq);

	double eigenValues[3];
	double eigenVectors[3][3];

	eigenDecompose3X3(normalEq, eigenVectors, eigenValues);

	// descending order
	std::swap(eigenValues[0], eigenValues[2]);	
	V[0] = eigenVectors[0][2]; V[3] = eigenVectors[1][2]; V[6] = eigenVectors[2][2];
	V[1] = eigenVectors[0][1]; V[4] = eigenVectors[1][1]; V[7] = eigenVectors[2][1];
	V[2] = eigenVectors[0][0]; V[5] = eigenVectors[1][0]; V[8] = eigenVectors[2][0];

	// Handle situation:
	// 1. det(V) == -1
	//    - multiply the first column of V by -1
	if (MATRIX_DETERMINANT3X3(V) < 0.0)
	{
		// convert V into a rotation (multiply column 1 by -1)
		V[0] *= -1.0;
		V[3] *= -1.0;
		V[6] *= -1.0;
	}

	Sigma[0] = (eigenValues[0] > 0.0) ? sqrt(eigenValues[0]) : 0.0;
	Sigma[1] = (eigenValues[1] > 0.0) ? sqrt(eigenValues[1]) : 0.0;
	Sigma[2] = (eigenValues[2] > 0.0) ? sqrt(eigenValues[2]) : 0.0;

	//printf("--- Sigma ---\n");
	//printf("%G %G %G\n", Sigma[0][0], Sigma[1][1], Sigma[2][2]);

	// compute inverse of singular values
	// also check if singular values are close to zero
	double SigmaInverse[3];
	SigmaInverse[0] = (Sigma[0] > singularValue_eps) ? (1.0 / Sigma[0]) : 0.0;
	SigmaInverse[1] = (Sigma[1] > singularValue_eps) ? (1.0 / Sigma[1]) : 0.0;
	SigmaInverse[2] = (Sigma[2] > singularValue_eps) ? (1.0 / Sigma[2]) : 0.0;

	// compute U using the formula:
	// U = F * V * diag(SigmaInverse)
	MATRIX_MULTIPLY3X3(F, V, U);
	MATRIX_MULTIPLY_DIAG_RIGHT3X3(U, SigmaInverse, U);

	// In theory, U is now orthonormal, U^T U = U U^T = I .. it may be a rotation or a reflection, depending on F.
	// But in practice, if singular values are small or zero, it may not be orthonormal, so we need to fix it.
	// Handle situation:
	// 2. An entry of Sigma is near zero
	// ---------------------------------------------------------

	/*
	printf("--- SigmaInverse ---\n");
	SigmaInverse.print();
	printf(" --- U ---\n");
	U.print();
	*/

	if ((Sigma[0] < singularValue_eps) && (Sigma[1] < singularValue_eps) && (Sigma[2] < singularValue_eps))
	{
		// extreme case, all singular values are small, material has collapsed almost to a point
		// see [Irving 04], p. 4
		MATRIX_SET_IDENTITY3X3(U);
	}
	else if ((Sigma[1] < singularValue_eps) && (Sigma[2] < singularValue_eps))
	{
		// handle the case where two singular values are small, but the third one is not
		// handle it by computing two (arbitrary) vectors orthogonal to the eigenvector for the large singular value
		// only the column dimA can be trusted, columns dimB and dimC correspond to tiny singular values
		double tmpVec1[3] = { U[0], U[3], U[6] }; // column dimA
		double tmpVec2[3], tmpVec3[3];
		VectorFindOrthoNormal3d(tmpVec1, tmpVec2);
		VECTOR_CROSS_PRODUCT3(tmpVec1, tmpVec2, tmpVec3);
		double len;
		VECTOR_NORMALIZE3(tmpVec3, tmpVec3, len);
		U[1] = tmpVec2[0];
		U[4] = tmpVec2[1];
		U[7] = tmpVec2[2];
		U[2] = tmpVec3[0];
		U[5] = tmpVec3[1];
		U[8] = tmpVec3[2];
		if (MATRIX_DETERMINANT3X3(U) < 0.0)
		{
			U[1] *= -1.0;
			U[4] *= -1.0;
			U[7] *= -1.0;
		}
	}
	else if (Sigma[2] < singularValue_eps)
	{
		// handle the case where one singular value is small, but the other two are not
		// handle it by computing the cross product of the two eigenvectors for the two large singular values
		double tmpVec1[3] = { U[0], U[3], U[6] };
		double tmpVec2[3] = { U[1], U[4], U[7] };
		double tmpVec3[3];
		VECTOR_CROSS_PRODUCT3(tmpVec1, tmpVec2, tmpVec3);
		double len;
		VECTOR_NORMALIZE3(tmpVec3, tmpVec3, len);
		U[2] = tmpVec3[0];
		U[5] = tmpVec3[1];
		U[8] = tmpVec3[2];
		if (MATRIX_DETERMINANT3X3(U) < 0.0)
		{
			U[2] *= -1.0;
			U[5] *= -1.0;
			U[8] *= -1.0;
		}
	}
	else
	{
		// Handle situation:
		// 3. negative determinant (Tet is inverted in solid mechanics)
		//    - check if det(U) == -1
		//    - If yes, then negate the minimal element of Sigma
		//      and the corresponding column of U
		if (MATRIX_DETERMINANT3X3(U) < 0.0)
		{
			// negate the smallest singular value
			Sigma[2] *= -1.0;
			U[2] *= -1.0;
			U[5] *= -1.0;
			U[8] *= -1.0;
		}
	}
}

// Computes the Polar Decomposition of a general 3x3 matrix M.
// M = Q * S
// M is 3x3 input matrix
// Q is 3x3 orthogonal output matrix, Q Q^T = Q^T Q = I 
// S is 3x3 symmetric positive-definite output matrix
// Note: det(Q)=sgn(det(M)); this sign can be 1 or -1, depending on M
// if forceRotation is 1, the matrix Q will be a rotation, S will be symmetric, but not necessarily positive-definite
// M is not modified
// All matrices are row-major
double PRINTX_LIB::polarDecompose3X3(double * M, double * Q, double * S, int forceRotation, double tolerance)
{
	double Mk[9];
	double Ek[9];
	double det, M_oneNorm, M_infNorm, E_oneNorm;
	int useSVD = 0;

	// Mk = M^T
	for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			Mk[3 * i + j] = M[3 * j + i];

	M_oneNorm = matrixOneNorm3X3(Mk);
	M_infNorm = matrixInfNorm3X3(Mk);

	do
	{
		double MadjTk[9];

		// row 2 x row 3
		VectorCrossProduct3d(Mk + 3, Mk + 6, MadjTk);
		// row 3 x row 1
		VectorCrossProduct3d(Mk + 6, Mk, MadjTk + 3);
		// row 1 x row 2
		VectorCrossProduct3d(Mk, Mk + 3, MadjTk + 6);

		det = Mk[0] * MadjTk[0] + Mk[1] * MadjTk[1] + Mk[2] * MadjTk[2];

		if ((det <= 1e-6) && forceRotation)
		{
			useSVD = 1;
			break;
		}

		if (det == 0.0)
		{
			printf("Warning (polarDecomposition) : zero determinant encountered.\n");
			break;
		}

		double MadjT_one = matrixOneNorm3X3(MadjTk);
		double MadjT_inf = matrixInfNorm3X3(MadjTk);

		double gamma = sqrt(sqrt((MadjT_one * MadjT_inf) / (M_oneNorm * M_infNorm * det * det)));
		double g1 = gamma * 0.5;
		double g2 = 0.5 / (gamma * det);

		for (int i = 0; i<9; i++)
		{
			Ek[i] = Mk[i];
			Mk[i] = g1 * Mk[i] + g2 * MadjTk[i];
			Ek[i] -= Mk[i];
		}

		E_oneNorm = matrixOneNorm3X3(Ek);
		M_oneNorm = matrixOneNorm3X3(Mk);
		M_infNorm = matrixInfNorm3X3(Mk);
	} while (E_oneNorm > M_oneNorm * tolerance);

	if (useSVD)
	{
#if 0	// do nothing in this case
		// use the SVD algorithm to compute Q
		double modifiedSVD_singularValue_eps = tolerance;
		double Um[9], Vm[9], Vm_[9];
		double Lambda[3];
		int modifiedSVD = 1;
		SVD3X3(M, Um, Lambda, Vm, modifiedSVD_singularValue_eps, modifiedSVD);

		MATRIX_TRANSPOSE3X3(Vm, Vm_);
		MATRIX_MULTIPLY3X3(Um, Vm_, Q);
#endif
	}
	else
	{
		// Q = Mk^T 
		for (int i = 0; i<3; i++)
			for (int j = 0; j<3; j++)
				Q[3 * i + j] = Mk[3 * j + i];

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				S[3 * i + j] = 0;
				for (int k = 0; k < 3; k++)
					S[3 * i + j] += Mk[3 * i + k] * M[3 * k + j];
			}

		// S must be symmetric; enforce the symmetry
		S[1] = S[3] = 0.5 * (S[1] + S[3]);
		S[2] = S[6] = 0.5 * (S[2] + S[6]);
		S[5] = S[7] = 0.5 * (S[5] + S[7]);
	}


	return det;
}