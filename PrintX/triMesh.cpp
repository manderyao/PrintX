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
//  triangle mesh
//	when read a new mesh, it will realloc the memory automatically
//**************************************************************************************

#include "triMesh.h"
#include <map>
#include "vec3.h"
#include "mat3.h"
#include "mat4.h"
#include <iostream>
#include <fstream>
#include <queue>
#include <set>
#include <list>
#include <unordered_map>
#include "opgl.h"
#include "colorManager.h"
#include <sstream>
#include <unordered_set>

using namespace PRINTX_LIB;

TriMesh::TriMesh(const TriMesh &mesh)
{
	n_verts = mesh.n_verts;
	n_faces = mesh.n_faces;
	verts = (double*)malloc(sizeof(double) * 3 * n_verts);
	faces = (int*)malloc(sizeof(int) * 3 * n_faces);
	memcpy(verts, mesh.verts, sizeof(double) * 3 * n_verts);
	memcpy(faces, mesh.faces, sizeof(int) * 3 * n_faces);

	hasSegments = false;
	face_ids = 0;

	hasFaceNormals = false;
	hasVertNormals = false;
	faceNormals=0;
	vertNormals=0;

	hasVertColors = false;
	vert_colors = 0;
}

TriMesh::TriMesh()
{
	n_verts = 0;
	n_faces = 0;
	verts = 0;
	faces = 0;

	hasSegments = false;
	face_ids = 0;

	hasFaceNormals = false;
	hasVertNormals = false;
	faceNormals = 0;
	vertNormals = 0;

	hasVertColors = false;
	vert_colors = 0;
}

TriMesh::~TriMesh()
{
	free(verts);
	free(faces);
	free(face_ids);
	free(faceNormals);
	free(vertNormals);
	free(vert_colors);
}

void TriMesh::clone(const TriMesh &mesh)
{
	n_verts = mesh.n_verts;
	n_faces = mesh.n_faces;
	verts = (double*)realloc(verts, sizeof(double) * 3 * n_verts);
	faces = (int*)realloc(faces, sizeof(int) * 3 * n_faces);
	memcpy(verts, mesh.verts, sizeof(double) * 3 * n_verts);
	memcpy(faces, mesh.faces, sizeof(int) * 3 * n_faces);
	patches.clear();

	hasSegments = false;
	hasFaceNormals = false;
	hasVertNormals = false;
	hasVertColors = false;
}

void TriMesh::createRecPyramid(double *p, double *p0, double *p1, double *p2, double *p3)
{
	n_verts = 5;
	n_faces = 6;
	verts = (double*)realloc(verts, sizeof(double) * 3 * n_verts);
	faces = (int*)realloc(faces, sizeof(int) * 3 * n_faces);

	VectorSet3d(verts, p);
	VectorSet3d(verts+3, p0);
	VectorAdd3d(verts + 3, p, verts + 3);
	VectorSet3d(verts + 6, p1);
	VectorAdd3d(verts + 6, p, verts + 6);
	VectorSet3d(verts + 9, p2);
	VectorAdd3d(verts + 9, p, verts + 9);
	VectorSet3d(verts + 12, p3);
	VectorAdd3d(verts + 12, p, verts + 12);

	faces[0] = 0; faces[1] = 1; faces[2] = 2;
	faces[3] = 0; faces[4] = 2; faces[5] = 3;
	faces[6] = 0; faces[7] = 3; faces[8] = 4;
	faces[9] = 0; faces[10] = 4; faces[11] = 1;
	faces[12] = 1; faces[13] = 4; faces[14] = 3;
	faces[15] = 1; faces[16] = 3; faces[17] = 2;

	hasSegments = false;
	hasFaceNormals = false;
	hasVertNormals = false;
	hasVertColors = false;
}



void TriMesh::createHexPrism(double lx, double ly, double lz)
{
	double cube_verts[14 * 3] = { 0.500000, -1.000000, -0.866025, -0.500000, -1.000000, -0.866026, -1.000000, -1.000000, -0.000000, -0.500000, - 1.000000, 0.866025,
		0.500000, -1.000000, 0.866025, 1.000000, -1.000000, 0.000000, 0.500000, 1.000000, -0.866025, -0.500000, 1.000000, -0.866026, -1.000000, 1.000000, - 0.000000,
		-0.500000, 1.000000, 0.866025, 0.500000, 1.000000, 0.866025, 1.000000, 1.000000, 0.000000, 0.000000, - 1.000000, 0.000000, 0.000000, 1.000000, 0.000000 };
	int n_cube_verts = 14;
	int cube_faces[24 * 3] = { 0, 1, 7, 0, 7, 6, 1, 2, 8, 1, 8, 7, 2, 3, 9, 2, 9, 8, 3, 4, 10, 3, 10, 9, 4, 5, 11, 4, 11, 10, 5, 0, 6, 5, 6, 11, 
		1, 0, 12, 2, 1, 12, 3, 2, 12, 4, 3, 12, 5, 4, 12, 0, 5, 12, 6, 7, 13, 7, 8, 13, 8, 9, 13, 9, 10, 13, 10, 11, 13, 11, 6, 13 };
	int n_cube_faces = 24;

	n_verts = n_cube_verts;
	n_faces = n_cube_faces;
	verts = (double*)realloc(verts, sizeof(double) * 3 * n_verts);
	faces = (int*)realloc(faces, sizeof(int) * 3 * n_faces);

	memcpy(verts, cube_verts, sizeof(cube_verts[0])*n_verts * 3);
	memcpy(faces, cube_faces, sizeof(cube_faces[0])*n_faces * 3);

	scale(lx / 2, ly / 2, lz / 2);

	hasSegments = false;
	hasFaceNormals = false;
	hasVertNormals = false;
	hasVertColors = false;
}

void TriMesh::createOBBCube(double *pt, double *offx, double *offy, double *offz)
{
	createCube(1, 1, 1);
	translate(0.5, 0.5, 0.5);

	for (int i = 0; i < n_verts; i++)
	{
		double x = verts[3 * i];
		double y = verts[3 * i+1];
		double z = verts[3 * i+2];
		verts[3 * i] = x*offx[0] + y*offy[0] + z*offz[0];
		verts[3 * i+1] = x*offx[1] + y*offy[1] + z*offz[1];
		verts[3 * i+2] = x*offx[2] + y*offy[2] + z*offz[2];
	}

	translate(pt[0], pt[1], pt[2]);

	//hasSegments = false;
	//hasFaceNormals = false;
	//hasVertNormals = false;
	//hasVertColors = false;
}

void TriMesh::createFromCurv(std::vector<double> &verts_, double *offset, int n_off)
{
	int np = verts_.size() / 3;
	n_verts = n_off * np;
	n_faces = 2*(n_off-2) + (np - 1) * 2 *n_off;

	verts = (double*)realloc(verts, sizeof(double) * 3 * n_verts);
	faces = (int*)realloc(faces, sizeof(int) * 3 * n_faces);

	n_verts = 0;
	for (int i = 0; i < np; i++)
	{
		for (int j = 0; j < n_off; j++)
		{
			verts[3 * n_verts] = verts_[3 * i] + offset[3*j];
			verts[3 * n_verts + 1] = verts_[3 * i + 1] + offset[3 * j+1];
			verts[3 * n_verts + 2] = verts_[3 * i + 2] + offset[3 * j+2];
			n_verts++;
		}
	}

	n_faces = 0;
	for (int j = n_off-1; j > 1; j--)
	{
		faces[3 * n_faces] = 0;
		faces[3 * n_faces + 1] = j;
		faces[3 * n_faces + 2] = j-1;
		n_faces++;
	}
	for (int j = 1; j < n_off - 1; j++)
	{
		faces[3 * n_faces] = n_verts - n_off;
		faces[3 * n_faces + 1] = n_verts - n_off + j;
		faces[3 * n_faces + 2] = n_verts - n_off + j+1;
		n_faces++;
	}

	for (int i = 0; i < np - 1; i++)
	{
		int base = i * n_off;
		for (int j = 0; j < n_off-1; j++)
		{
			faces[3 * n_faces] = base+j;
			faces[3 * n_faces + 1] = base + j+1;
			faces[3 * n_faces + 2] = base + n_off+j + 1;
			n_faces++;

			faces[3 * n_faces] = base + j;
			faces[3 * n_faces + 1] = base + n_off + j + 1;
			faces[3 * n_faces + 2] = base + n_off + j;
			n_faces++;
		}
		faces[3 * n_faces] = base + n_off - 1;
		faces[3 * n_faces + 1] = base;
		faces[3 * n_faces + 2] = base + n_off;
		n_faces++;

		faces[3 * n_faces] = base + n_off - 1;
		faces[3 * n_faces + 1] = base + n_off;
		faces[3 * n_faces + 2] = base + 2*n_off -1;
		n_faces++;
	}

	hasSegments = false;
	hasFaceNormals = false;
	hasVertNormals = false;
	hasVertColors = false;
}

void TriMesh::createFromCurvVary(std::vector<double> &verts_, std::vector<double> & offset1, std::vector<double> &offset2, std::vector<double> &offset3, std::vector<double> &offset4)
{
	int np = verts_.size() / 3;
	n_verts = 4 * np;
	n_faces = 4 + (np - 1) * 8;

	verts = (double*)realloc(verts, sizeof(double) * 3 * n_verts);
	faces = (int*)realloc(faces, sizeof(int) * 3 * n_faces);

	n_verts = 0;
	for (int i = 0; i < np; i++)
	{
		verts[3 * n_verts] = verts_[3 * i] + offset1[3*i];
		verts[3 * n_verts + 1] = verts_[3 * i + 1] + offset1[3 * i+1];
		verts[3 * n_verts + 2] = verts_[3 * i + 2] + offset1[3 * i+2];
		n_verts++;

		verts[3 * n_verts] = verts_[3 * i] + offset2[3 * i];
		verts[3 * n_verts + 1] = verts_[3 * i + 1] + offset2[3 * i + 1];
		verts[3 * n_verts + 2] = verts_[3 * i + 2] + offset2[3 * i + 2];
		n_verts++;

		verts[3 * n_verts] = verts_[3 * i] + offset3[3 * i];
		verts[3 * n_verts + 1] = verts_[3 * i + 1] + offset3[3 * i + 1];
		verts[3 * n_verts + 2] = verts_[3 * i + 2] + offset3[3 * i + 2];
		n_verts++;

		verts[3 * n_verts] = verts_[3 * i] + offset4[3 * i];
		verts[3 * n_verts + 1] = verts_[3 * i + 1] + offset4[3 * i + 1];
		verts[3 * n_verts + 2] = verts_[3 * i + 2] + offset4[3 * i + 2];
		n_verts++;
	}

	n_faces = 0;
	faces[3 * n_faces] = 0;
	faces[3 * n_faces + 1] = 3;
	faces[3 * n_faces + 2] = 2;
	n_faces++;

	faces[3 * n_faces] = 0;
	faces[3 * n_faces + 1] = 2;
	faces[3 * n_faces + 2] = 1;
	n_faces++;

	faces[3 * n_faces] = n_verts - 4;
	faces[3 * n_faces + 1] = n_verts - 3;
	faces[3 * n_faces + 2] = n_verts - 2;
	n_faces++;

	faces[3 * n_faces] = n_verts - 4;
	faces[3 * n_faces + 1] = n_verts - 2;
	faces[3 * n_faces + 2] = n_verts - 1;
	n_faces++;

	for (int i = 0; i < np - 1; i++)
	{
		faces[3 * n_faces] = 4 * i;
		faces[3 * n_faces + 1] = 4 * i + 4;
		faces[3 * n_faces + 2] = 4 * i + 7;
		n_faces++;

		faces[3 * n_faces] = 4 * i;
		faces[3 * n_faces + 1] = 4 * i + 7;
		faces[3 * n_faces + 2] = 4 * i + 3;
		n_faces++;

		faces[3 * n_faces] = 4 * i + 1;
		faces[3 * n_faces + 1] = 4 * i + 5;
		faces[3 * n_faces + 2] = 4 * i + 4;
		n_faces++;

		faces[3 * n_faces] = 4 * i + 1;
		faces[3 * n_faces + 1] = 4 * i + 4;
		faces[3 * n_faces + 2] = 4 * i;
		n_faces++;

		faces[3 * n_faces] = 4 * i + 2;
		faces[3 * n_faces + 1] = 4 * i + 6;
		faces[3 * n_faces + 2] = 4 * i + 5;
		n_faces++;

		faces[3 * n_faces] = 4 * i + 2;
		faces[3 * n_faces + 1] = 4 * i + 5;
		faces[3 * n_faces + 2] = 4 * i + 1;
		n_faces++;

		faces[3 * n_faces] = 4 * i + 3;
		faces[3 * n_faces + 1] = 4 * i + 7;
		faces[3 * n_faces + 2] = 4 * i + 6;
		n_faces++;

		faces[3 * n_faces] = 4 * i + 3;
		faces[3 * n_faces + 1] = 4 * i + 6;
		faces[3 * n_faces + 2] = 4 * i + 2;
		n_faces++;
	}

	hasSegments = false;
	hasFaceNormals = false;
	hasVertNormals = false;
	hasVertColors = false;
}

void TriMesh::createFromCurv(std::vector<double> &verts_, double *offset1, double *offset2, double *offset3, double *offset4)
{
	int np = verts_.size() / 3;
	n_verts = 4*np;
	n_faces = 4 + (np - 1) * 8;

	verts = (double*)realloc(verts, sizeof(double) * 3 * n_verts);
	faces = (int*)realloc(faces, sizeof(int) * 3 * n_faces);

	n_verts = 0;
	for (int i = 0; i < np;i++)
	{
		verts[3 * n_verts] = verts_[3 * i] + offset1[0];
		verts[3 * n_verts + 1] = verts_[3 * i + 1] + offset1[1];
		verts[3 * n_verts + 2] = verts_[3 * i + 2] + offset1[2];
		n_verts++;

		verts[3 * n_verts] = verts_[3 * i] + offset2[0];
		verts[3 * n_verts + 1] = verts_[3 * i + 1] + offset2[1];
		verts[3 * n_verts + 2] = verts_[3 * i + 2] + offset2[2];
		n_verts++;

		verts[3 * n_verts] = verts_[3 * i] + offset3[0];
		verts[3 * n_verts + 1] = verts_[3 * i + 1] + offset3[1];
		verts[3 * n_verts + 2] = verts_[3 * i + 2] + offset3[2];
		n_verts++;

		verts[3 * n_verts] = verts_[3 * i] + offset4[0];
		verts[3 * n_verts + 1] = verts_[3 * i + 1] + offset4[1];
		verts[3 * n_verts + 2] = verts_[3 * i + 2] + offset4[2];
		n_verts++;
	}

	n_faces = 0;
	faces[3 * n_faces] = 0;
	faces[3 * n_faces+1] = 3;
	faces[3 * n_faces + 2] = 2;
	n_faces++;

	faces[3 * n_faces] = 0;
	faces[3 * n_faces + 1] = 2;
	faces[3 * n_faces + 2] = 1;
	n_faces++;

	faces[3 * n_faces] = n_verts-4;
	faces[3 * n_faces + 1] = n_verts - 3;
	faces[3 * n_faces + 2] = n_verts - 2;
	n_faces++;

	faces[3 * n_faces] = n_verts - 4;
	faces[3 * n_faces + 1] = n_verts - 2;
	faces[3 * n_faces + 2] = n_verts - 1;
	n_faces++;

	for (int i = 0; i < np-1; i++)
	{
		faces[3 * n_faces] = 4*i;
		faces[3 * n_faces + 1] = 4 * i+4;
		faces[3 * n_faces + 2] = 4 * i+7;
		n_faces++;

		faces[3 * n_faces] = 4 * i;
		faces[3 * n_faces + 1] = 4 * i + 7;
		faces[3 * n_faces + 2] = 4 * i + 3;
		n_faces++;

		faces[3 * n_faces] = 4 * i+1;
		faces[3 * n_faces + 1] = 4 * i + 5;
		faces[3 * n_faces + 2] = 4 * i + 4;
		n_faces++;

		faces[3 * n_faces] = 4 * i + 1;
		faces[3 * n_faces + 1] = 4 * i + 4;
		faces[3 * n_faces + 2] = 4 * i;
		n_faces++;

		faces[3 * n_faces] = 4 * i + 2;
		faces[3 * n_faces + 1] = 4 * i + 6;
		faces[3 * n_faces + 2] = 4 * i + 5;
		n_faces++;

		faces[3 * n_faces] = 4 * i + 2;
		faces[3 * n_faces + 1] = 4 * i + 5;
		faces[3 * n_faces + 2] = 4 * i+1;
		n_faces++;

		faces[3 * n_faces] = 4 * i + 3;
		faces[3 * n_faces + 1] = 4 * i + 7;
		faces[3 * n_faces + 2] = 4 * i + 6;
		n_faces++;

		faces[3 * n_faces] = 4 * i + 3;
		faces[3 * n_faces + 1] = 4 * i + 6;
		faces[3 * n_faces + 2] = 4 * i + 2;
		n_faces++;
	}

	hasSegments = false;
	hasFaceNormals = false;
	hasVertNormals = false;
	hasVertColors = false;
}

void TriMesh::createCube(double lx, double ly, double lz)
{
	double cube_verts[8 * 3] = { -1, -1, -1, -1, -1, 1, 1, -1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, 1, 1, 1, 1, 1, -1 };
	int n_cube_verts = 8;
	int cube_faces[12 * 3] = { 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6, 3, 0, 4, 3, 4, 7, 5, 6, 7, 5, 7, 4, 1, 0, 3, 1, 3, 2 };
	int n_cube_faces = 12;

	n_verts = n_cube_verts;
	n_faces = n_cube_faces;
	verts = (double*)realloc(verts, sizeof(double) * 3 * n_verts);
	faces = (int*)realloc(faces, sizeof(int) * 3 * n_faces);

	memcpy(verts, cube_verts, sizeof(cube_verts[0])*n_verts * 3);
	memcpy(faces, cube_faces, sizeof(cube_faces[0])*n_faces * 3);

	scale(lx / 2, ly / 2, lz / 2);

	hasSegments = false;
	hasFaceNormals = false;
	hasVertNormals = false;
	hasVertColors = false;
}

bool TriMesh::loadSeg(const std::string &filename, bool do_bound)
{
	std::ifstream fin;
	fin.open(filename);

	if (fin.fail())
		return false;

	fin >> n_verts;
	verts = (double*)realloc(verts, sizeof(double) * 3 * n_verts);
	for (int i = 0; i < n_verts; i++)
	{
		fin >> verts[3 * i] >> verts[3 * i + 1] >> verts[3 * i + 2];
	}

	int n_patches;
	fin >> n_faces>>n_patches;

	faces = (int*)realloc(faces, sizeof(int) * 3 * n_faces);
	face_ids = (int*)realloc(face_ids, sizeof(int) * n_faces);
	
	if (do_bound)
	{
		patches.resize(n_patches, Patch());
		for (int i = 0; i < n_faces; i++)
		{
			fin >> faces[3 * i] >> faces[3 * i + 1] >> faces[3 * i + 2] >> face_ids[i];
			patches[face_ids[i]].faces.push_back(i);
		}

		int edge_len, p1, p2, v1, v2;
		while (fin >> edge_len)
		{
			Edge edge;
			fin >> p1 >> p2;
			edge.sides.first = p1;
			edge.sides.second = p2;
			for (int i = 0; i < edge_len; i++)
			{
				fin >> v1 >> v2;
				if (edge.verts.empty())
					edge.verts.push_back(v1);
				edge.verts.push_back(v2);
			}
			patches[p1].edges.push_back(edge);

			Edge edge2;
			edge2.sides.first = p2;
			edge2.sides.second = p1;
			for (int i = edge.verts.size() - 1; i >= 0; i--)
			{
				edge2.verts.push_back(edge.verts[i]);
			}
			patches[p2].edges.push_back(edge2);
		}
	}
	else
	{
		for (int i = 0; i < n_faces; i++)
		{
			fin >> faces[3 * i] >> faces[3 * i + 1] >> faces[3 * i + 2] >> face_ids[i];
		}
		buildPatches();
	}

	fin.close();
	hasSegments = true;

	return true;
}

void TriMesh::loadMesh(double *verts_, int n_verts_, int *faces_, int n_faces_)
{
	n_verts = n_verts_;
	n_faces = n_faces_;

	verts = (double*)realloc(verts, sizeof(double) * 3 * n_verts);
	faces = (int*)realloc(faces, sizeof(int) * 3 * n_faces);

	memcpy(verts, verts_, sizeof(double) * 3 * n_verts);
	memcpy(faces, faces_, sizeof(int) * 3 * n_faces);

	hasSegments = false;
	hasFaceNormals = false;
	hasVertNormals = false;
	hasVertColors = false;
}

void TriMesh::loadMesh(std::vector<double> &vert_pos, std::vector<int> &face_index)
{
	n_verts = vert_pos.size()/3;
	n_faces = face_index.size()/3;

	verts = (double*)realloc(verts, sizeof(double) * 3 * n_verts);
	faces = (int*)realloc(faces, sizeof(int) * 3 * n_faces);


	for (int v = 0; v < vert_pos.size(); v++)
	{
		verts[v] = vert_pos[v];
	}
	for (int v = 0; v < face_index.size(); v++)
	{
		faces[v] = face_index[v];
	}

	hasSegments = false;
	hasFaceNormals = false;
	hasVertNormals = false;
	hasVertColors = false;
}


#ifndef PAIR
#define PAIR std::pair<int, int>
#endif

void TriMesh::writeV(const std::string &filename, float *vs)
{
	std::ofstream ofs(filename);
	for (int i = 0; i < n_verts; i++)
	{
		ofs << vs[i] << "\n";
	}
	ofs.close();
}

void TriMesh::setV(float *vs)
{
	for (int i = 0; i < n_verts; i++)
	{
		vs[i] = 0.5 - vert_colors[3 * i];
		if (vs[i] < 0)vs[i] = 0;
	}
}

void TriMesh::loadV(const std::string &filename, float *vs)
{
	std::ifstream ifs(filename);
	for (int i = 0; i < n_verts; i++)
	{
		ifs >> vs[i];
	}
	ifs.close();

}
void  TriMesh::smoothV(float *vs, int iters)
{
	//build topology
	std::set<int> *adjVerts = new std::set<int>[n_verts];	//edge->face
	for (int i = 0; i < n_faces; i++)
	{
		adjVerts[faces[3 * i]].insert(faces[3 * i + 1]);
		adjVerts[faces[3 * i]].insert(faces[3 * i + 2]);
		adjVerts[faces[3 * i + 1]].insert(faces[3 * i]);
		adjVerts[faces[3 * i + 1]].insert(faces[3 * i + 2]);
		adjVerts[faces[3 * i + 2]].insert(faces[3 * i + 1]);
		adjVerts[faces[3 * i + 2]].insert(faces[3 * i]);
	}

	float *vs_tmp = new float[n_verts];

	for (int k = 0; k < iters;k++)
	{
		for (int i = 0; i < n_verts; i++)
		{
			vs_tmp[i] = vs[i];
			for (std::set<int>::iterator iter = adjVerts[i].begin(); iter != adjVerts[i].end(); iter++)
			{
				int v_adj = *iter;
				vs_tmp[i] += vs[v_adj];
			}
			vs[i] = vs_tmp[i] / (adjVerts[i].size() + 1);
		}
	}
	
	
	delete[]vs_tmp;
	delete[]adjVerts;

}


void TriMesh::setHeatColor(float *vs, float sc, bool edgeOnly)
{
	float v_max=0;
	for (int i = 0; i < n_verts; i++)
	{
		v_max = std::max(v_max, vs[i]);
	}

	vert_colors = (float*)realloc(vert_colors, sizeof(float) * 3 * n_verts);
	
	if (edgeOnly)
	{
		for (int i = 0; i < n_verts; i++)
		{
			VectorSetConst3f(vert_colors + 3 * i, 0.7f);
		}

		int n_patches = patches.size();
		for (int i = 0; i < n_patches; i++)
		{
			for (int j = 0; j < patches[i].edges.size(); j++)
			{
				Edge &edge = patches[i].edges[j];
				if (edge.otherSide()>i)
				{
					for (int k = 0; k < edge.len() + 1; k++)
					{
						int v_idx = edge.verts[k];
						float v_ = vs[v_idx] / v_max *sc;
						getHeatColor(v_, vert_colors + 3 * v_idx, vert_colors + 3 * v_idx + 1, vert_colors + 3 * v_idx + 2);
					}
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < n_verts; i++)
		{
			float v_ = vs[i] / v_max *sc;
			getHeatColor(v_, vert_colors + 3 * i, vert_colors + 3 * i + 1, vert_colors + 3 * i + 2);
		}

	}

	hasVertColors = true;
}

void TriMesh::writeOff(const std::string &filename)
{
	std::ofstream ofs(filename);

	if (hasVertColors)
		ofs << "COFF\n";
	else
		ofs << "OFF\n";

	ofs << n_verts << " " << n_faces << " 0\n";

	// write all of the points.
	for (int i = 0; i < n_verts; i++)
	{
		ofs << verts[3 * i] << " " << verts[3 * i + 1] << " " << verts[3 * i + 2];
		if (hasVertColors)
		{
			ofs << " " << vert_colors[3 * i] * 255.0f << " " << vert_colors[3 * i + 1] * 255.0f << " " << vert_colors[3 * i + 1] * 255.0f;
		}
		ofs << "\n";
	}

	// write all of the faces.
	for (int i = 0; i < n_faces; i++)
	{
		ofs << faces[3 * i] << " " << faces[3 * i + 1] << " " << faces[3 * i + 2]<<"\n";
	}

	ofs.close();
}

bool TriMesh::loadOff(const std::string &filename)
{
	std::ifstream ifs(filename);
	//std::string tempBuf;
	char tempBuf[1024];
	bool goodLoad = true;
	int numEdges;
	int i;
	bool loadColor;
	// Grab the first string. If it's "OFF", we think this
	// is an OFF file and continue. Otherwise we give up.
	ifs >> tempBuf;
	//if (tempBuf.compare("OFF")!=0) 
	if (strcmp(tempBuf, "OFF") == 0)
	{
		loadColor = false;
	}
	else if (strcmp(tempBuf, "COFF") == 0)
	{
		loadColor = true;
	}
	else
		goodLoad = false;

	if (goodLoad) {
		ifs >> n_verts >> n_faces >> numEdges;
		if (n_verts < 1 ||
			n_faces < 1) {
			goodLoad = false;
			n_verts = 0;
			n_faces = 0;
		}
		else {
			verts = (double*)realloc(verts, sizeof(double) * 3 * n_verts);
			if (loadColor)
				vert_colors = (float*)realloc(vert_colors, sizeof(float) * 3 * n_verts);
			faces = (int*)realloc(faces, sizeof(int) * 3 * n_faces);
		}
	}

	float alpha;
	if (goodLoad) 
	{
		// Load all of the points.
		for (i = 0; i < n_verts; i++) 
		{
			ifs >> verts[3 * i] >> verts[3 * i + 1] >> verts[3 * i + 2];
			if (loadColor)
			{
				ifs >> vert_colors[3 * i] >> vert_colors[3 * i + 1] >> vert_colors[3 * i + 2] >> alpha;
				vert_colors[3 * i] /= 255.0;
				vert_colors[3 * i+1] /= 255.0;
				vert_colors[3 * i+2] /= 255.0;
			}
		}

		int facesize;
		// Load all of the faces.
		for (i = 0; i < n_faces; i++) 
		{
			// This tells us how many points make up
			// this face.
			ifs >> facesize;
			if (facesize != 3)
			{
				goodLoad = false;
				std::cout << "ERROR: " << filename << " is not a valid triangle mesh\n";
				break;
			}

			// So we declare a new array of that size
			// And load its elements with the vertex indices.
			ifs >> faces[3 * i] >> faces[3 * i+1] >> faces[3 * i+2];
		}
	}

	hasSegments = false;
	hasFaceNormals = false;
	hasVertNormals = false;
	hasVertColors = loadColor;
	ifs.close();
	return goodLoad;

}

static bool naive_load_obj(const std::string & filename, std::vector<double> &vertices, std::vector<int> &elements)
{
	std::ifstream in(filename);
	if (!in)
	{
		std::cout << "Cannot open " << filename << "\n";
		return false;
	}

	vertices.clear();
	elements.clear();

	std::string line;
	while (getline(in, line))
	{
		if (line.substr(0, 2) == "v ")
		{
			std::istringstream s(line.substr(2));
			double v[3]; 
			s >> v[0] >> v[1] >> v[2];
			vertices.push_back(v[0]);
			vertices.push_back(v[1]);
			vertices.push_back(v[2]);
		}
		else if (line.substr(0, 2) == "f ")
		{
			std::istringstream s(line.substr(2));
			std::string buf;
			int v[3];

			int cnt = 0;
			while (s >> buf)
			{
				std::istringstream ss(buf);
				ss >> v[cnt];
				cnt++;
			}

			if (cnt != 3)
			{
				std::cout << "ERROR: " << filename << " is not a valid triangle mesh\n";
				exit(-1);
			}

			elements.push_back(v[0] - 1);
			elements.push_back(v[1] - 1);
			elements.push_back(v[2] - 1);
		}
		else if (line[0] == '#')
		{
			/* ignoring this line */
		}
		else
		{
			/* ignoring this line */
		}
	}
	
	in.close();
	return true;
}

bool TriMesh::loadObj(const std::string &filename)
{
	std::vector<double> vertices;
	std::vector<int> elements;
	if (!naive_load_obj(filename, vertices, elements))
		return false;

	n_verts = vertices.size() / 3;
	n_faces = elements.size() / 3;
	verts = (double*)realloc(verts, sizeof(double) * n_verts*3);
	faces = (int*)realloc(faces, sizeof(int) * n_faces*3);

	for (int i = 0; i < n_verts * 3; i++)
	{
		verts[i] = vertices[i];
	}
	for (int i = 0; i < n_faces * 3; i++)
	{
		faces[i] = elements[i];
	}

	hasSegments = false;
	hasFaceNormals = false;
	hasVertNormals = false;
	hasVertColors = false;
	return true;
}

void TriMesh::writeObj(const std::string &filename, double *verts, int num_verts, int *tris, int num_tris, double *offset, double scale)
{
	std::ofstream fout(filename);

	for (int i = 0; i < num_verts; i++)
	{
		fout << "v " << offset[0] + scale*verts[3 * i] << " " << offset[0] + scale*verts[3 * i + 1] << " " << offset[0] + scale*verts[3 * i + 2] << "\n";
	}

	for (int i = 0; i < num_tris; i++)
	{
		fout << "f " << tris[3 * i] + 1 << " " << tris[3 * i + 1] + 1 << " " << tris[3 * i + 2] + 1 << "\n";
	}

	fout.close();
}

void TriMesh::writeObj(const std::string &filename, double *verts, int num_verts, int *tris, int num_tris)
{
	std::ofstream fout(filename);

	for (int i = 0; i < num_verts; i++)
	{
		fout << "v " << verts[3 * i] << " " << verts[3 * i + 1] << " " << verts[3 * i + 2] << "\n";
	}

	for (int i = 0; i < num_tris; i++)
	{
		fout << "f " << tris[3 * i] + 1 << " " << tris[3 * i + 1] + 1 << " " << tris[3 * i + 2] + 1 << "\n";
	}

	fout.close();
}

void TriMesh::writeMesh(TriMesh &mesh)
{
	std::vector<double> vert_pos;
	std::vector<int> face_index;
	writeMesh(vert_pos, face_index);
	mesh.loadMesh(vert_pos, face_index);
}

void TriMesh::writeMesh(std::vector<double> &vert_pos, std::vector<int> &face_index)
{
	vert_pos.clear();
	face_index.clear();

	for (int i = 0; i < n_verts; i++)
	{
		vert_pos.push_back(verts[3 * i]);
		vert_pos.push_back(verts[3 * i+1]);
		vert_pos.push_back(verts[3 * i+2]);
	}

	for (int i = 0; i < n_faces; i++)
	{
		face_index.push_back(faces[3 * i]);
		face_index.push_back(faces[3 * i+1]);
		face_index.push_back(faces[3 * i+2]);
	}
}
void TriMesh::writeObj(const std::string &filename)
{
	std::ofstream fout(filename);

	for (int i = 0; i < n_verts; i++)
	{
		fout << "v " << verts[3 * i] << " " << verts[3 * i + 1] << " " << verts[3 * i + 2] << "\n";
	}

	for (int i = 0; i < n_faces; i++)
	{
		fout << "f " << faces[3 * i] + 1 << " " << faces[3 * i + 1] + 1 << " " << faces[3 * i + 2] + 1 << "\n";
	}

	fout.close();

}

// -1: scale 3 dimensional; 0-2: x,y,z
void TriMesh::scale(double scx, double scy, double scz)
{
		for (int i = 0; i < n_verts; i++)
		{
			verts[3 * i] *= scx;
			verts[3 * i+1] *= scy;
			verts[3 * i+2] *= scz;
		}
}

void TriMesh::rotate(Quaternion &q)
{
	double v[3];
	for (int i = 0; i < n_verts; i++)
	{
		q.rotateVec(verts + 3 * i, v);
		VectorSet3d(verts + 3 * i, v);
	}
	hasFaceNormals = false;
	hasVertNormals = false;
}

// 0-2: x,y,z, in radians
void TriMesh::rotate(int axis, double rot)
{
	Quaternion q;
	switch (axis)
	{
	case 0:
		q.init(0, 0, rot);
		break;
	case 1:
		q.init(0, rot, 0);
		break;
	case 2:
		q.init(rot, 0, 0);
		break;		
	}

	//double mat[9];
	//q.toMatrix3(mat);
	//double v[3];
	//for (int i = 0; i < n_verts; i++)
	//{
	//	MATRIX_VECTOR_MULTIPLY3X3(mat, verts + 3 * i, v);
	//	VectorSet3d(verts + 3 * i, v);
	//}
	double v[3];
	for (int i = 0; i < n_verts; i++)
	{
		q.rotateVec(verts + 3 * i, v);
		VectorSet3d(verts + 3 * i, v);
	}

	hasFaceNormals = false;
	hasVertNormals = false;
}

void TriMesh::translate(double x, double y, double z)
{
	for (int i = 0; i < n_verts; i++)
	{
		verts[3 * i] += x;
		verts[3 * i+1] += y;
		verts[3 * i+2] += z;
	}
}

void TriMesh::buildBbox(double *bmin, double *bmax)
{
	VectorSetConst3d(bmin, 10000);
	VectorSetConst3d(bmax, -10000);

	for (int i = 0; i < n_verts; i++)
	{
		VectorMax3d(bmax, verts + 3 * i, bmax);
		VectorMin3d(bmin, verts + 3 * i, bmin);
	}
}


typedef std::pair<int, int> Pair;


void TriMesh::buildPatches()
{
	// count # of segments and construct
	int n_patches = 0;

	for (int i = 0; i < n_faces; i++)
	{
		if (face_ids[i] + 1>n_patches)
			n_patches = face_ids[i] + 1;
	}

	patches.resize(n_patches, Patch());
	for (int i = 0; i < n_faces; i++)
	{
		patches[face_ids[i]].faces.push_back(i);
	}


	//	build topology
	std::map<Pair, int> adjFaces;	//edge->face
	std::set<int> *adjVerts = new std::set<int>[n_verts];	//edge->face
	for (int i = 0; i < n_faces; i++)
	{
		adjFaces[Pair(faces[3 * i], faces[3 * i + 1])] = i;
		adjFaces[Pair(faces[3 * i + 1], faces[3 * i + 2])] = i;
		adjFaces[Pair(faces[3 * i + 2], faces[3 * i])] = i;
		adjVerts[faces[3 * i]].insert(faces[3 * i + 1]);
		adjVerts[faces[3 * i]].insert(faces[3 * i + 2]);
		adjVerts[faces[3 * i + 1]].insert(faces[3 * i]);
		adjVerts[faces[3 * i + 1]].insert(faces[3 * i + 2]);
		adjVerts[faces[3 * i + 2]].insert(faces[3 * i + 1]);
		adjVerts[faces[3 * i + 2]].insert(faces[3 * i]);
	}

	/////////////Extract boundary for each patch
	std::set<Pair> visited;
	for (std::map<Pair, int>::iterator iter = adjFaces.begin(); iter != adjFaces.end(); iter++)
	{
		if (visited.find(iter->first) != visited.end())
			continue;

		int v1 = iter->first.first;
		int v2 = iter->first.second;	
		int f = face_ids[iter->second];
		int f_ = face_ids[adjFaces[Pair(v2, v1)]];

		if (f != f_)
		{
			std::list<int> bd;
			bd.push_back(v1);
			bd.push_back(v2);

			while (v1 != v2)
			{
				int next = -1;
				for (std::set<int>::iterator its = adjVerts[v1].begin(); its != adjVerts[v1].end(); its++)	// search vertices before v1
				{
					if (face_ids[adjFaces[Pair(*its, v1)]] == f && face_ids[adjFaces[Pair(v1, *its)]] != f)
					{
						next = *its;
						break;
					}
				}
				if (next < 0) //not a complete boundary
				{
					printf("ERROR: not a complete boundary!\n");
					exit(-1);
				}
				else
				{
					if (face_ids[adjFaces[Pair(v1, next)]] == f_)
					{
						bd.push_front(next);
						visited.insert(Pair(next, v1));
						v1 = next;
					}
					else
						break;
				}
			}


			if (v1 == v2)	// finished
			{
				// create a new edge
				Edge edge;
				edge.sides.first = f;
				edge.sides.second = f_;
				edge.verts = std::vector<int>(bd.begin(), bd.end());
				patches[f].edges.push_back(edge);
				bd.clear();
			}
			else
			{
				while (v1 != v2)
				{
					int next = -1;
					for (std::set<int>::iterator its = adjVerts[v2].begin(); its != adjVerts[v2].end(); its++)	// search vertices after v2
					{
						if (face_ids[adjFaces[Pair(v2, *its)]] == f && face_ids[adjFaces[Pair(*its, v2)]] != f)
						{
							next = *its;
							break;
						}
					}
					if (next < 0) //not a complete boundary
					{
						printf("ERROR: not a complete boundary!\n");
						exit(-1);
					}
					else
					{
						if (face_ids[adjFaces[Pair(next, v2)]] == f_)
						{
							bd.push_back(next);
							visited.insert(Pair(v2, next));
							v2 = next;
						}
						else
						{
							// create a new edge
							Edge edge;
							edge.sides.first = f;
							edge.sides.second = f_;
							edge.verts = std::vector<int>(bd.begin(), bd.end());
							patches[f].edges.push_back(edge);

							// go on search next segment with a different neighbor patch
							bd.clear();
							f_ = face_ids[adjFaces[Pair(next, v2)]];
							bd.push_back(v2);
							bd.push_back(next);
							visited.insert(Pair(v2, next));
							v2 = next;
						}
					}
				}

				if (!bd.empty())
				{
					// create a new edge
					Edge edge;
					edge.sides.first = f;
					edge.sides.second = f_;
					edge.verts = std::vector<int>(bd.begin(), bd.end());
					patches[f].edges.push_back(edge);
					bd.clear();
				}
			}
			
		}
	}


	delete[]adjVerts;
}

void TriMesh::loadMeshSeg(std::vector<double> &vert_pos, std::vector<int> &face_index, std::vector<int> &face_id)
{
	// copy verts and faces
	n_verts = vert_pos.size() / 3;
	n_faces = face_index.size() / 3;
	verts = (double*)realloc(verts, sizeof(double) * 3 * n_verts);
	faces = (int*)realloc(faces, sizeof(int) * 3 * n_faces);
	face_ids = (int*)realloc(face_ids, sizeof(int) * n_faces);
	for (int v = 0; v < vert_pos.size(); v++)
	{
		verts[v] = vert_pos[v];
	}
	for (int v = 0; v < face_index.size(); v++)
	{
		faces[v] = face_index[v];
		face_ids[v] = face_id[v];
	}
	buildPatches();

	hasSegments = true;
	hasFaceNormals = false;
	hasVertNormals = false;
	hasVertColors = false;
}


void TriMesh::segment()
{
	if (hasSegments)
		return;

	if (!hasVertColors)
	{
		printf("ERROR: segment with no colors assigned!\n");
		exit(-1);
	}

	// build topology
	std::map<Pair, int> adjFaces;	//edge->face
	for (int i = 0; i < n_faces; i++)
	{
		adjFaces[Pair(faces[3 * i], faces[3 * i + 1])] = i;
		adjFaces[Pair(faces[3 * i + 1], faces[3 * i + 2])] = i;
		adjFaces[Pair(faces[3 * i + 2], faces[3 * i])] = i;
	}


	// BFS to build face ids
	face_ids = (int*)realloc(face_ids, sizeof(int) * n_faces);
	for (int i = 0; i < n_faces; i++)
	{
		face_ids[i] = -1;
	}

	int n_patches = 0;
	for (int i = 0; i < n_faces; i++)
	{
		if (face_ids[i] == -1)
		{
			int seed = i;

			std::queue<int> q;
			q.push(seed);
			face_ids[seed] = n_patches;

			while (!q.empty())
			{
				int current = q.front();
				q.pop();

				for (int j = 0; j < 3; j++)
				{
					int v1 = faces[3 * current + j];
					int v2 = faces[3 * current + (j + 1) % 3];
					if (!isSeg(v1) || !isSeg(v2))
					{
						int f = adjFaces[Pair(v2, v1)];
						if (face_ids[f] == -1)
						{
							q.push(f);
							face_ids[f] = n_patches;
						}
					}
				}
			}
			n_patches++;
		}
	}

	buildPatches();

	printf("Segment finished with %d segments\n", patches.size());

	hasSegments = true;
}



void TriMesh::exportSeg(const std::string &filename, bool do_bound)
{
	if (!hasSegments)
		segment();
	std::ofstream fout(filename);

	int n_patches = patches.size();

	fout << n_verts <<"\n";
	for (int i = 0; i < n_verts; i++)
	{
		fout << verts[3 * i] << " " << verts[3 * i + 1] << " " << verts[3 * i + 2] << "\n";
	}
	fout << n_faces << " " << n_patches << "\n";
	for (int i = 0; i < n_faces; i++)
	{
		fout << faces[3 * i] << " " << faces[3 * i + 1] << " " << faces[3 * i + 2] << " " << face_ids[i] << "\n";
	}

	if (do_bound)	// export boundary of each patch
	{
		for (int i = 0; i < n_patches-1; i++)
		{
			for (int j = 0; j < patches[i].edges.size(); j++)
			{
				Edge &edge = patches[i].edges[j];
				if (edge.otherSide()>i)
				{
					fout << edge.len() << " " << i << " " << edge.otherSide() << "\n";
					for (int k = 0; k < edge.len(); k++)
					{
						fout << edge.verts[k] << " " << edge.verts[k + 1] << "\n";
					}
				}
			}
		}
	}
	fout.close();
}

void TriMesh::reorder(int *order)
{
	int n_patches = patches.size();
	int *rank = new int[n_patches];

	for (int i = 0; i < n_patches; i++)
	{
		rank[order[i]] = i;
	}

	Patch *tmp_patches = new Patch[n_patches];
	for (int i = 0; i < n_patches;i++)
	{
		tmp_patches[i] = patches[order[i]];
	}
	for (int i = 0; i < n_patches; i++)
	{
		patches[i] = tmp_patches[i];
		for (std::vector<int>::iterator iter = patches[i].faces.begin(); iter != patches[i].faces.end(); iter++)
		{
			face_ids[*iter] = i;
		}
		for (std::vector<Edge>::iterator iter = patches[i].edges.begin(); iter != patches[i].edges.end(); iter++)
		{
			iter->sides.first = i;
			iter->sides.second = rank[iter->sides.second];
		}
	}
	delete[]rank;
	delete[]tmp_patches;
}

void TriMesh::renderFlatAllPatchesTrans(float alpha)
{
	if (!hasFaceNormals)
		buildFaceNormals();
	if (!hasSegments)
		segment();

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	int p1, p2, p3;
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < n_faces; i++)
	{
		glNormal3fv(faceNormals + 3 * i);
		int id = face_ids[i];
		glColor4f(renderColor[3 * id], renderColor[3 * id + 1], renderColor[3 * id + 2], alpha);
		p1 = faces[3 * i]; p2 = faces[3 * i + 1]; p3 = faces[3 * i + 2];

		glVertex3dv(verts + 3 * p1);
		glVertex3dv(verts + 3 * p2);
		glVertex3dv(verts + 3 * p3);
	}
	glEnd();
	glDisable(GL_BLEND);
}
void TriMesh::renderFlatAllPatches()
{
	if (!hasFaceNormals)
		buildFaceNormals();
	if (!hasSegments)
		segment();

	int p1, p2, p3;

	glBegin(GL_TRIANGLES);
	for (int i = 0; i < n_faces; i++)
	{
		glNormal3fv(faceNormals + 3 * i);
		int id = face_ids[i];
		glColor3f(renderColor[3 * id], renderColor[3 * id + 1], renderColor[3 * id + 2]);
		p1 = faces[3 * i]; p2 = faces[3 * i + 1]; p3 = faces[3 * i + 2];

		glVertex3dv(verts + 3 * p1);
		glVertex3dv(verts + 3 * p2);
		glVertex3dv(verts + 3 * p3);
	}
	glEnd();
}


void TriMesh::renderFlatPatch(int id)
{

	if (!hasFaceNormals)
		buildFaceNormals();
	if (!hasSegments)
		segment();

	glDisable(GL_CULL_FACE);
	int p1, p2, p3;


	glBegin(GL_TRIANGLES);
	glColor3f(renderColor[3 * id], renderColor[3 * id + 1], renderColor[3 * id + 2]);
	for (std::vector<int>::iterator iter = patches[id].faces.begin(); iter != patches[id].faces.end(); iter++)
	{
		glNormal3fv(faceNormals + 3 * (*iter));

		p1 = faces[3 * (*iter)]; p2 = faces[3 * (*iter) + 1]; p3 = faces[3 * (*iter) + 2];

		glVertex3dv(verts + 3 * p1);
		glVertex3dv(verts + 3 * p2);
		glVertex3dv(verts + 3 * p3);
	}
	glEnd();

	glEnable(GL_CULL_FACE);
}

void TriMesh::renderEdge()
{
	if (!hasSegments)
		segment();

	if (!hasVertColors)
		glColor3f(0.0, 0.0, 0.0);
	
	glLineWidth(2.0);

	glBegin(GL_LINES);
	int n_patches = patches.size();
	for (int i = 0; i <n_patches; i++)
	{
		for (int j = 0; j < patches[i].edges.size(); j++)
		{
			Edge &edge = patches[i].edges[j];
			
			if (edge.otherSide()>edge.thisSide())
			{
				for (int k = 0; k < edge.len(); k++)
				{
					if (hasVertColors)
					{
						glColor3fv(vert_colors + 3 * edge.verts[k]);
					}
					glVertex3dv(verts + 3 * edge.verts[k]);
					if (hasVertColors)
					{
						glColor3fv(vert_colors + 3 * edge.verts[k+1]);
					}
					glVertex3dv(verts + 3 * edge.verts[k + 1]);
				}
			}
		}
	}

	glEnd();

}

void TriMesh::renderFlatTrans(float *override_color, float alpha)
{
	if (!hasFaceNormals)
		buildFaceNormals();

	int p1, p2, p3;
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	if (!override_color && !hasVertColors)
	{
		glColor4f(0.5, 0.5, 0.5, alpha);
	}
	else if (override_color)
	{
		glColor4f(override_color[0], override_color[1], override_color[2], alpha);
	}

	glBegin(GL_TRIANGLES);
	for (int i = 0; i < n_faces; i++)
	{
		glNormal3fv(faceNormals + 3 * i);

		p1 = faces[3 * i]; p2 = faces[3 * i + 1]; p3 = faces[3 * i + 2];

		if (!override_color)
		{
			if (hasVertColors)
			{
				glColor4f(vert_colors[3 * p1], vert_colors[3 * p1 + 1], vert_colors[3 * p1+2], alpha);
				glVertex3dv(verts + 3 * p1);
				glColor4f(vert_colors[3 * p2], vert_colors[3 * p2 + 1], vert_colors[3 * p2 + 2], alpha);
				glVertex3dv(verts + 3 * p2);
				glColor4f(vert_colors[3 * p3], vert_colors[3 * p3 + 1], vert_colors[3 * p3 + 2], alpha);
				glVertex3dv(verts + 3 * p3);
			}
		}
		else
		{
			glVertex3dv(verts + 3 * p1);
			glVertex3dv(verts + 3 * p2);
			glVertex3dv(verts + 3 * p3);
		}
	}
	glEnd();

	glDisable(GL_BLEND);
}

void TriMesh::renderFlat(float *override_color)
{
	if (!hasFaceNormals)
		buildFaceNormals();

	int p1, p2, p3;

	if (!override_color && !hasVertColors)
	{
		glColor3f(0.5,0.5,0.5);
	}
	else if (override_color)
	{
		glColor3fv(override_color);
	}

	glBegin(GL_TRIANGLES);
	for (int i = 0; i < n_faces; i++)
	{
		glNormal3fv(faceNormals + 3 * i);

		p1 = faces[3 * i]; p2 = faces[3 * i + 1]; p3 = faces[3 * i + 2];

		if (!override_color)
		{
			if (hasVertColors)
			{
				glColor3fv(vert_colors + 3 * p1);
				glVertex3dv(verts + 3 * p1);
				glColor3fv(vert_colors + 3 * p2);
				glVertex3dv(verts + 3 * p2);
				glColor3fv(vert_colors + 3 * p3);
				glVertex3dv(verts + 3 * p3);
			}
		}
		else
		{
			glVertex3dv(verts + 3 * p1);
			glVertex3dv(verts + 3 * p2);
			glVertex3dv(verts + 3 * p3);
		}
	}
	glEnd();
}

void TriMesh::renderSmooth(float *override_color)
{
	if (!hasVertNormals)
		buildVertNormals();

	int p1, p2, p3;

	if (!override_color && !hasVertColors)
	{
		glColor3f(0.5, 0.5, 0.5);
	}
	else if (override_color)
	{
		glColor3fv(override_color);
	}

	glBegin(GL_TRIANGLES);
	for (int i = 0; i < n_faces; i++)
	{
		p1 = faces[3 * i]; p2 = faces[3 * i + 1]; p3 = faces[3 * i + 2];

		if (!override_color)
		{
			if (hasVertColors)
			{
				glColor3fv(vert_colors + 3 * p1);
				glNormal3fv(vertNormals + 3 * p1);
				glVertex3dv(verts + 3 * p1);
				glColor3fv(vert_colors + 3 * p2);
				glNormal3fv(vertNormals + 3 * p2);
				glVertex3dv(verts + 3 * p2);
				glColor3fv(vert_colors + 3 * p3);
				glNormal3fv(vertNormals + 3 * p3);
				glVertex3dv(verts + 3 * p3);
			}
		}
		else
		{
			glNormal3fv(vertNormals + 3 * p1);
			glVertex3dv(verts + 3 * p1);
			glNormal3fv(vertNormals + 3 * p2);
			glVertex3dv(verts + 3 * p2);
			glNormal3fv(vertNormals + 3 * p3);
			glVertex3dv(verts + 3 * p3);
		}
	}
	glEnd();
}


void TriMesh::buildFaceNormals()
{
	if (!hasFaceNormals)
	{
		faceNormals = (float*)realloc(faceNormals, sizeof(float) * n_faces * 3);

		int n1, n2, n3;
		double norm[3];
		for (int i = 0; i < n_faces; i++)
		{
			n1 = faces[3 * i];
			n2 = faces[3 * i + 1];
			n3 = faces[3 * i + 2];

			TriangleNormal3d(verts + 3 * n1, verts + 3 * n2, verts + 3 * n3, norm);
			VECTOR_SET3(faceNormals + 3 * i, norm);
		}

		hasFaceNormals = true;
	}
}


void TriMesh::buildVertNormals()
{
	if (!hasFaceNormals)
		buildFaceNormals();

	if (!hasVertNormals)
	{
		vertNormals = (float*)realloc(vertNormals, sizeof(float) * n_verts * 3);

		std::set<int> *adjFaces = new std::set<int>[n_verts];	//edge->face
		for (int i = 0; i < n_faces; i++)
		{
			adjFaces[faces[3 * i]].insert(i);
			adjFaces[faces[3 * i + 1]].insert(i);
			adjFaces[faces[3 * i + 2]].insert(i);
		}
		for (int i = 0; i < n_verts; i++)
		{
			VectorSetConst3f(vertNormals + 3 * i, 0);
			for (std::set<int>::iterator iter = adjFaces[i].begin(); iter != adjFaces[i].end(); iter++)
			{
				VectorAdd3f(vertNormals + 3 * i, faceNormals + 3 * (*iter), vertNormals + 3 * i);
			}
			VectorNormalize3f(vertNormals + 3 * i, vertNormals + 3 * i);
		}

		hasVertNormals = true;
		delete[]adjFaces;
	}
}


void TriMesh::setDefaultColor()
{
	vert_colors = (float*)realloc(vert_colors, sizeof(float) * 3 * n_verts);
	for (int i = 0; i < 3 * n_verts; i++)
	{
		vert_colors[i] = 0.5;
	}
	hasVertColors = true;
}

void TriMesh::cleanUnrefVerts()
{
	int *newid = new int[n_verts];
	for (int i = 0; i < n_verts;i++)
	{
		newid[i] = -1;
	}

	
	for (int i = 0; i < 3*n_faces; i++)
	{
		int v = faces[i];
		newid[v] = 0;
	}

	int cnt = 0;
	for (int i = 0; i < n_verts; i++)
	{
		if (newid[i] != -1)
		{
			newid[i] = cnt;
			VectorSet3d(verts + 3 * cnt, verts + 3 * i);
			cnt++;
		}
	}

	for (int i = 0; i < 3 * n_faces; i++)
	{
		int v = faces[i];
		faces[i] = newid[v];
	}

	n_verts = cnt;
	verts = (double*)realloc(verts, sizeof(double) * 3 * cnt);
	delete[]newid;

	hasSegments = false;
	//hasFaceNormals = false;
	hasVertNormals = false;
	hasVertColors = false;
}