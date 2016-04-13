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
//  Provide a feature to segment a mesh by vetex color (black as boundary)
//**************************************************************************************

#ifndef __PRINTX_TRI_MESH_H__
#define __PRINTX_TRI_MESH_H__


#include <vector>
#include <string>
#include <set>
#include "quaternion.h"

namespace PRINTX_LIB
{
class TriMesh
{
public:
	double *verts;
	int n_verts;
	int *faces;
	int n_faces;

	// constructor
	TriMesh(const TriMesh &mesh);
	TriMesh();
	~TriMesh();
	void clone(const TriMesh &mesh);

	//	create a rectangle pyramid
	void createRecPyramid(double *p, double *p0, double *p1, double *p2, double *p3);

	//	create a hexagonal prism given length on each dimension
	void createHexPrism(double lx, double ly, double lz);

	//	create a cube given length on each dimension
	void createCube(double lx, double ly, double lz);

	//	create a OBB cube given center and XYZ axis
	void createOBBCube(double *pt_, double *offx, double *offy, double *offz);	// careful on the order
	
	//void createSphere(double x, double y, double z,  double rx, double ry, double rz);

	//	expand a curve to a rectangular-stick-like mesh, with various thickness on the curve points
	void createFromCurvVary(std::vector<double> &verts_, std::vector<double> & offset1, std::vector<double> &offset2, std::vector<double> &offset3, std::vector<double> &offset4);
	
	//	expand a curve to a rectangular-stick-like mesh, with uniform thickness on the curve points
	void createFromCurv(std::vector<double> &verts, double *offset1, double *offset2, double *offset3, double *offset4);	// careful on the order
	
	//	expand a curve to a polygonal-stick-like mesh, with uniform thickness on the curve points
	void createFromCurv(std::vector<double> &verts, double *offset, int n_off);	// careful on the order
	
	//	load/write mesh
	void writeMesh(TriMesh &mesh);
	void writeMesh(std::vector<double> &vert_pos, std::vector<int> &face_index);
	void loadMesh(std::vector<double> &vert_pos, std::vector<int> &face_index);
	void loadMesh(double *verts_, int n_verts_, int *faces_, int n_faces_);
	//	load a mesh with face ids and segment
	void loadMeshSeg(std::vector<double> &vert_pos, std::vector<int> &face_index, std::vector<int> &face_id);
	
	//	load/write Wavefront .obj, only verts and faces
	bool loadObj(const std::string &filename);
	void writeObj(const std::string &filename);
	static void writeObj(const std::string &filename, double *verts, int num_verts, int *tris, int num_tris);
	static void writeObj(const std::string &filename, double *verts, int num_verts, int *tris, int num_tris, double *offset, double scale);

	//	load/write .off, load/write vert colors if COFF
	bool loadOff(const std::string &filename);
	void writeOff(const std::string &filename);

	// build heatmap given input intensity values, sc for scale
	void setHeatColor(float *vs, float sc, bool edgeOnly);
	void loadV(const std::string &filename, float *vs);
	void writeV(const std::string &filename, float *vs);
	void setV(float *vs);	// build intensity values by vertex colors

	//	smooth value 
	void smoothV(float *vs, int iters);

	//	clean unrefernced vertices
	void cleanUnrefVerts();

	// transformation
	void scale(double scx, double scy, double scz);
	void rotate(int axis, double rot);// 0-2: x,y,z, in radians
	void rotate(Quaternion &q);
	void translate(double x, double y, double z);

	// build bounding box
	void buildBbox(double *bmin, double *bmax);

	
	/////////////////////////////for mesh segmentation
	// an edge consists of multiple vertices and 2 sides(ids)
	class Edge
	{
	public:
		std::vector<int> verts;
		std::pair<int, int> sides;
		Edge()
		{
			verts.clear();
			sides.first = -1;
			sides.second = -1;
		}
		int len()	// TODO: real physical length
		{
			return (int)verts.size() - 1;
		}
		int thisSide()
		{
			return sides.first;
		}
		int otherSide()
		{
			return sides.second;
		}
	};
	// a patch consists of multiple faces and multiple edges
	class Patch
	{
	public:
		std::vector<int> faces;
		std::vector<Edge> edges;
		Patch()
		{
			faces.clear();
			edges.clear();
		}
		int len()
		{
			int total = 0;
			for (std::vector<Edge>::iterator iter = edges.begin(); iter != edges.end(); iter++)
			{
				total += iter->len();
			}
			return total;
		}
	};
	bool hasSegments;
	int *face_ids;
	std::vector<Patch> patches;
	bool loadSeg(const std::string &filename, bool do_bound);	// do_bound: read patch info
	void exportSeg(const std::string &filename, bool do_bound);	// do_bound: export patch info
	void reorder(int *order);
	void segment();			// segment by vertex colors
	void buildPatches();	// build patches by face ids
	bool isSeg(int v)	// black color as separator
	{
		return vert_colors[3 * v] < 0.001 && vert_colors[3 * v + 1] < 0.001 && vert_colors[3 * v + 2] < 0.001;
	}


	/////////////////////////////normals for render
	bool hasFaceNormals;
	bool hasVertNormals;
	float* faceNormals;
	float* vertNormals;
	void buildFaceNormals();
	void buildVertNormals();
	void renderFlatTrans(float *override_color, float alpha);
	void renderFlat(float *override_color);
	void renderFlatPatch(int id);
	void renderFlatAllPatches();
	void renderFlatAllPatchesTrans(float alpha);
	void renderSmooth(float *override_color);
	void renderEdge();

	/////////////////////////////colors
	bool hasVertColors;
	float *vert_colors;	//0-255
	void setDefaultColor();
};
 }
#endif