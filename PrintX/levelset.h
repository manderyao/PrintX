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
***********************************************************************
//	level set
//	all level set values are positive and IDs are used to indicate inside or outside (ID=0 inside, ID=-1 outside)
//	Fedkiw, Stanley Osher Ronald, and S. Osher. "Level set methods and dynamic implicit surfaces." Surfaces 44 (2002): 77.
*/


#ifndef __PRINTX_LEVEL_SET_H__
#define __PRINTX_LEVEL_SET_H__

#include "triMesh.h"
#include <algorithm> 
#include "makeLevelset3.h"
#include <list>
#include "priorityQueueFM.h"
#include "vec3.h"
#include "cube.h"

namespace PRINTX_LIB
{
	class Levelset
	{
	public:
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////special cells for acceleration{
		enum{ UNPROCESSED = 0, PROCESSING = 1, PROCESSED = 2 };
		enum{ ID_INSIDE = 0, ID_OUTSIDE = -1 };

		// fast marching interface cells
		int *interface_cells;
		int num_interface_cells;

		inline void allocSpecialCells()
		{
			num_interface_cells = 0;
			interface_cells = (int*)malloc(sizeof(int)*n);
		}
		inline void deallocSpecialCells()
		{
			free(interface_cells);
		}

		void buildInterfaceCells();
		//////////////////special cells for acceleration}
		////////////////////////////////////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////helper function{
		inline int getIndex(int i, int j, int k)
		{
			return k*nxy + j*gn[0] + i;
		}
		inline void getAdjCells(int index, int*qindex)				// neighbors: -x, +x, -y, +y, -z, +z; if out bound, =-1
		{
			int *p = pos + 3 * index;
			qindex[0] = (p[0] - 1) >= 0 ? (index - 1) : -1;
			qindex[1] = (p[0] + 1) < gn[0] ? (index + 1) : -1;
			qindex[2] = (p[1] - 1) >= 0 ? (index - gn[0]) : -1;
			qindex[3] = (p[1] + 1) < gn[1] ? (index + gn[0]) : -1;
			qindex[4] = (p[2] - 1) >= 0 ? (index - nxy) : -1;
			qindex[5] = (p[2] + 1) < gn[2] ? (index + nxy) : -1;
		}
		inline void getCubeCells(int x, int y, int z, int*qindex)	// the 8 corners of this cube
		{
			for (int i = 0; i < 8; i++)
			{
				qindex[i] = getIndex(x + CUBE_STRUCT::vert_to_idx[i][0], y + CUBE_STRUCT::vert_to_idx[i][1], z + CUBE_STRUCT::vert_to_idx[i][2]);
			}
		}

		//	merge connected component to remove fragments
		void mergeCntCmpnt();

		//////////////////helper function}
		////////////////////////////////////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////Cell parameters{
		int gn[3];				//number of grids each dimension
		int n, nxy;				//total number, total number on xy slice
		double origin[3];
		double cell_size;
		double DIST_LIMIT;		// for those unprocessed cells

		//	level sets
		int *ids;				// -1 for outside, 0 for inside
		double *dists;			//distance to nearest interface point, always positive
		double *dist_tmp;
		unsigned char *dist_tags;		//UNPROCESSED = 0, PROCESSING = 1, PROCESSED = 2
		double *curvs;					//curvature
		int *pos;						//3d indices

		inline double getDist(int index)
		{
			return (ids[index] == ID_INSIDE) ? ((dist_tags[index] == PROCESSED) ? -dists[index] : -DIST_LIMIT) : ((dist_tags[index] == PROCESSED) ? dists[index] : DIST_LIMIT);
		}
		inline void allocCells()
		{
			ids = (int*)malloc(sizeof(int)*n);
			dists = (double*)malloc(sizeof(double)*n);
			dist_tmp = (double*)malloc(sizeof(double)*n);
			dist_tags = (unsigned char *)malloc(sizeof(unsigned char)*n);
			curvs = (double*)malloc(sizeof(double)*n);
			pos = (int*)malloc(sizeof(int) * 3 * n);
		}
		inline void deallocCells()
		{
			free(ids);
			free(dists);
			free(dist_tmp);
			free(dist_tags);
			free(curvs);
			free(pos);
		}
		void setAllDists(double dist)
		{
			for (int i = 0; i < n; i++)
			{
				dists[i] = dist;
			}
		}
		void setAllDistTags(unsigned char dist_tag)
		{
			for (int i = 0; i < n; i++)
			{
				dist_tags[i] = dist_tag;
			}
		}
		Levelset(TriMesh &mesh, double cellSize, bool enforceConnected, SDFGen::Array3i *closest_tri_ = NULL);
		~Levelset();

		double getDist(double x, double y, double z);	// get level set
		double getDistW(double x, double y, double z);	// input is world space coordinate
		//////////////////Cell parameters}
		////////////////////////////////////////////////////////////////////////////////////////////////////////



		////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////Redistance{
		PriorityQueueFM<double> *fm_queue = 0;
		void redist(double band_width=2.8);		//1+sqrt(3.0), enough for gradient/curvature calculation 
		void redistProcess();					// normalize level sets on the interface
		void fastMarch(double band_width);		// fast march
		inline void allocRedist()
		{
			fm_queue = new PriorityQueueFM<double>(n);
		}
		inline void deallocRedist()
		{
			delete fm_queue;
		}
		double compute_extended_distance(int index);
		void smooth(double curv_thres = 0.1, double speed = 0.2);
		void buildCurvature(int index);
		inline void advect(int idx, double v)
		{
			dists[idx] = getDist(idx) + v;
			ids[idx] = ID_OUTSIDE;
			if (dists[idx] < 0)
			{
				dists[idx] = -dists[idx];
				ids[idx] = ID_INSIDE;
			}
		}
		void advectAll(double v)
		{
			for (int i = 0; i < n; i++)
				advect(i, v);
		}
		void advectAllW(double v)	// world space
		{
			double vv = v / cell_size;
			for (int i = 0; i < n; i++)
				advect(i, vv);
		}
		//////////////////Redistance}
		////////////////////////////////////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////marching cube variables{
		int *mc_interface_cells;
		int num_mc_interface_cells;
		int* marching_mid_pts;		// middle point id

		void marchingCubes(std::vector<double> &mc_verts, std::vector<int> &mc_faces);

		void allocMarchingCubes()
		{
			num_mc_interface_cells = 0;
			mc_interface_cells = (int*)malloc(sizeof(int)*n);
			marching_mid_pts = (int*)malloc(sizeof(int)*n * 3);		//+x,+y,+z 3 for every grid
		}

		void deallocMarchingCubes()
		{
			free(mc_interface_cells);
			free(marching_mid_pts);
		}
		//////////////////marching cube variables}
		////////////////////////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////Render{

		//////////////////Render}
		////////////////////////////////////////////////////////////////////////////////////////////////////////
	};
}
#endif