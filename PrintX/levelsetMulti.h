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
//	Multiple level sets. All level set values are positive and IDs are used to indicate inside or outside
//	Implement the correction algorithm from the following paper:
//	Frank Losasso, Tamar Shinar, Andrew Selle, and Ronald Fedkiw. 2006. Multiple interacting liquids. ACM Trans. Graph. 25, 3 (July 2006), 812-819.
//	Based on the assumption that each time advection no more than one layer, trivial tricks are used to accelerate the level sets computations



#ifndef __PRINTX_LEVEL_SET_MULTI_H__
#define __PRINTX_LEVEL_SET_MULTI_H__

#include "triMesh.h"
#include <algorithm> 
#include "makeLevelset3.h"
#include <list>
#include "priorityQueueFM.h"
#include "vec3.h"


namespace PRINTX_LIB
{
	class LevelsetMulti
	{
	public:
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////special cells for acceleration{
		enum{ UNPROCESSED = 0, PROCESSING = 1, PROCESSED = 2 };
		enum{CELL_CHANGED=1};	// special markers
		// ID start from 0

		// keep a list of cells on interface
		int **interface_cells;
		int *num_interface_cells;
		// keep a list of changed cells for average
		int *changed_cells;
		int num_changed_cells;

		//	TODO: save space
		void allocSpecialCells()
		{
			interface_cells = (int**)malloc(sizeof(int*)*num_clusters);
			num_interface_cells = (int*)malloc(sizeof(int)*num_clusters);
			for (int i = 0; i < num_clusters; i++)
			{
				num_interface_cells[i] = 0;
				interface_cells[i] = (int*)malloc(sizeof(int)*(n));
			}

			changed_cells = (int*)malloc(sizeof(int)*(n));
			num_changed_cells = 0;
		}
		void deallocSpecialCells()
		{
			for (int i = 0; i < num_clusters; i++)
			{
				free(interface_cells[i]);
			}
			free(num_interface_cells);
			free(changed_cells);
		}

		inline void addChangedCell(int index)
		{
			mark[index] |= CELL_CHANGED;
			changed_cells[num_changed_cells++] = index;
		}
		inline bool isCellChanged(int index)
		{
			return (mark[index] & CELL_CHANGED) != 0;
		}
		void clearChangedCells()
		{
			for (int i = 0; i < num_changed_cells; i++)
			{
				mark[changed_cells[i]] &= ~CELL_CHANGED;
			}
			num_changed_cells = 0;
		}

		void buildInterfaceCells();	
		void clearInterfaceCells()
		{
			for (int i = 0; i < num_clusters; i++)
			{
				num_interface_cells[i] = 0;
			}
		}
		///////////////////special cells for acceleration}
		////////////////////////////////////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////help function{
		inline int getIndex(int i, int j, int k)
		{
			return k*nxy + j*gn[0] + i;
		}
		inline bool isInside(int index)
		{
			return ids[index] >= 0;
		}

		void getAdjCells(int index, int*qindex);				// neighbors: -x, +x, -y, +y, -z, +z; if out bound, =-1
		void getAdjCells2(int index, int*qindex);				// neighbors: -x, +x, -y, +y, -z, +z; if out bound, =index
		void getCubeCells(int x, int y, int z, int*qindex);		// the 8 corners of this cube

		//	merge connected component to remove fragments
		void mergeCntCmpnt();
		void mergeToNeighbor(int x, int y, int z, SDFGen::Array3uc &region_id_decided);	// merge to majority of neighbors
		///////////////////help function}
		////////////////////////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////Cell parameters{
		int gn[3];				//	number of cells each dimension
		int n;					//	total number 
		int nxy;				//	total number on the xy face, for acceleration
		double origin[3];
		double cell_size;
		int num_clusters;
		char *clusterMark;

		//	level sets
		double DIST_LIMIT;		// for those unprocessed cells
		int *ids;				// -1 for outside, 0,1,2... for piece domains
		int *id_tmp;
		double *dist_surf;		//distance to nearest surf point, not changed over process, always positive
		double **dists;			//distance to nearest interface point, one for each cluster, always positive
		double *dist_tmp;
		unsigned char **dist_tags;		//UNPROCESSED = 0, PROCESSING = 1, PROCESSED = 2
		int *mark;						//special marks
		double *curvs;					//curvature
		int *pos;						//3d indices

		// construct from a segmented mesh
		// enforceConnected: each piece is guaranteed to be a connected component by removing small fragments
		LevelsetMulti(TriMesh &mesh, double cellSize, bool enforceConnected, SDFGen::Array3i *closest_tri_ = NULL);
		~LevelsetMulti();

		inline double getDist(int index, int c)
		{
			return (ids[index] == c) ? ((dist_tags[c][index] == PROCESSED) ? -dists[c][index] : -DIST_LIMIT) : ((dist_tags[c][index] == PROCESSED) ? dists[c][index] : DIST_LIMIT);
		}
		void allocCells()
		{
			ids = (int*)malloc(sizeof(int)*n);
			id_tmp = (int*)malloc(sizeof(int)*n);
			dist_surf = (double*)malloc(sizeof(double)*n);
			mark = (int*)malloc(sizeof(int)*n);
			curvs = (double*)malloc(sizeof(double)*n);
			dist_tmp = (double*)malloc(sizeof(double)*n);
			dists = (double**)malloc(sizeof(double*)*num_clusters);
			dist_tags = (unsigned char**)malloc(sizeof(unsigned char*)*num_clusters);
			for (int i = 0; i < num_clusters; i++)
			{
				dists[i] = (double*)malloc(sizeof(double)*n);
				dist_tags[i] = (unsigned char *)malloc(sizeof(unsigned char)*n);
			}
			pos = (int*)malloc(sizeof(int) * 3 * n);

			clusterMark = (char*)malloc(sizeof(char)*num_clusters);;
		}
		void deallocCells()
		{
			for (int i = 0; i < num_clusters; i++)
			{
				free(dists[i]);
				free(dist_tags[i]);
			}
			free(dists);
			free(dist_tags);
			free(ids);
			free(id_tmp);
			free(dist_surf);
			free(mark);
			free(dist_tmp);
			free(curvs);
			free(pos);

			free(clusterMark);
		}
		void setAllDist(int idx, double dist)
		{
			for (int j = 0; j < num_clusters; j++)
			{
				dists[j][idx] = dist;
			}
		}
		void setAllDist(double dist)
		{
			for (int j = 0; j < num_clusters; j++)
			{
				for (int i = 0; i < n; i++)
				{
					dists[j][i] = dist;
				}
			}
		}
		void setAllDistTag(unsigned char dist_tag)
		{
			for (int j = 0; j < num_clusters; j++)
			{
				for (int i = 0; i < n; i++)
				{
					dist_tags[j][i] = dist_tag;
				}
			}
		}
		void getGrad(int index, int c, double *grad);
		void getApproxGrad(double x, double y, double z, int c, double *grad);
		double getDistSurf(double x, double y, double z);	// get surface level set
		double getDistSurfW(double x, double y, double z);	// input is world space coordinate
		double getDist(double x, double y, double z, int c);	// get level set
		double getDistW(double x, double y, double z, int c);	// input is world space coordinate
		int getID(double x, double y, double z);	// get ID
		int getIDW(double x, double y, double z);	// input is world space coordinate

		//Frank Losasso, Tamar Shinar, Andrew Selle, and Ronald Fedkiw. 2006. Multiple interacting liquids.In ACM SIGGRAPH 2006 Papers(SIGGRAPH '06)
		//return if id changed
		bool average(int index);
		///////////////////Cell parameters}
		////////////////////////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////Redistance{
		PriorityQueueFM<double> *fm_queue = 0;
		inline void allocRedist()
		{
			fm_queue = new PriorityQueueFM<double>(n);
		}
		inline void deallocRedist()
		{
			delete fm_queue;
		}
		void buildCurvature(int index, int c);
		double compute_extended_distance(int index, int c);

		// redist all cells until reaching bandwidth
		void redist(double band_width = 2.8);
		void redistProcess();
		void fastMarch(double band_width, int c);	// fast march

		void smooth(double curv_thres = 0.1, double speed = 0.2);
		void smoothCluster(int c, double curv_thres, double speed);
		///////////////////Redistance}
		////////////////////////////////////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////marching cube variables{
		int *mc_interface_cells;	// first:index; second: if on surface
		int num_mc_interface_cells;
		int* marching_mid_pts;		// middle point id

		void marchingCubes(std::vector<double> *mc_verts, std::vector<int> *mc_faces);
		void marchingCubes(int c, std::vector<double> &mc_verts, std::vector<int> &mc_faces);

		inline void allocMarchingCubes()
		{
			mc_interface_cells = (int*)malloc(sizeof(int)*n);
			marching_mid_pts = (int*)malloc(sizeof(int)*n * 3);		//+x,+y,+z 3 for every cell
		}

		inline void deallocMarchingCubes()
		{
			free(mc_interface_cells);
			free(marching_mid_pts);
		}
		///////////////////marching cube variables}
		////////////////////////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Render{
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		void renderInterfaceCells(int c);
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Render}
		////////////////////////////////////////////////////////////////////////////////////////////////////////
	};
}
#endif