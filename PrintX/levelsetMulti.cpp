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
//	Assume each time advection no more than one layer
//**************************************************************************************

#include "levelsetMulti.h"
#include "cube.h"
#include "opgl.h"
#include "colorManager.h"
#include <queue>
#include <unordered_map>
#include <algorithm>

using namespace PRINTX_LIB;

// find the smallest phi1 and phi2, shift all phis until phi1+phi2=0
bool LevelsetMulti::average(int index)
{
	int minid1 = 0, minid2 = 1;
	double mindist1 = getDist(index, 0), mindist2 = getDist(index, 1), dist;
	if (mindist2 < mindist1)
	{
		std::swap(minid1, minid2);
		std::swap(mindist1, mindist2);
	}

	for (int i = 2; i < this->num_clusters; i++)
	{
		dist = getDist(index, i);
		if (dist < mindist1)
		{
			minid2 = minid1;
			minid1 = i;
			mindist2 = mindist1;
			mindist1 = dist;
		}
		else if (dist < mindist2)
		{
			minid2 = i;
			mindist2 = dist;
		}
	}

	bool change = (ids[index] != minid1);
	double anchor = 0.5*(mindist1 + mindist2);

	// shift
	for (int i = 0; i < this->num_clusters; i++)
	{
		dists[i][index] = getDist(index, i) - anchor;
	}


	dists[minid1][index] = -dists[minid1][index];
	ids[index] = minid1;
	return change;
}

LevelsetMulti::LevelsetMulti(TriMesh &mesh, double cellSize, bool enforceConnected, SDFGen::Array3i *closest_tri_)
{
	if (!mesh.hasSegments)
	{
		printf("Error Constructing LevelsetMulti: Mesh is not segmented!!!\n");
		exit(-1);
	}

	std::vector<SDFGen::Vec3f> vertList;
	std::vector<SDFGen::Vec3ui> faceList;
	for (int v = 0; v < mesh.n_verts; v++)
	{
		vertList.push_back(SDFGen::Vec3f(mesh.verts[3 * v], mesh.verts[3 * v + 1], mesh.verts[3 * v + 2]));
	}
	for (int f = 0; f < mesh.n_faces; f++)
	{
		faceList.push_back(SDFGen::Vec3ui(mesh.faces[3 * f], mesh.faces[3 * f + 1], mesh.faces[3 * f + 2]));
	}

	double bmax[3], bmin[3];
	mesh.buildBbox(bmin, bmax);
	SDFGen::Vec3f minBox(bmin), maxBox(bmax);


	//Add padding around the box.
	int padding = 2;
	float dx = cellSize;
	SDFGen::Vec3f unit(1, 1, 1);
	minBox -= padding*dx*unit;
	maxBox += padding*dx*unit;
	SDFGen::Vec3ui sizes = SDFGen::Vec3ui((maxBox - minBox) / dx);
	sizes[0]++; sizes[1]++; sizes[2]++;


	gn[0] = sizes[0]; gn[1] = sizes[1]; gn[2] = sizes[2]; n = gn[0] * gn[1] * gn[2]; nxy = gn[0] * gn[1];
	origin[0] = minBox[0]; origin[1] = minBox[1]; origin[2] = minBox[2];
	cell_size = cellSize;

	num_clusters = mesh.patches.size();


	//generate distance field
	SDFGen::Array3f phi(sizes[0], sizes[1], sizes[2]);
	SDFGen::Array3vec3f pt(sizes[0], sizes[1], sizes[2]);
	SDFGen::Array3i closest_tri;
	make_level_set3(faceList, vertList, minBox, dx, sizes[0], sizes[1], sizes[2], phi, closest_tri, pt);

	// return closest_tri if requested
	if (closest_tri_)
		*closest_tri_ = closest_tri;

	allocCells();

	//	initialize level sets
	for (int k = 0; k < gn[2]; k++)
		for (int j = 0; j < gn[1]; j++)
			for (int i = 0; i < gn[0]; i++)
			{
		int index = getIndex(i, j, k);
		ids[index] = mesh.face_ids[closest_tri(i, j, k)];
		dist_surf[index] = (phi(i, j, k) / cell_size);

		setAllDist(index, fabs(dist_surf[index]));	// this dist doesn't matter 
		mark[index] = 0;

		pos[3 * index] = i; pos[3 * index + 1] = j; pos[3 * index + 2] = k;
			}

	DIST_LIMIT = 10000;

	if (enforceConnected)
		mergeCntCmpnt();

	setAllDistTag(UNPROCESSED);

	allocSpecialCells();
	allocRedist();
	allocMarchingCubes();

	redist();
}

LevelsetMulti::~LevelsetMulti()
{
	deallocSpecialCells();
	deallocCells();
	deallocRedist();
	deallocMarchingCubes();
}


void LevelsetMulti::buildInterfaceCells()
{
	clearInterfaceCells();

	int qidx[6];

	for (int idx = 0; idx < n; idx++)
	{
		memset(clusterMark, 0, sizeof(clusterMark[0])*num_clusters);

		int id0 = ids[idx];
		getAdjCells(idx, qidx);
		for (int ii = 0; ii < 6; ii++)
		{
			if (qidx[ii] >= 0 && (ids[qidx[ii]]) != id0)
			{
				if (!clusterMark[id0])
				{
					clusterMark[id0] = 1;
					interface_cells[id0][num_interface_cells[id0]++] = idx;
				}
				if (!clusterMark[ids[qidx[ii]]])
				{
					clusterMark[ids[qidx[ii]]] = 1;
					interface_cells[ids[qidx[ii]]][num_interface_cells[ids[qidx[ii]]]++] = idx;
				}
			}
		}
	}
}

int LevelsetMulti::getIDW(double x, double y, double z)
{
	return getID((x - origin[0]) / cell_size, (y - origin[1]) / cell_size, (z - origin[2]) / cell_size);
}


int LevelsetMulti::getID(double x, double y, double z)
{
	double min_dist;
	int min_id = -1;

	for (int i = 0; i < num_clusters; i++)
	{
		double d_ = getDist(x, y, z, i);
		if (min_id < 0 || d_ < min_dist)
		{
			min_id = i;
			min_dist = d_;
		}
	}
	return min_id;
}


double LevelsetMulti::getDistW(double x, double y, double z, int cluster)
{

	return getDist((x - origin[0]) / cell_size, (y - origin[1]) / cell_size, (z - origin[2]) / cell_size, cluster);
}

double LevelsetMulti::getDist(double x, double y, double z, int cluster)
{
	int i = std::floor(x), j = std::floor(y), k = std::floor(z);

	double a = x - i, b = y - j, c = z - k;

	int qindex[8];
	getCubeCells(i, j, k, qindex);

	return	(1 - a)*(1 - b)*(1 - c)*getDist(qindex[0], cluster) + a*(1 - b)*(1 - c)*getDist(qindex[3], cluster) +
		(1 - a)*b*(1 - c)*getDist(qindex[4], cluster) + a*b*(1 - c)*getDist(qindex[7], cluster) +
		(1 - a)*(1 - b)*c*getDist(qindex[1], cluster) + a*(1 - b)*c*getDist(qindex[2], cluster) +
		(1 - a)*b*c*getDist(qindex[5], cluster) + a*b*c*getDist(qindex[6], cluster);
}

double LevelsetMulti::getDistSurfW(double x, double y, double z)
{
	return getDistSurf((x - origin[0]) / cell_size, (y - origin[1]) / cell_size, (z - origin[2]) / cell_size);
}

double LevelsetMulti::getDistSurf(double x, double y, double z)
{
	int i = std::floor(x), j = std::floor(y), k = std::floor(z);

	double a = x - i, b = y - j, c = z - k;

	int qindex[8];
	getCubeCells(i, j, k, qindex);

	return	(1 - a)*(1 - b)*(1 - c)*dist_surf[qindex[0]] + a*(1 - b)*(1 - c)*dist_surf[qindex[3]] +
		(1 - a)*b*(1 - c)*dist_surf[qindex[4]] + a*b*(1 - c)*dist_surf[qindex[7]] +
		(1 - a)*(1 - b)*c*dist_surf[qindex[1]] + a*(1 - b)*c*dist_surf[qindex[2]] +
		(1 - a)*b*c*dist_surf[qindex[5]] + a*b*c*dist_surf[qindex[6]];
}

void LevelsetMulti::getGrad(int index, int c, double *grad)
{
	int qindex[6];
	getAdjCells2(index, qindex);

	grad[0] = (getDist(qindex[1], c) - getDist(qindex[0], c)) / 2;
	grad[1] = (getDist(qindex[3], c) - getDist(qindex[2], c)) / 2;
	grad[2] = (getDist(qindex[5], c) - getDist(qindex[4], c)) / 2;
}

void LevelsetMulti::getApproxGrad(double x, double y, double z, int c, double *grad)
{
	int i = std::round(x), j = std::round(y), k = std::round(z);
	getGrad(getIndex(i, j, k), c, grad);
}

void LevelsetMulti::getAdjCells2(int index, int*qindex)
{
	int *p = pos + 3 * index;
	qindex[0] = (p[0] - 1) >= 0 ? (index - 1) : index;
	qindex[1] = (p[0] + 1) < gn[0] ? (index + 1) : index;
	qindex[2] = (p[1] - 1) >= 0 ? (index - gn[0]) : index;
	qindex[3] = (p[1] + 1) < gn[1] ? (index + gn[0]) : index;
	qindex[4] = (p[2] - 1) >= 0 ? (index - nxy) : index;
	qindex[5] = (p[2] + 1) < gn[2] ? (index + nxy) : index;
}

void LevelsetMulti::getAdjCells(int index, int*qindex)
{
	int *p = pos + 3 * index;
	qindex[0] = (p[0] - 1) >= 0 ? (index - 1) : -1;
	qindex[1] = (p[0] + 1) < gn[0] ? (index + 1) : -1;
	qindex[2] = (p[1] - 1) >= 0 ? (index - gn[0]) : -1;
	qindex[3] = (p[1] + 1) < gn[1] ? (index + gn[0]) : -1;
	qindex[4] = (p[2] - 1) >= 0 ? (index - nxy) : -1;
	qindex[5] = (p[2] + 1) < gn[2] ? (index + nxy) : -1;
}


// the 8 corners of this cube
void LevelsetMulti::getCubeCells(int x, int y, int z, int*qindex)
{
	for (int i = 0; i < 8; i++)
	{
		qindex[i] = getIndex(x + CUBE_STRUCT::vert_to_idx[i][0], y + CUBE_STRUCT::vert_to_idx[i][1], z + CUBE_STRUCT::vert_to_idx[i][2]);
	}
}


// normalize level sets on the interface
void LevelsetMulti::redistProcess()
{
	int index, qindex[6];
	double diff;
	int id, qid;


	for (int c = 0; c < num_clusters; c++)
	{
		// initialize
		for (int i = 0; i < num_interface_cells[c]; i++)
		{
			index = interface_cells[c][i];
			dist_tmp[index] = dists[c][index];
		}

		for (int i = 0; i < num_interface_cells[c]; i++)
		{
			index = interface_cells[c][i];
			if (ids[index] == c)
			{
				getAdjCells(index, qindex);
				for (int ii = 0; ii < 6; ii++)
				{
					if (qindex[ii] >= 0 && ids[qindex[ii]] != c && (diff = dists[c][index] + dists[c][qindex[ii]]) > 1)
					{
						double dist_ = dists[c][index] / diff;
						dist_tmp[index] = std::min(dist_tmp[index], dist_);

						dist_ = dists[c][qindex[ii]] / diff;
						dist_tmp[qindex[ii]] = std::min(dist_tmp[qindex[ii]], dist_);
					}
				}
			}
		}

		// restore
		for (int i = 0; i < num_interface_cells[c]; i++)
		{
			index = interface_cells[c][i];
			dists[c][index] = dist_tmp[index];
		}
	}

}


static double extended_distance_2_value(double theta1, double theta2)
{
	if (theta1 > 9999) return 1 + theta2;
	if (theta2 > 9999) return 1 + theta1;

	if (abs(theta1 - theta2) > 1)
	{
		if (theta1 > theta2) return 1 + theta2;
		else return 1 + theta1;
	}
	double a = 2;
	double b = -2 * (theta1 + theta2);
	double c = theta1*theta1 + theta2*theta2 - 1;
	return (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
}


static double extended_distance_3_value(double theta1, double theta2, double theta3)
{
	if (theta1 > 9999) return extended_distance_2_value(theta2, theta3);
	if (theta2 > 9999) return extended_distance_2_value(theta1, theta3);
	if (theta3 > 9999) return extended_distance_2_value(theta1, theta2);

	int max_index = VECTOR_MAX_INDEX3S(theta1, theta2, theta3);
	double max_theta = std::max(theta1, std::max(theta2, theta3));
	if ((max_theta - theta1)*(max_theta - theta1) + (max_theta - theta2)*(max_theta - theta2) + (max_theta - theta3)*(max_theta - theta3) > 1)
	{
		if (max_index == 0)		return extended_distance_2_value(theta2, theta3);
		else if (max_index == 1)		return extended_distance_2_value(theta1, theta3);
		else		return extended_distance_2_value(theta1, theta2);
	}
	double a = 3;
	double b = -2 * (theta1 + theta2 + theta3);
	double c = theta1*theta1 + theta2*theta2 + theta3*theta3 - 1;
	return (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
}

double LevelsetMulti::compute_extended_distance(int index, int c)
{
	double theta1 = 10000, theta2 = 10000, theta3 = 10000;
	int qindex[6];

	getAdjCells(index, qindex);

	// X direction	
	if (qindex[0] >= 0 && dist_tags[c][qindex[0]] == PROCESSED)
	{
		theta1 = dists[c][qindex[0]];
	}
	if (qindex[1] >= 0 && dist_tags[c][qindex[1]] == PROCESSED)
	{
		theta1 = std::min(theta1, dists[c][qindex[1]]);
	}

	// Y direction	
	if (qindex[2] >= 0 && dist_tags[c][qindex[2]] == PROCESSED)
	{
		theta2 = dists[c][qindex[2]];
	}
	if (qindex[3] >= 0 && dist_tags[c][qindex[3]] == PROCESSED)
	{
		theta2 = std::min(theta2, dists[c][qindex[3]]);
	}

	// Z direction	
	if (qindex[4] >= 0 && dist_tags[c][qindex[4]] == PROCESSED)
	{
		theta3 = dists[c][qindex[4]];
	}
	if (qindex[5] >= 0 && dist_tags[c][qindex[5]] == PROCESSED)
	{
		theta3 = std::min(theta2, dists[c][qindex[5]]);
	}

#if 0
	if (theta1 < 0 || theta2 < 0 || theta3 < 0)
		printf("ERROR: in compute_extended_distance, theta < 0");
#endif
	return extended_distance_3_value(theta1, theta2, theta3);

}

void LevelsetMulti::fastMarch(double band_width, int c)
{
	fm_queue->setDistsPointer(dists[c]);

	int index;
	int qindex[6];
	for (int i = 0; i < num_interface_cells[c]; i++)
	{
		index = interface_cells[c][i];
		fm_queue->push(index);
		dist_tags[c][index] = PROCESSED;
	}

	while (!fm_queue->empty())
	{
		index = fm_queue->pop();
		dist_tags[c][index] = PROCESSED;

		if (band_width <= dists[c][index])	//exceed bandwidth, terminate
			continue;

		getAdjCells(index, qindex);

		for (int i = 0; i < 6; i++)
		{
			if (qindex[i] >= 0)
			{
				// calculate extended distance
				if (dist_tags[c][qindex[i]] == UNPROCESSED)
				{
					double dist = compute_extended_distance(qindex[i], c);
					dists[c][qindex[i]] = dist;
					dist_tags[c][qindex[i]] = PROCESSING;
					fm_queue->push(qindex[i]);
				}
				else if (dist_tags[c][qindex[i]] == PROCESSING)
				{
					double dist = compute_extended_distance(qindex[i], c);

					if (dists[c][qindex[i]] - dist > 1e-6)
					{
						fm_queue->update_dist(qindex[i], dist);
					}
				}
			}
		}
	}
}

void LevelsetMulti::redist(double band_width)
{
	// clear process marks
	for (int  c= 0; c < num_clusters; c++)
	{
		memset(dist_tags[c], 0, sizeof(dist_tags[0][0])*n);// set to UNPROCESSED
	}

	buildInterfaceCells();

	redistProcess();	// normalize level sets on the interface
	for (int c = 0; c < num_clusters; c++)
	{
		fastMarch(band_width, c);
	}
}


void LevelsetMulti::buildCurvature(int index, int c)
{
	int *p = pos + 3 * index;
	int i = p[0], j = p[1], k = p[2];

	int i_0 = i > 0 ? i - 1 : 0;
	int i_1 = i < gn[0] - 1 ? i + 1 : gn[0] - 1;
	int j_0 = j > 0 ? j - 1 : 0;
	int j_1 = j < gn[1] - 1 ? j + 1 : gn[1] - 1;
	int k_0 = k > 0 ? k - 1 : 0;
	int k_1 = k < gn[2] - 1 ? k + 1 : gn[2] - 1;

	double center = getDist(index, c);
	double x0 = getDist(getIndex(i_0, j, k), c);
	double x1 = getDist(getIndex(i_1, j, k), c);
	double y0 = getDist(getIndex(i, j_0, k), c);
	double y1 = getDist(getIndex(i, j_1, k), c);
	double z0 = getDist(getIndex(i, j, k_0), c);
	double z1 = getDist(getIndex(i, j, k_1), c);

	double theta_x = (x1 - x0)*0.5;
	double theta_y = (y1 - y0)*0.5;
	double theta_z = (z1 - z0)*0.5;

	double theta_x2 = theta_x*theta_x;
	double theta_y2 = theta_y*theta_y;
	double theta_z2 = theta_z*theta_z;

	double theta_xx = (x0 + x1 - 2 * center);
	double theta_yy = (y0 + y1 - 2 * center);
	double theta_zz = (z0 + z1 - 2 * center);

	double theta_xy = (getDist(getIndex(i_1, j_1, k), c) - getDist(getIndex(i_0, j_1, k), c) - getDist(getIndex(i_1, j_0, k), c) + getDist(getIndex(i_0, j_0, k), c))*0.25;
	double theta_xz = (getDist(getIndex(i_1, j, k_1), c) - getDist(getIndex(i_0, j, k_1), c) - getDist(getIndex(i_1, j, k_0), c) + getDist(getIndex(i_0, j, k_0), c))*0.25;
	double theta_yz = (getDist(getIndex(i, j_1, k_1), c) - getDist(getIndex(i, j_0, k_1), c) - getDist(getIndex(i, j_1, k_0), c) + getDist(getIndex(i, j_0, k_0), c))*0.25;

#if 1 //mean curvature
	double result =
		theta_x2*theta_yy + theta_y2*theta_xx - 2 * theta_x*theta_y*theta_xy
		+ theta_x2*theta_zz + theta_z2*theta_xx - 2 * theta_x*theta_z*theta_xz
		+ theta_y2*theta_zz + theta_z2*theta_yy - 2 * theta_y*theta_z*theta_yz;

	float gradient_squared = theta_x2 + theta_y2 + theta_z2 + 1e-16;
	float gradient_length = sqrt(gradient_squared);

	curvs[index] = result / (gradient_squared*gradient_length);
#else	//Surface diffusion
	float N = (theta_x2*theta_xx + theta_y2*theta_yy + theta_z2*theta_zz + 2 * theta_x*theta_y*theta_xy + 2 * theta_x*theta_z*theta_xz + 2 * theta_z*theta_y*theta_yz) / (theta_x2 + theta_y2 + theta_z2 + 1e-8);
	float du = theta_xx + theta_yy + theta_zz;
	return du - N;
#endif
}

void LevelsetMulti::smoothCluster(int c, double curv_thres, double speed)
{
	int index;
	double diff;
	for (int i = 0; i < num_interface_cells[c]; i++)
	{
		index = interface_cells[c][i];
		buildCurvature(interface_cells[c][i], c);

		if (curvs[index] > curv_thres)
		{
			diff = curvs[index] - curv_thres;

			if (c == ids[index])
			{
				dists[c][index] = dists[c][index] - speed*diff;
			}
			else
				dists[c][index] = dists[c][index] + speed*diff;

			if (!isCellChanged(index))
				addChangedCell(index);
		}
		else if (curvs[index] < -curv_thres)
		{
			diff = curvs[index] + curv_thres;

			if (c == ids[index])
			{
				dists[c][index] = dists[c][index] - speed*diff;
			}
			else
				dists[c][index] = dists[c][index] + speed*diff;

			if (!isCellChanged(index))
				addChangedCell(index);
		}
	}
}


void LevelsetMulti::smooth(double curv_thres, double speed)
{
	clearChangedCells();

	for (int c = 0; c < num_clusters; c++)
	{
		smoothCluster(c, curv_thres, speed);
	}

	for (int i = 0; i < num_changed_cells; i++)
	{
		average(changed_cells[i]);	// correction
	}

	clearChangedCells();

	redist();
}

// connected component with a specific ID
class ID_Grup_MULTI
{
public:
	int id;
	std::vector<SDFGen::Vec3i> grup;
	ID_Grup_MULTI()
	{
		id = -1;
	}
	ID_Grup_MULTI(int id_, std::vector<SDFGen::Vec3i> &grup_)
	{
		id = id_;
		grup = grup_;
	}
	bool operator<(const ID_Grup_MULTI& other) const
	{
		return id < other.id || (id == other.id && grup.size() < other.grup.size());
	}
};

void LevelsetMulti::mergeCntCmpnt()
{
	SDFGen::Array3uc region_id_decided(gn[0], gn[1], gn[2], (unsigned char)0);

	std::set<ID_Grup_MULTI> id_grups;
	std::vector<SDFGen::Vec3i> id_grup;
	std::queue<SDFGen::Vec3i> bfs;
	int id_;

	for (int i = 0; i < gn[0]; i++)
	{
		for (int j = 0; j < gn[1]; j++)
		{
			for (int k = 0; k < gn[2]; k++)
			{
				if (!region_id_decided(i, j, k))
				{
					// BFS to extract connected component with this ID

					bfs = std::queue<SDFGen::Vec3i>();
					bfs.push(SDFGen::Vec3i(i, j, k));
					id_ = ids[getIndex(i, j, k)];
					region_id_decided(i, j, k) = 1;
					id_grup.clear();
					id_grup.push_back(SDFGen::Vec3i(i, j, k));

					int cnt = 1;
					while (!bfs.empty())
					{
						SDFGen::Vec3i seed = bfs.front();
						bfs.pop();

						if (seed[0] - 1 >= 0 && !region_id_decided(seed[0] - 1, seed[1], seed[2]) && ids[getIndex(seed[0] - 1, seed[1], seed[2])] == id_)
						{
							bfs.push(SDFGen::Vec3i(seed[0] - 1, seed[1], seed[2]));
							region_id_decided(seed[0] - 1, seed[1], seed[2]) = 1;
							cnt++;
							id_grup.push_back(SDFGen::Vec3i(seed[0] - 1, seed[1], seed[2]));
						}
						if (seed[0] + 1 < gn[0] && !region_id_decided(seed[0] + 1, seed[1], seed[2]) && ids[getIndex(seed[0] + 1, seed[1], seed[2])] == id_)
						{
							bfs.push(SDFGen::Vec3i(seed[0] + 1, seed[1], seed[2]));
							region_id_decided(seed[0] + 1, seed[1], seed[2]) = 1;
							cnt++;
							id_grup.push_back(SDFGen::Vec3i(seed[0] + 1, seed[1], seed[2]));
						}
						if (seed[1] - 1 >= 0 && !region_id_decided(seed[0], seed[1] - 1, seed[2]) && ids[getIndex(seed[0], seed[1] - 1, seed[2])] == id_)
						{
							bfs.push(SDFGen::Vec3i(seed[0], seed[1] - 1, seed[2]));
							region_id_decided(seed[0], seed[1] - 1, seed[2]) = 1;
							cnt++;
							id_grup.push_back(SDFGen::Vec3i(seed[0], seed[1] - 1, seed[2]));
						}
						if (seed[1] + 1 < gn[1] && !region_id_decided(seed[0], seed[1] + 1, seed[2]) && ids[getIndex(seed[0], seed[1] + 1, seed[2])] == id_)
						{
							bfs.push(SDFGen::Vec3i(seed[0], seed[1] + 1, seed[2]));
							region_id_decided(seed[0], seed[1] + 1, seed[2]) = 1;
							cnt++;
							id_grup.push_back(SDFGen::Vec3i(seed[0], seed[1] + 1, seed[2]));
						}
						if (seed[2] - 1 >= 0 && !region_id_decided(seed[0], seed[1], seed[2] - 1) && ids[getIndex(seed[0], seed[1], seed[2] - 1)] == id_)
						{
							bfs.push(SDFGen::Vec3i(seed[0], seed[1], seed[2] - 1));
							region_id_decided(seed[0], seed[1], seed[2] - 1) = 1;
							cnt++;
							id_grup.push_back(SDFGen::Vec3i(seed[0], seed[1], seed[2] - 1));
						}
						if (seed[2] + 1 < gn[2] && !region_id_decided(seed[0], seed[1], seed[2] + 1) && ids[getIndex(seed[0], seed[1], seed[2] + 1)] == id_)
						{
							bfs.push(SDFGen::Vec3i(seed[0], seed[1], seed[2] + 1));
							region_id_decided(seed[0], seed[1], seed[2] + 1) = 1;
							cnt++;
							id_grup.push_back(SDFGen::Vec3i(seed[0], seed[1], seed[2] + 1));
						}
					}
					id_grups.insert(ID_Grup_MULTI(id_, id_grup));
				}
			}
		}
	}

	// Only keep the largest component for each ID
	region_id_decided.assign(1);
	std::set<ID_Grup_MULTI>::const_iterator iter, iter_next;
	for (iter = id_grups.begin(); iter != id_grups.end(); iter++)
	{
		iter_next = std::next(iter);
		if (iter_next != id_grups.end() && iter->id == iter_next->id)
		{
			// find largest connected component for each id
			for (std::vector<SDFGen::Vec3i>::const_iterator vec_iter = iter->grup.begin(); vec_iter != iter->grup.end(); vec_iter++)
			{
				region_id_decided((*vec_iter)[0], (*vec_iter)[1], (*vec_iter)[2]) = 0;
			}
		}
	}

	//	flood filled for the rest cells
	int cnt;
	while (true)
	{
		cnt = 0;
		for (int i = 0; i < gn[0]; i++)
		{
			for (int j = 0; j < gn[1]; j++)
			{
				for (int k = 0; k < gn[2]; k++)
				{
					if (!region_id_decided(i, j, k))
					{
						mergeToNeighbor(i, j, k, region_id_decided);
						cnt++;
					}

				}
			}
		}
		if (cnt == 0)
			break;
	}

}

void LevelsetMulti::mergeToNeighbor(int x, int y, int z, SDFGen::Array3uc &region_id_decided)
{
	typedef std::pair<int, int> Pair;


	std::unordered_map<int, int> map;
	int id_;
	int max_id = -1;
	int max_n = 0;
	if (x - 1 >= 0 && region_id_decided(x - 1, y, z))
	{
		id_ = ids[getIndex(x - 1, y, z)];
		if (map.find(id_) == map.end())
			map.insert(Pair(id_, 1));
		else
			map[id_]++;

		if (map[id_] > max_n)
		{
			max_n = map[id_];
			max_id = id_;
		}
	}
	if (x + 1 < gn[0] && region_id_decided(x + 1, y, z))
	{
		id_ = ids[getIndex(x + 1, y, z)];
		if (map.find(id_) == map.end())
			map.insert(Pair(id_, 1));
		else
			map[id_]++;

		if (map[id_] > max_n)
		{
			max_n = map[id_];
			max_id = id_;
		}
	}
	if (y - 1 >= 0 && region_id_decided(x, y - 1, z))
	{
		id_ = ids[getIndex(x, y - 1, z)];
		if (map.find(id_) == map.end())
			map.insert(Pair(id_, 1));
		else
			map[id_]++;

		if (map[id_] > max_n)
		{
			max_n = map[id_];
			max_id = id_;
		}
	}
	if (y + 1 < gn[1] && region_id_decided(x, y + 1, z))
	{
		id_ = ids[getIndex(x, y + 1, z)];
		if (map.find(id_) == map.end())
			map.insert(Pair(id_, 1));
		else
			map[id_]++;

		if (map[id_] > max_n)
		{
			max_n = map[id_];
			max_id = id_;
		}
	}
	if (z - 1 >= 0 && region_id_decided(x, y, z - 1))
	{
		id_ = ids[getIndex(x, y, z - 1)];
		if (map.find(id_) == map.end())
			map.insert(Pair(id_, 1));
		else
			map[id_]++;

		if (map[id_] > max_n)
		{
			max_n = map[id_];
			max_id = id_;
		}
	}
	if (z + 1 < gn[2] && region_id_decided(x, y, z + 1))
	{
		id_ = ids[getIndex(x, y, z + 1)];
		if (map.find(id_) == map.end())
			map.insert(Pair(id_, 1));
		else
			map[id_]++;

		if (map[id_] > max_n)
		{
			max_n = map[id_];
			max_id = id_;
		}
	}

	ids[getIndex(x, y, z)] = max_id;
	region_id_decided(x, y, z) = 1;
}


void LevelsetMulti::marchingCubes(int c, std::vector<double> &mc_verts, std::vector<int> &mc_faces)
{
	// set mc interface cells
	num_mc_interface_cells = 0;
	int qindex[8];
	for (int i = 0; i < gn[0] - 1; i++)
	{
		for (int j = 0; j < gn[1] - 1; j++)
		{
			for (int k = 0; k < gn[2] - 1; k++)
			{
				int index = getIndex(i, j, k);
				getCubeCells(i, j, k, qindex);
				if (ids[index] != c)
				{
					for (int t = 1; t < 8; t++)
					{
						if (ids[qindex[t]] == c)
						{
							mc_interface_cells[num_mc_interface_cells++] = index;
							break;
						}
					}
				}
				else
				{
					for (int t = 1; t < 8; t++)
					{
						if (ids[qindex[t]] != c)
						{
							mc_interface_cells[num_mc_interface_cells++] = index;
							break;
						}
					}
				}
			}
		}
	}

	// init
	for (int i = 0; i < num_mc_interface_cells; i++)
	{
		int index = mc_interface_cells[i];
		marching_mid_pts[3 * index] = -1;
		marching_mid_pts[3 * index + 1] = -1;
		marching_mid_pts[3 * index + 2] = -1;		//	middle points not found yet
	}


	mc_verts.clear();
	mc_faces.clear();
	for (int i = 0; i < num_mc_interface_cells; i++)
	{
		int index = mc_interface_cells[i];
		int *gp = pos + 3 * index, idx[3], mid_idx[12];
		getCubeCells(gp[0], gp[1], gp[2], qindex);

		for (int t = 0; t < 12; t++)	// traverse each edge
		{
			idx[0] = qindex[CUBE_STRUCT::edge_to_vert[t][0]];
			idx[1] = qindex[CUBE_STRUCT::edge_to_vert[t][1]];
			idx[2] = CUBE_STRUCT::edge_to_vert[t][2];			// direction
			int id1 = ids[idx[0]];
			int id2 = ids[idx[1]];
			if ((id1 == c && id2 != c) || (id1 != c && id2 == c))
			{
				if (marching_mid_pts[3 * idx[0] + idx[2]] >= 0)
				{
					mid_idx[t] = marching_mid_pts[3 * idx[0] + idx[2]];	// calculated already
				}
				else
				{
					double p0 = dists[c][idx[0]];
					double p1 = dists[c][idx[1]];
					double p01 = p0 / (p0 + p1 + 1e-16);
					int *gp1 = pos + 3 * idx[0], *gp2 = pos + 3 * idx[1];
					double p[3];
					VECTOR_SCALE2_ADD3(1 - p01, gp1, p01, gp2, p);
					mc_verts.push_back(origin[0] + cell_size*p[0]);
					mc_verts.push_back(origin[1] + cell_size*p[1]);
					mc_verts.push_back(origin[2] + cell_size*p[2]);
					marching_mid_pts[3 * idx[0] + idx[2]] = mc_verts.size() / 3 - 1;
					mid_idx[t] = mc_verts.size() / 3 - 1;
				}
			}
			else
				mid_idx[t] = -1;
		}

		int type = 0;
		for (int t = 0; t < 8; t++)
		{
			type += (ids[qindex[t]] == c) ? 1 << t : 0;
		}

		double norm[3];
		for (int j = 0; j < CUBE_STRUCT::case_to_numpolys[type]; j++)
		{
			VectorSet3i(idx, CUBE_STRUCT::edge_connect_list[type][j]);
			mc_faces.push_back(mid_idx[idx[2]]);
			mc_faces.push_back(mid_idx[idx[1]]);
			mc_faces.push_back(mid_idx[idx[0]]);
		}
	}
}

void LevelsetMulti::marchingCubes(std::vector<double> *mc_verts, std::vector<int> *mc_faces)
{
	int index;

	// build a boundary
	for (int i = 0; i < gn[0]; i++)for (int j = 0; j < gn[1]; j++)
	{
		index = getIndex(i,j,0);
		id_tmp[index] = ids[index];
		ids[index] = -1;
		index = getIndex(i, j, gn[2]-1);
		id_tmp[index] = ids[index];
		ids[index] = -1;
	}
	for (int i = 0; i < gn[0]; i++)for (int k = 1; k < gn[2]-1; k++)
	{
		index = getIndex(i, 0, k);
		id_tmp[index] = ids[index];
		ids[index] = -1;
		index = getIndex(i, gn[1]-1, k);
		id_tmp[index] = ids[index];
		ids[index] = -1;
	}
	for (int j = 1; j < gn[1]-1; j++)for (int k = 1; k < gn[2]-1; k++)
	{
		index = getIndex(0, j, k);
		id_tmp[index] = ids[index];
		ids[index] = -1;
		index = getIndex(gn[0]-1, j, k);
		id_tmp[index] = ids[index];
		ids[index] = -1;
	}


	// marching cubes
	for (int i = 0; i < num_clusters; i++)
	{
		marchingCubes(i, mc_verts[i], mc_faces[i]);
	}

	//restore boundary
	for (int i = 0; i < gn[0]; i++)for (int j = 0; j < gn[1]; j++)
	{
		index = getIndex(i, j, 0);
		ids[index] = id_tmp[index];
		index = getIndex(i, j, gn[2] - 1);
		ids[index] = id_tmp[index];
	}
	for (int i = 0; i < gn[0]; i++)for (int k = 1; k < gn[2] - 1; k++)
	{
		index = getIndex(i, 0, k);
		ids[index] = id_tmp[index];
		index = getIndex(i, gn[1] - 1, k);
		ids[index] = id_tmp[index];
	}
	for (int j = 1; j < gn[1] - 1; j++)for (int k = 1; k < gn[2] - 1; k++)
	{
		index = getIndex(0, j, k);
		ids[index] = id_tmp[index];
		index = getIndex(gn[0] - 1, j, k);
		ids[index] = id_tmp[index];
	}
}


void LevelsetMulti::renderInterfaceCells(int cluster)
{
	glPointSize(1.5f);
	glBegin(GL_POINTS);


	if (cluster<0)
	{
		for (int c = 0; c < num_clusters; c++)
		{
			glColor3fv(renderColor + 3 * c);
			for (int i = 0; i < num_interface_cells[c]; i++)
			{
				int index = interface_cells[c][i];
				glVertex3f((float)pos[3 * index], (float)pos[3 * index + 1], (float)pos[3 * index + 2]);
			}
		}
	}
	else
	{
		int c = cluster;
		glColor3fv(renderColor + 3 * c);
		for (int i = 0; i < num_interface_cells[c]; i++)
		{
			int index = interface_cells[c][i];
			glVertex3f((float)pos[3 * index], (float)pos[3 * index + 1], (float)pos[3 * index + 2]);
		}
	}
	glEnd();
}