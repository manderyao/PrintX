#include "TriMesh.h"
#include "Levelset.h"
#include "LevelsetMulti.h"
#include "meshOp.h"
using namespace PRINTX_LIB;

int main(int argc, char **argv)
{
	TriMesh mesh, mesh_in;

	bool success = mesh.loadSeg("../model/armadillo.seg", false);
	if (!success)
	{
		printf("Error: please run Demo_MeshSeg first to generate a .seg file or provide yours!\n");
		return -1;
	}

	printf("Successfully loaded ../model/armadillo.seg.\n");

	Levelset lvs_in(mesh, 2.0, true);
	lvs_in.redist(4.0);
	lvs_in.advectAllW(2.0);
	for (int i = 0; i < 30; i++)
	{
		lvs_in.smooth(0.3);
	}
	std::vector<double> mc_verts_in;
	std::vector<int> mc_faces_in;
	lvs_in.marchingCubes(mc_verts_in, mc_faces_in);
	mesh_in.loadMesh(mc_verts_in, mc_faces_in);
	printf("Built offset surface\n");

	LevelsetMulti lvs(mesh, 2.0, false);	
	for (int i = 0; i < 30; i++)
	{
		lvs.smooth(0.2);
	}
	std::vector<double> *mc_verts = new std::vector<double>[lvs.num_clusters];
	std::vector<int> *mc_faces = new std::vector<int>[lvs.num_clusters];
	lvs.marchingCubes(mc_verts, mc_faces);
	printf("Built piece enclosures\n");

	Nef_polyhedron nefpoly, nefpoly_in;
	loadNefPoly(nefpoly, mesh.verts, mesh.n_verts, mesh.faces, mesh.n_faces, true, true, true);
	printf("Loaded the surf mesh as a NefPoly\n");
	loadNefPoly(nefpoly_in, mesh_in.verts, mesh_in.n_verts, mesh_in.faces, mesh_in.n_faces, true, true, true);
	printf("Loaded the offset surf mesh as a NefPoly\n");
	for (int i = 0; i < lvs.num_clusters; i++)
	{
		Nef_polyhedron nefpoly1, nefpoly2;
		loadNefPoly(nefpoly1, &mc_verts[i][0], mc_verts[i].size() / 3, &mc_faces[i][0], mc_faces[i].size() / 3, true, true, true);
		printf("Loaded the piece enclosure as a NefPoly\n");
		boolOperator(nefpoly1, nefpoly, nefpoly2, BOOL_INTERSECT);
		boolOperator(nefpoly2, nefpoly_in, nefpoly1, BOOL_DIFF);

		if (cleanNefPoly(nefpoly1))
		{
			std::string dest = std::string("../model/armadillo_piece") + std::to_string(i) + ".obj";
			writeNefPoly2OBJ(nefpoly1, dest);
			printf("Write piece %d to %s\n", i, dest.c_str());
		}
		else
			printf("Fail to clean piece %d!!!\n", i);
	}

	delete[]mc_verts;
	delete[]mc_faces;
	return 0;
}
