#include "TriMesh.h"
#include "LevelsetMulti.h"

using namespace PRINTX_LIB;

TriMesh mesh;

int main(int argc, char **argv)
{
	bool success = mesh.loadSeg("../model/armadillo.seg", false);
	if (!success)
	{
		printf("Error: please run Demo_MeshSeg first to generate a .seg file or provide yours!\n");
		return -1;
	}

	LevelsetMulti lvs(mesh, 1.5, false);	
	printf("Level set initialized\n");
	for (int i = 0; i < 30; i++)
	{	
		lvs.smooth(0.2);
		printf("Level set smooth Iter: %d\n", i);
	}

	std::vector<double> *mc_verts = new std::vector<double>[lvs.num_clusters];
	std::vector<int> *mc_faces = new std::vector<int>[lvs.num_clusters];
	lvs.marchingCubes(mc_verts, mc_faces);
	for (int i = 0; i < lvs.num_clusters; i++)
	{
		printf("Write cluster %d to .OBJ\n", i);
		TriMesh::writeObj(std::string("../model/armadillo_cluster")+std::to_string(i)+".obj", &mc_verts[i][0], mc_verts[i].size() / 3, &mc_faces[i][0], mc_faces[i].size() / 3);
	}

	return 0;
}