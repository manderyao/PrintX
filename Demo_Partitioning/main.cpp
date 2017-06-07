#include "TriMesh.h"
#include "LevelsetMulti.h"
#include "meshOp.h"
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
	
	printf("Successfully loaded ../model/armadillo.seg.\n");

	LevelsetMulti lvs(mesh, 1.5, false);	
	printf("Level set initialized\n");
	for (int i = 0; i < 30; i++)
	{	
		lvs.smooth(0.2);
		printf("Level set smooth Iter: %d\n", i);
	}

	printf("Start marching cubes\n");
	std::vector<double> *mc_verts = new std::vector<double>[lvs.num_clusters];
	std::vector<int> *mc_faces = new std::vector<int>[lvs.num_clusters];
	lvs.marchingCubes(mc_verts, mc_faces);

	Nef_polyhedron nefpoly;
	loadNefPoly(nefpoly, mesh.verts, mesh.n_verts, mesh.faces, mesh.n_faces);
	for (int i = 0; i < lvs.num_clusters; i++)
	{
		printf("Build piece %d\n", i);
		Nef_polyhedron nefpoly1, nefpoly2;
		loadNefPoly(nefpoly1, &mc_verts[i][0], mc_verts[i].size() / 3, &mc_faces[i][0], mc_faces[i].size() / 3);		
		boolOperator(nefpoly1, nefpoly, nefpoly2, BOOL_INTERSECT);

		printf("Write piece %d to .OBJ\n", i);
		if (cleanNefPoly(nefpoly2))
			writeNefPoly2OBJ(nefpoly2, std::string("../model/armadillo_piece") + std::to_string(i) + ".obj");
		else
			printf("Fail to clean piece %d!!!\n", i);
	}

	delete[]mc_verts;
	delete[]mc_faces;
	return 0;
}



//int main()
//{
//	TriMesh mesh, mesh2;
//	mesh.createCube(1, 1, 1);
//	mesh2.createCube(1, 1, 1);
//	mesh2.translate(0.3, 0.3, 0.3);
//	Nef_polyhedron nefpoly, nefpoly1, nefpoly2, nefpoly3;
//
//	Polyhedron poly;
//	PolyhedronBuilder<Polyhedron::HalfedgeDS> builder_in(mesh.verts, mesh.n_verts, mesh.faces, mesh.n_faces);
//	poly.delegate(builder_in);
//	nefpoly = Nef_polyhedron(poly);
//	printf("Built test!");
//
//	loadNefPoly(nefpoly1, mesh.verts, mesh.n_verts, mesh.faces, mesh.n_faces);
//	printf("Built first!");
//	loadNefPoly(nefpoly2, mesh2.verts, mesh2.n_verts, mesh2.faces, mesh2.n_faces);
//	boolOperator(nefpoly1, nefpoly2, nefpoly3, BOOL_INTERSECT);
//	writeNefPoly2OBJ(nefpoly3, std::string("../model/cube.obj"));
//}