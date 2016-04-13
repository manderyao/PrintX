#include "TriMesh.h"
#include "Levelset.h"
#include "opgl.h"

using namespace PRINTX_LIB;

TriMesh mesh, mesh_in;

SphericalCamera *camera = 0;
float focus[3], up[3] = { 0, 1, 0 };
float cameraRadius = 100;

int g_iLeftMouseButton = 0, g_iMiddleMouseButton = 0, g_iRightMouseButton = 0;
int g_vMousePos[2] = { 0, 0 };
int shiftPressed = 0;
int altPressed = 0;
int ctrlPressed = 0;

int windowWidth = 1024;
int windowHeight = 768;
int windowID;



void mouse(int button, int state, int x, int y)
{
	switch (button)
	{
	case GLUT_LEFT_BUTTON:
		g_iLeftMouseButton = (state == GLUT_DOWN);
		shiftPressed = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);
		altPressed = (glutGetModifiers() == GLUT_ACTIVE_ALT);
		ctrlPressed = (glutGetModifiers() == GLUT_ACTIVE_CTRL);
		break;

	case GLUT_MIDDLE_BUTTON:
		g_iMiddleMouseButton = (state == GLUT_DOWN);
		shiftPressed = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);
		ctrlPressed = (glutGetModifiers() == GLUT_ACTIVE_CTRL);
		break;

	case GLUT_RIGHT_BUTTON:
		g_iRightMouseButton = (state == GLUT_DOWN);
		shiftPressed = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);
		ctrlPressed = (glutGetModifiers() == GLUT_ACTIVE_CTRL);
		break;
	}

	g_vMousePos[0] = x;
	g_vMousePos[1] = y;
}


void motion(int x, int y)
{
	int mouseDeltaX = x - g_vMousePos[0];
	int mouseDeltaY = y - g_vMousePos[1];

	g_vMousePos[0] = x;
	g_vMousePos[1] = y;

	if (g_iLeftMouseButton) // handles camera rotations
	{
		camera->MoveRight(-mouseDeltaX);
		camera->MoveUp(mouseDeltaY);
		//std::cout << camera->Phi << "           " << camera->Theta << "\n";
	}

	if (g_iRightMouseButton) // handle zoom in/out
	{
		camera->ZoomIn(cameraRadius * mouseDeltaY);
	}

	if ((g_iMiddleMouseButton)) // handle focus translation
	{
		camera->MoveFocusRight(-mouseDeltaX);
		camera->MoveFocusUp(mouseDeltaY);
	}

	glutSetWindow(windowID);
	glutPostRedisplay();
}


void reshape(int w, int h)
{
	windowWidth = w;
	windowHeight = h;

	glutSetWindow(windowID);
	glViewport(0, 0, w, h);
	//glutPositionWindow(0, 0);
	glutReshapeWindow(w, h);

	glEnable(GL_LIGHT0);
	GLfloat aGa[] = { 0.1, 0.1, 0.1, 1.0 };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, aGa);
	glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

	GLfloat lKa[] = { 0.65, 0.65, 0.65, 1.0 };
	GLfloat lKd[] = { 0.35, 0.35, 0.35, 1.0 };
	GLfloat lKs[] = { 0.1, 0.1, 0.1, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, lKa);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lKd);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lKs);
	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	glutPostRedisplay();
}

void idle()
{
	//glutPostRedisplay();
}

void keyboard(unsigned char Key, int x, int y)
{
}


void display()
{
	glEnable(GL_DEPTH_TEST);
	glClearColor(1, 1, 1, 0);
	glClearDepth(1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);


	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	camera->setOrtho(-cameraRadius*windowWidth / (float)windowHeight, cameraRadius*windowWidth / (float)windowHeight, -cameraRadius, cameraRadius, 0.1f, 1000);
	camera->project();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	camera->Look();

	float green[3] = { 0, 1, 0 };
	float red[3] = { 1, 0.5, 0 };

	mesh_in.renderFlat(red);
	mesh.renderFlatTrans(green, 0.4);
	

	glutSwapBuffers();
}


int main(int argc, char **argv)
{
	mesh.loadOff("../model/armadillo.off");
	double b_min[3], b_max[3];
	mesh.buildBbox(b_min, b_max);

	Levelset lvs(mesh, 1.5, true);
	printf("Level set initialized\n");
	lvs.redist(4.0);
	lvs.advectAllW(2.0);
	for (int i = 0; i < 30; i++)
	{
		lvs.smooth(0.3);
		printf("Level set smooth Iter: %d\n", i);
	}
	//lvs.mergeCntCmpnt();	// use this to remove island and keep only one piece

	std::vector<double> mc_verts;
	std::vector<int> mc_faces;
	lvs.marchingCubes(mc_verts, mc_faces);
	mesh_in.loadMesh(mc_verts, mc_faces);
	printf("Write offset surface to armadillo_in.obj\n");
	mesh_in.writeObj("../model/armadillo_in.obj");

	opgl_glut_init(argc, argv, "MeshSeg", windowWidth, windowHeight, &windowID, display, reshape, mouse, motion, keyboard, idle);
	focus[0] = (b_min[0] + b_max[0]) / 2.0;
	focus[1] = (b_min[1] + b_max[1]) / 2.0;
	focus[2] = (b_min[2] + b_max[2]) / 2.0;
	camera = new SphericalCamera(cameraRadius, 0.0f, 0.0f, focus, up, 0.02);
	glutMainLoop();

	return 0;
}