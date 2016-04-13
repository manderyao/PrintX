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
//  traditional opengl
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef __PRINTX_OPGL_H__
#define __PRINTX_OPGL_H__

#include "camera.h"
#include "material.h"

namespace PRINTX_LIB
{
	static bool opgl_glew_init()
	{
		GLenum err = glewInit();
		if (GLEW_OK != err)
		{
			printf("GLEW init failed: %s!\n", glewGetErrorString(err));
			return false;
		}
		else
		{
			//printf("GLEW init success!\n");
			return true;
		}
	}

	static void opgl_glut_init(int argc, char* argv[], char * windowTitle, int windowWidth, int windowHeight, int * windowID,
		void(*displayFunc)(), void(*reshape)(int, int), void(*mouse)(int button, int state, int x, int y),
		void(*motion)(int, int), void(*keyboard)(unsigned char key, int x, int y), void(*idle)(void))
	{
		// Initialize GLUT.
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE);
		glutInitWindowSize(windowWidth, windowHeight);
		*windowID = glutCreateWindow(windowTitle);

		// Setup GLUT callbacks.
		glutDisplayFunc(displayFunc);
		glutReshapeFunc(reshape);
		glutMouseFunc(mouse);
		glutMotionFunc(motion);;
		glutKeyboardFunc(keyboard);
		glutIdleFunc(idle);
		//GLUI_Master.set_glutIdleFunc (idle);
		//GLUI_Master.set_glutIdleFunc(idleFunction);
		//GLUI_Master.set_glutKeyboardFunc(keyboardFunction);
		//GLUI_Master.set_glutReshapeFunc(reshape);
		//GLUI_Master.set_glutMouseFunc(mouseButtonActivityFunction);
	}
}
#endif