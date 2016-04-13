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



#ifndef __PRINTX_MATERIAL_H__
#define __PRINTX_MATERIAL_H__

#include "opgl_header.h"
namespace PRINTX_LIB
{
class Material
{
public:
	Material(float *Ka_, float *Kd_, float *Ks_, float alpha_, float shiness_, int use_specular_)
		: alpha(alpha_), shiness(shiness_), use_specular(use_specular_)
	{
		memcpy(Ka, Ka_, sizeof(float) * 3);
		memcpy(Kd, Kd_, sizeof(float) * 3);
		memcpy(Ks, Ks_, sizeof(float) * 3);
	}

	inline void materialize(GLenum face)
	{
		glMaterialfv(face, GL_AMBIENT, Ka);
		glMaterialfv(face, GL_DIFFUSE, Kd);
		glMaterialfv(face, GL_SPECULAR, Ks);
		glMaterialf(face, GL_SHININESS, shiness);
	}

public:
		float Ka[3], Kd[3], Ks[3];
		float alpha;
		float shiness;
		int use_specular;
};
}
#endif