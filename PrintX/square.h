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
// square indices look up and marching squares look up
//**************************************************************************************

#ifndef __PRINTX_SQUARE_H__
#define __PRINTX_SQUARE_H__

namespace PRINTX_LIB
{
	namespace SQUARE_STRUCT
	{
		static int vert_to_idx[4][2] =
		{
			{ 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 },
		};

		//last element for direction: 0:+x; 1:+y; 2:+z;
		static int edge_to_vert[4][3] =
		{
			{ 0, 1, 0 }, { 1, 2, 1 }, { 3, 2, 0 }, { 0, 3, 1 }
		};


		// for idx=[0-3], if idx is inside, type += 1<<idx;

		// find num of triangles according to the type above
		static int case_to_numpolys[16] =
		{
			0, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 0
		};

		// find triangles according to the type above
		static int edge_connect_list[16][2][2] =
		{
			{ { -1, -1 }, { -1, -1 } },
			{ { 0, 3 }, { -1, -1 } },
			{ { 1, 0 }, { -1, -1 } },
			{ { 1, 3 }, { -1, -1 } },

			{ { 2, 1 }, { -1, -1 } },
			{ { 2, 3 }, { 0, 1 } },
			{ { 2, 0 }, { -1, -1 } },
			{ { 2, 3 }, { -1, -1 } },

			{ { 3, 2 }, { -1, -1 } },
			{ { 0, 2 }, { -1, -1 } },
			{ { 1, 2 }, { 3, 0 } },
			{ { 1, 2 }, { -1, -1 } },

			{ { 3, 1 }, { -1, -1 } },
			{ { 0, 1 }, { -1, -1 } },
			{ { 3, 0 }, { -1, -1 } },
			{ { -1, -1 }, { -1, -1 } },
		};
	}
}
#endif