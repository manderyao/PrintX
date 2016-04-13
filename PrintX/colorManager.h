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
* color manager
*/


#ifndef __PRINTX_COLOR_MANAGER_H__
#define __PRINTX_COLOR_MANAGER_H__

namespace PRINTX_LIB
{
	static float renderColor[60] = { 0xce / 255.0, 0x66 / 255.0, 0x66 / 255.0,
		0x66 / 255.0, 0xab / 255.0, 0xce / 255.0,
		0xce / 255.0, 0x9a / 255.0, 0x66 / 255.0,
		0x66 / 255.0, 0xce / 255.0, 0x66 / 255.0,
		0x66 / 255.0, 0x78 / 255.0, 0xce / 255.0,
		0xce / 255.0, 0x66 / 255.0, 0x89 / 255.0,
		0xcd / 255.0, 0xff / 255.0, 0x66 / 255.0,
		0x66 / 255.0, 0xce / 255.0, 0x89 / 255.0,
		0x66 / 255.0, 0xff / 255.0, 0x66 / 255.0,
		0x66 / 255.0, 0xcd / 255.0, 0xff / 255.0,

		0xce / 255.0, 0x66 / 255.0, 0x66 / 255.0,
		0x66 / 255.0, 0xab / 255.0, 0xce / 255.0,
		0xce / 255.0, 0x9a / 255.0, 0x66 / 255.0,
		0x66 / 255.0, 0xce / 255.0, 0x66 / 255.0,
		0x66 / 255.0, 0x78 / 255.0, 0xce / 255.0,
		0xce / 255.0, 0x66 / 255.0, 0x89 / 255.0,
		0xcd / 255.0, 0xff / 255.0, 0x66 / 255.0,
		0x66 / 255.0, 0xce / 255.0, 0x89 / 255.0,
		0x66 / 255.0, 0xff / 255.0, 0x66 / 255.0,
		0x66 / 255.0, 0xcd / 255.0, 0xff / 255.0
	};

	template <class T>
	static void getHeatColor(T value, float *red, float *green, float *blue)
	{
#if 1
		const int NUM_COLORS = 4;
		static float color[NUM_COLORS][3] = { { 0, 0, 1 }, { 0, 1, 0 }, { 1, 1, 0 }, { 1, 0, 0 } };
		// A static array of 4 colors:  (blue,   green,  yellow,  red) using {r,g,b} for each.
#else
		const int NUM_COLORS = 3;
		static float color[NUM_COLORS][3] = { { 0, 0, 1 }, { 1, 1, 0 }, { 1, 0, 0 } };
		// A static array of 4 colors:  (blue,   green,  yellow,  red) using {r,g,b} for each.
#endif
		int idx1;        // |-- Our desired color will be between these two indexes in "color".
		int idx2;        // |
		float fractBetween = 0;  // Fraction between "idx1" and "idx2" where our value is.

		if (value <= 0)      { idx1 = idx2 = 0; }    // accounts for an input <=0
		else if (value >= 1)  { idx1 = idx2 = NUM_COLORS - 1; }    // accounts for an input >=0
		else
		{
			value = value * (NUM_COLORS - 1);        // Will multiply value by 3.
			idx1 = floor(value);                  // Our desired color will be after this index.
			idx2 = idx1 + 1;                        // ... and before this index (inclusive).
			fractBetween = value - float(idx1);    // Distance between the two indexes (0-1).
		}

		*red = (color[idx2][0] - color[idx1][0])*fractBetween + color[idx1][0];
		*green = (color[idx2][1] - color[idx1][1])*fractBetween + color[idx1][1];
		*blue = (color[idx2][2] - color[idx1][2])*fractBetween + color[idx1][2];
	}
}
#endif