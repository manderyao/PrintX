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
//	priority queue for fast marching, values (level sets) are supposed to be all non-negative
//**************************************************************************************


#ifndef __PRINTX_PRIORITY_QUEUE_FM_H__
#define __PRINTX_PRIORITY_QUEUE_FM_H__

#include "stdlib.h"

namespace PRINTX_LIB
{
	template <class TYPE>
	class PriorityQueueFM
	{
	public:
		TYPE *dists;	// reference to the comparing values

		//	a min heap
		int *index;
		int *p;
		int num, num_max;

		PriorityQueueFM(int N)
		{
			num_max = N;
			index = (int*)malloc(sizeof(int*)*num_max);
			p = (int*)malloc(sizeof(int)*num_max);
			num = 0;
		}
		~PriorityQueueFM()
		{
			free(index);
			free(p);
		}
		void setDistsPointer(TYPE *dists_)
		{
			dists = dists_;
			num = 0;
		}
		void push(int idx)
		{
			index[num] = idx;
			p[idx] = num;
			num++;
			shift_up(num - 1);
		}
		bool empty()
		{
			return num <= 0;
		}
		int top()
		{
			if (num > 0)
				return index[0];
			else
				return -1;
		}
		int pop()
		{
			int ret = -1;
			if (num > 0)
			{
				ret = index[0];
				index[0] = index[num - 1];
				p[index[0]] = 0;
				num--;
				shift_down(0);
			}
			return ret;
		}
		void clear()
		{
			num = 0;
		}
		void swap(int i, int j)
		{
			int tmp = index[i];
			index[i] = index[j];
			index[j] = tmp;

			p[index[i]] = i;
			p[index[j]] = j;
		}

		void shift_up(int idx)
		{
			int current = idx, parent;
			while (current > 0)
			{
				parent = (current - 1) / 2;
				if (dists[index[current]] < dists[index[parent]])
				{
					swap(current, parent);
				}
				else
					break;

				current = parent;
			}
		}
		void shift_down(int idx)
		{
			int current = idx, child = 2 * current + 1;

			while (child < num)
			{
				if (child + 1 < num && dists[index[child + 1]] < dists[index[child]])
				{
					child = child + 1;
				}

				if (dists[index[child]] < dists[index[current]])
				{
					swap(current, child);
				}
				current = child;
				child = 2 * current + 1;
			}
		}
		void update_dist(int idx, TYPE dist)
		{
			int pp = p[idx];
			if (pp >= 0 && pp < num && index[pp] == idx)
			{
				if (dist < dists[idx])
				{
					dists[idx] = dist;
					shift_up(pp);
				}
			}
			else
				printf("Error: update a grid that is not in the priority queue!!!!!!!!!!\n");
		}
	};
}

#endif