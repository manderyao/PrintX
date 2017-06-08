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
//	mesh operation using CGAL (boolean, clean, ...)

#ifndef __PRINTX_Mesh_OP_H__
#define __PRINTX_Mesh_OP_H__

#include <iostream> 
#include <fstream>
#include <CGAL/Gmpq.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/print_wavefront.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include "ManifoldFix.h"

namespace PRINTX_LIB
{
	static enum{ BOOL_INTERSECT = 0, BOOL_UNION = 1, BOOL_DIFF = 2 };

	// A modifier creating a triangle with the incremental builder.
	template<class HDS>
	class PolyhedronBuilder : public CGAL::Modifier_base < HDS >
	{
	public:
		double *verts;	// shallow copy, careful when the original data changes
		int n_verts;
		int *tris;
		int n_tris;
		PolyhedronBuilder(double *verts_, int n_verts_, int *tris_, int n_tris_) : verts(verts_), n_verts(n_verts_), tris(tris_), n_tris(n_tris_) {}
		void operator()(HDS& hds)
		{
			typedef typename HDS::Vertex   Vertex;
			typedef typename Vertex::Point Point;
			// create a cgal incremental builder
			CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
			B.begin_surface(n_verts, n_tris);
			// add the polyhedron vertices
			for (int i = 0; i < n_verts; i += 1)
			{
				B.add_vertex(Point(verts[3 * i], verts[3 * i + 1], verts[3 * i + 2]));
			}
			// add the polyhedron triangles
			for (int i = 0; i < n_tris; i += 1)
			{
				B.begin_facet();
				B.add_vertex_to_facet(tris[3 * i]);
				B.add_vertex_to_facet(tris[3 * i + 1]);
				B.add_vertex_to_facet(tris[3 * i + 2]);
				B.end_facet();
			}
			// finish up the surface
			B.end_surface();
		}
	};


	// build NefPoly from vertices and faces
	static bool loadNefPoly(Nef_polyhedron &nef_poly, double *verts, int n_verts, int *tris, int n_tris, bool regular, bool normal, bool useDefault)
	{
		Polyhedron poly;
		PolyhedronBuilder<Polyhedron::HalfedgeDS> builder_in(verts, n_verts, tris, n_tris);
		poly.delegate(builder_in);
		
		//nef_poly = Nef_polyhedron(poly);
		return polyhe2nef(poly, nef_poly, regular, normal, useDefault);	 // safe
	}

	// export NefPoly to vertices and faces
	static bool exportNefPoly(Nef_polyhedron &nef_poly, std::vector<double> &verts, std::vector<int> &tris)
	{
		Polyhedron poly;
		if (nef2polyhe(nef_poly, poly))// safe
		{
			verts.clear();
			tris.clear();
			for (Polyhedron::Vertex_iterator v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
			{
				verts.push_back(CGAL::to_double(v->point().x()));
				verts.push_back(CGAL::to_double(v->point().y()));
				verts.push_back(CGAL::to_double(v->point().z()));
			}
			for (Polyhedron::Facet_iterator i = poly.facets_begin(); i != poly.facets_end(); ++i)
			{
				Polyhedron::Halfedge_around_facet_circulator j = i->facet_begin();
				do
				{
					tris.push_back(std::distance(poly.vertices_begin(), j->vertex()));
				} while (++j != i->facet_begin());
			}
			return true;
		}
		else
			return false;
	}

	// write to obj directly
	static bool writeNefPoly2OBJ(Nef_polyhedron &nef_poly, std::string &filename)
	{
		Polyhedron poly;
		if (nef2polyhe(nef_poly, poly))// safe
		{
			std::ofstream ofs(filename);
			CGAL::print_polyhedron_wavefront(ofs, poly);
			return true;
		}
		else
			return false;
	}

	// clean Nef_polyhedron, return  if successful
	static bool cleanNefPoly(Nef_polyhedron &nef_poly)
	{
		if (!nef_poly.is_simple())
		{
			return make2Manifold(nef_poly);
		}
		else
			return true;
	}


	// boolean on two Nef_polyhedron
	static void boolOperator(Nef_polyhedron &nef_mesh1, Nef_polyhedron &nef_mesh2, Nef_polyhedron &nef_result, int op)
	{
		switch (op)
		{
		case BOOL_INTERSECT:
			nef_result = nef_mesh1*nef_mesh2;
			break;
		case BOOL_UNION:
			nef_result = nef_mesh1 + nef_mesh2;
			break;
		case BOOL_DIFF:
			nef_result = nef_mesh1 - nef_mesh2;
			break;
		}
	}
}

#endif