#pragma once
//
//  Geometry_2.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 13/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#include "Constraint.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
//#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>


namespace Geometry_2 {
using namespace Constraint;

const int Dimension = 2;

namespace CGAL_Traits {
// CGAL types
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<InfoType, K> Vb0;
typedef CGAL::Regular_triangulation_vertex_base_2<K,Vb0> Vb;
typedef CGAL::Regular_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Regular_triangulation_2<K, Tds> RT;

typedef RT::Vertex_handle Vertex_handle;
typedef RT::Weighted_point Weighted_point;
typedef std::pair<Weighted_point, InfoType> Full_point;
typedef RT::Bare_point Point; 
typedef CGAL::Vector_2<K> Vector;

typedef RT::Finite_faces_iterator Finite_faces_iterator;
typedef RT::Finite_vertices_iterator Finite_vertices_iterator;
typedef RT::Vertex_circulator Vertex_circulator;
typedef RT::Edge_circulator Edge_circulator;
typedef RT::Edge Edge; // Face handle, int
typedef RT::Face_handle Face_handle;

ScalarType Parabola(const Point & p){return p[0]*p[0]+p[1]*p[1];}

}


 std::ostream & operator << (std::ostream &, const CGAL_Traits::RT &);
 std::ostream & operator << (std::ostream &, const CGAL_Traits::Weighted_point &);
 
 namespace CGT=CGAL_Traits;
 
 #include "Geometry_2.hxx"
 
}
