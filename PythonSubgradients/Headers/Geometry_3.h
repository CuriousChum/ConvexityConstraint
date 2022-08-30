#pragma once
//
//  Geometry_3.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 13/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#include "Constraint.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
//#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
//#include <CGAL/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>


namespace Geometry_3 {
    using namespace Constraint;
    
    const int Dimension=3;
    
    namespace CGT {
        // CGAL types
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
		typedef CGAL::Triangulation_vertex_base_with_info_3<InfoType, K> Vb0;
		typedef CGAL::Regular_triangulation_vertex_base_3<K,Vb0> Vb;
		typedef CGAL::Regular_triangulation_cell_base_3<K> Fb;
		typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
		typedef CGAL::Regular_triangulation_3<K, Tds> RT;

        typedef RT::Weighted_point Weighted_point;
        typedef RT::Bare_point Point;
        typedef CGAL::Vector_3<K> Vector;
        typedef RT::Vertex_handle Vertex_handle;
        typedef RT::Cell_handle Cell_handle;
        
        typedef RT::Cell_circulator Cell_circulator; // cells incident to an edge.
        typedef RT::Finite_cells_iterator Finite_cells_iterator;
        typedef RT::Finite_edges_iterator Finite_edges_iterator;
        typedef RT::Finite_vertices_iterator Finite_vertices_iterator;
        
        typedef std::pair<Weighted_point, InfoType> Full_point;
        
        
        ScalarType Parabola(const Point &);        
//        IndexType OnSameBoundary(Cell_handle);
        
        // Computes a basic indicator of the degeneracy of a simplex, (small is bad)
        // and returns volume and gradients of basis functions by reference.
        ScalarType Degeneracy(Cell_handle, ScalarType &, std::array<Vector,Dimension+1> &);
    }
    
    std::ostream & operator << (std::ostream &, const CGT::Weighted_point &);
    std::ostream & operator << (std::ostream &, const CGT::RT &);
    
#include "Geometry_3.hxx"

}
