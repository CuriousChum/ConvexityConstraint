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
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
//#include <CGAL/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>


namespace Geometry_3 {
    using namespace Constraint;
    
    const int Dimension=3;
    
    namespace CGT { // CGAL_Traits
        struct InfoType {
            IndexType index=BadIndex;
            IndexType boundary=0;
            InfoType(IndexType index_, IndexType boundary_):index(index_), boundary(boundary_){};
            InfoType(){};
            bool OnBoundary() const {return boundary;}
            friend std::ostream & operator << (std::ostream & os, const InfoType & p){
                return os << "{" << p.index << "," << p.boundary << "}";}
        };
        
        // CGAL types
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

	// Modified 18/3/2022
//	https://doc.cgal.org/4.7/Triangulation_3/classCGAL_1_1Regular__triangulation__filtered__traits__3.html
//	        typedef CGAL::Regular_triangulation_filtered_traits_3<K> RT_Traits;
		typedef CGAL::Regular_triangulation_euclidean_traits_3<K> RT_Traits;
	
        typedef CGAL::Triangulation_vertex_base_with_info_3
        <InfoType, RT_Traits> Vb;
        typedef CGAL::Regular_triangulation_cell_base_3<RT_Traits> Cb;
        typedef CGAL::Triangulation_data_structure_3<Vb,Cb> Tds;
        typedef CGAL::Regular_triangulation_3<RT_Traits, Tds> RT;
        
        typedef RT::Weighted_point Weighted_point;
        typedef RT::Point Point; //typename CGAL::Point_2<K>
        typedef CGAL::Vector_3<K> Vector;
        

        typedef RT::Vertex_handle Vertex_handle;
        typedef RT::Cell_handle Cell_handle;
        
        typedef RT::Cell_circulator Cell_circulator; // cells incident to an edge.
        typedef RT::Finite_cells_iterator Finite_cells_iterator;
        typedef RT::Finite_edges_iterator Finite_edges_iterator;
        typedef RT::Finite_vertices_iterator Finite_vertices_iterator;
        
        typedef std::pair<Weighted_point, InfoType> Full_point;
        
        
        ScalarType Parabola(const Point &);
        ScalarType Height(const Weighted_point &);
        
        IndexType OnSameBoundary(Cell_handle);
        
        // Computes a basic indicator of the degeneracy of a simplex, (small is bad)
        // and returns volume and gradients of basis functions by reference.
        ScalarType Degeneracy(Cell_handle, ScalarType &, std::array<Vector,Dimension+1> &);
    }
    
    std::ostream & operator << (std::ostream &, const CGT::Weighted_point &);
    std::ostream & operator << (std::ostream &, const CGT::RT &);
    
#include "Geometry_3.hxx"

}
