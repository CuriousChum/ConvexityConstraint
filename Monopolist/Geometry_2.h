//
//  Geometry_2.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 13/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_Geometry_2_h
#define CGalTest_Geometry_2_h

#include "Constraint.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>


namespace Geometry_2 {
    using namespace Constraint;
    
    const int Dimension = 2;
    struct PointType : public std::array<ScalarType,Dimension> {
        PointType(std::initializer_list<ScalarType> c) {std::copy(c.begin(), c.end(), this->begin());};
        PointType(){};
        friend std::ostream & operator << (std::ostream & os, const PointType & p){
            print_range(os,p.begin(),p.end()); return os;}
    };
    
    
    namespace CGAL_Traits {
        //        typedef IndexType InfoType;
        
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
        typedef CGAL::Regular_triangulation_filtered_traits_2<K> RT_Traits;
        typedef CGAL::Regular_triangulation_vertex_base_2<RT_Traits> Vbase;
        typedef CGAL::Triangulation_vertex_base_with_info_2
        <InfoType, RT_Traits, Vbase> Vb;
        typedef CGAL::Regular_triangulation_face_base_2<RT_Traits> Cb;
        typedef CGAL::Triangulation_data_structure_2<Vb,Cb> Tds;
        typedef CGAL::Regular_triangulation_2<RT_Traits, Tds> RT;
        
        typedef RT::Vertex_handle Vertex_handle;
        typedef RT::Weighted_point Weighted_point;
        typedef std::pair<Weighted_point, InfoType> Full_point;
        typedef RT::Point Point; //typename CGAL::Point_2<K>
        typedef CGAL::Vector_2<K> Vector;
        
        typedef RT::Finite_faces_iterator Finite_faces_iterator;
        typedef RT::Finite_vertices_iterator Finite_vertices_iterator;
        typedef RT::Vertex_circulator Vertex_circulator;
        typedef RT::Edge_circulator Edge_circulator;
        typedef RT::Edge Edge; // Face handle, int
        
        typedef RT::Face_handle Face_handle;
        
        
        struct EdgeRefinementPriority;
        Edge EdgeToRefine(Face_handle);
        void RefineFaces(RT & rt, const std::vector<Face_handle> &);
        
        ScalarType Parabola(const Point & p){return p[0]*p[0]+p[1]*p[1];}
    }
    
    
    std::ostream & operator << (std::ostream &, const CGAL_Traits::RT &);
    std::ostream & operator << (std::ostream &, const CGAL_Traits::Weighted_point &);

    namespace CGT=CGAL_Traits;
    
#include "Geometry_2.hxx"
    
}


#endif
