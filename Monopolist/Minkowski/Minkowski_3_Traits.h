//
//  Minkowski_3_Traits.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 18/08/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_Minkowski_3_Traits_h
#define CGalTest_Minkowski_3_Traits_h

#include "Constraint.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_3.h>
//#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>


namespace Minkowski_3 {
    using namespace Constraint;
    const int Dimension=3;
    namespace CGT {
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef K::Point_3 Point;
        typedef K::Vector_3 Vector;
        
        // A vertex type with an index
        template <class Refs>
        struct My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> {
            typedef CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> superclass;
            IndexType index=BadIndex;
            My_vertex(const Point & p):superclass(p){};
            My_vertex(){};
        };
        
        struct My_items : public CGAL::Polyhedron_items_3 {
            template <class Refs, class Traits>
            struct Vertex_wrapper {
                typedef My_vertex<Refs> Vertex;
            };
            
        };
        
        typedef CGAL::Polyhedron_3<K, My_items>  Polyhedron;
    }
    
}


#endif
