//
//  Refinement_2_Traits.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 21/08/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_Refinement_2_Traits_h
#define CGalTest_Refinement_2_Traits_h

#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "Geometry_2.h"

namespace Refinement_2_Test {
    using namespace Geometry_2::CGAL_Traits;
    
    
    void ShowFiniteFaces(const RT & rt, std::ostream & os=std::cout){
        std::cout << "Showing finite faces\n";
        for(auto fit=rt.finite_faces_begin(); fit!=rt.finite_faces_end(); ++fit){
            os << "Face\n";
            for(int i=0; i<3; ++i)
                os << fit->vertex(i)->point().point() << "\n";
        }

    }
    

    void Test(){
        typedef std::pair<Weighted_point,InfoType> FullPoint;
        std::vector<FullPoint> pts =
        {
            { Weighted_point(Point(0,0),0),InfoType(0,1) },
            { Weighted_point(Point(1,0),0),InfoType(1,1) },
            { Weighted_point(Point(1,1),0),InfoType(2,1) },
            { Weighted_point(Point(0,1),0),InfoType(3,1) },
        };
        
        RT rt = RT(pts.begin(),pts.end());
        rt.infinite_vertex()->info() = {Constraint::InfiniteIndex, -1};
        int error = (int)rt.number_of_hidden_vertices();
        
        auto fit = rt.finite_faces_begin();
        for(int i=0; i<3; ++i)
            std::cout << fit->vertex(i)->point() << "\n";
        rt.flip(fit,0);
        std::cout << "flipping\n";
        for(int i=0; i<3; ++i)
            std::cout << fit->vertex(i)->point() << "\n";
        std::cout << "Inserting\n";
        auto vit = rt.insert_in_edge(  Weighted_point(Point(1./2,1./2),0), fit,1 );
        vit->info().index=4;
        
        ShowFiniteFaces(rt);
        
        RT::Face_handle h0, h1;
        {
            auto ftest0 = rt.finite_faces_begin(); ++ftest0;
            h0=ftest0;
        }
        
        fit = rt.finite_faces_begin();
        vit=rt.insert_in_edge( Weighted_point(Point(1.,1./2),0),fit,0);
        vit->info().index=5;
        
        ShowFiniteFaces(rt);
        
        h1 = rt.finite_faces_begin();
        std::cout << "Handles preserved : " << (h0==h1) << "\n";
    }
    
    std::vector<Face_handle>
    FaceHandlesFromIndices(const RT & rt, std::vector<int> indices){
        std::sort(indices.begin(),indices.end());
        std::vector<Face_handle> handles;
        auto fit = rt.finite_faces_begin();
        auto iit=indices.begin();
        int counter=0;
        for(;fit!=rt.finite_faces_end() && iit!=indices.end();
            ++fit, ++counter){
            if(counter==*iit){
                handles.push_back(fit);
                ++iit;
            }
        }
        return handles;
    }
    
    
    
    void Test2(){
        typedef std::pair<Weighted_point,InfoType> FullPoint;
        std::vector<FullPoint> pts =
        {
            { Weighted_point(Point(0,0),0),InfoType(0,1) },
            { Weighted_point(Point(1,0),0),InfoType(1,1) },
            { Weighted_point(Point(1,1),0),InfoType(2,1) },
            { Weighted_point(Point(0,1),0),InfoType(3,1) },
        };
        
        RT rt = RT(pts.begin(),pts.end());
        
        ShowFiniteFaces(rt);
        RefineFaces(rt,FaceHandlesFromIndices(rt, {1,2}));
        ShowFiniteFaces(rt);
        RefineFaces(rt,FaceHandlesFromIndices(rt, {3}));
        ShowFiniteFaces(rt);
        RefineFaces(rt,FaceHandlesFromIndices(rt, {4}));
        ShowFiniteFaces(rt);
        
        
    }
    
    
}



/*
 struct Traits {
 struct Point_2 {K::Point_2 position; IndexType index;};
 };
 typedef CGAL::HalfedgeDS_default<Traits> HDS;
 typedef CGAL::HalfedgeDS_decorator<HDS>  Decorator;
 typedef HDS::Halfedge_iterator           Iterator;
 
 
 
 int Test() {
 HDS hds;
 Decorator decorator(hds);
 decorator.create_loop();
 CGAL_assertion( decorator.is_valid());
 int n = 0;
 for ( Iterator i = hds.halfedges_begin(); i != hds.halfedges_end(); ++i )
 ++n;
 std::cout << n << "\n";
 
 auto edge =hds.halfedges_begin();
 decorator.make_hole(edge->opposite());
 CGAL_assertion( decorator.is_valid());
 
 auto e=decorator.create_loop();
 CGAL_assertion( decorator.is_valid());
 
 
 // I do not see how to create and insert new (half)edges !!??!!
 //        decorator.insert_halfedge(e, edge);
 //        CGAL_assertion( decorator.is_valid());
 
 return 0;
 }*/


#endif
