//
//  Regular_Test.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 28/01/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_Regular_Test_h
#define CGalTest_Regular_Test_h


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_filtered_traits_2<K>  Traits;
typedef CGAL::Regular_triangulation_2<Traits> Regular_triangulation;

typedef Regular_triangulation::Weighted_point Weighted_Point;
typedef Regular_triangulation::Finite_vertices_iterator Finite_vertices_iterator;
typedef Regular_triangulation::Vertex_circulator Vertex_circulator;
typedef Regular_triangulation::Edge_circulator Edge_circulator;

void CGalTest_Manual()
{
    std::ifstream in("data/regular.cin");
    Regular_triangulation::Weighted_point wp;
    int count = 0;
    std::vector<Regular_triangulation::Weighted_point> wpoints;
    while(in >> wp){
        count++;
        wpoints.push_back(wp);
    }
    Regular_triangulation rt(wpoints.begin(), wpoints.end());
    rt.is_valid();
    std::cout << "number of inserted points : " << count << std::endl;
    std::cout << "number of vertices :  " ;
    std::cout << rt.number_of_vertices() << std::endl;
    std::cout << "number of hidden vertices :  " ;
    std::cout << rt.number_of_hidden_vertices() << std::endl;
}

void CGal_MyTest(){
    
    std::vector<Weighted_Point> wpoints;
    
    int n=5;
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            wpoints.push_back(Weighted_Point(i,j,0));
    Regular_triangulation rt(wpoints.begin(), wpoints.end());
    rt.is_valid();

    std::cout << rt.number_of_vertices() << std::endl;
    
    // Iterate over vertices.
    for(Finite_vertices_iterator vIt = rt.finite_vertices_begin(), vEnd = rt.finite_vertices_end();
        vIt!=vEnd; ++vIt){
        std::cout << "Vertex : " << vIt->point() << std::endl;
        
        // Iterate over neighbor vertices
        Vertex_circulator vc = rt.incident_vertices(vIt), vDone(vc);
        do {
            std::cout << "Neighbor : " << vc->point() << std::endl;
        } while(++vc != vDone);
        
        Edge_circulator eC = rt.incident_edges(vIt), eDone(eC);
        do {
            if(rt.is_infinite(eC)){
                std::cout << "Infinite edge\n";
                break;
            }
            std::cout << "Finite edge, source point : " <<
            eC->first->vertex( rt.cw( eC->second ) )->point() <<
            " opposite point : " <<
            eC->first->vertex( rt.ccw( eC->second ) )->point() <<
            "\n";
            
            // dual segment
            CGAL::Object o = rt.dual(eC);
            auto dualSeg = CGAL::object_cast<K::Segment_2>(&o);
            if(!dualSeg) std::cout << "Infinite dual ray\n";
            else std::cout << "dual segment : source " << dualSeg->source() << ", target " << dualSeg->target() << "\n";
//            std::cout << "dual segment : source " << dualSeg.source() << "\n";
            
            //std::cout << "Edge endpoint indices : " << eC->first-> << std::endl;
            
            //<< eC.first << "," << eC.second << "\n";
        } while(++eC!=eDone);
        
        // Compute the region areas, gradient and hessian wrt weights.
        
        // Get the dual edge lengths for gradients
        
        
        
    }
    
}


#endif
