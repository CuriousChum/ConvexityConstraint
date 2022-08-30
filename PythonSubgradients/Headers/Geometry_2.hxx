#pragma once
//
//  Geometry_2.hxx
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 13/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

// %%%%%%%%%%%%%%% Printing %%%%%%%%%%%%%%%%%

std::ostream & operator << (std::ostream & os, const CGAL_Traits::Weighted_point & p){
	os << "{" << p.weight() << "," <<
	"{" << p.point().x() << "," << p.point().y() << "}}";
	return os;
}


std::ostream & operator << (std::ostream & os, const CGAL_Traits::RT & rt){
	using namespace CGAL_Traits;
	
	os << "{" << '"' << "points" << '"' << "-> {";
	
	for(Finite_vertices_iterator vIt = rt.finite_vertices_begin(), vEnd = rt.finite_vertices_end();
        vIt!=vEnd; ++vIt){
        os << "{" <<
        vIt->info() << "," <<
        vIt->point().weight() << "," <<
        "{" << vIt->point().x() << "," << vIt->point().y() << "}" <<
        "},";
    }
    
    os << "}," << '"' << "faces" << '"' << "-> {";
    
    for(auto fIt=rt.finite_faces_begin(), fEnd=rt.finite_faces_end();
        fIt!=fEnd; ++fIt){
        
/*        // Eliminate flat triangles
        if(fIt->vertex(0)->info().boundary &
           fIt->vertex(1)->info().boundary &
           fIt->vertex(2)->info().boundary) continue;
        */
	
        std::array<Vector,3> g; // gradients of basis functions
        std::array<ScalarType,3> h; // heights
        for(int i=0; i<3; ++i){
            Vertex_handle
            p = fIt->vertex(rt.cw(i)),
            q = fIt->vertex(rt.ccw(i)),
            r = fIt->vertex(i);
            g[i]=p->point().point() - q->point().point();
            
            const ScalarType x = r->point().x(), y=r->point().y();
            h[i] = x*x+y*y - r->point().weight();
        }
        const ScalarType area2 = CGAL::determinant(g[0],g[1]);
        Vector gradient(0,0);
        for(int i=0; i<3; ++i){
            g[i] = g[i].perpendicular(CGAL::POSITIVE) / area2;
            gradient = gradient + h[i]*g[i];
        }
	
	os << "{" <<
	fIt->vertex(0)->info() << "," <<
	fIt->vertex(1)->info() << "," <<
	fIt->vertex(2)->info() << "},{" <<
	gradient[0] << "," << gradient[1] << "}},";

    }
    
    os << "}}";
    return os;
}

