//
//  Geometry_2.hxx
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 13/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_Geometry_2_hxx
#define CGalTest_Geometry_2_hxx

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
        
        // Eliminate flat triangles
        if(fIt->vertex(0)->info().boundary &
           fIt->vertex(1)->info().boundary &
           fIt->vertex(2)->info().boundary) continue;
        
        std::array<Vector,3> g; // gradients of basis functions
        std::array<ScalarType,3> h; // heights
        for(int i=0; i<3; ++i){
            Vertex_handle
            p = fIt->vertex(rt.cw(i)),
            q = fIt->vertex(rt.ccw(i)),
            r = fIt->vertex(i);
            g[i]=p->point() - q->point();
            
            const ScalarType x = r->point().x(), y=r->point().y();
            h[i] = x*x+y*y - r->point().weight();
        }
        const ScalarType area2 = CGAL::determinant(g[0],g[1]);
        Vector gradient(0,0);
        for(int i=0; i<3; ++i){
            g[i] = g[i].perpendicular(CGAL::POSITIVE) / area2;
            gradient = gradient + h[i]*g[i];
        }
        
        
        os << "{{" <<
        fIt->vertex(0)->info().index << "," <<
        fIt->vertex(1)->info().index << "," <<
        fIt->vertex(2)->info().index << "},{" <<
        fIt->vertex(0)->info().boundary << "," <<
        fIt->vertex(1)->info().boundary << "," <<
        fIt->vertex(2)->info().boundary << "},{" <<
        gradient[0] << "," << gradient[1] << "}},";
        
    }
    
    os << "}}";
    return os;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%% Refinement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

namespace CGAL_Traits {
    
    struct EdgeRefinementPriority {
        bool operator() (Edge a, Edge b) const {
            const Face_handle af=a.first, bf=b.first;
            const int ai=a.second, bi=b.second;
            const ScalarType
            ad = CGAL::squared_distance(af->vertex(RT::cw(ai))->point(), af->vertex(RT::ccw(ai))->point()),
            bd = CGAL::squared_distance(bf->vertex(RT::cw(bi))->point(), bf->vertex(RT::ccw(bi))->point());
            
            typedef std::tuple<ScalarType,Face_handle,int> TupleType;
            return TupleType(ad,af,ai) > TupleType(bd,bf,bi); // Refine longuest edges first
        }
    };
    
    Edge EdgeToRefine(Face_handle f) {
        Edge e0(f,0);
        for(int i=1; i<3; ++i){
            Edge e1(f,i);
            if( !EdgeRefinementPriority()(e0,e1) ) e0=e1;
        }
        return e0;
    }
    
    void RefineFaces(RT & rt, const std::vector<Face_handle> & faces){
        size_t counter = rt.number_of_vertices();
        
        std::set<Edge,EdgeRefinementPriority> queue;
        for(Face_handle f : faces)
            queue.insert(EdgeToRefine(f));
        
        while(!queue.empty()){
            Edge e0 = *queue.begin();
            Edge e1 = rt.mirror_edge(e0);
            Face_handle f1 = e1.first; int i1=e1.second;
            if(!rt.is_infinite(f1) && e1 != EdgeToRefine(f1)){
                queue.insert( EdgeToRefine(f1) );
                continue;
            }
            queue.erase(e0);
            queue.erase(e1);
            const Vertex_handle
            vp = f1->vertex(RT::cw(i1)),
            vq = f1->vertex(RT::ccw(i1));
            
            const Weighted_point &
            p = vp->point(),
            q = vq->point();
            
            Vertex_handle newVertex =
            rt.insert_in_edge(Weighted_point( CGAL::midpoint(p.point(),q.point()), (p.weight()+q.weight())/2. ),
                              f1,i1);
            newVertex->info().index = (int)counter++;
            newVertex->info().boundary = vp->info().boundary & vq->info().boundary;
        }
        
    }
    
}


#endif
