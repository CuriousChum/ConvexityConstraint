//
//  Geometry_3.hxx
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 13/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_Geometry_3_hxx
#define CGalTest_Geometry_3_hxx

ScalarType CGT::Parabola(const Point & p){
    return p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
}


ScalarType CGT::Height(const Weighted_point & p){
    return Parabola(p.point()) - p.weight();
}

IndexType CGT::OnSameBoundary(Cell_handle fh){
    IndexType flag = -1;
    for(int i=0; i<4; ++i)
        flag = flag & fh->vertex(i)->info().boundary;
    return flag;
}

ScalarType CGT::Degeneracy(Cell_handle fh, ScalarType & volume, std::array<Vector,Dimension+1> & g){
    std::array<Point,Dimension+1> p;
    for(int i=0; i<4; ++i)
        p[i] = fh->vertex(i)->point();
    
    const Vector
    e01 = p[1]-p[0],
    e02 = p[2]-p[0],
    e03 = p[3]-p[0];
    
    const ScalarType volume6 = CGAL::determinant(e01,e02,e03);
    
    volume = volume6/6;
    
    const ScalarType sqlSum = e01.squared_length()+e02.squared_length()+e03.squared_length();
    const ScalarType degeneracy = volume6 / (sqlSum*sqrt(sqlSum));
    
    g[1] = CGAL::cross_product(e02,e03)/volume6;
    g[2] = CGAL::cross_product(e03,e01)/volume6;
    g[3] = CGAL::cross_product(e01,e02)/volume6;
    
    g[0] = -(g[1]+g[2]+g[3]);
    
    return degeneracy;
}

// %%%%%%%%%%%%%%% Export %%%%%%%%%%%%%%


std::ostream & operator << (std::ostream & os , const CGT::Weighted_point & p){
    os << "{" << p.weight() << ",{" <<
    p[0] << "," << p[1] << "," << p[2] << "}}";
    return os;
}
std::ostream & operator << (std::ostream & os, const CGT::RT & rt){
    using namespace CGT;
    const ScalarType degen = 1e-8; // TO DO : put elsewhere.
    
    
    os << "{" << '"' << "points" << '"' << "-> {";

    for(Finite_vertices_iterator vIt = rt.finite_vertices_begin(), vEnd = rt.finite_vertices_end();
        vIt!=vEnd; ++vIt)
        os << *vIt;

    os << "}," << '"' << "cells" << '"' << "-> {";

    for(Finite_cells_iterator cIt = rt.finite_cells_begin(), cEnd=rt.finite_cells_end();
        cIt!=cEnd; ++cIt){
        
        if(OnSameBoundary(cIt)) continue;
        
        std::array<Vector,Dimension+1> g;
        ScalarType volume;
        if(Degeneracy(cIt,volume,g)<=degen) continue;
        
        Vector grad(0,0,0);
        for(int i=0; i<=Dimension; ++i)
            grad = grad + Height(cIt->vertex(i)->point()) * g[i];
        
        std::array<IndexType,Dimension+1> indices;
        for(int i=0; i<=Dimension; ++i)
            indices[i]=cIt->vertex(i)->info().index;
        
        /*{
            std::array<Point,Dimension+1> pts;
            for(int i=0; i<=Dimension; ++i) pts[i]=cIt->vertex(i)->point().point();
            
            std::ostream & os = std::cout << "\n\n"
            ExportArrayArrow(indices)
            ExportArrayArrow(pts)
            ExportArrayArrow(g);
            print_range(os,grad.cartesian_begin(),grad.cartesian_end());
            os << "\n\n";
        }*/
        
	
        os << "{";
		for(auto it=indices.begin(); it!=indices.end(); ++it) os<<*it;
//        print_range(os,indices.begin(),indices.end());
        os << ",";
		for(auto it=grad.cartesian_begin(); it!=grad.cartesian_end(); ++it) os<<*it;
//        print_range(os,grad.cartesian_begin(),grad.cartesian_end());
        os << "},";
    }
    
    os << "}}";
    return os;
}


#endif
