//
//  Minkowski_3_Test.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 18/08/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_Minkowski_3_Test_h
#define CGalTest_Minkowski_3_Test_h

#include "Minkowski_3_Traits.h"
#include "MinkowskiConstraint_3.h"
#include <CGAL/point_generators_3.h>
#include <list>


namespace Minkowski_3 {
    
    // compute the tangent plane of a point
    template <typename K>
    typename K::Plane_3 tangent_plane (typename K::Point_3 const& p) {
        typename K::Vector_3 v(p.x(), p.y(), p.z());
        v = v / sqrt(v.squared_length());
        typename K::Plane_3 plane(v.x(), v.y(), v.z(), -(p - CGAL::ORIGIN) * v);
        return plane;
    }
    
    /*
    void Test(){
        using namespace CGT;
        
        // number of generated planes
        int N = 50;
        // generates random planes on a sphere
        std::list<K::Plane_3> planes;
        CGAL::Random_points_on_sphere_3<Point> g;
        for (int i = 0; i < N; i++) {
            planes.push_back(tangent_plane<K>(*g++));
        }
        // define polyhedron to hold the intersection
        Polyhedron P;
        // compute the intersection
        // if no point inside the intersection is provided, one
        // will be automatically found using linear programming
        CGAL::halfspace_intersection_3(planes.begin(),
                                       planes.end(),
                                       P,
                                       boost::make_optional(Point(0, 0, 0)) );
        
        for(auto plane : planes){
            std::cout << plane.a() << "," << plane.b() << "\n";
        }
        
        for(auto fit=P.facets_begin(); fit!=P.facets_end(); ++fit){
            std::cout <<
            fit->plane().a() << "," <<
            fit->plane().b() << "," <<
            fit->plane().c() << "," <<
            fit->plane().d() << "," <<
            "\n";
            
            auto eit = fit->facet_begin();
            std::cout << eit->vertex()->point().x() << "\n";
        }
    }*/
    
    void Test2(){
        using namespace CGT;
        typedef CGAL::Creator_uniform_3<double, Point>  PointCreator;
        
        CGAL::Random_points_in_sphere_3<Point, PointCreator> gen(100.0);
        // generate 250 points randomly on a sphere of radius 100.0
        // and copy them to a vector
        std::vector<Point> points;
        CGAL::cpp11::copy_n( gen, 250, std::back_inserter(points) );
        // define polyhedron to hold convex hull
        Polyhedron poly;
        
        // compute convex hull of non-collinear points
        CGAL::convex_hull_3(points.begin(), points.end(), poly);
        std::cout << "The convex hull contains " << poly.size_of_vertices() << " vertices" << std::endl;
        
        for(auto vit=poly.vertices_begin(); vit!=poly.vertices_end(); ++vit){
            std::cout << vit->point().x() <<
            "," << " index : " << vit->index << "\n";
        }
    }
    
    std::vector<std::array<ScalarType,4> >
    cubeEqns =
    {{0.,0.,1.,0.5},{1.,0.,0.,0.5},{0.,1.,0.,0.5},{-1.,0.,0.,0.5},{0.,0.,-1.,0.5},{0.,-1.,0.,0.5}},
    octahedronEqns =
    {{0.57735,-0.57735,0.57735,0.408248},{0.57735,0.57735,0.57735,0.408248},{-0.57735,0.57735,0.57735,0.408248},{-0.57735,-0.57735,0.57735,0.408248},{-0.57735,-0.57735,-0.57735,0.408248},{0.57735,-0.57735,-0.57735,0.408248},{-0.57735,0.57735,-0.57735,0.408248},{0.57735,0.57735,-0.57735,0.408248}},
    tetraEqns =
    {{0.,0.,-1.,0.204124},{-0.942809,0.,0.333333,0.204124},{0.471405,-0.816497,0.333333,0.204124},{0.471405,0.816497,0.333333,0.204124}},
    dodecaEqns =
    {{-0.894427,0.,-0.447214,1.11352},{0.894427,0,0.447214,1.11352},{0.276393,-0.850651,0.447214,1.11352},{0.,0.,1.,1.11352},{0.276393,0.850651,0.447214,1.11352},{0.723607,0.525731,-0.447214,1.11352},{0.723607,-0.525731,-0.447214,1.11352},{-0.276393,0.850651,-0.447214,1.11352},{0.,0.,-1.,1.11352},{-0.276393,-0.850651,-0.447214,1.11352},{-0.723607,-0.525731,0.447214,1.11352},{-0.723607,0.525731,0.447214,1.11352}};
    
    void Test3(){ // Volumes of platonician solids are correctly recovered
        using namespace CGT;
        const auto & eqns = dodecaEqns; //tetraEqns; //octahedronEqns; //cubeEqns;
        
        std::vector<CGT::Vector> pts;
        std::vector<ScalarType> val;
        for(auto eq : eqns) {pts.push_back({eq[0],eq[1],eq[2]}); val.push_back(eq[3]);}
        ConvexityConstraint cstr(pts);
        cstr.SetValues(val);
        cstr.Compute(31);
        
        std::cout
        ExportVarArrow(cstr.volume)
        ExportVarArrow(cstr);
    }

}

#endif
