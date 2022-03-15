//
//  Minkowski_3.h
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 03/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef Minkowski_Minkowski_3_h
#define Minkowski_Minkowski_3_h

#include "QuotientedNewton.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>

namespace Minkowski_3 {
    namespace N = QuotientedNewton;
    typedef N::ScalarType ScalarType;
    typedef N::IndexType IndexType;
    const int Dimension=3;

    namespace G {
        typedef Eigen::Matrix<ScalarType, Dimension, 1> Vector;
        typedef Eigen::Matrix<ScalarType, Dimension, Eigen::Dynamic> VectorList;
    }
    
    namespace CGT {
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef K::Point_3 Point;
        typedef K::Vector_3 Vector;
        
        // A vertex type with an index
        template <class Refs>
        struct My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> {
            typedef CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> superclass;
            IndexType index=N::BadIndex;
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
    
    struct Functional : N::FunctionalBase {
        enum class ErrCode {NoError=0,InteriorPoints,NonPositiveEntries};
        
        virtual int Compute(N::VectorType &, N::VectorType &);
        virtual int DescentDirection(N::VectorType &);
        
        G::VectorList pts;
        std::array<int,3> referenceIndices={N::BadIndex,N::BadIndex,N::BadIndex};
        //bool CheckData(){return true;}
        int Size() const {return (int)pts.cols();}
        
        //Accessors
        G::Vector Barycenter() const {return moment/volume;}
        ScalarType Volume() const {return volume;}
        const N::VectorType & FacetAreas() const {return facetAreas;}
        const std::vector<N::Triplet> Jacobian() const {return jacobian;}
    protected:
        CGT::Polyhedron poly;
        CGT::Vector Point(int i) const;
        
        ScalarType volume;
        N::VectorType facetAreas; // area gradient
        std::vector<N::Triplet> jacobian; // area jacobian
        G::Vector moment;
        G::VectorList facetMoments; // moment gradient
        
        ErrCode MakeHull(const N::VectorType &);
        void Compute(const N::VectorType &);
        void Normalize(N::VectorType &);
        
        void SetReferenceIndices();
        ScalarType ReferenceIndicesConditionNumber();
        ScalarType referenceIndicesConditionNumberThreshold = 0.1;
        
        
    };
    
    G::Vector GFromCGT(const CGT::Vector & u){return G::Vector{u.x(),u.y(),u.z()};}
    CGT::Vector CGTFromG(const G::Vector & u){return CGT::Vector(u[0],u[1],u[2]);}
    
#include "Minkowski_3.hxx"
}


#endif
