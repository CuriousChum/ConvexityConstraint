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
	/**CGAL points by default do not carry an index, but this  is needed in a number of settings.
	 There are two workarounds (both used here actually) :
	  - Use an std::map<Point,int>. This works since CGAL never rounds or otherwise modifies the input points. However, there is a small computational overhead.
	  - Define your own class that carries an index, which is the approach used here. In fact this approach is described in the CGAL manual. The syntax is a bit obscure but it works.
	 */
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef K::Point_3 Point;
        typedef K::Vector_3 Vector;
        
        // A CGAL-compatible vertex type with an index
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
		/**
		 Computes the lengths of the edges of the three dimensional polyhedron defined by the linear inequality constraints
		 H(a) = {z in R^3; < z, pts[i] >  <= a[i] for all indices i}
		 The "pts" are fixed, and the heights "a" are meant to be optimized.
		 When these areas are equal to a given value, the Minkowski problem is solved.
		 
		 We rely on CGAL to compute the dual convex K(a) body to H(a), using a convex hull procedure.
		 K(a) = Hull {pts[i]/a[i] for all indices i}
		 
		 Note the translation invariance H(a) + u = H( a[i] + <u,pts[i]>; 0<= i < I ) for any vector u in R^3,
		 hence Area (H(a)) = Area( H( a[i] + <u,pts[i]>; 0<= i < I ) ).
		 
		 See Minkowski_2:Functional for a similar but simpler class.
		 */

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
        
        ScalarType volume; // Lebesgue volume of H(a)
        N::VectorType facetAreas; // area gradient
        std::vector<N::Triplet> jacobian; // area jacobian
		
        G::Vector moment; // Moment of H(a), i.e. the integral of position.
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
