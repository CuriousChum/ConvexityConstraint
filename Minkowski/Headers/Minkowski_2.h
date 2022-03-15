//
//  Minkowski_2.h
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 02/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef Minkowski_Minkowski_2_h
#define Minkowski_Minkowski_2_h

#include "QuotientedNewton.h"

namespace Minkowski_2 {

    namespace N = QuotientedNewton;
    typedef N::ScalarType ScalarType;
    const int Dimension=2;

    namespace G { // geometry
        typedef Eigen::Matrix<ScalarType, Dimension, 1> Vector;
        typedef Eigen::Matrix<ScalarType, Dimension, Eigen::Dynamic> VectorList;

        ScalarType ScalarProduct(const Vector & u, const Vector & v){return u[0]*v[0]+u[1]*v[1];}
        ScalarType Determinant(const Vector & u, const Vector & v){return u[0]*v[1]-u[1]*v[0];}
        ScalarType Orientation(const Vector & u, const Vector & v, const Vector & w){return Determinant(v-u,w-u);}
        Vector Perp(const Vector & v){return Vector{-v[1],v[0]};}
    }
    
    struct Functional : N::FunctionalBase {
		/**
		 Computes the lengths of the edges of the two dimensional polyhedron defined by the linear inequality constraints
		 H(a) = {z in R^2; < z, pts[i] >  <= a[i] for all indices i}
		 The "pts" are fixed, and the heights "a" are meant to be optimized.
		 When these areas are equal to a given value, the Minkowski problem is solved.
		 
		 Note the translation invariance H(a) + u = H( a[i] + <u,pts[i]>; 0<= i < I ) for any vector u in R^2,
		 hence Area (H(a)) = Area( H( a[i] + <u,pts[i]>; 0<= i < I ) ).
		 */
		
        enum class ErrCode {NoError=0,Orientation=1,Sign=2,Size=3};

        virtual int Compute(N::VectorType &, N::VectorType &);
        virtual int DescentDirection(N::VectorType &);
        
        G::VectorList pts; // Assumed to be ordered trigonometrically
        std::array<int,2> referenceIndices={N::BadIndex,N::BadIndex}; // Should be automatically set
        void CheckData();
        int Size() const {return (int)pts.cols();}
        G::Vector Barycenter() const {return moment/area;}
        
        // Accessors
        ScalarType Area() const {return area;}
        const N::VectorType & FacetLengths() const {return facetLengths;}
        const std::vector<N::Triplet> Jacobian() const {return jacobian;}
    protected:
        
        ScalarType area; // Lebesgue area of the convex body
        N::VectorType facetLengths; // area gradient
        std::vector<N::Triplet> jacobian; // area jacobian
		
        G::Vector moment; // First order moment of the convex body (integral of position).
        G::VectorList facetMoments; // moment gradient
        
        ErrCode CheckInput(const N::VectorType &) const;
        void Compute(const N::VectorType &);
		/// Modifies a in such way that the convex body H(a) is centered on the origin.
        void Normalize(N::VectorType & a);

        void SetReferenceIndices();
        ScalarType ReferenceIndicesConditionNumber();
        ScalarType referenceIndicesConditionNumberThreshold = 0.2;
    };
    
#include "Minkowski_2.hxx"
}


#endif
