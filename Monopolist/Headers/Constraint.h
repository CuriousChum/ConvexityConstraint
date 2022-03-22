#pragma once
//
//  Constraint.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 29/01/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#include <vector>
#include <limits>
#include <iostream>

#include <numeric>
#include "NewtonSolvers.h"

namespace Constraint {
    namespace NS = NewtonSolvers;
    
    typedef NS::ScalarType ScalarType;
    typedef int IndexType;
    typedef unsigned int FlagType; // Conjunction of requests

    const ScalarType Infinity = std::numeric_limits<ScalarType>::infinity();
    const IndexType BadIndex = std::numeric_limits<IndexType>::max();
    const IndexType InfiniteIndex = BadIndex-1;

    typedef std::pair<IndexType, ScalarType> VecCoef;
    struct MatCoef;
    struct TensorCoef; // Trivariant tensor
    
    static void Sqrt(ScalarType &, std::vector<VecCoef> &, std::vector<MatCoef> &); // true hessian
    
    template<typename Coef> struct SameIndices;

    template<typename Coef, typename Traits = SameIndices<Coef> >
    void Simplify(std::vector<Coef> &,ScalarType tol=0);
    void Sqrt(ScalarType & val, std::vector<VecCoef> & g, std::vector<MatCoef> & h);
    void NLog(ScalarType & val, std::vector<VecCoef> & g, std::vector<MatCoef> & h); //-log
    
    
    // Barrier function - Sum_i log(c_i(x)), for the constraints c_i(x)>0 for all i.
    struct ConstraintType : NS::Functionnal {
        enum Request {
            RLogSum = 1<<0, RLogGrad = 1<<1, RLogHessian = 1<<2,
            RValues = 1<<3, RJacobian = 1<<4, RHessian = 1<<5
        };
        virtual void SetValues(const std::vector<ScalarType> &)
            {throw "ConstraintType::Compute error : must be specialized";};
        virtual void Compute(FlagType);
        
        IndexType numberOfConstraints=BadIndex, numberOfUnknowns=BadIndex;
        ScalarType logSum = Infinity;
        std::vector<ScalarType> logGrad;
        std::vector<MatCoef> logHessian; // Upper triangular part only
        
        std::vector<ScalarType> values; // The values of the constraints
        std::vector<MatCoef> jacobian;
        std::vector<TensorCoef> hessian;
        
        // Glue code : Barrier for the constraint
        virtual ScalarType Value(const NS::VectorType &) override;
        virtual const NS::VectorType & Gradient() override; // At latest position.
        virtual const NS::SparseMatrixType & Hessian() override;

        virtual std::string Name() const override {return "Unspecified constraint name";}
        virtual void PrintSelf(std::ostream & os) const;
        friend std::ostream & operator << (std::ostream & os, const ConstraintType & c){
            c.PrintSelf(os); return os;}
        
        // nlopt interface
        static ScalarType GeometricMean(const std::vector<ScalarType> &,
										std::vector<ScalarType> &, void*);
    protected:
        virtual void Clean(FlagType);
        virtual void Check(FlagType);
        
        virtual void ComputeValJacHess(FlagType)
            {throw "ConstraintType::ComputeValJacHess error : must be specialized";}
        virtual void ComputeLogarithms(FlagType);
        
        NS::VectorType grad_; NS::SparseMatrixType hess_;
    };
    
#include "Constraint.hxx"
}
