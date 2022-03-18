//
//  Constraint.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 29/01/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_Constraint_h
#define CGalTest_Constraint_h

#include <vector>
#include <limits>
#include <iostream>

#include <numeric>
//#include "Macros.h"
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
        
        int error=1; // error==0 -> constraint satisfied, else constraint violated.
        int numberOfConstraints=BadIndex;
        ScalarType logSum = Infinity;
        std::vector<ScalarType> logGrad;
        std::vector<MatCoef> logHessian; // Upper triangular part only
        
        std::vector<ScalarType> values;
        std::vector<MatCoef> jacobian;
        std::vector<TensorCoef> hessian;
        
        
        // Glue code : Barrier for the constraint
        virtual ScalarType Value(const NS::VectorType &);
        virtual const NS::VectorType & Gradient(); // At latest position.
        virtual const NS::SparseMatrixType & Hessian();

        virtual std::string Name() const {return "Unspecified constraint name";}
        virtual void PrintSelf(std::ostream & os) const;
        friend std::ostream & operator << (std::ostream & os, const ConstraintType & c){
            c.PrintSelf(os); return os;}
        
        // nlopt interface
        static ScalarType GeometricMean(const std::vector<ScalarType> &, std::vector<ScalarType> &, void*);
    protected:
        virtual void Clean(FlagType);
        virtual void Check(FlagType);
        
        virtual void ComputeValJacHess(FlagType)
            {throw "ConstraintType::ComputeValJacHess error : must be specialized";}
        virtual void ComputeLogarithms(FlagType);
        
        NS::VectorType grad_; NS::SparseMatrixType hess_;
        virtual size_t InputSize() const {return values.size();} //!! incorrect in general: numCons != numVar
    };
    
    struct PositivityConstraint : ConstraintType {
        virtual void SetValues(const std::vector<ScalarType> &);
        virtual void ComputeValJacHess(FlagType);
        virtual std::string Name() const {return "PositivityConstraint";}
    protected:
        std::vector<ScalarType> val;
    };

#include "Constraint.hxx"
}

#endif
