//
//  NewtonSolvers.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 11/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_NewtonSolvers_h
#define CGalTest_NewtonSolvers_h

// Optimizing functionnals with and without constraints, with a basic Newton algorithm and logarithmic barriers.
// (I know that I am reinventing the wheel, but it is surprisingly difficult to find optimizers,
// which accept constraints with changing hessian sparsity pattern.)


//#include "ConvexityConstraint_Traits.h"

#include <iostream>
#include <limits>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include "JMM_CPPLibs/Macros/ExportArrow.h"

namespace NewtonSolvers {
//    using namespace ConvexityConstraint_Traits; // Just for Scalar, Index, MatCoef
    
    // Linear systems solved with Eigen. (Basic conjugate gradient.)
    typedef double ScalarType;
    const ScalarType Infinity = std::numeric_limits<ScalarType>::infinity();
    typedef Eigen::SparseMatrix<ScalarType> SparseMatrixType;
    typedef Eigen::Matrix<ScalarType,Eigen::Dynamic,1> VectorType;

    
    template<int Dimension>
    std::vector<ScalarType> StdFromEigen(const Eigen::Matrix<ScalarType,Dimension,1> & u){
        std::vector<ScalarType> u_(u.size());
        for(int i=0; i<u.size(); ++i) u_[i]=u[i];
        return u_;}
    
    template<int Dimension = Eigen::Dynamic>
    Eigen::Matrix<ScalarType,Dimension,1> EigenFromStd(const std::vector<ScalarType> & u){
        Eigen::Matrix<ScalarType,Dimension,1> u_(u.size());
        for(int i=0; i<u.size(); ++i) u_[i]=u[i];
        return u_;}
    
struct Functionnal {
    // Return +infty if out of domain.
    virtual ScalarType Value(const VectorType &){throw "Functionnal must be specialized";}
    virtual const VectorType & Gradient()=0; // At latest position.
    virtual const SparseMatrixType & Hessian()=0;
    static ScalarType Evaluate(const std::vector<ScalarType> &, std::vector<ScalarType> &, void *);
};

struct StoppingCriterion {
    ScalarType stopBelow=-Infinity;
    ScalarType stopAbove=Infinity;
    int persistence=1; // Activate if inserted value violates constraints persistence successive times.
    enum Priority {Top,Medium,Ignore} priority=Top;
    
    void Insert(ScalarType);
    bool Abort(bool &) const;
    std::vector<ScalarType> values;
    bool active=false;
    
    void Clear();
    void PrintSelf(std::ostream & os) const;
    friend std::ostream & operator << (std::ostream & os, const StoppingCriterion & a){
        a.PrintSelf(os); return os;}
};
    
struct NewtonUnconstrained {
    int verbose=0;
    int maxIter=20;
    int iter=0;
    ScalarType dampingRatio = 0.5;
    enum class DampingType {Valid, Increasing} dampingStrategy = DampingType::Valid;
    enum class SolverType {ConjugateGradient, SimplicialLLT, Other} solverStrategy = SolverType::SimplicialLLT;
    
    StoppingCriterion sObjective, sGHNorm, sGHDecay, sDelta, sGNorm, sDirNorm;
    
    NewtonUnconstrained(){
        sObjective.stopAbove = std::numeric_limits<ScalarType>::max();
        sDelta.stopBelow = 0.01;
        sDelta.persistence = 3;
    }
    
    void Solve(Functionnal & pb, VectorType & x);

    std::ostream * runtimeOut = &std::cout;
    
    void Clear();
    void PrintSelf(std::ostream & os) const;
    friend std::ostream & operator << (std::ostream & os, const NewtonUnconstrained & a){
        a.PrintSelf(os); return os;}
};
    
struct NewtonConstrained : Functionnal {
    NewtonUnconstrained opt;
    
    int maxIter=20; // Adds sub-iterations in opt
    int iter=0;
    bool verbose=false;
    
    ScalarType multiplier = 1;
    ScalarType multiplierBound = 1e-4;
    ScalarType multiplierDamping = 0.5;
    ScalarType ghNormBase = 1;
    /*
    struct ConstraintType {
        ScalarType initialMultiplier=1;
        bool divideByNumberOfConstraints = true;
        Functionnal c;
    };
     */
    
    virtual ScalarType Value(const VectorType &);
    virtual const VectorType & Gradient();
    virtual const SparseMatrixType & Hessian();
    
    void Solve(Functionnal &, VectorType &, std::vector<Functionnal*> ={});
    
    std::ostream * runtimeOut = &std::cout;
    void PrintSelf(std::ostream & os) const;
    friend std::ostream & operator << (std::ostream & os, const NewtonConstrained & a){
        a.PrintSelf(os); return os;}

    
protected:
    ScalarType value_; VectorType gradient_; SparseMatrixType hessian_;
    
    Functionnal * objective;
    std::vector<Functionnal*> barriers;

};
    
    
#include "NewtonSolvers.hxx"
}

#endif