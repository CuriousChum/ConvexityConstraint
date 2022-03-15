//
//  QuotientedNewton.h
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 02/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

// This class is designed to solve an inverse problem F(x)=y using the newton method.
// The functional of interest may have some invariance, thus


#ifndef Minkowski_QuotientedNewton_h
#define Minkowski_QuotientedNewton_h

#include <numeric>
#include <array>
#include <algorithm>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "JMM_CPPLibs/Macros/ExportArrow.h"
#include "JMM_CPPLibs/Output/EnumToString.h"


//#include <Eigen/SparseCholesky>


namespace QuotientedNewton {
    
    typedef int IndexType;
    const IndexType BadIndex = std::numeric_limits<IndexType>::max();

    typedef double ScalarType;
    const ScalarType Infinity = std::numeric_limits<ScalarType>::infinity();
    typedef Eigen::Matrix<ScalarType,Eigen::Dynamic,1> VectorType;
    typedef Eigen::Triplet<ScalarType> Triplet;
    typedef Eigen::SparseMatrix<ScalarType> SparseMatrixType;

    
    struct FunctionalBase {
        virtual int Compute(VectorType &, VectorType &) {assert(false);};
        virtual int DescentDirection(VectorType &) = 0;
        ScalarType constraintMultiplier = 1; // For constraint penalization
    };
    
    
    struct NewtonScheme {
        ScalarType delta=1;
        ScalarType deltaGranularity = 1.4;
        ScalarType deltaMin = 1e-5;
        
        // Stopping criteria. Leave on first reached.
        bool forceResidueDecrease=true;
        int maxIterations=200;
        
        struct TestType {
            ScalarType lowerBound;
            int delay;
            std::vector<ScalarType> values;
            bool operator()(ScalarType); // returns true if should terminate
            bool Active() const;
            void PrintSelf(std::ostream & os) const;
            friend std::ostream & operator << (std::ostream & os, const TestType & u){
                u.PrintSelf(os);return os;}
        };
        
        TestType stepSizeTest {1e-3,5};
        TestType residueDecayTest{0.001,5};
        TestType residueMinTest{0.00001,4};
        void Clear();
        
        enum class StoppingCriterion {MaxIterations, ResidueDecayTest, ResidueMinTest, StepSizeTest, InfiniteResidue, EvaluationError, InversionError};
        
        StoppingCriterion Solve(FunctionalBase & op, const VectorType & target, VectorType & x);
        
        void PrintSelf(std::ostream & os) const;
        friend std::ostream & operator << (std::ostream & os, const NewtonScheme & u){u.PrintSelf(os);return os;}
        int GetNumberOfEffectiveIterations() const {return numberOfEffectiveIterations;}
        
        // profiling
        struct TimingsType {
            std::vector<ScalarType> apply, applyDifferentiate, solve;
            void Clear();
            void PrintSelf(std::ostream & os) const;
            friend std::ostream & operator << (std::ostream & os, const TimingsType & u){
                u.PrintSelf(os);return os;}
        } timings;
    protected:
        int numberOfEffectiveIterations=0;
        std::vector<ScalarType> symmetryDefects;
    };

    
    
    // Eigen to Std conversion
    
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
    
#include "QuotientedNewton.hxx"
}

template<> char const* enumStrings<QuotientedNewton::NewtonScheme::StoppingCriterion>::data[] =
{"MaxIterations", "ResidueDecayTest", "ResidueMinTest", "StepSizeTest", "InfiniteResidue",
    "EvaluationError","InversionError"};


#endif
