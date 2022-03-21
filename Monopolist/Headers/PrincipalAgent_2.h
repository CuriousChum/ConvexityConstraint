//
//  PrincipalAgent.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 28/01/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_PrincipalAgent_h
#define CGalTest_PrincipalAgent_h

#include "Geometry_2.h"

// Minimize int 0.5*\|grad u-x\|^2 + u, subject to u>=0, u convex.
// I may want to use CGAL, and hand made finite elements -> adaptive mesh refinement

namespace Geometry_2 {

struct PrincipalAgent : NS::Functionnal {
    CGT::RT rt;
    bool hasQuadraticCost = true;
    struct ObjectiveType {
        ScalarType constant;
        std::vector<ScalarType> linear;
        std::vector<MatCoef> quadratic; // hessian -> factor 1/2 when computing objective
        
        void PrintSelf(std::ostream & os) const;
        friend std::ostream & operator << (std::ostream & os, const ObjectiveType & a){
            a.PrintSelf(os); return os;}
    } terms;
    
    
    // Glue code : Barrier for the constraint
    virtual ScalarType Value(const NS::VectorType &);
    virtual const NS::VectorType & Gradient(){return grad_;} // At latest position.
    virtual const NS::SparseMatrixType & Hessian(){return hess_;}

    void PrintSelf(std::ostream & os) const;
    friend std::ostream & operator << (std::ostream & os, const PrincipalAgent & a){
        a.PrintSelf(os); return os;}
    
    std::vector<CGT::Full_point> Pts() const;
    void Refine(ScalarType, const NS::VectorType &); // Refine based on a quantile for error indicator.
    ScalarType ErrorIndicator(CGT::Face_handle, const NS::VectorType &) const;
    CGT::Vector Gradient(CGT::Face_handle, const NS::VectorType &) const;
protected:
    void MakeObjective();
    NS::VectorType grad_; NS::SparseMatrixType hess_;
    
    // sqrt(h_T) \| [grad u] \|_{L^2(\partial T)} . Estimator for H1 approx error.
};
    
#include "PrincipalAgent_2.hxx"
    
}


#endif
