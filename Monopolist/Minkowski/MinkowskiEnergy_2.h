//
//  MinkowskiEnergy_2.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 21/08/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_MinkowskiEnergy_2_h
#define CGalTest_MinkowskiEnergy_2_h

#include "MinkowskiConstraint_2.h"

// The Minkowski energy is Sum lambda_i alpha_i - log(Vol(alpha)). It is convex.
// Hessian is not sparse as is, hence useless.


namespace Minkowski_2  {
    struct MinkowskiEnergy : NS::Functionnal {
        
        ConvexityConstraint cstr;
        std::vector<ScalarType> target;
        
        ScalarType val; // value, gradient and hessian (at latest position)
        std::vector<VecCoef> g; NS::VectorType g_;
        std::vector<MatCoef> h; NS::SparseMatrixType h_;
        
        ScalarType Value(const NS::VectorType &);
        virtual const NS::VectorType & Gradient() {return g_;}
        virtual const NS::SparseMatrixType & Hessian() {return h_;}
        int size() const {return (int)target.size();}
        MinkowskiEnergy(const std::vector<Vector> & pts, const std::vector<ScalarType> & target)
        :cstr(pts),target(target){};
        
        static ScalarType Evaluate(const std::vector<ScalarType> &, std::vector<ScalarType> &, void *);
    };
    
    
    ScalarType MinkowskiEnergy::Value(const NS::VectorType & v_){
        assert(v_.size()==size());
        std::vector<ScalarType> v; v.reserve(size());
        for(int i=0; i<size(); ++i) v.push_back(v_[i]);
        cstr.SetValues(v);
        if(cstr.error) return Infinity;
        cstr.Compute(ConstraintType::RValues | ConstraintType::RJacobian);
        
        // Compute logarithms
        val=cstr.area; 
        g.clear();  for(int i=0; i<size(); ++i) g.push_back({i,cstr.values[i]});
        h=cstr.jacobian;
        
        //std::ostream & os = std::cout ExportVarArrow(val) ExportArrayArrow(v) << "\n";
        NLog(val,g,h);
        // Issue : the hessian is NOT sparse anymore !!
        
        // Add the linear term
        for(int i=0; i<size(); ++i){
            val += v[i]*target[i];
            g.push_back({i,target[i]});
        }
        Simplify(g);

        
        // Convert to eigen format
        g_.resize(size()); g_.fill(0.); for(auto c : g) g_[c.first] += c.second;
        
        std::vector<Eigen::Triplet<ScalarType> > triplets;
        triplets.reserve(h.size());
        for(auto c : h) triplets.push_back({c.i,c.j,c.v});
        h_.resize(size(),size()); h_.setZero(); h_.setFromTriplets(triplets.begin(),triplets.end());
        
        return val;
    }
    
    ScalarType MinkowskiEnergy::Evaluate(const std::vector<ScalarType> & x,
                                         std::vector<ScalarType> & grad,
                                         void * data) {
        assert(data!=nullptr);
        MinkowskiEnergy & energy = *static_cast<MinkowskiEnergy*>(data);
        NS::VectorType x_(x.size()); for(int i=0; i<x.size(); ++i) x_[i]=x[i];
        const ScalarType result = energy.Value(x_);
        if(!grad.empty()){
            auto grad_ = energy.Gradient();
            for(int i=0; i<x.size(); ++i) grad[i]=grad_[i];
        }
        
        static int nIterFinite; if(result<INFINITY) nIterFinite++;
        static int nIter=0; nIter++;
        std::ostream & os =std::cout
        << "Evaluating"
        ExportVarArrow(result)
        << "\n"
        ExportArrayArrow(x)
        ExportVarArrow(nIter)
        ExportVarArrow(nIterFinite)
        << "\n\n";
        
        return result;
    }

}



#endif
