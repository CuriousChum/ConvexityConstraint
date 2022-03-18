//
//  PrincipalAgent_3.hxx
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 13/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_PrincipalAgent_3_hxx
#define CGalTest_PrincipalAgent_3_hxx

PrincipalAgent::PrincipalAgent(const CGT::RT & rt, ScalarType degen)
//:minDegeneracy(degen_)
{
    using namespace CGT;
    const size_t n = rt.number_of_vertices();
    std::vector<Eigen::Triplet<ScalarType> > triplets;
    
    NS::VectorType integral(n);
    integral.fill(0.);
    
    for(Finite_cells_iterator fIt = rt.finite_cells_begin(), fEnd = rt.finite_cells_end();
        fIt!=fEnd; ++fIt){
        
        std::array<Vector,Dimension+1> g;
        ScalarType volume;
        if(OnSameBoundary(fIt)) continue;
        if(Degeneracy(fIt,volume,g) <= degen) continue; // omit
        
        std::array<IndexType,Dimension+1> indices;
        for(int i=0; i<=Dimension; ++i) {
            const IndexType index = fIt->vertex(i)->info().index;
            assert(0<=index && index < n);
            indices[i]=index;
        }
        
        
        // Hessian
        for(int i=0; i<=Dimension; ++i)
            for(int j=0; j<=Dimension; ++j)
                triplets.push_back({indices[i],indices[j], volume*g[i]*g[j]});
        
        // Integral
        for(int i=0; i<=Dimension; ++i)
            integral[indices[i]] += volume/(Dimension+1);
    }
    
    hess_.resize((int)n,(int)n);
    hess_.setFromTriplets(triplets.begin(), triplets.end());

    NS::VectorType halfParabola(n);
    halfParabola.fill(std::numeric_limits<ScalarType>::quiet_NaN());
    
    // Parabola
    for(Finite_vertices_iterator vIt = rt.finite_vertices_begin(), vEnd = rt.finite_vertices_end();
        vIt!=vEnd; ++vIt)
        halfParabola(vIt->info().index) = 0.5*Parabola(vIt->point());
    
    if(halfParabola.hasNaN()) throw "Principal agent error : invalid vertices indices";
    
    grad0_ = integral - hess_*halfParabola;
}

ScalarType PrincipalAgent::Value(const NS::VectorType & x){
    const NS::VectorType hx = hess_*x;
    grad_ = hx + grad0_;
    return x.dot(grad0_+0.5*hx);
}

void PrincipalAgent::PrintSelf(std::ostream & os) const {
    std::vector<ScalarType> grad0;
    for(int i=0; i<grad0_.size(); ++i)
        grad0.push_back(grad0_[i]);
    
    std::vector<MatCoef> hess;
    for (int k=0; k<hess_.outerSize(); ++k)
        for (NS::SparseMatrixType::InnerIterator it(hess_,k); it; ++it)
            hess.push_back(MatCoef{(IndexType)it.row(), (IndexType)it.col(), it.value()});
    
    os << "{"
    ExportArrayArrow(grad0)
    ExportArrayArrow(hess)
    << "}";
}

#endif
