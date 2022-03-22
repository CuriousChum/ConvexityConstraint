#pragma once
//
//  ConvexityConstraint_Traits.hxx
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 29/01/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

// %%%%%%%%%%%%%% Coefficient types %%%%%%%%%%%%

struct MatCoef {
    IndexType i,j;
    ScalarType v;
    void PrintSelf(std::ostream & os) const {os << "{" << i << "," << j << "," << v << "}";}
    friend std::ostream & operator << (std::ostream & os, const MatCoef & c){
        c.PrintSelf(os); return os;}
    typedef std::tuple<IndexType,IndexType,ScalarType> TupleType;
    bool operator < (const MatCoef & c) const {return TupleType(i,j,v) < TupleType(c.i,c.j,c.v);}
};

struct TensorCoef {
    IndexType i,j,k;
    ScalarType v;
    void PrintSelf(std::ostream & os) const {
        os << "{" << i << "," << j << "," << k << "," << v << "}";}
    friend std::ostream & operator << (std::ostream & os, const TensorCoef & c){
        c.PrintSelf(os); return os;}
    typedef std::tuple<IndexType,IndexType,IndexType,ScalarType> TupleType;
    bool operator < (const TensorCoef & c) const {
        return TupleType(i,j,k,v) < TupleType(c.i,c.j,c.k,c.v);}
};


// %%%%%%%%%%%%%%% Eliminate redundancies %%%%%%%%%%%%%%

template<> struct SameIndices<VecCoef> {
    bool Same(const VecCoef & a, const VecCoef & b){
        return std::get<0>(a)==std::get<0>(b);}
    ScalarType & Val(VecCoef & a){return std::get<1>(a);}
};


template<> struct SameIndices<MatCoef> {
    bool Same(const MatCoef & a, const MatCoef & b){
        return a.i==b.i && a.j==b.j;}
    ScalarType & Val(MatCoef & a){return a.v;}
};

template<> struct SameIndices<TensorCoef> {
    bool Same(const TensorCoef & a, const TensorCoef & b){
        return a.i==b.i && a.j==b.j && a.k==b.k;}
    ScalarType & Val(TensorCoef & a){return a.v;}
};


template<typename Coef, typename Traits>
void Simplify(std::vector<Coef> & v,ScalarType tol){
    if(v.empty()) return;
    std::sort(v.begin(),v.end());
    Traits traits;
    
    auto i=v.begin(), j=v.begin();
    ++i;
    for(; i!=v.end(); ++i){
        if(fabs(traits.Val(*i))<=tol)
            continue;
        else if(traits.Same(*i,*j))
            traits.Val(*j)+=traits.Val(*i);
        else {
            ++j;
            if(i!=j)
                *j=*i;
        }
    }
    v.erase(++j,v.end());
}

void Sqrt(ScalarType & val, std::vector<VecCoef> & g, std::vector<MatCoef> & h){
    
    const ScalarType sVal = sqrt(val);
    if(sVal==0) throw "Sqrt error : differentiation at 0";
    const ScalarType a1 =  1./(2.*sVal);
    const ScalarType a2 = -1./(4.*sVal*val);
    
    for(MatCoef & c : h)
        c.v *= a1;
    for(const VecCoef v1 : g)
        for(const VecCoef v2 : g)
            h.push_back(MatCoef{v1.first,v2.first, a2 * v1.second * v2.second});
    
    for(VecCoef & c : g)
        c.second *= a1;
    
    val=sVal;
}

void NLog(ScalarType & val, std::vector<VecCoef> & g, std::vector<MatCoef> & h){
    
    /*
    std::ostream & os = std::cout << "Begin "
    ExportVarArrow(val)
    ExportArrayArrow(g)
    ExportArrayArrow(h)
    << "\n";*/
    
    if(val==0) throw "Log error: evaluation at 0";
    const ScalarType iNV =  -1./val;
    const ScalarType iV2 = iNV*iNV;
    
    for(MatCoef & c : h)
        c.v *= iNV;
    for(const VecCoef v1 : g)
        for(const VecCoef v2 : g)
            h.push_back(MatCoef{v1.first,v2.first, iV2 * v1.second * v2.second});
    
    for(VecCoef & c : g)
        c.second *= iNV;
    
    /*
    os << "End "
    ExportVarArrow(val)
    ExportArrayArrow(g)
    ExportArrayArrow(h)
    << "\n";*/
    
    val=-log(val);
}

// %%%%%%%%%%%%%%% Constraint Type %%%%%%%%%%%%%%

void ConstraintType::Clean(FlagType r){
//    assert(!error);
    if(r & RLogSum) logSum=0;
    if(r & RLogGrad) logGrad.clear();
    if(r & RLogHessian) logHessian.clear();
    if(r & RValues) values.clear();
    if(r & RJacobian) jacobian.clear();
    if(r & RHessian) hessian.clear();
}

void ConstraintType::Check(FlagType r){
    assert(!error);
    if(r & RLogSum) assert(logSum<Infinity);
    if(r & RLogGrad) assert(!logGrad.empty());
    if(r & RLogHessian) assert(!logHessian.empty());
    if(r & RValues) assert(!values.empty());
    if(r & RJacobian) assert(!jacobian.empty());
}


ScalarType ConstraintType::Value(const NS::VectorType & x){
    Clean(RLogGrad | RLogHessian);
    std::vector<ScalarType> y(x.size());
    for(int i=0; i<x.size(); ++i) y[i]=x[i];
    SetValues(y);
    if(error!=0) return Infinity;
    Compute(RLogSum);
    return logSum;
}

const NS::VectorType & ConstraintType::Gradient(){
    if(logGrad.empty()) {
        Compute(RLogGrad | RLogHessian);
        const int n=(int)logGrad.size();
        grad_.resize(n);
        hess_.resize(0,0); hess_.resize(n,n); // should clear values
        for(int i=0; i<logGrad.size(); ++i) grad_[i]=logGrad[i];
        
        std::vector<Eigen::Triplet<ScalarType> > triplets;
        triplets.reserve(logHessian.size());
        
        ScalarType maxv=0;
        for(auto c : logHessian){
            triplets.push_back({c.i,c.j,c.v});
            maxv=std::max(maxv,c.v);
        }
        /*
        std::cout
        ExportVarArrow(maxv)
        ExportVarArrow(logHessian.size())
        ExportVarArrow(Name())
        ExportVarArrow(logGrad.size())
        << "\n";*/
        
        hess_.setFromTriplets(triplets.begin(),triplets.end());
    }
    return grad_;
}

const NS::SparseMatrixType & ConstraintType::Hessian(){
    if(logGrad.empty()) Gradient();
    return hess_;
}


void ConstraintType::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(error)
    ExportVarArrow(numberOfConstraints)
    ExportVarArrow(logSum)
    ExportArrayArrow(logGrad)
    ExportArrayArrow(logHessian)
    ExportArrayArrow(values)
    ExportArrayArrow(jacobian)
    ExportArrayArrow(hessian)
    << "}";
}

void ConstraintType::Compute(FlagType r0){
    FlagType r=r0;
    if(r & (RLogSum | RLogGrad | RLogHessian) ) r = r | RValues;
    if(r & (RLogGrad | RLogHessian) ) r = r | RJacobian;
    if(r & RLogHessian) r = r | RHessian;
    
    Clean(r);
    ComputeValJacHess(r);
    
    Simplify(jacobian);
    Simplify(hessian);
    ComputeLogarithms(r0);
};

// %%%%%%%%%%%%%%%%%%%%%%%%%%% NLOPT interface %%%%%%%%%%%%%%%%%

ScalarType ConstraintType::GeometricMean(const std::vector<ScalarType> & input, std::vector<ScalarType> & grad, void*data){
    ConstraintType & cons = *static_cast<ConstraintType*>(data);
    std::ostream & os = std::cout ExportArrayArrow(input) << "\n";
    cons.SetValues(input);
    if(cons.error) {std::cout << cons.Name() << ", error in SetValues\n"; return 1.;}
    cons.Compute(RValues | RJacobian);
    if(cons.error) {std::cout << cons.Name() << ", error in Compute\n"; return 0.;};
    
    const std::vector<ScalarType> & x = cons.values;
    const std::vector<MatCoef> & jac = cons.jacobian;
    
    std::fill(grad.begin(),grad.end(),0.);
    if(*std::min_element(x.begin(),x.end())<=0)
        {std::cout << cons.Name() << ", error in Logarithms\n"; return 0.;};
    
    // Taylor development for the sum of logarithms
    if(!grad.empty())
        for(const MatCoef & t : jac)
            grad[t.j] += t.v/x[t.i];
    ScalarType logSum = 0;
    int n=0;
    for(ScalarType v : x) if(v!=Infinity) {logSum += log(v); ++n;}
    if(n==0) return 1;
    
    // Normalization
    for(ScalarType & g : grad) g/=n;
    logSum/=n;
    
    // Exponentiation
    ScalarType result = exp(logSum);
    for(ScalarType & g : grad) g*=result;
    
    
    /*        std::ostream & os = std::cout
     ExportVarArrow(result)
     << "\n";*/
    
    std::cout << "ConstraintType::GeometricMean, value : " << result << ", " << cons.Name() << "\n";
    
    // Negating, as required by NLOPT convention
    result*=-1;
    for(ScalarType & g : grad) g*=-1.;

    return result;
}



// %%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute Logarithms %%%%%%%%%%%%%%%%%%%%%%%%

void ConstraintType::ComputeLogarithms(FlagType r){
    // Jacobian and hessian have been computed and sorted. Cleanup done.
    if(r & RLogGrad) logGrad.resize(InputSize(),0);
    
    ScalarType tValue;
    std::vector<VecCoef> tJacobian;
    std::vector<MatCoef> tHessian;
    
    std::vector<MatCoef>::iterator jacIt = jacobian.begin();
    std::vector<TensorCoef>::iterator hessIt = hessian.begin();
    
    for(int i=0; i<values.size(); ++i){
        tValue = values[i];
        if(tValue<0) {error=1; logSum=Infinity; return;}
        if(tValue==Infinity)
            continue;
        
        tJacobian.clear();
        tHessian.clear();
        
        for(;jacIt!=jacobian.end() && jacIt->i == i; ++jacIt)
            tJacobian.push_back({jacIt->j,jacIt->v});
        
        for(;hessIt!=hessian.end() && hessIt->i == i; ++hessIt)
            tHessian.push_back({hessIt->j,hessIt->k,hessIt->v});
        
        // Compute logarithm
        
        NLog(tValue,tJacobian,tHessian);
        
        Simplify(tJacobian);
        Simplify(tHessian);
        
        // Save computations
        logSum += tValue;
        if(r & RLogGrad)
            for(const VecCoef & c : tJacobian)
                logGrad[c.first] += c.second;
        if(r & RLogHessian)
            logHessian.insert(logHessian.end(), tHessian.begin(), tHessian.end());
    }
    
    assert(jacIt==jacobian.end());
    assert(hessIt==hessian.end());
    
    Simplify(logHessian);
}
