#pragma once
//
//  NewtonConstrained.hxx
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 11/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

// %%%%%%%%%%%%%%% Functionnal interface to NLOPT %%%%%%%%%%%%%%%

ScalarType Functionnal::Evaluate(const std::vector<ScalarType> & x,
								 std::vector<ScalarType> & grad, void* data){
    Functionnal & func = *static_cast<Functionnal*>(data);
    const ScalarType result = func.Value(EigenFromStd(x));
    if(!grad.empty()) grad = StdFromEigen(func.Gradient());
    std::cout << "Functional::Evaluate, got value : " << result << "\n";
    return result;
}

// %%%%%%%%%%%%%%%%%%%% Stopping Criterion %%%%%%%%%%%%%%%%%%%%%%%

void StoppingCriterion::Insert(ScalarType v){
    values.push_back(v);
    if(values.size()>=persistence){
        const size_t
        nBelow = std::count_if(values.end()-persistence, values.end(),
                               [this](ScalarType s){return s<=stopBelow;}),
        nAbove = std::count_if(values.end()-persistence, values.end(),
                               [this](ScalarType s){return s>=stopAbove;});
        if(nBelow==persistence||nAbove==persistence) active=true;
    }
}

bool StoppingCriterion::Abort(bool & abort) const {
    if(priority==Top && active) return true;
    if(priority==Medium || priority==Top) abort = abort && active;
    return false;
}

void StoppingCriterion::Clear(){
    values.clear();
    active=false;
}

void StoppingCriterion::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(stopBelow)
    ExportVarArrow(stopAbove)
    ExportVarArrow(persistence)
    ExportVarArrow(priority)
    << "\n"
    ExportArrayArrow(values)
    ExportVarArrow(active)
    << "}";
}

// %%%%%%%%%%%%%%%%%%%%% Newton Unconstrained %%%%%%%%%%%%%%%%%%%%%%

void NewtonUnconstrained::Solve(Functionnal & pb, VectorType & x){
    VectorType dir; // descent direction
    while(maxIter > iter++){
        sObjective.Insert(pb.Value(x));
        const VectorType & grad = pb.Gradient();
        const SparseMatrixType & hess = pb.Hessian();
        
        if(solverStrategy==SolverType::ConjugateGradient){
            Eigen::ConjugateGradient<SparseMatrixType> cg;
            cg.compute(hess); 
            dir = cg.solve(grad);
            if(verbose) std::cout << "ConjugateGradient\n";
        } else if(solverStrategy==SolverType::SimplicialLLT) {
            Eigen::SimplicialLDLT<SparseMatrixType> solver(hess);
            dir = solver.solve(grad);
            if(verbose) std::cout << "SimplicialLDLT\n";
        } else {
//            Eigen::PastixLDLT<SparseMatrixType> solver(hess); // Too much to import.
            /*
            Eigen::BiCGSTAB<SparseMatrixType> solver;
            solver.compute(hess);
            dir = solver.solve(grad);
            std::cout << "BiCGSTAB\n";
             */ // Does not work at all ???
            assert(hess.isCompressed());
            Eigen::SparseQR<SparseMatrixType, Eigen::COLAMDOrdering<int> > solver;
            solver.compute(hess.selfadjointView<Eigen::Upper>());
            dir = solver.solve(grad);
            if(verbose) std::cout << "SparseQR\n";
        }
        
        ScalarType delta=1;
        ScalarType val = Infinity;
        
        if(dampingStrategy==DampingType::Increasing){
            
            val = SafeValue(pb,x-dir);
            while(delta >= sDelta.stopBelow){
                delta/=2;
				
                const ScalarType cmp = SafeValue(pb,x-delta*dir);
                
                if(runtimeOut)
                    (*runtimeOut) << " Tried " ExportVarArrow(delta)
					<< " got " ExportVarArrow(cmp) << "\n";
				
                if(cmp!=Infinity && cmp > val){
                    delta *= 2;
                    x = x-delta*dir;
                    val = pb.Value(x); // No safety net here (value already tried)
                    break;
				} else {
					val = cmp;}
            }
        } else { // Valid
            
            while(delta >= sDelta.stopBelow &&
                  (val = SafeValue(pb,x-delta*dir) ) > sObjective.values.back()){
                if(runtimeOut)
                    (*runtimeOut) << " Tried " ExportVarArrow(delta) << " got " ExportVarArrow(val) << "\n";
                delta *= dampingRatio;
            }
            x= x-delta*dir;
        }
        
        
        
        // Handling stopping criteria
        sDelta.Insert(delta);
        sObjective.Insert(val);
        const ScalarType ghNorm = grad.dot(dir);
        if(! sGHNorm.values.empty() )
            sGHDecay.Insert(ghNorm/sGHNorm.values.back());
        sGHNorm.Insert( ghNorm );
        
        sGNorm.Insert( grad.norm() );
        sDirNorm.Insert( dir.norm() );
        
        if(runtimeOut){
            (*runtimeOut) << "Newton unconstrained : "
            ExportVarArrow(iter)
            ExportVarArrow(delta)
            ExportVarArrow(val)
            ExportVarArrow(ghNorm)
            << "\n";
        }
        
        // Check if should stop.
        bool q=true;
        if(sObjective.Abort(q) || sGHNorm.Abort(q) || sGHDecay.Abort(q)
           || sDelta.Abort(q) || sGNorm.Abort(q) || sDirNorm.Abort(q)
           || q) {
/*            std::cout << "Aborting "
            ExportVarArrow(q)
            << "\n";*/
            break;}
    }
}

void NewtonUnconstrained::Clear(){
    iter=0;
    sObjective.Clear(); sGHNorm.Clear(); sGHDecay.Clear();
    sDelta.Clear(); sGNorm.Clear(); sDirNorm.Clear();
}

void NewtonUnconstrained::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(maxIter)
    ExportVarArrow(iter)
    ExportVarArrow(dampingRatio) << "\n"
    ExportVarArrow(sObjective) << "\n"
    ExportVarArrow(sGHNorm) << "\n"
    ExportVarArrow(sGHDecay) << "\n"
    ExportVarArrow(sDelta) << "\n"
    ExportVarArrow(sGNorm) << "\n"
    ExportVarArrow(sDirNorm) << "\n"
    << "}";
}

// %%%%%%%%%%%%%%%%%% Newton Constrained %%%%%%%%%%%%%%%%%%%

ScalarType NewtonConstrained::Value(const VectorType & x){
    if(!verbose){
        value_=objective->Value(x);
        for(auto c : barriers)
            value_ += multiplier * c->Value(x);
    } else {
		std::ostream & os = std::cout;
//		os << "Operator evaluation at " ExportArrayArrow(StdFromEigen(x));
        ScalarType v = objective->Value(x);
        value_=v; os << "Objective value : " << v;
        for(auto c : barriers) {
            v=c->Value(x);
            value_ += multiplier * v;
            os << ", Constraint value : " << v;}
        os << "\n";
    }
    return value_;
}

const VectorType & NewtonConstrained::Gradient(){
    gradient_ = objective->Gradient();
    for(auto c : barriers)
        gradient_ += multiplier * c->Gradient();
    return gradient_;
}

const SparseMatrixType & NewtonConstrained::Hessian(){
    hessian_ = objective->Hessian();
    for(auto c : barriers)
        hessian_ += multiplier * c->Hessian();
    return hessian_;
}

void NewtonConstrained::Solve(Functionnal & objective_,
                              VectorType & x,
                              std::vector<Functionnal*> barriers_){
    objective=&objective_;
    barriers=barriers_;
    
    while (multiplier>=multiplierBound) {
        if(runtimeOut){
            (*runtimeOut) << "\nNewton constrained, "
            ExportVarArrow(multiplier)
//            ExportVarArrow(opt)
            << " beginning unconstrained inner optimization\n";
        }

        opt.Clear();
        opt.sGHNorm.stopBelow = multiplier*ghNormBase;
        const VectorType x_save = x;
        opt.Solve(*this,x);
        
        iter += opt.iter;
        multiplier*=multiplierDamping;
        
        if(opt.sObjective.values.back()==Infinity){
            std::cout << "Aborting optimization (infinite objective).\n";
            x=x_save;
            break;}
        else if(!opt.sGHNorm.active){
            std::cout << "Aborting optimization (sGHNorm) .\n";
            break;}
        else if(iter>maxIter){
            std::cout << "Aborting (maxIter).\n";
            break;}
    }
    
    std::cout << "Finishing optimization.\n";
}

void NewtonConstrained::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(maxIter)
    ExportVarArrow(iter)
    ExportVarArrow(multiplier)
    ExportVarArrow(multiplierBound)
    ExportVarArrow(multiplierDamping)
    ExportVarArrow(ghNormBase)
    ExportVarArrow(opt)
    << "}";
    
}
