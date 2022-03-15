//
//  QuotientedNewton.hxx
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 02/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef Minkowski_QuotientedNewton_hxx
#define Minkowski_QuotientedNewton_hxx

bool
NewtonScheme::TestType::operator()(ScalarType val){
    values.push_back(val);
    return Active();
}

bool NewtonScheme::TestType::Active() const {
    if(values.size()<delay) return false;
    return *::std::max_element(values.end()-delay, values.end()) <  lowerBound;
}

void
NewtonScheme::Clear(){
    stepSizeTest.values.clear();
    residueDecayTest.values.clear();
    residueMinTest.values.clear();
    timings.Clear();
    symmetryDefects.clear();
    numberOfEffectiveIterations=0;
}

void NewtonScheme::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(delta)
    ExportVarArrow(deltaGranularity)
    ExportVarArrow(deltaMin)
    ExportVarArrow(forceResidueDecrease)
    ExportVarArrow(maxIterations)
    ExportVarArrow(stepSizeTest)
    ExportVarArrow(residueDecayTest)
    ExportVarArrow(residueMinTest)
    ExportVarArrow(numberOfEffectiveIterations)
    ExportVarArrow(timings)
    ExportArrayArrow(symmetryDefects)
    << "}";
}

void NewtonScheme::TestType::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(lowerBound)
    ExportVarArrow(delay)
    ExportArrayArrow(values)
    << '"' << "Active" << '"' << "->" << Active()
    << "}";
}

NewtonScheme::StoppingCriterion
NewtonScheme::Solve(FunctionalBase & op,
                    const VectorType & target,
                    VectorType & x)
{
    assert(target.size()==x.size());
    assert(!target.hasNaN());
    assert(!x.hasNaN());
    Clear();

    VectorType diff;
    
    const auto s = x.size();
    VectorType last_x;      last_x.resize(s); // assign x ??
    VectorType next_x;      next_x.resize(s);
    VectorType last_diff;   last_diff.resize(s);
    VectorType next_diff;   next_diff.resize(s);
    
    for(numberOfEffectiveIterations=0; numberOfEffectiveIterations<maxIterations; ++numberOfEffectiveIterations){

        if( op.Compute(x, diff) ) return StoppingCriterion::EvaluationError;
        assert(!diff.hasNaN());
        diff-=target;
        
        const ScalarType residue = diff.norm();
        if( residueMinTest(residue) ) return StoppingCriterion::ResidueMinTest;
        
        const ScalarType residueDecay = (numberOfEffectiveIterations==0) ? 1. :
        1. - residue / *(residueMinTest.values.end()-2);
        if( residueDecayTest(residueDecay) ) return StoppingCriterion::ResidueDecayTest;
        
        std::cout
        << "residue : " << residue
        << ", residue decay : " << residueDecayTest.values.back()
        << std::endl;
        
        // Compute the descent direction
        if( op.DescentDirection(diff) ) return StoppingCriterion::InversionError;
        const VectorType & direction = diff;
        assert(!direction.hasNaN());
        std::cout << "direction found" << std::endl;
        
        ScalarType last_delta = delta*deltaGranularity;
        ScalarType last_residue = std::numeric_limits<ScalarType>::infinity();
        
        while(true){
            ScalarType next_delta = last_delta/deltaGranularity, next_residue;
            
            next_x=x+next_delta*direction;
            if( op.Compute(next_x,next_diff) ) next_residue = Infinity;
            else { next_diff -= target; next_residue=next_diff.norm();}
            
            std::cout
            <<"delta trial : " << next_delta
            << ", residue trial : " << next_residue
            << "\n";
            
            if(next_residue < Infinity &&
               next_residue >= last_residue &&
               (!forceResidueDecrease || last_residue <= residue))
                break;
            
            std::swap(last_delta, next_delta);
            std::swap(last_residue, next_residue);
            last_x.swap(next_x);
            last_diff.swap(next_diff);
            
            if(last_delta < deltaMin || deltaGranularity <= 1) break;
        }
        
        last_x.swap(x);
        //        std::swap(last_diff, diff); // recomputed anyway in order to get the gradient
        if( stepSizeTest(last_delta) ) return StoppingCriterion::StepSizeTest;
        if(last_residue==std::numeric_limits<ScalarType>::infinity()) return StoppingCriterion::InfiniteResidue;
        std::cout << "next iter\n\n";
    }
    return StoppingCriterion::MaxIterations;
}


// Timings

void NewtonScheme::TimingsType::Clear(){apply.clear();applyDifferentiate.clear();solve.clear();}
void NewtonScheme::TimingsType::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportArrayArrow(solve)
    ExportArrayArrow(apply)
    ExportArrayArrow(applyDifferentiate)
    << "}";
}


#endif
