//
//  Meissner_3.hxx
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 04/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef Minkowski_Meissner_3_hxx
#define Minkowski_Meissner_3_hxx


ScalarType EvaluateVolume(const std::vector<ScalarType> & x, std::vector<ScalarType> & grad, void * data){
    Minkowski_3::Functional & mink = *static_cast<Minkowski_3::Functional*>(data);
    const int n=mink.Size()/2;
    N::VectorType input = N::EigenFromStd(x);
    N::VectorType minkInput(2*n), minkOutput(2*n);
    minkInput << input, (N::VectorType::Constant(n,1.)-input);
    
    const int error = mink.Compute(minkInput,minkOutput); if(error) return N::Infinity;
    if(mink.FacetAreas().minCoeff()<=0.) return N::Infinity;
    if(!grad.empty())
        grad = N::StdFromEigen(N::VectorType( minkOutput.head(n) - minkOutput.tail(n) ));
    
    static int counter=0;
    //    std::ofstream os; os.open("EvaluateVolume_"+std::to_string(++counter)+".txt");
    std::cout << "EvaluateVolume "
    ExportVarArrow(mink.Volume())
    ExportVarArrow(grad.empty())
    ExportVarArrow(++counter)
    << "\n\n";
    
    return mink.Volume();
}

ScalarType EvaluateConstraint(const std::vector<ScalarType> & x, std::vector<ScalarType> & grad, void * data){
    Minkowski_3::Functional & mink = *static_cast<Minkowski_3::Functional*>(data);
    const int n=mink.Size()/2;
    N::VectorType input = N::EigenFromStd(x);
    N::VectorType minkInput(2*n), minkOutput(2*n);
    minkInput << input, (N::VectorType::Constant(n,1.)-input);
    
    const int error = mink.Compute(minkInput,minkOutput); if(error) return 0;
    
    // Note the minus sign: constraints in nlopt must be <= 0
    const ScalarType constraintValue = - N::ConstraintsProduct(mink.FacetAreas(), mink.Jacobian(), minkOutput);
    
    if(!grad.empty())
        grad = N::StdFromEigen(N::VectorType( -(minkOutput.head(n) - minkOutput.tail(n)) ));
    
    std::ostream & os = std::cout << "EvaluateConstraint "
    ExportVarArrow(constraintValue)
    //    ExportArrayArrow(N::StdFromEigen(mink.FacetLengths()))
    << "\n\n";
    
    return constraintValue;
    
}



#endif
