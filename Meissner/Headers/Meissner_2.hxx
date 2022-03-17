//
//  Meissner_2.hxx
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 03/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef Minkowski_Meissner_2_hxx
#define Minkowski_Meissner_2_hxx


/*
Functional::Functional(const G::VectorList & halfPts){
    G::VectorList & pts = mink.pts;
    pts.resize(Dimension, 2*halfPts.size());
    for(int i=0; i<Size(); ++i)
        pts.col(i) = halfPts.col(i);
    for(int i=0; i<Size(); ++i)
        pts.col(i) = halfPts.col(Size()+i);
}*/

/*
int Functional::Compute(N::VectorType & input, N::VectorType & output){
    const int n=Size();
    N::VectorType minkInput(2*n), minkOutput(2*n);
    minkInput << input, (Eigen::Constant(n,1.)-input);
    
    const int error = mink.Compute(minkInput,minkOutput); if(error) return error;
    input = minkInput.head(n); // mink.Compute adjusts weights to centers the polytope
    output = minkOutput.head(n) - minkOutput.tail(n);
    
    // Now add the penalty for convexity
    
}*/

/**
 Evaluates the area of the two dimensional convex polyhedron
 H(a) with a = [x,1-x],
 where brackets denote vector concatenation.
 This polyhedron has width one, provided 0<= x <= 1, and the points involved in the polyhedron definition have unit norm
 and obey pts[n+i] = -pts[i]
 
 See the class Minkowski_2::Functional for the definition of H(a). (An instance  is passed as an opaque third argument.)
 
 Input :
  - x, vector defining the convex polyhedron shape.
  - data, an appropriate Minkowski_2::Functional  instance.
 Output :
  - grad, the gradient of Area( H( [x,1-x] ) ) with respect to x.
 Return value : Area( H( [x,1-x] ) )
 */
ScalarType EvaluateVolume(const std::vector<ScalarType> & x,
						  std::vector<ScalarType> & grad, void * data){
    Minkowski_2::Functional & mink = *static_cast<Minkowski_2::Functional*>(data);
    const int n=mink.Size()/2;
    N::VectorType input = N::EigenFromStd(x);
    N::VectorType minkInput(2*n), minkOutput(2*n);
    minkInput << input, (N::VectorType::Constant(n,1.)-input);
    
    const int error = mink.Compute(minkInput,minkOutput); if(error) return N::Infinity;
    if(!grad.empty())
        grad = N::StdFromEigen(N::VectorType( minkOutput.head(n) - minkOutput.tail(n) ));
    
    static int counter=0;
//    std::ofstream os; os.open("EvaluateVolume_"+std::to_string(++counter)+".txt");
    std::cout << "EvaluateVolume "
    ExportVarArrow(mink.Area())
    ExportVarArrow(grad.empty())
    ExportVarArrow(++counter)
    << "\n\n";
    
    return mink.Area();
}

/**
 Evaluates the opposite of the geometric mean of the edge lengths of the two dimensional convex polyhedron
 H( [x,1-x] ),
 see the definition of EvaluateVolume for details. Also returns the gradient.
 
 We thus encode a large number of (strict) positivity constraints into a single (strict) positivity constraint, using the geometric
 mean. Note that this may be regarded Bad Practice from the optimization standpoint, but it simplifies the interface, and the
 nlopt optimizer appears to be robust enough to handle it. 
 */

ScalarType EvaluateConstraint(const std::vector<ScalarType> & x, std::vector<ScalarType> & grad, void * data){
    Minkowski_2::Functional & mink = *static_cast<Minkowski_2::Functional*>(data);
    const int n=mink.Size()/2;
    N::VectorType input = N::EigenFromStd(x);
    N::VectorType minkInput(2*n), minkOutput(2*n);
    minkInput << input, (N::VectorType::Constant(n,1.)-input);
    
    const int error = mink.Compute(minkInput,minkOutput); if(error) return 0;
    
    // Note the minus sign: constraints in nlopt must be <= 0
    const ScalarType constraintValue = - N::ConstraintsProduct(mink.FacetLengths(), mink.Jacobian(), minkOutput);
    
    if(!grad.empty())
        grad = N::StdFromEigen(N::VectorType( -(minkOutput.head(n) - minkOutput.tail(n)) ));
    
    std::ostream & os = std::cout << "EvaluateConstraint "
    ExportVarArrow(constraintValue)
//    ExportArrayArrow(N::StdFromEigen(mink.FacetLengths()))
    << "\n\n";
    
    return constraintValue;

}

#endif
