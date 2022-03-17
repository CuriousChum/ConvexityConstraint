//
//  Meissner_3_Test.cpp
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 04/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#include <fstream>
#include "Headers/Meissner_3.h"
#include "nlopt.hpp"


namespace Meissner_3_Test {
    using namespace Meissner_3;
    typedef N::VectorType NV;
    
    
    void OptimizeMeissner3(int subdivision, ScalarType stopval = 0.4197){
		std::cout << "Numerical optimization of Meissner problem "
		<< " using point set number " << subdivision << std::endl;
		
		std::string filename = "HalfPts_"+std::to_string(subdivision)+".txt";
        std::ifstream is;
        is.open(filename);
        
		if(!is.is_open()){
			std::cout << "Error !! Could not open file " << filename << ".\n"
			<< "This file is likely in the ConvexityConstraint/Minkowski/data folder of the code repo. (Otherwise, it can be generated with the ConvexityConstraint/Minkowski/HalfPts.nb notebook). The file must be copied to the same folder as the compiled executable.\n";
			return;
		}
		
        int n,checkDim;
        is >> n;
        is >> checkDim;
        assert(checkDim==Dimension);
        
        
        Minkowski_3::Functional mink;
        G::VectorList & pts = mink.pts;
        pts.resize(Dimension,2*n);
        for(int i=0; i<n; ++i)
            for(int j=0; j<Dimension; ++j)
                is >> pts(j,i);
            
        
        for(int i=0; i<n; ++i)
            pts.col(n+i) = -pts.col(i);
        
        {
            std::ostream & os = std::cout
            ExportVarArrow(n)
            ExportArrayArrow(N::StdFromEigen(NV(pts.col(0))))
            << "\n";
        }
        
		/* Initialization of the optimization.
		 The initialization must be such that all the facets are non-empty.
		 Depending on the initialization, the optimization may be attracted by one or the other of Meissner's (conjecturedly optimal) bodies.
		 */
		
		// Modifying this parameter (e.g. flipping the sign) is sometimes enough so that
		// the numerica optimizer is attracted to one or the other of Meissner bodies.
		const ScalarType delta=-0.02;
		
        Eigen::Matrix<ScalarType,4,3> tetra;
        tetra <<
        0.,0.,-1.,
        -0.942809,0.,0.333333,
        0.471405,-0.816497,0.333333,
        0.471405,0.816497,0.333333;
        
        
        std::vector<ScalarType> x(n,0.5),init;
        for(int i=0; i<n; ++i){
            for(int j=0; j<4; ++j){
                ScalarType s = pts.col(i).dot(tetra.row(j));
                x[i] += delta*(s*s*s);
            }
        }
        init=x;
        
        
        {
            NV input(2*n), output(2*n);
            for(int i=0; i<n; ++i){
                input[i]=x[i];
                input[n+i]=1-x[i];
                }
            const int result = mink.Compute(input,output);
            std::ostream & os = std::cout
            ExportVarArrow(result)
            ExportArrayArrow(N::StdFromEigen(input))
            ExportArrayArrow(N::StdFromEigen(output))
            << "\n";
        }
        
        
        nlopt::opt algo=nlopt::opt(nlopt::LD_SLSQP, n); //LD_MMA
        algo.set_min_objective(EvaluateVolume, static_cast<void*>(&mink));
        algo.add_inequality_constraint(EvaluateConstraint, static_cast<void*>(&mink));
        algo.set_lower_bounds(0);
        algo.set_upper_bounds(1);
        algo.set_maxeval(10000);
        
        algo.set_stopval(stopval);
        
        ScalarType minimum;
        std::cout << "Opt begin\n";
		nlopt::result result;
		try {result = algo.optimize(x,minimum);}
		catch(nlopt::roundoff_limited & e){
			std::cout << "Caught roundoff_limited exception\n";}
        std::cout << "Opt end\n";
        
        std::ofstream os; os.open("Meissner_3_"+std::to_string(subdivision)+".txt");
        os << "{"
        ExportVarArrow(minimum)
        ExportVarArrow(result)
        ExportArrayArrow(x)
        ExportArrayArrow(init)
        ExportArrayArrow(N::StdFromEigen(NV(pts.row(0))))
        ExportArrayArrow(N::StdFromEigen(NV(pts.row(1))))
        ExportArrayArrow(N::StdFromEigen(NV(pts.row(2))))
        << "}";

    }
}

int main(int argc, const char * argv[]) {
    
        
	int i=4;
	if(argc>1) i=atoi(argv[1]);
	
	double stopVal = 0.4197;
	if(argc>2) stopVal = atof(argv[2]);
	
	if(i>0) Meissner_3_Test::OptimizeMeissner3(i,stopVal);
	else {for(i=2;i<=7;++i) Meissner_3_Test::OptimizeMeissner3(i,stopVal);}
	
    return 0;
}
