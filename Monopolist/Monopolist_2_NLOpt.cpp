//
//  PrincipalAgent_NLOpt.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 07/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_PrincipalAgent_NLOpt_h
#define CGalTest_PrincipalAgent_NLOpt_h


#include "PrincipalAgent_Test.h"
#include "nlopt.hpp"

namespace PrincipalAgent_Test {
    
    
    void PA_NLOpt_0(int n=30, ShapeType shape=ShapeType::Square){
        std::ofstream os;
        
        ConvexityConstraint cvx(MakeShape(n,shape));
        std::vector<ScalarType> x;
        for(auto p : cvx.GetPts()) x.push_back(CGT::Parabola(p.first));
        
        cvx.SetValues(x);
        cvx.Compute(31);
        std::cout ExportVarArrow(cvx.error) ExportVarArrow(cvx.logSum) << "\n";
        
        os.open("cvxTest.txt");
        os << "{" ExportVarArrow(cvx) << "}";
        os.close();
        
        
        PrincipalAgent pa;
        pa.rt = cvx.rt; // Hope this is a true deep copy.
        
        
        nlopt::opt algo=nlopt::opt(nlopt::LD_MMA , (int)x.size()); //LD_SLSQP, LD_MMA
        algo.set_min_objective(PrincipalAgent::Evaluate, static_cast<void*>(&pa));
        algo.add_inequality_constraint(ConstraintType::GeometricMean, static_cast<void*>(&cvx), 0.);
        //            algo.add_inequality_constraint(ConstraintType::GeometricMean, static_cast<void*>(&lip));
        algo.set_lower_bounds(0);
        algo.set_maxeval(1000);

        ScalarType minimum;
        nlopt::result result = algo.optimize(x,minimum);

        cvx.Compute(31);
        
        const std::string filename =
        "PA_NLOpt_0_"+enumToRealString(shape)+".txt";
        os.open(filename);
        os << "{"
        ExportVarArrow(pa)
        ExportVarArrow(cvx)
        ExportArrayArrow(x)
        ExportVarArrow(result)
        ExportVarArrow(minimum)
        << "}";
    }

    
    
    void PALinear2(std::vector<CGT::Full_point> pts, int nRefine, std::string file = "PALinear2"){
        PrincipalAgent pa;
        pa.rt = CGT::RT(pts.begin(),pts.end());
        pa.rt.infinite_vertex()->info() = {InfiniteIndex, -1};
        std::ofstream os;
        //        pa.hasQuadraticCost=false;
        
        for(int iRefine = 0; iRefine < nRefine; ++iRefine){
            const size_t n=pts.size();
            ConvexityConstraint cvx(pts);
            LipschitzConstraint lip(pts);
            lip.lastOfLineOnly = false;

            nlopt::opt algo=nlopt::opt(nlopt::LD_MMA , (int)n); //LD_SLSQP
            algo.set_min_objective(PrincipalAgent::Evaluate, static_cast<void*>(&pa));
            algo.add_inequality_constraint(ConstraintType::GeometricMean, static_cast<void*>(&cvx), 0);
//            algo.add_inequality_constraint(ConstraintType::GeometricMean, static_cast<void*>(&lip));
            algo.set_lower_bounds(0);
            algo.set_maxeval(50);
            
            ScalarType minimum;
            std::vector<ScalarType> x(n),init;
            for(int i=0; i<n; ++i) x[i]=CGT::Parabola(pts[i].first)/10.;
            init=x;
            
            {
                cvx.SetValues(x);
                lip.SetValues(x);
                cvx.Compute(31);
                lip.Compute(31);
                os.open("cvxTest.txt");
                os << "{" ExportVarArrow(cvx) ExportVarArrow(lip) << "}";
                os.close();
            }
            
            std::cout << "Opt begin\n";
            nlopt::result result = algo.optimize(x,minimum);
            std::cout << "Opt end\n";

            {
                cvx.SetValues(x);
                lip.SetValues(x);
                cvx.Compute(31);
                lip.Compute(31);
                
                const std::string filename = file+"_iter"+std::to_string(iRefine)+".txt";
                os.open(filename);
                os << "{"
                ExportVarArrow(pa)
                ExportVarArrow(cvx)
                ExportVarArrow(lip)
                ExportArrayArrow(x)
                ExportVarArrow(result)
                ExportVarArrow(minimum)
                << "}";
                os.close();
            }
            
            std::cout << "Points before refinement : " ExportVarArrow(pa.Pts().size());
            const ScalarType refinementQuantile = 0.3;
            pa.Refine(refinementQuantile,NS::EigenFromStd(x));
            pts = pa.Pts();
            std::cout ExportVarArrow(pts.size()) << "\n";
        }
    }

}


int main(int argc, const char * argv[]){

	PrincipalAgent_Test::PALinear2(PrincipalAgent_Test::MakeShape(10,PrincipalAgent_Test::ShapeType::Square),
								  1, "PARef2_Square_Linear");
	
    PrincipalAgent_Test::PA_NLOpt_0(10);
}


#endif
