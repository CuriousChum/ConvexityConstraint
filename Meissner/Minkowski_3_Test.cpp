//
//  Minkowski_3_Test.cpp
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 03/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#include <fstream>
#include "Headers/Minkowski_3.h"

/**
This file contains some undocumented validation tests for the Minkowski problem in three dimensions :
find a convex body whose facets areas have a given distribution.
 
Please contact the author for details in case of need.
*/


namespace Minkowski_3_Test {
    using namespace Minkowski_3;
    std::ostream & os = std::cout;
    typedef N::VectorType NV;
    
    typedef std::vector<std::array<ScalarType,4> > EqnsType;
    EqnsType
    cubeEqns =
    {{0.,0.,1.,0.5},{1.,0.,0.,0.5},{0.,1.,0.,0.5},{-1.,0.,0.,0.5},{0.,0.,-1.,0.5},{0.,-1.,0.,0.5}},
    octahedronEqns =
    {{0.57735,-0.57735,0.57735,0.408248},{0.57735,0.57735,0.57735,0.408248},{-0.57735,0.57735,0.57735,0.408248},{-0.57735,-0.57735,0.57735,0.408248},{-0.57735,-0.57735,-0.57735,0.408248},{0.57735,-0.57735,-0.57735,0.408248},{-0.57735,0.57735,-0.57735,0.408248},{0.57735,0.57735,-0.57735,0.408248}},
    tetraEqns =
    {{0.,0.,-1.,0.204124},{-0.942809,0.,0.333333,0.204124},{0.471405,-0.816497,0.333333,0.204124},{0.471405,0.816497,0.333333,0.204124}},
    dodecaEqns =
    {{-0.894427,0.,-0.447214,1.11352},{0.894427,0,0.447214,1.11352},{0.276393,-0.850651,0.447214,1.11352},{0.,0.,1.,1.11352},{0.276393,0.850651,0.447214,1.11352},{0.723607,0.525731,-0.447214,1.11352},{0.723607,-0.525731,-0.447214,1.11352},{-0.276393,0.850651,-0.447214,1.11352},{0.,0.,-1.,1.11352},{-0.276393,-0.850651,-0.447214,1.11352},{-0.723607,-0.525731,0.447214,1.11352},{-0.723607,0.525731,0.447214,1.11352}};
    
    
    // Check volumes, barycenters
    void Test1(){
        const auto & eqns = dodecaEqns;
        const int n=(int)eqns.size();
        Functional op;
        G::VectorList & pts = op.pts;
        pts.resize(Dimension,n);
        N::VectorType input(n), target(n);
        
        for(int i=0; i<n; ++i){
            for(int j=0; j<Dimension; ++j)
                pts(j,i) = eqns[i][j];
            input[i]=eqns[i][Dimension];
        }
        
        G::Vector bary{0.1,0.2,0.3};
        input+=pts.transpose()*bary;
        
        
        op.Compute(input, target);
        os 
        ExportVarArrow(op.Volume())
        ExportArrayArrow(N::StdFromEigen(op.Barycenter()))
        ExportArrayArrow(N::StdFromEigen(target))
        << "\n";
    }
    
    struct {
        ScalarType operator () (const G::Vector & v) const {
            const ScalarType x=v[0], y=v[1], z=v[2];
            return sqrt(x*x+2*y*y-2*x*y+3*z*z) + std::max(3*x+2*y+z,0.)/3.;
            //            return std::max(x+2*y,0.)+sqrt(x*x+y*y);
        }
    } norm;
    
    
    void Test2(){
        const auto & eqns = dodecaEqns;
        const int n=(int)eqns.size();
        Functional op;
        G::VectorList & pts = op.pts;
        pts.resize(Dimension,n);
        NV sol(n), target(n);
        
        for(int i=0; i<n; ++i){
            for(int j=0; j<Dimension; ++j)
                pts(j,i) = eqns[i][j];
            sol[i]=norm(pts.col(i));
        }
        op.referenceIndices={3,7,9};
//        op.CheckData();
        op.Compute(sol, target);
                
        N::NewtonScheme solver;
        N::VectorType x = N::VectorType::Constant(n, 1);
        
        solver.Solve(op,target,x);
        
        std::ofstream os;
        os.open("Minkowski_3_Test2.txt");
        
        os << "{"
        ExportVarArrow((sol-x).norm())
        ExportArrayArrow(N::StdFromEigen(sol))
//        ExportArrayArrow(N::StdFromEigen(x))
        ExportArrayArrow(N::StdFromEigen(target))
        << "\n"
        ExportArrayArrow(N::StdFromEigen(NV(pts.row(0))))
        ExportArrayArrow(N::StdFromEigen(NV(pts.row(1))))
        ExportArrayArrow(N::StdFromEigen(NV(pts.row(2))))
        << "}";
    }
    
    // Assumes that vertex normals are distinct
    void MinkowskiInterpolate(std::array<EqnsType*,2> eqns, int steps=5, std::string fileName="InterpolationDefault.txt"){
        Functional ops[2];
        N::VectorType targets[2];
        
        for(int e=0; e<=1; ++e){
            Functional & op = ops[e];
            G::VectorList & pts = op.pts;
            EqnsType eqn = *eqns[e];
            
            pts.resize(Dimension,eqn.size());
            const int n=op.Size();
            N::VectorType x(n);
            for(int i=0; i<n; ++i){
                for(int j=0; j<Dimension; ++j){
                    pts(j,i) = eqn[i][j];
                    x(i)=eqn[i][Dimension];
                }
            }
            op.Compute(x,targets[e]);
        }
        
        //const int n = op[0].Size()+op[1].Size();
        Functional op;
        G::VectorList & pts = op.pts;
        pts.resize(Dimension,ops[0].Size()+ops[1].Size());
        pts << ops[0].pts , ops[1].pts;
        
        std::ofstream os;
        os.open(fileName);
        os << "{"
        ExportVarArrow(steps)
        ExportArrayArrow(N::StdFromEigen(NV(pts.row(0))))
        ExportArrayArrow(N::StdFromEigen(NV(pts.row(1))))
        ExportArrayArrow(N::StdFromEigen(NV(pts.row(2))));

        
        const int n=op.Size();
        N::VectorType target(n), x(n);
        N::NewtonScheme solver;
        x.setConstant(1);
        std::cout << op.Compute(x,target) << op.pts;
        
        for(ScalarType i=1; i<steps; ++i){
            target << (i/steps)*targets[0], (1-i/steps)*targets[1];
//            std::cout << target << "\n\n";
            x.setConstant(1);
            solver.Solve(op,target,x);
			os << i << "->{"
			ExportArrayArrow(N::StdFromEigen(x))
			<< "},";
			/*
            os << i << "->";
            print_container(os,N::StdFromEigen(x));
            os << ",";*/

        }
        os << "}";
    }
    
}

int main(int argc, const char * argv[]) {

	std::cout << "Use --help as first argument to learn about usage.\n";
	if(argc>=2 && !strcmp(argv[1],"--help")){
		std::cout <<
		"This executable is currently untested and in development.\n"
		;
	};

    Minkowski_3_Test::Test1();
    
//    Minkowski_3_Test::MinkowskiInterpolate({&Minkowski_3_Test::octahedronEqns, &Minkowski_3_Test::tetraEqns},30);
}
