//
//  PrincipalAgent_Test.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 11/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_PrincipalAgent_Test_h
#define CGalTest_PrincipalAgent_Test_h

#include <fstream>

#include "JMM_CPPLibs/Output/EnumToString.h"

#include "PrincipalAgent.h"
#include "ConvexityConstraint.h"
#include "LipschitzConstraint.h"

namespace PrincipalAgent_Test {enum class ShapeType {Square, Triangle, Circle}; }
template<> char const* enumStrings<PrincipalAgent_Test::ShapeType>::data[] = {"Square", "Triangle", "Circle"};

namespace PrincipalAgent_Test {
    using namespace Geometry_2;

    std::vector<CGT::Full_point> MakeShape(int n, ShapeType shape){
        const ScalarType h=1./n;
        const ScalarType pi = 4.*atan(1.);
        std::vector<CGT::Full_point> pts;
        int counter=0;
        
        typedef CGT::Weighted_point WP;
        typedef CGT::InfoType IT;
        
        switch (shape) {
            case ShapeType::Square:
                for(int i=0; i<n; ++i)
                    for(int j=0; j<n; ++j){
                        IndexType flag =
                        (i==0 ? 1 : 0) | (j==0 ? 2 : 0) |
                        (i==n-1 ? 4 : 0) | (j==n-1 ? 8 : 0);
                        pts.push_back({WP({1.+(i+0.5)*h,1.+(j+0.5)*h},0),
                            IT(counter++,flag)});
                    }
                
                break;
                
            case ShapeType::Circle:{
                const ScalarType c0=1.5, r0=0.5;
                pts.push_back({WP({c0,c0},0),IT(counter++,0)});
                for(int i=2; i<=n; i+=2){
                    const ScalarType r = i*r0*h;
                    const ScalarType k = (int)floor(2.*pi*r/h);
                    for(int j=0; j<k; ++j){
                        const ScalarType t = (2.*pi*j)/k;
                        pts.push_back({WP({c0+r*cos(t), c0+r*sin(t)},0), IT(counter++,0)});
                    }
                }
                break;}
                
            case ShapeType::Triangle:{
                const ScalarType c0=1.5;
                int counter=0;
                for(int i=0; i<n; ++i)
                    for(int j=0; j<n-i; ++j){
                        const ScalarType c=cos(pi/12), s=sin(pi/12);
                        
                        IndexType flag = 0;
                        if(i==0) flag = flag | 1;
                        if(j==0) flag = flag | 1<<1;
                        if(i+j==n-1) flag = flag | 1<<2;
                        
                        pts.push_back({WP({c0-h*(i*c+j*s),c0-h*(i*s+j*c)},0),IT(counter++,flag)});

                    }
                break;}
        }

        return pts;
    }
    
    // %%%%%%%%%%%%%%%% Standard Principal Agent %%%%%%%%%%%%%%%%%
    
    
    void PA0(int n=30, ShapeType shape=ShapeType::Square){
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
        
        
        PositivityConstraint pos;
        
        std::vector<NS::Functionnal*> constraints = {&cvx, &pos}; //
        
        PrincipalAgent pa;
        pa.rt = cvx.rt; // Hope this is a true deep copy.
        
        NS::NewtonConstrained newton;
        newton.maxIter=50;
        
        /*
        newton.multiplier = 2.e-5;
        newton.opt.sDelta.stopBelow = 1e-4;
         */
        
        newton.multiplier = 1./x.size();
        newton.multiplierBound *= newton.multiplier; // Vanishing constraint penalisation
//        newton.multiplierBound = newton.multiplier; // Fixed penality

        NS::VectorType xx(x.size());
        for(int i=0; i<x.size(); ++i) xx[i]=x[i];
                
        newton.Solve(pa,xx,constraints);
        
        cvx.Compute(31);
        pos.Compute(31);
        
        const std::string filename =
        "PA0_"+enumToRealString(shape)+".txt";
        os.open(filename);
        os << "{"
        ExportVarArrow(newton)
        ExportVarArrow(pa)
        ExportVarArrow(cvx)
        ExportVarArrow(pos)
        ExportArrayArrow(x)
        << "}";
        
        std::cout
        ExportVarArrow(newton)
        << "\n";
    }
    
    
    // %%%%%%%%%%%%%%%% A basic opt test %%%%%%%%%%%%%%%%
    
    struct TestFunction : ConstraintType {
        virtual void SetValues(const std::vector<ScalarType> & x_){
            x = x_[0];
            error=0;
        }
        
        virtual void Compute(FlagType r) {
            if(r & RLogSum) logSum = 0.5*x*x;
            if(r & RLogGrad) {logGrad.resize(1); logGrad[0] = x;}
            if(r & RLogHessian){logHessian.resize(1); logHessian[0] = {0,0,1};}
        };
        ScalarType x;
    };
    
    void Opt0(){
        const int n=1;
        TestFunction pb;
        PositivityConstraint pos;
        std::vector<NS::Functionnal*> constraints = {&pos};
        NS::NewtonConstrained newton;
        newton.multiplier=1;
        newton.multiplierBound=1e-4; // Only one step
        
        NS::VectorType xx(n);
        xx[0]=2;
        newton.Solve(pb,xx,constraints);
        
        std::ofstream os;
        os.open("Opt0.txt");
        os << "{"
        ExportVarArrow(newton)
        ExportVarArrow(pos)
        ExportVarArrow(xx[0])
        << "}";
    }

    
    
    
    // %%%%%%%%%%%%%%%%%%%% With refinement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    void PA1(std::vector<CGT::Full_point> pts, int nRefine, std::string file = "PA1"){
        PrincipalAgent pa;
        pa.rt = CGT::RT(pts.begin(),pts.end());
        pa.rt.infinite_vertex()->info() = {InfiniteIndex, -1};
        std::ofstream os;
        
        for(int iRefine = 0; iRefine < nRefine; ++iRefine){
            
            ConvexityConstraint cvx(pts);
            PositivityConstraint pos;
            std::vector<NS::Functionnal*> constraints = {&cvx, &pos}; //
            
            NS::NewtonConstrained newton;
            newton.maxIter=50;

            newton.multiplier = 1./pts.size();
            newton.multiplierBound *= newton.multiplier;
            
            NS::VectorType x(pts.size());
            for(int i=0; i<x.size(); ++i) x[i]=CGT::Parabola(pts[i].first);
            
            cvx.SetValues(NS::StdFromEigen(x));
            cvx.Compute(31);
            os.open("cvxTest.txt");
            os << "{" ExportVarArrow(cvx) << "}";
            os.close();

            
            newton.Solve(pa,x,constraints);
            
            cvx.Compute(31);
            pos.Compute(31);
            
            const std::string filename = file+"_iter"+std::to_string(iRefine)+".txt";
            os.open(filename);
            os << "{"
            ExportVarArrow(newton)
            ExportVarArrow(pa)
            ExportVarArrow(cvx)
            ExportVarArrow(pos)
            ExportArrayArrow(NS::StdFromEigen(x))
            << "}";
            os.close();
            
            std::cout
            ExportVarArrow(newton)
            << "\n";
            
            // Checking error indicators, gradients
            
            
            
            std::cout << "Points before refinement : " ExportVarArrow(pa.Pts().size());
            const ScalarType refinementQuantile = 0.3;
            pa.Refine(refinementQuantile,x);
            pts = pa.Pts();
            std::cout ExportVarArrow(pts.size()) << "\n";
        }
    }
    
    struct {
        ScalarType operator()(const CGT::Point & p){
            const CGT::Vector v = p-CGT::Point(1,1);
            return std::max( 0., v.squared_length()-0.8);
        }
    } testFunc;
    
    void Ref0(int nRefine, int startSize=10){
        std::vector<CGT::Full_point> pts = MakeShape(startSize,ShapeType::Square);
        PrincipalAgent pa;
        pa.rt = CGT::RT(pts.begin(),pts.end());
        pa.rt.infinite_vertex()->info() = {InfiniteIndex, -1};
        std::ofstream os;
        os.open("Ref0.txt"); os << "{";
        
        for(int iRefine=0; iRefine<nRefine; ++iRefine){
            NS::VectorType x(pts.size());
            for(int i=0; i<pts.size(); ++i)
                x[i] = testFunc(pts[i].first.point());
            pa.Refine(0.3,x);
            pts = pa.Pts();
            
            // Also check Lipschitz regularity
            LipschitzConstraint lip(pts);
            
            os << '"' << "iter" << iRefine << '"' << "-> {"
            ExportVarArrow(pa)
            ExportVarArrow(lip)
            << "},";
        }
        os << "}";
        
    }
    
    void PALinear(std::vector<CGT::Full_point> pts, int nRefine, std::string file = "PALinear"){
        PrincipalAgent pa;
        pa.rt = CGT::RT(pts.begin(),pts.end());
        pa.rt.infinite_vertex()->info() = {InfiniteIndex, -1};
        std::ofstream os;
//        pa.hasQuadraticCost=false;
        
        for(int iRefine = 0; iRefine < nRefine; ++iRefine){
            
            ConvexityConstraint cvx(pts);
            PositivityConstraint pos;
            LipschitzConstraint lip(pts);
            lip.lastOfLineOnly = false;
            std::vector<NS::Functionnal*> constraints = {&cvx, &pos, &lip}; //, &lip
            
            NS::NewtonConstrained newton; //newton.verbose=true;
            newton.maxIter=50;
//            newton.opt.solverStrategy = NS::NewtonUnconstrained::SolverType::ConjugateGradient;
            newton.opt.solverStrategy = NS::NewtonUnconstrained::SolverType::Other;
            
            newton.multiplier = 1./pts.size();
            newton.multiplierBound *= newton.multiplier;
            newton.multiplierDamping = 0.7;
            
            NS::VectorType x(pts.size());
            for(int i=0; i<x.size(); ++i) x[i]=CGT::Parabola(pts[i].first)/10.;
            
            cvx.SetValues(NS::StdFromEigen(x));
            lip.SetValues(NS::StdFromEigen(x));
            cvx.Compute(31);
            lip.Compute(31);
            os.open("cvxTest.txt");
            os << "{" ExportVarArrow(cvx) ExportVarArrow(lip) << "}";
            os.close();
            
            
            newton.Solve(pa,x,constraints);
            
            cvx.Compute(31);
            pos.Compute(31);
            lip.Compute(31);
            
            const std::string filename = file+"_iter"+std::to_string(iRefine)+".txt";
            os.open(filename);
            os << "{"
            ExportVarArrow(newton)
            ExportVarArrow(pa)
            ExportVarArrow(cvx)
            ExportVarArrow(lip)
            ExportVarArrow(pos)
            ExportArrayArrow(NS::StdFromEigen(x))
            << "}";
            os.close();
            
            std::cout
            ExportVarArrow(newton)
            << "\n";
            
            // Checking error indicators, gradients
            
            
            
            std::cout << "Points before refinement : " ExportVarArrow(pa.Pts().size());
            const ScalarType refinementQuantile = 0.3;
            pa.Refine(refinementQuantile,x);
            pts = pa.Pts();
            std::cout ExportVarArrow(pts.size()) << "\n";
        }
    }

}

#endif
