//
//  PrincipalAgent_3_Test.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 13/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_PrincipalAgent_3_Test_h
#define CGalTest_PrincipalAgent_3_Test_h

#include "Geometry_3.h"
#include <fstream>

namespace PrincipalAgent_3_Test {
    using namespace Geometry_3;
    
    
    void PA0(int n=4, ScalarType eps=1){
        
        //Principal agent problem in 3D.
        std::ofstream os;
        
        // Create cube
        const ScalarType h = 1./n;
        const CGT::Vector ei(h,0,0), ej(0,h,0), ek(0,0,h);
        const CGT::Point origin(1.+0.5*h, 1.+0.5*h, 1.+0.5*h);
        
        std::vector<CGT::Full_point> pts;
        int counter=0;
        for(int i=0; i<n; ++i)
            for(int j=0; j<n; ++j)
                for(int k=0; k<n; ++k)
                    pts.push_back({
                        CGT::Weighted_point(origin + i*ei+j*ej+k*ek, 0.),
                        CGT::InfoType(counter++,
                                      (i==0 ?  1 : i==n-1 ?  2 : 0) |
                                      (j==0 ?  4 : j==n-1 ?  8 : 0) |
                                      (k==0 ? 16 : k==n-1 ? 32 : 0) )
                    });
        
        NS::VectorType xx(pts.size());
        for(int i=0; i<pts.size(); ++i)
            xx[i]=CGT::Parabola( pts[i].first.point() );
        
        // Setup constraints and objective
        
        PositivityConstraint pos;
        ConvexityConstraint cvx(pts);
        std::vector<NS::Functionnal*> constraints = {&cvx,&pos};
        
        ScalarType val = cvx.Value(xx);
        PrincipalAgent pa(cvx.rt);

        /*{
            std::ostream & os = std::cout
            ExportVarArrow(val)
            ExportArrayArrow(cvx.values)
            ExportVarArrow(pa)
//            ExportVarArrow(cvx)
            << "\n";
        }*/
        
//        return;
        

        // Setup Newton
        
        NS::NewtonConstrained newton;
        newton.maxIter=50;
        
        newton.multiplier = 1./pts.size();
        //newton.multiplierBound = eps*newton.multiplier;
        newton.multiplierBound = 4e-7; //1.5e-05;
        
        
        newton.Solve(pa,xx,constraints);
        
        cvx.Compute(31);
        
        const std::string filename = "PAO_3D_Cube.txt";
        os.open(filename);
        os << "{"
        ExportVarArrow(pa)
        ExportVarArrow(newton)
        ExportVarArrow(cvx)
        << "}";
        
    }
    
    
    void Tetra0(){
        
        using CGT::Weighted_point;
        using CGT::Point;
        using CGT::InfoType;
        // Test basis function gradients on a simplex.
        int counter=0;
        std::vector<CGT::Full_point> pts =
        {
            {Weighted_point(Point(1,0,0),0),InfoType(counter++,1)},
            {Weighted_point(Point(0,1,0),0),InfoType(counter++,1)},
            {Weighted_point(Point(0,0,1),0),InfoType(counter++,1)},
            {Weighted_point(Point(0,0,0),0),InfoType(counter++,1)}
        };
        
        CGT::RT rt(pts.begin(), pts.end());
        std::array<CGT::Vector,4 > grads;
        std::array<int,4> indices;
        ScalarType volume;
        
        CGT::Cell_handle fh = rt.finite_cells_begin();
        for(int i=0; i<4; ++i) indices[i]=fh->vertex(i)->info().index;
        
        const ScalarType degeneracy =
        CGT::Degeneracy(fh, volume, grads);
        
        std::ostream & os = std::cout
        ExportArrayArrow(grads)
        ExportVarArrow(volume)
        ExportVarArrow(degeneracy)
        ExportArrayArrow(indices)
        << "\n";
    }
    
    void Octa0(){
        using CGT::Weighted_point;
        using CGT::Point;
        using CGT::InfoType;
        // Test basis function gradients on a simplex.
        int counter=0;
        std::vector<CGT::Full_point> pts =
        {
            {Weighted_point(Point(0,0,0),0),InfoType(counter++,0)},
            {Weighted_point(Point(1,0,0),0),InfoType(counter++,1)},
            {Weighted_point(Point(0,1,0),0),InfoType(counter++,1)},
            {Weighted_point(Point(0,0,1),0),InfoType(counter++,1)},
            {Weighted_point(Point(-1,0,0),0),InfoType(counter++,1)},
            {Weighted_point(Point(0,-1,0),0),InfoType(counter++,1)},
            {Weighted_point(Point(0,0,-1),0),InfoType(counter++,1)}
        };
        
        CGT::RT rt(pts.begin(), pts.end());

        
        NS::VectorType xx(pts.size());
        for(int i=0; i<pts.size(); ++i)
            xx[i]=CGT::Parabola( pts[i].first.point() );
        
        ConvexityConstraint cvx(pts);
        
        ScalarType val = cvx.Value(xx);
        cvx.Compute(-1);
        {
            std::ostream & os = std::cout
            ExportVarArrow(val)
            ExportArrayArrow(cvx.values)
            ExportVarArrow(cvx)
            ExportVarArrow(rt)
            << "\n";
        }
        
        return;
    }

}

#endif
