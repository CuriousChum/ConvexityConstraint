//
//  LipschitzConstraint.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 09/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_LipschitzConstraint_h
#define CGalTest_LipschitzConstraint_h

// Implement the constraint D_x u \leq 1, D_y u \leq 1, as required by the monopolist problem.

#include "Geometry_2.h"
namespace Geometry_2 {
    
    struct LipschitzConstraint : ConstraintType {
        LipschitzConstraint(const std::vector<CGT::Full_point> &);
        
        virtual void SetValues(const std::vector<ScalarType> & x){input=x;error=0;}
        virtual void ComputeValJacHess(FlagType);
        
        virtual void PrintSelf(std::ostream &) const;
        virtual std::string Name() const {return "LipschitzConstraint";}
        bool lastOfLineOnly = true;
    protected:
        std::vector<MatCoef> constraints;
        std::vector<ScalarType> input;
        virtual size_t InputSize() const {return input.size();}
    };

    LipschitzConstraint::LipschitzConstraint(const std::vector<CGT::Full_point> & pts){
        if(pts.empty()) return;
        struct OrderType {
            int i;         typedef std::tuple<ScalarType, ScalarType> TupleType;
            bool operator()(const CGT::Point & p, const CGT::Point & q) const {
                return TupleType{p[i],p[1-i]} < TupleType{q[i],q[1-i]};
            }
        } order;
        
        for(int i=0; i<2; ++i){
            order.i=i;
            std::map<CGT::Point,IndexType,OrderType > sorter(order);
            for(const auto & p : pts)
                sorter.emplace(p.first.point(), p.second.index);
            typedef std::pair<CGT::Point,IndexType> PI;
            
            std::array<PI,2> lastP; lastP[0] = *sorter.begin(); // donc...
            for(const PI & pi : sorter){
                if((!lastOfLineOnly || pi.first[i] != lastP[0].first[i]) && lastP[0].first[i]==lastP[1].first[i])
                    constraints.push_back({lastP[0].second,lastP[1].second,lastP[0].first[1-i]-lastP[1].first[1-i]});
                lastP[1]=lastP[0];
                lastP[0]=pi;
            }
            if(lastP[0].first[i]==lastP[1].first[i]) // Last line
                constraints.push_back({lastP[0].second,lastP[1].second,lastP[0].first[1-i]-lastP[1].first[1-i]});
        }
        numberOfConstraints = (int)constraints.size();
    }
    
    void LipschitzConstraint::ComputeValJacHess(FlagType r){
        
        /*
         numberOfConstraints=(int)input.size();
         if(r & RValues) values.reserve(numberOfConstraints);
         if(r & RJacobian) jacobian.reserve(numberOfConstraints);
         for(int i=0 ; i<numberOfConstraints; ++i) {
         if(r & RValues) values.push_back(input[i]);
         if(r & RJacobian) jacobian.push_back({i,i,1});
         }
         */

         
        /*
        values.push_back(input[0]);
        jacobian.push_back({0,0,1});
        
        error=0;
        */
        
        int counter=0;
        for(auto c : constraints){
            
            if(r & RValues) values.push_back(c.v-input[c.i]+input[c.j]);
            if(r & RJacobian) jacobian.push_back({counter,c.i,-1.});
            if(r & RJacobian) jacobian.push_back({counter,c.j, 1.});
            
            /* // Just bound constraints for testing
            if(r & RValues) values.push_back(1.-input[c.i]);
            if(r & RJacobian) jacobian.push_back({counter,c.i,-1.});
            */
            
            
/*            std::cout
            ExportVarArrow(c)
            ExportVarArrow(counter)
            ExportVarArrow(input[c.i])
            << "\n";*/
/*
            if(r & RValues) values.push_back(c.v - (input[c.i]-input[c.j]) );
            if(r & RJacobian){
                jacobian.push_back({counter, c.i, -1});
                jacobian.push_back({counter, c.j, 1});
            }
 */
            ++counter;
        }
//        std::cout << "\n";
    }
    
    void LipschitzConstraint::PrintSelf(std::ostream & os) const {
        os << "{"
        ExportArrayArrow(constraints)
        << '"' << "superclass" << '"' << "->";
        ConstraintType::PrintSelf(os);
        os << "}";
    }
}




/*void PositivityConstraint::SetValues(const std::vector<ScalarType> & val_){
 val = val_;
 error = std::accumulate(val.begin(), val.end(), 0, [](int k, ScalarType x) {return k+int(x<=0);});
 };
 
 
 void PositivityConstraint::ComputeValJacHess(FlagType r){
 numberOfConstraints=(int)val.size();
 if(r & RValues) values.reserve(numberOfConstraints);
 if(r & RJacobian) jacobian.reserve(numberOfConstraints);
 for(int i=0 ; i<numberOfConstraints; ++i) {
 if(r & RValues) values.push_back(val[i]);
 if(r & RJacobian) jacobian.push_back({i,i,1});
 }
 }*/


#endif
