#pragma once
//
//  ConvexityConstraint_2.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 28/01/2015. Some from Q. Merigot.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#include <vector>

#include "Geometry_2.h"
#include "ConvexityConstraint_1.h"

namespace Geometry_2 {

/**
 This functional evaluates the area of the subgradient cells associated with all the interior discretization points.
 Taking their logarithms, we obtain a barrier function for convexity in the domain interior.
 
 The boundary points are ignored however.
 */
struct ConvexityConstraint : ConstraintType {
    typedef CGAL_Traits::Full_point Full_point;
    ConvexityConstraint(std::vector<Full_point> pts_):pts(pts_){};
    
    virtual void SetValues(const std::vector<ScalarType> &);
    virtual void Compute(FlagType);
    // To do : work with ComputeValJacHess as others ?

    CGAL_Traits::RT rt; // To do : read only
    virtual void PrintSelf(std::ostream & os) const;
    virtual std::string Name() const {return "ConvexityConstraint";}
    const std::vector<Full_point> & GetPts() const {return pts;}
    size_t Size() const {return pts.size();}    
protected:
    std::vector<Full_point> pts;
    ScalarType Height(IndexType index) const;
};
    
    
    
    
    /* // Later
struct LipschitzRegularityConstraint : ConstraintType {
    LipschitzRegularityConstraint(const CGAL_Traits::RT & rt);
    
    
protected:
    std::vector<MatCoef> edges;
};
*/
    
    /*
     std::vector<ScalarType> heights;
     const std::vector<PointType> & points;
     std::vector<InfoType> info;
     */
    
    /*
     typedef CGAL_Traits::InfoType InfoType;
     
     ConvexityConstraint(const std::vector<PointType> & points_, std::vector<InfoType> info_ = {})
     :points(points_), nPts(points_.size()), info(info_){
     if(info.empty()) for(int i=0; i<nPts; ++i) info.push_back(InfoType(i,0));};
     */


    
#include "ConvexityConstraint_2.hxx"

}
