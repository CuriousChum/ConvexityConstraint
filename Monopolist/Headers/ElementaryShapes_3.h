#pragma once
//
//  ElementaryShapes_3.h
//  Monopolist
//
//  Created by Jean-Marie Mirebeau on 21/03/2022.
//

#include "Geometry_3.h"
#include "JMM_CPPLibs/Output/EnumToString.h"

namespace Geometry_3 {
enum class ShapeType {Cube};

/**
 Generate a list of points, which regularly sample an elementary shape.
 These point sets usually serve as input to compute a Delaunay triangulation.
 
 The shapes are :
  - The cube [1,2]^3
 
 The points come with
  - a null weight
  - an index
  - a flag, if they belong to a flat side of the shape.
 
 Input :
  - number of points on each side (approx)
  - shape : class of shape
  - theta (optional) : rotation of the shape (around its barycenter}, given as Euler angles
  - bary (optional) : barycenter of the shape
 */
#pragma message("Feature incomplete function (add more shapes)")
std::vector<CGT::Full_point>
MakeShape(int n, ShapeType shape){
	//,PointType theta={0.,0.,0.}, PointType bary={Infinity,Infinity,Infinity}
	const ScalarType h=1./(n-1);
	std::vector<CGT::Full_point> pts;
	int counter=0;
	
	typedef CGT::Weighted_point WP;
	typedef CGT::InfoType IT;
	ScalarType bary_;
	
	switch (shape) {
		case ShapeType::Cube: {

			// Create cube
			const ScalarType h = 1./n;
			const CGT::Vector ei(h,0,0), ej(0,h,0), ek(0,0,h);
			const CGT::Point origin(1.+0.5*h, 1.+0.5*h, 1.+0.5*h);
			
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
			break;}
	
	} // switch

	return pts;
} // MakeShape
} // namespace Geometry_3

template<> char const* enumStrings<Geometry_3::ShapeType>::data[] = {"Cube"};
