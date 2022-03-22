#pragma once
//
//  ElementaryShapes_2.h
//  Monopolist
//
//  Created by Jean-Marie Mirebeau on 18/03/2022.
//

#include "JMM_CPPLibs/Output/EnumToString.h"
namespace Geometry_2 {

enum class ShapeType {Square, Triangle, Circle};
/**
 Generate a list of points, which regularly sample an elementary shape.
 These point sets usually serve as input to compute a Delaunay triangulation.
 
 The shapes are :
  - The square [1,2]^2
  - The disk centered at (1.5,1.5) with radius 0.5
  - An equilateral triangle.
 
 The points come with
  - a null weight
  - an index
  - a flag, if they belong to a flat side of the shape.
 
 Input :
  - number of points on each side (approx)
  - shape : class of shape
  - theta (optional) : rotation of the shape (around its barycenter}
  - bary (optional) : barycenter of the shape
 
 */
std::vector<CGT::Full_point>
MakeShape(int n, ShapeType shape,ScalarType theta=0., PointType bary={Infinity,Infinity}){
	const ScalarType h=1./(n-1);
	const ScalarType pi = 4.*atan(1.);
	std::vector<CGT::Full_point> pts;
	int counter=0;
	
	typedef CGT::Weighted_point WP;
	typedef CGT::InfoType IT;
	ScalarType bary_;
	
	switch (shape) {
		case ShapeType::Square:
			for(int i=0; i<n; ++i)
				for(int j=0; j<n; ++j){
					IndexType flag =
					(i==0 ? 1 : 0) | (j==0 ? 2 : 0) |
					(i==n-1 ? 4 : 0) | (j==n-1 ? 8 : 0);
					pts.push_back({WP({1.+i*h,1.+j*h},0),
						IT(counter++,flag)});
				}
			bary_=1.5;
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
			bary_=1.5;
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
			bary_=1.1053600322580646;
			break;}
	}
	
	// Translate and rotate the shape
	if(bary[0]!=Infinity || bary[1]!=Infinity || theta!=0){
		for(ScalarType & x : bary) {if(x==Infinity) x=bary_;}
		ScalarType c=cos(theta), s=sin(theta);
		for(CGT::Full_point & p : pts){
			auto & q = p.first.point();
			ScalarType x=q.x(), y=q.y();
			x-=bary_; y-=bary_;
			ScalarType xo=x; x = c*x-s*y; y=s*xo+c*y;
			x+=bary[0]; y+=bary[1];
			p.first = WP({x,y},0); // CGAL coordinates must remain constant, so a new point is made.
		}
	}
	return pts;
}

} // namespace Geometry_2

template<> char const* enumStrings<Geometry_2::ShapeType>::data[] = {"Square", "Triangle", "Circle"};
