#pragma once

//
//  DispatchAndRun.h
//  Monopolist
//
//  Created by Jean-Marie Mirebeau on 09/06/2022.
//

#include "JMM_CPPLibs/Output/FileIO.h"
typedef IO_<FileIO> IO;
typedef typename IO::Msg Msg;
typedef typename IO::WarnMsg WarnMsg;


//#include "ConvexityConstraint_1.h"
#include "ConvexityConstraint_2.h"
#include "ConvexityConstraint_3.h"

namespace {
using Constraint::ScalarType;
using Constraint::FlagType;
using Constraint::ConstraintType;

void Export(IO & io, FlagType flag, ConstraintType & c){
	io.currentSetter = TraitsIO::SetterTag::Compute;
	using CT = ConstraintType;
	std::vector<ScalarType> v,i,j,k; // IO only accepts scalars
	if(flag & CT::RLogSum){io.Set<ScalarType>("logSum",c.logSum);}
	if(flag & CT::RLogGrad){io.SetVector<ScalarType>("logGrad",c.logGrad);}
	if(flag & CT::RLogHessian){
		const size_t size = c.logHessian.size();
		v.resize(size); i.resize(size); j.resize(size);
		for(size_t s=0; s<size; ++s) {
			const Constraint::MatCoef & a =  c.logHessian[s];
			v[s] = a.v; i[s] = a.i; j[s] = a.j;
		}
		io.SetVector<ScalarType>("logHessian_v",v);
		io.SetVector<ScalarType>("logHessian_i",i);
		io.SetVector<ScalarType>("logHessian_j",j);
	}
	if(flag & CT::RValues){io.SetVector<ScalarType>("values",c.values);}
	if(flag & CT::RJacobian){
		const size_t size = c.jacobian.size();
		v.resize(size); i.resize(size); j.resize(size);
		for(size_t s=0; s<size; ++s) {
			const Constraint::MatCoef & a =  c.jacobian[s];
			v[s] = a.v; i[s] = a.i; j[s] = a.j;
		}
		io.SetVector<ScalarType>("jacobian_v",v);
		io.SetVector<ScalarType>("jacobian_i",i);
		io.SetVector<ScalarType>("jacobian_j",j);
	}
	if(flag & CT::RHessian){
		const size_t size = c.hessian.size();
		v.resize(size); i.resize(size); j.resize(size); k.resize(size);
		for(size_t s=0; s<size; ++s) {
			const Constraint::TensorCoef & a =  c.hessian[s];
			v[s] = a.v; i[s] = a.i; j[s] = a.j; k[s] = a.k;
		}
		io.SetVector<ScalarType>("hessian_v",v);
		io.SetVector<ScalarType>("hessian_i",i);
		io.SetVector<ScalarType>("hessian_j",j);
		io.SetVector<ScalarType>("hessian_k",k);
	}
}


void Run1(IO & io,FlagType flag,std::vector<ScalarType> heights,std::vector<bool> exclude){}

void Run2(IO & io,FlagType flag,std::vector<ScalarType> heights,std::vector<bool> exclude){
	using namespace Geometry_2;
	// Stupid conversion
	const auto points_ = io.GetVector<std::array<ScalarType,Dimension> >("points");
	std::vector<CGT::Point> points; points.reserve(points_.size());
	for(const auto & p : points_) points.push_back({p[0],p[1]});
	
	// Compute and export
	ConvexityConstraint cvx(points,heights,exclude,flag);
	Export(io,flag,cvx);
}

void Run3(IO & io,FlagType flag,std::vector<ScalarType> heights,std::vector<bool> exclude){
	using namespace Geometry_3;
	// Stupid conversion
	const auto points_ = io.GetVector<std::array<ScalarType,Dimension> >("points");
	std::vector<CGT::Point> points; points.reserve(points_.size());
	for(const auto & p : points_) points.push_back({p[0],p[1],p[2]});
	
	// Compute and export
	ConvexityConstraint cvx(points,heights,exclude,flag);
	Export(io,flag,cvx);
}


void Run(IO & io){
	// Input : points, values, and a flag about what to compute
	// Output : the required subgradient areas, their jacobians, etc
	io.arrayOrdering = TraitsIO::ArrayOrdering::RowMajor;
	typedef IO::ScalarType ScalarType;
	const auto dims = io.GetDimensions<ScalarType>("points");
	if(dims.size()!=2) {ExceptionMacro("Error: graph must be a depth two array");}
	const auto flag = Constraint::FlagType(io.Get<ScalarType>("flag"));
	
	const auto heights = io.GetVector<ScalarType>("heights");

	// Excluded points table ... TODO : tranfer more efficiently.
	const auto exclude_ = io.GetVector<ScalarType>("exclude");
	std::vector<bool> exclude; exclude.reserve(exclude_.size());
	for(const ScalarType & b : exclude_) exclude.push_back(bool(b));

/*	using CT = Constraint::ConstraintType;
	if(io.Get<ScalarType>("LogSum",0))     flag |= CT::RLogSum;
	if(io.Get<ScalarType>("LogGrad",0))    flag |= CT::RLogGrad;
	if(io.Get<ScalarType>("LogHessian",0)) flag |= CT::RLogHessian;
	if(io.Get<ScalarType>("Values",0))     flag |= CT::RValues;
	if(io.Get<ScalarType>("Jacobian",0))   flag |= CT::RJacobian;
	if(io.Get<ScalarType>("Hessian",0))    flag |= CT::RHessian;*/
			
	switch (dims.back()) {
		case 1: return Run1(io,flag,heights,exclude);
		case 2: return Run2(io,flag,heights,exclude);
		case 3: return Run3(io,flag,heights,exclude);
	}
}

}


