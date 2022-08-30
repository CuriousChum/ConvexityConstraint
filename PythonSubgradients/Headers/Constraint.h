#pragma once
//
//  Constraint.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 29/01/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#include <vector>
#include <limits>
#include <iostream>

#include <numeric>
//#include "NewtonSolvers.h"

namespace Constraint {

typedef double ScalarType;
typedef int IndexType;
typedef unsigned int FlagType; // Conjunction of requests


const ScalarType Infinity = std::numeric_limits<ScalarType>::infinity();
const IndexType BadIndex = std::numeric_limits<IndexType>::max();
const IndexType InfiniteIndex = BadIndex-1;
//const FlagType RoundBoundary = 1<<31, ExcludedBoundary = 1<<30;


struct InfoType { // Additional info for each vertex
	IndexType index=BadIndex;
	bool boundary=false;
	// FlagType boundary=0;
	// InfoType(IndexType index_, IndexType boundary_):index(index_), boundary(boundary_){};
	// InfoType(){};
	bool OnBoundary() const {return boundary;}
	friend std::ostream & operator << (std::ostream & os, const InfoType & p){
		return os << "{" << p.index << "," << p.OnBoundary() << "}";}
 };

typedef std::pair<IndexType, ScalarType> VecCoef;
struct MatCoef;
struct TensorCoef; // Trivariant tensor

static void Sqrt(ScalarType &, std::vector<VecCoef> &, std::vector<MatCoef> &); // true hessian

template<typename Coef> struct SameIndices;

template<typename Coef, typename Traits = SameIndices<Coef> >
void Simplify(std::vector<Coef> &,ScalarType tol=0);
void Sqrt(ScalarType & val, std::vector<VecCoef> & g, std::vector<MatCoef> & h);
void NLog(ScalarType & val, std::vector<VecCoef> & g, std::vector<MatCoef> & h); //-log


/**
 Barrier function - Sum_i log(c_i(x)), for the constraints c_i(x)>0 for all i.
 
 Only SetValues and ComputeValJacHess need to be specialized in the subclass.
 */

struct ConstraintType {
	enum Request {
		RLogSum = 1<<0, RLogGrad = 1<<1, RLogHessian = 1<<2,
		RValues = 1<<3, RJacobian = 1<<4, RHessian = 1<<5
	};
	
	/// Set the values where the constraint is to be evaluated
	virtual void SetValues(const std::vector<ScalarType> &)
	{throw "ConstraintType::Compute error : must be specialized";};
	virtual void Compute(FlagType);
	
	// Read only
	IndexType numberOfConstraints=BadIndex, numberOfUnknowns=BadIndex;
	ScalarType logSum = Infinity;
	std::vector<ScalarType> logGrad;
	std::vector<MatCoef> logHessian; // Upper triangular part only
	
	std::vector<ScalarType> values; // The values of the constraints
	std::vector<MatCoef> jacobian;
	std::vector<TensorCoef> hessian;
	
	virtual std::string Name() const {return "Unspecified constraint name";}
	virtual void PrintSelf(std::ostream & os) const;
	friend std::ostream & operator << (std::ostream & os, const ConstraintType & c){
		c.PrintSelf(os); return os;}
	
	// nlopt interface
	static ScalarType GeometricMean(const std::vector<ScalarType> &,
									std::vector<ScalarType> &, void*);
protected:
	virtual void Clean(FlagType);
	virtual void Check(FlagType);
	
	/// Compute the value, jacobian, and hessian of the constraint
	virtual void ComputeValJacHess(FlagType)
	{throw "ConstraintType::ComputeValJacHess error : must be specialized";}
	
	/// Compute the value, jacobian and hessian of the barrier, defined as the sum of the logarithms of the constraint
	virtual void ComputeLogarithms(FlagType);	
};

#include "Constraint.hxx"
}
