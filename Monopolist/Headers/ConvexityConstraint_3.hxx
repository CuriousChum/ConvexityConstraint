#pragma once
//
//  ConvexityConstraint_3.hxx
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 13/02/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

ConvexityConstraint::ConvexityConstraint(const std::vector<CGT::Full_point> & pts_, ScalarType degen_)
:pts(pts_), degen(degen_) {
    // Check that indices are valid.
    NS::VectorType check(pts.size());
    check.fill(std::numeric_limits<ScalarType>::quiet_NaN());
    for(const CGT::Full_point & p : pts){
        const IndexType index = p.second.index;
        if(index<0 || index>=pts.size())
			throw NS::DataError("Convexity constraint 3D: Invalid indices");
        check[index]=0;
    }
    if(check.hasNaN())
		throw NS::DataError("Convexity constraint error: Invalid indices");
}

void ConvexityConstraint::SetValues(const std::vector<ScalarType> & x){
	numberOfUnknowns = (IndexType) x.size();
    using namespace CGT;
    if(x.size()!=pts.size())
		throw DataError("Convexity constraint error: Invalid vector size");
    for(Full_point & p : pts) {
        const ScalarType w = Parabola(p.first.point()) - x[p.second.index];
        p.first = Weighted_point(p.first.point(),w);
    }
    rt = RT(pts.begin(), pts.end());
	if(rt.number_of_vertices()!=pts.size())
		throw NS::DomainError("ConvexityConstraint 3D : non-convex data");
//    std::cout << "Computing0 " << "\n";

}

// %%%%%%%%%%%%%%%%%%%%% Geometry %%%%%%%%%%%%%%%%%%%%%%

void ConvexityConstraint::ComputeValJacHess(FlagType r){
    using namespace CGT;
    
    if(r & RValues) values.resize(pts.size(),ScalarType(0));
    
//    std::cout << "Computing " << r << "\n";
    std::ostream & os = std::cout;
    
    std::vector< VecCoef > dualFacetGradient;
    
    // Iterate over edges. Each edge defines an orthogonal facet of a subgradient cell.
    for(Finite_edges_iterator eIt = rt.finite_edges_begin(), eEnd = rt.finite_edges_end();
        eIt!=eEnd; ++eIt){
        
        // Get the endpoints
        const Cell_handle containingCell = eIt->first;
        Vertex_handle
        source = containingCell->vertex(eIt->second),
        target = containingCell->vertex(eIt->third);
        
        if(source->info().boundary && target->info().boundary) continue;
        const Vector eVecRaw = target->point() - source->point();
        const Vector eVec = eVecRaw / eVecRaw.squared_length();
        
        // This normalization of eVec yields dual facets area and gradient divided by |eVecRaw|,
        // which corresponds to the subgradient measure differentiation.

        /*
        os
        ExportVarArrow(eVec)
        ExportVarArrow(source->point())
        ExportVarArrow(target->point())
        << "\n";*/
        
        ScalarType dualFacetArea=0;
        dualFacetGradient.clear();
        
        Cell_circulator fIt = rt.incident_cells(*eIt), fEnd = fIt;
        do{
            std::array<Vector,Dimension+1> fGrad, gGrad; // gradients of basis functions
            ScalarType fVol, gVol;
            
            Cell_circulator gIt = fIt;
            ++gIt;
            
            if(Degeneracy(fIt, fVol, fGrad) <= degen || // does compute fVol, fGrad
               Degeneracy(gIt, gVol, gGrad) <= degen )
                throw "ConvexityConstraint error : Degenerate cell."; // Fix me !!
            
            Vector fG(0,0,0), gG(0,0,0); // true gradients
            std::array<ScalarType,Dimension+1> fIndices, gIndices;
            
            for(int i=0; i<=Dimension; ++i){
                const Vertex_handle fP = fIt->vertex(i), gP = gIt->vertex(i);
                
                fG = fG + fGrad[i] * Height(fP->point());
                gG = gG + gGrad[i] * Height(gP->point());
                
                fIndices[i] = fP->info().index;
                gIndices[i] = gP->info().index;
            }
            
            dualFacetArea += CGAL::determinant(eVec,fG,gG);
            
            const Vector fCross =   CGAL::cross_product(eVec, fG);
            const Vector gCross = - CGAL::cross_product(eVec, gG);
            
            for(int i=0; i<=Dimension; ++i){
                dualFacetGradient.push_back({fIndices[i], gCross*fGrad[i]});
                dualFacetGradient.push_back({gIndices[i], fCross*gGrad[i]});
            }
            
        } while (++fIt != fEnd);
        
        Simplify(dualFacetGradient);
        
//        os << "\n" ExportArrayArrow(dualFacetGradient) << "\n";
        
        // The (normalized) dual facet area is the differential of the voronoi cell measure,
        // w.r.t. heightDiff
        
        for(int orientation=0; orientation<2; ++orientation, std::swap(source,target)){
        
            // Upon reversing orientation, eVec changes sign, but circulation around the edge also changes. So dualFacetArea and dualFacetGradient are unchanged.
            
            const ScalarType heightDiff = Height(target->point()) - Height(source->point());
            
            /*
            os
            ExportVarArrow(dualFacetArea)
            ExportVarArrow(heightDiff)
            << "\n";*/

            const IndexType
            sourceIndex = source->info().index,
            targetIndex = target->info().index;
            
            if( source->info().boundary ==0 ){
                if(r & RValues) values[sourceIndex] += heightDiff * dualFacetArea / 3.;
                if(r & RJacobian) {
                    jacobian.push_back({sourceIndex,targetIndex,dualFacetArea});
                    jacobian.push_back({sourceIndex,sourceIndex,-dualFacetArea});
                }
                if(r & RHessian)
                    for(const VecCoef & c : dualFacetGradient){
                        hessian.push_back({sourceIndex,targetIndex,c.first,   c.second}); 
                        hessian.push_back({sourceIndex,sourceIndex,c.first, - c.second});
                    }
            } // if not on boundary
        } // for orientation
    }
    
    // By convention, put infinite values for boundary indices.
    for(Finite_vertices_iterator vIt = rt.finite_vertices_begin(), vEnd = rt.finite_vertices_end();
        vIt!=vEnd; ++vIt){
        const InfoType info = vIt->info();
        if(info.boundary)
            values[info.index] = Infinity;
    }
}

// %%%%%%%%%%%%%%%%%%%%%%%% Export %%%%%%%%%%%%%%%%%%%%%%%

void ConvexityConstraint::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(rt) // pts come with rt
    ExportVarArrow(degen)
	<< "\"superclass\"->";
    ConstraintType::PrintSelf(os);
    os << "}";
}

// %%%%%%%%%%%%%%%%% Boundary constraints %%%%%%%%%%%%%%%
std::map<FlagType, std::unique_ptr<ConstraintType> >
BoundaryConvexityConstraints(const std::vector<CGT::Full_point> & pts){
	typedef std::bitset<8*sizeof(FlagType)> FlagSet;
	
	std::map<FlagType, std::unique_ptr<ConstraintType> > result;
	std::set<FlagType> edges, faces;
	for(const CGT::Full_point & p : pts){
		const FlagType flag = p.second.boundary & ~RoundBoundary;
		switch(FlagSet(flag).count()){
			case 1:faces.insert(flag); break;
			case 2:edges.insert(flag); break;
		}
	}
	
	for(FlagType edge : edges){
		std::vector<CGT::Full_point> bpts;
		for(const CGT::Full_point & p : pts){
			if((p.second.boundary & edge) == edge) bpts.push_back(p);}
		
		// By assumption, these points are aligned. Take two to find the direction.
		if(bpts.size()<=1) continue;
		const CGT::Point p0 = bpts[0].first.point();
		const CGT::Vector _v = bpts[1].first.point() - p0;
		const CGT::Vector v = _v/sqrt(_v.squared_length());
		
		// Get the abcissa and indices associated to the boundary points
		std::vector<ScalarType> abcissa;
		std::vector<IndexType> indices;
		for(const CGT::Full_point & p : bpts){
			abcissa.push_back( v*(p.first.point()-p0) );
			indices.push_back( p.second.index );
		}
		
		assert(std::is_sorted(abcissa.begin(),abcissa.end())); // TODO : sort if needed
		
		auto cvx1 = std::make_unique<Geometry_1::ConvexityConstraint>(abcissa);
		auto resampled = std::make_unique<ResampledConstraint>(std::move(cvx1),indices);
		result.insert({edge,std::move(resampled)});
	}
	
	for(FlagType face : faces){
		std::vector<CGT::Full_point> bpts;
		for(const CGT::Full_point & p : pts){
			if(p.second.boundary & face) bpts.push_back(p);}

		// TODO : Construct an orthogonal basis of the plane spanned
		assert(false);
//		const CGT::Vector v0,v1;
		
//		std::vector<Geometry_2::Cgal_Traits::Point>
		
	}
	
}
