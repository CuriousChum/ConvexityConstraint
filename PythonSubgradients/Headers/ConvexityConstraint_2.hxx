#pragma once
//
//  ConvexityConstraint_2.hxx
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 28/01/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

// %%%%%%%%%%%%%%%%%%%%%%%%%% Convexity %%%%%%%%%%%%%%%%%%%%%%%%%%%

ConvexityConstraint::
ConvexityConstraint(const std::vector<Point> & pts,
					const std::vector<ScalarType> & heights_,
					const std::vector<bool> & exclude, FlagType flag):
heights(heights_){
	// Check sizes
	const size_t size = heights.size();
	if(pts.size()!=size || exclude.size()!=size)
		ExceptionMacro("Inconsistent data sizes");
	numberOfUnknowns=size; numberOfConstraints=size;
	
	// Generate set of full points
	using Full_point = CGT::Full_point;
	std::vector<Full_point> fpts; fpts.reserve(size);
	for(IndexType i=0; i<size; ++i)
		fpts.push_back({CGT::Weighted_point(pts[i],CGT::Parabola(pts[i])-heights[i]),
			InfoType{i,exclude[i]}});

	rt.insert(fpts.begin(),fpts.end());
	rt.infinite_vertex()->info() = {InfiniteIndex,true};
	if(rt.number_of_hidden_vertices() > 0)
		ExceptionMacro("ConvexityConstraint 2D : non-convex data (hidden vertices)");
	if(rt.number_of_vertices() != size)
		ExceptionMacro("Error : insertion not successful (likely duplicate point)")
	Compute(flag);
}

void ConvexityConstraint::ComputeValJacHess(FlagType r){
	using namespace CGAL_Traits;
    
    std::vector<Vector> edges; // outward edges from current point
    std::vector<IndexType> indices; // neighbor indices
    std::vector<ScalarType> diffs; // height differences
    
    ScalarType tValue;
    std::vector<VecCoef> tJacobian;
    std::vector<MatCoef> tHessian;
	
    if(r & RValues) values.resize(heights.size(),ScalarType(0));
	
    // Iterate over vertices.
    for(Finite_vertices_iterator vIt = rt.finite_vertices_begin(), vEnd = rt.finite_vertices_end();
        vIt!=vEnd; ++vIt){
        
        edges.clear();
        indices.clear();
        diffs.clear();
        tValue=0;
        tJacobian.clear();
        tHessian.clear();
        
        const IndexType baseIndex = vIt->info().index;
        const ScalarType baseHeight = Height(baseIndex);
		
        // Iterate over neighbor vertices, in clockwise order
        Vertex_circulator vC = rt.incident_vertices(vIt), vDone(vC);
        bool onBoundary=false;
        do {
            const IndexType index = vC->info().index;
            if(index==InfiniteIndex){onBoundary=true; break;}
            edges.push_back(vC->point().point() - vIt->point().point());
            indices.push_back(index);
            diffs.push_back(Height(index) - baseHeight);
        } while(++vC != vDone);
        
        // Abort if on boundary
        if(onBoundary || vIt->info().OnBoundary()){
            if(r & RValues) values[baseIndex]=Infinity;
            continue;
        }
        
        // Construct value, Jacobian, Hessian
        const IndexType n= (IndexType)edges.size();
        
        for(IndexType j=0; j<n; ++j){
            const IndexType i=(j+n-1)%n, k=(j+1)%n;
            // evaluate determinants
            
            const ScalarType
            Dij = CGAL::determinant(edges[i],edges[j]),
            Dik = CGAL::determinant(edges[i],edges[k]),
            Djk = CGAL::determinant(edges[j],edges[k]);
            
            const ScalarType
            cij = 1./Dij,
            cj  = Dik/(Dij*Djk);
            
            const ScalarType di=diffs[i], dj=diffs[j];
            
            tValue+=cij*di*dj-0.5*cj*dj*dj;
            
            if(r & RJacobian) {
                tJacobian.push_back({i, cij*dj});
                tJacobian.push_back({j, cij*di-cj*dj});
            }
			if(r & RHessian) {
				tHessian.push_back({i,j, cij});
				tHessian.push_back({j,i, cij}); // true symmetric hessian
				tHessian.push_back({j,j, -cj});
			}
        }
        
		assert(tValue>=0);
        
        Simplify(tJacobian);
        if(r & RValues) values[baseIndex]=tValue;
        
        if(r & RJacobian){
            ScalarType sum=0;
            for(const VecCoef c : tJacobian){
                const IndexType i=indices[c.first];
                jacobian.push_back(MatCoef{baseIndex, i, c.second});
                sum+=c.second;
            }
            jacobian.push_back({baseIndex, baseIndex,-sum});
        }
        
		if(r & RHessian){
			for(const MatCoef c : tHessian){
				const IndexType i = indices[c.i];
				const IndexType j = indices[c.j];
				const ScalarType val = c.v;
				
				hessian.push_back(TensorCoef{baseIndex,i,j,val});
				hessian.push_back(TensorCoef{baseIndex,baseIndex,baseIndex,val});
				hessian.push_back(TensorCoef{baseIndex,i,baseIndex,-val});
				hessian.push_back(TensorCoef{baseIndex,baseIndex,j,-val});
				// Fine since tHessian is symmetric w.r.t. i,j
			}
		}
    } // for finite vertices
}

void ConvexityConstraint::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(rt)
    << '"' << "superclass" << '"' << "->";
    ConstraintType::PrintSelf(os);
    os << "}";
}
