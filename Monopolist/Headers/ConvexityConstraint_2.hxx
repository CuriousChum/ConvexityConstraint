#pragma once
//
//  ConvexityConstraint_2.hxx
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 28/01/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

// %%%%%%%%%%%%%%%%%%%%%%%%%% Convexity %%%%%%%%%%%%%%%%%%%%%%%%%%%

void ConvexityConstraint::SetValues(const std::vector<ScalarType> & heights)
{
    assert(heights.size()==Size());
    for(int i=0; i<Size(); ++i){
        CGT::Weighted_point & p = pts[i].first;
        p = CGT::Weighted_point(p.point(), CGT::Parabola(p) - heights[i]);
    }
    
    rt = CGT::RT(pts.begin(),pts.end());
    rt.infinite_vertex()->info() = {InfiniteIndex, -1};
    error = (int)rt.number_of_hidden_vertices();
}

ScalarType ConvexityConstraint::Height(IndexType index) const {
    assert(0<=index && index<Size());
    const CGT::Weighted_point & p = pts[index].first;
    return CGT::Parabola(p) - p.weight();
}

void ConvexityConstraint::Compute(FlagType r){
    Clean(r);
    using namespace CGAL_Traits;
    
    /*
    std::cout << "\n\n";
    print_range(std::cout, heights.begin(), heights.end());
    std::cout << "\n\n";*/

    std::vector<Vector> edges; // outward edges from current point
    std::vector<IndexType> indices; // neighbor indices
    std::vector<ScalarType> diffs; // height differences
    
    ScalarType tValue;
    std::vector<VecCoef> tJacobian;
    std::vector<MatCoef> tHessian;

    
    if(r & RValues) values.resize(Size(),ScalarType(0));
    if(r & RLogGrad) logGrad.resize(Size(),ScalarType(0));
    
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
            edges.push_back(vC->point() - vIt->point());
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
            
            if(r & (RLogGrad | RLogHessian | RJacobian) ) {
                tJacobian.push_back({i, cij*dj});
                tJacobian.push_back({j, cij*di-cj*dj});
                
                if (r & RLogHessian) {
                    tHessian.push_back({i,j, cij});
                    tHessian.push_back({j,i, cij}); // true symmetric hessian
                    tHessian.push_back({j,j, -cj});
                }
            }
        }
        
        if(tValue<0){
            std::ostream & os = std::cout
            ExportArrayArrow(edges)
            ExportVarArrow(vIt->point())
            ExportVarArrow(baseIndex)
            ExportArrayArrow(diffs)
            ExportArrayArrow(pts)
            << "\n";
            assert(false);
        }
        
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
        
        if(r & (RLogSum | RLogGrad | RLogHessian)){
            NLog(tValue,tJacobian,tHessian);
            Simplify(tHessian);
        }
        
        ++numberOfConstraints;
        if(r & RLogSum) logSum+=tValue;
        
        if(r & RLogGrad) {
            ScalarType sum=0;
            for(const VecCoef c : tJacobian){
                const IndexType i=indices[c.first];
                logGrad[i]+=c.second;
                sum+=c.second;
            }
            logGrad[baseIndex]-=sum;
        }
        
        if(r & RLogHessian){
            for(const MatCoef c : tHessian){
                const IndexType i = indices[c.i];
                const IndexType j = indices[c.j];
                const ScalarType val = c.v;
                
                logHessian.push_back(MatCoef{i,j,val});
                logHessian.push_back(MatCoef{baseIndex,baseIndex,val});
                logHessian.push_back(MatCoef{i,baseIndex,-val});
                logHessian.push_back(MatCoef{baseIndex,j,-val}); // Fine since tHessian is symmetric
            }
        }
    } // for finite vertices
    
    Simplify(jacobian);
    Simplify(logHessian);
}

void ConvexityConstraint::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(rt)
    ExportVarArrow(Size())
    ExportArrayArrow(pts)
    << '"' << "superclass" << '"' << "->";
    ConstraintType::PrintSelf(os);
    os << "}";
}


// %%%%%%%%%%%%%%%%%% Convexity constraints associated with the boundary %%%%%%%%%%%%%%%%%

std::vector<std::pair<int, std::unique_ptr<ConstraintType> > >
BoundaryConvexityConstraints(const std::vector<CGAL_Traits::Full_point> & pts){
	namespace CT = CGAL_Traits;
	std::vector<std::pair<int, std::unique_ptr<ConstraintType> > > result;
	
	// Find the possible boundary tags of the domain. (Each corresponds to an edge.)
	IndexType boundaries = 0;
	for(const CT::Full_point & p : pts) boundaries |= p.second.boundary;
	
	for(int iBoundary=0; iBoundary<8*sizeof(IndexType); ++iBoundary){
		if( !(boundaries & 1<<iBoundary) ) continue;
		
		// Select the points associated to this boundary tag.
		std::vector<CT::Full_point> bpts;
		for(const CT::Full_point & p : pts){
			if(p.second.boundary & 1<<iBoundary) bpts.push_back(p);}
		
		// By assumption, these points are aligned. Take two to find the direction.
		if(bpts.size()<=1) continue;
		const CT::Point p0 = bpts[0].first.point();
		const CT::Vector _v = bpts[1].first.point() - p0;
		const CT::Vector v = _v/sqrt(_v.squared_length());
		
		// Get the abcissa and indices associated to the boundary points
		std::vector<ScalarType> abcissa;
		std::vector<IndexType> indices;
		for(const CT::Full_point & p : bpts){
			abcissa.push_back( v*(p.first.point()-p0) );
			indices.push_back( p.second.index );
		}
		
		std::cout << "Boundary convexity constraints "
		ExportArrayArrow(abcissa)
		ExportVarArrow(v)
		ExportVarArrow(p0)
		<< std::endl;
		
		auto cvx1 = std::make_unique<Geometry_1::ConvexityConstraint>(abcissa);
		auto resampled = std::make_unique<ResampledConstraint>(std::move(cvx1),indices);
		result.push_back({iBoundary,std::move(resampled)});
	}
	return result;
}
