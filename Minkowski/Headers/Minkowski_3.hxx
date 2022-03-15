//
//  Minkowski_3.hxx
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 03/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef Minkowski_Minkowski_3_hxx
#define Minkowski_Minkowski_3_hxx

int Functional::Compute(N::VectorType & input, N::VectorType & output){
    ErrCode error;
    error=MakeHull(input);
	if(error!=ErrCode::NoError) return (int)error;
	
    Compute(input);
    output = facetAreas;
    Normalize(input);
    return (int)error;
}


CGT::Vector Functional::Point(int i) const {
    assert(0<=i && i<Size());
    G::Vector pt = pts.col(i);
    return CGT::Vector( pt[0], pt[1], pt[2]);
}

/// Computes the polytope K(a), see class description.
Functional::ErrCode Functional::MakeHull(const N::VectorType & input) {
    assert(input.size()==Size());
    
    std::vector<CGT::Point> scaledPts;
    scaledPts.reserve(Size());
    for(int i=0; i<Size(); ++i){
        const ScalarType val = input[i];
        const CGT::Vector p = Point(i);
        if(val<=0) return ErrCode::NonPositiveEntries;
        scaledPts.push_back( CGT::Point( p.x()/val, p.y()/val, p.z()/val ) );
    }
    
    CGAL::convex_hull_3(scaledPts.begin(), scaledPts.end(), poly);
    
	/* If some vertices are interior, then the dual polytope has empty facets,
	which is not acceptable for our optimization purposes. */
    if( poly.size_of_vertices() < Size() ) return ErrCode::InteriorPoints;
    
	// In the rest of this function, we assign indices to the vertices of the polyhedron.
    struct OrderingType {
        typedef std::array<ScalarType,3> ArrayType;
        bool operator()(CGT::Point a,CGT::Point b) const {
            return ArrayType{a.x(),a.y(),a.z()} < ArrayType{b.x(),b.y(),b.z()};}
    };
    
    std::map<CGT::Point,IndexType,OrderingType> sorter;
    for(int i=0; i<scaledPts.size(); ++i)
        sorter.emplace(scaledPts[i],i);
    
    for(auto it=poly.vertices_begin(); it!=poly.vertices_end(); ++it)
        it->index = sorter.find(it->point())->second;
    return ErrCode::NoError;
}

void Functional::Compute(const N::VectorType & input){
    jacobian.clear();
    facetAreas.resize(Size()); facetAreas.setZero();
    facetMoments.resize(Dimension,Size()); facetMoments.setZero();
    
    // The facet area derivatives are determined by the edges
    for(auto eit = poly.halfedges_begin(); eit != poly.halfedges_end(); ++eit){
        auto opp = eit->opposite();
        
        const IndexType
        si= opp->vertex()->index,
        ti= eit->vertex()->index,
        ui= opp->next()->vertex()->index,
        vi= eit->next()->vertex()->index;
        
        if(si>ti) continue; // -> Count only once each edge
        
        const CGT::Vector
        s=Point(si),
        t=Point(ti),
        u=Point(ui),
        v=Point(vi);
        
        const ScalarType
        sa=input[si],
        ta=input[ti],
        ua=input[ui],
        va=input[vi];
		
		/* Geometric computation here.
		 The triangles 
		 
		 */
        
        
        typedef CGT::K::Vector_2 Vector_2;
        typedef CGT::K::Aff_transformation_2 Aff_transformation_2;
        
        const auto Gram = Aff_transformation_2(s*s,s*t,s*t,t*t);
        const auto GramInv = Gram.inverse();
        
        const Vector_2
        ul = Vector_2(u*s,u*t).transform(GramInv),
        vl = Vector_2(v*s,v*t).transform(GramInv);
        
        const ScalarType
        dstv = CGAL::determinant(s,t,v),
        dtsu = CGAL::determinant(t,s,u); // both positive usually
        
        const ScalarType // Normalized edge length in dual polytope
        deltaV = (va-vl*Vector_2(sa,ta))/dstv,
        deltaU = (ua - ul*Vector_2(sa,ta))/dtsu;
        const ScalarType
        delta = deltaV + deltaU,
        st = s*t;
        
        // Contributions to facet area derivatives
        jacobian.push_back({si,si,-delta*st});
        jacobian.push_back({si,ti,delta});
        jacobian.push_back({ti,si,delta});
        jacobian.push_back({ti,ti,-delta*st});
        
        /*
        // Second derivatives of facet areas. // Not needed here
        if(r & RLogHessian){
            const std::array<IndexType,2> sti = {si,ti};
            for(IndexType ai : sti){
                for(IndexType bi : sti){
                    const ScalarType m = (ai==bi ? -st : 1);
                    
                    hessian.push_back({ai,bi,si, m*(-vl.x()/dstv -ul.x()/dtsu)});
                    hessian.push_back({ai,bi,ti, m*(-vl.y()/dstv -ul.y()/dtsu)});
                    hessian.push_back({ai,bi,ui, m/dtsu});
                    hessian.push_back({ai,bi,vi, m/dstv});
                }
            }
        }
        */
        
        
        // Now computing face moments
        // First get point p inside the edge, and edge tips pu, pv.
        const Vector_2 pCoef = Vector_2(sa,ta).transform(GramInv);
        
        const CGT::Vector
        p = pCoef.x()*s + pCoef.y()*t,
        crossST = CGAL::cross_product(s,t);
        
        const CGT::Vector // common edge tips
        pv = p+deltaV*crossST,
        pu = p-deltaU*crossST;
        
        const CGT::Vector // barycenters of facet sub-triangles
        bS = (s*sa+pu+pv)/3.,
        bT = (t*ta+pu+pv)/3.;

        facetMoments.col(si) += GFromCGT(bS*(CGAL::determinant(s,pu,pv)/2.));
        facetMoments.col(ti) += GFromCGT(bT*(CGAL::determinant(t,pv,pu)/2.));
        
        // TO DO : export pu,pv for visualization ?
    }
    
    
    //Simplify(jacobian);
    //Simplify(hessian);
    
    // Compute the facet areas, from their partial derivatives, using Euler Formula
    facetAreas.resize(Size()); facetAreas.setZero();
    for(auto c : jacobian)          facetAreas[c.row()]+=input[c.col()]*c.value();
    facetAreas/=2;
    
    // Compute the total volume, for the facet areas, using Euler formula.
    volume = facetAreas.dot(input)/3.;
    moment = (facetMoments*input)/4.;
}

void Functional::Normalize(N::VectorType & input) {
    input -= pts.transpose() * Barycenter();
}

int Functional::DescentDirection(N::VectorType & diff) {
    if(referenceIndices[0]==N::BadIndex) SetReferenceIndices();
    for(int i: referenceIndices)
        jacobian.push_back({i,i,1});
    N::SparseMatrixType m;
    m.resize(Size(), Size());
    m.setFromTriplets(jacobian.begin(), jacobian.end());
    
    Eigen::SimplicialLDLT<N::SparseMatrixType> solver(m);
    diff = solver.solve(diff);
    N::VectorType & dir = diff;
    dir=-dir;
    
    // Correction
    const G::Vector dirMoment = facetMoments*dir;
    const G::Vector correction = (facetMoments*pts.transpose()).inverse()*dirMoment;
    dir -= pts.transpose()*correction;
    return (int)ErrCode::NoError;
}

/*
bool Functional::CheckData(){
    // Linear independence of reference indices
    Eigen::Matrix<ScalarType,Dimension,Dimension> m;
    for(int i=0; i<Dimension; ++i)
        m.col(i) = pts.col(referenceIndices[i]);
    const ScalarType mu = fabs(m.determinant())/ std::pow(m.norm(),Dimension);
    std::cout ExportVarArrow(mu) << "\n";
    return mu>1e-4;
}
*/

void Functional::SetReferenceIndices(){
    const int maxIt = 100;
    int i;
    for(i=0; i<maxIt; ++i){
        for(int j=0; j<Dimension; ++j)
            referenceIndices[j] = std::rand()%Size();
        if(ReferenceIndicesConditionNumber()>referenceIndicesConditionNumberThreshold)
            break;
    }
    if(i==maxIt){ assert(false);
        throw "Could not set reference indices, check point set for coplanarity.";
    }
}

ScalarType Functional::ReferenceIndicesConditionNumber(){
    Eigen::Matrix<ScalarType,Dimension,Dimension> m;
    for(int i=0; i<Dimension; ++i)
        m.col(i) = pts.col(referenceIndices[i]);
    return fabs(m.determinant())/ std::pow(m.norm(),Dimension);
}


#endif
