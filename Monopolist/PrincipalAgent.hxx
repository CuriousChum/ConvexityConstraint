//
//  PrincipalAgent.hxx
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 29/01/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_PrincipalAgent_hxx
#define CGalTest_PrincipalAgent_hxx

void PrincipalAgent::MakeObjective(){
    terms.constant=0;
    terms.linear.clear();
    terms.linear.resize(rt.number_of_vertices(),0.);
    terms.quadratic.clear();
    hess_.resize(0,0);
    
    using namespace CGAL_Traits;
    
    // Enumerate triangles
    
    for(auto fIt=rt.finite_faces_begin(), fEnd=rt.finite_faces_end();
        fIt!=fEnd; ++fIt){
        
        // Eliminate flat triangles
        if(fIt->vertex(0)->info().boundary &
           fIt->vertex(1)->info().boundary &
           fIt->vertex(2)->info().boundary) continue;
        
        const Point barycenter = CGAL::barycenter(fIt->vertex(0)->point(), 1./3.,
                                                  fIt->vertex(1)->point(), 1./3.,
                                                  fIt->vertex(2)->point());
        const Point origin(0.,0.);
        std::array<Vector,3> g; // gradients of basis functions
        for(int i=0; i<3; ++i){
            Vertex_handle
            p = fIt->vertex(rt.cw(i)),
            q = fIt->vertex(rt.ccw(i));
            
            g[i]=p->point() - q->point();
//            g[i] = fIt->vextex(rt.ccw(i))->point()-fIt->vextex(rt.cw(i))->point();
        }
        const ScalarType area2 = CGAL::determinant(g[0],g[1]);
        for(int i=0; i<3; ++i)
            g[i] = g[i].perpendicular(CGAL::POSITIVE) / area2;
        
        const ScalarType area = area2/2;
        
        
        for(int i=0; i<3; ++i){
            terms.linear[fIt->vertex(i)->info().index] +=
            area*(1./3. - g[i]*(barycenter-origin)); // int u - <grad u, x> ??
        }
        
        if(hasQuadraticCost){
            for(int i=0; i<3; ++i)
                for(int j=0; j<3; ++j)
                    terms.quadratic.push_back({
                        fIt->vertex(i)->info().index,
                        fIt->vertex(j)->info().index,
                        area*g[i]*g[j]}); // int (1/2) * \|grad u\|^2
        }
    }
    Simplify(terms.quadratic);
}

ScalarType PrincipalAgent::Value(const NS::VectorType & x){
    
    size_t n=terms.linear.size();
    
    if(hess_.size()==0){
        if(terms.quadratic.empty()) MakeObjective();
        n=terms.linear.size();
        hess_.resize((int)n,(int)n);
        std::vector<Eigen::Triplet<ScalarType> > triplets;
        triplets.reserve(terms.quadratic.size());
        for(auto c : terms.quadratic)
            triplets.push_back({c.i,c.j,c.v});
        hess_.setFromTriplets(triplets.begin(),triplets.end());
    }
    
    NS::VectorType l_(n);
    for(int i=0; i<n; ++i) l_[i]=terms.linear[i];
    
    const NS::VectorType hx = hess_*x;
    grad_=l_ + hx;
    return terms.constant+x.dot(l_+0.5*hx);
}

void PrincipalAgent::ObjectiveType::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(constant)
    ExportArrayArrow(linear)
    ExportArrayArrow(quadratic)
    << "}";
}

void PrincipalAgent::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(rt)
    ExportVarArrow(terms)
    << "}";
}

// %%%%%%%%%%%%%%%%%%%% Error indicators, refinement %%%%%%%%%%%%%%%%%%
CGT::Vector PrincipalAgent::Gradient(CGT::Face_handle f, const NS::VectorType & x) const {
    std::array<CGT::Vector,3> edges;
    for(int i=0; i<3; ++i)
        edges[i] = f->vertex(CGT::RT::ccw(i))->point() - f->vertex(CGT::RT::cw(i))->point();
    
    CGT::Vector result(0,0);
    for(int i=0; i<3; ++i)
        result = result+ x[f->vertex(i)->info().index] * edges[i];
    
    return result.perpendicular(CGAL::POSITIVE) / CGAL::determinant(edges[0],edges[1]);
}


ScalarType PrincipalAgent::ErrorIndicator(CGT::Face_handle f, const NS::VectorType & x) const {
    std::array<CGT::Vector,3> edges;
    for(int i=0; i<3; ++i)
        edges[i] = f->vertex(CGT::RT::ccw(i))->point() - f->vertex(CGT::RT::cw(i))->point();
    std::array<ScalarType,3> edgeLengths;
    for(int i=0; i<3; ++i)
        edgeLengths[i] = sqrt(edges[i].squared_length());
    
    // triangle diameter
    const ScalarType h = std::max(edgeLengths[0], std::max(edgeLengths[1],edgeLengths[2]));
    
    const CGT::Vector g = Gradient(f,x);
    ScalarType result=0;
    
    for(int i=0; i<3; ++i)
        if(f->mirror_vertex(i)->info().index != InfiniteIndex)
            result += edgeLengths[i]*(Gradient(f->neighbor(i),x)-g).squared_length();
    return sqrt(h*result);
}

void PrincipalAgent::Refine(ScalarType quantile, const NS::VectorType & x) {
    std::map<ScalarType,CGT::Face_handle,std::greater<ScalarType> > indicators;
    for(CGT::Finite_faces_iterator fIt = rt.finite_faces_begin();
        fIt != rt.finite_faces_end(); ++fIt)
        indicators.insert({ErrorIndicator(fIt,x),fIt});
    
    int counter=0, counterMax = 1+int(indicators.size()*quantile);
    std::vector<CGT::Face_handle> facesToRefine; facesToRefine.reserve(counterMax);
    for(auto it = indicators.begin();
        it!=indicators.end() && counter<counterMax;
        ++it, ++counter)
        facesToRefine.push_back(it->second);
    
    RefineFaces(rt,facesToRefine);
    MakeObjective();
}


std::vector<CGT::Full_point> PrincipalAgent::Pts() const {
    std::vector<CGT::Full_point> pts; pts.reserve(rt.number_of_vertices());
    for(CGT::Finite_vertices_iterator vIt=rt.finite_vertices_begin();
        vIt!=rt.finite_vertices_end(); ++vIt)
        pts.push_back({vIt->point(),vIt->info()});
    std::sort(pts.begin(),pts.end(),[](CGT::Full_point p, CGT::Full_point q){
        return p.second.index < q.second.index;});
    return pts;
}

#endif
