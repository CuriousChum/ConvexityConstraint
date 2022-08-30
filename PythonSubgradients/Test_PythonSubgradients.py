import sys; sys.path.insert(0,"../../AdaptiveGridDiscretizations")
import numpy as np
import PythonSubgradients
from pysdot.domain_types import ConvexPolyhedraAssembly
from pysdot import PowerDiagram
from scipy.sparse import csr_matrix



aX = np.linspace(0,1,5)
Xint = np.meshgrid(aX,aX,indexing='ij')
Xint = np.reshape(Xint,(2,-1))

def BoundingBox(X,delta=1):
	"""A bounding box with some extra space"""
	Xmax = np.max(X,axis=1)
	Xmin = np.min(X,axis=1)
	Xdiff = Xmax-Xmin
	return [Xmin-delta*Xdiff,Xmax+delta*Xdiff]

(x0min,x1min),(x0max,x1max) = BoundingBox(Xint)
Xbd = np.array([(x0min,x1min),(x0min,x1max),(x0max,x1min),(x0max,x1max)]).T
X = np.concatenate([Xint,Xbd],axis=-1)

heights = 0.5*np.sum(X**2,axis=0)
exclude = np.full(X.shape[1],False)
exclude[-Xbd.shape[1]:]=True

print(X.T.shape,heights.shape,exclude.shape)
print(X.T)
out = PythonSubgradients.Measures(X.T,heights,exclude,order=1)

#----------------- Pysdot version -----------

K = ConvexPolyhedraAssembly()
r=10
K.add_box([-r,-r],[r,r])

def subgradients_diagram(φ,X,K):
    Y = np.moveaxis(X,0,-1).reshape(-1,2)
    ψ = (np.sum(X**2,axis=0) - 2*φ).reshape(-1)
    Y,ψ = (np.ascontiguousarray(e) for e in (Y,ψ)) # Important !
    return PowerDiagram(Y,ψ,K)

pd = subgradients_diagram(heights,X,K)
assert np.allclose(pd.integrals()[:-4],out['values'][:-4])



mvs = pd.der_integrals_wrt_weights()
Jac_pysdot = csr_matrix((-2*mvs.m_values,mvs.m_columns, mvs.m_offsets), shape=(mvs.v_values.size,)*2)

print(out['jacobian_i'])
Jac_cgal = csr_matrix((out['jacobian_v'],(out['jacobian_i'],out['jacobian_j'])))

print(Jac_cgal.shape)
print(Jac_cgal - Jac_pysdot[:-4])
