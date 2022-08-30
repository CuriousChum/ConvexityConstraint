from agd.Eikonal import FileIO,GetBinaryDir
import numpy as np

Flags = {
	'LogSum':2**0,
	'LogGrad':2**1,
	'LogHessian':2**2,
	'Values':2**3,
	'Jacobian':2**4,
	'Hessian':2**5
}

def Measures(points,heights,exclude,order=1):
	"""
	Computes the lebesgue measure of the subgradient cells of the non-excluded points, 
	and their derivatives w.r.t. the heights.
	Input:
	 - points (array of shape (N,d)) : d=1,2,3 is the dimension, N the number of points
	 - heights (array of shape (N,)) : the values of the function whose subgradients are requested.
	 - exclude (boolean array of shape (N,)) : points where subgradients are not to be computed.
	 - order : differentiation order

	Output : 
	- values : the subgradient areas
	- jacobian_v, jacobian_i, jacobian_j : if order>=1, 
		the jacobian of the subgradients as sparse matrix
	- hessian_v, hessian_i, hessian_j, hessian_k : if order>=2, 
		the hessian of the subgradients, as a sparse tensor
	"""
	assert 0<=order<=2
	flag=Flags['Values']
	if order>=1: flag+=Flags['Jacobian']
	if order>=2: flag+=Flags['Hessian']
	return _Compute(points,heights,exclude,flag)

def Barrier(points,heights,exclude,order=2):
	"""
	Computes the barrier function defined as the sum of the logarithms of the measures 
	of the subgradient cells of the non-excluded points.
	Input : see Measures function
	Output : 
	- LogSum : the value of the barrier function
	- logGrad : the gradient of the barrier function
	- logHessian_v, logHessian_i, logHessian_j : the sparse hessian of the barrier function
	Failure : 
	- retcode : 1, is returned
	- an explanatory error message is displayed
	"""
	assert 0<=order<=2
	flag = Flags['LogSum']
	if order>=1: flag+=Flags['LogGrad']
	if order>=2: flag+=Flags['LogHessian']
	return _Compute(points,heights,exclude,flag)

def _Compute(points,heights,exclude,flag,retnan=True):
	"""Calls a suitable binary to compute the data required by Measures and Barrier"""
	bin_dir = GetBinaryDir("FilePSG",None) 
	if points.ndim==1: points=points[None,:] # 1D special case
	size,dim = points.shape
	assert heights.shape == (size,)
	assert exclude.shape == (size,)
	psgIn = {"points":points,"heights":heights,"exclude":exclude,"flag":flag}
	psgOut = FileIO.WriteCallRead(psgIn, "FilePSG", bin_dir)
	for key,value in psgOut.items():
		if key[-2:] in ('_i','_j','_k'): psgOut[key]=value.astype(int)

	if psgOut.get('retcode',0):
		raise ValueError('Failed computing the sugradients areas')

	# Check for error, and return NaN in that case
	# if psgOut.get("retcode",0):
	# 	anan = np.full_like(heights,fill_value=np.nan)
	# 	afloat = np.array([],dtype=np.float)
	# 	aint = np.array([],dtype=np.int)

		# if flag & Flags['LogSum']: psgOut['logSum'] =np.nan
		# if flag & Flags['LogGrad']:psgOut['logGrad']=anan
		# if flag & Flags['LogHessian']:
		# 	psgOut.update({'logHessian_v':afloat, 'logHessian_i':aint, 'logHessian_j':aint})
		# if flag & Flags['Values']: psgOut['values'] =anan
		# if flag & Flags['Jacobian']:
		# 	psgOut.update({'jacobian_v':afloat, 'jacobian_i':aint, 'jacobian_j':aint})
		# if flag & Flags['Hessian']:
		# 	psgOut.update({'hessian_v':afloat, 'hessian_i':aint, 'hessian_j':aint, 'hessian_k':aint})

	return psgOut
