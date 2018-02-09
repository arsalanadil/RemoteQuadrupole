Here's a rundown of what each file contains/does:

i) RandomClusterGenerator.py: Generates random locations (in spherical co-ordinates) of a given
number of clusters. You want to change the parameters:

	n: specifies the number of clusters you want
	(line 23): np.save("filename", rand"n"clusters)

I've included one such file that I've been using for testing purposes, "rand50Clusters.npy"
You can generate more of these pretty easily.

ii) ClebGord.py: No surprises here. This just calculates the Clebsch-Gordan co-efficients;

iii)TransferFunction.py: Calculates transfer function; Note that this is not the generalized version, 
i.e. it only works for l =2 (I actually do have the generalized version for this but I haven't tested it 
extensively). This depends on the files "nlist.csv" and "qlist.csv" from Dave's code for the conversion
factors.

iv) VarMatrix_Generalized.py : This is the heart of this code. It calculates, first the
corelations between two galaxy clusters in their respective axes, and then applies a 
Wigner rotation to rotate them in the frame of our z-axis.
The final product is a matrix element (actually a matrix but we end up using only one
of the elements of this matrix) <a2m1(r1) a2m2*(r2)>, i.e. the corelation in the CMB 
quadrupole.
This routine depends on (ii) and (iii) above - for obvious reasons.

ii) VarAndCorMat.py : This is the code that you would actually be running; that is to
say that this is the routine which is dependent on everything else in the directory and
you just want to run this to get the output. Based on VarMatrixGeneralized.py, this 
routine calculates the Variance and Corelation Matrices and then combines them 
(see lines 97-100 in code) to give you the real valued Covariance Matrix. 

A couple of things about this routine -- 
This is where all your specifications (input parameters go). Specifically you want to
give this:
	a) the power spectrum (def powerSpec). Currently set to 1/k^3.... presumablyt that's the one
	you'd want to use for the standard mode of operation.
	
	b)the variable "randClusters" (line 24): import an nx3 matrix containing a randomly generated list of 
	clusters using the routine "RandomClusterGenerator.py"
	
	c) the variable "n" (line 25): the number of clusters you actually want to use from the file
	you just loaded. This is useful in cases where say I have a file containing 50 randomly generated
	clusters but I want to test something so I can just change this value to consider, say 10 out of 
	the 50 clusters to save computational time. 

	d) line 105: "np.save("filnename", VariableToSave)"; this is pretty self explanatory - just the name
	of the file you want to save your real valued Covariance Matrix as. You probably want to use this
	considering it takes a while for all of this to come up with an output and you presumably want to
	save it somewhere to play around with it

In summary, if you want to compare just two galaxy clusters, you would primarily be concerned with 
VarMatrixGeneralized.py; Note that in this case you would want to specify the power spectrum to this
routine (as it is normally supplied via VarAndCorMat.py)

If you want to compare several - you would want to generate a random list using RandomClusterGenerator.py,
save the output from that, and then load it in to VarAndCorMat.py and make 
the appropriate changes to the variables I mentioned above.

I've also included some figures that you might find helpful; I think those are pretty self explanatory.

"Beauty is the first test.... there is no permanent place for ugly mathematics in this world"
-G.H. Hardy
