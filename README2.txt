This gives a rundown of the code written thus far, up to 5th May 2017:

In the first README(please).txt document you will find information on running the piece of code that calculates
the Corelation and Variance matrices for the corelation for the remote quadrupole between some n number of 
galaxies across the skies. 

We are now interested in writing routines that let us also calculate, (a) the corelation between local multipoles,
i.e. the Cl's, and (b) the cross-corelation between a local multipole and a remote quadrupole.

As of now, these routines have not been subjected to extensive testing (they have passed some, but not all). 

The code:

(1) Creating corelations between remote quadrupoles: See README(please).txt

(2) TransFuncGen: The previous transfer function that was being used had the issue that it only calculated
values for l=2. This is the generalized version that let's us calculate it for any arbitrary l. Should test this.

	nlist and qlist: Generated using OmegaM = 0.32
	nlist1 and qlist1: Generated using OmegaM=0.30

(3) LocalPS.py: Generates the corelation between the local multipoles, i.e. the Cl's, for n clusters. 
Note that it does so analytically and not as a limiting case of (4) below. (although that is a good test). Also note that this routine only generates the diagonal elements for the matrix containing Cl’s. They’re actually put together in the matrix form in COVandREL.py (see method “def LocalCMB”). 

(3) VarMatrix_Gen_Local.py: Generates the cross corelation between the local multipole and remote quadrupole
for n clusters. 

(4) VarMatrix_Generalized.py: Generates remote quadrupole corelations (see README(please).txt)

(5) COVandREL.py: Father code of all. Combines all the above pieces of information to return the Variance and 
Corelation matrices. 

The last piece of code is what we want to mess around with. All parameters are passed to COVandREL.py
which then calls all the other classes as per the parameters specified. 

Specifically, we want to give this:
	(i) "ranClusters": an nx3 matrix containing cluster locations
	(ii) "n" : The number of clusters from randClusters we actually want to consider. n=>2 (for remote quadrupole)
	(iii) "l" : This is really lmax, i.e. the maximum l value upto which we want the local power spectrum
		and the cross-corelations. l=>2

Note: The CrossCMB part is extremely slow and takes a long time to calculate the desired matrices for COVandREL.py.
It is therefore desireable to first generate the matrix from VarMatrix_Gen_Local.py and store that so that
it can be reused in the future. 

Testing: 
Thus far we have passed most of the tests tried, but failed a crucial one. Specifically, an important test that
this routine has passed is that the Variance Matrix should be hermitian and indeed it is. Unfortunately, the 
relation matrix is not symmetric as we would expect it to be..... this is quite inconsistent with the last
result though…..

	Update: The above is no longer true. For some reason, we are actually getting a symmetric relation matrix. 
	There is however one problem: The diagonals are not all zero…. Namely there are three diagonal elements which are non-	zero. Wait, why do we even need that?? Shouldn’t we expect elements of Cl’s to appear in the diagonals? That’s exactly 	what happens actually…..

	OH I REMEMBER THE PROBLEM NOW. IT WAS NEVER SYMMETRICITY. The problem was with the real valued gamma matrix which 	wasn’t hermitian…..My assumption is that this problem is being caused by the three non-zero diagonal elements in the 	(symmetric) relation matrix. Note that GammaReal is still symmetric, just not positive definite. The diagonal elements 	are, in general, non-zero.	
