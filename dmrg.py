import numpy as N
import scipy as S
import pylab

############################################################################
def Local_Spin_Ops(TwoTimesSpin):

        #Natural Basis = {|Sz=S>=|0>, |Sz=S-1>=|1>,.... |Sz=S-i>=|i>...,|Sz=-S>=|2S>}
        #Here Operators in Matrix form for single site are defined for arbitary Spin
        #For Spin Model, that are Sz, S_minus=TransposeConjugate(S_plus)

        Spin=TwoTimesSpin*0.5
        N_Basis=int(TwoTimesSpin + 1)
        Sz_Fresh=N.zeros( (N_Basis, N_Basis), dtype=complex )
        Splus_Fresh=N.zeros( (N_Basis, N_Basis), dtype=complex )
        Sminus_Fresh=N.zeros( (N_Basis, N_Basis), dtype=complex ) #Reduntant Matrix


        #Sz is a diagonal Matrix :)
        for i in range(N_Basis):
                Sz_Fresh[i,i]=(Spin)-i


        #<nz|Splus|mz> = delta_{nz,mz+1} sqrt(S(S+1) - nz*mz); nz=S-n,mz=S-m
        #---> <n|Splus|m> = delta_{n,m-1}sqrt (S(S+1) - (S-n)*(S-m))
        for n in range(0,TwoTimesSpin):
                m=n+1
                Splus_Fresh[n,m]=N.sqrt( Spin*(Spin+1) - (Spin-n)*(Spin-m))
                Sminus_Fresh[m,n]=N.conjugate(Splus_Fresh[n,m])


        return Splus_Fresh,Sminus_Fresh,Sz_Fresh
############################################################################

def make_ham_of_superblock_given_sp_sm_sz(m,sp,sm,sz,spin):
	spsite, smsite, szsite = Local_Spin_Ops(int(2*spin + 1.e-6))
	# Term L-site1
	for ml in range(m):      # Left block 
		for mlprime in range(m): 
				cind1=
				H[cind1,cind2]+=
	
	# Term site1-site2
				cind1=
				H[cind1,cind2]+=
	
	# Term site2-R
	for mr in range(m):      # Right block 
		for mrprime in range(m):

	return H
############################################################################

def get_gs_of_superblock_H_and_make_its_dm(H):
	# print GS and Excited state energy
	eigs,vecs=N.linalg.eigh(H)

	return eigs[0],eigs[1]

############################################################################
def diagonalize_dm_truncate_and_find_new_spin_matrices(dm):		

	return m,sp,sm,sz

#########################################################
# User specifies spin
spin=float(sys.argv[1]) 
niter=int(sys.argv[2]) 
maxm=int(sys.argv[3]) 
#########################################################

########################################################
#Initialize sp,sm,sz
# Irrespective of spin 1/2 or spin 1 Hamiltonian, the boundary spins should be spin 1/2 to start off with
m=2 
sp, sm, sz = Local_Spin_Ops(int(2*0.5 + 1.e-6))
########################################################

# Do DMRG
for i in range(niter):
	H=make_ham_of_superblock_given_sp_sm_sz(m,sp,sm,sz,spin)	
	e0,e1,dm=get_gs_of_superblock_H_and_make_its_dm(H)
	m,sp,sm,sz=diagonalize_dm_truncate_and_find_new_spin_matrices(dm,maxm)
