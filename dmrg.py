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
	d=int(2*spin + 1 + 1.0e-6)  # Size of the free spin Hilbert space 
	
	# Term L-site1
	for ml in range(m):      # Left block 
		for mlprime in range(m): 
			for sl in range(d): 
				for slprime in range(d): 
					for mr in range(m):
						for sr in range(d): 
							cind1=ml*m*d*d + sl*d*m + mr*d + sr
							mrprime=mr  # Delta function
							srprime=sr  # Delta function    
							cind2=mlprime*m*d*d + slprime*d*m + mrprime*d + srprime
							H[cind1,cind2]+=0.5*sp[ml,mlprime]*smsite[sl,slprime]
							H[cind1,cind2]+=0.5*sm[ml,mlprime]*spsite[sl,slprime]
							H[cind1,cind2]+=sz[ml,mlprime]*szsite[sl,slprime]
	
	# Term site1-site2
	for sl in range(d):      # Left spin  
		for slprime in range(d): 
			for sr in range(d):  # Right spin
				for srprime in range(d): 
					for ml in range(m):
						for mr in range(m): 
							cind1=ml*m*d*d + sl*d*m + mr*d + sr
							mlprime=ml  # Delta function
							mrprime=mr # Delta function    
							cind2=mlprime*m*d*d + slprime*d*m + mrprime*d + srprime
							H[cind1,cind2]+=0.5*spsite[sl,slprime]*smsite[sr,srprime]
							H[cind1,cind2]+=0.5*smsite[sl,slprime]*spsite[sr,srprime]
							H[cind1,cind2]+=szsite[sl,slprime]*szsite[sr,srprime]
	
	# Term site2-R
	for mr in range(m):      # Right block 
		for mrprime in range(m): 
			for sr in range(d): 
				for srprime in range(d): 
					for ml in range(m):
						for sl in range(d): 
							cind1=ml*m*d*d + sl*d*m + mr*d + sr
							mlprime=ml  # Delta function
							slprime=sl  # Delta function    
							cind2=mlprime*m*d*d + slprime*d*m + mrprime*d + srprime
							H[cind1,cind2]+=0.5*sp[mr,mrprime]*smsite[sr,srprime]
							H[cind1,cind2]+=0.5*sm[mr,mrprime]*spsite[sr,srprime]
							H[cind1,cind2]+=sz[mr,mrprime]*szsite[sr,srprime]

	return H
############################################################################

def get_gs_of_superblock_H_and_make_its_dm(m,spin,H):
	# print GS and Excited state energy
	d=int(2*spin + 1 + 1.0e-6)  # Size of the free spin Hilbert space 
	eigs,vecs=N.linalg.eigh(H) # Doing full diag here, should be replaced by Lanczos
	gs=vecs[:,0] # Gets the zeroth column which is the ground state
	dm=N.zeros((m*d,m*d),dtype=float)
	
	for ml in range(m):
		for sl in range(d):
			ind1=ml*d + sl
			for mlprime in range(m):
				for slprime in range(d):
					ind2=mlprime*d + slprime
					for mr in range(m):
						for sr in range(d):
							cind1=ml*m*d*d + sl*d*m + mr*d + sr
							cind2=mlprime*m*d*d + slprime*d*m + mr*d + sr
							dm[ind1,ind2]+=gs[cind1]*gs[cind2]
	
	return eigs[0],eigs[1],dm

############################################################################
def diagonalize_dm_truncate_and_find_new_spin_matrices(dm,maxm):		

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
	e0,e1,dm=get_gs_of_superblock_H_and_make_its_dm(m,spin,H)
	m,sp,sm,sz=diagonalize_dm_truncate_and_find_new_spin_matrices(dm,maxm)
