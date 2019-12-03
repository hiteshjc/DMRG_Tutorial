import numpy as N
import scipy as S
import sys
import pylab
import copy 

#############################################################################
# This DMRG code is for educational purposes only.
# It was written by Nitin Kaushal (U Tennessee and IISER K) and Hitesh Changlani (FSU and MagLab) 
# on December 2,2019 as a demonstration for the main ideas
# of the original DMRG algorithm by Steve White (PRL 1992)
# The code does DMRG for an open spin chain - edge spin 1/2 are placed irrespective
# of the bulk spin specified by the user 

############################################################################
def Local_Spin_Ops(TwoTimesSpin):

        #Natural Basis = {|Sz=S>=|0>, |Sz=S-1>=|1>,.... |Sz=S-i>=|i>...,|Sz=-S>=|2S>}
        #Here Operators in Matrix form for single site are defined for arbitary Spin
        #For Spin Model, that are Sz, S_minus=TransposeConjugate(S_plus)

        Spin=TwoTimesSpin*0.5
        N_Basis=int(TwoTimesSpin + 1)
        Sz_Fresh=N.zeros( (N_Basis, N_Basis), dtype=float )
        Splus_Fresh=N.zeros( (N_Basis, N_Basis), dtype=float )
        Sminus_Fresh=N.zeros( (N_Basis, N_Basis), dtype=float ) #Reduntant Matrix


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

def make_ham_of_superblock_given_sp_sm_sz(m,sp,sm,sz,spin, H_lb, H_rb):

        #SUPERBLOCK = [LEFT BLOCK]  +  [RIGHTBLOCK] + [L-site1] + [site2-R] + [site1-site2]
        #SYSTEM = [LEFT BLOCK] + [L-site1]
        #ENVIROMENT = [site2-R] + [RIGHTBLOCK]
        #[SITE1], and [SITE2] are fresh sites

        d=int(2*spin + 1 + 1.0e-6)  # Size of the free spin Hilbert space 
        H=N.zeros( (m*d*d*m, m*d*d*m), dtype=float )

        spsite, smsite, szsite = Local_Spin_Ops(int(2*spin + 1.e-6))


        # Term Leftblock = (H_lb)mxm (X) (Iden)dxd (X) (Iden)mxm (X) (Iden)dxd ; (X) is direct product 
        for ml in range(m):      # Left block
                for mlprime in range(m):
                        for sl in range(d): # site1
                                slprime=sl
                                for sr in range(d): #site2
                                        for mr in range(m): # Right block
                                                cind1=ml*m*d*d + sl*d*m + mr*d + sr
                                                mrprime=mr  # Delta function
                                                srprime=sr  # Delta function
                                                cind2=mlprime*m*d*d + slprime*d*m + mrprime*d + srprime
                                                H[cind1,cind2]+=H_lb[ml,mlprime]
                                                    

        # Term RightBlock = (Iden)mxm (X) (Iden)dxd (X) (H_rb)mxm (X) (Iden)dxd
        for mr in range(m):      # Right block 
                for mrprime in range(m):
                        for sr in range(d): # Right spin
                                srprime=sr
                                for ml in range(m):
                                        for sl in range(d):
                                                cind1=ml*m*d*d + sl*d*m + mr*d + sr
                                                mlprime=ml  # Delta function
                                                slprime=sl  # Delta function    
                                                cind2=mlprime*m*d*d + slprime*d*m + mrprime*d + srprime
                                                H[cind1,cind2]+=H_rb[mr,mrprime] #Reflection Symmetry can be used by using H_lb instead of H_rb


	
	# Term L-site1 = (Opr1)mxm (X) (Opr2)dxd (X) (Iden)mxm (X) (Iden)dxd 
	for ml in range(m):      # Left block 
		for mlprime in range(m): 
			for sl in range(d): # site1 
				for slprime in range(d):
                                        for sr in range(d): #site2
					        for mr in range(m): # Right block 
							cind1=ml*m*d*d + sl*d*m + mr*d + sr
							mrprime=mr  # Delta function
							srprime=sr  # Delta function    
							cind2=mlprime*m*d*d + slprime*d*m + mrprime*d + srprime
							H[cind1,cind2]+=0.5*sp[ml,mlprime]*smsite[sl,slprime]
							H[cind1,cind2]+=0.5*sm[ml,mlprime]*spsite[sl,slprime]
							H[cind1,cind2]+=sz[ml,mlprime]*szsite[sl,slprime]
	
	# Term site1-site2 = (Iden)mxm (X) (Opr1)dxd (X) (Iden)mxm (X) (Opr2)dxd 
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
	
	# Term site2-R = (Iden)mxm (X) (Iden)dxd (X) (Opr2)mxm (X) (Opr1)dxd
	for mr in range(m):      # Right block 
		for mrprime in range(m): 
			for sr in range(d): # Right spin
				for srprime in range(d): 
					for ml in range(m):
						for sl in range(d): 
							cind1=ml*m*d*d + sl*d*m + mr*d + sr
							mlprime=ml  # Delta function
							slprime=sl  # Delta function    
							cind2=mlprime*m*d*d + slprime*d*m + mrprime*d + srprime

                                                        #Reflection symmetry is used here
							H[cind1,cind2]+=0.5*sp[mr,mrprime]*smsite[sr,srprime] 
							H[cind1,cind2]+=0.5*sm[mr,mrprime]*spsite[sr,srprime]
							H[cind1,cind2]+=sz[mr,mrprime]*szsite[sr,srprime]




        #Creating System, Enviroment Hamil--##############################################################################

        H_system=N.zeros( (m*d, m*d), dtype=float )

        # Term Leftblock = (H_lb)mxm (X) (Iden)dxd
        for ml in range(m):      # Left block
                for mlprime in range(m):
                        for sl in range(d): # site1
                                slprime=sl
                                cind1=ml*d + sl
                                cind2=mlprime*d + slprime
                                H_system[cind1,cind2]+=H_lb[ml,mlprime]




        # Term L-site1 = (Opr1)mxm (X) (Opr2)dxd
        for ml in range(m):      # Left block
                for mlprime in range(m):
                        for sl in range(d): # site1
                                for slprime in range(d):
                                        cind1=ml*d + sl
                                        cind2=mlprime*d + slprime
                                        H_system[cind1,cind2]+=0.5*sp[ml,mlprime]*smsite[sl,slprime]
                                        H_system[cind1,cind2]+=0.5*sm[ml,mlprime]*spsite[sl,slprime]
                                        H_system[cind1,cind2]+=sz[ml,mlprime]*szsite[sl,slprime]

        H_enviroment=copy.deepcopy(H_system)
        ###########################################################################################################
        

        return H, H_system, H_enviroment
############################################################################

def get_gs_of_superblock_H_and_make_its_dm(m,spin,H):
	# print GS and Excited state energy
	d=int(2*spin + 1 + 1.0e-6)  # Size of the free spin Hilbert space 
	eigs,vecs=N.linalg.eigh(H) # Doing full diag here, should be replaced by Lanczos
	gs=vecs[:,0] # Gets the zeroth column which is the ground state
        exc1=vecs[:,1] #1st excited state
        weight=1.0 # weight=1.0 will only use GS vector i.e. pure ensemble, weight<1 corresponds to mixed ensemble of GS and 1st Excited state

	dm=N.zeros((m*d,m*d),dtype=float)
	
        #Reduced_Density_Matrix=dm=Trace_{sr,mr}(|GS><GS|)

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
							dm[ind1,ind2]+=weight*gs[cind1]*gs[cind2] + (1.0 - weight)*exc1[cind1]*exc1[cind2]
	
	return eigs[0],eigs[1],dm

############################################################################
def diagonalize_dm_truncate_and_find_new_matrices(dm,maxm,m_old,H_system):		

        eigs_dm,vecs_dm=N.linalg.eigh(dm) #Eigenvalues are in increasing order
       #print (eigs_dm)

        spsite, smsite, szsite = Local_Spin_Ops(int(2*spin + 1.e-6))
        d=int(2*spin + 1 + 1.0e-6)
        m=min(maxm,d*m_old) #Exact until maxm>d*m_old, m is the new no. of states kept

        sz=N.zeros( (m, m), dtype=float )
        sp=N.zeros( (m, m), dtype=float )
        sm=N.zeros( (m, m), dtype=float )
        H_lb=N.zeros( (m, m), dtype=float )



        #<ml|O|mlprime>=\sum <ml|ml_old sl > <s|O|sl prime> <ml_old sprime|mlprime> 
        for ml in range(m):
                for mlprime in range(m):
                        for ml_old in range(m_old):
                                for sl in range(d):
                                        for slprime in range(d):
                                                mlprime_old=ml_old 
                                                ind1=ml_old*d + sl
                                                ind2=mlprime_old*d + slprime
                                                sz[ml,mlprime]+=vecs_dm[ind1,m_old*d-1-ml]*vecs_dm[ind2,m_old*d-1-mlprime]*szsite[sl,slprime]
                                                sp[ml,mlprime]+=vecs_dm[ind1,m_old*d-1-ml]*vecs_dm[ind2,m_old*d-1-mlprime]*spsite[sl,slprime]
                                                sm[ml,mlprime]+=vecs_dm[ind1,m_old*d-1-ml]*vecs_dm[ind2,m_old*d-1-mlprime]*smsite[sl,slprime]

           
        #<ml|H|mlprime>=\sum <ml|ml_old sl > <ml_old sl |H| ml_oldprime sl prime> <ml_oldprime sprime|mlprime> 
        for ml in range(m):
                for mlprime in range(m):
                        for ml_old in range(m_old):
                                for mlprime_old in range(m_old):
                                        for sl in range(d):
                                        	for slprime in range(d):
							#slprime=sl
							ind1=ml_old*d + sl
							ind2=mlprime_old*d + slprime
							#H_lb[ml,mlprime]+=vecs_dm[ind1,m_old*d-1-ml]*vecs_dm[ind2,m_old*d-1-mlprime]*H_system[ml_old,mlprime_old]
							H_lb[ml,mlprime]+=vecs_dm[ind1,m_old*d-1-ml]*vecs_dm[ind2,m_old*d-1-mlprime]*H_system[ind1,ind2]




        H_rb=copy.deepcopy(H_lb) #Reflection Symmetry

        truncerror=1.0
        for n in range(m):
                truncerror-=eigs_dm[m_old*d-1-n]

	return m,truncerror,sp,sm,sz,H_lb,H_rb

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
#Initialize H_lb, H_rb
H_lb=N.zeros( (m, m), dtype=float )
H_rb=N.zeros( (m, m), dtype=float )
########################################################

# Do DMRG
print "#iterno   SytemLength  E0  E1  m_kept  Truncerror  E0/site  Gap=E1-E0"
for i in range(niter):
	H,H_system,H_enviroment=make_ham_of_superblock_given_sp_sm_sz(m,sp,sm,sz,spin, H_lb, H_rb)	
	e0,e1,dm=get_gs_of_superblock_H_and_make_its_dm(m,spin,H)
        print (i),(4+2*i),(e0),(e1),
	m_old=m
	m,truncerror,sp,sm,sz,H_lb,H_rb=diagonalize_dm_truncate_and_find_new_matrices(dm,maxm,m_old,H_system)
        print (m),(truncerror),((e0)/float(4+2*i)),(e1-e0)
        
        
