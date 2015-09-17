
def left_lead(Energy, H00, VLR, VRL):
	g00 = invert((Energy+i*eta)*I - H00) 
	t0 = np.dot(g00, VRL)
	s0 = np.dot(g00, VLR)
	temp = np.identity((2*size+1))
	s_sum = np.identity(2*size+1)
	T = np.zeros(2*size+1)
	while np.sum(np.abs(T - temp)) > tolerance:
		temp = T
		T = T + np.dot(s_sum,t0)
		M = invert(np.identity(2*size+1) - np.dot(t0,s0) - np.dot(s0,t0))
		t0 = np.dot(M , np.dot(t0,t0))
		s_sum = np.dot(s_sum,s0) 
		s0 = np.dot(M, np.dot(s0,s0))
	return np.dot(invert(np.identity(2*size+1) - np.dot(g00,np.dot(VLR,T))),g00)

def right_lead(Energy, H00, VLR, VRL):
	g00 = invert((Energy+i*eta)*I - H00) 
	t0 = np.dot(g00, VLR)
	s0 = np.dot(g00, VRL)
	temp = np.identity((2*size+1))
	s_sum = np.identity(2*size+1)
	T = np.zeros(2*size+1)
	while np.sum(np.abs(T - temp)) > tolerance:
		temp = T
		T = T + np.dot(s_sum,t0)
		M = invert(np.identity(2*size+1) - np.dot(t0,s0) - np.dot(s0,t0))
		t0 = np.dot(M , np.dot(t0,t0))
		s_sum = np.dot(s_sum,s0) 
		s0 = np.dot(M, np.dot(s0,s0))
	return np.dot(invert(np.identity(2*size+1) - np.dot(g00,np.dot(VRL,T))),g00)


def left_lead_recursive(Energy, H00):
	g00 = invert((Energy+i*eta)*I - H00) 
	SL = g00
	SLold = np.zeros(shape=(2*size+1,2*size+1))
	#Build Left Lead
	while np.abs(np.sum(SL - SLold)) > tolerance:
		SLold = np.array(SL)
		SL = np.dot( invert(I - np.dot(np.dot(g00,VRL),np.dot(SLold,VLR) )) , g00)
	return SL

def right_lead_recursive(Energy, H00):
	g00 = invert((Energy+i*eta)*I - H00) 
	SR = g00
	SRold = np.zeros(shape=(2*size+1,2*size+1))
	#Right Lead	
	while np.abs(np.sum(SR - SRold)) > tolerance:
		SRold = np.array(SR)
		SR = np.dot( invert(I - np.dot(np.dot(g00,VLR),np.dot(SRold,VRL) )) , g00)
	return SR

def H_build(t1,t2):
    """
    8AGNR Labelling scheme:

    0 \_/ 8
      / \
      \_/
      / \
    7 \_/ 15

    """
    H = np.zeros(shape=(2*size,2*size))
    for x in xrange(size-1): H[x][x+1] = t2
    for x in xrange(size,2*size-1): H[x][x+1] = t2
    for x in xrange(1,size,2): H[x][x+size] = t1
    H = H + np.transpose(H)
    Hbigger = np.zeros(shape=(2*size+1,2*size+1)) 
    Hbigger[:2*size,:2*size] = H
    return Hbigger
	
def VRL_build(t1,t2):
	VRL = np.zeros(shape=(2*size,2*size))
	for x in xrange(0,size,2): VRL[x][x+size] = t1
	VRLbigger = np.zeros(shape=(2*size+1,2*size+1))
	VRLbigger[:2*size,:2*size] = VRL
	return VRLbigger

import numpy as np
from scipy.linalg import inv as invert
import shelve
import os

global size, tolerance, eta, i , t, ea, tau, lambda_a, lambda_c , I

def setValues( strain = 0, angle = 0):
    """  
    Sets Global values to be used in all calculations.
    Strain = 0.0 - 0.2 for realistic strain.
    Angle: 0 = Zigzag direction, pi/2 = Armchair direction, only works for these directions!!
    connex is the hopping parameter tau between the impurity and the substrate.
    Fermi energy allows the fermi energy to be set.
    pt prints some of these values to the screen.
    No return value.
    """  
    # Strain matrix
    sig = 0.165
    E11 = strain*(np.cos(angle)**2 - sig*np.sin(angle)**2)
    E12 = strain*((1+sig)*np.cos(angle)*np.sin(angle))
    E22 = strain*(np.sin(angle)**2 - sig*np.cos(angle)**2)
    # New distances
    l1 = 1 + E22
    l2 = 1 + 0.75*E11 - (np.sqrt(3.0)/2)*E12 + 0.25*E22
    # New hoppings
    t1 = -np.exp(-3.37*(l1-1)) 
    t2 = -np.exp(-3.37*(l2-1))
    return t1, t2

def check(M):
    for a in range(len(M)):
        for b in range(len(M)):
            if M[a][b] != 0.0: print a , b , M[a][b]



tolerance = 1.0e-8
eta = 1.0e-4
i = 1j

t = -1.0

#f = open("ldos_6agnr.dat","w")

def conductance_pristine(Energy, H00,VLR,VRL,SL,SR):

#       Giir = np.dot(invert( I - np.dot(np.dot(SL,VLR),np.dot(SR,VRL))) , SL)          
#       Gjjr = np.dot(invert( I - np.dot(np.dot(SR,VRL),np.dot(SL,VLR))) , SR)
#       Gjir = np.dot(np.dot(SR,VRL),Giir)

        Giir = np.dot(invert( I - np.dot(np.dot(SL,VLR),np.dot(SR,VRL)), True, False) , SL)
        Gjjr = np.dot(invert( I - np.dot(np.dot(SR,VRL),np.dot(SL,VLR)), True, False) , SR)
        Gjir = np.dot(np.dot(SR,VRL),Giir)

        Giitilde = (1.0/(2.0*i))*(Giir - np.conj(Giir))
        Gjjtilde = (1.0/(2.0*i))*(Gjjr - np.conj(Gjjr))
        Gjitilde = (1.0/(2.0*i))*(Gjir - np.conj(Gjir))

        return 4.0*np.trace(np.dot(np.dot(Giitilde, VLR), np.dot(Gjjtilde,VRL)) - np.dot(np.dot(VLR,Gjitilde),np.dot(VLR,Gjitilde)) ).real

size = 8



#for size in xrange(8,9,4):
for strain in np.arange(0.00,0.001,0.01):
#for strain in [0.0]:
    t1, t2 = setValues(strain,np.pi/2.0)
    I = np.identity(2*size+1)
    H00 = H_build(t1,t2)
    VRL = VRL_build(t1,t2)
    VLR = np.transpose(VRL)
#    os.makedirs( str(strain))
#    SLdic = shelve.open(str(strain)+'/SL.db', writeback = True)
#    SRdic = shelve.open(str(strain)+'/SR.db', writeback = True)
#    fldos = open('cond_and_ldos_'+str(strain)+'.dat','w')    
    SLdic = {}
    SRdic = {}
    for Energy in np.arange(-3.0,3.01,0.01):
    	SLoriginal = right_lead(Energy, H00, VLR, VRL)
    	SRoriginal = left_lead(Energy, H00, VLR, VRL)
    
    	SLdic[str(Energy)] = SLoriginal
    	SRdic[str(Energy)] = SRoriginal
    
    	G = np.dot(invert(I - np.dot(np.dot(SLoriginal,VLR),np.dot(SRoriginal, VRL))),SLoriginal)
    	ldos = (-2.0/np.pi)*G[size/2,size/2].imag
        cond = conductance_pristine(Energy, H00, VLR, VRL, SLoriginal, SRoriginal)
#    	print Energy, ldos , cond, strain
        print Energy, G[0][0].real, G[0][0].imag
 #       fldos.write( str(Energy) + '\t' + str(ldos) + '\t' + str(cond) + '\n')
        
