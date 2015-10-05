from GF import *
from math import exp as expRe
from scipy.integrate import dblquad

Seps = 0.2	# The epsilon for strain
Ssigma = 0.165		# Poisson's ration in graphene
Salpha = 3.37		# A strain constant, taken from the literature


def SHoppingGen(ratio1,ratio2):
  """Calculates the hopping integrals t1, t2 given the ratio of the relevant bond lengths to the pristine length."""
  t1 = t*expRe(-Salpha*(ratio1-1.0))
  t2 = t*expRe(-Salpha*(ratio2-1.0))
  return t1, t2


def SHoppingZ(Seps,Ssigma):
  """Gets the correct bond length ratios for zigzag strain"""
  ratio1 = 1.0 + 3.0*Seps/4.0 - Seps*Ssigma/4.0	# R1/R0
  ratio2 = 1.0 - Seps*Ssigma	# R2/R0
  return SHoppingGen(ratio1,ratio2)


def SHoppingA(Seps,Ssigma):
  """Gets the correct bond length ratios for armchair strain"""
  ratio1 = 1.0 + Seps/4.0 - 3.0/4.0*Seps*Ssigma	# R1/R0
  ratio2 = 1.0 + Seps	# R2/R0
  return SHoppingGen(ratio1,ratio2)


def gSBulkDouble(m,n,s,E):
  """Calculates the graphene GF under strain. Uses the double integral. Tested only for pristine values"""
  def gI(kA,kZ):  
    if s == 0:
      N = E
    elif s == 1:
      N = t2 + 2*t1*exp(1j*kA)*cos(kZ)
    else:
      N = t2 + 2*t1*exp(-1j*kA)*cos(kZ)
    Beps2 = t2**2+4*t1*t2*cos(kA)*cos(kZ)+4*t1**2 *cos(kZ)**2
    
    g = (1.0/(2*pi**2))*N*exp(1j*kA*(m+n)+1j*kZ*(m-n))/(E**2-Beps2)
    return g
  
  def gIRe(kA,kZ):
    return gI(kA,kZ).real
  def gIIm(kA,kZ):
    return gI(kA,kZ).imag
  
  gr = dblquad(gIRe, -pi/2.0, pi/2.0, lambda kZ: -pi, lambda kZ: pi)[0]
  gim = dblquad(gIIm, -pi/2.0, pi/2.0, lambda kZ: -pi, lambda kZ: pi)[0]
  return gr + 1j*gim


def gSBulk(m,n,s,E):
  """Calculates the graphene GF with strain. Uses only one integral"""
  def gI(kZ):
    q = acos( (E**2 - t2**2 - 4.0*t1**2 *cos(kZ)**2)/(4.0*t1*t2 *cos(kZ) ) )

    if q.imag < 0.0: q = -q

    Const = 1j/(4*pi*t1*t2)
    Den = cos(kZ)*sin(q)

    if s == 0:
      sig = copysign(1,m+n)
      return Const*E*exp( 1j*(sig*q*(m+n) + kZ*(m-n) ) )/ Den 
    elif s == 1:
      sig = copysign(1,m+n)
      f = t2 + 2*t1*exp(sig*1j*q)*cos(kZ)
      return Const*f*exp( 1j*(sig*q*(m+n) + kZ*(m-n) ) )/Den  
    elif s == -1:
      sig = copysign(1,m+n-1)
      ft = t2 + 2*t1*exp(-sig*1j*q)*cos(kZ)
      return Const*ft*exp( 1j*(sig*q*(m+n) + kZ*(m-n) ) )/Den 
    else:
      print "Sublattice error in gSBulk"
      
  return C_int(gI,-pi/2,pi/2)


def gSGNRIntegral(nE,m1,n1,m2,n2,s,E):
  """Calculates the GF of a GNR for a strained system, without the integral being performed analytically"""
  def gI(kA):  
    if s == 0:
      N = E
    elif s == 1:
      N = t2 + 2*t1*exp(1j*kA)*cos(pi*j/nE)
    else:
      N = t2 + 2*t1*exp(-1j*kA)*cos(pi*j/nE)
    Beps2 = t2**2+4*t1*t2*cos(kA)*cos(pi*j/nE)+4*t1**2 *cos(pi*j/nE)**2
    
    g = N*exp(1j*kA*(m2+n2-m1-n1))*sin(pi*j*(m1-n1)/nE)*sin(pi*j*(m2-n2)/nE)/(E**2-Beps2)
    return g
  
  g = 0
  for j in range(1,nE):
    g += C_int(gI, -pi, pi)
  return g/(pi*nE)


def gSGNR(nE,m1,n1,m2,n2,s,E):
  """Calculates the GF for a nanoribbon, want to incorporate strain at some point"""
  def gterm(ky):
    q = acos( (E**2 - t2**2 - 4.0*t1**2 *cos(ky)**2)/(4.0*t1*t2 *cos(ky) ) )
    if q.imag < 0.0: q = -q

    Const = 1j/(2.0*nE*t1*t2)
    Den = cos(ky)*sin(q)

    if s == 0:
      sig = copysign(1,m2+n2-m1-n1)
      return Const*E*exp( 1j*sig*q*(m2+n2-m1-n1) )*sin(ky*(m2-n2))*sin(ky*(m1-n1))/Den 
    elif s == 1:
      sig = copysign(1,m2+n2-m1-n1)
      f = t2 + 2.0*t1*cos(ky)*exp(sig*1j*q)
      return Const*f*exp( 1j*sig*q*(m2+n2-m1-n1) )*sin(ky*(m2-n2))*sin(ky*(m1-n1))/Den 
    elif s == -1:
      sig = copysign(1,m2+n2-m1-n1-1)
      ft = t2 + 2.0*t1*cos(ky)*exp(-sig*1j*q)
      return Const*ft*exp( 1j*sig*q*(m2+n2-m1-n1) )*sin(ky*(m2-n2))*sin(ky*(m1-n1))/Den 
    else:
      print 'Sublattice error in gSGNR'
  
  def limit_term(ky):
    if s == 0:
      N_ab = E
    elif (s == 1) or (s == -1):
      N_ab = t2
    else:
      print 'Sublattice error in gSGNR'
      
    return 2.0*N_ab*sin(ky*(m2-n2))*sin(ky*(m1-n1))/( nE*( E**2 - t2**2 ) )
  
  g = 0.0
  if nE % 2 == 0:
    for j in range(1,nE):
      if j == nE/2: continue		# Avoid singularities
      g += gterm(pi*j/nE)
    if m2+n2-m1-n1 == 0:
      g += limit_term(pi/2)
  else: 
    for j in range(1,nE):
      g += gterm(pi*j/nE)
      
  return g



if __name__ == "__main__":
  t1,t2 = SHoppingA(Seps,Ssigma)
  N = 5
  gL,gR,VLR,VRL = Leads(N,E)
  print RecAdd(gL,gR,VLR,VRL)

