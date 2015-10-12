from GFRoutines import *
from numpy.linalg import norm
from Recursive import *


def gGNRSubs1(nE,m,n,E):
  g = gRib_Arm(nE,m,n,m,n,0,E)
  return Dyson1(g,eps_imp)


def gGNRTop1(nE,m,n,E):      
  """Calculates the onsite GF of a top adsorbed impurity in a GNR"""
  # Introduce the connecting GFs
  g = np.zeros((2,2),dtype=complex)
  g[0,0] = gRib_Arm(nE,m,n,m,n,0,E)	# Connecting site
  g[1,1] = 1.0/(E-eps_imp)		# Impurity site
  
  V = np.zeros([2,2],dtype=complex)	# Connecting perturbation
  V[0,1] = V[1,0] = tau
  
  G = Dyson(g,V)
  
  # Return the part of the matrix that governs the impurity behaviour. 
  return G[1,1]


def GMx2Subs(nE,mI,nI,mP,nP,s,E):
  """Creates the GF matrix for a single subsitutional impurity. Order is Probe -> subsitutional"""
  g = gMx2GNR(nE,mP,nP,mI,nI,s,E)
  V = [[0,0],[0,eps_imp]]
  return Dyson(g,V)


def GMx1CenterProbe(nE,mC,nC,mP,nP,sP,E):  
  """Returns the GF Mx of the center adosrbed impurity and also a single probe site.
  The site order is probe,hexagon,impurity site."""
  n = 8		# Total number of sites (hexagon + impurity + probe)
  
  rC = np.array([mC,nC,0])	# Position of center adsorbed impurity (bottom left site)
  rHex = np.array([[0,0,0],[0,0,1],[1,0,0],[1,-1,1],[1,-1,0],[0,-1,1]])		# All of the sites of a hexagon (w.r.t bottom left)
  r = np.concatenate(([[mP,nP,sP]],rHex + rC))	# Our complete list of positions, probe site, hexagon
  
  gMx = np.zeros((n,n),dtype=complex)
  gMx[:n-1,:n-1] = gMxnGNR(nE,r,E)
  gMx[n-1,n-1] = 1.0/(E-eps_imp)
    
  V = np.zeros([n,n],dtype=complex)
  V[:n-1,n-1] = tau
  V[n-1,:n-1] = tau
  
  GMx = Dyson(gMx,V)
  
  return GMx


def GMx1Center(nE,mC,nC,E):  
  """Gets the GF a Center Adsorbed impurity."""
  n = 7		# Total number of sites (hexagon + impurity)
  
  rC = np.array([mC,nC,0])	# Position of center adsorbed impurity (bottom left site)
  rHex = np.array([[0,0,0],[0,0,1],[1,0,0],[1,-1,1],[1,-1,0],[0,-1,1]])		# All of the sites of a hexagon (w.r.t bottom left)
  r = rHex + rC		# Our complete list of positions, probe site, hexagon
  
  gMx = np.zeros((n,n),dtype=complex)
  gMx[:n-1,:n-1] = gMxnGNR(nE,r,E)	# First n-1 elements are the hexagon sites
  gMx[n-1,n-1] = 1.0/(E-eps_imp)	# Final site is that of the impurity
    
  V = np.zeros([n,n],dtype=complex)
  V[0,n-1] = V[1,n-1] = V[3,n-1] = V[4,n-1] = tau1
  V[n-1,0] = V[n-1,1] = V[n-1,3] = V[n-1,4] = tau1
  
  V[2,n-1] = V[5,n-1] = tau2
  V[n-1,2] = V[n-1,5] = tau2
  
  GMx = Dyson(gMx,V)
  
  return GMx



def GMxCenterRec2(N,p,ImpList,E):
  """Calculates the GF of a strip in an AGNR in the presence of center adsorbed impurities"""
  nimp = len(ImpList)
  gL,gR,VLR,VRL = Leads(N,E)		# Get Leads
  H = HArmStrip(N,p,CenterList=ImpList)	# Hamiltonian with Center adsorbed impurities
  gC = gGen(E,H)
  sizeP = gL.shape[0]
  sizeI = gC.shape[0]
  VLRb, VRLb = PadZeros(VLR,(sizeP,sizeI)), PadZeros(VRL,(sizeI,sizeP))	# To get the full GF, VLR and VRL must be padded to match left and right cells.
  gL = RecAdd(gL,gC,VLRb,VRLb)
  VLRb, VRLb = PadZeros(VLR,(sizeI,sizeP)), PadZeros(VRL,(sizeP,sizeI))
  g = RecAdd(gR,gL,VRLb,VLRb)
  return g


def GMxSubsRec(N,ImpList,E):
  """Gets the GF of a strip containing a single substitutional impurity"""
  gL,gR,VLR,VRL = Leads(N,E)
  HM = HArmStripSubs(N,ImpList)
  gM = gGen(E,HM)
  gM = RecAdd(gL,gM,VLR,VRL)
  gM = RecAdd(gR,gM,VRL,VLR)
  return gM


def GMxTopRec(N,ImpList,E):
  """Gets the GF of a strip containing a single top-adsorbed impurity"""
  gL,gR,VLR,VRL = Leads(N,E)
  gL = PadZeros(gL,(11,11))
  gR = PadZeros(gR,(11,11))
  VLR = PadZeros(VLR,(11,11))
  VRL = PadZeros(VRL,(11,11))
  HM = HArmStripTop(N,ImpList)
  gM = gGen(E,HM)
  gM = RecAdd(gL,gM,VLR,VRL)
  gM = RecAdd(gR,gM,VRL,VLR)
  return gM


def GMxCenterRec(N,ImpList,E):
  """Gets the GF of a strip containing a single center adsorbed impurity, cannot deal with interstitial sites"""
  gL,gR,VLR,VRL = Leads(N,E)
  gL = PadZeros(gL,(11,11))
  gR = PadZeros(gR,(11,11))
  VLR = PadZeros(VLR,(11,11))
  VRL = PadZeros(VRL,(11,11))
  HM = HArmStripCenter(N,ImpList)
  gM = gGen(E,HM)
  gM = RecAdd(gL,gM,VLR,VRL)
  gM = RecAdd(gR,gM,VRL,VLR)
  return gM


if __name__ == "__main__":  
  nE = 6
  m,n = 1,0
  Elist = np.linspace(-3.0+1j*eta,3.0+1j*eta,201)
  DOSSlist = [-gGNRSubs1(nE,m,n,E).imag/pi for E in Elist]
  DOSTlist = [-gGNRTop1(nE,m,n,E).imag/pi for E in Elist]
  DOSClist = [-GMx1Center(nE,1,0,E)[6,6].imag/pi for E in Elist]
  pl.plot(Elist.real,DOSSlist,label='Subs')
  pl.plot(Elist.real,DOSTlist,label='Top')
  pl.plot(Elist.real,DOSClist,label='Center')
  pl.legend()
  pl.show()

  

