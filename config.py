from math import exp as expRe

# Functions for determining paramters
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

# math parameters
hbar = 1.0
eta = 1.0e-4
# material parameters
t = -1.0
EF = 0.0

# Impurity parameters
eps_imp = 1.0
tau = -1.0

# Dynamic parameters
U = 10.0
hw0 = 1.0e-3	# A default value for the Zeeman field. 
dtol = 1.0e-6		# Stands for "default tolerance". Could be better
wf=EF	# This needs a lot of thought/work

# Strain parameters
Seps = 0.1	# The epsilon for strain
Ssigma = 0.165		# Poisson's ration in graphene
Salpha = 3.37		# A strain constant, taken from the literature

t1,t2 = SHoppingA(Seps,Ssigma)
tau1,tau2 = SHoppingA(Seps,Ssigma)

