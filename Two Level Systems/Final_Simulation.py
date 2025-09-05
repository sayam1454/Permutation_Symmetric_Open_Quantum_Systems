from Basis_states import *
from J_op import *
from Initialization import *
from Dynamics import *

# Number of qubits
N=10  

# Initial state parameters |theta,phi>
theta = np.pi
phi = np.pi/4

# Lindblad parameters
gamma=2           #Decay rate
k=1             #Pumping rate
gamma_phi=1       #Dephasing rate
gamma_col=3       #Collective Decay rate
k_col=0             #Collective Pumping rate
gamma_phi_col=0       #Collective Dephasing rate


t_final=8


# Calculating required parameters
Numberbasis=SymmetricLiouvilleQubitNumberBasis(N)
P_values=Initiate_rho(Numberbasis,theta,phi)
Liouvillian=evolution(Numberbasis,gamma, k, gamma_phi, gamma_col, k_col, gamma_phi_col)
results=result(t_final,Liouvillian,Numberbasis,N,P_values)

output(results,N,t_final)