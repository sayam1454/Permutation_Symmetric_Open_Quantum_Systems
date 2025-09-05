from Liouville_basis import *
from State_initializtion import *
from J_operators import *
from Evolution import *


# Number of qutrits
N=10

# Lindblad parameters
gamma20=0             #Decay rate between 2->0
k02=0                 #Pumping rate between 0->2
gamma_phi02=0         #Dephasing rate between 0 and 2 states
gamma_col20=0        #Collective Decay rate between 2->0
k_col02=0.5             #Collective Pumping rate between 0->2
gamma_phi_col02=0     #Collective Dephasing rate between 0 and 2 states

gamma21=0.3             #Decay rate between 2->1
k12=0                #Pumping rate between 1->2
gamma_phi12=0         #Dephasing rate between 1 and 2 states
gamma_col21=0         #Collective Decay rate between 2->1
k_col12=0             #Collective Pumping rate between 1->2
gamma_phi_col12=0     #Collective Dephasing rate between 1 and 2 states


gamma10=0             #Decay rate between 1->0 
k01=0                 #Pumping rate between 0->1
gamma_phi01=0         #Dephasing rate between 0 and 1 states
gamma_col10=0.8         #Collective Decay rate between 1->0 
k_col01=0             #Collective Pumping rate between 0->1
gamma_phi_col01=0     #Collective Dephasing rate between 0 and 1 states


t_final=8


# Calculating required parameters
Numberbasis=SymmetricLiouvilleQubitNumberBasis(N)
P_values=Initiate_rho(Numberbasis)
Liouvillian=Evolution(Numberbasis,gamma20, k02, gamma_phi02, gamma_col20, k_col02, gamma_phi_col02,gamma21, k12, gamma_phi12, gamma_col21, k_col12, gamma_phi_col12,gamma10, k01, gamma_phi01, gamma_col10, k_col01, gamma_phi_col01)
results=result(t_final,Liouvillian,Numberbasis,N,P_values)

output(results,N,t_final)