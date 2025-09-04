import qutip.piqs as piqs
from qutip import mesolve, expect
import time
import numpy as np

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

t_eval=np.linspace(0,5,100)

#PIQS Simulation
start_time_piqs=time.time()

[jx , jy , jz] = piqs.jspin (N)
piqs_sys = piqs.Dicke(N=N, pumping=k,emission=gamma, dephasing=gamma_phi, collective_emission=gamma_col, collective_pumping=k_col)
D=piqs_sys.liouvillian()
#rho0_piqs = piqs.ghz(N)
#rho0_piqs = piqs.dicke (N,N/2,-N/2)
rho0_piqs = piqs.css (N,np.pi-theta,-phi,basis="dicke",coordinates="polar")

result = mesolve (D , rho0_piqs , t_eval,[])
rhot = result.states
jzt = expect (rhot , jz)
piqs_pop_n11 = [i+N/2 for i in jzt]

jxt = expect (rhot , jx)

end_time_piqs=time.time()
print("Time taken for PIQS to run: ",end_time_piqs-start_time_piqs)
