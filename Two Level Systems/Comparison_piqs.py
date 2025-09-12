import qutip.piqs as piqs
from qutip import mesolve, expect
from Basis_states import *
from J_op import *
from Initialization import *
from Dynamics import *
import time

# Number of qubits
N=40  


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
t_eval=np.linspace(0,t_final,100)





#PIQS Simulation
start_time_piqs=time.time()
[jx , jy , jz] = piqs.jspin (N)
piqs_sys = piqs.Dicke(N=N, pumping=k,emission=gamma, dephasing=gamma_phi, collective_emission=gamma_col, collective_pumping=k_col)
D=piqs_sys.liouvillian()
#rho0_piqs = piqs.ghz(N)
#rho0_piqs = piqs.dicke (N,N/2,-N/2)
rho0_piqs = piqs.css (N,np.pi-theta,-phi,basis="dicke",coordinates="polar")

final_val = mesolve (D , rho0_piqs , t_eval,[])
rhot = final_val.states
jzt = expect (rhot , jz)
piqs_pop_n11 = [i+N/2 for i in jzt]
piqs_pop_n00 = [(N-i) for i in piqs_pop_n11]
end_time_piqs=time.time()
print("Time taken for PIQS to run: ",end_time_piqs-start_time_piqs)

# Calculating required parameters
start_time_code=time.time()
Numberbasis=SymmetricLiouvilleQubitNumberBasis(N)
P_values=Initiate_rho(Numberbasis,theta,phi)
Liouvillian=evolution(Numberbasis,gamma, k, gamma_phi, gamma_col, k_col, gamma_phi_col)
results=result(t_final,Liouvillian,Numberbasis,N,P_values)
t_step=100
lab=['Ground state population','1st excited state population','2nd excited state population']
end_time_code=time.time()
print("Time taken for the code to run: ",end_time_code-start_time_code)


plt.figure(figsize=(10,6))
plt.xlabel("Time")
plt.ylabel("Population")
plt.title(r"Evolution of population of various states for N="+str(N)+" qubits")
plt.grid(True)


plt.plot(t_eval, piqs_pop_n00, '*', label="P0 from piqs")
plt.plot(t_eval, piqs_pop_n11, '*', label="P1 from piqs")
for i in range(len(results)):
    plt.plot(t_eval,results[i],label=lab[i])

plt.legend()
plt.show()