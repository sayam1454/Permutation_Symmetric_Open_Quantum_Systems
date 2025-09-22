from Liouville_basis import *
from State_initializtion import *
from J_operators import *
from Evolution import *
from qutip import *


# Number of qutrits
N=2

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
t_eval = np.linspace(0, t_final, 100)

# define basis states

f = basis(3,0)  #excited state
e = basis(3,1)  #intermediate state
g = basis(3,2)  #ground state

# initial state
#psi0 = tensor([g]*N)     #ground state
psi0 = tensor([f]*N)     #Fully excited state
#psi0 = tensor([e]*N) 
rho0 = ket2dm(psi0)

# Hamiltonian (zero for now)
H = 0 * tensor([qeye(3)]*N)

# Takes from i to j state
def sigma_operator(i,j):
    mat = np.zeros((3,3))
    mat[abs(j-2),abs(i-2)]=1
    return Qobj(mat)

# Local operator acting on site k
def sigma_operator_local(i,j,N,k):
    ident = qeye(3)
    op_list = [ident]*N
    op_list[k] = sigma_operator(i,j)
    return tensor(op_list)

# Collective operator (sum over all sites)
def sigma_operator_col(i,j,N):
    return sum(sigma_operator_local(i,j,N,k) for k in range(N))

c_ops = []

# local channels
for k in range(N):
    c_ops += [
        np.sqrt(gamma20) * sigma_operator_local(2,0,N,k),
        np.sqrt(k02) * sigma_operator_local(0,2,N,k),
        np.sqrt(gamma_phi02) * (sigma_operator_local(2,2,N,k) - sigma_operator_local(0,0,N,k)),

        np.sqrt(gamma21) * sigma_operator_local(2,1,N,k),
        np.sqrt(k12) * sigma_operator_local(1,2,N,k),
        np.sqrt(gamma_phi12) * (sigma_operator_local(2,2,N,k) - sigma_operator_local(1,1,N,k)),

        np.sqrt(gamma10) * sigma_operator_local(1,0,N,k),
        np.sqrt(k01) * sigma_operator_local(0,1,N,k),
        np.sqrt(gamma_phi01) * (sigma_operator_local(1,1,N,k) - sigma_operator_local(0,0,N,k)),
    ]

# collective channels
c_ops += [
    np.sqrt(gamma_col20) * sigma_operator_col(2,0,N),
    np.sqrt(k_col02) * sigma_operator_col(0,2,N),
    np.sqrt(gamma_phi_col02) * (sigma_operator_col(2,2,N) - sigma_operator_col(0,0,N)),

    np.sqrt(gamma_col21) * sigma_operator_col(2,1,N),
    np.sqrt(k_col12) * sigma_operator_col(1,2,N),
    np.sqrt(gamma_phi_col12) * (sigma_operator_col(2,2,N) - sigma_operator_col(1,1,N)),

    np.sqrt(gamma_col10) * sigma_operator_col(1,0,N),
    np.sqrt(k_col01) * sigma_operator_col(0,1,N),
    np.sqrt(gamma_phi_col01) * (sigma_operator_col(1,1,N) - sigma_operator_col(0,0,N)),
]

# observables: collective populations
proj_g = sigma_operator_col(0,0,N)
proj_e = sigma_operator_col(1,1,N)
proj_f = sigma_operator_col(2,2,N)

e_ops = [proj_g, proj_e, proj_f]

final_val = mesolve(H, rho0, t_eval, c_ops, e_ops)


# Calculating required parameters
Numberbasis=SymmetricLiouvilleQutritNumberBasis(N)
P_values=Initiate_rho(Numberbasis)
Liouvillian=Evolution(Numberbasis,gamma20, k02, gamma_phi02, gamma_col20, k_col02, gamma_phi_col02,gamma21, k12, gamma_phi12, gamma_col21, k_col12, gamma_phi_col12,gamma10, k01, gamma_phi01, gamma_col10, k_col01, gamma_phi_col01)
results=final_result(t_final,Liouvillian,Numberbasis,N,P_values)




#Plotting
lab=['Ground state population','1st excited state population','2nd excited state population']

plt.figure(figsize=(10,6))
for i in range(len(results)):
    plt.plot(t_eval,results[i],label=lab[i])
plt.plot(t_eval, final_val.expect[0],"*", label="Population of ground state from qutip")
plt.plot(t_eval, final_val.expect[1],'*', label="Population of 1st excited state from qutip")
plt.plot(t_eval, final_val.expect[2],'*', label="Population of 2nd excited state from qutip")

plt.xlabel("Time")
plt.ylabel("Population")
plt.title(r"Evolution of population of various states for N="+str(N)+" qutrits")
plt.grid(True)
plt.legend()
plt.show()
