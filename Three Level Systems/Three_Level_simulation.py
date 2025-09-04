#Imports
import numpy as np
from itertools import product
import math
from scipy.sparse import kron, csr_matrix
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
import qutip.piqs as piqs
from qutip import mesolve, expect
from scipy.sparse.linalg import expm_multiply
from scipy.sparse import csc_matrix,csr_matrix,lil_matrix


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


#Creating the number basis class
class SymmetricLiouvilleQubitNumberBasis:
    def __init__(self, N):
        self.N = N
        self.q = 3
        self.d = self.q ** 2
        self.basis = self._generate_basis()
        self.index_map = {tuple(state): i for i, state in enumerate(self.basis)}

    def _generate_basis(self):
        def compositions(n, k):
            if k == 1:
                yield (n,)
            else:
                for i in range(n + 1):
                    for rest in compositions(n - i, k - 1):
                        yield (i,) + rest

        basis = list(compositions(self.N, self.d))
        basis.reverse()  
        return basis

    def dim(self):
        return len(self.basis)

    def get_index(self, n_tuple):
        return self.index_map.get(tuple(n_tuple), -1)

    def show_basis(self, compact=True):
        for i, state in enumerate(self.basis):
            if compact:
                print(f"{i}: {state}")
            else:
                mat = np.array(state).reshape(self.q, self.q)
                print(f"{i}:\n{mat}\n")

Numberbasis=SymmetricLiouvilleQubitNumberBasis(N)


#Defining the collective operators
def collective_J_right(x,y):
    L= lil_matrix((Numberbasis.dim(),Numberbasis.dim()))

    for idx,n in enumerate(Numberbasis.basis):  
        nij=np.array(n).reshape(3,3)

        for k in range(3):
            n_copy = nij.copy()
            n_copy[x,k] += 1
            n_copy[y,k] -= 1
           
            new_index = Numberbasis.get_index(n_copy.flatten())
           
            if new_index != -1:
                L[idx,new_index]+=(nij[x,k]+1) if nij[y,k] > 0 else 0
    return L

def collective_J_left(x,y):
    M= lil_matrix((Numberbasis.dim(),Numberbasis.dim()))
    for idx,n in enumerate(Numberbasis.basis):  
        nij=np.array(n).reshape(3,3)

        for k in range(3):
            n_copy = nij.copy()
            n_copy[k,x] -= 1
            n_copy[k,y] += 1
           
            new_index = Numberbasis.get_index(n_copy.flatten())
           
            if new_index != -1:
                M[idx,new_index]+=(nij[k,y]+1) if nij[k,x] > 0 else 0

    return M

J_right22=collective_J_right(2,2)
J_right11=collective_J_right(1,1)
J_right00=collective_J_right(0,0)
J_right20=collective_J_right(2,0)
J_right02=collective_J_right(0,2)
J_right21=collective_J_right(2,1)
J_right12=collective_J_right(1,2)
J_right10=collective_J_right(1,0)
J_right01=collective_J_right(0,1)  

J_left22=collective_J_left(2,2)
J_left11=collective_J_left(1,1)
J_left00=collective_J_left(0,0)
J_left20=collective_J_left(2,0)
J_left02=collective_J_left(0,2)
J_left21=collective_J_left(2,1)
J_left12=collective_J_left(1,2)
J_left10=collective_J_left(1,0)
J_left01=collective_J_left(0,1)


P_values0=np.zeros(Numberbasis.dim())
#P_values0[0]=1  #Ground state
P_values0[-1]=1 #Excited state

P_values=P_values0.copy()

#Constructing the Entire Liouvillian

Hamiltonian=lil_matrix((Numberbasis.dim(),Numberbasis.dim()))              #Hamiltonian


#Defining the rate of change of P values
dP_dt=lil_matrix((Numberbasis.dim(),Numberbasis.dim()))
for idx, n in enumerate(Numberbasis.basis):
    nij=np.array(n).reshape(3,3)


    #For decay 2->0
    if gamma20!=0:
        nij_copy=nij.copy()
        nij_copy[2,2]+=1
        nij_copy[0,0]-=1

        i=Numberbasis.get_index(nij_copy.flatten())
   
        dP_dt[idx,idx]+=-gamma20/2*((2*nij[2,2]+nij[2,0]+nij[0,2]+nij[2,1]+nij[1,2]))

        if i!=-1:
            dP_dt[idx,i]+=gamma20/2*(2*(nij[2,2]+1))

    #For pumping 0->2
    if k02!=0:
        nij_copy=nij.copy()
        nij_copy[2,2]-=1
        nij_copy[0,0]+=1

        j=Numberbasis.get_index(nij_copy.flatten())

        dP_dt[idx,idx]+=-k02/2*((2*nij[0,0]+nij[1,0]+nij[0,1]+nij[2,0]+nij[0,2]))

        if j!=-1:
            dP_dt[idx,j]+=k02/2*(2*(nij[0,0]+1))

    #For dephasing 0,2
    if gamma_phi02!=0:
        dP_dt[idx,idx]+=-gamma_phi02*(nij[2,0]+nij[0,2])



    #For decay 1->0
    if gamma10!=0:
        nij_copy=nij.copy()
        nij_copy[1,1]+=1
        nij_copy[0,0]-=1

        i=Numberbasis.get_index(nij_copy.flatten())
   
        dP_dt[idx,idx]+=-gamma10/2*((2*nij[1,1]+nij[1,0]+nij[0,1]+nij[1,2]+nij[2,1]))

        if i!=-1:
            dP_dt[idx,i]+=gamma10/2*(2*(nij[1,1]+1))

    #For pumping 0->1
    if k01!=0:
        nij_copy=nij.copy()
        nij_copy[1,1]-=1
        nij_copy[0,0]+=1

        j=Numberbasis.get_index(nij_copy.flatten())

        dP_dt[idx,idx]+=-k01/2*((2*nij[0,0]+nij[1,0]+nij[0,1]+nij[0,2]+nij[2,0]))

        if j!=-1:
            dP_dt[idx,j]+=k01/2*(2*(nij[0,0]+1))

    #For dephasing 0,1
    if gamma_phi01!=0:
        dP_dt[idx,idx]+=-gamma_phi01*(nij[1,0]+nij[0,1])



    #For decay 2->1
    if gamma21!=0:
        nij_copy=nij.copy()
        nij_copy[2,2]+=1
        nij_copy[1,1]-=1

        i=Numberbasis.get_index(nij_copy.flatten())
   
        dP_dt[idx,idx]+=-gamma21/2*((2*nij[2,2]+nij[2,1]+nij[1,2]+nij[0,2]+nij[2,0]))

        if i!=-1:
            dP_dt[idx,i]+=gamma21/2*(2*(nij[2,2]+1))

    #For pumping 1->2
    if k12!=0:
        nij_copy=nij.copy()
        nij_copy[2,2]-=1
        nij_copy[1,1]+=1

        j=Numberbasis.get_index(nij_copy.flatten())

        dP_dt[idx,idx]+=-k12/2*((2*nij[1,1]+nij[2,1]+nij[1,2]+nij[1,0]+nij[0,1]))

        if j!=-1:
            dP_dt[idx,j]+=k12/2*(2*(nij[1,1]+1))

    #For dephasing 1,2
    if gamma_phi12!=0:
        dP_dt[idx,idx]+=-gamma_phi12*(nij[2,1]+nij[1,2])


    
#For collective Decay 2->0
if gamma_col20!=0:
    dP_dt+=gamma_col20/2*(2*J_left02@J_right20-J_left20@J_left02-J_right02@J_right20)

#For collective Pumping  0->2
if k_col02!=0:
    dP_dt+=k_col02/2*(2*J_left20@J_right02-J_left02@J_left20-J_right20@J_right02)


#For collective Decay 2->1
if gamma_col21!=0:
    dP_dt+=gamma_col21/2*(2*J_left12@J_right21-J_left21@J_left12-J_right12@J_right21)

#For collective Pumping  1->2
if k_col12!=0:
    dP_dt+=k_col12/2*(2*J_left21@J_right12-J_left12@J_left21-J_right21@J_right12)


#For collective Decay 1->0
if gamma_col10!=0:
    dP_dt+=gamma_col10/2*(2*J_left01@J_right10-J_left10@J_left01-J_right01@J_right10)

#For collective Pumping  0->1
if k_col01!=0:
    dP_dt+=k_col01/2*(2*J_left10@J_right01-J_left01@J_left10-J_right10@J_right01)


#Time parameters
t_final=50
t_steps=100
t_eval = np.linspace(0, t_final, t_steps)

P_values0=P_values.copy()
Pt=expm_multiply(dP_dt, P_values0, 0, t_final, t_steps).transpose()

#Defining the observables
#Extracting <J11>
index_list_1=[]
for idx,n in enumerate(Numberbasis.basis):
    nij=np.array(n).reshape(3,3)
   

    if nij[0,0]+nij[1,1]==N:
        index_list_1.append(idx)


J11=0*Pt[0]

for i in range(1,len(index_list_1)):
    J11+=i*Pt[index_list_1[i]]    #Population of the |1><1|


#Extracting <J22>
index_list_2=[]
for idx,n in enumerate(Numberbasis.basis):
    nij=np.array(n).reshape(3,3)
   

    if nij[0,0]+nij[2,2]==N:
        index_list_2.append(idx)


J22=0*Pt[0]

for i in range(1,len(index_list_2)):
    J22+=i*Pt[index_list_2[i]]    #Population of the |2><2|

J00=np.zeros(len(J22))

for i in range(len(J00)):
    J00[i]=N-J11[i]-J22[i]


plt.figure(figsize=(10,6))

plt.plot(t_eval, J11, label="Population of 1st excited state")
plt.plot(t_eval, J22, label="Population of 2nd excited state")
plt.plot(t_eval, J00, label="Population of ground state")


plt.xlabel("Time")
plt.ylabel("Population")
plt.title(r"Evolution of population of various states for N="+str(N)+" qubits")
plt.grid(True)
plt.legend()
plt.show()