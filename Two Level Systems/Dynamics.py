from Basis_states import *
from J_op import *
from Initialization import *

def evolution(Numberbasis,gamma, k, gamma_phi, gamma_col, k_col, gamma_phi_col):
    J_right11=collective_J_right(1,1,Numberbasis)
    J_right00=collective_J_right(0,0,Numberbasis)
    J_right10=collective_J_right(1,0,Numberbasis)
    J_right01=collective_J_right(0,1,Numberbasis)  
    J_left11=collective_J_left(1,1,Numberbasis)
    J_left00=collective_J_left(0,0,Numberbasis)
    J_left10=collective_J_left(1,0,Numberbasis)
    J_left01=collective_J_left(0,1,Numberbasis)

    collective_Jz_right=0.5*(J_right11 - J_right00)
    collective_Jz_left=0.5*(J_left11 - J_left00)
    collective_Jy_right=0.5j*(J_right10 - J_right01)
    collective_Jy_left=0.5j*(J_left10 - J_left01)


    #Defining the rate of change of P values
    dP_dt=lil_matrix((Numberbasis.dim(),Numberbasis.dim()))
    for idx, n in enumerate(Numberbasis.basis):
        nij=np.array(n).reshape(2,2)


        #For decay
        if gamma!=0:
            nij_copy=nij.copy()
            nij_copy[1,1]+=1
            nij_copy[0,0]-=1

            i=Numberbasis.get_index(nij_copy.flatten())
   
            dP_dt[idx,idx]+=-gamma/2*((2*nij[1,1]+nij[1,0]+nij[0,1]))

            if i!=-1:
                dP_dt[idx,i]+=gamma/2*(2*(nij[1,1]+1))


        #For pumping
        if k!=0:
            nij_copy=nij.copy()
            nij_copy[1,1]-=1
            nij_copy[0,0]+=1

            j=Numberbasis.get_index(nij_copy.flatten())

            dP_dt[idx,idx]+=-k/2*((2*nij[0,0]+nij[1,0]+nij[0,1]))

            if j!=-1:
                dP_dt[idx,j]+=k/2*(2*(nij[0,0]+1))


        #For dephasing
        if gamma_phi!=0:
            dP_dt[idx,idx]+=-gamma_phi*(nij[1,0]+nij[0,1])


    
    #For collective Decay
    if gamma_col!=0:
        dP_dt+=gamma_col/2*(2*J_left01@J_right10-J_left10@J_left01-J_right01@J_right10)

    #For collective Pumping
    if k_col!=0:
        dP_dt+=k_col/2*(2*J_left10@J_right01-J_left01@J_left10-J_right10@J_right01)

    return dP_dt


def result(t_final,dP_dt,Numberbasis,N,P_values):
    t_steps=100

    P_values0=P_values.copy()
    Pt=expm_multiply(dP_dt, P_values0, 0, t_final, t_steps).transpose()

    #Defining the observables
    #Extracting <J11>
    index_list_1=[]
    for idx,n in enumerate(Numberbasis.basis):
        nij=np.array(n).reshape(2,2)

        if nij[0,0]+nij[1,1]==N:
            index_list_1.append(idx)


    J11=0*Pt[0]

    for i in range(1,len(index_list_1)):
        J11+=i*Pt[index_list_1[i]]    #Population of the |1><1|


    J00=np.zeros(len(J11))

    for i in range(len(J00)):
        J00[i]=N-J11[i]
    results=[J00,J11]

    return results


def output(results,N,final_time):
    t_step=100
    t_eval=np.linspace(0,final_time,t_step)

    lab=['Ground state population','1st excited state population','2nd excited state population']

    plt.figure(figsize=(10,6))
    for i in range(len(results)):
        plt.plot(t_eval,results[i],label=lab[i])

    plt.xlabel("Time")
    plt.ylabel("Population")
    plt.title(r"Evolution of population of various states for N="+str(N)+" qubits")
    plt.grid(True)
    plt.legend()
    plt.show()
