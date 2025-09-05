from J_operators import *
from Liouville_basis import *
from State_initializtion import *

def Evolution(Numberbasis,gamma20, k02, gamma_phi02, gamma_col20, k_col02, gamma_phi_col02,gamma21, k12, gamma_phi12, gamma_col21, k_col12, gamma_phi_col12,gamma10, k01, gamma_phi01, gamma_col10, k_col01, gamma_phi_col01):
    J_right22=collective_J_right(2,2,Numberbasis)
    J_right11=collective_J_right(1,1,Numberbasis)
    J_right00=collective_J_right(0,0,Numberbasis)
    J_right20=collective_J_right(2,0,Numberbasis)
    J_right02=collective_J_right(0,2,Numberbasis)
    J_right21=collective_J_right(2,1,Numberbasis)
    J_right12=collective_J_right(1,2,Numberbasis)
    J_right10=collective_J_right(1,0,Numberbasis)
    J_right01=collective_J_right(0,1,Numberbasis)  
    J_left22=collective_J_left(2,2,Numberbasis)
    J_left11=collective_J_left(1,1,Numberbasis)
    J_left00=collective_J_left(0,0,Numberbasis)
    J_left20=collective_J_left(2,0,Numberbasis)
    J_left02=collective_J_left(0,2,Numberbasis)
    J_left21=collective_J_left(2,1,Numberbasis)
    J_left12=collective_J_left(1,2,Numberbasis)
    J_left10=collective_J_left(1,0,Numberbasis)
    J_left01=collective_J_left(0,1,Numberbasis)

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

    return dP_dt


def result(t_final,dP_dt,Numberbasis,N,P_values):
    t_steps=100

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

    results=[J00,J11,J22]

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
    plt.title(r"Evolution of population of various states for N="+str(N)+" qutrits")
    plt.grid(True)
    plt.legend()
    plt.show()