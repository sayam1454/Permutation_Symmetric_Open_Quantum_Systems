from Basis_states import *
from J_op import *

def Initiate_rho(Numberbasis,theta,phi):
    P_values0=np.zeros(Numberbasis.dim())
    P_values0[0]=1  #Ground state
    
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

    P_values=P_values0.copy()

     #Taking to arbitary theta

    if theta != 0:
        H=-(collective_Jy_right-collective_Jy_left)*1j
        P_values=expm_multiply(H, P_values0, 0, theta, 2)[-1]


     #Taking to arbitary phi
    if phi != 0:
        H=(collective_Jz_right-collective_Jz_left)*1j
        P_values0=P_values.copy()
        P_values=expm_multiply(H, P_values0, 0, phi, 2)[-1]

    return P_values