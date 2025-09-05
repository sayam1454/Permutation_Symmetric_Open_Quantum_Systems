from J_operators import *
from Liouville_basis import *

def Initiate_rho(Numberbasis):
    P_values0=np.zeros(Numberbasis.dim())
    #P_values0[0]=1  #Ground state
    P_values0[-1]=1 #Excited state

    P_values=P_values0.copy()

    return P_values