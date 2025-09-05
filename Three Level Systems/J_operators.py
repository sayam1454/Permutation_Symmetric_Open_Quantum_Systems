from Liouville_basis import *

#Defining the collective operators
def collective_J_right(x,y,Numberbasis):
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

def collective_J_left(x,y,Numberbasis):
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
