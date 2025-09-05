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

#Creating the number basis class
class SymmetricLiouvilleQutritNumberBasis:
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