import qutip
import numpy as np
from scipy.linalg import expm
from functools import reduce
import matplotlib.pyplot as plt

def hamiltonian(V, C, E, coefficient=5):
    H = 0
    for i in range(V):
        create_power = qutip.create(3)**len(E[i])
        h = [qutip.identity(3)] * (V+C)
        h[C+i] = create_power
        for j in E[i]:
            h[j] = qutip.destroy(3)
        if H == 0:
            H = qutip.tensor(h)
        else:
            H += qutip.tensor(h)
    H = coefficient * H
    return H + H.dag()
    

def drive_graph(V, C, E, vecs, time, sample_t, coefficient=5, psi0=[], eigenstates=False):
    H = hamiltonian(V, C, E, coefficient=coefficient)
    print("generated H")

    step = 0.01

    t = np.arange(0, time, step)
    results = qutip.krylovsolve(H, psi0, t, krylov_dim=20, sparse=True)

    print("exponentiated")
    fidelities = [[] for _ in range(len(vecs))]

    for i, v in enumerate(vecs):
        for s in results.states:
            fidelities[i].append(np.abs(s.overlap(v))**2)
        
    print("computed")
    for i in range(len(fidelities)):
        plt.plot(t, fidelities[i], label=i)
    plt.legend()
    plt.show()

N=15
levels=3
psi0 = qutip.ket([1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 3)
sol1 = qutip.ket([0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 2], 3)
sol2 = qutip.ket([0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 2, 0], 3)
sol3 = qutip.ket([0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 2], 3)
sol4 = qutip.ket([0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 2, 0, 0], 3)
sol5 = qutip.ket([0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 2, 0], 3)
sol6 = qutip.ket([0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 2, 0, 0], 3)

target = (sol1 + sol2 + sol3 + sol4 + sol5 + sol6).unit()
drive_graph(9, 6, [[0, 3], [0, 4], [0, 5], [1, 3], [1, 4], [1, 5], [2, 3], [2, 4], [2, 5]], [psi0, target], time=5, sample_t=3, psi0=psi0)