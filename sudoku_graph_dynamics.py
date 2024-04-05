from qutip import ket, create, destroy, identity, tensor, krylovsolve, sigmax, projection, propagator, state_number_index, state_index_number, Qobj, basis
import numpy as np
from scipy.linalg import expm
from scipy.sparse import coo_matrix
from functools import reduce
import matplotlib.pyplot as plt
import itertools
import networkx as nx
import time

# scarred
def H_scar(V, C, E, coefficient=5):
    H = 0
    for i in range(V):
        print("hi")
        create_power = create(len(E[i])+1)**len(E[i])
        h = [identity(2)] * C + [identity(len(E[i])+1) for i in range(V)]
        h[C+i] = create_power
        for j in E[i]:
            h[j] = destroy(2)
        if H == 0:
            H = tensor(h)
        else:
            H += tensor(h)
    H = coefficient * H
    return H + H.dag()

# thermalized (summing over all terms no products to ensure constraints)
def H_thermalized(V, C, E, coefficient = 5):
    H = 0
    for i in range(V):
        for j in E[i]:
            h = [identity(2)] * C + [identity(len(E[i])+1) for i in range(V)]
            h[C+i] = create(len(E[i])+1)
            h[j] = destroy(2)
            if H == 0:
                H = tensor(h)
            else:
                H += tensor(h)
    H = coefficient * H
    return H + H.dag()
    

def drive_graph(H, vecs, time, psi0=[], name=None, eigenstates=False, plot=True):
    step = 0.01
    t = np.arange(0, time, step)
    results = krylovsolve(H, psi0, t, krylov_dim=20, sparse=True)

    # print("exponentiated")
    fidelities = []
    for s in results.states:
        # print(s)
        fidelities.append(sum([np.abs(s.overlap(v))**2 for v in vecs]))

    # print("computed")
    m = 0
    if name != None:
        if plot == True:
            plt.plot(t, fidelities, label=name)
        m = max(fidelities)
    else:
        if plot == True:
            plt.plot(t, fidelities)

    # plt.legend()
    # plt.show()
    return results.states[-1]

def propagate_graph(H, vecs, time, psi0):
    U = propagator(H, time)
    final = U * psi0

    print(np.abs(final.overlap(vecs[0]))**2)
    print(np.abs(final.overlap(vecs[1]))**2)

adj_mat = []

def sudoku_recursion(state, cliques, N, depth):
    # dims = [2]*(3*N**2) + [4]*(N**3)
    basisv = [state]
    # print(depth)
    for i in set(cliques.flatten()):
        site = 3*N**2 + i
        vert = i%(N**2)
        col = N**2 + (i%N) + N*(i//(N**2))
        row = 2*N**2 + (i//N)

        if state[vert] == 1 and state[col] == 1 and state[row] == 1:
            s = state.copy()
            s[site] += 3
            s[vert] -= 1
            s[col] -= 1
            s[row] -= 1

            sol = [sum([(1+j)*s[(3+j)*N**2 + i] for j in range(N)])//3 for i in range(N**2)]
            if sol not in sols:
                sols.append(sol)
                print(len(sols))
                print(sol)

            cliques_copy = cliques.copy()

            cliques_copy = np.delete(cliques_copy, vert)
            cliques_copy = np.delete(cliques_copy, col)
            cliques_copy = np.delete(cliques_copy, row)

            basisv += sudoku_recursion(s, cliques_copy, N, depth+1)
    # print(len(basis))
    return basisv

def sudoku_basis(N, inp):
    # find solutions for each occupation number
    initial = {}
    for i in inp:
        print(i, (inp[i]-1)*(N**2))
        initial[i + (inp[i]-1)*(N**2)] = inp[i]
    print(initial)
    basisv = []
    state = [1]*3*N**2 + [0]*N**3
    cliques = [[] for _ in range(3*N**2)]
    for i in range(N**3):
        vert = i%(N**2)
        col = N**2 + (i%N) + N*(i//(N**2))
        row = 2*N**2 + (i//N)
        cliques[vert].append(i)
        cliques[col].append(i)
        cliques[row].append(i)

    to_remove = []
    for i in initial:
        vert = i%(N**2)
        col = N**2 + (i%N) + N*(i//(N**2))
        row = 2*N**2 + (i//N)
        to_remove.append(vert)
        to_remove.append(col)
        to_remove.append(row)

        state[i + 3*N**2] = initial[i]
        state[vert] = 0
        state[col] = 0
        state[row] = 0

    cliques = [cliques[i] for i in range(len(cliques)) if i not in to_remove]

    cliques = np.array(cliques)
    basisv = sudoku_recursion(state, cliques, N, 0)
    return basisv
    
sols = []
G = nx.Graph()

def sudoku_direct(N, grid, constraints, depth, parent, sols, G, subgrid=True, localize_term=False):
    print(grid)
    print(parent)
    # print(constraints)
    if not any(0 in s for s in grid):
        return
    
    for i in range(N):
        for j in range(N):
            if grid[i,j] != 0:
                continue
            for k in range(1, N+1):
                if constraints[k-1,0,i] or constraints[k-1,1,j] or (subgrid==True and constraints[k-1,2,(i//2)*2 + (j//2)]):
                    continue
                grid2 = grid.copy()
                grid2[i,j] = k
                index = len(sols)
                for l in range(len(sols)):
                    if np.array_equal(grid2, sols[l]):
                        index = l

                G.add_edge(parent, index)
                
                if index != len(sols):
                    continue

                if np.count_nonzero(grid2) == N*N:
                    G.add_edge(index, index)

                const = constraints.copy()
                const[k-1,0,i] = True
                const[k-1,1,j] = True
                if subgrid == True:
                    const[k-1,2,(i//2)*2 + (j//2)] = True
                sols.append(grid2)
                sudoku_direct(N, grid2, const, depth+1, index, sols, G, subgrid=subgrid)

def get_target_states(basisv, count):
    target = []
    for i in range(len(basisv)):
        # if 0 not in basisv[i]:
        #     print(basisv[i])
        #     target.append(basis(len(basisv), i))

        # uncomment this to include dead-end states in our calculations
        # elif basisv[i].count(0) <=2:
        #     zero_ = basisv[i].index(0)
        #     flag = False
        #     b = np.array(basisv[i])
        #     b = np.reshape(b, (3,3))
        #     for j in range(1,3+1):
        #         if j not in b[zero_//3,:] and j not in b[:,zero_%3]:
        #             flag = True
        #     if flag == False:
        #         print(b[zero_//3,:])
        #         print(b[:,zero_%3])
        #         print(basisv[i])
        #         target.append(basis(len(basisv), i))
        # print(list(basisv[i]))
        if list(basisv[i]).count(0) == count:
            target.append(basis(len(basisv), i))
    return target


def sudoku_hamiltonian(initial={}, psi0vec=None, localized=False, data="sudoku_subspace.txt"):
    # read in basis states, convert to index number
    # compute edges by identifying nearest neighbors in scattering graph through hamming distance of 1
    # graph adjacency
    txt = open(data, "r")
    lines = txt.readlines()
    txt.close()

    lines = [lines[i] for i in range(len(lines)) if i%2 == 1]
    basisv = [[0] * 9] + [[int(j) for j in i[1:-2].replace(" ", "").split(",")] for i in lines]

    basis_start = []
    for i in range(len(basisv)):
        flag = False
        for j in range(9):
            if j in initial and basisv[i][j] != initial[j]:
                flag = True
        if flag == False:
            basis_start.append(basisv[i])

    basisv = basis_start
    psi0 = 0

    if psi0vec == None:
        psi0vec = [0]*9
        for i in range(9):
            if i in initial:
                psi0vec[i] = initial[i]
        psi0vec = [psi0vec]

    for i in range(len(basisv)):
        if basisv[i] in psi0vec:
            if psi0 == 0:
                psi0 = basis(len(basisv), i)
            else:
                psi0 += basis(len(basisv), i)
    psi0 = psi0.unit()
    # print(basisv)
    target = get_target_states(basisv, 0)
    H = []
    for i in range(len(basisv)):
        H.append([])
        H[i] = [0] * len(basisv)
        # tilted potential
        c = basisv[i].count(0)
        if c == 0:
            c = 1
        if localized and basisv[i].count(0) == 0:
            print("hi")
            H[i][i] = 100
        # H[i][i] = c
        for j in range(len(basisv)):
            hamming = 0
            for k in range(9):
                if basisv[i][k] != basisv[j][k]:
                    hamming+=1
            
            if hamming == 1:
                H[i][j] = 1
        # print(i)
    print(len(basisv))
    dims = [[len(basisv)], [len(basisv)]]
    H = Qobj(H, dims)
    return H + H.dag(), target, psi0, basisv

def project_total_number_cube(V, C, E, H):
    # construct constraint satisfying states
    # find Hamiltonian action graph on these states
    # this can be done by directly computing excitation and backreaction in a loop
    pass

def project_total_number(V, C, E, H):
    # need to generate projectors for every constraint satisfying configuration
    # note in general this is.. hard but we can customize some solutions

    # generate permutations of 1, 2, 3
    # then generate all results where 1, 2, or all 3 are in the permutation
    # then clean up (or don't)
    # at most N * N! constraint satisfying configurations
    # any configuration of the graph fully determines occupation of constraint modes

    # or a more efficient way: for each choice of a vertex, now choose every other valid vertex and add those, so on..
    # can do this by pruning the graph one step at a time and looping over remaining configurations
    # for filled in range(1,3):
    #     for i in itertools.permutations(range(3), filled):
    adj = {0: [1, 2, 4], 1: [0, 5, 3], 2: [0, 3, 6], 3: [1, 2, 7], 4: [0, 5, 6], 5: [1, 4, 7], 6: [2, 4, 7], 7: [3, 5, 6]}
    configs = []
    # for a 2x2 cube
    for i in range(8):
        config = [0] * 8
        config[i] = 3
        configs.append(config)
        for j in range(8):
            if j in adj[i]:
                continue

            second = config.copy()
            second[j] = 3
            configs.append(second)
            for k in range(8):
                if k in adj[i] or k in adj[j]:
                    continue
                third = second.copy()
                third[k] = 3
                configs.append(third)

                for l in range(8):
                    if l in adj[i] or l in adj[j] or l in adj[k]:
                        continue
                    fourth = third.copy()
                    fourth[l] = 3
                    configs.append(fourth)
    for j in range(len(configs)):
        c = configs[j]
        cliques = [1] * C
        for i in range(len(c)):
            if c[i] == 3:
                for k in E[i]:
                    cliques[k] = 0
        configs[j] = cliques + c
        configs[j] = state_number_index([2] * C + [4] * V, configs[j])

    configs.append(state_number_index([2] * C + [4] * V, [1]*C + [0]*V))

    configs = list(set(configs))
    dims = [[2]*C + [4]*V, [2]*C + [4]*V]
    labels = [state_index_number(dims[0], i) for i in sorted(configs)]
    row = np.array(configs)
    col = np.array(configs)
    data = np.ones_like(row)
    dim = (2 ** C) * (4 ** V)
    projection_matrix = coo_matrix((data, (row, col)), shape=(dim, dim))

    projection = Qobj(projection_matrix, dims = dims)
    # print(configs)

    # configs = [tensor(ket(c[:C], 2), ket(c[C:], 4)) for c in configs]

    # print("working")

    # projectors = sum([projection(2**(V+C), c, c) for c in configs])
    H_reduced = projection*H*projection
    H_csr = H_reduced.data.tocsr()

    print("H")

    non_zero_rows = np.where(H_csr.getnnz(axis=1) > 0)[0]
    H_reduced = H_csr[non_zero_rows][:, non_zero_rows]
    H_reduced = Qobj(H_reduced, dims=dims)

    return H_reduced, labels

def draw_custom_node(ax, pos, node_data, n=2, flag=False):
    rect_size = 0.2
    text_size = 8
    if flag == False:
        rect_size = 0.02
        text_size = 2
        rect = plt.Rectangle((pos[0] - rect_size/2, pos[1] - rect_size/2), rect_size, rect_size, color='white', ec='black', zorder=2, clip_on=False)
        ax.add_patch(rect)
        return
    rect = plt.Rectangle((pos[0] - rect_size/2, pos[1] - rect_size/2), rect_size, rect_size, color='white', ec='black', zorder=2, clip_on=False)
    ax.add_patch(rect)

    # Split the node data into a nxn grid
    grid_data = np.array(node_data).reshape(n, n)

    # Calculate the size of each square in the grid
    square_size = rect_size / n

    # Draw the dividing lines for the grid
    # ax.plot([pos[0] - rect_size/2, pos[0] + rect_size/2], [pos[1], pos[1]], color='black', zorder=2)  # Horizontal line
    # ax.plot([pos[0], pos[0]], [pos[1] - rect_size/2, pos[1] + rect_size/2], color='black', zorder=2)  # Vertical line

    # Draw the numbers in the grid
    for i in range(n):
        for j in range(n):
            ax.text(pos[0] - rect_size/2 + square_size/2 + j*square_size, pos[1] + rect_size/2 - square_size/2 - i*square_size, grid_data[i, j], ha='center', va='center', fontsize=text_size, zorder=3)


def draw_hamiltonian(H, labels):
    print("Hello")
    H_mat = np.abs(H.full())
    print("Hello")
    G = nx.from_numpy_array(H_mat)
    G.remove_nodes_from(list(nx.isolates(G)))
    # labels = {i: f"{str(labels[i][:12])}\n{str(labels[i][12:])}" for i in G.nodes()}
    # nx.draw(G, with_labels=True, labels=labels, font_size=6)

    pos = nx.spring_layout(G, scale=1)
    pos = {node: (x*2, y*2) for node, (x, y) in pos.items()}
    plt.figure(figsize=(9,9))
    ax = plt.gca()
    ax.set_aspect('equal')
    nx.draw_networkx_edges(G, pos, ax=ax)

    for i in range(len(labels)):
        grid = labels[i][12:]
        print(grid)
        data_grid = [(grid[:4][i] + 2*grid[4:][i])//3 for i in range(4)]
        draw_custom_node(ax, pos[i], data_grid)

    plt.axis('off')
    plt.show()

def draw_sudoku(H, labels):
    print("Hello")
    H_mat = np.abs(H.full())
    print("Hello")
    G = nx.from_numpy_array(H_mat)
    G.remove_nodes_from(list(nx.isolates(G)))
    # labels = {i: f"{str(labels[i][:12])}\n{str(labels[i][12:])}" for i in G.nodes()}
    # nx.draw(G, with_labels=True, labels=labels, font_size=6)

    pos = nx.spring_layout(G, scale=1)
    pos = {node: (x*2, y*2) for node, (x, y) in pos.items()}
    plt.figure(figsize=(15,15))
    ax = plt.gca()
    ax.set_aspect('equal')
    nx.draw_networkx_edges(G, pos, ax=ax)

    for i in range(len(labels)):
        flag = False
        if 0 not in labels[i]:
            flag = True

        # uncomment to include dead-end states (kernel of forward scattering Hamiltonian)
        elif labels[i].count(0) <=2:
            zero_ = labels[i].index(0)
            f = False
            b = np.array(labels[i])
            b = np.reshape(b, (3,3))
            for j in range(1,3+1):
                if j not in b[zero_//3,:] and j not in b[:,zero_%3]:
                    f = True
            if f == False:
                flag = True
        
        draw_custom_node(ax, pos[i], labels[i], n=3, flag=flag)


    plt.axis('off')
    plt.show()

def stabilized_scar():
    H, target, psi0, basisv = sudoku_hamiltonian(initial={0:1, 2:2}, psi0vec=[[1,0,2,0,0,0,0,0,0]])
    final_ = drive_graph(H, target, 3, psi0, f"Hamming Distance")
    # draw_sudoku(H, labels=basisv)

    H, target, psi0, basisv = sudoku_hamiltonian(initial={0:1, 2:2}, psi0vec=[[1,0,2,0,0,0,0,0,0]])
    psi0 = final_
    for i in range(8):
        target = get_target_states(basisv, i)
        final_ = drive_graph(H, target, 3, psi0, f"Hamming Distance {i}")
    
    plt.title("3x3 Sudoku Target States Fidelity")
    plt.legend()
    plt.show()

def sudoku_4():
    # grid = np.array([[2, 4, 0, 3], [3, 0, 4, 2], [4, 2, 0, 1], [0, 0, 0, 0]])
    grid = np.array([[2, 4, 0, 3], [3, 0, 0, 2], [4, 0, 0, 1], [0, 0, 0, 0]])
    constraints = np.array([[[False] * 4] * 3] * 4)

    for k in range(4):
        for i in range(4):
            if (k+1) in grid[i]:
                constraints[k, 0, i] = True
            for j in range(4):
                if (k+1) in grid[:, j]:
                    constraints[k, 1, j] = True
                subgrid = grid[(i//2)*2:(i//2+1)*2, (j//2)*2:(j//2+1)*2]
                if any((k+1) in s for s in subgrid):
                    constraints[k,2,(i//2)*2 + (j//2)] = True


    sols = [grid]
    G = nx.Graph()
    sudoku_direct(4, grid=grid, constraints=constraints, depth=0, parent=0, sols=sols, G=G)
    nx.write_adjlist(G, "sudoku_graph_2.adjlist")
        
    G = nx.read_adjlist("sudoku_graph_2.adjlist")
    H = nx.adjacency_matrix(G)
    # dim = 7297
    dim = len(sols)
    # print(H)
    H = Qobj(H.todense(), [[dim], [dim]])
    H = H + H.dag()
    # 1163
    # max_v = 0
    # max_i = 0
    for i in range((grid == 0).sum()):
        drive_graph(H, get_target_states([j.flatten() for j in sols], i), 10, basis(dim, 0), name=i)
        # if m > max_v and i > 2:
        #     max_v = m
        #     max_i = i
        # print(m)

    # print(max_i, max_v)
    plt.legend()
    plt.show()

def sudoku_scaling_test(N):
    for n in range(3,N+1):
        grid = np.zeros((n, n))
        if n == 3:
            grid = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
        if n == 4:
            grid = np.array([[2, 4, 0, 3], [3, 0, 0, 2], [4, 0, 0, 1], [0, 0, 0, 0]])
        if n == 5:
            grid = np.array([[0, 2, 3, 4, 0], [2, 0, 1, 5, 3], [3, 0, 4, 2, 0], [0, 0, 0, 3, 0], [0, 3, 0, 1, 0]])

        constraints = np.array([[[False] * n] * 2] * n)

        for k in range(n):
            for i in range(n):
                if (k+1) in grid[i]:
                    constraints[k, 0, i] = True
                for j in range(n):
                    if (k+1) in grid[:, j]:
                        constraints[k, 1, j] = True

        sols = [grid]
        G = nx.Graph()
        sudoku_direct(n, grid=grid, constraints=constraints, depth=0, parent=0, sols=sols, G=G, subgrid=False)
        nx.write_adjlist(G, "sudoku_graph_5.adjlist")

        G = nx.read_adjlist("sudoku_graph_5.adjlist")
        H = nx.adjacency_matrix(G)

        dim = len(sols)
        H = Qobj(H.todense(), [[dim], [dim]])
        H = H + H.dag()

        drive_graph(H, get_target_states([j.flatten() for j in sols], 0), 10, basis(dim, 0), name=n)

    plt.legend()
    plt.show()

def draw_sudoku_4():
    grid = np.array([[2, 4, 0, 3], [3, 0, 0, 2], [4, 0, 0, 1], [0, 0, 0, 0]])
    constraints = np.array([[[False] * 4] * 3] * 4)

    for k in range(4):
        for i in range(4):
            if (k+1) in grid[i]:
                constraints[k, 0, i] = True
            for j in range(4):
                if (k+1) in grid[:, j]:
                    constraints[k, 1, j] = True
                subgrid = grid[(i//2)*2:(i//2+1)*2, (j//2)*2:(j//2+1)*2]
                if any((k+1) in s for s in subgrid):
                    constraints[k,2,(i//2)*2 + (j//2)] = True


    sols = [grid]
    G = nx.Graph()
    sudoku_direct(4, grid=grid, constraints=constraints, depth=0, parent=0, sols=sols, G=G)
    nx.write_adjlist(G, "sudoku_graph_3.adjlist")

    G = nx.read_adjlist("sudoku_graph_3.adjlist")
    H = nx.adjacency_matrix(G)
    # dim = 7297
    dim = len(sols)
    # print(H)
    H = Qobj(H.todense(), [[dim], [dim]])
    H = H + H.dag()

    # draw_sudoku(H, labels=)

# sudoku_4()
# start = time.time()
# V, E, cliques, psi0, target = rooks(3, initial=[1, 3])
# V, E, cliques, psi0, target = broken_cube()
# H = H_scar(V, E, cliques)
# print("H generated")
# drive_graph(H, [psi0, target], time=5, psi0=psi0)
# propagate_graph(H, [psi0, target], 1, psi0)
# print(time.time()-start)

# H, psi0, target = sudoku_basis(10)
# drive_graph(H, [psi0, target], time=5, psi0=psi0)

# PXP_rooks(3, 50)
    
# test_circle()
    
# test_initial()
    
# test_geometries()
# sudoku_basis(3)
# print(sudoku_basis(2))
# test_1d()

# test_rooks()
    
# V, E, cliques, psi0, target = rooks(3)
# H = H_scar(V, E, cliques)
# draw_hamiltonian(H)

# V, E, cliques, psi0, target = cube()
# H = H_scar(V, E, cliques)
# reduced_H, labels = project_total_number(V, E, cliques, H)
# print("Hi")
# draw_hamiltonian(reduced_H, labels)

# H, target, psi0, basisv = sudoku_hamiltonian(initial={0:1, 4: 2, 8: 3}, psi0vec=[[1,3,2,3,2,1,2,1,3]])
# drive_graph(H, target, 20, psi0, "3 Initial Conditions")


# H, target, psi0, basisv = sudoku_hamiltonian(initial={0:1}, psi0vec=[[1,3,2,3,2,1,2,1,3]])
# drive_graph(H, target, 20, psi0, "1 Initial Condition")
# H, target, psi0, basisv = sudoku_hamiltonian(psi0vec=[[1, 3, 2, 3, 2, 1, 2, 1, 3],
# [1, 3, 2, 3, 2, 1, 2, 1, 3],
# [1, 2, 3, 3, 1, 2, 2, 3, 1],
# [1, 3, 2, 2, 1, 3, 3, 2, 1],
# [1, 2, 3, 2, 3, 1, 3, 1, 2],
# [1, 3, 2, 3, 2, 1, 2, 1, 3]])
# drive_graph(H, target, 20, psi0, "0 Initial Condition")
# plt.title("3x3 Sudoku Target States Fidelity")
# plt.legend()
# plt.show()

# sudoku_basis(9, {10:1, 11:2, 13:3, 14:4, 15:5, 16:6, 17:7, 
#                  19:3, 20:4, 21:5, 23: 6, 24: 1, 25: 8, 26: 2, 
#                  29: 1, 31: 5, 32: 8, 33: 2, 35: 6, 
#                  38: 8, 39: 6, 44: 1, 
#                  46: 2, 50: 7, 52: 5, 
#                  56: 3, 57: 7, 59: 5, 61: 2, 62: 8, 
#                  64: 8, 67: 6, 69: 7, 
#                  72: 2, 74: 7, 76: 8, 77: 3, 78: 6, 79: 1, 80: 5})
    
# grid = np.array([[0]*9,
#            [0, 1, 2, 0, 3, 4, 5, 6, 7],
#            [0, 3, 4, 5, 0, 6, 1, 8, 2],
#            [0, 0, 1, 0, 5, 8, 2, 0, 6],
#            [0, 0, 8, 6, 0, 0, 0, 0, 1],
#            [0, 2, 0, 0, 0, 7, 0, 5, 0],
#            [0, 0, 3, 7, 0, 5, 0, 2, 8],
#            [0, 8, 0, 0, 6, 0, 7, 0, 0],
#            [2, 0, 7, 0, 8, 3, 6, 1, 5]])
# grid = np.array([[8, 0, 0, 2, 0, 0, 4, 0, 0], 
#                 [0, 0, 0, 0, 3, 0, 0, 9, 0],
#                 [4, 0, 0, 5, 0, 0, 0, 2, 0],
#                 [0, 2, 0, 3, 0, 0, 0, 5, 0],
#                 [0, 9, 0, 0, 0, 6, 0, 0, 1],
#                 [0, 0, 0, 0, 0, 5, 7, 0, 0],
#                 [0, 5, 0, 4, 0, 8, 0, 0, 0],
#                 [1, 0, 0, 0, 0, 0, 8, 0, 0],
#                 [0, 0, 0, 6, 7, 0, 0, 0, 5]])
# grid = np.array([[0, 0, 3, 0, 4, 0, 6, 0, 5], 
#                 [0, 1, 7, 8, 3, 0, 0, 0, 0],
#                 [5, 0, 4, 0, 7, 0, 1, 0, 0],
#                 [8, 0, 0, 7, 0, 4, 0, 9, 0],
#                 [0, 3, 1, 9, 0, 6, 0, 8, 0],
#                 [9, 7, 0, 5, 8, 0, 4, 0, 0],
#                 [0, 4, 0, 0, 0, 1, 7, 0, 6],
#                 [1, 6, 0, 0, 5, 7, 0, 0, 0],
#                 [0, 2, 0, 0, 0, 8, 0, 1, 3]])
    

# pos = nx.spring_layout(G, scale=1)
# pos = {node: (x*2, y*2) for node, (x, y) in pos.items()}
# plt.figure(figsize=(12,12))
# ax = plt.gca()
# ax.set_aspect('equal')
# nx.draw_networkx_edges(G, pos, ax=ax)
# print(pos)
# dim = len(pos.items())
# for i in range(dim):
#     flag = True
    
#     draw_custom_node(ax, pos[str(i)], [0]*9, n=3, flag=flag)

# plt.axis('off')
# plt.show()
    
sudoku_scaling_test(5)