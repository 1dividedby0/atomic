
def cube():
    cliques = [[0, 4, 8], [0, 5, 9], [1, 4, 10], [1, 5, 11], [2, 6, 8], [2, 7, 9], [3, 6, 10], [3, 7, 11]]
    psi0 = tensor(ket([1] * 12, 2), ket([0]*8, 4))
    sols = [[3, 0, 0, 3, 0, 3, 3, 0], [0, 3, 3, 0, 3, 0, 0, 3]]
    sols = [tensor(ket([0]*12, 2), ket(s, 4)) for s in sols]
    return 8, 12, cliques, psi0, sols

def cubeplus():
    cliques = [[0, 2, 5], [1, 2, 6], [0, 3, 7], [1, 3, 8], [0, 4, 9], [1, 4, 10], 
               [11, 13, 5], [12, 13, 6], [11, 14, 7], [12, 14, 8], [11, 15, 9], [12, 15, 10]]
    # psi0 = tensor(ket([1] * 12, 2), ket([0]*8, 4))
    # sols = [[3, 0, 0, 3, 0, 3, 3, 0], [0, 3, 3, 0, 3, 0, 0, 3]]
    # sols = [tensor(ket([0]*12, 2), ket(s, 4)) for s in sols]
    # target = sum(sols).unit()
    return 12, 15, cliques, psi0, target

def broken_cube(initial=None):
    cliques = [[0, 3, 6], [0, 4, 7], [1, 3, 8], [1, 4], [2, 5, 6], [2, 7], [5, 8]]
    psi0 = 0
    if initial == None:   
        psi0 = tensor(ket([1] * 9, 2), ket([0]*3, 4), ket([0], 3), ket([0], 4), ket([0]*2, 3))
    else:
        bank = [1]*9

        for ic in initial:
            for j in cliques[ic]:
                bank[j] = 0
        
        psi0 = tensor(ket(bank, 2), ket([0]*3, 4), ket([0], 3), ket([0], 4), ket([0]*2, 3))

    sols = [[3, 0, 0, 2, 0, 2, 2], [0, 3, 3, 0, 3, 0, 0]]
    target = []
    for j in range(len(sols)):
        if initial:
            skip = False
            for ic in initial:
                # if initial condition not satisfied for this solution, then don't check overlap
                if sols[j][ic] == 0:
                    skip = True
                
                # make sure initial conditions are set to 0 since we have fewer particles
                elif sols[j][ic] > 0:
                    sols[j][ic] = 0

            if skip == True:
                continue

        prod = tensor(ket([0]*9, 2), ket(sols[j][:3], 4), ket([sols[j][3]], 3), ket([sols[j][4]], 4), ket(sols[j][5:7], 3))
        target.append(prod)

    return 7, 9, cliques, psi0, target

def ladder():
    # first row, then second row
    cliques = [[0, 4], [0, 1, 5], [1, 6], [4, 2], [2, 3, 5], [3, 6]]
    psi0 = tensor(ket([1] * 7, 2), ket([0], 3), ket([0],4), ket([0], 3), ket([0], 3), ket([0], 4), ket([0], 3))
    sols = [[2, 0, 2, 0, 3, 0], [0, 3, 0, 2, 0, 2]]
    sols = [tensor(ket([0]*7, 2), ket([s[0]], 3), ket([s[1]],4), ket([s[2]], 3), ket([s[3]], 3), ket([s[4]], 4), ket([s[5]], 3)) for s in sols]
    return 6, 7, cliques, psi0, sols

# def rooks():
#     cliques = [[0, 3], [0, 4], [0, 5], [1, 3], [1, 4], [1, 5], [2, 3], [2, 4], [2, 5]]
#     psi0 = tensor(ket([1, 1, 1, 1, 1, 1], 2), ket([0, 0, 0, 0, 0, 0, 0, 0, 0], 3))
#     sol1 = tensor(ket([0, 0, 0, 0, 0, 0], 2), ket([2, 0, 0, 0, 2, 0, 0, 0, 2], 3))
#     sol2 = tensor(ket([0, 0, 0, 0, 0, 0], 2), ket([2, 0, 0, 0, 0, 2, 0, 2, 0], 3))
#     sol3 = tensor(ket([0, 0, 0, 0, 0, 0], 2), ket([0, 2, 0, 2, 0, 0, 0, 0, 2], 3))
#     sol4 = tensor(ket([0, 0, 0, 0, 0, 0], 2), ket([0, 2, 0, 0, 0, 2, 2, 0, 0], 3))
#     sol5 = tensor(ket([0, 0, 0, 0, 0, 0], 2), ket([0, 0, 2, 2, 0, 0, 0, 2, 0], 3))
#     sol6 = tensor(ket([0, 0, 0, 0, 0, 0], 2), ket([0, 0, 2, 0, 2, 0, 2, 0, 0], 3))
#     target = (sol1 + sol2 + sol3 + sol4 + sol5 + sol6).unit()
#     return 9, 6, cliques, psi0, target

# torus
def rooks(N, initial=None):
    cliques = [[i%N, N + i//N] for i in range(N**2)]
    print(cliques)
    psi0 = 0
    if initial == None:
        psi0 = tensor(ket([1]*2*N, 2), ket([0]*(N**2), 3))
    else:
        bank = [1]*2*N
        sites = [0]*(N**2)

        for ic in initial:
            bank[ic%N] = 0
            bank[N + ic//N] = 0
        
        psi0 = tensor(ket(bank, 2), ket(sites, 3))

    target = []
    # print(itertools.permutations(list(range(N))))
    for i in itertools.permutations(list(range(N))):
        vec = []
        for j in range(N):
            row = [0]*N
            row[i[j]] = 2
            vec += row
        
        if initial:
            skip = False
            for ic in initial:
                # if initial condition not satisfied for this solution, then don't check overlap
                if vec[ic] == 0:
                    skip = True
                
                # make sure initial conditions are set to 0 since we have fewer particles
                elif vec[ic] == 2:
                    vec[ic] = 0

            if skip == True:
                continue
        
        print(vec)
        # print(vec)
        prod = tensor(ket([0]*2*N, 2), ket(vec, 3))
        target.append(prod)

    return N**2, 2*N, cliques, psi0, target

def PXP_rooks(N, time):
    H = 0
    for i in range(N):
        h = [identity(2)] * (N**2)
        h[i] = sigmax()
        neighbors = []
        if i%N == N-1:
            neighbors.append(i-1)
        elif i%N == 0:
            neighbors.append(i+1)
        else:
            neighbors.append(i-1)
            neighbors.append(i+1)
        
        if i//N == N-1:
            neighbors.append(i-N)
        elif i//N == 0:
            neighbors.append(i+N)
        else:
            neighbors.append(i-N)
            neighbors.append(i+N)

        for n in neighbors:
            h[n] = projection(2, 0, 0)

        if H == 0:
            H = tensor(h)
        else:
            H += tensor(h)

    target = 0

    for i in itertools.permutations(list(range(N))):
        vec = []
        for j in range(N):
            row = [0]*N
            row[i[j]] = 1
            vec += row

        prod = ket(vec, 2)
        if target == 0:
            target = prod
        else:
            target += prod

    target = target.unit()

    # psi0 = ket([0]*(N**2), 2)
    psi0 = target

    step = 0.01
    t = np.arange(0, time, step)
    results = krylovsolve(H, psi0, t, krylov_dim=20, sparse=True)

    vecs = [psi0, target]
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


def even_chain(N):
    cliques = [[0]] + [[i-1, i] for i in range(1,N-1)] + [[N-2]]
    psi0 = tensor(ket([1]*(N-1), 2), ket([0], 2), ket([0]*(N-2), 3), ket([0], 2))
    neel1 = 0
    neel2 = 0

    if N%2 == 0:
        neel1 = [1] + [0, 2] * (N//2 - 1) + [0]
        neel2 = [0] + [2, 0] * (N//2 - 1) + [1]

    neel1 = tensor(ket([0]*(N-1), 2), ket([neel1[0]], 2), ket(neel1[1:N-1], 3), ket([neel1[N-1]], 2))
    neel2 = tensor(ket([0]*(N-1), 2), ket([neel2[0]], 2), ket(neel2[1:N-1], 3), ket([neel2[N-1]], 2))
    
    target = [neel1, neel2]

    return N, N-1, cliques, psi0, target

def circle(N):
    cliques = [[N-1, 0]] + [[i-1, i] for i in range(1,N)]
    psi0 = tensor(ket([1]*N, 2), ket([0]*N, 3))
    neel1 = 0
    neel2 = 0

    if N%2 == 0:
        neel1 = [2, 0] * (N//2)
        neel2 = [0, 2] * (N//2)

    elif N%2 == 1:
        neel1 = [2, 0] * (N//2) + [0]
        neel2 = [0, 2] * (N//2) + [0]

    neel1 = tensor(ket([0]*N, 2), ket(neel1, 3))
    neel2 = tensor(ket([0]*N, 2), ket(neel2, 3))
    
    target = [neel1, neel2]

    return N, N, cliques, psi0, target

def onehot(N):
    cliques = [[0]]*N
    psi0 = tensor(ket([1], 2), ket([0]*N, 2))
    target = [[0]*i + [1] + [0]*(N-i-1) for i in range(N)]

    return N, 1, cliques, psi0, [tensor(ket([0], 2), tensor(ket(t, 2))) for t in target]

def kagome(initial_condition=False):
    cliques = [[0], [0, 1, 2], [0, 2, 3], [1], [1, 2, 3], [3]]
    psi0 = 0
    target = 0
    if initial_condition:
        psi0 = tensor(ket([0, 1, 1, 1],2), ket([0],2), ket([0, 0],4), ket([0], 2), ket([0],4), ket([0],2))
        target = tensor(ket([0, 0, 1, 0], 2), ket([0],2), ket([0, 0],4), ket([1], 2), ket([0],4), ket([1],2))
    else:
        psi0 = tensor(ket([1]*4,2), ket([0],2), ket([0, 0],4), ket([0], 2), ket([0],4), ket([0],2))
        target = tensor(ket([0, 0, 1, 0], 2), ket([1],2), ket([0, 0],4), ket([1], 2), ket([0],4), ket([1],2))
    

    return 6, 4, cliques, psi0, [target]

def extended_kagome(initial_condition=False):
    cliques = [[0, 1], [1], [0, 2, 3], [0, 1, 3, 4], [2], [2, 3, 4], [4]]
    psi0 = 0
    target = 0
    if initial_condition:
        psi0 = tensor(ket([0, 1, 1, 1, 1],2), ket([0],3), ket([0], 2), ket([0],4), ket([0], 5), ket([0], 2), ket([0],4), ket([0],2))
        target = [tensor(ket([0, 0, 1, 0, 0], 2), ket([0],2), ket([0, 0],4), ket([1], 2), ket([0],4), ket([1],2))]
    else:
        psi0 = tensor(ket([1]*5,2), ket([0],3), ket([0], 2), ket([0],4), ket([0], 5), ket([0], 2), ket([0],4), ket([0],2))
        first = tensor(ket([0, 0, 0, 1, 0], 2), ket([2],3), ket([0], 2), ket([0],4), ket([0], 5), ket([1], 2), ket([0],4), ket([1],2))
        second = tensor(ket([0, 0, 0, 0, 0], 2), ket([0],3), ket([1], 2), ket([3],4), ket([0], 5), ket([0], 2), ket([0],4), ket([1],2))
        third = tensor(ket([1, 0, 0, 1, 0], 2), ket([0],3), ket([1], 2), ket([0],4), ket([0], 5), ket([1], 2), ket([0],4), ket([1],2))
        # fourth = tensor(ket([0, 0, 0, 0, 0], 2), ket([2],3), ket([0], 2), ket([0],4), ket([0], 5), ket([0], 2), ket([3],4), ket([0],2))
        # fifth = tensor(ket([1, 0, 0, 0, 0], 2), ket([0],3), ket([1], 2), ket([0],4), ket([0], 5), ket([0], 2), ket([3],4), ket([0],2))
        target = [first, second, third]
    
    return 7, 5, cliques, psi0, target

def test_geometries():
    drive_time = 1

    V, E, cliques, psi0, target = kagome()
    H = H_scar(V, E, cliques)
    m6 = drive_graph(H, target, time=drive_time, psi0=psi0, name="kagome (2/3)")

    # V, E, cliques, psi0, target = extended_kagome()
    # H = H_scar(V, E, cliques)
    # m6 = drive_graph(H, target, time=drive_time, psi0=psi0, name="extended kagome (5/7)")

    V, E, cliques, psi0, target = ladder()
    H = H_scar(V, E, cliques)
    m1 = drive_graph(H, target, time=drive_time, psi0=psi0, name="ladder (7/6)")

    V, E, cliques, psi0, target = circle(6)
    H = H_scar(V, E, cliques)
    m2 = drive_graph(H, target, time=drive_time, psi0=psi0, name="circle (1)")

    V, E, cliques, psi0, target = rooks(3)
    H = H_scar(V, E, cliques)
    m3 = drive_graph(H, target, time=drive_time, psi0=psi0, name="rooks (2/3)")

    V, E, cliques, psi0, target = rooks(2)
    H = H_scar(V, E, cliques)
    m3 = drive_graph(H, target, time=drive_time, psi0=psi0, name="rooks (1)")

    V, E, cliques, psi0, target = broken_cube()
    H = H_scar(V, E, cliques)
    m4 = drive_graph(H, target, time=drive_time, psi0=psi0, name="broken cube (9/7)")

    V, E, cliques, psi0, target = onehot(3)
    H = H_scar(V, E, cliques)
    m5 = drive_graph(H, target, time=drive_time, psi0=psi0, name="cluster (1/3)")

    # plt.plot([1/3, 2/3, 1, 7/6, 9/7], [m5, m3, m2, m1, m4])

    plt.title("Target State Fidelity for Different Geometries")
    plt.xlabel("Time")
    plt.ylabel("Target State Fidelity")

    plt.legend()
    plt.show()

def test_circle():
    drive_time = 1
    V, E, cliques, psi0, target = circle(8)
    H = H_scar(V, E, cliques)
    m8 = drive_graph(H, [target], time=drive_time, psi0=psi0, name="circle (8)", plot=False)
    V, E, cliques, psi0, target = circle(6)
    H = H_scar(V, E, cliques)
    m6 = drive_graph(H, [target], time=drive_time, psi0=psi0, name="circle (6)", plot=False)
    V, E, cliques, psi0, target = circle(4)
    H = H_scar(V, E, cliques)
    m4 = drive_graph(H, [target], time=drive_time, psi0=psi0, name="circle (4)", plot=False)
    plt.plot([4, 6, 8], [m4, m6, m8])
    plt.title("Maximum Overlap with Target State for 1D Chain with PBC")
    plt.xlabel("# Vertices")
    plt.ylabel("Max Fidelity")
    plt.legend()
    plt.show()

def test_initial():
    drive_time = 3
    V, E, cliques, psi0, target = kagome(initial_condition=True)
    H = H_scar(V, E, cliques)
    drive_graph(H, [target], time=drive_time, psi0=psi0, name="Initial Condition")

    V, E, cliques, psi0, target = kagome(initial_condition=False)
    H = H_scar(V, E, cliques)
    drive_graph(H, [target], time=drive_time, psi0=psi0, name="No Initial Condition")

    plt.title("Target State Fidelity for Nearest Neighbor Triangle Lattice")
    plt.xlabel("Time")
    plt.ylabel("Target State Fidelity")

    plt.legend()
    plt.show()

def test_rooks():
    drive_time = 2
    V, E, cliques, psi0, target = rooks(3)
    H = H_scar(V, E, cliques)
    drive_graph(H, [target], time=drive_time, psi0=psi0, name="0 Initial Conditions")

    V, E, cliques, psi0, target = rooks(3, initial=[1])
    H = H_scar(V, E, cliques)
    drive_graph(H, [target], time=drive_time, psi0=psi0, name="1 Initial Condition")

    V, E, cliques, psi0, target = rooks(3, initial=[1, 3])
    H = H_scar(V, E, cliques)
    drive_graph(H, [target], time=drive_time, psi0=psi0, name="2 Initial Conditions")

    plt.title("Target State Fidelity for 3x3 Rooks Puzzle")
    plt.ylabel("Target State Fidelity")
    plt.xlabel("Time")
    plt.legend()
    plt.show()

def test_1d():
    drive_time = 2
    V, E, cliques, psi0, target = even_chain(8)
    H = H_scar(V, E, cliques)
    drive_graph(H, [target], time=drive_time, psi0=psi0, name="Nearest Neighbor")

    V, E, cliques, psi0, target = circle(8)
    H = H_scar(V, E, cliques)
    drive_graph(H, [target], time=drive_time, psi0=psi0, name="Nearest Neighbor with PBC")

    V, E, cliques, psi0, target = onehot(8)
    H = H_scar(V, E, cliques)
    drive_graph(H, [target], time=drive_time, psi0=psi0, name="Handshakes at a Party")

    plt.title("Target State Fidelity for 1D Constraint Configurations")
    plt.ylabel("Target State Fidelity")
    plt.xlabel("Time")
    plt.legend()
    plt.show()
