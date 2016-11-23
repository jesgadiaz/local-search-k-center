import random
import datetime
import timeit
import math
import heapq
import os

################### FUNCTIONS ##############################################################

# Update Distance. Get the distance from every vertex to its nearest center in the partial solution
def update_distance(centers, max_index):
    global distance
    for i in range(0, n):
        if matrix[i][centers[max_index]] < distance[i]:
            distance[i] = matrix[i][centers[max_index]]

# Get the farthest node
def farthest_node():
    max_dist = 0
    max_dist_node = 0
    global distance
    for i in range(0,n):
        if distance[i] > max_dist:
            max_dist = distance[i]
            max_dist_node = i
    return max_dist_node


# Get the set of critical neighbors
def get_critical_neighbors(node, s):
    l = []
    for i in range(0, n):
        if matrix[i][node] < s:
            l.append(i)
    return l

# Get the solution size of the currently constructed solution
def solution_size():
    max_dist = 0
    global distance
    for i in range(0,n):
        if distance[i] > max_dist:
            max_dist = distance[i]
    return max_dist

def sol_size(sol):
    max_dist = 0
    for i in range(0, n):
        min_dist = float("inf")
        for j in range(0, k):
            if matrix[i][sol[j]] < min_dist:
                min_dist = matrix[i][sol[j]]
        if max_dist < min_dist:
            max_dist = min_dist
    return max_dist


# This function generates a 2-approximated solution with the Gon algorithm
def sdr_gon():
    C = []
    global distance, prob
    distance = []
    for i in range(0,n):
        distance.append(float("inf"))
    f = random.randint(0, n-1)
    C.append(f)
    update_distance(C, 0)
    for i in range(1, k):
        f = farthest_node()
        if random.random() <= prob:
            C.append(f)
        else:
            C.append(random.randint(0, n-1))
        update_distance(C, i)
    size = solution_size()
    out = [size, C]
    return out

# This function generates a 2-approximated solution with the Gon algorithm
def random_solution():
    C = []
    #for i in range(0, k):
    #    # e = random.randint(0, n-1)
    #    # while e in x_opt:
    #    #    e = random.randint(0, n-1)
    #    C.append(i)
    #    # x_curr.append(i)
    #random.shuffle(C)

    for i in range(0, k):
        e = random.randint(0, n - 1)
        while e in C:
            e = random.randint(0, n - 1)
        C.append(e)
    #for i in range(0, k):
    #    C.append(random.randint(0, n-1))
    size = sol_size(C)
    out = [size, C]
    return out

# Get the critical (farthest) node
def criticalNode(C):
    max_dist = 0
    farthest_node = 0
    for i in range(0, n):
        min_dist = float("inf")
        for j in range(0, k):
            if matrix[i][C[j]] <= min_dist:
                min_dist = matrix[i][C[j]]
        if min_dist > max_dist:
            max_dist = min_dist
            farthest_node = i
    return farthest_node

# Make the boolean version of the input solution
def getCurrentBoolean(in_sol):
    l = []
    for i in range(0, n):
        l.append(False)
    for i in range(0, len(in_sol)):
        l[in_sol[i]] = True
    return l

# Get First and Second closest centers
def firstAndSecond(in_sol):
    global F0
    global F1
    global D0
    global D1
    global M
    F0 = []
    F1 = []
    D0 = []
    D1 = []
    M = []
    for i in range(0, n):
        F0.append(0)
        F1.append(0)
        D0.append(0)
        D1.append(0)
        M.append(0)
    for i in range(0, n):
        min_dist = float("inf")
        for j in range(0, k):
            if matrix[i][in_sol[j]] < min_dist:
                min_dist = matrix[i][in_sol[j]]
                F0[i] = in_sol[j]
                D0[i] = min_dist
        min_dist = float("inf")
        for j in range(0, k):
            if matrix[i][in_sol[j]] < min_dist and in_sol[j] != F0[i]:
                min_dist = matrix[i][in_sol[j]]
                F1[i] = in_sol[j]
                D1[i] = min_dist

# Add facility to current solution
def addFacility(f):
    global F0
    global F1
    global D0
    global D1
    Sc = 0
    for i in range(0, n):
        if matrix[i][f] < D0[i]:
            D1[i] = D0[i]
            F1[i] = F0[i]
            D0[i] = matrix[i][f]
            F0[i] = f
        elif matrix[i][f] < D1[i]:
            D1[i] = matrix[i][f]
            F1[i] = f
        if D0[i] > Sc:
            Sc = D0[i]
    return Sc

# Remove facility to current solution
def removeFacility(f, current_sol):
    global F0
    global F1
    global D0
    global D1
    Sc = 0
    for i in range(0, n):
        if F0[i] == f:
            D0[i] = D1[i]
            F0[i] = F1[i]
            findNext(i, current_sol)
        elif F1[i] == f:
            findNext(i, current_sol)
        if D0[i] > Sc:
            Sc = D0[i]
    return Sc

# Find Next function
def findNext(f, current_sol):
    min_dist = float("inf")
    global D1
    global F1
    for i in range(0, k):
        if current_sol[i] != F0[f]:
            if matrix[current_sol[i]][f] < min_dist:
                F1[f] = current_sol[i]
                min_dist = matrix[current_sol[i]][f]
                D1[f] = min_dist

# Clean M array
def cleanMarray():
    global M
    M = []
    for i in range(0, n):
        M.append(0)


def remove_redundancies(L):
    out = []
    for i in range(0, len(L)):
        if L[i][0] != L[i][1]:
            out.append(L[i])
    return out



# Critical Interchange heuristic
def interchange(in_sol, current_size):
    x_opt = []
    x_curr = []
    temp_list = []
    for i in range(0, n):
        temp_list.append(i)
    temp_list = [x for x in temp_list if x not in in_sol]
    for i in range(0, k):
        x_opt.append(in_sol[i])
        x_curr.append(in_sol[i])
    for i in range(k, n):
        x_opt.append(temp_list[i-k])
        x_curr.append(temp_list[i-k])
    c1 = []
    c1_curr = []
    for i in range(0, n):
        min_dist = float("inf")
        nearest_center = 0
        for j in range(0, k):
            if matrix[i][x_opt[j]] < min_dist:
                min_dist = matrix[i][x_opt[j]]
                nearest_center = x_opt[j]
        c1.append(nearest_center)
        c1_curr.append(nearest_center)
    c2 = []
    c2_curr = []
    for i in range(0, n):
        min_dist = float("inf")
        second_nearest_center = 0
        for j in range(0, k):
            if x_opt[j] != c1[i]:
                if matrix[i][x_opt[j]] < min_dist:
                    min_dist = matrix[i][x_opt[j]]
                    second_nearest_center = x_opt[j]
        c2.append(second_nearest_center)
        c2_curr.append(second_nearest_center)
    f_opt = 0
    f_curr = 0
    i_star = 0
    i_star_curr = 0
    for i in range(0, n):
        if matrix[i][c1[i]] > f_opt:
            f_opt = matrix[i][c1[i]]
            i_star = i
            f_curr = matrix[i][c1[i]]
            i_star_curr = i
##################### BÜSQUEDA LOCAL ####################
    improve_locally = True
    while improve_locally:
        f_star = float("inf")
        # REDUCE THE NEIGHBORHOOD
        J = []
        for i in x_curr[k:n]:
            if matrix[i][i_star_curr] < f_curr:
                J.append(i)
        if len(J) == 0:
            for i in x_curr[k:n]:
                J.append(i)
        for new in J:
            # Initialization
            f = 0
            r = []
            z = []
            for i in range(0, n):
                r.append(0)
                z.append(0)

            # Add facility
            for i in range(0, n):
                if matrix[i][new] < matrix[i][c1_curr[i]]:
                    if matrix[i][new] > f:
                        f = matrix[i][new]
                else:
                    if matrix[i][c1_curr[i]] > r[c1_curr[i]]:
                        r[c1_curr[i]] = matrix[i][c1_curr[i]]
                    if matrix[i][new] < matrix[i][c2_curr[i]]:
                        temp = matrix[i][new]
                    else:
                        temp = matrix[i][c2_curr[i]]
                    if temp > z[c1_curr[i]]:
                        z[c1_curr[i]] = temp

            # Best deletion
            g1 = 0
            g2 = 0
            l_star = 0
            for i in range(0, k):
                if r[x_curr[i]] > g1:
                    g1 = r[x_curr[i]]
                    l_star = i
            for i in range(0, k):
                if r[x_curr[i]] > g2 and i != l_star:
                    g2 = r[x_curr[i]]

            min_g = float("inf")
            for i in range(0, k):
                if i != l_star:
                    temp = max([f, z[x_curr[i]], g1])
                    if temp < min_g:
                        min_g = temp
                        old = x_curr[i]
                else:
                    temp = max([f, z[x_curr[i]], g2])
                    if temp < min_g:
                        min_g = temp
                        old = x_curr[i]
            f = min_g
            if f < f_star:
                f_star = f
                new_star = new
                old_star = old
        if f_star >= f_curr:
            improve_locally = False
            #print(f_curr)
        else:
            f_curr = f_star
            index_old = x_curr.index(old_star)
            index_new = x_curr.index(new_star)
            x_curr[index_old] = new_star
            x_curr[index_new] = old_star
            for i in range(0, n):
                if c1_curr[i] == old_star:
                    if matrix[i][new_star] <= matrix[i][c2_curr[i]]:
                        c1_curr[i] = new_star
                    else:
                        c1_curr[i] = c2_curr[i]
                        min_dist = float("inf")
                        for j in range(0, k):
                            if matrix[i][x_curr[j]] < min_dist and x_curr[j] != c1_curr[i]:
                                min_dist = matrix[i][x_curr[j]]
                                c2_curr[i] = x_curr[j]
                else:
                    if matrix[i][c1_curr[i]] > matrix[i][new_star]:
                        c2_curr[i] = c1_curr[i]
                        c1_curr[i] = new_star
                    else:
                        if matrix[i][new_star] < matrix[i][c2_curr[i]]:
                            c2_curr[i] = new_star
                        else:
                            if c2_curr[i] == old_star:
                                min_dist = float("inf")
                                for j in range(0, k):
                                    if matrix[i][x_curr[j]] < min_dist and x_curr[j] != c1_curr[i]:
                                        min_dist = matrix[i][x_curr[j]]
                                        c2_curr[i] = x_curr[j]
            max_distance = 0
            for j in range(0, n):
                if matrix[j][c1_curr[j]] > max_distance:
                    max_distance = matrix[j][c1_curr[j]]
                    i_star_curr = j
    local_opt = []
    for i in range(0, k):
        local_opt.append(x_curr[i])
    return [local_opt, f_curr]

def critical_interchange_sa(in_sol, bks, current_size, T, k_sa, sa):
    #cardinality()
    #rep = 2*n
    #rep = int((n/10) * T)
    #max_rep = T
    T = 1
    max_rep = 2 * n * random.random()
    #rep = math.floor(random.random() * 2 * n)
    current_sol = []
    prev_mov_bad = True
    global F0
    global F1
    global D0
    global D1
    global M
    global card
    global best_known_size
    best_size = current_size
    tabu_movement = 0
    #tabu_list = []
    #for i in range(0, n):
    #    tabu_list.append(False)
    for i in range(0, n):
        list.append(False)
    for i in range(0, k):
        current_sol.append(in_sol[i])
    rep = 0
    while rep < max_rep:
    #for i in range(0, max_rep):
        rep =  rep + 1
        L = []
        # The critical node is also the farthest node
        current_sol_boolean = getCurrentBoolean(current_sol)
        critical_neighbors = get_critical_neighbors(criticalNode(current_sol), bks)
        firstAndSecond(current_sol)
        C = float("inf")
        for j in range(0, len(critical_neighbors)):
            Sc = addFacility(critical_neighbors[j])
            if j > 0:
                current_sol_boolean[critical_neighbors[j-1]] = False # Erase the previous added facility
            current_sol_boolean[critical_neighbors[j]] = True # Check with true the new added facility
            cleanMarray()
            for m in range(0, n):
                #if card[m] > 1:
                if matrix[critical_neighbors[j]][m] <= D1[m]:
                    if matrix[critical_neighbors[j]][m] > M[F0[m]]:
                        M[F0[m]] = matrix[critical_neighbors[j]][m]
                elif D1[m] > M[F0[m]]:
                    M[F0[m]] = D1[m]
                if Sc > M[F0[m]]:
                    M[F0[m]] = Sc
                #else:
                #    M[F0[m]] = float("inf")
            for m in range(0, k):
                if M[current_sol[m]] == 0:
                    M[current_sol[m]] = current_size
                #if not tabu_list[critical_neighbors[j]] or M[current_sol[m]] < best_size:
                if tabu_movement != critical_neighbors[j] or M[current_sol[m]] < best_size:
                    if M[current_sol[m]] == C:
                        L.append([current_sol[m], critical_neighbors[j]])
                    elif M[current_sol[m]] < C:
                        L = []
                        L.append([current_sol[m], critical_neighbors[j]])
                        C = M[current_sol[m]]
            Sc = removeFacility(critical_neighbors[j], current_sol)
            #remove redundant interchanges
            #L = remove_redundancies(L)
        #if len(L) == 0:
        #    tabu_list = []
        #    for i in range(0, n):
        #        tabu_list.append(False)
        #else:
        # perform the best change
        if C < current_size:
            selection = L[random.randint(0, len(L)-1)]
            for j in range(0, k):
                if current_sol[j] == selection[0]:
                    current_sol[j] = selection[1]
            current_size = C
            # If in the previous movement it got worst, then clean the tabu list
            #if prev_mov_bad:
            #    tabu_list = []
            #    for i in range(0, n):
            #        tabu_list.append(False)
            #    prev_mov_bad = False
        # perform any change
        else:
            selection = [current_sol[random.randint(0, k-1)], critical_neighbors[random.randint(0, len(critical_neighbors)-1)]]
            #selection = [current_sol[random.randint(0, k-1)], random.randint(0, n-1)]
            for j in range(0, k):
                if current_sol[j] == selection[0]:
                    current_sol[j] = selection[1]
            new_sol_size = sol_size(current_sol)
            if new_sol_size < current_size:
                current_size = new_sol_size
            else:
                if sa:
                    p_sa = math.pow(math.e, -((new_sol_size-current_size) / (k_sa * T)))
                else:
                    p_sa = -1
                if random.random() <= p_sa:
                    current_size = new_sol_size #sol_size(current_sol)
                    #tabu_list[selection[0]] = True
                    tabu_movement = selection[0]
                else:
                    if len(L) == 0:
                        current_size = new_sol_size
                    else:
                        for j in range(0, k):
                            if current_sol[j] == selection[1]:
                                current_sol[j] = selection[0]
                        selection = L[random.randint(0, len(L)-1)]
                        #tabu_list[selection[0]] = True
                        tabu_movement = selection[0]
                        for j in range(0, k):
                            if current_sol[j] == selection[0]:
                                current_sol[j] = selection[1]
                        current_size = C
                #prev_mov_bad = True
        if current_size < best_size:
            best_size = current_size
        if best_size < best_known_size:
            best_known_size = best_size
            #print(best_known_size)
            #cardinality()
        if best_known_size <= target[v-1]:
            break
    return best_known_size


    # Get First closest centers
def assignedCenter(in_sol):
    global F0
    F0 = []
    for i in range(0, n):
        F0.append(0)
        D0.append(0)
    for i in range(0, n):
        min_dist = float("inf")
        for j in range(0, k):
            if matrix[i][in_sol[j]] < min_dist:
                min_dist = matrix[i][in_sol[j]]
                F0[i] = in_sol[j]
                D0[i] = min_dist
        min_dist = float("inf")



################################ BEGINNING OF THE CODE ################################################################
cwd = os.getcwd()
f_out = open(cwd + '/out.txt', 'w')
# We text the algorithm over the set of 40 instances pmed from OR-Lib
#for v in range(9, 41):
#for v in [22, 11, 16, 21, 26, 31, 35, 38, 39]:
for v in [1, 6, 11, 16, 21, 26, 31, 35, 38]:
#for v in [1, 2, 3, 6, 7, 11 ,12, 16, 17, 21, 22, 26, 27, 31, 32, 35, 36, 38, 39]:
#for v in [17, 21, 22, 26, 27, 31, 32, 35, 36, 38, 39]:
#for v in [27]:
    # The path of every instance is placed here (pmed1, pmed2,..., pmed40)
    f = open(cwd + '/Lib/pmed' + repr(v) + '.txt', 'r')
    # These target values are the optimal solutions' size
    target = [127, 98, 93, 74, 48, 84, 64, 55, 37, 20, 59, 51, 35, 26, 18, 47, 39, 28, 18, 13, 40, 38, 22, 15, 11, 38,
              32, 18, 13, 9, 30, 29, 15, 11, 30, 27, 15, 29, 23, 13]
    #target = [127, 98, 93, 74, 48, 84, 64, 55, 37, 20, 59, 51, 35, 26, 18, 47, 39, 29, 19, 14, 40, 38, 23, 15, 11, 38,
    #          32, 19, 13, 10, 30, 30, 16, 11, 30, 28, 16, 29, 24, 14]

    # Setting up the instance (BEGIN)***************************************************************
    # Get the values of n, m and k
    str = f.readline()
    str = str.split()
    n = int(str[0])
    m = int(str[1])
    k = int(str[2])
    # Fill the Adjacencies Matrix
    matrix = []
    for i in range(0,n):
        list = []
        for j in range(0,n):
            list.append(float("inf"))
        matrix.append(list)
    for i in range(0, m):
        str = f.readline()
        str = str.split()
        v1 = int(str[0]) - 1
        v2 = int(str[1]) - 1
        weight = float(str[2])
        matrix[v1][v2] = weight
        matrix[v2][v1] = weight
    # Apply the Floyd-Marshal algorithm
    for i in range(0, n):
        matrix[i][i] = 0
    for i in range(0, n):
        for j in range(0, n):
            for l in range(0, n):
                if matrix[i][j] == float("inf") or matrix[i][l] == float("inf"):
                    cost = float("inf")
                else:
                    cost = matrix[i][j] + matrix[i][l]
                if cost < matrix[j][l]:
                    matrix[j][l] = cost
    f.close()
    # Setting up the instance (END)***************************************************************


    # The algorithm can be executed over every instance any number of times with a different seed
    #print("*************************************************************************************************")
    print("100 executions of pmed" + repr(v) + ": ")
    f_out.write("100 executions of pmed" + repr(v) + ": " + '\n')
    for h in range(0, 10):
    #for h in range(17, 40):
        random.seed(h)
        # These variables are used by the Local Search method "Critical Interchange"
        F0 = []
        F1 = []
        D0 = []
        D1 = []
        M = []
        ds_time = 0
        ls_time = 0
        sa = True
        done = False

        interchange_size = float("inf")

        k_sa = 100
        T = 1

        is_local_opt = False

        init_ds_time = timeit.default_timer()
        ds = random_solution()
        size = ds[0]
        #C = ds[1]
        best_known_size = size
        best_sol = []
        for i in range(0, k):
            best_sol.append(ds[1][i])

        init_ls_time = timeit.default_timer()
        current_sol = []
        prev_mov_bad = True
        best_size = best_known_size
        current_size = best_known_size
        tabu_movement = 0
        list = []
        for i in range(0, n):
            list.append(False)
        for i in range(0, k):
            current_sol.append(best_sol[i])
        while not done and timeit.default_timer() - init_ds_time < 1000:
            neighbors = []
            for j in range(0, n):
                if j not in current_sol:
                    neighbors.append(j)

            selection = [current_sol[random.randint(0, k - 1)],
                         neighbors[random.randint(0, len(neighbors) - 1)]]

            for j in range(0, k):
                if current_sol[j] == selection[0]:
                    current_sol[j] = selection[1]
            new_sol_size = sol_size(current_sol)

            if new_sol_size < current_size:
                current_size = new_sol_size
            else:
                p_sa = math.pow(math.e, -((new_sol_size - current_size) / (k_sa * T)))
                if random.random() <= p_sa:
                    current_size = new_sol_size
                else:
                    current_sol = current_sol

#                    interchange_size = interchange(current_sol, current_size)

            if current_size < best_size:
                best_size = current_size
            #if interchange_size < best_size:
            #    best_size = interchange_size
            if best_size < best_known_size:
                best_known_size = best_size
                #print(best_known_size)
                # cardinality()
            if best_known_size <= target[v - 1]:
                done = True

        # Here we are printing the execution time until reaching the target
        #time = ds_time + ls_time
        time = timeit.default_timer() - init_ds_time
        if time >= 1000:
            time = 1000
        print(repr(time))
        f_out.write(repr(time) + '\n')
        #print(best_known_size)
f_out.close()