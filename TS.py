import random
import datetime
import timeit
import math
import os

################### FUNCTIONS ##############################################################




################################ BEGINNING OF THE CODE ################################################################
cwd = os.getcwd()
#for v in range(1, 41):
#for v in [1, 2, 3, 6, 7, 11, 12, 16, 17, 21, 22, 26, 27, 31, 32, 35, 36, 38, 39]:
for v in [1, 6, 11, 16, 21, 26, 31, 35, 38]:
#for v in [38]:

    #print("Loading input....................................................................")
    # Open the selected file
    #f = open('/home/jgd/ORLib/u1817_5.txt', 'r')
    #f = open('/home/jgd/ORLib/pmed11.txt', 'r')
    f = open(cwd + '/Lib/pmed' + repr(v) + '.txt', 'r')
    #f = open('/home/jgd/ORLib/pmed' + repr(v) + '.txt', 'r')
    #f = open('C:/Users/netcomlab1/Desktop/PAPER/ORLib/pmed' + repr(v) + '.txt', 'r')
    #f = open('/home/jgd/ORLib/test.txt', 'r')
    #target = 93
    target = [127, 98, 93, 74, 48, 84, 64, 55, 37, 20, 59, 51, 35, 26, 18, 47, 39, 28, 18, 13, 40, 38, 22, 15, 11, 38,
              32, 18, 13, 9, 30, 29, 15, 11, 30, 27, 15, 29, 23, 13]
    #target = [127, 98, 93, 74, 48, 84, 64, 55, 37, 20, 59, 51, 36, 26, 18, 47, 39, 29, 19, 14, 40, 38, 23, 15, 11, 38,
    #          32, 19, 13, 10, 30, 30, 16, 11, 30, 28, 16, 29, 24, 14]

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

    f.close()

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
    #print("Input loaded")



    # The output will be written in a text file
    #out = open("output", "w")

    n_tl = 2

    max_time = 1000

    print("*************************************************************************************************")
    print("100 executions of pmed" + repr(v) + ": ")
    # SdR-DS
    for h in range(0, 100):
        # INITIALIZATION
        random.seed(h)
        init_time = timeit.default_timer()
        x_opt = []
        x_curr = []
        for i in range(0, n):
            e = random.randint(0, n-1)
            while e in x_opt:
                e = random.randint(0, n-1)
            x_opt.append(e)
            x_curr.append(e)
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
        cn1 = k-1
        cn2 = n-1
        t1 = math.floor(1 + (k-1) * random.random())
        t2 = math.floor(1 + (n-k-1) * random.random())
        improve = True

        old_star = 0
        new_star = 0

        max_time_reached = False

        while f_opt > target[v-1] and not max_time_reached:
        # MAIN STEP
            # REDUCE THE NEIGHBORHOOD
            J = []
            for i in x_curr[k:n-t2]:
                if matrix[i][i_star_curr] < f_curr:
                    J.append(i)
            if len(J) == 0 or not improve:
                for i in x_curr[k:n]:
                    if i not in x_curr[n-t2:n]:
                        J.append(i)

            # FIND THE BEST IN THE NEIGHBORHOOD
            f_star = float("inf")
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
                    if matrix[i][new] < matrix[i][c1[i]]:
                        if matrix[i][new] > f:
                            f = matrix[i][new]
                    else:
                        if matrix[i][c1[i]] > r[c1[i]]:
                            r[c1[i]] = matrix[i][c1[i]]
                        if matrix[i][new] < matrix[i][c2[i]]:
                            temp = matrix[i][new]
                        else:
                            temp = matrix[i][c2[i]]
                        if temp > z[c1[i]]:
                            z[c1[i]] = temp

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

            # Improvement
            if f_star < f_opt:
                f_opt = f_star
                improve = True
                x_opt = x_curr
                if f_opt <= target[v-1]:
                    exec_time = timeit.default_timer() - init_time
                    print(exec_time)
            else:
                improve = False
                exec_time = timeit.default_timer() - init_time
                if exec_time > max_time:
                    print(max_time)
                    max_time_reached = True

            # Updating
            f_curr = f_opt
            max_dist = 0
            for i in range(0, n):
                min_dist = float("inf")
                for j in range(0, k):
                    if matrix[i][x_curr[j]] < min_dist:
                        min_dist = matrix[i][x_curr[j]]
                if min_dist > max_dist:
                    max_dist = min_dist
                    i_star = i
                    i_star_curr = i_star
            if n_tl == 1:
                a = old_star
                index_old = x_curr.index(old_star)
                index_new = x_curr.index(new_star)
                x_curr[index_old] = new_star
                x_curr[index_new] = x_curr[cn2]
                x_curr[cn2] = a
            elif n_tl == 2:
                a = old_star
                index_old = x_curr.index(old_star)
                index_new = x_curr.index(new_star)
                x_curr[index_old] = x_curr[cn1]
                x_curr[cn1] = new_star
                x_curr[index_new] = x_curr[cn2]
                x_curr[cn2] = a
            cn1 = cn1 - 1
            if cn1 == (k-1) -t1:
            #if cn1 == k - t1:
                cn1 = k-1
            cn2 = cn2 - 1
            #if cn2 == n - t2:
            if cn2 == n-1 - t2:
                cn2 = n-1

            for i in range(0, n):
                if c1[i] == old_star:
                    if matrix[i][new_star] <= matrix[i][c2[i]]:
                        c1[i] = new_star
                    else:
                        c1[i] = c2[i]
                        min_dist = float("inf")
                        for j in range(0, k):
                            if matrix[i][x_curr[j]] < min_dist and x_curr[j] != c1[i]:
                                min_dist = matrix[i][x_curr[j]]
                                c2[i] = x_curr[j]
                else:
                    if matrix[i][c1[i]] > matrix[i][new_star]:
                        c2[i] = c1[i]
                        c1[i] = new_star
                    else:
                        if matrix[i][new_star] < matrix[i][c2[i]]:
                            c2[i] = new_star
                        else:
                            if c2[i] == old_star:
                                min_dist = float("inf")
                                for j in range(0, k):
                                    if matrix[i][x_curr[j]] < min_dist and x_curr[j] != c1[i]:
                                        min_dist = matrix[i][x_curr[j]]
                                        c2[i] = x_curr[j]



