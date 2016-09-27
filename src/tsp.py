import math

# Nearest Neighbor法
def nn(lst):
    print("Start Nearest Neighbor Method.")
    route,remaining_points = [lst[0]],lst[1:]
    i = 0
    while (len(remaining_points) > 0):
        if i % 1000 ==0:
            print("{} points are connected.".format(i))
        last = route[-1]
        near = find_nearest(last, remaining_points)
        route.append(near)
        remaining_points.remove(near)
        i+=1
    print("{} points are connected.".format(i))
    print("Completed Nearest Neghbor Method.")
    return route

# 距離の計算
def dist(pt1, pt2):
    return math.sqrt((pt1[0] - pt2[0])**2 + (pt1[1]-pt2[1])**2)

# 最も近い地点を探すメソッド
def find_nearest(pt, lst):
    best = float('inf')
    best_index = -1
    for i in range(len(lst)):
        if dist(lst[i], pt) < best:
            best = dist(lst[i], pt)
            best_index = i
    return lst[best_index]

def calc_totaldist(lst):
    tdist = 0
    edges = [[lst[i],lst[i+1]] for i in range(len(lst)-1)] + [(lst[-1], lst[0])]
    for e in edges:
        tdist += dist(e[0],e[1])
    return tdist

# 2-opt法
def opt2(route):
    print("Start 2-Opt Method.")
    size = len(route)
    print(size)
    total_dist = calc_totaldist(route)
    print("Current total distance is {}".format(total_dist))
    improve = 21
    itr = 1
    while improve > 20:
        improve = 0
        for i in range(0,size - 2):
            nexti = i+1
            for j in range(i + 2, size):
                if j == size - 1:
                    nextj = 0
                else:
                    nextj = j + 1
                change = dist(route[i],route[j]) + dist(route[nexti],route[nextj]) - dist(route[i],route[nexti]) - dist(route[j],route[nextj])
                if change < 0:
                    swap_edges(nexti,j+1,route)
                    improve += -1*change
        total_dist = total_dist + -1*improve
        print("Finished {} loop. Current total distance is {}.({} reduced)".format(itr,total_dist,improve))
        itr += 1
    print("Completed 2-Opt.")
    return route

def swap_edges(ii,jj,route):
    new_subroute = route[ii:jj]        
    route[ii:jj] = new_subroute[::-1]
