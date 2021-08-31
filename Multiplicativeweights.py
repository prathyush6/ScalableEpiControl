import networkx as nx
import random
import math

def MWU(G, S, epsilon, delta_val, lambda_val, M, foldername):
    #initialize ell values
    ell = {} 
    for u in G.nodes():
        if u in S:
           ell[u] = 0
        else:
           ell[u] = delta_val
    #define capacity of nodes
    c = {}
    for u in G.nodes():
        if u in S:
           c[u] = 0
        else:
           c[u] = lambda_val
 
    A = set() #copies of nodes. eg. a(u, j) is a copy of u attached to it in sample H_j
  
    s = "-1"
    #Iterations
    max_itr = math.log((1+epsilon)/delta_val, 1+epsilon)
    max_itr = math.floor(max_itr)
    spcomp = 0
    r = 1
    print("Max no. of iterations "+str(max_itr), flush = True)
    #for r in range(1, max_itr+1):
    H = [] 
    while r < max_itr+1:
        #load samples into memoery
        if r % 2 == 1:
           H.clear()
           rs = r
           for j in range(rs, rs+2):
               fp = open(foldername+"/sample"+str(j)+".txt", 'r')
               lines = fp.readlines()
               Hj = nx.Graph()
               for line in lines:
                   cols = line.strip().split("\t")
                   u = cols[0]
                   v = cols[1]
                   Hj.add_edge(u, v)
               H.append(Hj)
               fp.close()
        rq = r % 2 -1
        Hr, Ar = generateDiGraph(H[rq], M, S, s, ell, c, delta_val)
        A = A.union(Ar)
        parent_t, dist_t = nx.dijkstra_predecessor_and_distance(Hr, s, weight = 'weight')
        spcomp += 1
        parent, dist, alpha, active_nodes = parent_dist_alpha_computation(Hr, delta_val, epsilon, r, parent_t, dist_t, S)
        while (alpha <= delta_val * math.pow((1+epsilon), r):
              if len(active_nodes) == 0:
                 print("No more active nodes")
                 break
              for u in active_nodes:
                  path = []
                  v = u
                  minc = math.inf
                  while (v != "-1"):
                        path.append(v)
                        v = parent[v]
                  #print(path, u)
                  for v in path:
                      if v != "-1" and v not in S:
                         if c[v] < minc:
                            minc = c[v]
                  d = 1 
                  while (True):
                        sum_path = 0
                        for v in path:
                            if v not in S:
                               sum_path += ell[v] * math.pow(1+(epsilon* minc)/c[v], d)
                        if sum_path > delta_val * math.pow((1+epsilon), r):
                           break
                        d += 1
                  for v in path:
                      if v not in S:
                         ell[v] = ell[v] * math.pow(1+(epsilon* minc)/c[v], d)
              #update weights on edges before next shortest path search
              for edge in Hr.edges():
                  u = edge[0] 
                  v = edge[1]
                  Hr[u][v]["weight"] = ell[v]
              parent_t, dist_t = nx.dijkstra_predecessor_and_distance(Hr, s, weight = 'weight')
              spcomp += 1
              parent, dist, alpha, active_nodes = parent_dist_alpha_computation(Hr, delta_val, epsilon, r, parent_t, dist_t, S)
        r = r+1
        print(r, max_itr, flush = True)
    noA = len(A)
    return ell, spcomp
        
def parent_dist_alpha_computation(Hr, delta_val, epsilon, r, parent_t, dist_t, S):
    parent = {}
    dist = {}
    alpha = math.inf
    active_nodes = set()   
    for u in Hr.nodes():
        if u == "-1":
           continue
        parent[u] = parent_t[u][0]
        dist[u] = dist_t[u]
        if u in S or "," not in u:
           continue
        if dist[u] < delta_val * math.pow((1+epsilon), r):
           active_nodes.add(u)
        if dist[u] < alpha:
           alpha = dist[u]
    return parent, dist, alpha, active_nodes


#generate a directed graph Hr from H   
def generateDiGraph(H, M, S, s, ell, c, delta_val):
    Hr = nx.DiGraph()
    #add bi-directional edges only when one of the nodes is in S. add uni-directional edge from s in S to u otherwise.
    for edge in H.edges():
        u = edge[0]
        v = edge[1]
        if v not in S:
           Hr.add_edge(u, v, weight = ell[v])
        if u not in S:
           Hr.add_edge(v, u, weight = ell[u])
    #add an edge from the super source s to each source in S
    for u in S:
        Hr.add_edge(s, u, weight = ell[u])
    Ar = set()
    #attach random copy a(v,j) to each v in H:
    for v in H.nodes():
        if v in S:
           continue
        j = random.randint(0, M-1)
        cv = v+","+str(j) #copy of v
        Ar.add(cv)
        c[cv] = 1.0/M
        if cv not in ell:
           ell[cv] = delta_val
        Hr.add_edge(v, cv, weight = ell[cv])
    return Hr, Ar     

#generate a random sample
def generateSamples(G, S, p):
    Hd = []
    for i in range(0, 10):
        H = nx.Graph()
        Hd.append(H)      
    for edge in G.edges():
        for i in range(0, 10):
            r = random.random()
            if r <= p:
               u = edge[0]
               v = edge[1]
               Hd[i].add_edge(u,v)
    #connect the super source s to every source node in Gp with an edge 
    s = "-1"
    for i in range(0, 10):
        Hd[i].add_node(s) 
    for node in S:
        for i in range(0, 10):
            Hd[i].add_edge(s, node)  
 
    #compute the BFS tree of Gp with s as the source
    Td = []
    for i in range(0, 10):
        Td.append(nx.bfs_tree(Hd[i], s))
        Td[i].remove_node(s)

    #compute number of infections in SIR process (percolation)
    infections = []
    infected_nodes = []
    for i in range(0, 10):
        infections.append(0)
        infected_nodes.append(set())
    for i in range(0, 10):
        for node in Td[i].nodes():
            infections[i] += 1
            infected_nodes[i].add(node)

    #H is the sampled graph induced by nodes in T, with vertex set V(H) equal to all nodes nodes reachable 
    #from sources and the edgeset E(H) is all edges in Gp whose both endpoints are in V(H)
    H = []
    for i in range(0, 10):
        H.append(nx.Graph())
    for i in range(0, 10):
        for edge in Hd[i].edges():
            v1 = edge[0]
            v2 = edge[1]
            if v1 in infected_nodes[i] and v2 in infected_nodes[i]:
               H[i].add_edge(v1, v2)
 
    #print("H_j^': sampled graph")
    #print("#infections ", len(H.nodes()))
    #print("#edges ", len(H.edges()))
    return H
            
#generate a random sample
def generateSample(G, S, p):
    Gp = nx.Graph()
    for edge in G.edges():
        r = random.random()
        if r <= p:
           u = edge[0]
           v = edge[1]
           Gp.add_edge(u,v)
    #connect the super source s to every source node in Gp with an edge
    s = "-1"
    Gp.add_node(s)

    for node in S:
        Gp.add_edge(s, node)

    #compute the BFS tree of Gp with s as the source
    T = nx.bfs_tree(Gp, s)
    T.remove_node(s)

    #compute number of infections in SIR process (percolation)
    infections = 0
    infected_nodes = set()
    for node in T.nodes():
        infections += 1
        infected_nodes.add(node)

    #H is the sampled graph induced by nodes in T, with vertex set V(H) equal to all nodes nodes reachable
    #from sources and the edgeset E(H) is all edges in Gp whose both endpoints are in V(H)
    H = nx.Graph()
    for edge in Gp.edges():
        v1 = edge[0]
        v2 = edge[1]
        if v1 in infected_nodes and v2 in infected_nodes:
           H.add_edge(v1, v2)

    #print("H_j^': sampled graph")
    #print("#infections ", len(H.nodes()))
    #print("#edges ", len(H.edges()))
    return H

