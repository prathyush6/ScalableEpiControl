import sys
import networkx as nx
import math
import Multiplicativeweights
import Simulation

def algo(G, S, p, B, epsilon, M, path, network, lambda_val):
    print("====Input details====", flush = True)
    n = len(G.nodes())
    m = len(G.edges())
    print("|V|, |E| of network are "+str(n)+" and "+str(m), flush = True)
    print("Set of sources: "+str(S), flush = True)
    print("Budget B: "+str(B), flush = True)
    print("error parameter \epsilon: "+str(epsilon), flush = True)
    print("output path: "+str(path), flush = True)
    print("=====================", flush = True)
    print("====Parameter values====", flush = True)
    #set delta value
    delta_val = (1+epsilon)/math.pow((1+epsilon)*6, 1.0/epsilon)
    print("delta: "+str(delta_val), flush = True)
    print("lambda: "+str(lambda_val), flush = True)
    total_spcomp = 0
    ell, spcomp = Multiplicativeweights.MWU(G, S, epsilon, delta_val, lambda_val, M, path)
    print("Done with given lambda", flush = True)
    total_spcomp += spcomp
    budget_used, avginfections = compute_budget(ell, M)
    print("lambda, budget_used, avg. infections: ", lambda_val, budget_used, avginfections, flush = True)
    outstring = ""
    outstring += str(M)+","+str(B)+","+str(budget_used)+","+str(avginfections)
    print("Total shortest path computations: "+str(total_spcomp), flush = True) 
    fw.close()

def generateGraph(network):
    G = nx.Graph()
    fp = open(network, 'r')
    lines = fp.readlines()
    for line in lines:
        cols = line.split()
        u = cols[0]
        v = cols[1]
        G.add_edge(u, v)
    return G

def generateSources(sourcefile):
    S = set()
    fp = open(sourcefile, 'r')
    lines = fp.readlines()
    for line in lines:
        u = line.strip()
        S.add(u)
    return S

def compute_budget(ell, M):
    budget_used = 0
    no_infections = 0
    for u in ell.keys():
        if "," not in u:
            budget_used += ell[u]
        else:
            no_infections += ell[u]
    avginfections = no_infections/M + len(S)
    return budget_used, avginfections

if __name__ == '__main__':
   network = sys.argv[1] #contact network file name
   sourcefile = sys.argv[2] #source nodes in a file; one source per line
   B = int(sys.argv[3]) #number of vaccinations available
   p = float(sys.argv[4]) #transmission probability p
   epsilon = float(sys.argv[5]) #error parameter epsilon value
   M = int(sys.argv[6]) #number of sampled graphs
   foldername = sys.argv[7] #name of folder containing sample dags
   lambda_init = float(sys.argv[8]) #initial value of lambda
   print("Input Details", flush = True)
   print("network, sourcefile, B, epsilon, M, foldername ", network, sourcefile, B, epsilon, M, foldername, flush = True)
   #Generate a graph G and a set of source nodes S from the input file
   G = generateGraph(network)
   S = generateSources(sourcefile)   
   algo(G, S, p, B, epsilon, M, foldername, network, lambda_init)

