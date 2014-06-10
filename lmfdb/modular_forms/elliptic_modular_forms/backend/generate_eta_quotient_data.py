from sage.all import DirichletGroup
from dirichlet_conrey import DirichletGroup_conrey
#import pymongo
# from lmfdb.website import dbport
from sage.all import ZZ, QQ
#from lmfdb.modular_forms.elliptic_modular_forms.backend.web_modforms import WebModFormSpace
import itertools

def generate_eta_quotient_data(maxN=100, maxk=12, maxr=None, minN=1):
# TODO: say which condition means what!
    data = dict()
    for N in range(minN, maxN+1):
        data[N] = dict()
        deltas = tuple(sorted(Integer(N).divisors()))
        delta_n = len(deltas)
        if maxr == None:
            maxrr = max(24 // delta_n, 1) # This allows for higher levels (N) without taking too much time
        else:
            maxrr = maxr
        for rdeltas in itertools.product(range(-1*maxrr,maxrr+1), repeat=delta_n): #This runs through the exponents r_d of the eta(d tau) for d a divisor of N uo to a given bound
            #print "step 1:", rdeltas
            if sum([deltas[j]*rdeltas[j] for j in range(0,delta_n)]) == 24: #This is to have the resulting eta quotient equal to q+O(q^2) 
                k = Integer(sum(rdeltas)) #This is double the weight of the eta quotient
                #print "step 2:", rdeltas, k
                if k % 2 == 0 and k <= 2*maxk and 4 <= k: #Make shure, we have an integer weight within given bounds
                    #print "step 3:", rdeltas
                    if sum([rdeltas[j]*deltas[j] for j in range(0,delta_n)]) % 24 == 0: #Condition from Borcherds Refl... Thm 6.2
                        #print "step 4:", rdeltas
                        if N * sum([Integer(rdeltas[j])/Integer(deltas[j]) for j in range(0,delta_n)]) % 24 == 0: # As above
                            #print "step 5:", rdeltas
                            k = k // 2
                            conrey_n, conrey_c, dirich_n, dirich_c = eta_quotient_conrey_number(N, k, deltas, rdeltas) #Get the chararcter and their numbers
                            chi = dirich_n #This depends on the definition of chi in the WebModFormSpace
                            #W = WebModFormSpace(N=int(N), k=int(k), chi=int(chi)) #Get the WebModFormSpace 
                            #print dirich_c
                            print N, rdeltas, k, conrey_n, dirich_n
                            NF = Newforms(dirich_c, k, base_ring=None, names='x')
                            #print NF
                            eta_q = eta_quotient(arguments=deltas, exponents=rdeltas, prec=10, ret=1)
                            #print eta_q
                            for nf in NF:
                                #print nf, nf.base_ring()
                                if nf.base_ring() == QQ:
                                    nf_q = nf.q_expansion(prec=10)
                                    if eta_q == nf_q:
                                        # Here we should compere the coefficients up to the Sturm bound
                                        #prec = W.sturm_bound()
                                        if k in data[N].keys():
                                            if conrey_n in data[N][k].keys():
                                                #data[N][k][conrey_n].append([deltas, rdeltas, eta_q])
                                                data[N][k][conrey_n].append((deltas, rdeltas))
                                            else:
                                                #data[N][k][conrey_n] = [[deltas, rdeltas, eta_q]]
                                                data[N][k][conrey_n] = [(deltas, rdeltas)]
                                        else:
                                            data[N][k] = dict()
                                            #data[N][k][conrey_n] = [[deltas, rdeltas, eta_q]]
                                            data[N][k][conrey_n] = [(deltas, rdeltas)]
                                            #print eta_q, chi
    return data

def eta_quotient_conrey_number(N, k, arguments, exponents):
    r"""
    Under the conditions given in Borcherds Reflection groups
    of Lorentzian lattices, 
    this should compute the eta-quotient and the conrey number
    of a dirichlet character such that the eta-quotient
    is a modular form for Gamma0(N) and the character.
    """
    A = prod([arguments[j]**(exponents[j] % 2) for j in range(0,len(arguments))])
    D = DirichletGroup(N)
    Dc = DirichletGroup_conrey(N)
    if N % 4 == 0:
        c1 = kronecker_character(-1)**(k+((kronecker(-1,A)-1)//2))
        c2 = kronecker_character_upside_down(2**(2*k) * A)
        # print "c1", c1
        # print "c2", c2
        # print "D", D
        dirich_c = D(c1)*D(c2)
    else:
        c3 = kronecker_character_upside_down(A)
        # print "c3", c3
        # print "D", D
        dirich_c = D(c3)

    #print "N", N
    #print "D", D
    #print "chi:", chi
    #print "Dc", Dc
    conrey_c = Dc.from_sage_character(dirich_c)
    conrey_n = conrey_c.number()
    dirich_n = D.list().index(dirich_c)
    #print "conrey_n", conrey_n
    return conrey_n, conrey_c, dirich_n, dirich_c
    
    
def eta_quotient(arguments=None,exponents=None,prec=100,ret=0):
    if arguments == None:
        arguments = []
    if exponents == None:
        arguments = []
    assert len(arguments)==len(exponents)
    eta =  qexp_eta(ZZ[['q']], prec)
    R = eta.parent()
    q = R.gens()[0]
    res = R(1)
    prefak = 0
    for i in range(len(arguments)):
        res = res*eta.subs({q:q**arguments[i]})**exponents[i] + O(q**prec)
        prefak = prefak+arguments[i]*exponents[i]
    r = prefak % 24
    if ret == 1 and r == 0:
        return res*q**(prefak/24)+O(q**prec)
    else:
        return res + O(q**prec),prefak/24

def eta_quotient_newforms():
    data = {
        1: {12: {0: [((1,), (24,))]}},
        2: {8: {1: [((1, 2), (8, 8))]}},
        3: {6: {1: [((1, 3), (6, 6))]}},
        4: {
            5: {3: [((1, 2, 4), (4, 2, 4))]},
            6: {1: [((1, 2, 4), (0, 12, 0))]}
        },
        5: {4: {1: [((1, 5), (4, 4))]}},
        6: {4: {1: [((1, 2, 3, 6), (2, 2, 2, 2))]}},
        7: {3: {6: [((1, 7), (3, 3))]}},
        8: {
            3: {3: [((1, 2, 4, 8), (2, 1, 1, 2))]},
            4: {1: [((1, 2, 4, 8), (0, 4, 4, 0))]}
        },
        9: {4: {1: [((1, 3, 9), (0, 8, 0))]}},
        10: {},
        11: {2: {1: [((1, 11), (2, 2))]}},
        12: {3: {5: [((1, 2, 3, 4, 6, 12), (0, 3, 0, 0, 3, 0))]}},
        13: {},
        14: {2: {1: [((1, 2, 7, 14), (1, 1, 1, 1))]}},
        15: {2: {1: [((1, 3, 5, 15), (1, 1, 1, 1))]}},
        16: {
            3: {15: [((1, 2, 4, 8, 16), (0, 0, 6, 0, 0))]},
            4: {1: [((1, 2, 4, 8, 16), (0, -4, 16, -4, 0))]},
            6: {1: [((1, 2, 4, 8, 16), (0, -12, 36, -12, 0))]} # Entered by hand
        },
        17: {},
        18: {},
        19: {},
        20: {2: {1: [((1, 2, 4, 5, 10, 20), (0, 2, 0, 0, 2, 0))]}},
        21: {},
        22: {},
        23: {},
        24: {2: {1: [((1, 2, 3, 4, 6, 8, 12, 24), (0, 1, 0, 1, 1, 0, 1, 0))]}},
        25: {},
        26: {},
        27: {2: {1: [((1, 3, 9, 27), (0, 2, 2, 0))]}},
        28: {},
        29: {},
        30: {},
        31: {},
        32: {
            2: {1: [((1, 2, 4, 8, 16, 32), (0, 0, 2, 2, 0, 0))]},
            3: {15: [((1, 2, 4, 8, 16, 32), (0, -2, 5, 5, -2, 0))]}
        },
        33: {},
        34: {},
        35: {},
        36: {2: {1: [((1, 2, 3, 4, 6, 9, 12, 18, 36), (0, 0, 0, 0, 4, 0, 0, 0, 0))]}},
        37: {},
        38: {},
        39: {}
    }
    return data

def check_eta_quotient_newforms():
    data = eta_quotient_newforms()
    for N in data.keys():
        for k in data[N].keys():
            for conrey_n in data[N][k].keys():
                for eta_qs in data[N][k][conrey_n]:
                    arguments = eta_qs[0]
                    exponents = eta_qs[1]
                    Dc = DirichletGroup_conrey(N)
                    conrey_nn, conrey_c, dirich_n, dirich_c = eta_quotient_conrey_number(N, k, arguments, exponents)
                    orbit_conrey_n = [Dc.from_sage_character(d).number() for d in dirich_c.galois_orbit()]
                    NF = Newforms(dirich_c, k, base_ring=None, names='x')
                    b = False
                    eta_q = eta_quotient(arguments, exponents, prec = 10, ret = 1)
                    for nf in NF:
                        #print nf, nf.base_ring()
                        if nf.base_ring() == QQ:
                            nf_q = nf.q_expansion(prec=10)
                            print "\n", nf_q, "\n", eta_q, "\n"
                            if eta_q == nf_q:
                                b = True
                    print N, k, conrey_n, orbit_conrey_n, arguments, exponents, conrey_n == conrey_nn, b#, eta_q
                    