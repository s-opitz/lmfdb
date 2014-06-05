from sage.all import DirichletGroup
from dirichlet_conrey import DirichletGroup_conrey
#import pymongo
# from lmfdb.website import dbport
from sage.all import ZZ, QQ
#from lmfdb.modular_forms.elliptic_modular_forms.backend.web_modforms import WebModFormSpace
import itertools

def generate_eta_quotient_data(maxN=100, maxk=12, maxr=40):

    data = dict()
    for N in range(1, maxN+1):
        data[N] = dict()
        deltas = tuple(sorted(Integer(N).divisors()))
        delta_n = len(deltas)
        for rdeltas in itertools.product(range(-1*maxr,maxr+1), repeat=delta_n):
            k = Integer(sum(rdeltas))
            if k % 2 == 0 and k <= 2*maxk and 0 < k:
                if sum([rdeltas[j]*deltas[j] for j in range(0,delta_n)]) % 24 == 0:
                    if N * sum([Integer(rdeltas[j])/Integer(deltas[j]) for j in range(0,delta_n)]) % 24 == 0:
                        k = k // 2
                        # TODO find A

                        pro = prod([deltas[j]**rdeltas[j] for j in range(0,delta_n)])
                        A1 = Integer(pro.denominator()*pro.numerator()).squarefree_part()
                        # print "pro", pro, "A1", A1
                        As = [delta for delta in deltas if (Integer(delta)/ A1).is_square()]
                        #print pro, A1, As, N
                        for A in As:
                            try:
                                eta_q,chi = eta_quotient_with_conrey_number(A, N, k, deltas, rdeltas)
                            except TypeError as e:
                                print type(e), e, "A", A, "N", N, "k", k
                            #W = WebModFormSpace(N=int(N), k=int(k), chi=int(chi))
                            else:
                                if k in data[N].keys():
                                    if (deltas, rdeltas) in data[N][k].keys():
                                        data[N][k][(deltas, rdeltas)].append((eta_q, chi, A))
                                    else:
                                        data[N][k][(deltas, rdeltas)] = [(eta_q, chi, A)]
                                else:
                                    data[N][k] = dict()
                                    data[N][k][(deltas, rdeltas)] = [(eta_q, chi, A)]
                                    #print eta_q, chi
    return data

def eta_quotient_with_conrey_number(A, N, k, arguments, exponents, prec=5):
    r"""
    Under the conditions given in Borcherds Reflection groups
    of Lorentzian lattices, 
    this should compute the eta-quotient and the conrey number
    of a dirichlet character such that the eta-quotient
    is a modular form for Gamma0(N) and the character.
    """

    D = DirichletGroup(N)
    Dc = DirichletGroup_conrey(N)
    if N % 4 == 0:
        c1 = kronecker_character(-1)**(k+((kronecker(-1,A)-1)//2))
        c2 = kronecker_character_upside_down(2**(2*k) * A)
        # print "c1", c1
        # print "c2", c2
        # print "D", D
        chi = D(c1)*D(c2)
    else:
        c3 = kronecker_character_upside_down(A)
        # print "c3", c3
        # print "D", D
        chi = D(c3)

    #print "N", N
    #print "D", D
    #print "chi:", chi
    #print "Dc", Dc
    conrey_n = Dc.from_sage_character(chi)
    #print "conrey_n", conrey_n
    eta_q = eta_quotient(arguments=arguments, exponents=exponents, prec=prec, ret=1)
    return eta_q, conrey_n
    
    
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
        res = res*eta.subs({q:q**arguments[i]})**exponents[i]
        prefak = prefak+arguments[i]*exponents[i]
    r = prefak % 24
    if ret == 1 and r == 0:
        return res*q**(prefak/24)
    else:
        return res,prefak/24