from sage.all import DirichletGroup
from dirichlet_conrey import DirichletGroup_conrey
#import pymongo
# from lmfdb.website import dbport
from sage.all import ZZ, QQ
from lmfdb.modular_forms.elliptic_modular_forms.backend.web_modforms import WebModFormSpace

def generate_eta_quotient_data(maxN=100, maxk=12):
    data = dict()
    for N in range(1, maxN):
        data[N] = dict()
        for k in range(2, maxk):
            for chiorbit in [chi for chi in range ]
            W = WebModFormSpace(N=int(N), k=int(k), chi=int(chi))
            finds = ms.find({'t': [int(N), int(k), int(0)]})
            print finds.count()
            in_db = finds.count() > 0
            data[N][k] = {'dimension': dim, 'in_db': in_db}
            # print N,k,data[N][k]
            
    
    return ms, data

def eta_quotient_with_conrey_number(A, N, arguments, exponents):
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
        c2 = kronecker_character(2**(2k) * A)
        chi = D(c1)*D(c2)
    else:
        chi = kronecker_character(A)

    conrey_n = Dc.from_sage_character(chi)
    eta_q = eta_quotient(arguments=arguments, exponents=exponents, N=N, ret=1)
    return eta_q, conrey_n
    
    
def eta_quotient(arguments=None,exponents=None,N=100,ret=0):
    if arguments == None:
        arguments = []
    if exponents == None:
        arguments = []
    assert len(arguments)==len(exponents)
    eta =  qexp_eta(ZZ[['q']], N)
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