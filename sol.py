from Crypto.Util.number import *
import itertools
from tqdm import tqdm

def generate_basis(n):
    basis = [True] * n
    for i in range(3, int(n**0.5)+1, 2):
        if basis[i]:
            basis[i*i::2*i] = [False]*((n-i*i-1)//(2*i)+1)
    return [2] + [i for i in range(3, n, 2) if basis[i]]


def miller_rabin(n, b):
    """
    Miller Rabin test testing over all
    prime basis < b
    """
    basis = generate_basis(b)
    if n == 2 or n == 3:
        return True

    if n % 2 == 0:
        return False

    r, s = 0, n - 1
    while s % 2 == 0:
        r += 1
        s //= 2
    for b in basis:
        x = pow(b, s, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True



def xgcd(a, b):
    s = 0
    t = 1
    r = b
    s1 = 1
    t1 = 0
    r1 = a
    while not (r == 0):
        q = r1//r
        r1, r = r, r1-q*r
        s1, s = s, s1-q*s
        t1, t = t, t1-q*t
        #print(r1, s1, t1)
    return (r1, s1, t1)

def crt1(residues, modulos):
    #print(modulos)
    rm = list(zip(residues, modulos))
    cur_res, cur_mod = rm[0]
    for r, m in rm[1:]:
        
        g = GCD(cur_mod, m)
        
        if not r % g == cur_res % g:
            return -1, -1
        r1, s, t = xgcd(m//g, cur_mod//g)
       # print(r, cur_res, r % g, cur_res%g, s, t)
        cur_res = cur_res * m//g * s + r * cur_mod//g * t
        cur_mod *= m//g
        cur_res %= cur_mod
    return cur_res, cur_mod

primes = generate_basis(64)
print(len(primes))
fool = []
h = 3
def legendre(a, p):
    # if a == 2:
        # return (-1)**((p**2-1)//8)
    return pow(a, (p-1)//2, p)

for p in primes:
    f = set()
    for i in generate_basis(200*p)[1:]:
        #print(legendre(p, i))
        if legendre(p, i) == i-1:
            f.add(i % (4*p))
            
    fool.append(list(f))
    
print(primes)
#print(fool)
ks = [1, 998244353, 233]   
fool2 = []
for p, f in enumerate(fool):
    prime = primes[p]
    m = prime*4
    cur_set = set(f)
    for i in range (1, h):
        new_set = set()
        for ff in f:
            if ((ff + ks[i] - 1)*inverse(ks[i], m)) % 4 == 3:
                new_set.add(((ff + ks[i] - 1)*inverse(ks[i], m)) % m)
        cur_set = cur_set.intersection(new_set)
    fool2.append(cur_set)
print(fool2)
mm = 1
for a in fool2:
    mm *= len(a)

print(mm)
pr = 0
for tup in itertools.product(*fool2):
    residues = []
    modulos = []
    #modul = 1
    for i, t in enumerate(tup):
        residues.append(t)
        modulos.append(primes[i]*4)
        # g = GCD(primes[i]*4, modul)
        # modul *= (primes[i]*4)//g
    residues.append(ks[1]-inverse(ks[2], ks[1]))
    modulos.append(ks[1])
    residues.append(ks[2]-inverse(ks[1], ks[2]))    
    modulos.append(ks[2])
    #print(modulos)

    sol, modul = crt1(residues, modulos)
    #rint(sol)
    #print(modul)
    #print(sol)
    found = False
    if not sol == -1:
        cur_t = sol
        #print(sol % modul)
        cur_t = 2**73*modul + cur_t
        for i in tqdm(range(100000)):
            if isPrime(cur_t):
                fin = cur_t
                facs = [cur_t]
                for ii in range(1, h):
                    facs.append(ks[ii]*(cur_t-1) + 1)
                    fin *= ks[ii]*(cur_t-1) + 1
                #print(fin.bit_length())
                if(miller_rabin(fin, 64)):
                    print(isPrime(fin))
                    print(fin)
                    print(facs)
                    if fin.bit_length() >= 600 and fin.bit_length() <= 900:
                        found = True
                        break
            cur_t += modul

    if found:
        break
