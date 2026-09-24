# verif_sage.sage
# A executer en local : sage verif_sage.sage    (ou magma verif.magma, section separee plus bas)
#
# But : produire une sortie brute verifiable pour etre citee dans l'article
# (Section 4), a la place de l'affirmation retiree faute de log retenu.
#
# Copiez-collez l'INTEGRALITE de la sortie de ce script dans la conversation --
# pas un resume -- pour qu'elle puisse etre citee avec le vrai log en annexe.

import time

print("=== A) Calcul direct (aucun raccourci), det M_p(sigma) == 0 ? p<150 ===")
def det_direct(p):
    n = (p - 1) // 2
    Fp = GF(p)
    Rxy.<x, s> = PolynomialRing(Fp, 2)
    f = x**7 + s*x + 1
    fn = f**n
    def coeff_ij(i, j):
        e = p*i - j
        return fn.coefficient({x: e})
    Sp = PolynomialRing(Fp, 's')
    M = matrix(Sp, 3, 3)
    for i in range(1, 4):
        for j in range(1, 4):
            M[i-1, j-1] = Sp(coeff_ij(i, j))
    return M.determinant()

t0 = time.time()
hits = []
for p in primes(3, 150):
    D = det_direct(p)
    if D == 0:
        hits.append(p)
print("  primes p<150 avec det identiquement nul :", hits)
print("  attendu : [3, 7, 11, 23]  -> ", "OK" if hits == [3,7,11,23] else "DESACCORD")
print("  temps :", round(time.time()-t0, 1), "s")

print("\n=== B) Table des R_a (recalculee independamment en Sage, arithmetique QQ exacte) ===")
def ratio(beta, j):
    alpha = QQ(-(j + beta)) / QQ(7)
    gamma = QQ(-1)/QQ(2) - alpha - beta
    if beta >= 1:
        return QQ(beta) / (gamma + 1)
    r = alpha
    for t in range(5):
        r *= (gamma - t)
    return r / 720

def Ra_value(a):
    beta = {(i, j): (a*i - j) % 7 for i in (1,2,3) for j in (1,2,3)}
    rho = {}
    for i in (1,2,3):
        rho[(i,1)] = QQ(1)
        for j in (1,2):
            rho[(i, j+1)] = rho[(i,j)] * ratio(beta[(i,j)], j)
    perms = Permutations([1,2,3])
    sums = {tuple(pm): sum(beta[(i, tuple(pm)[i-1])] for i in (1,2,3)) for pm in perms}
    mn = min(sums.values())
    T = [pm for pm in perms if sums[tuple(pm)] == mn]
    tot = QQ(0)
    for pm in T:
        pm = tuple(pm)
        inv_count = sum(1 for x in range(3) for y in range(x+1,3) if pm[x] > pm[y])
        sgn = -1 if inv_count % 2 else 1
        t = QQ(sgn)
        for i in (1,2,3):
            t *= rho[(i, pm[i-1])]
        tot += t
    return mn, T, tot

Ra_table = {}
for a in range(1, 7):
    mn, T, val = Ra_value(a)
    Ra_table[a] = val
    print(f"  a={a}: min_beta={mn}, |T|={len(T)}, R_a = {val}")

print("\n=== C) Formule (prod m_i1)*R_a vs extraction directe multinomiale, 101<=p<PMAX_C ===")
PMAX_C = 20000
def entries_low(p):
    n = (p - 1)//2
    fact = [1]*(n+1)
    for k in range(1, n+1):
        fact[k] = fact[k-1]*k % p
    inv = [pow(f, -1, p) for f in fact]
    low = {}
    for i in (1,2,3):
        for j in (1,2,3):
            N = p*i - j
            al = N // 7
            be = N - 7*al
            ga = n - al - be
            if ga < 0 or al > n:
                return None
            low[(i,j)] = (be, fact[n]*inv[al] % p * inv[be] % p * inv[ga] % p)
    return low

def fr_mod(x, p):
    return int(x.numerator()) % p * pow(int(x.denominator()) % p, -1, p) % p

tested = 0
bad = []
for p in primes(101, PMAX_C):
    a = p % 7
    low = entries_low(p)
    if low is None:
        continue
    best = None
    Tp = []
    for pm in Permutations([1,2,3]):
        pm = tuple(pm)
        s = sum(low[(i, pm[i-1])][0] for i in (1,2,3))
        if best is None or s < best:
            best = s; Tp = [pm]
        elif s == best:
            Tp.append(pm)
    tot = 0
    for pm in Tp:
        inv_count = sum(1 for x in range(3) for y in range(x+1,3) if pm[x] > pm[y])
        sgn = -1 if inv_count % 2 else 1
        t = sgn
        for i in (1,2,3):
            t = t * low[(i, pm[i-1])][1] % p
        tot = (tot + t) % p
    prod = 1
    for i in (1,2,3):
        prod = prod * low[(i,1)][1] % p
    tested += 1
    if tot != prod * fr_mod(Ra_table[a], p) % p:
        bad.append(p)
print(f"  premiers testes: {tested} | desaccords: {bad}")

print("\nTermine.")
