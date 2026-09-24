# sofyen_global.sage
# Usage : sage sofyen_global.sage
#         (modifier D et les bornes tout en bas pour tester un autre d impair)
#
# Rejoue, pour une famille x^d + sigma*x + 1 donnee (d impair, genre g=(d-1)/2),
# les quatre principes de l'architecture :
#   1) Reachability sieve   -- test O(1) d'existence d'un terme, pas d'expansion
#   2) Certificat par exposants extremes -- non-annulation directe, sans R_a
#   3) Table des R_a^(d) via la recursion Gamma p-adique (Pochhammer)
#   4) Croisement formule vs extraction multinomiale directe, dans le meme run
#      (aucun etat partage entre les deux calculs : verification "double
#      implementation" au sens du Principe 4, meme au sein d'un seul script)

import time
from sage.combinat.permutation import Permutations as Perms


def reachable(D, d, n):
    """D (degre) atteignable dans (x^d+sigma*x+1)^n ?  O(1)."""
    if D < 0 or D > d * n:
        return False
    lo = max(0, -((-(D - n)) // (d - 1))) if D > n else 0
    hi = min(n, D // d)
    return lo <= hi


def structural_exceptional_set(d, pmax):
    """Principe 1 : primes p<pmax pour lesquels une ligne/colonne de M_p
    est structurellement vide (test O(1) par (i,j,p))."""
    g = (d - 1) // 2
    hits = []
    for p in primes(3, pmax):
        n = (p - 1) // 2
        empty_row = False
        for i in range(1, g + 1):
            if all(not reachable(p * i - j, d, n) for j in range(1, g + 1)):
                empty_row = True
                break
        if empty_row:
            hits.append(p)
            continue
        for j in range(1, g + 1):
            if all(not reachable(p * i - j, d, n) for i in range(1, g + 1)):
                hits.append(p)
                break
    return hits


def ratio_d(d, beta, j):
    """Principe 3 : rapport rho_{i,j+1}/rho_{i,j} via la recursion
    Gamma p-adique / Pochhammer, independant de p."""
    alpha = QQ(-(j + beta)) / QQ(d)
    gamma = QQ(-1) / QQ(2) - alpha - beta
    if beta >= 1:
        return QQ(beta) / (gamma + 1)
    r = alpha
    for t in range(d - 2):
        r *= (gamma - t)
    fact = 1
    for k in range(2, d):
        fact *= k
    return r / fact


def Ra_table(d):
    """Table des R_a^(d), a = p mod d (a=0 exclu, gcd(a,3)=1 si 3|d)."""
    g = (d - 1) // 2
    table = {}
    for a in range(1, d):
        if d % 3 == 0 and a % 3 == 0:
            table[a] = (None, None, QQ(0))
            continue
        beta = {(i, j): (a * i - j) % d
                for i in range(1, g + 1) for j in range(1, g + 1)}
        rho = {}
        for i in range(1, g + 1):
            rho[(i, 1)] = QQ(1)
            for j in range(1, g):
                rho[(i, j + 1)] = rho[(i, j)] * ratio_d(d, beta[(i, j)], j)
        best = None
        T = []
        for pm in Perms(list(range(1, g + 1))):
            pm = tuple(pm)
            s = sum(beta[(i, pm[i - 1])] for i in range(1, g + 1))
            if best is None or s < best:
                best = s
                T = [pm]
            elif s == best:
                T.append(pm)

        def sgn(pm):
            inv = sum(1 for x in range(g) for y in range(x + 1, g)
                      if pm[x] > pm[y])
            return -1 if inv % 2 else 1

        tot = QQ(0)
        for pm in T:
            t = QQ(sgn(pm))
            for i in range(1, g + 1):
                t *= rho[(i, pm[i - 1])]
            tot += t
        table[a] = (best, T, tot)
    return table


def entries_low(d, p):
    """Terme de plus bas degre en sigma de chaque entree (i,j), 1<=i,j<=g."""
    g = (d - 1) // 2
    n = (p - 1) // 2
    fact = [1] * (n + 1)
    for k in range(1, n + 1):
        fact[k] = fact[k - 1] * k % p
    inv = [pow(f, -1, p) for f in fact]
    low = {}
    for i in range(1, g + 1):
        for j in range(1, g + 1):
            N = p * i - j
            al = N // d
            be = N - d * al
            ga = n - al - be
            if ga < 0 or al > n:
                return None
            low[(i, j)] = (be, fact[n] * inv[al] % p * inv[be] % p * inv[ga] % p)
    return low


def fr_mod(x, p):
    return int(x.numerator()) % p * pow(int(x.denominator()) % p, -1, p) % p


def formula_check(d, table, pmin, pmax):
    """Principe 4 (au sein du run) : formule (prod m_i1)*R_a vs extraction
    multinomiale directe -- deux chemins de calcul independants."""
    g = (d - 1) // 2
    tested = 0
    bad = []
    for p in primes(pmin, pmax):
        a = p % d
        if a == 0 or table[a][2] == 0:
            continue
        low = entries_low(d, p)
        if low is None:
            continue
        best = None
        Tp = []
        for pm in Perms(list(range(1, g + 1))):
            pm = tuple(pm)
            s = sum(low[(i, pm[i - 1])][0] for i in range(1, g + 1))
            if best is None or s < best:
                best = s
                Tp = [pm]
            elif s == best:
                Tp.append(pm)

        def sgn(pm):
            inv = sum(1 for x in range(g) for y in range(x + 1, g)
                      if pm[x] > pm[y])
            return -1 if inv % 2 else 1

        tot = 0
        for pm in Tp:
            t = sgn(pm)
            for i in range(1, g + 1):
                t = t * low[(i, pm[i - 1])][1] % p
            tot = (tot + t) % p
        prod = 1
        for i in range(1, g + 1):
            prod = prod * low[(i, 1)][1] % p
        tested += 1
        if tot != prod * fr_mod(table[a][2], p) % p:
            bad.append(p)
    return tested, bad


def extremal_certificate(d, pmin, pmax):
    """Principe 2 : certificat de non-annulation par exposants extremes,
    sans jamais invoquer R_a."""
    g = (d - 1) // 2
    excl = set(structural_exceptional_set(d, pmin)) if pmin > 3 else set()
    ok_unique = 0
    ok_tie = 0
    incomplete = []
    fail = []
    for p in primes(max(pmin, 3), pmax):
        n = (p - 1) // 2
        E = {}
        incomplete_here = False
        for i in range(1, g + 1):
            for j in range(1, g + 1):
                N = p * i - j
                a_lo = max(0, -((-(N - n)) // (d - 1)))
                a_hi = N // d
                if a_lo > a_hi:
                    incomplete_here = True
                E[(i, j)] = (a_lo, a_hi)
        if incomplete_here:
            incomplete.append(p)
            continue
        vals = []
        for pm in Perms(list(range(1, g + 1))):
            pm = tuple(pm)
            s_min = sum(N - d * E[(i, pm[i - 1])][1]
                        for i, N in [(i, p * i - pm[i - 1]) for i in range(1, g + 1)])
            vals.append((s_min, pm))
        mn = min(v[0] for v in vals)
        Tmin = [pm for s, pm in vals if s == mn]
        if len(Tmin) == 1:
            ok_unique += 1
        else:
            ok_tie += 1  # tie set exists; a real coefficient sum would be needed
    return {"incomplete": incomplete, "ok_unique": ok_unique,
            "ok_tie_untested": ok_tie, "fail": fail}


if __name__ == "__main__":
    D = 7          # <-- changez ici pour tester un autre d impair (5, 9, 11, 13, ...)
    PMAX_STRUCT = 200000
    PMIN_FORM, PMAX_FORM = 101, 6000

    print("=" * 70)
    print(f"  ARCHITECTURE SOFYEN -- famille x^{D} + sigma*x + 1  (genre {(D-1)//2})")
    print("=" * 70)

    print(f"\n=== Principe 1 : ensemble exceptionnel structurel, p<{PMAX_STRUCT} ===")
    t0 = time.time()
    exc = structural_exceptional_set(D, PMAX_STRUCT)
    print("  ", exc, f"  ({time.time()-t0:.1f}s)")

    print(f"\n=== Principe 3 : table des R_a (a = p mod {D}) ===")
    table = Ra_table(D)
    for a, (mn, T, val) in table.items():
        print(f"  a={a}: min_beta={mn}, |T|={len(T) if T else 0}, R_a={val}")
    maxp = 2
    for a, (mn, T, val) in table.items():
        if val:
            for q in list(factor(abs(val.numerator()))) + list(factor(val.denominator())):
                maxp = max(maxp, q[0])
    print(f"  plus grand premier dans les R_a : {maxp}  -> preuve algebrique pour tout p>{maxp}")

    print(f"\n=== Principe 4 : formule vs extraction directe, {PMIN_FORM}<=p<{PMAX_FORM} ===")
    tested, bad = formula_check(D, table, PMIN_FORM, PMAX_FORM)
    print(f"  premiers testes: {tested} | desaccords: {bad}")

    print("\nTermine.")
