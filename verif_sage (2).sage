"""
Verification independante (Sage) du "Sofyen Theorem" (ex-Theoreme C) :
  y^2 = x^7 + sigma*x + 1, genre 3, Cartier-Manin M_p(sigma).
  det M_p(sigma) = 0 identiquement dans F_p[sigma]  <=>  p in {3,7,11,23}.

Reimplementation COMPLETEMENT INDEPENDANTE du code Python original (aucune
ligne recopiee) :

  A) calcul DIRECT de (x^7+sigma*x+1)^n mod p (n=(p-1)/2), sans AUCUN
     raccourci combinatoire, puis determinant symbolique en sigma.
     -> confirme Regime I / le pont fini / le brute-force p<100.
     Cout : croit vite avec p, reserve a p<~150.

  B) reimplementation independante de la formule du coefficient extremal
     (prod m_i1)*R_a, en arithmetique rationnelle exacte Sage (QQ).

  C) comparaison de B) contre une extraction VRAIMENT directe (multinomiale,
     sans passer par la formule Ra) -> teste Regime II a grande echelle
     sans jamais faire confiance a la formule elle-meme.

Usage :
  sage verif_sage.sage
(ou, dans un notebook / la console sage : load("verif_sage.sage"))
"""

# ============================================================
# A) Calcul direct, aucun raccourci : det M_p(sigma) == 0 ?
# ============================================================

def cartier_manin_direct(p):
    """Matrice 3x3 M(sigma) sur F_p[sigma], obtenue en developpant
    f(x)^n = (x^7+sigma*x+1)^n TEL QUEL (aucune formule de coefficient),
    n=(p-1)//2, g=3. Renvoie une matrice a coefficients dans F_p[sigma]."""
    Fp = GF(p)
    S = PolynomialRing(Fp, 's'); s = S.gen()
    R = PolynomialRing(S, 'x'); x = R.gen()
    n = (p - 1) // 2
    f = x**7 + s*x + 1
    fn = f**n                      # puissance directe (carre-et-multiplie de Sage)
    M = matrix(S, 3, 3)
    for i in range(1, 4):
        for j in range(1, 4):
            N = p*i - j
            M[i-1, j-1] = fn[N] if 0 <= N <= fn.degree() else S(0)
    return M

def det_identiquement_nul(p):
    M = cartier_manin_direct(p)
    return M.determinant() == 0

print("=== A) Calcul direct (aucun raccourci) : det M_p(sigma) == 0 ? (p<PMAX_A) ===")
PMAX_A = 100   # augmenter prudemment : le cout croit vite (n=(p-1)/2 termes)
exceptionnels_direct = []
for p in prime_range(3, PMAX_A):
    if det_identiquement_nul(p):
        exceptionnels_direct.append(p)
print("  primes p<%d avec det identiquement nul :" % PMAX_A, exceptionnels_direct)
attendu = [q for q in [3, 7, 11, 23] if q < PMAX_A]
print("  attendu :", attendu, " -> ", "OK" if exceptionnels_direct == attendu else "ECART !")


# ============================================================
# B) Formule du coefficient extremal, reimplementee independamment
# ============================================================

def Ra_value(a):
    """R_a (rationnel exact), calcule via la definition mathematique
    directe -- reecrit ici depuis zero, en Sage/QQ, pas copie du Python."""
    def ratio(beta, j):
        alpha = QQ(-(j + beta)) / 7
        gamma = QQ(-1)/QQ(2) - alpha - beta
        if beta >= 1:
            return QQ(beta) / (gamma + 1)
        r = alpha
        for t in range(5):
            r *= (gamma - t)
        return r / 720
    beta = {(i, j): (a*i - j) % 7 for i in (1, 2, 3) for j in (1, 2, 3)}
    rho = {}
    for i in (1, 2, 3):
        rho[(i, 1)] = QQ(1)
        for j in (1, 2):
            rho[(i, j+1)] = rho[(i, j)] * ratio(beta[(i, j)], j)
    perms = Permutations([1, 2, 3])
    sums = {tuple(pm): sum(beta[(i, tuple(pm)[i-1])] for i in (1, 2, 3)) for pm in perms}
    mn = min(sums.values())
    T = [pm for pm in perms if sums[tuple(pm)] == mn]
    total = QQ(0)
    for pm in T:
        t = pm.sign()
        for i in (1, 2, 3):
            t *= rho[(i, tuple(pm)[i-1])]
        total += t
    return mn, T, total

print("\n=== B) Table des R_a (recalculee independamment en Sage) ===")
Ra_table = {}
for a in range(1, 7):
    mn, T, val = Ra_value(a)
    Ra_table[a] = val
    print(f"  a={a}: min_beta={mn}, |T|={len(T)}, R_a = {val}")


def m_ij(p, i, j):
    """Coefficient de plus bas degre en sigma de l'entree (i,j), mod p."""
    n = (p - 1) // 2
    N = p*i - j
    beta = N % 7
    alpha = (N - beta) // 7
    gamma = n - alpha - beta
    if gamma < 0 or alpha > n:
        return None
    val = (factorial(n) * inverse_mod(factorial(alpha), p)
           * inverse_mod(factorial(beta), p) * inverse_mod(factorial(gamma), p)) % p
    return beta, val

def coeff_via_formule(p):
    """(prod m_i1) * R_a mod p -- utilise la table Ra_table."""
    a = p % 7
    prod_m1 = 1
    for i in (1, 2, 3):
        r = m_ij(p, i, 1)
        if r is None:
            return None
        prod_m1 = (prod_m1 * r[1]) % p
    Ra = Ra_table[a]
    num = Ra.numerator() % p
    den = Ra.denominator() % p
    if den % p == 0:
        return None
    return (prod_m1 * num * inverse_mod(den, p)) % p

def coeff_direct_multinomial(p):
    """Coefficient de sigma^Bmin dans det M_p(sigma), calcule en extrayant
    VRAIMENT les termes multinomiaux extremaux (aucun recours a R_a) --
    verification independante de la formule elle-meme."""
    n = (p - 1) // 2
    entries = {}
    for i in (1, 2, 3):
        for j in (1, 2, 3):
            r = m_ij(p, i, j)
            if r is None:
                return None
            entries[(i, j)] = r  # (beta, valeur)
    perms = Permutations([1, 2, 3])
    sums = {tuple(pm): sum(entries[(i, tuple(pm)[i-1])][0] for i in (1, 2, 3)) for pm in perms}
    mn = min(sums.values())
    T = [pm for pm in perms if sums[tuple(pm)] == mn]
    total = 0
    for pm in T:
        t = 1 if pm.sign() == 1 else -1
        for i in (1, 2, 3):
            t = (t * entries[(i, tuple(pm)[i-1])][1]) % p
        total = (total + t) % p
    return total % p

print("\n=== C) Formule (B) vs extraction directe multinomiale, comparaison sur 101<=p<PMAX_C ===")
PMAX_C = 20000   # rapide (arithmetique modulaire uniquement) ; augmenter librement
tested = 0
mism = []
for p in prime_range(101, PMAX_C):
    lhs = coeff_direct_multinomial(p)
    rhs = coeff_via_formule(p)
    if lhs is None or rhs is None:
        continue
    tested += 1
    if lhs != rhs:
        mism.append(p)
print(f"  premiers testes: {tested} | desaccords: {mism}")
print("\nTermine.")
