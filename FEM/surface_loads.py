import numpy as np
from utils import kos


# --- 1. Допоміжні функції для NUSHAB ---

def _first_valid_nf(nf, nodes):
    """Шукає перше значення NF > 10 серед перелічених вузлів."""
    for node in nodes:
        if node > 0 and nf[node - 1] > 10:
            return nf[node - 1]
    return 0


def _calc_j_flag(i1, i2, i3, i4):
    """Обчислює прапорець J (J1, J2, J3) за станом граничних умов."""
    if i2 > 10 and i3 > 10:
        return 1
    if i2 < 10 and i3 > 10 and i4 > 10:
        return 2
    if i3 < 10 and i1 > 10 and i2 > 10:
        return 3
    if i1 < 10 and i2 > 10 and i3 < 10:
        return 4
    if i2 < 10 and i3 > 10 and i4 < 10:
        return 5
    if i2 < 10 and i3 < 10:
        return 6
    return 0


def nushab(nf, n, state):
    """Формування шаблону оточення вузла N та визначення прапорців J1, J2, J3."""
    m1, m2, m3 = state.m1, state.m2, state.m3

    # 1. Координати вузла N у 3D-сітці
    ks1 = kos(n, 1, state)
    ks2 = kos(n, 2, state)
    ks3 = kos(n, 3, state)

    # 2. Обчислення осьових суміжних вузлів N19..N34
    n19 = n - 2 if ks1 > 2 else 0
    n20 = n - 1 if ks1 > 1 else 0
    n21 = n + 1 if ks1 + 1 <= m1 else 0
    n22 = n + 2 if ks1 + 2 <= m1 else 0

    n14 = n - 2 * m1 if ks2 > 2 else 0
    n17 = n - m1 if ks2 > 1 else 0
    n24 = n + m1 if ks2 + 1 <= m2 else 0
    n25 = n + 2 * m1 if ks2 + 2 <= m2 else 0

    n31 = n - 2 * m1 * m2 if ks3 > 2 else 0
    n9  = n - m1 * m2 if ks3 > 1 else 0
    n33 = n + m1 * m2 if ks3 + 1 <= m3 else 0
    n34 = n + 2 * m1 * m2 if ks3 + 2 <= m3 else 0

    # 3. Обчислення кутових та діагональних вузлів
    n1 = n2 = n3 = n4 = n5 = n6 = n7 = n8 = n10 = n11 = n12 = 0
    if ks3 > 1:
        n1  = n - 1 - 2 * m1 - m1 * m2 if (ks1 > 1 and ks2 > 2) else 0
        n2  = n - 2 * m1 - m1 * m2 if ks2 > 2 else 0
        n3  = n - 2 - m1 - m1 * m2 if (ks1 > 2 and ks2 > 1) else 0
        n4  = n - 1 - m1 - m1 * m2 if (ks1 > 1 and ks2 > 1) else 0
        n5  = n - m1 - m1 * m2 if ks2 > 1 else 0
        n6  = n + 1 - m1 - m1 * m2 if (ks1 + 1 <= m1 and ks2 > 1) else 0
        n7  = n - 2 - m1 * m2 if ks1 > 2 else 0
        n8  = n - 1 - m1 * m2 if ks1 > 1 else 0
        n10 = n + 1 - m1 * m2 if ks1 + 1 <= m1 else 0
        n11 = n - 1 + m1 - m1 * m2 if (ks1 > 1 and ks2 + 1 <= m2) else 0
        n12 = n + m1 - m1 * m2 if ks2 + 1 <= m2 else 0

    n13 = n - 1 - 2 * m1 if (ks1 > 1 and ks2 > 2) else 0
    n15 = n - 2 - m1 if (ks1 > 2 and ks2 > 1) else 0
    n16 = n - 1 - m1 if (ks1 > 1 and ks2 > 1) else 0
    n18 = n + 1 - m1 if (ks1 + 1 <= m1 and ks2 > 1) else 0
    n23 = n - 1 + m1 if (ks1 > 1 and ks2 + 1 <= m2) else 0
    n30 = n - m1 - 2 * m1 * m2 if (ks3 > 2 and ks2 > 1) else 0
    n32 = n - m1 + m1 * m2 if (ks3 + 1 <= m3 and ks2 > 1) else 0

    n26 = n27 = n28 = n29 = 0
    if ks1 > 1:
        n26 = n - 1 - m1 - 2 * m1 * m2 if (ks3 > 2 and ks2 > 1) else 0
        n27 = n - 1 - 2 * m1 * m2 if ks3 > 2 else 0
        n28 = n - 1 - m1 + m1 * m2 if (ks3 + 1 <= m3 and ks2 > 1) else 0
        n29 = n - 1 + m1 * m2 if ks3 + 1 <= m3 else 0

    # 4. Обчислення прапорців IA, IB, IT
    ia1 = _first_valid_nf(nf, [n3, n7, n15, n19])
    ia2 = _first_valid_nf(nf, [n4, n8, n16, n20])
    ia3 = _first_valid_nf(nf, [n5, n9, n17, n])
    ia4 = _first_valid_nf(nf, [n6, n10, n18, n21])

    ib1 = _first_valid_nf(nf, [n1, n2, n13, n14])
    ib2 = _first_valid_nf(nf, [n4, n5, n16, n17])
    ib3 = _first_valid_nf(nf, [n8, n9, n20, n])
    ib4 = _first_valid_nf(nf, [n11, n12, n23, n24])

    it1 = _first_valid_nf(nf, [n26, n27, n30, n31])
    it2 = _first_valid_nf(nf, [n4, n8, n5, n9])
    it3 = _first_valid_nf(nf, [n16, n20, n17, n])
    it4 = _first_valid_nf(nf, [n28, n29, n32, n33])

    # 5. Розрахунок J1, J2, J3
    state.j1 = _calc_j_flag(ia1, ia2, ia3, ia4)
    state.j2 = _calc_j_flag(ib1, ib2, ib3, ib4)
    state.j3 = _calc_j_flag(it1, it2, it3, it4)

    # 6. Збереження шаблону в стан
    state.n19, state.n20, state.n21, state.n22 = n19, n20, n21, n22
    state.n14, state.n17, state.n24, state.n25 = n14, n17, n24, n25
    state.n31, state.n9,  state.n33, state.n34 = n31, n9, n33, n34


# --- 2. Геометричні та матричні підпрограми ---

def r_dist(n1, n2, X):
    """Евклідова відстань між вузлами n1 та n2 (1-based indices)."""
    if n1 <= 0 or n2 <= 0:
        return 0.0
    return float(np.linalg.norm(X[n1 - 1] - X[n2 - 1]))


def prs(j, a, b, f1, fo, f2):
    """Апроксимація похідної за 3-точковим або 2-точковим шаблоном."""
    if j == 1:
        return (a**2 * f2 + (b**2 - a**2) * fo - b**2 * f1) / (a * b * (a + b))
    elif j == 2:
        return (2.0 * a * (-a * f2 + (a + b) * fo - b * f1) + (a**2 * f2 + (b**2 - a**2) * fo - b**2 * f1)) / (a * b * (a + b))
    elif j == 3:
        return (2.0 * b * (a * f2 - (a + b) * fo + b * f1) + (a**2 * f2 + (b**2 - a**2) * fo - b**2 * f1)) / (a * b * (a + b))
    elif j == 4:
        return (fo - f1) / a if a != 0 else 0.0
    elif j == 5:
        return (f2 - fo) / b if b != 0 else 0.0
    return 0.0


def obrm33(A, state=None):
    """Обернення матриці 3x3 та збереження детермінанта в state."""
    det = float(np.linalg.det(A))
    if state is not None:
        state.det = det
    inv_A = np.linalg.inv(A)
    return inv_A, det


# --- 3. Тензори перетворення координат ---

def tipck(n, x, state):
    """Обчислення матриці TIPCK (3x3)."""
    CK = np.zeros((3, 3), dtype=float)
    stencils = [
        (state.j1, state.n19, state.n20, state.n21, state.n22),
        (state.j2, state.n14, state.n17, state.n24, state.n25),
        (state.j3, state.n31, state.n9,  state.n33, state.n34)
    ]

    for col, (j, n_p2, n_p1, n_n1, n_n2) in enumerate(stencils):
        if j == 6:
            continue

        a, b = 0.0, 0.0
        idx_f1, idx_fo, idx_f2 = -1, n - 1, -1

        if j == 1:
            a, b = r_dist(n_p1, n, x), r_dist(n, n_n1, x)
            idx_f1, idx_f2 = n_p1 - 1, n_n1 - 1
        elif j == 2:
            a, b = r_dist(n, n_n1, x), r_dist(n_n1, n_n2, x)
            idx_f1, idx_fo, idx_f2 = n - 1, n_n1 - 1, n_n2 - 1
        elif j == 3:
            a, b = r_dist(n_p2, n_p1, x), r_dist(n_p1, n, x)
            idx_f1, idx_fo, idx_f2 = n_p2 - 1, n_p1 - 1, n - 1
        elif j == 4:
            a = r_dist(n_p1, n, x)
            idx_f1 = n_p1 - 1
        elif j == 5:
            b = r_dist(n, n_n1, x)
            idx_f2 = n_n1 - 1

        for k in range(3):
            f1 = x[idx_f1, k] if idx_f1 >= 0 and j in (1, 2, 3, 4) else 0.0
            fo = x[idx_fo, k]
            f2 = x[idx_f2, k] if idx_f2 >= 0 and j in (1, 2, 3, 5) else 0.0
            CK[k, col] = prs(j, a, b, f1, fo, f2)

    return CK


def ckcn(x, n, state):
    """Обчислення тензорів CK та CN."""
    CK = tipck(n, x, state)
    CN, _ = obrm33(CK, state)
    return CK, CN


def gkn(CK, state=None):
    """Обчислює метричний тензор GK = CK^T * CK та його обернений GN."""
    GK = CK.T @ CK
    GN, _ = obrm33(GK, state)
    return GK, GN


# --- 4. Головний розрахунок поверхневого потоку ---

def _process_node_flux(n, kg, pp, nux, nf, x, qr, state, povna3_func):
    """Обробка одного вузла всередині HAPOS3."""
    nushab(nf, n, state)

    CK, _ = ckcn(x, n, state)
    GK, _ = gkn(CK, state)

    A = np.zeros(3, dtype=float)
    for i in range(3):
        A[i] = povna3_func(n, kg, x, nux, i + 1, state) * pp / np.sqrt(GK[i, i])

    qr[n - 1] += CK @ A


def hapos3(n1, n2, n3, n4, kg, nux, nf, X, QR, state, povna3_func):
    """Основна функція обчислення поверхневого потоку на елементі."""
    r12, r24, r14 = r_dist(n1, n2, X), r_dist(n2, n4, X), r_dist(n1, n4, X)
    r13, r34 = r_dist(n1, n3, X), r_dist(n3, n4, X)

    p1 = 0.5 * (r12 + r24 + r14)
    area1 = np.sqrt(abs(p1 * (p1 - r12) * (p1 - r24) * (p1 - r14)))

    p2 = 0.5 * (r13 + r14 + r34)
    area2 = np.sqrt(abs(p2 * (p2 - r14) * (p2 - r13) * (p2 - r34)))

    pp = (area1 + area2) * 0.25

    for node in (n1, n2, n3, n4):
        _process_node_flux(node, kg, pp, nux, nf, X, QR, state, povna3_func)