import numpy as np
from utils import (
    coor, wpro, ckpro, ckpr2s, ckpr3s, obrm33, tran33, gkn,
    defcuo, defc2s, defc3s, fun, fun2, fun3, alfat, sigmaterm,
    prfmk, prfm2s, prfm3s, defor1, defcue, defcu2, defcu3,
    sigdek, rmnli, kos, stx
)

def grprck(nu, ck, state):
    """
    Підпрограма GRPRCK[cite: 12].
    Формування матриці перетворення координат для вузла.
    """
    m3 = state.m3
    ks3 = kos(nu, 3, state)
    state.sver = False

    if ks3 == m3:
        state.sver = True
        alf = np.arctan(1.0)
        ca = np.cos(alf)
        sa = np.sin(alf)

        ck[0, 0] = ca
        ck[1, 0] = sa
        ck[2, 0] = 0.0

        ck[0, 1] = -sa
        ck[1, 1] = ca
        ck[2, 1] = 0.0

        ck[0, 2] = 0.0
        ck[1, 2] = 0.0
        ck[2, 2] = 1.0


def qtep3(v, mnl, state):
    """
    Підпрограма QTEP3[cite: 11].
    Перетворення компонентів у глобальну систему координат.
    """
    c = np.zeros((3, 3), dtype=float)

    for i in range(1, mnl + 1):
        nu = state.nn[i - 1]
        grprck(nu, c, state)

        if not state.sver:
            continue

        n1 = i - 1
        n2 = i - 1 + mnl
        n3 = i - 1 + 2 * mnl

        w1 = v[n1] * c[0, 0] + v[n2] * c[1, 0] + v[n3] * c[2, 0]
        w2 = v[n1] * c[0, 1] + v[n2] * c[1, 1] + v[n3] * c[2, 1]
        w3 = v[n1] * c[0, 2] + v[n2] * c[1, 2] + v[n3] * c[2, 2]

        v[n1] = w1
        v[n2] = w2
        v[n3] = w3


def brafm3(mn, mnl, al, w, xe, vt, bt, v, mnl3, state):
    """
    Обчислення теплового напруження.
    """
    kt1 = state.kt1
    kt2 = state.kt2
    kt3 = state.kt3
    mske = state.mske
    mcx = state.mcx
    corect = state.corect
    det = state.det

    # Ініціалізація локальних масивів та нумерації
    et = np.zeros(9, dtype=float)
    st = np.zeros(9, dtype=float)
    bt_arr = np.array(bt, dtype=float)

    gk = np.zeros((3, 3), dtype=float)
    gn = np.zeros((3, 3), dtype=float)
    ck = np.zeros((3, 3), dtype=float)
    cn = np.zeros((3, 3), dtype=float)
    cnt = np.zeros((3, 3), dtype=float)
    fi = np.zeros((3, 3), dtype=float)
    f = np.zeros(6, dtype=float)
    fii = np.zeros(9, dtype=float)
    stt = np.zeros(6, dtype=float)
    gkdop = np.zeros((3, 3), dtype=float)
    pdop = np.zeros((3, 3), dtype=float)
    p = np.zeros(3, dtype=float)

    # Очищення векторів V
    for i in range(mnl3):
        v[i] = 0.0

    coor(xe, al, 0.0, 0.0, 0.0)
    wpro(w[0:], mnl, xe[0:], al)
    wpro(w[mnl:], mnl, xe[mnl:], al)
    wpro(w[2 * mnl:], mnl, xe[2 * mnl:], al)
    wpro(vt, mnl, bt_arr, al)

    # Цикли інтегрування за Гауссом
    for i in range(1, kt1 + 1):
        o1 = state.csi[0, i - 1]
        q1 = state.as_arr[0, i - 1]
        for j in range(1, kt2 + 1):
            o2 = state.csi[1, j - 1]
            q2 = state.as_arr[1, j - 1]
            for k in range(1, kt3 + 1):
                o3 = state.csi[2, k - 1]
                q3 = state.as_arr[2, k - 1]

                if mcx < 11:
                    ckpro(ck, mn, mnl, w, o1, o2, o3)
                elif mcx == 11:
                    ckpr2s(ck, mn, mnl, w, o1, o2, o3)
                elif mcx == 12:
                    ckpr3s(ck, mn, mnl, w, o1, o2, o3)

                obrm33(ck, cn)
                tran33(cn, cnt)
                gkn(gk, gn, ck)

                if mske == 2:
                    if mcx < 11:
                        defcuo(gkdop, w, w, mnl, o1, o2, o3)
                    elif mcx == 11:
                        defc2s(gkdop, w, w, mnl, o1, o2, o3)
                    elif mcx == 12:
                        defc3s(gkdop, w, w, mnl, o1, o2, o3)

                    gkdop[1, 0] = gkdop[0, 1]
                    gkdop[2, 0] = gkdop[0, 2]
                    gkdop[2, 1] = gkdop[1, 2]

                    for ick in range(3):
                        for jck in range(3):
                            gk[ick, jck] = gkdop[ick, jck] + (gk[ick, jck] - gkdop[ick, jck]) * corect
                    obrm33(gk, gn)

                aq = q1 * q2 * q3 * np.sqrt(det)

                if mcx < 11:
                    t_val = fun(vt, o1, o2, o3)
                elif mcx == 11:
                    t_val = fun2(vt, o1, o2, o3)
                elif mcx == 12:
                    t_val = fun3(vt, o1, o2, o3)

                at = alfat(t_val) * t_val
                et[0] = at
                et[4] = at
                et[8] = at

                sigmaterm(et, st, ck, t_val)
                stt[0] = st[0]
                stt[1] = st[3]
                stt[2] = st[6]
                stt[3] = st[4]
                stt[4] = st[7]
                stt[5] = st[8]

                for ii in range(1, mnl + 1):
                    itek = mnl * (ii - 1)
                    if mcx < 11:
                        prfmk(p, mn, al[itek], o1, o2, o3)
                    elif mcx == 11:
                        prfm2s(p, mn, al[itek], o1, o2, o3)
                    elif mcx == 12:
                        prfm3s(p, mn, al[itek], o1, o2, o3)

                    for jj in range(1, 4):
                        jrr = mnl * (jj - 1)
                        defor1(fi, ck, p, jj)

                        if mske == 2:
                            if mcx < 11:
                                defcue(pdop, w[jrr:], al[itek], o1, o2, o3)
                            elif mcx == 11:
                                defcu2(pdop, w[jrr:], al[itek], o1, o2, o3)
                            elif mcx == 12:
                                defcu3(pdop, w[jrr:], al[itek], o1, o2, o3)

                            for ick in range(3):
                                for jck in range(3):
                                    fi[ick, jck] = pdop[ick, jck] + (fi[ick, jck] - pdop[ick, jck]) * corect

                        sigdek(fi, fii, cnt)
                        tran33(fii, fi)
                        f[0] = fi[0, 0]
                        f[1] = fi[0, 1]
                        f[2] = fi[0, 2]
                        f[3] = fi[1, 1]
                        f[4] = fi[1, 2]
                        f[5] = fi[2, 2]

                        agp = rmnli(f, stt) * aq
                        ij_idx = mnl * (jj - 1) + (ii - 1)
                        v[ij_idx] += agp


def pntpro(numb, i, state):
    """
    Формування глобальних номерів вузлів для поточного скінченного елемента.
    """
    m1 = state.m1
    m2 = state.m2

    m = state.mm
    n = state.nn
    l = state.ll

    kn = 0
    for k in range(1, l + 1):
        for j in range(1, n + 1):
            for ii in range(1, m + 1):
                kn += 1
                # 0-based індексація для Python
                # Оригінальна формула Fortran: NUMB + II - 1 + M1*(J - 1) + M1*M2*(K - 1)
                state.nn[kn - 1] = numb + (ii - 1) + m1 * (j - 1) + m1 * m2 * (k - 1)
                

def obratm(a, n, l, m):
    """
    Інвертування матриці методом Гаусса-Жордана з вибором головного елемента.
    """
    d = 1.0

    # Головний цикл по K від 0 до N-1 (0-based)
    for k in range(n):
        kk = k * n + k
        biga = a[kk]
        l_idx = k
        m_idx = k

        # Пошук максимального елемента в підматриці
        for j in range(k, n):
            iz = j * n
            for i in range(k, n):
                ij = iz + i
                if abs(biga) < abs(a[ij]):
                    biga = a[ij]
                    l_idx = i
                    m_idx = j

        l[k] = l_idx
        m[k] = m_idx

        # Перестановка рядків
        j_val = l[k]
        if j_val > k:
            for i in range(n):
                ki = i * n + k
                ji = i * n + j_val
                hold = -a[ki]
                a[ki] = a[ji]
                a[ji] = hold

        # Перестановка стовпців
        i_val = m[k]
        if i_val > k:
            jp = i_val * n
            for j in range(n):
                jk = k * n + j
                ji = jp + j
                hold = -a[jk]
                a[jk] = a[ji]
                a[ji] = hold

        if biga == 0.0:
            d = 0.0
            return d

        # Ділення елементів рядка/стовпця
        for i in range(n):
            if i != k:
                ik = k * n + i
                a[ik] = a[ik] / (-biga)

        for i in range(n):
            ij = i * n
            hold = a[k * n + i] if i != k else 0.0 # тимчасове збереження
            # перевірки індексів для модифікації залишку матриці
            for j in range(n):
                if i != k and j != k:
                    # точне відображення формули KJ = IJ - I + K з Fortran
                    # в 0-based: ij = i*n + j
                    pass

        # Повне перетворення матриці (жорданові виключення)
        for i in range(n):
            ik = k * n + i
            hold = a[ik]
            for j in range(n):
                ij = i * n + j
                if i != k and j != k:
                    kj = ij - i + k
                    a[ij] = hold * a[kj] + a[ij]

        for j in range(n):
            kj = k * n + j
            if j != k:
                a[kj] = a[kj] / biga

        d = d * biga
        a[kk] = 1.0 / biga

    # Зворотні перестановки (кроки 100-150 у Fortran)
    k = n
    while k > 0:
        k -= 1
        i_val = l[k]
        if i_val > k:
            jq = k * n
            jr = i_val * n
            for j in range(n):
                jk = jq + j
                hold = a[jk]
                ji = jr + j
                a[jk] = -a[ji]
                a[ji] = hold

        j_val = m[k]
        if j_val > k:
            for i in range(n):
                ki = i * n + k
                ji = i * n + j_val
                hold = a[ki]
                a[ki] = -a[ji]
                a[ji] = hold

    return d


def alpro(al, jn, la, ma, state):
    """
    Формування матриці апроксимації та її звернення.
    """
    mm = state.mm
    nn = state.nn
    ll = state.ll

    mmh = mm - 1
    nnh = nn - 1
    llh = ll - 1
    mn = mm * nn

    # Цикл по LL
    for k in range(1, ll + 1):
        a3 = float(k)
        if llh == 1:
            o3 = 2.0 * a3 - 3.0
        elif llh == 2:
            o3 = a3 - 2.0
        elif llh == 3:
            o3 = (2.0 * a3 - 5.0) / 3.0
        else:
            o3 = a3

        for j in range(1, nn + 1):
            a2 = float(j)
            if nnh == 1:
                o2 = 2.0 * a2 - 3.0
            elif nnh == 2:
                o2 = a2 - 2.0
            elif nnh == 3:
                o2 = (2.0 * a2 - 5.0) / 3.0
            else:
                o2 = a2

            for i in range(1, mm + 1):
                a1 = float(i)
                if mmh == 1:
                    o1 = 2.0 * a1 - 3.0
                elif mmh == 2:
                    o1 = a1 - 2.0
                elif mmh == 3:
                    o1 = (2.0 * a1 - 5.0) / 3.0
                else:
                    o1 = a1

                n_idx = i + mm * (j - 1) + mn * (k - 1)

                for kk in range(1, ll + 1):
                    lt3 = kk - 1
                    for jj in range(1, nn + 1):
                        lt2 = jj - 1
                        for ii in range(1, mm + 1):
                            lt1 = ii - 1
                            m_idx = ii + mm * (jj - 1) + mn * (kk - 1)

                            x1 = stx(lt1, o1)
                            x2 = stx(lt2, o2)
                            x3 = stx(lt3, o3)

                            # 0-based індексація для матриці AL (вхідна розмірність JN x JN)
                            al[(n_idx - 1), (m_idx - 1)] = x1 * x2 * x3

    obratm(al, jn, 0.0, la, ma)


def fotepq(nf, x, t, q, nux, al, mnl3, w, xe, vt, bt, v, mn, mnl, state):
    """
    Формування правих частин та зв'язків для кінцевих елементів.
    """
    m1 = state.m1
    m2 = state.m2
    m3 = state.m3

    i1 = state.i1
    i2 = state.i2
    i3 = state.i3

    i1d = i1 - 1
    i2d = i2 - 1
    i3d = i3 - 1

    if i1d != 0: i1b = (m1 - 1) // i1d
    else: i1b = 1

    if i2d != 0: i2b = (m2 - 1) // i2d
    else: i2b = 1

    if i3d != 0: i3b = (m3 - 1) // i3d
    else: i3b = 1

    for ii3c in range(1, i3b + 1):
        iks3 = 1 + i3d * (ii3c - 1)
        for ii2c in range(1, i2b + 1):
            iks2 = 1 + i2d * (ii2c - 1)
            for ii1c in range(1, i1b + 1):
                iks1 = 1 + i1d * (ii1c - 1)

                j = iks1 + m1 * (iks2 - 1) + m1 * m2 * (iks3 - 1)
                inkp = nf[j - 1]

                if inkp < 10:
                    continue

                pntpro(j, 1)

                for i1c in range(1, mnl + 1):
                    kks1 = state.nn[i1c - 1]
                    ks2 = i1c + mnl
                    ks3 = i1c + 2 * mnl

                    bt[i1c - 1] = t[kks1 - 1]
                    xe[i1c - 1] = x[kks1 - 1, 0]
                    xe[ks2 - 1] = x[kks1 - 1, 1]
                    xe[ks3 - 1] = x[kks1 - 1, 2]

                brafm3(mn, mnl, al, w, xe, vt, bt, v, mnl3)

                if state.utem and state.preo:
                    qtep3(v, mnl)

                for i in range(1, mnl + 1):
                    kks1 = state.nn[i - 1]
                    ks2 = i + mnl
                    ks3 = ks2 + mnl

                    q[kks1 - 1, 0] += v[i - 1]
                    q[kks1 - 1, 1] += v[ks2 - 1]
                    q[kks1 - 1, 2] += v[ks3 - 1]


# ==========================================
# РЕАЛІЗОВАНА ПІДПРОГРАМА NMAQ1
# ==========================================
def nmaq1(nf, ng, x, t, q, nux, state):
    """
    Ініціалізує локальні масиви для кінцевого елемента та викликає подальші процедури
    """
    m = state.ldk[0]
    n = state.ldk[1]
    l = state.ldk[2]

    al = np.zeros(64, dtype=float)
    w = np.zeros(24, dtype=float)
    xe = np.zeros(24, dtype=float)
    vt = np.zeros(8, dtype=float)
    bt = np.zeros(8, dtype=float)
    v = np.zeros(24, dtype=float)
    la = np.zeros(8, dtype=int)
    ma = np.zeros(8, dtype=int)

    mn = m * n
    mnl = mn * l
    mnl3 = mnl * 3

    # Передаємо state у підпрограми для доступу до спільних змінних
    alpro(al, mnl, la, ma, state)
    fotepq(nf, x, t, q, nux, al, mnl3, w, xe, vt, bt, v, mn, mnl, state)