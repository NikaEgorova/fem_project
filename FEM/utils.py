import numpy as np
from surface_loads import hapos3


def telos(i1, i2, i3, k1, k2, k3, nf, state):
    """
    Формування ознак елементів/вузлів сітки в об'ємі.
    """
    m1, m2, m3 = state.m1, state.m2, state.m3

    for l1 in range(i1, k1 + 1):
        for l2 in range(i2, k2 + 1):
            for l3 in range(i3, k3 + 1):
               # Обчислення 0-based індексу в одновимірному масиві
                nu = (l1 - 1) + m1 * (l2 - 1) + m1 * m2 * (l3 - 1)

                if l1 < k1 and l2 < k2 and l3 < k3:
                    nf[nu] = 71
                else:
                    nf[nu] = 7


def zakrep(i1, i2, i3, k1, k2, k3, iz, nf, state):
    """
    Задавання граничних умов (закріплень) для вузлів сітки.
    """
    m1, m2, m3 = state.m1, state.m2, state.m3

    for l1 in range(i1, k1 + 1):
        for l2 in range(i2, k2 + 1):
            for l3 in range(i3, k3 + 1):
                # Обчислення 0-based індексу в одновимірному масиві
                nu = (l1 - 1) + m1 * (l2 - 1) + m1 * m2 * (l3 - 1)
                
                # Модифікація прапорця вузла на місці
                nf[nu] -= iz


def hagpo3(nx1, nx2, nx3, kx1, kx2, kx3, kg, nux, nf, x, qr, povna_func, state):
    """
    Обчислення та прикладення поверхневого навантаження на гранях елементів.
    """
    m1, m2, m3 = state.m1, state.m2, state.m3

    # Перехід GOTO(1,1,2,2,3,3), KG реалізовано через умовні оператори
    if kg in (1, 2):
        # Грані паралельні площині YZ (фіксований X)
        for j in range(nx2, kx2):
            for k in range(nx3, kx3):
                n1 = nx1 + m1 * (j - 1) + m1 * m2 * (k - 1)
                n2 = n1 + m1
                n3 = n1 + m1 * m2
                n4 = n3 + m1

                # Передаємо в HAPOS3 вузли у 0-based індексації (-1)
                hapos3(n1 - 1, n2 - 1, n3 - 1, n4 - 1, kg, nux, nf, x, qr, povna_func, state)

    elif kg in (3, 4):
        # Грані паралельні площині XZ (фіксований Y)
        for i in range(nx1, kx1):
            for k in range(nx3, kx3):
                n1 = i + m1 * (nx2 - 1) + m1 * m2 * (k - 1)
                n2 = n1 + 1
                n3 = n1 + m1 * m2
                n4 = n3 + 1

                hapos3(n1 - 1, n2 - 1, n3 - 1, n4 - 1, kg, nux, nf, x, qr, povna_func, state)

    elif kg in (5, 6):
        # Грані паралельні площині XY (фіксований Z)
        for i in range(nx1, kx1):
            for j in range(nx2, kx2):
                n1 = i + m1 * (j - 1) + m1 * m2 * (nx3 - 1)
                n2 = n1 + 1
                n3 = n1 + m1
                n4 = n3 + 1

                hapos3(n1 - 1, n2 - 1, n3 - 1, n4 - 1, kg, nux, nf, x, qr, povna_func, state)

def nud(ng, i, state=None):
    """
    Повертає номер рівняння з масиву NG за 1-based індексом i.
    """
    return int(ng[i - 1])

def pntpro(numb, i, state):
    """
    Заповнює масив state.nn індексами вузлів локального блоку.
    
    numb: 1-based індекс опорного вузла
    i: резервний параметр Fortran (збережено для сумісності сигнатури)
    """
    # Параметри блоку KSO (M, N, L) із збережених у darmir значень state.ldk
    m = state.ldk[0]  # M
    n = state.ldk[1]  # N
    l = state.ldk[2]  # L

    m1 = state.m1
    m2 = state.m2

    kn = 0
    for k in range(1, l + 1):
        for j in range(1, n + 1):
            for ii in range(1, m + 1):
                # Формування 1-based індексу вузла для зберігання в state.nn
                node_idx = numb + (ii - 1) + m1 * (j - 1) + m1 * m2 * (k - 1)
                state.nn[kn] = node_idx
                kn += 1

def nnnst3(nf, ng, l, state, data_txt=None):
    """
    Обчислення напівширини стрічки NSTR.
    """
    # Зчитуємо M, N, LL із збережених у darmir значень state.ldk
    m = state.ldk[0]  # 2
    n = state.ldk[1]  # 2
    ll_kso = state.ldk[2]  # 2

    mnl = m * n * ll_kso  # 2 * 2 * 2 = 8
    nnstr = -1000000

    i1d = m - 1
    i2d = n - 1
    i3d = ll_kso - 1

    i1b = (state.m1 - 1) // i1d if i1d != 0 else 1
    i2b = (state.m2 - 1) // i2d if i2d != 0 else 1
    i3b = (state.m3 - 1) // i3d if i3d != 0 else 1

    for ii3c in range(1, i3b + 1):
        iks3 = 1 + i3d * (ii3c - 1)
        for ii2c in range(1, i2b + 1):
            iks2 = 1 + i2d * (ii2c - 1)
            for ii1c in range(1, i1b + 1):
                iks1 = 1 + i1d * (ii1c - 1)

                # 1-based індекс опорного вузла J
                j = iks1 + state.m1 * (iks2 - 1) + state.m1 * state.m2 * (iks3 - 1)

                if nf[j - 1] < 10:
                    continue

                # Заповнення state.nn (вузли макроелемента)
                pntpro(j, 1, state)

                ngmin = 1000000
                ngmax = -1000000

                for i1c in range(mnl):
                    kks1 = state.nn[i1c]
                    if kks1 <= 0 or nf[kks1 - 1] < 0:
                        continue

                    ngg = nud(ng, kks1, state)

                    if ngg < ngmin:
                        ngmin = ngg
                    if ngg > ngmax:
                        ngmax = ngg

                maxmin = ngmax - ngmin
                if maxmin > nnstr:
                    nnstr = maxmin

    nstr = nnstr + l  # додаємо ll = 3

    if data_txt:
        data_txt.write(f"          NSTR={nstr:5d}\n")

    return nstr


def ipamat(nar, state):
    """
    Відновлена підпрограма IPAMAT.
    Ініціалізує номер матеріалу (MAT) у стані за кодом армування NAR.
    """
    state.nar = nar
    state.mat = nar if (nar is not None and nar > 0) else 1

    if not hasattr(state, "mre"):
        state.mre = 1
    if not hasattr(state, "ire"):
        state.ire = 0


def fgaus(kti):
    """
    Обчислює вузли (csi) та ваги (as) квадратур Гаусса 1..5 порядку.
    """
    csi = np.zeros(10, dtype=float)
    as_w = np.zeros(10, dtype=float)

    if kti == 1:
        csi[0] = 0.0
        as_w[0] = 2.0

    elif kti == 2:
        csi[0] = -np.sqrt(1.0 / 3.0)
        csi[1] = -csi[0]
        for i in range(2):
            pe = 3.0 * csi[i]
            as_w[i] = 2.0 / ((1.0 - csi[i] ** 2) * pe**2)

    elif kti == 3:
        csi[0] = 0.0
        csi[1] = np.sqrt(3.0 / 5.0)
        csi[2] = -csi[1]
        for i in range(3):
            pe = 0.5 * (15.0 * csi[i] ** 2 - 3.0)
            as_w[i] = 2.0 / ((1.0 - csi[i] ** 2) * pe**2)

    elif kti == 4:
        csi[0] = np.sqrt((15.0 - 2.0 * np.sqrt(30.0)) / 35.0)
        csi[1] = -csi[0]
        csi[2] = np.sqrt((15.0 + 2.0 * np.sqrt(30.0)) / 35.0)
        csi[3] = -csi[2]
        for i in range(4):
            pe = 0.5 * (35.0 * csi[i] ** 3 - 15.0 * csi[i])
            as_w[i] = 2.0 / ((1.0 - csi[i] ** 2) * pe**2)

    elif kti == 5:
        csi[0] = 0.0
        csi[1] = np.sqrt((35.0 + 2.0 * np.sqrt(70.0)) / 63.0)
        csi[2] = -csi[1]
        csi[3] = np.sqrt((35.0 - 2.0 * np.sqrt(70.0)) / 63.0)
        csi[4] = -csi[3]
        for i in range(5):
            pe = (
                315.0 * csi[i] ** 4 - 210.0 * csi[i] ** 2 + 15.0
            ) * 0.125
            as_w[i] = 2.0 / ((1.0 - csi[i] ** 2) * pe**2)

    return csi, as_w


def parke(m, n, l, state):
    """
    Налаштовує розмірності макроелемента, схему інтегрування MCX
    та заповнює масиви квадратур Гаусса по 3 осях.
    """
    state.ldk[0] = m
    state.ldk[1] = n
    state.ldk[2] = l
    state.ldk[3] = 1

    mnl = m * n * l
    state.kt1 = m
    state.kt2 = n
    state.kt3 = l

    # Визначення конфігурації MCX за таблицею MNL
    mcx_map = {
        12: 2,
        16: 3,
        18: 4,
        24: 5,
        27: 6,
        32: 7,
        36: 8,
        48: 9,
        64: 10,
    }
    state.mcx = mcx_map.get(mnl, 1)

    # Виклик відновленої IPAMAT
    nar = getattr(state, "nar", 1)
    ipamat(nar, state)

    # Ініціалізація масивів точок та ваг Гаусса (ACSI та AAS)
    if not hasattr(state, "acsi"):
        state.acsi = np.zeros((3, 10), dtype=float)
        state.aas = np.zeros((3, 10), dtype=float)

    kts = [state.kt1, state.kt2, state.kt3]
    for axis_idx, kt in enumerate(kts):
        csi, as_w = fgaus(kt)
        state.acsi[axis_idx, :kt] = csi[:kt]
        state.aas[axis_idx, :kt] = as_w[:kt]

def nmak1(nf, ng, x, t, q, nux, state):
    pass

def grprck(i, c):
    # Заглушка для функції отримання матриці перетворення
    pass


def sigke3(u, x, t, nf, i, nux, nms, d, sig, at, en, sn):
    # Заглушка для обчислення напружень
    pass


def coor(xe, al, v1, v2, v3):
    """Заглушка для COOR."""
    pass

def wpro(w_target, mnl, xe_source, al):
    """Заглушка для WPRO."""
    pass

def ckpro(ck, mn, mnl, w, o1, o2, o3):
    """Заглушка для CKPRO."""
    pass

def ckpr2s(ck, mn, mnl, w, o1, o2, o3):
    """Заглушка для CKPR2S."""
    pass

def ckpr3s(ck, mn, mnl, w, o1, o2, o3):
    """Заглушка для CKPR3S."""
    pass

def obrm33(a, b):
    """Заглушка для OBRM33."""
    pass

def tran33(a, b):
    """Заглушка для TRAN33."""
    pass



def defcuo(gkdop, w1, w2, mnl, o1, o2, o3):
    """Заглушка для DEFCUO."""
    pass

def defc2s(gkdop, w1, w2, mnl, o1, o2, o3):
    """Заглушка для DEFC2S."""
    pass

def defc3s(gkdop, w1, w2, mnl, o1, o2, o3):
    """Заглушка для DEFC3S."""
    pass

def fun(vt, o1, o2, o3):
    """Заглушка для FUN."""
    return 0.0

def fun2(vt, o1, o2, o3):
    """Заглушка для FUN2."""
    return 0.0

def fun3(vt, o1, o2, o3):
    """Заглушка для FUN3."""
    return 0.0

def alfat(t):
    """Заглушка для ALFAT."""
    return 0.0

def sigmaterm(et, st, ck, t):
    """Заглушка для SIGMATERM."""
    pass

def prfmk(p, mn, al_val, o1, o2, o3):
    """Заглушка для PRFMK."""
    pass

def prfm2s(p, mn, al_val, o1, o2, o3):
    """Заглушка для PRFM2S."""
    pass

def prfm3s(p, mn, al_val, o1, o2, o3):
    """Заглушка для PRFM3S."""
    pass

def defor1(fi, ck, p, jj):
    """Заглушка для DEFOR1."""
    pass

def defcue(pdop, w_val, al_val, o1, o2, o3):
    """Заглушка для DEFCUE."""
    pass

def defcu2(pdop, w_val, al_val, o1, o2, o3):
    """Заглушка для DEFCU2."""
    pass

def defcu3(pdop, w_val, al_val, o1, o2, o3):
    """Заглушка для DEFCU3."""
    pass

def sigdek(fi, fii, cnt):
    """Заглушка для SIGDEK."""
    pass

def rmnli(f, stt):
    """Заглушка для RMNLI."""
    return 0.0

def uvis3(nu, j, state):
    return 0.0

def stx(i, a):
    """
    Обчислення степеня параметра a залежно від індексу i.
    """
    if i <= 0:
        return 1.0
    return a ** i


def kos(n, j, state):
    """
    Визначення індексу вузла за заданим напрямком (1, 2 або 3) на основі розмірностей з state.
    """
    m1 = state.m1
    m2 = state.m2
    mn = m1 * m2

    if j == 1:
        i3 = (n - 1) // mn + 1
        m = n - (i3 - 1) * mn
        i2 = (m - 1) // m1 + 1
        k = m - (i2 - 1) * m1
    elif j == 2:
        i3 = (n - 1) // mn + 1
        m = n - (i3 - 1) * mn
        k = (m - 1) // m1 + 1
    elif j == 3:
        k = (n - 1) // mn + 1
    else:
        k = 0

    return k


def pecda3(x, q, t, nf, ng, nux, nms, state, data_txt):
    """
    Друк табличних даних (координати, навантаження, температура, ознаки).
    """
    # Відновлений та об'єднаний текстовий заголовок з масиву A1
    a1 = (
        '  ВУЗОЛ       КООРДИНАТИ                 '
        'НАВАНТАЖЕННЯ           ТЕМП-РА       ОЗНАКА '
    )

    # Текстовий заголовок з масиву A2
    a2 = (
        '    NU         Z1          Z2          Z3          '
        'Q1          Q2          Q3          TT          NF     NG   '
    )

    # FORMAT 1: Запис заголовків (відповідає FORMAT(/,32A4/))
    data_txt.write('\n' + a1 + '\n')
    data_txt.write('\n' + a2 + '\n')

    # FORMAT 2: 3X,I4, 3X,D12.5 (x7), 3X,I4, 3X,I4
    for i in range(nms):
        node_idx = i + 1

        line = (
            f"   {node_idx:4d}   "
            f"{x[i, 0]:12.5E}   {x[i, 1]:12.5E}   {x[i, 2]:12.5E}   "
            f"{q[i, 0]:12.5E}   {q[i, 1]:12.5E}   {q[i, 2]:12.5E}   "
            f"{t[i]:12.5E}   "
            f"{int(nf[i]):4d}   {int(ng[i]):4d}\n"
        )
        data_txt.write(line)



def preuci(u, nux, nms, state):
    """
    Перетворення переміщень у новій системі координат.
    """
    # Доступ до прапорців з COMMON/PRCILC/ через state
    sver = state.sver

    for i in range(nms):
        c = np.zeros((3, 3), dtype=float)
        grprck(i, c)

        if not sver:
            continue

        u1 = u[i, 0] * c[0, 0] + u[i, 1] * c[1, 0] + u[i, 2] * c[2, 0]
        u2 = u[i, 0] * c[0, 1] + u[i, 1] * c[1, 1] + u[i, 2] * c[2, 1]
        u3 = u[i, 0] * c[0, 2] + u[i, 1] * c[1, 2] + u[i, 2] * c[2, 2]

        u[i, 0] = u1
        u[i, 1] = u2
        u[i, 2] = u3


def prinke(u, x, t, nf, nux, nms, state, data_txt):
    """
    Друк результатів розрахунку: переміщень та напружень.
    """
    # Доступ до констант матеріалу з COMMON через state
    data_txt.write(f'E11={state.e11} E22={state.e22} E33={state.e33}\n')
    data_txt.write(f'G12={state.g21} G13={state.g31} G23={state.g23}\n')
    data_txt.write(f'G21={state.g21} G31={state.g31} G32={state.g23}\n')
    data_txt.write(f'V12={state.v12} V13={state.v13} V23={state.v23}\n')
    data_txt.write(f'V21={state.v21} V31={state.v31} V32={state.v32}\n')

    state.sih = 0.0

    # Вивід заголовків формату 10 та 21
    data_txt.write(
        '\n     ЗНАЧЕННЯ ПЕРЕМІЩЕНЬ ТА НАВАНТАЖЕНЬ В ЦЕНТРАХ КІНЦЕВИХ ЕЛЕМЕНТІВ У ДЕКАРТОВІЙ СИСТЕМІ КООРДИНАТ\n\n')
    data_txt.write(
        "    NU       U1           U2           U3           SG11         SG12         SG13         SG22         SG23         SG33\n")

    d = np.zeros(9, dtype=float)
    sig = np.zeros(9, dtype=float)
    en = np.zeros(9, dtype=float)
    sn = np.zeros(9, dtype=float)

    for i_idx in range(nms):
        # 1-based індекс для сумісності з логікою Фортрану
        i_fortran = i_idx + 1
        ks2 = kos(i_fortran, 2, state)

        state.a1 = 0.0

        if nf[i_idx] < 10:
            # Мітка 2
            sig.fill(0.0)
        elif ks2 == state.m2:
            # Мітка 2
            sig.fill(0.0)
        else:
            # Виклик SIGKE3
            sigke3(u, x, t, nf, i_fortran, nux, nms, d, sig, 0.0, en, sn)
            for j in range(9):
                sig[j] += sn[j] * state.dn

        # Розрахунок довжини вектора переміщень у площині та формування рядка виводу
        u1_val = u[i_idx]  # відповідає U(I)
        u2_val = u[i_idx + nux]  # U(I+NUX)
        u3_val = u[i_idx + nux * 2]  # U(I+NUX*2)
        u_norm = np.sqrt(u1_val ** 2 + u2_val ** 2)

        line = (
            f"  {i_fortran:4d}  "
            f"{u1_val:12.5e} {u2_val:12.5e} {u3_val:12.5e} {u_norm:12.5e} "
            f"{sig[1]:12.5e} {sig[2]:12.5e} {sig[4]:12.5e} {sig[5]:12.5e} {sig[8]:12.5e}\n"
        )
        data_txt.write(line)


def nuglob(x, nux, j, nf, ng, l, state):
    """
    Глобальна нумерація вузлів тривимірної конструкції.
    """
    m1 = state.m1
    m2 = state.m2
    m3 = state.m3

    nms = m1 * m2 * m3
    lg = 1 - l
    b_val = 0.000001

    for i in range(nms):
        ip = nf[i]
        if ip < 0:
            continue

        ngg = ng[i]
        if ngg > 0:
            continue

        lg += l
        ng[i] = lg

        ipp = i + 1
        if ipp >= nms:
            continue

        for ii in range(ipp, nms):
            irr = nf[ii]
            if irr < 0:
                continue

            ngg_ii = ng[ii]
            if ngg_ii > 0:
                continue

            jj = 1
            for kj in range(j):
                diff = abs(x[i, kj] - x[ii, kj])
                if diff > b_val:
                    jj = 0
                    break

            if jj == 1:
                ng[ii] = lg

    state.neq = lg + l - 1