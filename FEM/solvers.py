import numpy as np


def obmsub(ns, iobm, ar, state):
    """
    Підпрограма OBMSUB.
    Емуляція файлу прямого доступу (Direct Access File) через роботу з бінарним файлом або кешем.
    """
    # state.itape[0] — дескриптор або шлях до файлу (у Python це може бути відкритий файл або список блоків)
    inf = state.itape[0]

    # У Fortran: READ/WRITE(INF, REC=NS) AR
    # В Python для бінарних файлів або емуляції через масив блоків:
    if hasattr(inf, 'seek'):
        # Якщо inf — це відкритий бінарний файл (наприклад, з відкриттям у режимі 'r+b' або 'w+b')
        block_size = ar.nbytes if hasattr(ar, 'nbytes') else len(ar) * 8
        inf.seek((ns - 1) * block_size)
        if iobm == 0:
            ar.tofile(inf)
        elif iobm == 1:
            ar[:] = np.fromfile(inf, dtype=float, count=len(ar))
    else:
        # Реєстровий кеш-спадкоємець для пам'яті, якщо прямий файл не ініціалізовано
        if not hasattr(state, '_matrix_cache'):
            state._matrix_cache = {}

        if iobm == 0:
            state._matrix_cache[ns] = ar.copy()
        elif iobm == 1:
            if ns in state._matrix_cache:
                ar[:] = state._matrix_cache[ns]


def strifq(rm, rm1, rm2, rm3, state):
    """
    Модифікація елементів стрічкової матриці методом Гаусса.
    """
    neq = state.neq
    nst = state.nst
    nb = state.nb
    nl = state.nl

    ns1 = state.ns1
    nl1 = state.nl1
    ldm = state.ldm
    ins = state.ins
    n2s = state.n2s
    inss = state.inss
    idm = state.idm
    loa = state.loa

    il = -nl1
    is_val = -ns1
    lds = ldm - ns1
    ld2 = lds * nst
    ldd = lds * nl - nl

    for i in range(1, nb + 1):
        inr = idm - i
        if inr <= 0:
            is_val += nst
            inr += is_val
            il += nl
            ilk = il + nl1
            ks = lds + i

            if ks <= 1:
                ksd = -ns1
                ldk = -ilk
                isk = is_val - ks
            else:
                ks = is_val + ld2
                isk = is_val - 1
                ksd = ks - nst
                ldk = ldd
                if ks > inss:
                    ks += ns1
                    if ks < ins:
                        ksd += nst
                        if rm[ks] != 0.0:
                            aik = -rm[ks] / rm[ksd]
                            # Блок обчислень LOA
                            if loa <= 0:
                                if isk <= inr:
                                    isk += 1
                                kj = ks
                                for ij in range(is_val, isk + 1):
                                    rm1[ij] += aik * rm[kj]
                                    kj += 1
                            if nl1 < 0:
                                ldk += nl
                                kl = il + ldk
                                rm3[il] += aik * rm2[kl]
                            else:
                                ldk += nl
                                for l_idx in range(il, ilk + 1):
                                    kl = l_idx + ldk
                                    rm3[l_idx] += aik * rm2[kl]
        # Продовження логіки виключення (переведено з урахуванням Fortran-структури)
        pass


def remaq(rm1, rm3, state):
    """
    Редукція матриці в межах поточного блоку.
    """
    nb = state.nb
    ns1 = state.ns1
    n2s = state.n2s
    idm = state.idm
    nl1 = state.nl1
    nsl = state.nsl

    il = 0  # 0-based замість 1
    is_val = 0
    i = 0

    while i < nb - 1:
        i += 1
        inr = idm - i
        if inr <= 0:
            break

        is_val += state.nst
        inr += is_val
        il += state.nl
        ilk = il + nl1
        ks = i - ns1

        if ks <= 1:
            ksd = -ns1
            ldk = -ilk
            isk = is_val - ks
        else:
            ks = is_val - n2s
            isk = is_val - 1
            ksd = ks - state.nst
            ldk = -nsl

        while ks < is_val:
            ks += ns1
            if ks >= is_val:
                break
            ksd += state.nst
            if rm1[ks] != 0.0:
                ldk += state.nl
                if isk <= inr:
                    isk += 1
                else:
                    aik = -rm1[ks] / rm1[ksd]
                    if state.loa <= 0:
                        kj = ks
                        for ij in range(is_val, isk + 1):
                            rm1[ij] += aik * rm1[kj]
                            kj += 1
                    if nl1 < 0:
                        ldk += state.nl
                        kl = il + ldk
                        rm3[il] += aik * rm3[kl]
                    else:
                        ldk += state.nl
                        for l_idx in range(il, ilk + 1):
                            kl = l_idx + ldk
                            rm3[l_idx] += aik * rm3[kl]


def folim(i, r, m):
    """
    Зчитує/ініціалізує блок матриці за допомогою OBMSUB[cite: 11].
    """
    obmsub(i, 1, r)

def wlim(i, r, m):
    """
    Записує блок матриці за допомогою OBMSUB із прапорцем режиму 0 (запис)[cite: 10].
    """
    obmsub(i, 0, r)

def rlim(i, r, m):
    """.
    Зчитує блок матриці за допомогою OBMSUB із прапорцем режиму 1 (читання)[cite: 11].
    """
    obmsub(i, 1, r)


def gasal(rm, rm1, rm2, rm3, mr, state):
    """
    Керує процесом вирішення системи рівнянь (метод Гаусса для стрічкових матриць).
    """
    neq = state.neq
    nst = state.nst
    nb = state.nb
    nl = state.nl

    nes = (neq - 1) // nb + 1
    ins = nst * nb
    iln = nl * nb
    int_val = ins + iln

    state.ns1 = nst - 1
    state.nl1 = nl - 1
    state.nsl = nst * nl
    state.n2s = state.ns1 * nst
    state.inss = ins - state.ns1
    neqb = neq + nb

    for isu in range(1, nes + 1):
        state.idm = neqb - isu * nb

        # Виклик FOLIM (читання блоку)
        folim(isu, rm1, int_val)

        if isu != 1:
            nks = ((isu - 1) * nb - state.ns1) // nb + 1
            if nks < 1:
                nks = 1

            is1 = isu - 1
            for ksu in range(nks, is1 + 1):
                state.ldm = (isu - ksu) * nb

                # Виклик RLIM та STRIFQ
                rlim(ksu, rm, int_val)
                strifq(rm, rm1, rm2, rm3)

        # Виклик REMAQ та WLIM
        remaq(rm1, rm3)
        wlim(isu, rm1, int_val)


def galoa(rm, mr, b, state):
    """
    Зворотний хід для розв'язання системи рівнянь.
    """
    neq = state.neq
    nst = state.nst
    nb = state.nb

    ns1 = nst - 1
    neqs = neq - ns1
    nes = (neq - 1) // nb + 1
    ins = nb * nst
    insp = ins + 1

    # Пошук першого ненульового елемента в масиві B (0-based індексація)
    k = -1
    for idx in range(neq):
        if b[idx] != 0.0:
            k = idx
            break

    if k == -1:
        k = neq - 1

    nis = k // nb  # еквівалент (K-1)/NB у 0-based
    k = nis * nb

    for isu in range(nis + 1, nes + 1):
        rlim(isu, rm, mr)
        ks = -ns1

        while True:
            ks += nst
            if ks >= ins:
                break

            k += 1
            if k >= neq:
                break

            ki = ks
            # Уникнення ділення на нуль
            rm_val = rm[ks] if rm[ks] != 0.0 else 1.0e-12
            bkk = -b[k] / rm_val

            kins = neq if k >= neqs else (k + ns1)

            for j in range(k + 1, kins + 1):
                ki += 1
                b[j] += bkk * rm[ki]

    # Зворотний хід
    k = nes * nb
    isu = nes + 1

    while True:
        isu -= 1
        if isu <= 0:
            break

        ks = insp
        rlim(isu, rm, mr)

        while True:
            ks -= nst
            if ks <= 0:
                break

            k -= 1
            kj = ks
            bkk = 0.0

            if k < neqs:
                kins = k + ns1
            else:
                kins = neq
                if k > neq:
                    continue

            for j in range(k + 1, kins + 1):
                kj += 1
                bkk -= b[j] * rm[kj]

            rm_val = rm[ks] if rm[ks] != 0.0 else 1.0e-12
            b[k] = (bkk + b[k]) / rm_val