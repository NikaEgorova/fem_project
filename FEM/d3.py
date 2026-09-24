import numpy as np
from utils import telos, zakrep, hagpo3, nuglob, nnnst3, parke


def darmir(x, nf, ng, t, q, nux, state, data_txt):
    """
    Варіант 3. "DARMIR" у 3-х мірній постановці
    """
    # Ініціалізація прапорців
    state.utem = False
    state.uvin = False
    state.upro = False
    state.preo = False
    state.ipec = 0

    state.smesi = True
    state.abolinsh = False
    state.vanin = False
    state.variaciam = False
    state.variaciap = False
    state.kasparov = False
    state.kristensen = False
    state.klastornir = False
    state.klastornio = False
    state.klastornis = False
    state.energ = False
    state.dlina = False

    # Розміри сіткової області
    state.m1 = 12
    m21 = 5
    state.m2 = 4 * m21 + 1
    state.m3 = 4

    data_txt.write(f"M1 x M2 x M3 {state.m1} x {state.m2} x {state.m3}\n")
    nms = state.m1 * state.m2 * state.m3

    state.ldk[0] = 2  # M
    state.ldk[1] = 2  # N
    state.ldk[2] = 2  # L
    state.ldk[6] = 1

    ma1 = state.m1 - 1
    ma2 = state.m2 - 1
    ma3 = state.m3 - 1

    state.ilis = 911
    state.nus = 1

    state.itape[0] = 8
    state.itape[1] = 1
    state.itape[2] = 2
    state.dn = 0.0

    state.mske = 2
    state.corect = 0.0

    pi = 4.0 * np.arctan(1.0)

    # Геометрія та матеріали
    state.dc = 0.07e-2
    state.ico = 1100.0

    state.ec11 = 1277.5
    state.ec22 = 1277.5
    state.vc23 = 0.3
    state.vc12 = 0.3
    state.gc12 = state.ec11 / (2.0 * (1.0 + state.vc12))

    state.er11 = 4.4
    state.er22 = state.er11
    state.vr12 = 0.49
    state.vr23 = state.vr12
    state.gr12 = state.er11 / (2.0 * (1.0 + state.vr12))

    state.ac11 = 1.0e-6
    state.ac22 = 1.0e-6
    state.am11 = 1.0e-6
    state.am22 = 1.0e-6

    ksloy = 1
    state.ha = 0.05
    state.rr = 0.1
    state.hl = 0.05
    state.pa2 = 1.0
    state.pa3 = getattr(state, 'pa3', 0.0)
    state.ph = getattr(state, 'ph', 0.0)
    state.ht = ksloy * state.ha

    # Інтенсивність навантаження (COMMON /QNAGR3/ DSPUS / QNAGR3)
    state.dspus = 1.0
    state.qnagr3 = state.dspus

    state.c1 = (pi * (state.dc ** 2.0) / (4.0 * state.ha)) * state.ico
    state.c1 = 0.2
    state.c2 = 1.0 - state.c1

    danpri(state, data_txt)
    telos(1, 1, 1, state.m1, state.m2, state.m3, nf, state)

    # Граничні умови
    zakrep(1, 1, 1, state.m1, state.m2, 1, 4, nf, state)
    zakrep(1, 1, state.m3, state.m1, state.m2, state.m3, 4, nf, state)
    zakrep(1, 1, 1, state.m1, 1, state.m3, 2, nf, state)
    zakrep(1, state.m2, 1, state.m1, state.m2, state.m3, 1, nf, state)

    korcil(x, nux, state)
    hagpo3(1, 1, 1, 1, state.m2, state.m3, 1, nux, nf, x, q, povna3, state)

    j = 3
    ll = 3
    nuglob(x, nux, j, nf, ng, ll, state)

    nstr = nnnst3(nf, ng, ll, state)
    state.ldk[5] = nstr  # LDK(6)

    parke(2, 2, 2, state)
    state.ldk[3] = 3  # LDK(4)

    return


def alfat(t_val):
    return 8.2e-6


def uvis3(nu, j, state):
    return 0.0


def upro3(nu, ng, iz, state):
    return 0.0


def povna3(nu, ng_val, x, nux, jm, state):
    a = 0.0
    if ng_val == 1 and jm == 1:
        a = state.qnagr3
    return a


def kasat(xt, x1, x2, x3, state):
    xt[0] = -x2
    xt[1] = x1
    xt[2] = 0.0
    return xt


def korcil(x, nux, state):
    pi2 = np.pi / 2.0

    for n1 in range(1, state.m1 + 1):
        for n2 in range(1, state.m2 + 1):
            for n3 in range(1, state.m3 + 1):
                nf_idx = n1 + state.m1 * (n2 - 1) + state.m1 * state.m2 * (n3 - 1)
                idx = nf_idx - 1

                x[idx, 0] = (state.rr + state.ht / (state.m1 - 1) * (n1 - 1)) * np.cos(pi2 / (state.m2 - 1) * (n2 - 1))
                x[idx, 1] = (state.rr + state.ht / (state.m1 - 1) * (n1 - 1)) * np.sin(pi2 / (state.m2 - 1) * (n2 - 1))
                x[idx, 2] = state.hl / (state.m3 - 1) * (n3 - 1)
    return


def danpri(state, file_out):
    file_out.write('***ДАНІ ПРО ВИРІШУВАНУ ЗАДАЧУ***\n')
    file_out.write('***ГРЕБЕНЮК СЕРГІЙ МИКОЛАЙОВИЧ***\n')
    file_out.write(f'ШИРИНА ПЛАСТИНИ PA2= {getattr(state, "pa2", 0.0)}\n')
    file_out.write(f'ШИРИНА ПЛАСТИНИ PA3= {getattr(state, "pa3", 0.0)}\n')
    file_out.write(f'HA - ТОВЩИНА АРМОВАНОГО ШАРУ= {state.ha}\n')
    file_out.write(f'ТОВЩИНА ПЛАСТИНИ PH= {getattr(state, "ph", 0.0)}\n')
    file_out.write(f'ІНТЕНСИВНІСТЬ НАВАНТАЖЕННЯ QNAGR3= {state.qnagr3}\n')
    file_out.write(f'СІТКА РОЗБИТТЯ M1= {state.m1}\n')
    file_out.write(f'СІТКА РОЗБИТТЯ M2= {state.m2}\n')
    file_out.write(f'СІТКА РОЗБИТТЯ M3= {state.m3}\n')

    if state.smesi: file_out.write('ТЕОРІЯ СУМІШЕЙ\n')
    if state.abolinsh: file_out.write('ТЕОРІЯ АБОЛІНЬША\n')
    if state.vanin: file_out.write('ТЕОРІЯ ВАНІНА\n')
    if state.dlina: file_out.write('ТЕОРІЯ З УРАХУВАННЯМ ДОВЖИНИ ВОЛОКОН\n')
    if state.variaciam: file_out.write('ВАРІАЦІЙНА ТЕОРІЯ (МЕТАЛИ)\n')
    if state.variaciap: file_out.write('ВАРІАЦІЙНА ТЕОРІЯ (ПОЛІМЕРИ)\n')
    if state.kasparov: file_out.write('ТЕОРІЯ КАСПАРОВА\n')
    if state.kristensen: file_out.write('ТЕОРІЯ КРІСТЕНСЕНА\n')
    if state.klastornir: file_out.write('ТЕОРІЯ КЛАСТОРНОЇ (РАДІАЛЬНІ)\n')
    if state.klastornio: file_out.write('ТЕОРІЯ КЛАСТОРНОЇ (ОКРУЖНІ)\n')
    if state.klastornis: file_out.write('ТЕОРІЯ КЛАСТОРНОЇ (ЛІНІЙНІ)\n')
    if state.energ: file_out.write('ЕНЕРГЕТИЧНА ТЕОРІЯ\n')

    if state.mske != 2:
        file_out.write('ТРАДИЦІЙНА СХЕМА СКІНЧЕННОГО ЕЛЕМЕНТА\n')
    else:
        file_out.write('МОМЕНТНА СХЕМА СКІНЧЕННОГО ЕЛЕМЕНТА\n')

    file_out.write(f'DC - ДІАМЕТР ВОЛОКОН АРМУВАННЯ= {state.dc}\n')
    file_out.write(f'ICO - ЧАСТОТА АРМУВАННЯ= {state.ico}\n')
    file_out.write(f'C1 - ОБ\'ЄМНИЙ ВМІСТ ВОЛОКОН= {state.c1}\n')
    file_out.write(f'EC11 - МОДУЛЬ ПРУЖНОСТІ ДЛЯ МАТЕРІАЛУ КОРДУ= {state.ec11}\n')
    file_out.write(f'EC22 - МОДУЛЬ ПРУЖНОСТІ ДЛЯ МАТЕРІАЛУ КОРДУ= {state.ec22}\n')
    file_out.write(f'GC12 - МОДУЛЬ ЗСУВУ ДЛЯ МАТЕРІАЛУ КОРДУ= {state.gc12}\n')
    file_out.write(f'VC12 - КОЕФІЦІЄНТ ПУАССОНА ДЛЯ МАТЕРІАЛУ КОРДУ= {state.vc12}\n')
    file_out.write(f'VC23 - КОЕФІЦІЄНТ ПУАССОНА ДЛЯ МАТЕРІАЛУ КОРДУ= {state.vc23}\n')
    file_out.write(f'ER11 - МОДУЛЬ ПРУЖНОСТІ МАТЕРІАЛУ МАТРИЦІ= {state.er11}\n')
    file_out.write(f'ER22 - МОДУЛЬ ПРУЖНОСТІ МАТЕРІАЛУ МАТРИЦІ= {state.er22}\n')
    file_out.write(f'GR12 - МОДУЛЬ ЗСУВУ МАТЕРІАЛУ МАТРИЦІ= {state.gr12}\n')
    file_out.write(f'VR12 - КОЕФІЦІЄНТ ПУАССОНА МАТЕРІАЛУ МАТРИЦІ= {state.vr12}\n')
    file_out.write(f'VR23 - КОЕФІЦІЄНТ ПУАССОНА МАТЕРІАЛУ МАТРИЦІ= {state.vr23}\n')
    file_out.write('*****************************************\n')