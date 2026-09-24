import numpy as np
from utils import telos, zakrep, hagpo3, nuglob, nnnst3, parke
from material_models import ukarpinos

def darmir(x, nf, ng, t, q, nux, state, data_txt):
    """
    Варіант 4. "DARMIR" В 3-X MEPНОЙ ПОСТАНОВКЕ
    """
    m, n, l = 2, 2, 2
    pi = 4.0 * np.arctan(1.0)

    # Признаки температури, зміщень, упр. осн., перетворення
    state.utem = False
    state.uvin = False
    state.upro = False
    state.preo = False
    state.ipec = 0

    state.smesi = False
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
    state.karpinos = True
    state.stolyarova = False

    # Розміри сіткової області
    state.m1 = 4
    m21 = 2
    state.m2 = 4 * m21 + 1
    state.m3 = 3

    data_txt.write(f"M1 x M2 x M3 {state.m1} x {state.m2} x {state.m3}\n")
    nms = state.m1 * state.m2 * state.m3

    state.ldk[0] = m
    state.ldk[1] = n
    state.ldk[2] = l
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

    # Облік моментної схеми
    state.mske = 2
    state.corect = 0.0

    # DC - діаметр волокон армування
    state.dc = 0.07e-2
    # ICO - частота армування
    state.ico = 1100.0

    # Пружні сталі матеріалу корда
    state.ec11 = 73.0
    state.ec22 = 73.0
    state.vc23 = 0.22
    state.vc12 = 0.22
    state.gc12 = state.ec11 / (2.0 * (1.0 + state.vc12))

    # Пружні сталі матеріалу матриці
    state.er11 = 3.24
    state.er22 = state.er11
    state.vr12 = 0.385
    state.vr23 = state.vr12
    state.gr12 = state.er11 / (2.0 * (1.0 + state.vr12))

    # Термопружні характеристики матеріалу корда
    state.ac11 = 1e-6
    state.ac22 = 1e-6

    # Термопружні характеристики матеріалу матриці
    state.am11 = 1e-6
    state.am22 = 1e-6

    # Кількість шарів
    ksloy = 1
    # HA - товщина армованого шару
    state.ha = 0.05
    state.rr = 0.1
    # Довжина конструкції
    state.hl = 0.01
    # Ширина конструкції
    state.pa2 = 1.0
    # Товщина конструкції
    state.ht = ksloy * state.ha
    # Інтенсивність навантаження
    state.dspus = 5.0

    state.c1 = 0.5
    state.c2 = 1.0 - state.c1
    state.cq = np.sqrt(0.5)

    # Задання топології тривимірної конструкції
    ukarpinos(state)

    data_txt.write(f"E1={state.e11} E2={state.e22}\n")
    data_txt.write(
        f"E3={state.e33} G12={state.g12} G23={state.g23} G13={state.g13}\n"
    )
    data_txt.write(
        f"V12={state.v12} V21={state.v21} V13={state.v13} V31={state.v31} V23={state.v23} V32={state.v32}\n"
    )

    arrr = -state.dspus * (1.0 - state.v12 * state.v21) * (
        state.rr ** (np.sqrt(state.e11 / state.e22) + 1.0)
    ) / (
        (np.sqrt(state.e22 * state.e11) - state.e11 * state.v21)
        * ((state.rr + state.ha) ** (2.0 * np.sqrt(state.e11 / state.e22)))
        + (np.sqrt(state.e11 * state.e22) + state.e11 * state.v21)
        * (state.rr ** (2.0 * np.sqrt(state.e11 / state.e22)))
    )
    urr = arrr * (
        state.rr ** (np.sqrt(state.e11 / state.e22))
        - ((state.rr + state.ha) ** (2.0 * np.sqrt(state.e11 / state.e22)))
        / (state.rr ** (np.sqrt(state.e11 / state.e22)))
    )
    data_txt.write(f"UR={urr} {arrr}\n")

    arrr = -state.dspus * (1.0 - state.v23 * state.v32) * (
        state.rr ** (np.sqrt(state.e33 / state.e22) + 1.0)
    ) / (
        (np.sqrt(state.e22 * state.e33) - state.e33 * state.v23)
        * ((state.rr + state.ha) ** (2.0 * np.sqrt(state.e33 / state.e22)))
        + (np.sqrt(state.e33 * state.e22) + state.e33 * state.v23)
        * (state.rr ** (2.0 * np.sqrt(state.e33 / state.e22)))
    )
    urr = arrr * (
        state.rr ** (np.sqrt(state.e33 / state.e22))
        - ((state.rr + state.ha) ** (2.0 * np.sqrt(state.e33 / state.e22)))
        / (state.rr ** (np.sqrt(state.e33 / state.e22)))
    )
    data_txt.write(f"UR={urr} {arrr}\n")

    danpri(state, data_txt)
    telos(1, 1, 1, state.m1, state.m2, state.m3, nf, state)

    # Граничні умови
    zakrep(1, 1, 2, state.m1 - 1, state.m2, 2, 4, nf, state)
    zakrep(1, 1, 1, state.m1 - 1, 1, state.m3, 2, nf, state)
    zakrep(1, state.m2, 1, state.m1 - 1, state.m2, state.m3, 1, nf, state)
    zakrep(state.m1, 1, 1, state.m1, state.m2, state.m3, 7, nf, state)

    korcil(x, nux, state)
    hagpo3(1, 1, 1, 1, state.m2, state.m3, 1, nux, nf, x, q, povna3, state)

    # Глобальна нумерація вузлів 3-х мер. кон-ції
    j = 3
    ll = 3
    state.neq = nuglob(x, nux, j, nf, ng, ll, state.neq)
    nstr = nnnst3(nf, ng, ll, state)
    state.ldk[5] = nstr  # LDK(6) у Fortran
    parke(m, n, l, state)
    state.ldk[3] = 3  # LDK(4) у Fortran
    return


def alfat(t):
    """КОЕФІЦІЄНТ ЛІНІЙНОГО ТЕМПЕРАТУРНОГО РОЗШИРЕННЯ ЗВ'ЯЗУЮЧОГО[cite: 14]"""
    return 8.2e-6


def uvis3(nu, j, state):
    return 0.0


def kasat(xt, x1, x2, x3, state):
    nu = state.nn[0]
    pi = 4.0 * np.arctan(1.0)
    xt[0] = 0.0
    xt[1] = 0.0
    xt[2] = 1.0
    return xt


def upro3(nu, ng, iz, state):
    return 0.0


def povna3(nu, ng, x, nux, jm, state):
    a = 0.0
    if ng == 1 and jm == 1:
        a = state.dspus
    return a


def korcil(x, nux, state):
    pi = 4.0 * np.arctan(1.0)
    pi2 = pi / 2.0
    for n1 in range(1, state.m1 + 1):
        for n2 in range(1, state.m2 + 1):
            for n3 in range(1, state.m3 + 1):
                nf_idx = n1 + state.m1 * (n2 - 1) + state.m1 * state.m2 * (n3 - 1)
                idx = nf_idx - 1  # 0-based індексація для Python
                x[idx, 0] = (
                    state.rr
                    + state.ht / (state.m1 - 1) * (n1 - 1)
                ) * np.cos(pi2 / (state.m2 - 1) * (n2 - 1))
                x[idx, 1] = (
                    state.rr
                    + state.ht / (state.m1 - 1) * (n1 - 1)
                ) * np.sin(pi2 / (state.m2 - 1) * (n2 - 1))
                x[idx, 2] = state.hl / (state.m3 - 1) * (n3 - 1)
    return x


def danpri(state, data_txt):
    """ДРУК ДАННИХ ДОСЛІДЖУВАНИХ ЗАДАЧ"""
    data_txt.write("***ДАННІ ПРО ДОСЛІДЖУВАНУ ЗАДАЧУ***\n")
    data_txt.write("***ГРЕБЕНЮК СЕРГІЙ МИКОЛАЙОВИЧ***\n")
    data_txt.write(f"ШИРИНА ПЛАСТИНИ PA2={state.pa2}\n")
    data_txt.write(f"ШИРИНА ПЛАСТИНИ PA3={state.pa3}\n")
    data_txt.write(f"HA - ТОВЩИНА АРМОВАНОГО ШАРУ {state.ha}\n")
    data_txt.write(f"ТОВЩИНА ПЛАСТИНИ PH={state.ph}\n")
    data_txt.write(f"ІНТЕНСИВНІСТЬ НАВАНТАЖЕННЯ QNAGR3={state.dspus}\n")
    data_txt.write(f"СІТКА РОЗБИТТЯ M1={state.m1}\n")
    data_txt.write(f"СІТКА РОЗБИТТЯ M2={state.m2}\n")
    data_txt.write(f"СІТКА РОЗБИТТЯ M3={state.m3}\n")

    if state.smesi:
        data_txt.write("ТЕОРІЯ СУМІШЕЙ\n")
    if state.abolinsh:
        data_txt.write("ТЕОРІЯ АБОЛІНЬША\n")
    if state.vanin:
        data_txt.write("ТЕОРІЯ ВАНІНА\n")
    if state.dlina:
        data_txt.write("ТЕОРІЯ З УРАХУВАННЯМ ДОВЖИНИ ВОЛОКОН\n")
    if state.variaciam:
        data_txt.write("ВАРІАЦІЙНА ТЕОРІЯ (МЕТАЛИ)\n")
    if state.variaciap:
        data_txt.write("ВАРІАЦІЙНА ТЕОРІЯ (ПОЛІМЕРИ)\n")
    if state.kasparov:
        data_txt.write("ТЕОРІЯ КАСПАРОВА\n")
    if state.kristensen:
        data_txt.write("ТЕОРІЯ КРІСТЕНСЕНА\n")
    if state.klastornir:
        data_txt.write("ТЕОРІЯ КЛАСТОРНИ (РАДІАЛЬНІ)\n")
    if state.klastornio:
        data_txt.write("ТЕОРІЯ КЛАСТОРНИ (ОКРУЖНІ)\n")
    if state.klastornis:
        data_txt.write("ТЕОРІЯ КЛАСТОРНИ (ЛІНІЙНІ)\n")
    if state.energ:
        data_txt.write("ЕНЕРГЕТИЧНА ТЕОРІЯ\n")

    if state.mske != 2:
        data_txt.write("ТРАДИЦІЙНА СХЕМА КІНЦЕВОГО ЕЛЕМЕНТА\n")
    if state.mske == 2:
        data_txt.write("МОМЕНТНА СХЕМА КІНЦЕВОГО ЕЛЕМЕНТА\n")

    data_txt.write(f"DC - ДІАМЕТР ВОЛОКОН АРМУВАННЯ {state.dc}\n")
    data_txt.write(f"ICO - ЧАСТОТА АРМУВАННЯ {state.ico}\n")
    data_txt.write(f"C1 - ОБ'ЄМНИЙ ВМІСТ ВОЛОКОН {state.c1}\n")
    data_txt.write(
        f"EC11 - МОДУЛЬ ПРУЖНОСТІ ДЛЯ МАТЕРІАЛУ КОРДА {state.ec11}\n"
    )
    data_txt.write(
        f"EC22 - МОДУЛЬ ПРУЖНОСТІ ДЛЯ МАТЕРІАЛУ КОРДА {state.ec22}\n"
    )
    data_txt.write(
        f"GC12 - МОДУЛЬ ЗСУВУ ДЛЯ МАТЕРІАЛУ КОРДА {state.gc12}\n"
    )
    data_txt.write(
        f"VC12 - КОЕФІЦІЄНТ ПУАССОНА ДЛЯ МАТЕРІАЛУ КОРДА {state.vc12}\n"
    )
    data_txt.write(
        f"VC23 - КОЕФІЦІЄНТ ПУАССОНА ДЛЯ МАТЕРІАЛУ КОРДА {state.vc23}\n"
    )
    data_txt.write(
        f"ER11 - МОДУЛЬ ПРУЖНОСТІ МАТЕРІАЛУ МАТРИЦІ {state.er11}\n"
    )
    data_txt.write(
        f"ER22 - МОДУЛЬ ПРУЖНОСТІ МАТЕРІАЛУ МАТРИЦІ {state.er22}\n"
    )
    data_txt.write(
        f"GR12 - МОДУЛЬ ЗСУВУ МАТЕРІАЛУ МАТРИЦІ {state.gr12}\n"
    )
    data_txt.write(
        f"VR12 - КОЕФІЦІЄНТ ПУАССОНА МАТЕРІАЛУ МАТРИЦІ {state.vr12}\n"
    )
    data_txt.write(
        f"VR23 - КОЕФІЦІЄНТ ПУАССОНА МАТЕРІАЛУ МАТРИЦІ {state.vr23}\n"
    )
    data_txt.write("*****************************************\n")