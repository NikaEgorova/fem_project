import numpy as np
from utils import nuglob, kos, telos, zakrep, hagpo3, nnnst3, parke


# ==========================================
# РЕАЛІЗОВАНІ ФУНКЦІЇ З ФАЙЛУ DARMIR
# ==========================================
def alfat(t):
    """КОЕФІЦІЄНТ ЛІНІЙНОГО ТЕМПЕРАТУРНОГО РОЗШИРЕННЯ"""
    return 8.2e-6


def kasat(xt, x1, x2, x3, state):
    nu = state.nn[0]
    ks1 = kos(nu, 1, state)

    if ks1 == 1:
        xt[0] = -x2
        xt[1] = x1
        xt[2] = 0.0
    elif ks1 == 2:
        xt[0] = 0.0
        xt[1] = 0.0
        xt[2] = 1.0
    return xt


def ugalf(x1, x2, x3):
    return 45.0 * np.pi / 180.0


def ugbet(x1, x2, x3):
    return 45.0 * np.pi / 180.0


def upro3(nu, ng, iz):
    return 0.0


def povna3(nu, ng, x, nux, jm, state):
    a = 0.0
    if ng == 1 and jm == 1:
        a = state.qnagr3
    return a


def korcil(x, nux, state):
    """Генерація циліндричних координат для вузлів"""
    pi2 = np.pi * 2.0

    # Уникнення ділення на нуль, якщо M-розмірності = 1
    m1_div = (state.m1 - 1) if state.m1 > 1 else 1
    m2_div = (state.m2 - 1) if state.m2 > 1 else 1
    m3_div = (state.m3 - 1) if state.m3 > 1 else 1

    # Використовуємо 0-based індексацію для циклів Python
    for n3 in range(state.m3):
        for n2 in range(state.m2):
            for n1 in range(state.m1):
                # Обчислення глобального індексу (1D), аналог Fortran: NF=N1+M1*(N2-1)+M1*M2*(N3-1)
                nf = n1 + state.m1 * n2 + state.m1 * state.m2 * n3

                r_val = state.rr + (state.ht / m1_div) * n1
                angle = (pi2 / m2_div) * n2
                z_val = (state.hl / m3_div) * n3

                x[nf, 0] = r_val * np.cos(angle)
                x[nf, 1] = r_val * np.sin(angle)
                x[nf, 2] = z_val


def danpri(state, data_txt):
    """ДРУК ДАНИХ ДОСЛІДЖУВАНИХ ЗАДАЧ"""
    data_txt.write('***ДАНІ ПРО РОЗВʼЯЗУВАНУ ЗАДАЧУ***\n')
    data_txt.write('***ГРЕБЕНЮК СЕРГІЙ МИКОЛАЙОВИЧ***\n')
    data_txt.write(f'ШИРИНА ПЛАСТИНИ PA2= {state.pa2}\n')
    data_txt.write(f'ШИРИНА ПЛАСТИНИ PA3= {state.pa3}\n')
    data_txt.write(f'HA - ТОВЩИНА АРМОВАНОГО ШАРУ {state.ha}\n')
    data_txt.write(f'ТОВЩИНА ПЛАСТИНИ PH= {state.ph}\n')
    data_txt.write(f'ІНТЕНСИВНІСТЬ НАВАНТАЖЕННЯ QNAGR3= {state.qnagr3}\n')
    data_txt.write(f'СІТКА РОЗБИТТЯ M1= {state.m1}\n')
    data_txt.write(f'СІТКА РОЗБИТТЯ M2= {state.m2}\n')
    data_txt.write(f'СІТКА РОЗБИТТЯ M3= {state.m3}\n')

    if state.smesi: data_txt.write('ТЕОРІЯ СУМІШЕЙ\n')
    if state.abolinsh: data_txt.write('ТЕОРІЯ АБОЛІНЬША\n')
    if state.vanin: data_txt.write('ТЕОРІЯ ВАНІНА\n')
    if state.dlina: data_txt.write('ТЕОРІЯ З УРАХУВАННЯМ ДОВЖИНИ ВОЛОКОН\n')
    if state.variaciam: data_txt.write('ВАРІАЦІЙНА ТЕОРІЯ (МЕТАЛИ)\n')
    if state.variaciap: data_txt.write('ВАРІАЦІЙНА ТЕОРІЯ (ПОЛІМЕРИ)\n')
    if state.kasparov: data_txt.write('ТЕОРІЯ КАСПАРОВА\n')
    if state.kristensen: data_txt.write('ТЕОРІЯ КРІСТЕНСЕНА\n')
    if state.klastornir: data_txt.write('ТЕОРІЯ КЛАСТОРНИ (РАДІАЛЬНІ)\n')
    if state.klastornio: data_txt.write('ТЕОРІЯ КЛАСТОРНИ (ОКРУЖНІ)\n')
    if state.klastornis: data_txt.write('ТЕОРІЯ КЛАСТОРНИ (ЛІНІЙНІ)\n')
    if state.energ: data_txt.write('ЕНЕРГЕТИЧНА ТЕОРІЯ\n')

    if state.mske != 2:
        data_txt.write('ТРАДИЦІЙНА СХЕМА СКІНЧЕНОГО ЕЛЕМЕНТА\n')
    elif state.mske == 2:
        data_txt.write('МОМЕНТНА СХЕМА СКІНЧЕНОГО ЕЛЕМЕНТА\n')

    data_txt.write(f'DC - ДІАМЕТР ВОЛОКОН АРМУВАННЯ {state.dc}\n')
    data_txt.write(f'ICO - ЧАСТОТА АРМУВАННЯ {state.ico}\n')
    data_txt.write(f'C1 - ОБʼЄМНИЙ ВМІСТ ВОЛОКОН {state.c1}\n')
    data_txt.write(f'EC11 - МОДУЛЬ ПРУЖНОСТІ ДЛЯ МАТЕРІАЛУ КОРДУ {state.ec11}\n')
    data_txt.write(f'EC22 - МОДУЛЬ ПРУЖНОСТІ ДЛЯ МАТЕРІАЛУ КОРДУ {state.ec22}\n')
    data_txt.write(f'GC12 - МОДУЛЬ ЗСУВУ ДЛЯ МАТЕРІАЛУ КОРДУ {state.gc12}\n')
    data_txt.write(f'VC12 - КОЕФІЦІЄНТ ПУАССОНА ДЛЯ МАТЕРІАЛУ КОРДУ {state.vc12}\n')
    data_txt.write(f'VC23 - КОЕФІЦІЄНТ ПУАССОНА ДЛЯ МАТЕРІАЛУ КОРДУ {state.vc23}\n')
    data_txt.write(f'ER11 - МОДУЛЬ ПРУЖНОСТІ МАТЕРІАЛУ МАТРИЦІ {state.er11}\n')
    data_txt.write(f'ER22 - МОДУЛЬ ПРУЖНОСТІ МАТЕРІАЛУ МАТРИЦІ {state.er22}\n')
    data_txt.write(f'GR12 - МОДУЛЬ ЗСУВУ МАТЕРІАЛУ МАТРИЦІ {state.gr12}\n')
    data_txt.write(f'VR12 - КОЕФІЦІЄНТ ПУАССОНА МАТЕРІАЛУ МАТРИЦІ {state.vr12}\n')
    data_txt.write(f'VR23 - КОЕФІЦІЄНТ ПУАССОНА МАТЕРІАЛУ МАТРИЦІ {state.vr23}\n')
    data_txt.write('*****************************************\n')

def darmir(x, nf, ng, t, q, nux, state, data_txt):
    """ГОЛОВНА ПІДПРОГРАМА ІНІЦІАЛІЗАЦІЇ ДАНИХ DARMIR"""
    # Значення DATA M,N,L/2,2,2/
    m, n, l = 2, 2, 2

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

    state.m1 = 3
    state.m2 = 9
    state.m3 = 4

    data_txt.write(f'M1 x M2 x M3 {state.m1} x {state.m2} x {state.m3}\n')

    state.ldk[0] = m
    state.ldk[1] = n
    state.ldk[2] = l
    state.ldk[6] = 1  # LDK(7) = 1

    state.ilis = 911
    state.nus = 1

    state.itape[0] = 8
    state.itape[1] = 1
    state.itape[2] = 2
    state.dn = 0.0

    state.mske = 2
    state.corect = 0.0

    state.dc = 0.07e-2
    state.ico = 1100.0

    state.ec11 = 1277.5
    state.ec22 = 1277.5
    state.vc23 = 0.3
    state.vc12 = 0.3
    state.gc12 = state.ec11 / (2.0 * (1.0 + state.vc12))

    state.er11 = 4.4
    state.er22 = state.er11
    state.vr12 = 0.4
    state.vr23 = state.vr12
    state.gr12 = state.er11 / (2.0 * (1.0 + state.vr12))

    state.ac11 = 1.0e-6
    state.ac22 = 1.0e-6
    state.am11 = 1.0e-6
    state.am22 = 1.0e-6

    ksloy = 1
    state.ha = 0.5e-2
    state.rr = 0.1
    state.hl = 0.2
    state.pa2 = 1.0
    state.ht = ksloy * state.ha
    state.dspus = 0.5

    state.c1 = (np.pi * (state.dc ** 2.0) / (4.0 * state.ha)) * state.ico
    state.c2 = 1.0 - state.c1

    danpri(state, data_txt)

    # Геометрія та граничні умови
    telos(1, 1, 1, state.m1, state.m2, state.m3, nf, state)

    # Окреме задання прапорця 71 для вузлів верхньої межі J=M2
    j = state.m2 - 1  # 0-based індекс для останнього M2
    for k in range(state.m3 - 1):
        for i in range(state.m1 - 1):
            nu = i + state.m1 * j + state.m1 * state.m2 * k
            nf[nu] = 71

    zakrep(1, 1, 1, state.m1, state.m2, 1, 7, nf, state)
    zakrep(1, 1, state.m3, state.m1, state.m2, state.m3, 7, nf, state)

    korcil(x, nux, state)

    hagpo3(1, 1, 1, 1, state.m2, state.m3, 1, nux, nf, x, q, povna3, state)

    j_val = 3
    ll = 3
    nuglob(x, nux, j_val, nf, ng, ll, state)

    nstr = nnnst3(nf, ng, ll, state)
    state.ldk[5] = nstr  # LDK(6) = NSTR

    parke(m, n, l, state)
    state.ldk[3] = 3  # LDK(4) = 3