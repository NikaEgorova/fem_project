import numpy as np
from d1 import darmir as darmir_var1
from d2 import darmir as darmir_var2
from d3 import darmir as darmir_var3
from d4 import darmir as darmir_var4
from nmaq import nmaq1
from solvers import gasal, iprzak, galoa
from utils import uvis3, pecda3, preuci, prinke
from mesh import nmakm, iprzak 


class FEMState:
    """Клас для зберігання глобальних змінних, які раніше були в COMMON блоках."""

    def __init__(self):
        # COMMON/IUPR/ та інші загальні змінні
        self.mre = np.zeros(4, dtype=int)
        self.ldk = np.zeros(8, dtype=int)
        self.itape = np.zeros(4, dtype=int)
        self.neq = 0
        self.m1 = 0
        self.m2 = 0
        self.m3 = 0
        self.ilis = 911
        self.nus = 1
        self.qnagr3 = 0.0

        # Логічні прапорці термопружності
        self.uvin = False
        self.utem = False
        self.preo = False
        self.sver = False
        self.upro = False

        # Загальні параметри системи
        self.is_val = 1
        self.loa = 0
        self.mcx = 0
        self.neqq = 0
        self.nst = 0
        self.nb = 0
        self.nzm = 0

        self.ipec = 0
        self.dn = 0.0
        self.nel = 0
        self.eps = 0.0
        self.mske = 2
        self.corect = 0.0

        # Геометрія
        self.pa2 = 0.0
        self.pa3 = 0.0
        self.ph = 0.0
        self.ro = 0.0
        self.rr = 0.0
        self.ht = 0.0
        self.hl = 0.0
        self.dc = 0.0
        self.ha = 0.0
        self.ico = 0.0

        # Упругі константи корду
        self.ec11 = 0.0
        self.ec22 = 0.0
        self.gc12 = 0.0
        self.vc12 = 0.0
        self.vc23 = 0.0

        # Упругі константи матриці
        self.er11 = 0.0
        self.er22 = 0.0
        self.gr12 = 0.0
        self.vr12 = 0.0
        self.vr23 = 0.0

        # Термоупругі характеристики
        self.ac11 = 0.0
        self.ac22 = 0.0
        self.am11 = 0.0
        self.am22 = 0.0
        self.ilis = 0.0
        self.c1 = 0.0
        self.c2 = 0.0
        self.qnagr3 = 0.0
        self.dspus = 0.0

        # Теорії та прапорці
        self.smesi = False
        self.abolinsh = False
        self.vanin = False
        self.variaciam = False
        self.variaciap = False
        self.kasparov = False
        self.kristensen = False
        self.klastornir = False
        self.klastornio = False
        self.klastornis = False
        self.energ = False
        self.dlina = False

        self.nn = np.zeros(64, dtype=int)  # Для KASAT

        # Додаткові прапорці для Варіанта 4 (теорія Карпіноса)
        self.karpinos = False

        self.rm_storage = None

# ==========================================
# ОСНОВНА ПРОГРАМА
# ==========================================
def main():
    state = FEMState()

    print("Оберіть варіант розрахунку (1-4):")
    print("1 - Варіант 1 (Циліндрична сітка, теорія сумішей, MSKE=2, M1=3, M2=9, M3=4)")
    print("2 - Варіант 2 (Прямокутна сітка, MSKE=1, C1=0.01)")
    print("3 - Варіант 3 (Деталізована сітка M1=12, M2=21, M3=4, VR12=0.49, C1=0.2)")
    print("4 - Варіант 4 (Теорія Карпіноса, SMESI=False, KARPINOS=True, C1=0.5)")

    choice = input("Введіть номер варіанта (1-4): ").strip()

    # Вибір відповідної функції darmir залежно від введеного значення
    if choice == '1':
        darmir_func = darmir_var1
        print("-> Обрано Варіант 1")
    elif choice == '2':
        darmir_func = darmir_var2
        print("-> Обрано Варіант 2")
    elif choice == '3':
        darmir_func = darmir_var3
        print("-> Обрано Варіант 3")
    elif choice == '4':
        darmir_func = darmir_var4
        print("-> Обрано Варіант 4 (Теорія Карпіноса)")
    else:
        print("Неправильний вибір! За замовчуванням використовується Варіант 1.")
        darmir_func = darmir_var1

    # Оголошення та ініціалізація робочих масивів (DATA оператори)
    nux = 1500

    # nf, ng заповнюються -1, t - нулями
    nf = np.full(nux, -1, dtype=int)
    ng = np.full(nux, -1, dtype=int)
    t = np.zeros(nux, dtype=float)

    # RM, RM1 ініціалізуються нулями
    rm = np.zeros(911, dtype=float)
    rm1 = np.zeros(911, dtype=float)

    # Q, X, UD ініціалізуються нулями (розмірність 1500x3)
    q = np.zeros((nux, 3), dtype=float)
    x = np.zeros((nux, 3), dtype=float)
    ud = np.zeros((nux, 3), dtype=float)

    # U ініціалізується нулями
    u = np.zeros(4500, dtype=float)

    # Задання параметрів файлів запису
    data_txt = open('DATA.TXT', 'w')

    # Задання загальних параметрів системи
    state.is_val = 1
    nl = 0
    state.loa = 0

    # Виклик обраної функції DARMIR
    darmir_func(x, nf, ng, t, q, nux, state, data_txt)

    # Призначення змінних після виконання DARMIR (зсув індексу масивів на -1)
    state.mcx = state.mre[1]  # У Fortran: MRE(2)
    state.neqq = state.neq
    state.nst = state.ldk[5]  # У Fortran: LDK(6)
    state.nb = state.ldk[4]  # У Fortran: LDK(5)

    nstn = state.nst * state.nb + nl * state.nb
    nsb = state.nst * state.nb + 1
    nms = state.m1 * state.m2 * state.m3

    # Уникнення ділення на нуль, якщо DARMIR ще не ініціалізував state.nb
    if state.nb != 0:
        nes = (state.neq - 1) // state.nb + 1
    else:
        nes = 0

    # СТВОРЕННЯ СХОВИЩА В RAM
    state.rm_storage = np.zeros((nes + 2, 911), dtype=np.float64)

    nmaq1(nf, ng, x, t, q, nux, state)
    state.utem = False

    # Друк табличних даних
    pecda3(x, q, t, nf, ng, nux, nms, state, data_txt)
    print('КОНТРОЛЬ PECDA3')

    # Заповнення сховища замість запису у бінарний файл
    for i in range(1, nes + 2):
        state.rm_storage[i] = rm.copy()

    # Формування матриці
    nmakm(nf, ng, x, t, q, nux, state)
    print('КОНТРОЛЬ NMAKM')

    print(f'РІШЕННЯ СИСТЕМИ {nes}')
    # Передаємо індекси початку масивів замість вказівників
    gasal(rm, rm1, 0, 0, nsb - 1, state)

    # Права частина термопружності
    for i in range(nms):
        ip = nf[i]
        if ip < 0:
            continue

        ipz, ii, jj, kk = iprzak(ip)

        # Індексація Python (0, 1, 2 замість 1, 2, 3)
        if ii == 1: q[i, 0] = 0.0
        if jj == 1: q[i, 1] = 0.0
        if kk == 1: q[i, 2] = 0.0

    for i in range(nms):
        ip = nf[i]
        if ip < 0:
            continue

        mu = ng[i] - 1  # Зсув до 0-based індексації
        for j in range(3):
            nu = mu + j
            u[nu] = u[nu] + q[i, j]

    print('ЗВОРОТНІЙ ХІД')
    galoa(rm, nstn, u, state)

    umax = 0.0
    umin = 0.0
    numax = 0
    namax = 0
    numin = 0
    namin = 0

    for i in range(nms):
        mu = ng[i] - 1
        ip = nf[i]
        if ip < 0:
            continue

        ipz, ii, jj, kk = iprzak(ip)

        for j in range(3):
            nu = mu + j
            ud[i, j] = u[nu]

            if ud[i, j] > umax:
                numax = i + 1  # Повертаємо 1-based індекси для логів
                namax = j + 1
                umax = ud[i, j]

            if ud[i, j] < umin:
                numin = i + 1
                namin = j + 1
                umin = ud[i, j]

            if not state.uvin:
                continue

            if j == 0 and ii == 1: ud[i, 0] = uvis3(i, 0, state)
            if j == 1 and jj == 1: ud[i, 1] = uvis3(i, 1, state)
            if j == 2 and kk == 1: ud[i, 2] = uvis3(i, 2, state)

    if state.preo:
        preuci(ud, nux, nms, state)

    # Друк результатів
    prinke(ud, x, t, nf, nux, nms, state, data_txt)

    # Запис мінімумів/максимумів у текстовий файл
    data_txt.write(f'umin= {umin} {numin} {namin}\n')
    data_txt.write(f'umax= {umax} {numax} {namax}\n')

    state.nzm = 1

    # Закриття файлів
    data_txt.close()


if __name__ == '__main__':
    main()