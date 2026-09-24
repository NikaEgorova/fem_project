import numpy as np
import scipy.linalg

def iprzak(iz):
    val = 7 - (int(iz) % 64)
    if val <= 0:
        return 0, 0, 0, 0

    return 0, val & 1, (val >> 1) & 1, (val >> 2) & 1

def obratm(a: np.ndarray) -> np.ndarray:
    """Обернення матриці in-place та повернення результату."""
    inv_a = scipy.linalg.inv(a, overwrite_a=True)
    if a.shape == inv_a.shape:
        a[:] = inv_a
    return inv_a


def stx(i: int, a: float) -> float:
    """Аналог FUNCTION STX (для інших модулів, де вона викликається точково)."""
    return float(a ** i) if i > 0 else 1.0

def build_al_vectorized(m1: int, m2: int, m3: int) -> np.ndarray:
    """Побудова матриці AL через Кронекерів добуток (замінює 6 циклів і STX)."""
    o1 = np.linspace(-1.0, 1.0, m1)
    o2 = np.linspace(-1.0, 1.0, m2)
    o3 = np.linspace(-1.0, 1.0, m3)

    v1 = o1[:, None] ** np.arange(m1)
    v2 = o2[:, None] ** np.arange(m2)
    v3 = o3[:, None] ** np.arange(m3)

    return np.kron(v3, np.kron(v2, v1))

def alpro(al: np.ndarray, jn: int, la, ma, state):
    """
    Ультра-швидка версія ALPRO.
    Працює коректно як з 1D (shape=(64,)), так і з 2D (shape=(8,8)) масивами al.
    """
    m1, m2, m3 = state.m1, state.m2, state.m3
    
# 1. Швидка генерація AL (матриця JN x JN)
    al_mat = build_al_vectorized(m1, m2, m3)

    # 2. Обернення матриці AL
    al_inv = scipy.linalg.inv(al_mat)

    # 3. Безаварійне записування в al (підтримує і 1D, і 2D масиви)
    al.flat[:] = al_inv.flat


def nmak1(nf, ng, x, t, q, nux, state):
    """
    Ініціалізація робочих масивів та виклик розрахункових
    блоків ALPRO і FOMAKR.
    
    :param nf: Масив прапорців граничних умов
    :param ng: Масив номерів рівнянь
    :param x: Масив координат вузлів [NUX, 3]
    :param t: Масив температур [NUX]
    :param q: Масив вузлових потоків/навантажень [NUX, 3]
    :param nux: Максимальна кількість вузлів
    :param state: Об'єкт стану (State), що містить параметри m1, m2, m3 (COMMON/KSO/)
    """
    # 1. Отримання розмірностей із COMMON/KSO/ (у state це m1, m2, m3)
    m = state.m1
    n = state.m2
    l = state.m3

    # 2. Обчислення загальних розмірностей
    mn = m * n
    mnl = mn * l
    mnl3 = mnl * 3
    mat = 300

    # 3. Виділення робочих буферних масивів (точні відповідники DIMENSION з Fortran)
    rmn = np.zeros(mat, dtype=np.float64)
    al = np.zeros(64, dtype=np.float64)
    la = np.zeros(8, dtype=np.int32)
    ma = np.zeros(8, dtype=np.int32)

    epo = np.zeros(144, dtype=np.float64)
    sio = np.zeros(144, dtype=np.float64)
    w = np.zeros(24, dtype=np.float64)
    xe = np.zeros(24, dtype=np.float64)
    vt = np.zeros(8, dtype=np.float64)
    bt = np.zeros(8, dtype=np.float64)
    v = np.zeros(24, dtype=np.float64)

    # 4. Виклик підпрограми формування геометрії/зв'язків ALPRO
    alpro(al, mnl, la, ma, state)

    # 5. Виклик основного блоку формування макроелементів FOMAKR
    fomakr(
        nf, ng, x, t, q, nux, rmn, mat,
        al, epo, sio, mnl3, w, xe, vt, bt, v, mn, mnl,
        state
    )

def nmakm(nf, ng, x, t, q, nux, state):
    """
    Обгортка для виклику підпрограми формування матриці NMAK1.
    Передає керування та масиви до основної підпрограми нумерації/формування сітки NMAK1.

    :param nf: Масив прапорців граничних умов
    :param ng: Масив номерів рівнянь / ступенів вільності
    :param x: Масив координат вузлів [NUX, 3]
    :param t: Масив узлових температур [NUX]
    :param q: Масив вузлових потоків / навантажень [NUX, 3]
    :param nux: Максимально припустима кількість вузлів
    :param state: Об'єкт стану (State)
    """
    return nmak1(nf, ng, x, t, q, nux, state)
