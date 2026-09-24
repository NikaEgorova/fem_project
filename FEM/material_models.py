


def ukarpinos(state):
    """
    Обчислення ефективних пружних характеристик односпрямовано 
    армованого волокнистого матеріалу за мікромеханічною моделлю Карпіноса.
    
    Результати записуються безпосередньо в об'єкт state (E11..E33, G12..G23, V12..V32).
    """
    # 1. Зчитування вихідних параметрів та концентрацій
    c1 = state.c1
    c2 = 1.0 - c1
    state.c2 = c2
    cq = getattr(state, "cq", 0.0)  # Відносна порожнистість волокна (CQ)

    ec = state.ec11
    vc = state.vc12
    er = state.er11
    vr = state.vr12

    cq2 = cq**2

    # 2. Модулі зсуву волокна (PARG1) та матриці (PARG2)
    parg1 = ec / (2.0 * (1.0 + vc))
    parg2 = er / (2.0 * (1.0 + vr))
    parg_ratio = parg2 / parg1

    qk1 = 3.0 - 4.0 * vr
    qk2 = 3.0 - 4.0 * vc

    # 3. Модуль Юнга вздовж волокон E11
    e11 = ec * c1 * (1.0 - cq2) + er * c2

    # 4. Поздовжній коефіцієнт Пуассона V12, V13
    v12_denom = (1.0 - cq2) * (1.0 + c2 + c1 * qk1) + c2 * (qk2 - 1.0 + 2.0 * cq2) * parg_ratio
    v12 = vr + c1 * (vc - vr) * (qk1 + 1.0) * (1.0 - cq2) / v12_denom
    v13 = v12

    # 5. Допоміжні структурні коефіцієнти EQ, EN, EL
    eq_num = c1 * (1.0 - cq2 - (1.0 + qk2 * cq2) * parg_ratio)
    eq_den = (1.0 - cq2) * qk1 + (1.0 + qk2 * cq2) * parg_ratio
    eq = eq_num / eq_den

    en_den = 2.0 * (1.0 - cq2) + (qk2 - 1.0 + 2.0 * cq2) * parg_ratio
    en = c2 + c1 * (qk1 + 1.0) * (1.0 - cq2) / en_den

    el_den = 1.0 - cq2 - (1.0 + qk2 * cq2) * parg_ratio
    el = 1.0 / (-c2 + (qk1 + 1.0) * (1.0 - cq2) / el_den)

    # 6. Поперечні модулі Юнга E22, E33
    e22_inv = (v12**2 / e11) + ((qk1 + 1.0) / (8.0 * parg2)) * (1.0 / en - 2.0 * eq / (1.0 + eq))
    e22 = 1.0 / e22_inv
    e33 = e22

    # 7. Поперечний коефіцієнт Пуассона V32, V23
    v32_term = -v12**2 / e11 + vr / (2.0 * parg2) + ((qk1 + 1.0) / (8.0 * parg2)) * ((en - 1.0) / en - 2.0 * eq / (1.0 + eq))
    v32 = e22 * v32_term
    v23 = v32

    # 8. Модулі зсуву G12, G13, G23
    g12_num = (1.0 + c1) * (1.0 - cq2) + c2 * (1.0 + cq2) * parg_ratio
    g12_den = c2 * (1.0 - cq2) + (1.0 + c1) * (1.0 + cq2) * parg_ratio
    g12 = parg2 * (g12_num / g12_den)
    g13 = g12

    g23 = parg2 / (1.0 - c1 * (qk1 + 1.0) * el)

    # 9. Зворотні коефіцієнти Пуассона V21, V31
    v21 = e22 * v12 / e11
    v31 = e33 * v13 / e11

    # 10. Запис обчислених характеристик у об'єкт стану
    state.e11, state.e22, state.e33 = e11, e22, e33
    state.g12, state.g13, state.g23 = g12, g13, g23
    state.v12, state.v13, state.v23 = v12, v13, v23
    state.v21, state.v31, state.v32 = v21, v31, v32