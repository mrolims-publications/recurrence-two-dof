import numpy as np

# ---------------------------------
#  Hénon-Heiles system
# ---------------------------------


def henon_heiles_T(p, parameters=None):
    p0, p1 = p

    return 0.5 * (p0**2 + p1**2)


def henon_heiles_gradT(p, parameters=None):
    return p


def henon_heiles_hessT(p, parameters=None):
    return np.eye(2)


def henon_heiles_V(q, parameters=None):
    q0, q1 = q
    return 0.5 * (q0**2 + q1**2) + q0**2 * q1 - (q1**3) / 3.0


def henon_heiles_gradV(q, parameters=None):
    q0, q1 = q
    return np.array([q0 * (1.0 + 2.0 * q1), q1 + q0**2 - q1**2])


def henon_heiles_hessV(q, parameters=None):
    q0, q1 = q
    return np.array([[1.0 + 2.0 * q1, 2.0 * q0], [2.0 * q0, 1.0 - 2.0 * q1]])


def px_from_E(E, q, p1):
    arg = 2.0 * (E - henon_heiles_V(q)) - p1**2

    if arg > 0:
        p0 = np.sqrt(arg)
        return 0, p0
    else:
        return -1, None


def initial_conditions(x, y_range, py_range, E, num_ic):
    q = np.zeros((num_ic, 2))
    p = np.zeros((num_ic, 2))

    np.random.seed(1312)
    for i in range(num_ic):
        while True:
            y = np.random.uniform(*y_range)
            py = np.random.uniform(*py_range)
            info, px = px_from_E(E, [x, y], py)
            if info == 0:
                q[i] = [x, y]
                p[i] = [px, py]
                break

    return q, p
