import os

import numpy as np
from pynamicalsys import TimeSeriesMetrics


def exec_single_IC(ds, q, p, total_time, n_cross, std_scale):

    full_pss = ds.poincare_section(q, p, n_cross)
    reduced_pss = np.empty((n_cross, 2), dtype=float)
    reduced_pss[:, 0] = full_pss[:, 2]
    reduced_pss[:, 1] = full_pss[:, 4]

    tsm = TimeSeriesMetrics(reduced_pss)
    rte = tsm.recurrence_time_entropy(
        threshold_mode="std",
        threshold=std_scale,
    )

    sali = ds.SALI(q, p, total_time, threshold=1e-16)

    return rte, sali


def main():
    path_to_data = "data"
    os.makedirs(path_to_data, exist_ok=True)

    energies = np.linspace(1 / 10, 1 / 6, 50)
