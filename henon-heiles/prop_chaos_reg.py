import multiprocessing
import os

import numpy as np
from model import px_from_E
from pynamicalsys import HamiltonianSystem, TimeSeriesMetrics


def compute_one_ic(energy, idx, total_time, total_crossing, std_scale=0.05):
    """Worker function: generate one random IC and compute RTE and SALI."""
    ds = HamiltonianSystem(model="henon heiles")
    q = np.zeros(2, dtype=np.float64)
    p = np.zeros(2, dtype=np.float64)
    q[0] = 0.0

    # Random initial conditions for q[1] and p[1] using idx as a seed to get reproducible but independent draws
    # We'll use a local RNG seeded by energy and idx to avoid collisions.
    rng = np.random.default_rng(seed=int(energy * 1e6) + idx)
    while True:
        q[1] = rng.uniform(-0.5, 1.0)
        p[1] = rng.uniform(-0.7, 0.7)
        info, p[0] = px_from_E(energy, q, p[1])
        if info == 0:
            break

    # Compute Poincaré section
    full_pss = ds.poincare_section(q, p, total_crossing)
    reduced_pss = np.empty((total_crossing, 2), dtype=float)
    reduced_pss[:, 0] = full_pss[:, 2]
    reduced_pss[:, 1] = full_pss[:, 4]

    tsm = TimeSeriesMetrics(reduced_pss)
    rte = tsm.recurrence_time_entropy(threshold_mode="std", threshold=std_scale)
    _, sali = ds.SALI(q, p, total_time, threshold=1e-16)
    return rte, sali


def main():
    path_to_data = "data"
    os.makedirs(path_to_data, exist_ok=True)

    energies = np.linspace(1 / 10, 1 / 6, 50)[::-1]
    total_crossing = 10000
    total_time = 1e4
    num_ic = 1000
    rte_threshold = 2.5
    sali_threshold = 1e-10

    datafile = f"{path_to_data}/prop_chaos_vs_E_N={total_crossing}_T={total_time:.1f}_num_ic={num_ic}_rte_thr={rte_threshold:.2f}_sali_thr={sali_threshold:.1e}.dat"

    # Use a fixed seed for the main process (not directly needed)
    np.random.seed(1312)

    with open(datafile, "w") as df:
        for energy in energies:
            # Prepare arguments: (energy, idx) for each IC
            args = [(energy, i, total_time, total_crossing) for i in range(num_ic)]

            # Parallel computation of RTE and SALI for all ICs of this energy
            with multiprocessing.Pool() as pool:
                results = pool.starmap(compute_one_ic, args)

            # Unpack results
            rte_vals = np.array([res[0] for res in results])
            sali_vals = np.array([res[1] for res in results])

            # Count chaotic orbits by each criterion
            num_chaos_rte = np.sum(rte_vals > rte_threshold)
            num_chaos_sali = np.sum(sali_vals < sali_threshold)

            num_reg_rte = num_ic - num_chaos_rte
            num_reg_sali = num_ic - num_chaos_sali

            prop_chaos_rte = num_chaos_rte / num_ic
            prop_reg_rte = num_reg_rte / num_ic
            prop_chaos_sali = num_chaos_sali / num_ic
            prop_reg_sali = num_reg_sali / num_ic

            df.write(
                f"{energy:.16f} {prop_chaos_rte:.16e} {prop_reg_rte:.16e} {prop_chaos_sali:.16e} {prop_reg_sali:.16e}\n"
            )
            df.flush()


if __name__ == "__main__":
    main()
