import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from model import initial_conditions
from pynamicalsys import HamiltonianSystem


def progress_bar(done, total, width=40):
    frac = done / total
    filled = int(width * frac)
    bar = "#" * filled + "-" * (width - filled)
    print(f"\r    [{bar}] {done}/{total}", end="", flush=True)


def compute_single_pss(qi, pi, n_points):
    ds = HamiltonianSystem(model="henon heiles")
    pss = ds.poincare_section(qi, pi, n_points)
    return pss[:, 2], pss[:, 4]


if __name__ == "__main__":
    path_to_data = "data"
    os.makedirs(path_to_data, exist_ok=True)

    E_vals = [1 / float(i) for i in range(6, 11)]
    num_ic = 200
    n_points = 5000
    y_range = (-0.6, 0.6)
    py_range = (-0.75, 0.75)

    max_workers = os.cpu_count()

    print("\n=== GLOBAL PARAMETERS ===")
    print("Model           : Henon-Heiles")
    print(f"E values        : {E_vals}")
    print(f"# ICs           : {num_ic}")
    print(f"# points / IC   : {n_points}")
    print(f"y range         : {y_range}")
    print(f"py range        : {py_range}")
    print(f"Output path     : {path_to_data}")
    print(f"Workers         : {max_workers}")
    print("=========================\n")

    times_per_E = {}

    for E in E_vals:
        print(f"\n>>> Starting E = {E:.5f}")
        start_E = time.perf_counter()

        q, p = initial_conditions(0.0, y_range, py_range, E, num_ic)
        datafile = f"{path_to_data}/pss_E={E:.5f}_n={n_points}.dat"
        pss_data = np.zeros((n_points * num_ic, 2), dtype=float)

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(compute_single_pss, q[i], p[i], n_points): i
                for i in range(num_ic)
            }

            done = 0
            for future in as_completed(futures):
                i = futures[future]
                col0, col1 = future.result()

                pss_data[i * n_points : (i + 1) * n_points, 0] = col0
                pss_data[i * n_points : (i + 1) * n_points, 1] = col1

                done += 1
                progress_bar(done, num_ic)

        print()
        np.savetxt(datafile, pss_data)

        elapsed = time.perf_counter() - start_E
        times_per_E[E] = elapsed

        print(f"<<< Finished E = {E:.5f} in {elapsed:.2f} s")
        print(f"Saved to: {datafile}")

    print("\n=== TIMING SUMMARY ===")
    for E, t in times_per_E.items():
        print(f"E = {E:.5f} -> {t:.2f} s")
    print("======================\n")
