import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from model import px_from_E
from pynamicalsys import HamiltonianSystem


def progress_bar(done, total, width=40):
    frac = done / total
    filled = int(width * frac)
    bar = "#" * filled + "-" * (width - filled)
    print(f"\r    [{bar}] {done}/{total}", end="", flush=True)


def compute_lle_for_y(i, y_val, energy, t_total):
    start = time.perf_counter()

    q = np.array([0.0, y_val], dtype=float)
    p = np.zeros(2, dtype=float)

    info, px = px_from_E(energy, q, p[1])
    if info != 0:
        elapsed = time.perf_counter() - start
        return i, y_val, np.nan, elapsed, info

    p[0] = px
    ds = HamiltonianSystem(model="henon heiles")
    lle_val = ds.lyapunov(q, p, t_total, num_exponents=1)[0]

    elapsed = time.perf_counter() - start
    return i, y_val, lle_val, elapsed, 0


def main():
    path_to_data = "data"
    os.makedirs(path_to_data, exist_ok=True)

    energy = 1 / 8
    t_total = np.float64(1e5)
    y_range = (-0.5, 1.0)
    num_y = 5000
    y = np.linspace(*y_range, num_y)

    lle = np.full(num_y, np.nan, dtype=float)
    times_y = np.full(num_y, np.nan, dtype=float)
    info_y = np.full(num_y, -999, dtype=int)

    max_workers = os.cpu_count() if len(sys.argv) == 1 else int(sys.argv[1])

    datafile = (
        f"{path_to_data}/lle_vs_y_E={energy:.5f}_T={t_total:.0f}_num_y={num_y}.dat"
    )

    print("\n=== GLOBAL PARAMETERS ===")
    print("Model            : Henon-Heiles")
    print(f"Energy           : {energy:.5f}")
    print(f"t_total          : {t_total:.0f}")
    print(f"y range          : {y_range}")
    print(f"num_y            : {num_y}")
    print(f"Workers          : {max_workers}")
    print(f"Output data      : {datafile}")
    print("=========================\n")

    total_start = time.perf_counter()
    print(">>> Starting scan over y")

    done = 0
    n_invalid = 0
    valid_times = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(compute_lle_for_y, i, y[i], energy, t_total): i
            for i in range(num_y)
        }

        for future in as_completed(futures):
            i, y_val, lle_val, elapsed, info = future.result()

            lle[i] = lle_val
            times_y[i] = elapsed
            info_y[i] = info

            if info != 0:
                n_invalid += 1
            else:
                valid_times.append(elapsed)

            done += 1
            progress_bar(done, num_y)

    print()
    total_elapsed = time.perf_counter() - total_start

    valid_times = np.array(valid_times, dtype=float)
    n_valid = valid_times.size

    if n_valid > 0:
        mean_valid = np.mean(valid_times)
        std_valid = np.std(valid_times)
        total_valid_time = np.sum(valid_times)
    else:
        mean_valid = np.nan
        std_valid = np.nan
        total_valid_time = np.nan

    np.savetxt(
        datafile,
        np.column_stack((y, lle)),
        fmt="%.16e",
    )

    print("\n=== RUN SUMMARY ===")
    print(f"Total y values    : {num_y}")
    print(f"Valid ICs         : {n_valid}")
    print(f"Invalid ICs       : {n_invalid}")
    print(f"Mean time (valid) : {mean_valid:.6f} s")
    print(f"Std  time (valid) : {std_valid:.6f} s")
    print(f"Total valid time  : {total_valid_time:.2f} s")
    print(f"Total wall time   : {total_elapsed:.2f} s")
    print(f"Saved data to     : {datafile}")
    print("===================\n")


if __name__ == "__main__":
    main()
