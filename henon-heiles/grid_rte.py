import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from model import px_from_E
from pynamicalsys import HamiltonianSystem, TimeSeriesMetrics


def progress_bar(done, total, width=40):
    frac = done / total
    filled = int(width * frac)
    bar = "#" * filled + "-" * (width - filled)
    print(f"\r    [{bar}] {done}/{total}", end="", flush=True)


def compute_rte_row(i, y_val, py_vals, energy, n_cross, std_scale):
    start = time.perf_counter()

    ds = HamiltonianSystem(model="henon heiles")

    n_py = py_vals.size
    rte_row = np.full(n_py, np.nan, dtype=float)
    info_row = np.full(n_py, -999, dtype=int)

    q = np.array([0.0, y_val], dtype=float)

    for j, py_val in enumerate(py_vals):
        p = np.array([0.0, py_val], dtype=float)

        info, px = px_from_E(energy, q, p[1])
        if info != 0:
            info_row[j] = info
            continue

        p[0] = px

        full_pss = ds.poincare_section(q, p, n_cross)

        reduced_pss = np.empty((n_cross, 2), dtype=float)
        reduced_pss[:, 0] = full_pss[:, 2]
        reduced_pss[:, 1] = full_pss[:, 4]

        tsm = TimeSeriesMetrics(reduced_pss)
        rte_row[j] = tsm.recurrence_time_entropy(
            threshold_mode="std",
            threshold=std_scale,
        )
        info_row[j] = 0

    elapsed = time.perf_counter() - start
    return i, rte_row, info_row, elapsed


def run_for_energy(
    energy, y, py, grid_size, n_cross, std_scale, max_workers, path_to_data
):
    rte = np.full((grid_size, grid_size), np.nan, dtype=float)
    info_grid = np.full((grid_size, grid_size), -999, dtype=int)
    times_y = np.full(grid_size, np.nan, dtype=float)

    datafile = (
        f"{path_to_data}/grid_rte_E={energy:.5f}_N={n_cross}"
        f"_grid_size={grid_size}_scale={std_scale:.3f}.dat"
    )

    print("\n====================================")
    print(f"Energy           : {energy:.5f}")
    print(f"Output data      : {datafile}")
    print("====================================\n")

    total_start = time.perf_counter()

    done = 0
    valid_times = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                compute_rte_row,
                i,
                y_val,
                py,
                energy,
                n_cross,
                std_scale,
            ): i
            for i, y_val in enumerate(y)
        }

        for future in as_completed(futures):
            i, rte_row, info_row, elapsed = future.result()

            rte[i, :] = rte_row
            info_grid[i, :] = info_row
            times_y[i] = elapsed

            if np.any(info_row == 0):
                valid_times.append(elapsed)

            done += 1
            progress_bar(done, grid_size)

    print()
    total_elapsed = time.perf_counter() - total_start

    n_valid = np.count_nonzero(info_grid == 0)
    n_invalid = np.count_nonzero(info_grid != 0)

    valid_times = np.array(valid_times, dtype=float)
    if valid_times.size > 0:
        mean_valid = np.mean(valid_times)
        std_valid = np.std(valid_times)
        total_valid_time = np.sum(valid_times)
    else:
        mean_valid = np.nan
        std_valid = np.nan
        total_valid_time = np.nan

    Y, PY = np.meshgrid(y, py, indexing="ij")

    np.savetxt(
        datafile,
        np.column_stack(
            (
                Y.ravel(),
                PY.ravel(),
                rte.ravel(),
                info_grid.ravel(),
            )
        ),
        fmt=["%.16e", "%.16e", "%.16e", "%d"],
    )

    print("\n=== RUN SUMMARY ===")
    print(f"Energy           : {energy:.5f}")
    print(f"Total ICs        : {grid_size * grid_size}")
    print(f"Valid ICs        : {n_valid}")
    print(f"Invalid ICs      : {n_invalid}")
    print(f"Mean row time    : {mean_valid:.6f} s")
    print(f"Std  row time    : {std_valid:.6f} s")
    print(f"Total row time   : {total_valid_time:.2f} s")
    print(f"Total wall time  : {total_elapsed:.2f} s")
    print(f"Saved data to    : {datafile}")
    print("===================\n")


def main():
    path_to_data = "data"
    os.makedirs(path_to_data, exist_ok=True)

    energies = [1 / 10, 1 / 9, 1 / 8, 1 / 7, 1 / 6]

    n_cross = 5000
    std_scale = 0.05

    y_range = (-0.5, 1.0)
    py_range = (-0.75, 0.75)
    grid_size = 1000

    y = np.linspace(*y_range, grid_size)
    py = np.linspace(*py_range, grid_size)

    max_workers = os.cpu_count() if len(sys.argv) == 1 else int(sys.argv[1])

    print("\n=== GLOBAL PARAMETERS ===")
    print("Model            : Henon-Heiles")
    print(f"Energies         : {[f'{E:.5f}' for E in energies]}")
    print(f"n_cross          : {n_cross}")
    print(f"std_scale        : {std_scale:.3f}")
    print(f"y range          : {y_range}")
    print(f"py range         : {py_range}")
    print(f"grid_size        : {grid_size}")
    print(f"Total ICs / E    : {grid_size * grid_size}")
    print(f"Workers          : {max_workers}")
    print("=========================\n")

    for energy in energies:
        run_for_energy(
            energy,
            y,
            py,
            grid_size,
            n_cross,
            std_scale,
            max_workers,
            path_to_data,
        )


if __name__ == "__main__":
    main()
