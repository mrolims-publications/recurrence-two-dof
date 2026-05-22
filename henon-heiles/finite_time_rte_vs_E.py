import contextlib
import multiprocessing
import os
import sys

import numpy as np
from model import px_from_E
from pynamicalsys import HamiltonianSystem, TimeSeriesMetrics


def progress_bar(done, total, width=40):
    frac = done / total
    filled = int(width * frac)
    bar = "#" * filled + "-" * (width - filled)
    print(f"    [{bar}] {done}/{total}", flush=True)


def flush_buffer(fh, buffer, nrows):
    """Write a 1D buffer of RTE values (one per line) to file."""
    if nrows == 0:
        return
    np.savetxt(fh, buffer[:nrows], fmt="%.16f")
    fh.flush()
    os.fsync(fh.fileno())


def run_single_energy(
    ds,
    energy,
    total_crossing,
    finite_crossing,
    datafile,
    chunk_size=1000,
    std_scale=0.05,
):
    q = np.zeros(2, dtype=np.float64)
    p = np.zeros(2, dtype=np.float64)
    q[0] = 0.0
    q[1] = -0.15
    p[1] = 0.0

    info, p[0] = px_from_E(energy, q, p[1])
    if info != 0:
        print("Invalid initial condition. Exiting...")
        sys.exit()

    num_windows = total_crossing // finite_crossing

    print("\n====================================")
    print(f"Energy           : {energy:.5f}")
    print(f"Total crossing   : {total_crossing}")
    print(f"Finite crossing  : {finite_crossing}")
    print(f"Num windows      : {num_windows}")
    print(f"Chunk size       : {chunk_size}")
    print(f"Output data      : {datafile}")
    print("====================================\n")

    # Buffer for RTE values (one per window)
    buffer = np.empty(chunk_size, dtype=np.float64)
    buffer_count = 0

    with open(datafile, "w", encoding="utf-8") as fh:
        for n in range(num_windows):
            # Get one window of Poincaré section points
            full_pss = ds.poincare_section(q, p, finite_crossing)

            # Extract the coordinates (x, y) for RTE calculation
            reduced_pss = np.empty((finite_crossing, 2), dtype=np.float64)
            reduced_pss[:, 0] = full_pss[:, 2]
            reduced_pss[:, 1] = full_pss[:, 4]

            # Compute RTE for this window
            tsm = TimeSeriesMetrics(reduced_pss)
            rte = tsm.recurrence_time_entropy(
                threshold_mode="std",
                threshold=std_scale,
            )

            # Store RTE in buffer
            buffer[buffer_count] = rte
            buffer_count += 1

            # Update state for next window
            q[:] = full_pss[-1, 1:3]
            p[:] = full_pss[-1, 3:]

            # Flush buffer when full
            if buffer_count == chunk_size:
                flush_buffer(fh, buffer, buffer_count)
                buffer_count = 0

            progress_bar(n + 1, num_windows)

        # Flush remaining RTE values
        if buffer_count > 0:
            flush_buffer(fh, buffer, buffer_count)

    print(f"\nSaved to: {datafile}")


def worker(energy, total_crossing, finite_crossing, path_to_data, path_to_logs):
    datafile = (
        f"{path_to_data}/finite_time_rte_E={energy:.5f}"
        f"_N={total_crossing}_n={finite_crossing}.dat"
    )
    logfile = (
        f"{path_to_logs}/finite_time_rte_E={energy:.5f}"
        f"_N={total_crossing}_n={finite_crossing}.log"
    )

    with open(logfile, "w", encoding="utf-8") as log_fh:
        with contextlib.redirect_stdout(log_fh), contextlib.redirect_stderr(log_fh):
            ds = HamiltonianSystem(model="henon heiles")
            run_single_energy(ds, energy, total_crossing, finite_crossing, datafile)


def main():
    path_to_data = "data"
    os.makedirs(path_to_data, exist_ok=True)
    path_to_logs = "logs"
    os.makedirs(path_to_logs, exist_ok=True)

    energies = np.linspace(1 / 10, 1 / 6, 50)

    total_crossing = int(1e10)
    finite_crossing = 400

    args = [
        (e, total_crossing, finite_crossing, path_to_data, path_to_logs)
        for e in energies
    ]

    with multiprocessing.Pool() as pool:
        pool.starmap(worker, args)


if __name__ == "__main__":
    main()
