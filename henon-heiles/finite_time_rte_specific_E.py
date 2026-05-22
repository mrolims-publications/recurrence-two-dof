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
    """Write buffer (1D or 2D) to file."""
    if nrows == 0:
        return
    # For 2D arrays, savetxt already uses space delimiter; for 1D, it writes one column.
    np.savetxt(fh, buffer[:nrows], fmt="%.16f", delimiter=" ")
    fh.flush()
    os.fsync(fh.fileno())


def run_single_energy(
    ds,
    energy,
    total_crossing,
    finite_crossing,
    datafile,
    save_pss,
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
    print(f"Save PSS         : {save_pss}")
    print(f"Output data      : {datafile}")
    print("====================================\n")

    if not save_pss:
        # Compact mode: one RTE per window
        buffer = np.empty(chunk_size, dtype=np.float64)
        buffer_count = 0

        with open(datafile, "w", encoding="utf-8") as fh:
            for n in range(num_windows):
                full_pss = ds.poincare_section(q, p, finite_crossing)

                reduced_pss = np.empty((finite_crossing, 2), dtype=np.float64)
                reduced_pss[:, 0] = full_pss[:, 2]
                reduced_pss[:, 1] = full_pss[:, 4]

                tsm = TimeSeriesMetrics(reduced_pss)
                rte = tsm.recurrence_time_entropy(
                    threshold_mode="std",
                    threshold=std_scale,
                )

                buffer[buffer_count] = rte
                buffer_count += 1

                q[:] = full_pss[-1, 1:3]
                p[:] = full_pss[-1, 3:]

                if buffer_count == chunk_size:
                    flush_buffer(fh, buffer, buffer_count)
                    buffer_count = 0

                progress_bar(n + 1, num_windows)

            if buffer_count > 0:
                flush_buffer(fh, buffer, buffer_count)

    else:
        # Full PSS mode: save (x, y, rte) for each point, repeating the same rte for the window
        rows_per_window = finite_crossing
        rows_per_chunk = chunk_size * rows_per_window
        buffer = np.empty((rows_per_chunk, 3), dtype=np.float64)
        buffer_count = 0

        with open(datafile, "w", encoding="utf-8") as fh:
            for n in range(num_windows):
                full_pss = ds.poincare_section(q, p, finite_crossing)

                # Extract x and y coordinates (columns 2 and 4, zero‑based)
                x_coords = full_pss[:, 2]
                y_coords = full_pss[:, 4]

                reduced_pss = np.column_stack((x_coords, y_coords))

                tsm = TimeSeriesMetrics(reduced_pss)
                rte = tsm.recurrence_time_entropy(
                    threshold_mode="std",
                    threshold=std_scale,
                )

                start = buffer_count
                end = start + rows_per_window

                buffer[start:end, 0] = x_coords
                buffer[start:end, 1] = y_coords
                buffer[start:end, 2] = rte  # same rte for all rows in this window

                buffer_count = end

                q[:] = full_pss[-1, 1:3]
                p[:] = full_pss[-1, 3:]

                if buffer_count == rows_per_chunk:
                    flush_buffer(fh, buffer, buffer_count)
                    buffer_count = 0

                progress_bar(n + 1, num_windows)

            if buffer_count > 0:
                flush_buffer(fh, buffer, buffer_count)

    print(f"\nSaved to: {datafile}")


def worker(
    energy, total_crossing, finite_crossing, save_pss, path_to_data, path_to_logs
):
    # Build file names: add "_pss" suffix if saving full PSS
    if save_pss:
        datafile = (
            f"{path_to_data}/finite_time_rte_E={energy:.5f}"
            f"_N={total_crossing}_n={finite_crossing}_pss.dat"
        )
        logfile = (
            f"{path_to_logs}/finite_time_rte_E={energy:.5f}"
            f"_N={total_crossing}_n={finite_crossing}_pss.log"
        )
    else:
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
            run_single_energy(
                ds, energy, total_crossing, finite_crossing, datafile, save_pss
            )


def main():
    path_to_data = "data"
    os.makedirs(path_to_data, exist_ok=True)
    path_to_logs = "logs"
    os.makedirs(path_to_logs, exist_ok=True)

    # User inputs
    try:
        total_crossing = int(input("Enter total number of crossings (e.g., 10000): "))
        finite_crossing = int(
            input("Enter number of crossings per window (e.g., 400): ")
        )
    except ValueError:
        print(
            "Invalid integer. Using defaults: total_crossing=10000, finite_crossing=400"
        )
        total_crossing = 10000
        finite_crossing = 400

    save_pss_input = (
        input("Save full Poincaré section points? (yes/no): ").strip().lower()
    )
    save_pss = save_pss_input in ("yes", "y", "true", "1")

    # Energies to compute (can be adjusted here or prompted as well)
    energies = [1 / 9, 1 / 8, 1 / 7]

    # Prepare arguments for each energy
    args = [
        (e, total_crossing, finite_crossing, save_pss, path_to_data, path_to_logs)
        for e in energies
    ]

    with multiprocessing.Pool() as pool:
        pool.starmap(worker, args)


if __name__ == "__main__":
    main()
