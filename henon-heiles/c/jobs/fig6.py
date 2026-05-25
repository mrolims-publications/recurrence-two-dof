import os
import sys


def main():
    cc = "icx"
    if len(sys.argv) > 1:
        cc = sys.argv[1]
    target = "finite_time_rte"
    rm_cmd = f"rm -rf bin/{target}.x"
    print(rm_cmd)
    os.system(rm_cmd)
    make_cmd = f"make {target} CC={cc}"
    print(make_cmd)
    os.system(make_cmd)

    if not os.path.isfile(f"bin/{target}.x"):
        print("Executable not found. Did compilation fail?")
        return 1

    energies = [1 / 9, 1 / 8, 1 / 7]
    total_time = int(5e5)
    finite_time = 400
    save_pss = 1

    for energy in energies:
        exe_cmd = (
            f"./bin/{target}.x {energy:.16f} {save_pss} {total_time} {finite_time} &"
        )
        print(exe_cmd)
        os.system(exe_cmd)


if __name__ == "__main__":
    main()
