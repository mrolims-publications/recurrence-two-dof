import os
import sys


def main():
    cc = "icx"
    if len(sys.argv) > 1:
        cc = sys.argv[1]
    target = "grid_rte"
    rm_cmd = f"rm -rf bin/{target}.x"
    print(rm_cmd)
    os.system(rm_cmd)
    make_cmd = f"make {target} CC={cc} USE_OMP=1"
    print(make_cmd)
    os.system(make_cmd)

    if not os.path.isfile(f"bin/{target}.x"):
        print("Executable not found. Did compilation fail?")
        return 1

    exe_cmd = f"./bin/{target}.x"
    print(exe_cmd)
    os.system(exe_cmd)


if __name__ == "__main__":
    main()
