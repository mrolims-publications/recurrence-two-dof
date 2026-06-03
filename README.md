# recurrence-two-dof

Code repository accompanying the publication entitled "Recurrence in two degrees of freedom Hamiltonian flows" by M. R. Sales, L. C. de Souza, I. L. Caldas, E. D. Leonel, and J. D. Szezech Jr.

This project contains the code to generate and plot all the data from all figures.

Currently, our manuscript is under review. In the meantime you can check out the [preprint].

## Requirements

- Python **3.8+**
- Dependencies listed in `requirements.txt`
- This project uses `pynamicalsys`
- A C compiler with **C99** support
  - GCC
  - Clang
  - Intel oneAPI `icx` (recommended for Linux)
- OpenMP
  - The codes use OpenMP for parallelization.

## Installation

### Python packages

Install the required Python packages with:

```bash
pip install -r requirements.txt
```

### Compiler

#### Linux

Install either GCC, Clang, or Intel oneAPI `icx`.

**GCC (Ubuntu/Debian):**

```bash
sudo apt install build-essential
```

**Intel oneAPI `icx` (recommended):**

See the official installation guide:

https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html

#### macOS (Apple Silicon)

Install OpenMP support:

```bash
brew install libomp
```

## Generating the data and figures

All the data will be stored in the following directories:

- `henon_heiles/data`
- `paradigm_hamiltonian/data`

These directories are created automatically.

All the figures will be stored in the directory `figures/`, which is created automatically by the `Plots.ipynb` notebook.

### Fig. 1

To generate the Poincaré sections shown in Fig. 1, run the script `fig1.py` inside the `henon_heiles` directory:

```bash
cd recurrence-two-dof/henon_heiles
python fig1.py
```

To generate the figure, run all cells under the `Fig. 1` heading inside the `Plots.ipynb` notebook.

### Fig. 2

The data of Fig. 2 is generated within the `Plots.ipynb` notebook under the `Fig. 2` cells.

### Fig. 3

To generate the largest Lyapunov exponent and RTE data shown in Fig. 3, run the scripts `fig3a.py` and `fig3b.py` inside the `henon_heiles` directory:

```bash
cd recurrence-two-dof/henon_heiles
python fig3a.py
python fig3b.py
```

To generate the figure, run all cells under the `Fig. 3` heading inside the `Plots.ipynb` notebook.

### Fig. 4

To generate the largest Lyapunov exponent and RTE data shown in Fig. 4, run the scripts `fig4_1.py` and `fig4_2.py` inside the `henon_heiles` directory:

```bash
cd recurrence-two-dof/henon_heiles
python fig4_1.py
python fig4_2.py
```

To generate the figure, run all cells under the `Fig. 4` heading inside the `Plots.ipynb` notebook.

### Fig. 5

To generate the finite-time RTE data shown in Fig. 5, run the script `fig5.py` inside the `henon_heiles/c/jobs/` directory:

```bash
cd recurrence-two-dof/henon_heiles/c
python jobs/fig5.py
```

It will build, compile, and execute the `henon_heiles/c/programs/finite_time_rte.c` program. By default, the compilation is done using `icx`. To change the compiler, pass it as a command-line argument:

```bash
cd recurrence-two-dof/henon_heiles/c
python jobs/fig5.py gcc
```

To generate the figure, run all cells under the `Fig. 5` heading inside the `Plots.ipynb` notebook.

### Fig. 6

To generate the phase space points shown in Fig. 6, run the script `fig6.py` inside the `henon_heiles/c/jobs/` directory:

```bash
cd recurrence-two-dof/henon_heiles/c
python jobs/fig6.py
```

It will build, compile, and execute the `henon_heiles/c/programs/finite_time_rte.c` program. By default, the compilation is done using `icx`. To change the compiler, pass it as a command-line argument:

```bash
cd recurrence-two-dof/henon_heiles/c
python jobs/fig6.py gcc
```

To generate the figure, run all cells under the `Fig. 6` heading inside the `Plots.ipynb` notebook.

### Fig. 7

The data of Fig. 7 is generated within the `Plots.ipynb` notebook under the `Fig. 7` cells.

### Fig. 8

To generate the phase space points shown in Fig. 8, run the script `fig8.py` inside the `henon_heiles/c/jobs/` directory:

```bash
cd recurrence-two-dof/henon_heiles/c
python jobs/fig8.py
```

It will build, compile, and execute the `henon_heiles/c/programs/finite_time_rte.c` program. By default, the compilation is done using `icx`. To change the compiler, pass it as a command-line argument:

```bash
cd recurrence-two-dof/henon_heiles/c
python jobs/fig8.py gcc
```

> **NOTE**: This will execute 50 programs simultaneously. Either run it on an HPC or remove `&` from the `exe_cmd` string inside the script.

To generate the figure, run all cells under the `Fig. 8` heading inside the `Plots.ipynb` notebook.

### Fig. 9

To generate the proportion of chaotic orbits using RTE and SALI shown in Fig. 9, run the script `fig9.py` inside the `henon_heiles` directory:

```bash
cd recurrence-two-dof/henon_heiles
python fig9.py
```

To generate the figure, run all cells under the `Fig. 9` heading inside the `Plots.ipynb` notebook.

### Fig. 10

To generate the stroboscopic maps and the RTE values shown in Fig. 10, run the scripts `fig10_1.py` and `fig10_2.py` inside the `paradigm_hamiltonian/jobs` directory:

```bash
cd recurrence-two-dof/paradigm_hamiltonian
python jobs/fig10_1.py
python jobs/fig10_2.py
```

To generate the figure, run all cells under the `Fig. 10` heading inside the `Plots.ipynb` notebook.

## Citation

Currently, our manuscript is under review. In the meantime, if you use this repository or parts of it in your work, please consider citing our preprint:

## Contact

For questions, suggestions, or collaboration, please reach out to:

📧 [rolim.sales.m[at]gmail.com](mailto:rolim.sales.m@gmail.com)  
🔗 [mrolims.github.io](https://mrolims.github.io)

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

This project was financed, in part, by the São Paulo Research Foundation (FAPESP, Brazil), under process number 2023/08698-9.

## Disclaimer

As opiniões, hipóteses e conclusões ou recomendações expressas neste material são de responsabilidade do(s) autor(es) e não necessariamente refletem a visão da Fundação de Amparo à Pesquisa do Estado de São Paulo (FAPESP, Brasil).

The opinions, hypotheses, and conclusions or recommendations expressed in this material are the sole responsibility of the author(s) and do not necessarily reflect the views of the São Paulo Research Foundation (FAPESP, Brazil).
