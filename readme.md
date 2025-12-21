# `Incompressible_Cavity`: Incompressible Flow Solver for Lid-Driven Cavity Flow

This repository contains the source code for Mini-Project 3 for MAE 557: Simulation and Modeling of Fluid Flows. The project consists of a C++ based numerical solver for the two-dimensional, unsteady, incompressible Navier-Stokes equations in a lid-driven cavity configuration.

The solver uses a finite volume method with a pressure-correction scheme and a Forward Euler time integration scheme. Output data is saved in the HDF5 format.

## Requirements

To compile and run this project, you will need the following dependencies installed:

*   **C++ Compiler:** A modern C++ compiler that supports OpenMP (e.g., `g++`).
*   **HDF5 Library (with C++ bindings):** Required for reading and writing simulation data files.
*   **Python 3:** Required for running the analysis and plotting scripts.
*   **Python Libraries:** `numpy`, `matplotlib`, `h5py`, `tqdm`.
    *   These can be installed via pip: `pip install numpy matplotlib h5py tqdm`

## Code Structure

The project is organized into the following directories and files:

```
.
├── data/                     # Recommended directory for simulation output (.h5 files)
├── convergence.py            # Script for spatial and temporal convergence plots
├── plotter.py                # General-purpose plotting script for single simulations
│
└── src/                      # Contains all C++ source code
    ├── main.cpp              # Main C++ driver, sets parameters and runs the simulation loop
    └── grid.hpp              # Header for the SimpleMesh class, data structures, and I/O
```

## How to Compile and Run

The workflow involves configuring simulation parameters directly in the source code, compiling, running the simulation, and then plotting the results.

### 1. Configure Simulation Parameters

Open `src/main.cpp` in a text editor. The primary simulation parameters are defined in the variables at the top of the `main` function. Adjust these values to set up your desired case.

| Variable | Description | Example Value |
| :--- | :--- | :--- |
| `nx1`, `nx2` | The number of grid points in the x and y directions. | `64` |
| `ng` | The number of ghost cells for boundary conditions. | `1` |
| `gamma` | The ratio of specific heats for the fluid. | `1.4` |
| `R` | The specific gas constant. | `287` |
| `CFL` | The Courant-Friedrichs-Lewy number for timestep stability. | `0.5` |
| `printfreq` | How often (in iterations) to print status to the console. | `100` |
| `Re` | The Reynolds number of the flow. | `100` |
| `T0` | The initial and reference temperature in Kelvin. | `300` |
| `Ma` | The Mach number of the flow (used to set `Uw`). | `0.8` |
| `P0` | The initial and reference pressure in Pascals. | `101330` |
| `omegat_max`| The final simulation time in non-dimensional units ($\omega t$). | `10.0` |
| `savedt` | The time interval (in $\omega t$) between saving output files. | `1.0` |
| `tol` | The convergence tolerance for the iterative pressure solver. | `1e-5` |

### 2. Compile the Code

From the project's **root directory**, run the following command:

```bash
g++ -O3 -fopenmp -I/path/to/hdf5/include -L/path/to/hdf5/lib -lhdf5_cpp -lhdf5 src/main.cpp -o fluid
```

**Important:** You must replace `/path/to/hdf5/include` and `/path/to/hdf5/lib` with the actual paths to your HDF5 installation. For example, on macOS with Homebrew:

```bash
g++ -O3 -fopenmp -I/opt/homebrew/include -L/opt/homebrew/lib -lhdf5_cpp -lhdf5 src/main.cpp -o fluid
```

### 3. Run the Simulation

The executable will save `.h5` files in the directory from which it is run. It is highly recommended to organize your runs inside the `data/` directory.

```bash
# Create a directory for your run
mkdir -p data/my_first_run

# Navigate into that directory
cd data/my_first_run

# Run the executable (use ../ to point to the root directory)
../fluid
```

## Analysis and Plotting

The repository includes Python scripts for visualizing results.

### 1. General Plotting (`plotter.py`)

This script generates PNG images for each output file from a single simulation run. It can also generate a multi-panel plot comparing different Mach numbers.

**Usage (Single Plot):**
Run the script from the root directory, pointing it to your data folder.

```bash
python plotter.py -d <directory> -f <field> -n <ngrid> --dt <save_dt>
```

| Argument | Description | Example Value |
| :--- | :--- | :--- |
| `-d` | Directory containing the `.h5` files. | `data/my_first_run` |
| `-f` | The physical field to plot. | `VelMag` |
| `-n` | The number of grid cells (e.g., 64 for 64x64). | `64` |
| `--dt` | The non-dimensional time between saves (matches `savedt`). | `1.0` |

**Available Fields:** `Density`, `VelX1`, `VelX2`, `Press`, `Temp`, `VelMag`, `dpdx`, `dpdy`.

**Usage (Multi-panel Mach Number Plot):**
The script has a special mode to generate a comparison plot across different Mach numbers.

**Important:** This mode contains **hardcoded paths** to the HDF5 files. You must edit the paths inside the `if args.multi:` block in `plotter.py` to point to your specific output files.

To run this mode:
```bash
python plotter.py --multi -f VelMag -n 64
```

### 2. Convergence Plots (`convergence.py`)

This script generates the spatial and temporal convergence plots.

**Important:** This script contains **hardcoded paths** and expects a specific directory structure inside the `data/` folder. You must place your simulation output in folders with these exact names:
*   `data/spat_conv_16/`
*   `data/spat_conv_32/`
*   `data/spat_conv_64/`
*   `data/spat_conv_128/`
*   `data/time_conv_01/`
*   `data/time_conv_005/`
*   `data/time_conv_001/`
*   `data/time_conv_0005/`

**Usage:**
Once your data is organized correctly, run the script from the root directory:
```bash
python convergence.py
```
It will generate `spatial_convergence.png` and `temporal_convergence.png`.
