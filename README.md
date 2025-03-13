# Reve package

Reve is a Python package for analyzing trajectories of atomic systems, such as molecules or crystals, in a periodic space.

## Installation

Clone the repository
```bash
git clone https://github.com/jperradin/reve.git
```

Create a virtual environment
```bash
python -m venv venv
source venv/bin/activate
```

Install the package and dependencies
```bash
cd reve
pip install -e .
```

Launch the package that will process the trajectory `tests/input-10_frames.xyz` and calculate and print the direct and periodic boundary condition (PBC) angle between atoms Si$_{22023}$, Si$_{1500}$ and Si$_{22398}$ for each frame 
```bash
python launch_reve.py
```

## Documentation

Coming soon (no)