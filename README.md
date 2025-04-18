# Reve package

Reve (verre en verlan) is a Python package for analyzing structural properties of trajectories of atomic systems.

## Installation

Clone the repository
```bash
git clone https://github.com/jperradin/reve.git
```

Deactivate the current virtual environment (Facultative)
```bash
deactivate
```
OR
```bash
conda deactivate
```

Create a virtual environment
```bash
cd reve
python -m venv venv
source venv/bin/activate
```

Install the package and dependencies
```bash
pip install -e .
```

Launch the package that will process the trajectory `tests/input-10_frames.xyz` and calculate and print the direct and periodic boundary condition (PBC) angle between atoms Si_22023, Si_1500 and Si_22398 for each frame (this is just a demo)

```bash
python launch_reve.py
```

## Documentation

Coming soon (no)
