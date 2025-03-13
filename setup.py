from setuptools import setup, find_packages
import os
import re

# Read version from __init__.py
with open(os.path.join('reve', '__init__.py'), 'r') as f:
    version_match = re.search(r"__version__ = ['\"]([^'\"]*)['\"]", f.read())
    version = version_match.group(1) if version_match else '0.0.0'

# Read long description from README.md
with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="reve",
    version=version,
    author="jperradin",
    author_email="jperradin@protonmail.com",
    description="Realistic Environment for Vitreous Exploration - A package for working with large-scale molecular dynamics trajectories",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jperradin/reve",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "numba>=0.53.0",
        "matplotlib>=3.4.0",
        "tqdm>=4.60.0",
        "scipy>=1.10.0"
    ],
    entry_points={
        "console_scripts": [
            "reve=reve.main:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)