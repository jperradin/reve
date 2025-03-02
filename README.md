# Reve

Reve is a C++ project designed to calculate the structural properties of glasses. It is optimized using the NVIDIA HPC SDK to ensure high performance and efficiency. This project is intended for researchers and developers who need to analyze the structural properties of glass materials.

## Features

- **Structural Property Calculation**: Reve calculates various structural properties of glasses, providing insights into their behavior and characteristics.
- **Multi-Frame Input Handling**: The project can handle input files containing multiple frames, making it versatile for different types of structural data.
- **Structural Analysis**: Supports structural analysis techniques such as mean squared displacement, allowing for in-depth examination of glass structures.
- **Customizable Settings**: Users can specify the number of atoms and species in the structure, tailoring the analysis to their specific needs.
- **Result Export**: Exports the results of the analysis to a specified directory, ensuring easy access and further processing.

## Getting Started

### Prerequisites

- C++ Compiler (e.g., GCC, Clang)
- [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk)
- CMake (for building the project)

### Installation

Clone the repository:

```bash
git clone https://github.com/jperradin/reve.git
cd reve
```

Build the project using CMake:

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

Prepare your input files with the structural data.
Run the Reve executable with the appropriate command-line arguments:

```bash
./reve --input <input_file> --output <output_directory> --atoms <number_of_atoms> --species <number_of_species>
```

### Example

```bash
./reve --input glass_structure.dat --output results/ --atoms 1000 --species 3
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request if you have any suggestions or improvements.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgments

- Thanks to the NVIDIA HPC SDK for providing the tools to optimize this project.
