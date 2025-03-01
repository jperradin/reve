#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <map>
#include "config.h"

struct Atom {
    int id;
    std::string type;
    double x, y, z;
    std::map<std::string, double> properties; // Additional properties from extended XYZ
};

struct Lattice {
    double xlo, xhi, ylo, yhi, zlo, zhi;
    double lx, ly, lz;
    // For non-orthogonal boxes (if needed)
    double xy, xz, yz;
    bool triclinic;
};

struct ClusterStats {
    int frameId;
    int numClusters;
    int largestClusterSize;
    bool percolates;
    std::vector<int> clusterSizes;
};

class TrajectoryAnalyzer {
private:
    Config config;
    std::vector<Atom> atoms;
    Lattice lattice;
    std::vector<ClusterStats> statistics;
    std::ifstream inputFile;
    
    // Parse extended XYZ header information
    bool parseExtendedXYZHeader(const std::string& line, int numAtoms) {
        std::istringstream iss(line);
        std::string key, value;
        
        // Initialize default lattice values
        lattice.triclinic = false;
        lattice.xy = lattice.xz = lattice.yz = 0.0;
        
        while (iss >> key) {
            if (key == "Lattice=") {
                // Parse lattice vectors in the format "Lattice=\"vx1 vx2 vx3 vy1 vy2 vy3 vz1 vz2 vz3\""
                std::string latticeStr;
                if (iss >> latticeStr) {
                    // Remove quotes if present
                    if (latticeStr.front() == '"' && latticeStr.back() == '"') {
                        latticeStr = latticeStr.substr(1, latticeStr.length() - 2);
                    }
                    
                    std::istringstream latticeStream(latticeStr);
                    double vx1, vx2, vx3, vy1, vy2, vy3, vz1, vz2, vz3;
                    if (latticeStream >> vx1 >> vx2 >> vx3 >> vy1 >> vy2 >> vy3 >> vz1 >> vz2 >> vz3) {
                        // Set box bounds based on lattice vectors
                        lattice.xlo = 0.0;
                        lattice.ylo = 0.0;
                        lattice.zlo = 0.0;
                        lattice.xhi = vx1;
                        lattice.yhi = vy2;
                        lattice.zhi = vz3;
                        lattice.lx = vx1;
                        lattice.ly = vy2;
                        lattice.lz = vz3;
                        
                        // Check if box is triclinic
                        if (vx2 != 0.0 || vx3 != 0.0 || vy1 != 0.0 || vy3 != 0.0 || vz1 != 0.0 || vz2 != 0.0) {
                            lattice.triclinic = true;
                            lattice.xy = vy1;
                            lattice.xz = vz1;
                            lattice.yz = vz2;
                        }
                    } else {
                        std::cerr << "Error: Failed to parse lattice vectors" << std::endl;
                        return false;
                    }
                }
            } else if (key == "Properties=") {
                // Parse properties information (will be used when reading atom data)
                // Format: Properties=species:S:1:pos:R:3:charge:R:1
                std::string propertiesStr;
                std::getline(iss, propertiesStr);
                // Store properties information for later use when reading atoms
                // This is just a placeholder - we'll implement the actual parsing when reading atoms
            }
        }
        
        return true;
    }
    
public:
    TrajectoryAnalyzer(const Config& cfg) : config(cfg) {
        inputFile.open(config.inputFile);
        if (!inputFile.is_open()) {
            std::cerr << "Error: Cannot open trajectory file " << config.inputFile << std::endl;
        }
    }
    
    ~TrajectoryAnalyzer() {
        if (inputFile.is_open()) {
            inputFile.close();
        }
    }
    
    bool readFrame(int frameId) {
        // Reset file position to beginning
        inputFile.clear();
        inputFile.seekg(0, std::ios::beg);
        
        std::string line;
        int currentFrame = 0;
        
        // Skip to the requested frame
        while (currentFrame < frameId && inputFile.good()) {
            if (std::getline(inputFile, line)) {
                int numAtoms = std::stoi(line);
                
                // Skip comment/header line and all atom lines for this frame
                std::getline(inputFile, line); // Skip comment/header
                for (int i = 0; i < numAtoms; i++) {
                    std::getline(inputFile, line); // Skip atom data
                }
                
                currentFrame++;
            }
        }
        
        if (currentFrame != frameId || !inputFile.good()) {
            std::cerr << "Error: Could not find frame " << frameId << std::endl;
            return false;
        }
        
        // Read the number of atoms
        std::getline(inputFile, line);
        int numAtoms = std::stoi(line);
        
        // Read the comment/header line which may contain extended XYZ metadata
        std::getline(inputFile, line);
        bool isExtendedXYZ = (line.find("Lattice=") != std::string::npos || 
                              line.find("Properties=") != std::string::npos);
        
        if (isExtendedXYZ) {
            if (!parseExtendedXYZHeader(line, numAtoms)) {
                std::cerr << "Error: Failed to parse extended XYZ header" << std::endl;
                return false;
            }
        } else {
            // For standard XYZ, use default box size or config values
            lattice.xlo = 0.0;
            lattice.ylo = 0.0;
            lattice.zlo = 0.0;
            lattice.xhi = config.defaultBoxSize;
            lattice.yhi = config.defaultBoxSize;
            lattice.zhi = config.defaultBoxSize;
            lattice.lx = config.defaultBoxSize;
            lattice.ly = config.defaultBoxSize;
            lattice.lz = config.defaultBoxSize;
            lattice.triclinic = false;
        }
        
        // Read atom positions
        atoms.clear();
        atoms.resize(numAtoms);
        
        for (int i = 0; i < numAtoms; i++) {
            if (!std::getline(inputFile, line)) {
                std::cerr << "Error: Unexpected end of file while reading atoms" << std::endl;
                return false;
            }
            
            std::istringstream iss(line);
            Atom& atom = atoms[i];
            atom.id = i + 1; // Default ID is the index + 1
            
            if (isExtendedXYZ) {
                // For extended XYZ, the format depends on the Properties field
                // For simplicity, we'll assume a common format: type x y z [additional properties]
                iss >> atom.type >> atom.x >> atom.y >> atom.z;
                
                // Read any additional properties
                double propValue;
                std::string propName;
                int propIndex = 0;
                while (iss >> propValue) {
                    propName = "prop" + std::to_string(propIndex++);
                    atom.properties[propName] = propValue;
                }
            } else {
                // Standard XYZ format: type x y z
                iss >> atom.type >> atom.x >> atom.y >> atom.z;
            }
        }
        
        std::cout << "Read frame " << frameId << " with " << numAtoms << " atoms";
        if (isExtendedXYZ) {
            std::cout << " (Extended XYZ format)";
        }
        std::cout << std::endl;
        
        return true;
    }
    // Support for LAMMPS trajectory format
    bool readLammpsFrame(int frameId) {
        // Reset file position to beginning
        inputFile.clear();
        inputFile.seekg(0, std::ios::beg);
        
        std::string line;
        int currentFrame = 0;
        
        // Skip to the requested frame
        while (currentFrame < frameId && inputFile.good()) {
            if (std::getline(inputFile, line) && line.find("ITEM: TIMESTEP") != std::string::npos) {
                currentFrame++;
                
                // Skip until next timestep or EOF
                while (std::getline(inputFile, line) && line.find("ITEM: TIMESTEP") == std::string::npos) {
                    // Just skip lines
                }
                
                // If we found another timestep, go back one line
                if (line.find("ITEM: TIMESTEP") != std::string::npos) {
                    inputFile.seekg(-static_cast<int>(line.length()) - 1, std::ios::cur);
                }
            }
        }
        
        if (currentFrame != frameId) {
            std::cerr << "Error: Could not find frame " << frameId << std::endl;
            return false;
        }
        
        // Read timestep
        std::getline(inputFile, line); // ITEM: TIMESTEP
        std::getline(inputFile, line); // actual timestep value
        int timestep = std::stoi(line);
        
        // Read number of atoms
        std::getline(inputFile, line); // ITEM: NUMBER OF ATOMS
        std::getline(inputFile, line); // actual number of atoms
        int numAtoms = std::stoi(line);
        
        // Read box bounds
        std::getline(inputFile, line); // ITEM: BOX BOUNDS
        bool triclinic = (line.find("xy xz yz") != std::string::npos);
        
        std::getline(inputFile, line);
        std::istringstream xbounds(line);
        xbounds >> lattice.xlo >> lattice.xhi;
        if (triclinic) xbounds >> lattice.xy;
        
        std::getline(inputFile, line);
        std::istringstream ybounds(line);
        ybounds >> lattice.ylo >> lattice.yhi;
        if (triclinic) ybounds >> lattice.xz;
        
        std::getline(inputFile, line);
        std::istringstream zbounds(line);
        zbounds >> lattice.zlo >> lattice.zhi;
        if (triclinic) zbounds >> lattice.yz;
        
        // Calculate box dimensions
        lattice.lx = lattice.xhi - lattice.xlo;
        lattice.ly = lattice.yhi - lattice.ylo;
        lattice.lz = lattice.zhi - lattice.zlo;
        lattice.triclinic = triclinic;
        
        // Read atom data
        std::getline(inputFile, line); // ITEM: ATOMS
        
        // Parse the atom properties from the header
        std::vector<std::string> atomProperties;
        size_t pos = line.find("ITEM: ATOMS") + 11; // Skip "ITEM: ATOMS"
        std::istringstream propStream(line.substr(pos));
        std::string prop;
        while (propStream >> prop) {
            atomProperties.push_back(prop);
        }
        
        // Read atom data
        atoms.clear();
        atoms.resize(numAtoms);
        
        for (int i = 0; i < numAtoms; i++) {
            if (!std::getline(inputFile, line)) {
                std::cerr << "Error: Unexpected end of file while reading atoms" << std::endl;
                return false;
            }
            
            std::istringstream iss(line);
            Atom& atom = atoms[i];
            
            // Read values based on the properties in the header
            for (size_t j = 0; j < atomProperties.size(); j++) {
                std::string value;
                iss >> value;
                
                if (atomProperties[j] == "id") {
                    atom.id = std::stoi(value);
                } else if (atomProperties[j] == "type") {
                    atom.type = value;
                } else if (atomProperties[j] == "x") {
                    atom.x = std::stod(value);
                } else if (atomProperties[j] == "y") {
                    atom.y = std::stod(value);
                } else if (atomProperties[j] == "z") {
                    atom.z = std::stod(value);
                } else {
                    // Store other properties
                    atom.properties[atomProperties[j]] = std::stod(value);
                }
            }
        }
        
        std::cout << "Read LAMMPS frame " << frameId << " (timestep " << timestep << ") with " 
                  << numAtoms << " atoms" << std::endl;
        return true;
    }
    
    void calculateLatticeProperties() {
        // Additional lattice calculations can be added here if needed
        std::cout << "Lattice dimensions: " << lattice.lx << " x " << lattice.ly << " x " << lattice.lz;
        if (lattice.triclinic) {
            std::cout << " (triclinic box with tilt factors: xy=" << lattice.xy 
                      << ", xz=" << lattice.xz << ", yz=" << lattice.yz << ")";
        }
        std::cout << std::endl;
    }
    
    void updateStatistics(const ClusterStats& stats) {
        statistics.push_back(stats);
    }
    
    void writeSummaryStatistics(const std::string& filename) {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Cannot open output file " << filename << std::endl;
            return;
        }
        
        outFile << "# Frame LargestCluster NumClusters Percolates\n";
        for (const auto& stat : statistics) {
            outFile << stat.frameId << " " 
                    << stat.largestClusterSize << " " 
                    << stat.numClusters << " " 
                    << (stat.percolates ? 1 : 0) << "\n";
        }
        
        outFile.close();
        std::cout << "Summary statistics written to " << filename << std::endl;
    }
    
    // Write the current frame in extended XYZ format
    void writeExtendedXYZ(const std::string& filename) {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Cannot open output file " << filename << std::endl;
            return;
        }
        
        // Write number of atoms
        outFile << atoms.size() << "\n";
        
        // Write header with lattice information
        outFile << "Lattice=\"" << lattice.lx << " 0.0 0.0 ";
        if (lattice.triclinic) {
            outFile << lattice.xy << " " << lattice.ly << " 0.0 ";
            outFile << lattice.xz << " " << lattice.yz << " " << lattice.lz << "\" ";
        } else {
            outFile << "0.0 " << lattice.ly << " 0.0 ";
            outFile << "0.0 0.0 " << lattice.lz << "\" ";
        }
        
        // Add properties information
        outFile << "Properties=species:S:1:pos:R:3";
        
        // Check if we have additional properties to write
        bool hasProperties = false;
        for (const auto& atom : atoms) {
            if (!atom.properties.empty()) {
                hasProperties = true;
                break;
            }
        }
        
        if (hasProperties) {
            // Get property names from the first atom that has properties
            for (const auto& atom : atoms) {
                if (!atom.properties.empty()) {
                    for (const auto& [propName, propValue] : atom.properties) {
                        outFile << ":" << propName << ":R:1";
                    }
                    break;
                }
            }
        }
        outFile << "\n";
        
        // Write atom data
        for (const auto& atom : atoms) {
            outFile << atom.type << " " << atom.x << " " << atom.y << " " << atom.z;
            
            // Write additional properties if any
            for (const auto& [propName, propValue] : atom.properties) {
                outFile << " " << propValue;
            }
            outFile << "\n";
        }
        
        outFile.close();
        std::cout << "Frame written to " << filename << " in extended XYZ format" << std::endl;
    }
    
    const std::vector<Atom>& getAtoms() const {
        return atoms;
    }
    
    const Lattice& getLattice() const {
        return lattice;
    }
};

#endif // TRAJECTORY_H
