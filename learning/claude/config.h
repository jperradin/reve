#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

class Config {
public:
    // Input/output settings
    std::string inputFile;
    std::string outputDir;
    
    // Trajectory settings
    int numAtoms;
    int startFrame;
    int endFrame;
    
    // Analysis settings
    double cutoff;
    bool periodicBoundary;
    int minClusterSize;
    
    // File format settings
    enum FileFormat {
        XYZ,
        EXTENDED_XYZ,
        LAMMPS
    };
    FileFormat fileFormat;
    
    // Default box size for standard XYZ files that don't specify box dimensions
    double defaultBoxSize;
    
    Config() : 
        inputFile("trajectory.xyz"),
        outputDir("output"),
        numAtoms(0),
        startFrame(1),
        endFrame(-1),
        cutoff(3.0),
        periodicBoundary(true),
        minClusterSize(3),
        fileFormat(EXTENDED_XYZ),
        defaultBoxSize(100.0) {}
    
    bool loadFromFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open config file " << filename << std::endl;
            return false;
        }
        
        std::string line, key, value;
        while (std::getline(file, line)) {
            // Skip comments and empty lines
            if (line.empty() || line[0] == '#') continue;
            
            std::istringstream iss(line);
            if (std::getline(iss, key, '=') && std::getline(iss, value)) {
                // Trim whitespace
                key.erase(0, key.find_first_not_of(" \t"));
                key.erase(key.find_last_not_of(" \t") + 1);
                value.erase(0, value.find_first_not_of(" \t"));
                value.erase(value.find_last_not_of(" \t") + 1);
                
                if (key == "input_file") inputFile = value;
                else if (key == "output_dir") outputDir = value;
                else if (key == "num_atoms") numAtoms = std::stoi(value);
                else if (key == "start_frame") startFrame = std::stoi(value);
                else if (key == "end_frame") endFrame = std::stoi(value);
                else if (key == "cutoff") cutoff = std::stod(value);
                else if (key == "periodic_boundary") periodicBoundary = (value == "true" || value == "1");
                else if (key == "min_cluster_size") minClusterSize = std::stoi(value);
                else if (key == "file_format") {
                    if (value == "xyz") fileFormat = XYZ;
                    else if (value == "extended_xyz") fileFormat = EXTENDED_XYZ;
                    else if (value == "lammps") fileFormat = LAMMPS;
                    else std::cerr << "Warning: Unknown file format: " << value << ", using extended_xyz" << std::endl;
                }
                else if (key == "default_box_size") defaultBoxSize = std::stod(value);
                else std::cerr << "Warning: Unknown configuration key: " << key << std::endl;
            }
        }
        
        std::cout << "Configuration loaded from " << filename << std::endl;
        return true;
    }
};

#endif // CONFIG_H
