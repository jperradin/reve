#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <filesystem>
#include <algorithm>
#include "config.h"
#include "trajectory.h"
#include "cluster.h"

namespace fs = std::filesystem;

void createDirectories(const std::string& outputDir) {
    if (!fs::exists(outputDir)) {
        fs::create_directories(outputDir);
        std::cout << "Created output directory: " << outputDir << std::endl;
    }
}

int countConfigurations(const std::string& inputFile, Config::FileFormat format) {
    std::ifstream file(inputFile);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open input file " << inputFile << std::endl;
        return 0;
    }
    
    std::string line;
    int count = 0;
    int numAtoms = 0;
    bool firstFrame = true;
    
    if (format == Config::LAMMPS) {
        // Count LAMMPS frames
        while (std::getline(file, line)) {
            if (line.find("ITEM: TIMESTEP") != std::string::npos) {
                count++;
                
                // If this is the first frame, try to get the number of atoms
                if (firstFrame) {
                    // Skip timestep value
                    std::getline(file, line);
                    
                    // Look for number of atoms
                    while (std::getline(file, line)) {
                        if (line.find("ITEM: NUMBER OF ATOMS") != std::string::npos) {
                            std::getline(file, line);
                            numAtoms = std::stoi(line);
                            firstFrame = false;
                            break;
                        }
                    }
                }
            }
        }
    } else {
        // Count XYZ or Extended XYZ frames
        while (std::getline(file, line)) {
            try {
                int atoms = std::stoi(line);
                if (firstFrame) {
                    numAtoms = atoms;
                    firstFrame = false;
                }
                count++;
                
                // Skip the comment line and atom data
                std::getline(file, line); // Skip comment/header
                for (int i = 0; i < atoms; i++) {
                    std::getline(file, line); // Skip atom data
                }
            } catch (const std::invalid_argument&) {
                // Not a number, skip this line
                continue;
            }
        }
    }
    
    std::cout << "Found " << count << " configurations with " << numAtoms << " atoms" << std::endl;
    return count;
}

int main(int argc, char* argv[]) {
    // Check if configuration file is provided
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file>" << std::endl;
        return 1;
    }
    
    // Initialize configuration
    Config config;
    if (!config.loadFromFile(argv[1])) {
        std::cerr << "Failed to load configuration from " << argv[1] << std::endl;
        return 1;
    }
    
    // Create output directories
    createDirectories(config.outputDir);
    
    // Count configurations in the input file
    int numConfigs = countConfigurations(config.inputFile, config.fileFormat);
    if (numConfigs == 0) {
        std::cerr << "No configurations found in input file" << std::endl;
        return 1;
    }
    
    // Adjust frame range if needed
    if (config.endFrame == -1 || config.endFrame > numConfigs) {
        config.endFrame = numConfigs;
    }
    
    // Initialize trajectory analyzer
    TrajectoryAnalyzer analyzer(config);
    
    // Process frames
    for (int frame = config.startFrame; frame <= config.endFrame; frame++) {
        std::cout << "Processing frame " << frame << " of " << config.endFrame << std::endl;
        
        // Read frame data based on file format
        bool frameRead = false;
        if (config.fileFormat == Config::LAMMPS) {
            frameRead = analyzer.readLammpsFrame(frame);
        } else {
            // Both XYZ and Extended XYZ use the same reader which auto-detects the format
            frameRead = analyzer.readFrame(frame);
        }
        
        if (!frameRead) {
            std::cerr << "Failed to read frame " << frame << std::endl;
            continue;
        }
        
        // Calculate lattice properties
        analyzer.calculateLatticeProperties();
        
        // Find clusters
        ClusterAnalyzer clusterAnalyzer(analyzer.getAtoms(), analyzer.getLattice(), config.cutoff);
        clusterAnalyzer.findClusters();
        
        // Calculate percolation properties
        clusterAnalyzer.calculatePercolationProperties();
        
        // Write results
        std::string frameOutputFile = config.outputDir + "/frame_" + std::to_string(frame) + ".dat";
        clusterAnalyzer.writeResults(frameOutputFile);
        
        // Also write the frame in extended XYZ format with cluster information
        std::string xyzOutputFile = config.outputDir + "/frame_" + std::to_string(frame) + ".xyz";
        analyzer.writeExtendedXYZ(xyzOutputFile);
        
        // Update summary statistics
        ClusterStats stats = clusterAnalyzer.getClusterStats();
        stats.frameId = frame;
        analyzer.updateStatistics(stats);
    }
    
    // Write summary statistics
    analyzer.writeSummaryStatistics(config.outputDir + "/summary.dat");
    
    std::cout << "Analysis complete!" << std::endl;
    return 0;
}
