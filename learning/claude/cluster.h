#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <fstream>
#include "trajectory.h"

class ClusterAnalyzer {
private:
    const std::vector<Atom>& atoms;
    const Lattice& lattice;
    double cutoff;
    double cutoffSquared;
    
    std::vector<std::vector<int>> adjacencyList;
    std::vector<int> clusterLabels;
    std::vector<int> clusterSizes;
    bool percolates;
    
    // Helper function to calculate distance with periodic boundary conditions
    double calculateDistance(const Atom& a1, const Atom& a2) {
        double dx = a1.x - a2.x;
        double dy = a1.y - a2.y;
        double dz = a1.z - a2.z;
        
        // Apply periodic boundary conditions
        if (dx > lattice.lx/2) dx -= lattice.lx;
        else if (dx < -lattice.lx/2) dx += lattice.lx;
        
        if (dy > lattice.ly/2) dy -= lattice.ly;
        else if (dy < -lattice.ly/2) dy += lattice.ly;
        
        if (dz > lattice.lz/2) dz -= lattice.lz;
        else if (dz < -lattice.lz/2) dz += lattice.lz;
        
        return dx*dx + dy*dy + dz*dz;
    }
    
    // Build adjacency list for atoms within cutoff distance
    void buildAdjacencyList() {
        int numAtoms = atoms.size();
        adjacencyList.resize(numAtoms);
        
        for (int i = 0; i < numAtoms; i++) {
            for (int j = i + 1; j < numAtoms; j++) {
                double distSq = calculateDistance(atoms[i], atoms[j]);
                if (distSq <= cutoffSquared) {
                    adjacencyList[i].push_back(j);
                    adjacencyList[j].push_back(i);
                }
            }
        }
    }
    
    // Breadth-first search to find connected components (clusters)
    void bfs(int startAtom, int clusterLabel) {
        std::queue<int> queue;
        queue.push(startAtom);
        clusterLabels[startAtom] = clusterLabel;
        int clusterSize = 1;
        
        // Track if this cluster spans the simulation box
        bool spansX = false, spansY = false, spansZ = false;
        double minX = atoms[startAtom].x, maxX = atoms[startAtom].x;
        double minY = atoms[startAtom].y, maxY = atoms[startAtom].y;
        double minZ = atoms[startAtom].z, maxZ = atoms[startAtom].z;
        
        while (!queue.empty()) {
            int current = queue.front();
            queue.pop();
            
            for (int neighbor : adjacencyList[current]) {
                if (clusterLabels[neighbor] == -1) {
                    clusterLabels[neighbor] = clusterLabel;
                    queue.push(neighbor);
                    clusterSize++;
                    
                    // Update cluster span
                    minX = std::min(minX, atoms[neighbor].x);
                    maxX = std::max(maxX, atoms[neighbor].x);
                    minY = std::min(minY, atoms[neighbor].y);
                    maxY = std::max(maxY, atoms[neighbor].y);
                    minZ = std::min(minZ, atoms[neighbor].z);
                    maxZ = std::max(maxZ, atoms[neighbor].z);
                }
            }
        }
        
        // Check if cluster spans the simulation box (percolates)
        spansX = (maxX - minX) > 0.9 * lattice.lx;
        spansY = (maxY - minY) > 0.9 * lattice.ly;
        spansZ = (maxZ - minZ) > 0.9 * lattice.lz;
        
        if (spansX || spansY || spansZ) {
            percolates = true;
        }
        
        clusterSizes.push_back(clusterSize);
    }
    
public:
    ClusterAnalyzer(const std::vector<Atom>& a, const Lattice& l, double c) 
        : atoms(a), lattice(l), cutoff(c), cutoffSquared(c*c), percolates(false) {}
    
    void findClusters() {
        // Build adjacency list
        buildAdjacencyList();
        
        // Initialize cluster labels
        int numAtoms = atoms.size();
        clusterLabels.assign(numAtoms, -1);
        clusterSizes.clear();
        percolates = false;
        
        // Find connected components (clusters)
        int clusterLabel = 0;
        for (int i = 0; i < numAtoms; i++) {
            if (clusterLabels[i] == -1) {
                bfs(i, clusterLabel);
                clusterLabel++;
            }
        }
        
        std::cout << "Found " << clusterSizes.size() << " clusters" << std::endl;
    }
    
    void calculatePercolationProperties() {
        // Sort cluster sizes in descending order
        std::sort(clusterSizes.begin(), clusterSizes.end(), std::greater<int>());
        
        std::cout << "Largest cluster size: " << (clusterSizes.empty() ? 0 : clusterSizes[0]) << std::endl;
        std::cout << "System percolates: " << (percolates ? "Yes" : "No") << std::endl;
    }
    
    void writeResults(const std::string& filename) {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Cannot open output file " << filename << std::endl;
            return;
        }
        
        // Write cluster sizes
        outFile << "# Cluster Sizes (in descending order)\n";
        for (size_t i = 0; i < clusterSizes.size(); i++) {
            outFile << i+1 << " " << clusterSizes[i] << "\n";
        }
        
        // Write atom cluster assignments
        outFile << "\n# Atom Cluster Assignments\n";
        outFile << "# AtomID ClusterID\n";
        for (size_t i = 0; i < atoms.size(); i++) {
            outFile << atoms[i].id << " " << clusterLabels[i] << "\n";
        }
        
        outFile.close();
        std::cout << "Results written to " << filename << std::endl;
    }
    
    ClusterStats getClusterStats() {
        ClusterStats stats;
        stats.numClusters = clusterSizes.size();
        stats.largestClusterSize = clusterSizes.empty() ? 0 : clusterSizes[0];
        stats.percolates = percolates;
        stats.clusterSizes = clusterSizes;
        return stats;
    }
};

#endif // CLUSTER_H
