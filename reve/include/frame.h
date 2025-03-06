#ifndef FRAME_H
#define FRAME_H

#include <vector>
#include <string>
#include "atom.h"

class Frame {
private:
    std::vector<Atom> atoms;        // Atoms in the frame
    int frameID;                    // Frame ID of the trajectory
    double lx, ly, lz;              // Box dimensions
    double volume;                  // Volume of the box
    std::vector<double> lattice;    // Lattice vectors
    int numAtoms;                   // Number of atoms in the frame

public:
    // Constructor
    Frame(int frameID, double boxX, double boxY, double boxZ);

    // Getters
    std::vector<Atom> getAtoms() const;
    int getFrameID() const;
    double getLX() const;
    double getLY() const;
    double getLZ() const;
    double getVolume() const;
    std::vector<double> getLattice() const;
    int getNumAtoms() const;

    // Setters
    void setAtoms(const std::vector<Atom>& atomList);
    void setFrameID(int frameID);
    void setLX(double boxX);
    void setLY(double boxY);
    void setLZ(double boxZ);
    void setVolume(double vol);
    void setLattice(const std::vector<double>& latticeVec);
    void setNumAtoms(int num);

    // Destructor
    ~Frame() {}

    // Display frame information
    void display() const;

    // Add atom to the frame
    void addAtom(const Atom& atom);

    // Remove atom from the frame
    void removeAtom(int atomID);

    // Get atom by element
    std::vector<Atom> getAtomsByElement(const std::string& elem) const;

    // Wrap atoms coordinates to the box
    void wrapAtoms();

    // Find nearest neighbors for each atom within periodic boundary conditions
    void findNeighbors(double cutoff);

    // Calculate the distance between two atoms
    double calculateDistance(const Atom& atom1, const Atom& atom2) const;

};

#endif // FRAME_H