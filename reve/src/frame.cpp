#include "frame.h"
#include <iostream>
#include <cmath>

// Constructor
Frame::Frame(int frameID, double boxX, double boxY, double boxZ)
    : frameID(frameID), lx(boxX), ly(boxY), lz(boxZ) {}

// Getters
std::vector<Atom> Frame::getAtoms() const {
    return atoms;
}

int Frame::getFrameID() const {
    return frameID;
}

double Frame::getLX() const {
    return lx;
}

double Frame::getLY() const {
    return ly;
}

double Frame::getLZ() const {
    return lz;
}

double Frame::getVolume() const {
    // Assuming box is cubic
    return lx * ly * lz;
}

std::vector<double> Frame::getLattice() const {
    // Lattice vectors are lx, ly, lz
    return {lx, ly, lz};
}

int Frame::getNumAtoms() const {
    return atoms.size();
}


// Setters
void Frame::setAtoms(const std::vector<Atom>& atomList) {
    atoms = atomList;
}

void Frame::setFrameID(int fID) {
    frameID = fID;
}

void Frame::setLX(double boxX) {
    lx = boxX;
}

void Frame::setLY(double boxY) {
    ly = boxY;
}

void Frame::setLZ(double boxZ) {
    lz = boxZ;
}

void Frame::setVolume(double vol) {
    volume = vol;
}

void Frame::setNumAtoms(int num) {
    numAtoms = num;
}

void Frame::setLattice(const std::vector<double>& latticeVec) {
    lx = latticeVec[0];
    ly = latticeVec[1];
    lz = latticeVec[2];
}

// Display frame information
void Frame::display() const {
    std::cout << "Frame ID: " << frameID << "\n"
              << "Box dimensions: (" << lx << ", " << ly << ", " << lz << ")\n"
              << "Volume: " << volume << "\n"
              << "Number of atoms: " << numAtoms << std::endl;
}

// Add atom to the frame
void Frame::addAtom(const Atom& atom) {
    atoms.push_back(atom); // push_back() adds the atom to the end of the vector
}

// Remove atom from the frame
void Frame::removeAtom(int atomID) {
    for (auto it = atoms.begin(); it != atoms.end(); ++it) {
        if (it->getID() == atomID) {
            atoms.erase(it);
            break;
        }
    }
}

// Get atom by element
std::vector<Atom> Frame::getAtomsByElement(const std::string& elem) const {
    std::vector<Atom> atomsByElement;
    for (const auto& atom : atoms) {
        if (atom.getElement() == elem) {
            atomsByElement.push_back(atom);
        }
    }
    return atomsByElement;
}

// Wrap atoms coordinates to the box
void Frame::wrapAtoms() {
    for (auto& atom : atoms) {
        double x = atom.getX();
        double y = atom.getY();
        double z = atom.getZ();

        // Calculate x + lx % lx
        x = fmod(x + lx, lx);
        y = fmod(y + ly, ly);
        z = fmod(z + lz, lz);

        atom.updateCoordinates(x, y, z);
    }
}

// Find nearest neighbors for each atom within periodic boundary conditions
// Using a CKD tree for efficient neighbor search
void Frame::findNeighbors(double cutoff) {
    // Create a CKD tree with the atom coordinates
    // Use the periodic boundary conditions for the box
    // Find the nearest neighbors within the cutoff distance
    // Update the neighbor list for each atom
    return;
}

// Calculate the distance between two atoms
double Frame::calculateDistance(const Atom& atom1, const Atom& atom2) const {
    double x1 = atom1.getX();
    double y1 = atom1.getY();
    double z1 = atom1.getZ();
    double x2 = atom2.getX();
    double y2 = atom2.getY();
    double z2 = atom2.getZ();

    // Calculate the distance between two atoms
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;

    // Apply periodic boundary conditions
    dx = fmod(dx + lx, lx);
    dy = fmod(dy + ly, ly);
    dz = fmod(dz + lz, lz);

    return sqrt(dx * dx + dy * dy + dz * dz);
}

int main() {
    // Create a frame
    Frame frame(1, 10.0, 10.0, 10.0);

    // Create an atom
    Atom atom("O", 5.0, 5.0, 5.0, 0.0, 0.0, 0.0, 1);

    // Add the atom to the frame
    frame.addAtom(atom);

    // Calculate the volume of the box
    double volume = frame.getVolume();

    // Display frame information
    frame.display();

    return 0;
}
// Output

