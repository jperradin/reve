#include "atom.h"
#include <iostream>

// Constructor
Atom::Atom(std::string elem, double xCoord, double yCoord, double zCoord, double vxComp, double vyComp, double vzComp, int atomID)
    : element(elem), x(xCoord), y(yCoord), z(zCoord), vx(vxComp), vy(vyComp), vz(vzComp), id(atomID) {}

// Getters
std::string Atom::getElement() const {
    return element;
}

double Atom::getX() const {
    return x;
}

double Atom::getY() const {
    return y;
}

double Atom::getZ() const {
    return z;
}

double Atom::getVX() const {
    return vx;
}

double Atom::getVY() const {
    return vy;
}

double Atom::getVZ() const {
    return vz;
}

int Atom::getID() const {
    return id;
}

std::vector<Atom> Atom::getNeighbors() const {
    return neighbors;
}

// Setters
void Atom::setElement(const std::string& elem) {
    element = elem;
}

void Atom::setX(double xCoord) {
    x = xCoord;
}

void Atom::setY(double yCoord) {
    y = yCoord;
}

void Atom::setZ(double zCoord) {
    z = zCoord;
}

void Atom::setVX(double vxComp) {
    vx = vxComp;
}

void Atom::setVY(double vyComp) {
    vy = vyComp;
}

void Atom::setVZ(double vzComp) {
    vz = vzComp;
}

void Atom::setID(int atomID) {
    id = atomID;
}

void Atom::setNeighbors(const std::vector<Atom>& neighborList) {
    neighbors = neighborList;
}

// Add neighbor to the list
void Atom::addNeighbor(const Atom& neighbor) {
    neighbors.push_back(neighbor);
}

// Remove neighbor from the list
void Atom::removeNeighbor(int neighborID) {
    for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
        if (it->getID() == neighborID) {
            neighbors.erase(it);
            break;
        }
    }
}

// Display atom information
void Atom::display() const {
    std::cout << "Element: " << element << "\n"
              << "Coordinates: (" << x << ", " << y << ", " << z << ")\n"
              << "Velocity: (" << vx << ", " << vy << ", " << vz << ")\n"
              << "ID: " << id << std::endl;
}

// Update coordinates
void Atom::updateCoordinates(double xCoord, double yCoord, double zCoord) {
    x = xCoord;
    y = yCoord;
    z = zCoord;
}

// Update velocity
void Atom::updateVelocity(double vxComp, double vyComp, double vzComp) {
    vx = vxComp;
    vy = vyComp;
    vz = vzComp;
}

/*
int main() {
    Atom hydrogen("Hydrogen", 1.23, 4.56, 7.89, 1.0, 1.0, 2.0, 1);
    hydrogen.display();

    return 0;
}
*/