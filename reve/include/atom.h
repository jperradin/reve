#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>

class Atom {
private:
    std::string element;    // Element symbol
    double x, y, z;         // Coordinates
    double vx, vy, vz;      // Velocity components
    int id;                 // Atom ID in the system
    std::vector<Atom> neighbors;    // Neighboring atoms

public:
    // Constructor
    Atom(std::string elem, double xCoord, double yCoord, double zCoord, double vxComp, double vyComp, double vzComp, int atomID);

    // Getters
    std::string getElement() const;
    double getX() const;
    double getY() const;
    double getZ() const;
    double getVX() const;
    double getVY() const;
    double getVZ() const;
    int getID() const;
    std::vector<Atom> getNeighbors() const;

    // Setters
    void setElement(const std::string& elem);
    void setX(double xCoord);
    void setY(double yCoord);
    void setZ(double zCoord);
    void setVX(double vxComp);
    void setVY(double vyComp);
    void setVZ(double vzComp);
    void setID(int atomID);
    void setNeighbors(const std::vector<Atom>& neighborList);

    // Destructor
    ~Atom() {}

    // Display atom information
    void display() const;

    // Add neighbor to the list
    void addNeighbor(const Atom& neighbor);

    // Remove neighbor from the list
    void removeNeighbor(int neighborID);

    // Update coordinates
    void updateCoordinates(double xCoord, double yCoord, double zCoord);

    // Update velocity
    void updateVelocity(double vxComp, double vyComp, double vzComp);

};

#endif // ATOM_H