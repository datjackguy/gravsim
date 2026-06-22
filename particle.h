#pragma once

#include "core.h"

//Particle Class
class Particle {
private:
    double mass;
    Vect3 pos;
    Vect3 vel;
    Vect3 force;
    Color colour = WHITE;
    double softening = 0.0;

public:
    Particle(const OneParticleSettings& settings);
    // Add (or subtract) force from total force acting on particle
    void addForce(Vect3 fnew, int sign = 1);
    // Calculate the new position
    void calcPos(double dt);
    // Calculate the new velocity
    void calcVel(double dt);
    // Reset force vector to 0
    void resetForce();

    //GETTERS/SETTERS
    void setColour(Color newColour) {colour = newColour;}
    Color getColour() const {return colour;}
    double getMass() const {return mass;}
    void setMass(double newMass) {mass = newMass;}
    Vect3 getForce() const {return force;}
    Vect3 getPos() const {return pos;}
    void setPos(Vect3 newPos) {pos = newPos;}
    Vect3 getVel() const {return vel;}
    void setVel(Vect3 newVel) {vel = newVel;}
    void setSoftening(double newSoft) {softening = newSoft;}
    double getSoftening() const {return softening;}
};