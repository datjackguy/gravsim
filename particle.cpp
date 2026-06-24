#include "core.h"
#include "particle.h"

//Particle Class
Particle::Particle(const OneParticleSettings &settings)
        : mass(settings.mass), pos(settings.position), vel(settings.velocity), colour(settings.colour) {
    }

// Add (or subtract) force from total force acting on particle
void Particle::addForce(Vect3 fnew, int sign) {
    force += sign * fnew;
}

// Calculate the new position
void Particle::calcPos(double dt) {
    pos += vel * dt;
}

// Calculate the new velocity
void Particle::calcVel(double dt) {
    vel += force * dt / mass;
}

// Reset force vector to 0
void Particle::resetForce() {
    force.set(0,0,0);
}
