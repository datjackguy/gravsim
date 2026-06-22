#pragma once

#include "core.h"
#include "particle.h"
#include "generation.h"

class Cluster {
private:
    Vect3 centrepos;
    Vect3 vel;
    Vect3 axis;
    double flattening;
    int particleCount;
    double radius;

    double genMass;
    bool rescaleMasses = false;

    std::vector<Particle> particles;
    bool centralMassExists = false;

public:
    Cluster(ClusterSettings settings);
    double enclosedMass(const Particle& particle) const;
    void AddCentralMass(double mass);
    void initialise_particles();
    std::vector<double> scaledMasses(std::vector<double> initialMasses);
    void assignVelocities();
    void reframeParticles();
    std::vector<Particle> getParticleList();
    void configureSoftening(double plotSize);
    double maxSoftening() const;
    void setClusterColour(Color colour);
    double calculateKE() const;
    double calculatePE() const;
    void removeInternalDrift();
    void assignRandomVelocities();
    void rescaleToVirialRatio(double target);
    void assignVirialVelocities(double target);
};