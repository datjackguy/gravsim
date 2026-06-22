#include "cluster.h"
#include <numeric>
#include <algorithm>

Cluster::Cluster(ClusterSettings settings) :
        centrepos(settings.position), vel(settings.systemVelocity), axis(settings.axis), flattening(settings.flattening), particleCount(settings.particleCount), radius(settings.size), genMass(settings.clusterGeneratedMass), rescaleMasses(settings.rescaleMasses) {
        
        initialise_particles();
        // assignVelocities();
        // reframeParticles();
    
    }

double Cluster::enclosedMass(const Particle& particle) const {
    Vect3 relativePos = particle.getPos();
    double radius = relativePos.magnitude();

    double mass = 0.0;

    for (size_t p = 0; p < particles.size(); p++) {
        Particle other = particles[p];
        Vect3 otherRelativePos = other.getPos();
        double otherRadius = otherRelativePos.magnitude();

        if (otherRadius < radius) {
            mass += other.getMass();
        }
    }

    return mass;
}

void Cluster::AddCentralMass(double mass) {
    if (!centralMassExists) {
        OneParticleSettings cmass;
        cmass.mass = mass;
        cmass.position = {0,0,0};
        cmass.velocity = {0,0,0};
        particles.push_back(Particle(cmass));
        centralMassExists=true;
    }
}

void Cluster::initialise_particles() {
    std::vector<double> masses = generate_masses(particleCount);
    if (rescaleMasses) {
        masses = scaledMasses(masses);
    }

    for (int i = 0; i < particleCount; i++) {
        Vect3 particlepos = generateDisk(radius, flattening, axis);
        Vect3 particlevel{0,0,0};
        OneParticleSettings psettings;
        psettings.mass = masses[i];
        psettings.velocity = particlevel;
        psettings.position = particlepos;
        particles.push_back(Particle(psettings));
    }

}

std::vector<double> Cluster::scaledMasses(std::vector<double> initialMasses) {
    double totalMass = std::accumulate(initialMasses.begin(), initialMasses.end(), 0.0);
    std::vector<double> newMasses;
    for (size_t i = 0; i < initialMasses.size(); i++) {
        double newMass = initialMasses[i]/totalMass*genMass;
        newMasses.push_back(newMass);
    }
    return newMasses;

}

void Cluster::assignVelocities() {
    double clusterSoftening = maxSoftening();

    for (int i = 0; i < particleCount; i++) {
        double encMass = enclosedMass(particles[i]);
        Vect3 currentPos = particles[i].getPos();

        Vect3 particlevel = diskVelocity(currentPos,encMass,axis,clusterSoftening);

        particles[i].setVel(particlevel);
    }
}

void Cluster::reframeParticles() {
    for (size_t i = 0; i < particles.size(); i++) {
        particles[i].setPos(particles[i].getPos()+centrepos);
        particles[i].setVel(particles[i].getVel()+vel);
    }
}

std::vector<Particle> Cluster::getParticleList() {
    return particles;
}

void Cluster::configureSoftening(double plotSize) {

    std::vector<double> masses;
    masses.reserve(particles.size());

    for (const Particle& particle : particles) {
        masses.push_back(particle.getMass());
    }

    std::sort(masses.begin(), masses.end());

    double referenceMass = masses[masses.size() / 2];

    double approximateSpacing = plotSize / std::cbrt(particles.size());
    double baseSoftening = 0.05 * approximateSpacing;
    // double baseSoftening = 0.1 * approximateSpacing;

    double minAllowedSoftening = 0.25 * baseSoftening;
    double maxAllowedSoftening = 5.0 * baseSoftening;

    for (Particle& particle : particles) {
        double massRatio = particle.getMass() / referenceMass;

        double massScale = std::cbrt(massRatio);
        massScale = clampDouble(massScale, 0.5, 10.0);

        double particleSoftening = baseSoftening * massScale;

        particleSoftening = clampDouble(
            particleSoftening,
            minAllowedSoftening,
            maxAllowedSoftening
        );

        particle.setSoftening(particleSoftening);
    }
}
double Cluster::maxSoftening() const {
    double maxSoft = 0.0;

    for (const Particle& particle : particles) {
        maxSoft = std::max(maxSoft, particle.getSoftening());
    }

    return maxSoft;
}
void Cluster::setClusterColour(Color colour) {
    for (Particle& particle : particles) {
        particle.setColour(colour);
    }
}
double Cluster::calculateKE() const {
    double kineticEnergy = 0.0;

    for (const Particle& particle : particles) {
        double speed = particle.getVel().magnitude();
        kineticEnergy += 0.5 * particle.getMass() * speed * speed;
    }

    return kineticEnergy;
}
double Cluster::calculatePE() const {
    double potentialEnergy = 0.0;

    for (size_t i = 0; i < particles.size(); i++) {
        for (size_t j = i + 1; j < particles.size(); j++) {
            Vect3 separation = particles[i].getPos() - particles[j].getPos();
            double r = separation.magnitude();

            double e1 = particles[i].getSoftening();
            double e2 = particles[j].getSoftening();
            double epsilonPair = std::sqrt(e1 * e1 + e2 * e2);

            potentialEnergy += -G * particles[i].getMass() * particles[j].getMass() / std::sqrt(r * r + epsilonPair * epsilonPair);
        }
    }

    return potentialEnergy;
}
void Cluster::removeInternalDrift() {
    Vect3 totalMomentum{0.0, 0.0, 0.0};
    double totalMass = 0.0;

    for (const Particle& particle : particles) {
        totalMomentum += particle.getVel() * particle.getMass();
        totalMass += particle.getMass();
    }

    Vect3 centreVelocity = totalMomentum / totalMass;

    for (Particle& particle : particles) {
        particle.setVel(particle.getVel() - centreVelocity);
    }
}
void Cluster::assignRandomVelocities() {
    for (Particle& particle : particles) {
        Vect3 velocity{randomDouble(-1.0, 1.0), randomDouble(-1.0, 1.0), randomDouble(-1.0, 1.0)};
        particle.setVel(velocity);
    }
    removeInternalDrift();
}
void Cluster::rescaleToVirialRatio(double target) {
    double KE = calculateKE();
    double PE = calculatePE();

    double targetKE = 0.5 * target * std::abs(PE);
    double velocityScale = std::sqrt(targetKE / KE);

    for (Particle& particle : particles) {
        particle.setVel(particle.getVel() * velocityScale);
    }

    removeInternalDrift();
}
void Cluster::assignVirialVelocities(double target) {
    assignRandomVelocities();
    rescaleToVirialRatio(target);
}