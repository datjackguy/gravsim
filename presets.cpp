#include "presets.h"

void buildSolarSystem(AllSettings &settings, InitialisationQueue &queue) {

    settings.user.plotSize = 5e12;
    settings.user.colourScheme = ColourOptions::None;
    settings.user.enableCentralMass = false;
    //Sun
    OneParticleSettings sun;
    sun.mass = M0;
    sun.position = Vect3{0,0,0};
    sun.velocity = Vect3{0,0,0};
    sun.particleColour = YELLOW;
    queue.addParticle(sun);
    //Earth
    OneParticleSettings earth;
    earth.mass = 5.972e24;
    earth.position = Vect3{1.496e11, 0, 0};
    earth.velocity = Vect3{0, 0, 29.78e3};
    earth.autoOrbit = true;
    // earth.orbitBody = sun;
    earth.axis = Vect3{0,1,0};
    earth.eccentricity = 0;
    earth.particleColour = SKYBLUE;
    queue.addParticle(earth);

}

void buildMilkyWay(AllSettings &settings, InitialisationQueue &queue) {

    settings.user.plotSize = 5e20;
    
    ClusterSettings mway;
    mway.centralMassEnabled = true;
    mway.centralMass = 6.0e10 * M0;
    // mway.centralMass = 5.0e5 * M0;
    mway.clusterGeneratedMass = 5.0e3 * M0;
    mway.axis = {0,1,0};
    mway.flattening = 0.02;
    mway.particleCount = settings.user.particleN;
    mway.IMF = MassDistribution::KroupaIMF;
    mway.rescaleMasses = true;
    mway.position = {0, 0, 0};
    mway.systemVelocity = {0, 0, 0};

    mway.size = 4.7e20;

    queue.addCluster(mway);

}