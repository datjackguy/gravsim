#include "presets.h"

enum SolarBodyIndex {
    Sun,
    Mercury,
    Venus,
    Earth,
    Mars,
};

// Epoch: 2026-06-21
const std::vector<BodyInfo> solarBodies = {
    { "Sun",     M0,     {},       -1, YELLOW },
    { "Mercury", 3.30110000e+23,
      { 5.7909027264243515e+10, 2.0563808222818331e-01, 1.2223228259191675e-01, 8.4295405569269877e-01, 5.0965744980885941e-01, 2.4701705436688464e+00, 0.0 },
      Sun, LIGHTGRAY },
    { "Venus", 4.86750000e+24,
      { 1.0820795156577223e+11, 6.7834459915174037e-03, 5.9242967355261945e-02, 1.3370365667415551e+00, 9.5705678571457831e-01, 1.0600514784568258e+00, 0.0 },
      Sun, ORANGE },
    { "Earth", 5.97220000e+24,
      { 1.4960372927198248e+11, 1.6533868724080230e-02, 6.7068254063374142e-05, 3.5628664897397693e+00, 4.4959740674350233e+00, 2.9344824688587621e+00, 0.0 },
      Sun, SKYBLUE },
    { "Mars", 6.41710000e+23,
      { 2.2792669727367676e+11, 9.3427669731034085e-02, 3.2244745041929798e-02, 8.6360951254386153e-01, 5.0022605015026658e+00, 8.0239174819756187e-01, 0.0 },
      Sun, RED },
    { "Jupiter", 1.89820000e+27,
      { 7.7826911911740234e+11, 4.8167932019439023e-02, 2.2751077610948151e-02, 1.7543537240096381e+00, 4.7734172623650615e+00, 1.8105242084012425e+00, 0.0 },
      Sun, BROWN },
    { "Saturn", 5.68340000e+26,
      { 1.4273837320023918e+12, 5.5383761041869178e-02, 4.3423749218298648e-02, 1.9833712965182662e+00, 5.9102237892424210e+00, 4.9113798319803132e+00, 0.0 },
      Sun, GOLD },
    { "Uranus", 8.68100000e+25,
      { 2.8828371650609946e+12, 4.7313949295762311e-02, 1.3494431918282992e-02, 1.2917141735101787e+00, 1.6072172444098847e+00, 4.5499756615824181e+00, 0.0 },
      Sun, SKYBLUE },
    { "Neptune", 1.02413000e+26,
      { 4.5004049470365693e+12, 1.0382839561587190e-02, 3.0839170232395356e-02, 2.2978858892954173e+00, 4.8929490944692304e+00, 5.4271357296076959e+00, 0.0 },
      Sun, BLUE },
    { "Pluto", 1.30300000e+22,
      { 5.9022717953897432e+12, 2.4637469178738050e-01, 2.9613346118579653e-01, 1.9226155090757482e+00, 1.9927106047001528e+00, 9.3086081514162089e-01, 0.0 },
      Sun, VIOLET },
};






void buildSolarSystem(AllSettings &settings, InitialisationQueue &queue) {

    double targetTime = 0;//randomDouble(0,secondsPerYear);

    settings.user.plotSize = 5e12;
    settings.user.colourScheme = ColourOptions::None;
    settings.user.enableCentralMass = false;

    std::vector<OneParticleSettings> allStates;

    for (BodyInfo body : solarBodies) {
        if (body.ParentBody == -1) {
            OneParticleSettings particleState;
            particleState.mass = body.mass;
            particleState.position = {0,0,0};
            particleState.velocity = {0,0,0};
            particleState.colour = body.colour;
            allStates.push_back(particleState);
        }
        else {
            BodyInfo parent = solarBodies[body.ParentBody];
            ParticleState relative = stateFromKeplerianOrbit(body.orbit,parent.mass,body.mass,targetTime);

            ParticleState absolute;
            absolute.position = allStates[body.ParentBody].position + relative.position;
            absolute.velocity = allStates[body.ParentBody].velocity + relative.velocity;
            OneParticleSettings particleState;
            particleState.mass = body.mass;
            particleState.position = absolute.position;
            particleState.velocity = absolute.velocity;
            particleState.colour = body.colour;

            allStates.push_back(particleState);
        }
    }

    // std::vector<double> asteroid_masses;
    for (int i = 0; i < 50; i++) {
        double mass = randomDouble(1e6,1e14);
        Vect3 pos = generateDisk(2.1/AUpermetre,3.3/AUpermetre, 0.2, {0,1,0});
        Vect3 vel = diskVelocity(pos,M0,{0,1,0});
        OneParticleSettings asteroid;
        asteroid.mass = mass;
        asteroid.colour = DARKGRAY;
        asteroid.position = pos;
        asteroid.velocity = vel;
        allStates.push_back(asteroid);
    }

    for (OneParticleSettings state : allStates) {
        queue.addParticle(state);
    }
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