//Compile Portable -- g++ -O3 gravsim.cpp presets.cpp particle.cpp generation.cpp cluster.cpp orbit.cpp -o gravsim.exe -I include -L lib -lraylib -lgdi32 -lwinmm -static-libgcc -static-libstdc++ -static

//Includes
#include <iostream>
// #include <cmath>
// #include <vector>
#include <random>
#include <thread>
#include <chrono>
#include <algorithm>
#include <deque>
#include <numeric>

#include "raylib.h"
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"
#include "raymath.h"

//Links
#include "core.h"
#include "generation.h"
#include "presets.h"
#include "particle.h"
#include "cluster.h"
#include "orbit.h"



//Utility Functions/Maths

Color randomColour(const std::vector<Color>&colourList) {
    int colouri = randomInteger(0,colourList.size()-1);
    return colourList[colouri];
}



//Convert a 3 value vector to Raylib's Vector3
Vector3 ToVector3(const std::vector<double>& values)
{
    return Vector3{
        static_cast<float>(values[0]),
        static_cast<float>(values[1]),
        static_cast<float>(values[2])
    };
}




Vect3 polarToCartesian(double r, double theta, double phi) {
    Vect3 vector;
    vector.x = r * std::sin(theta)*std::cos(phi);
    vector.y = r * std::sin(theta)*std::sin(phi);
    vector.z = r * std::cos(theta);
    return vector;
}



//Assigns particle colour by distance from origin
Color distanceColour(const Particle& particle, double plotSize) {
    double distance = particle.getPos().magnitude();
    //plotSize/2 <1/3, <2/3, <3/3, else
    if (distance < plotSize/6) {
        return VIOLET;
    }
    else if (distance < plotSize/3) {
        return RED;
    }
    else if (distance < plotSize/2) {
        return ORANGE;
    }
    else {
        return GREEN;
    }
}


//Initial Condition Generation

//Initialise all particles in queue
std::vector<Particle> InitialiseAll(InitialisationQueue queue) {
    std::vector<OneParticleSettings> particleinit = queue.getParticles();
    std::vector<ClusterSettings> clusterinit = queue.getClusters();

    std::vector<Particle> allParticles;

    //Add all inidividual particles
    for (OneParticleSettings p : particleinit) {
        allParticles.push_back(Particle(p));
    }
    //if (p.autoOrbit) {Calculate velocity}

    //Add all clusters
    for (ClusterSettings c : clusterinit) {
        Cluster cluster(c);
        if (c.centralMassEnabled) {
            cluster.AddCentralMass(c.centralMass);
        }
        cluster.setClusterColour(c.clusterColour);
        
        cluster.configureSoftening(c.size);

        if (c.centralMassEnabled) {
            cluster.assignVelocities();
        }
        else {
            cluster.assignVirialVelocities(0.6);
        }

        cluster.reframeParticles();
        std::vector<Particle> clusterParticles = cluster.getParticleList();
        allParticles.insert(allParticles.end(),clusterParticles.begin(),clusterParticles.end());

    }

    return allParticles;
}

std::vector<int> GenerateClusterSizes(int clusters, int particles) {
    int average = particles / clusters;

    int minSize = std::max(1, average / 2);
    int maxSize = average * 2;

    std::vector<int> sizes(clusters, minSize);

    int remaining = particles - minSize * clusters;

    int maxChunk = std::max(1, average / 2);

    while (remaining > 0) {
        int i = randomInteger(0, clusters - 1);

        int spaceLeft = maxSize - sizes[i];

        if (spaceLeft <= 0) {
            continue;
        }

        int chunkLimit = std::min({ maxChunk, remaining, spaceLeft });

        int amount = randomInteger(1, chunkLimit);

        sizes[i] += amount;
        remaining -= amount;
    }

    return sizes;
}


void autoClusters(const AllSettings &settings, InitialisationQueue &queue) {
    if (settings.init.clusterCount == 1) {
        ClusterSettings clustersettings;
        clustersettings.axis = Vect3{0,1,0};
        clustersettings.position = Vect3{0,0,0};
        clustersettings.systemVelocity = Vect3{0,0,0};
        clustersettings.particleCount = settings.user.particleN;
        clustersettings.size = settings.user.plotSize/2;
        clustersettings.flattening = 0.5;
        clustersettings.rescaleMasses = false;
        clustersettings.centralMassEnabled = settings.user.enableCentralMass;
        clustersettings.centralMass = 1e34;

        queue.addCluster(clustersettings);
    }
    else {
        std::vector<int> split = GenerateClusterSizes(settings.init.clusterCount,settings.user.particleN);
        for (int c = 0; c < settings.init.clusterCount; c++) {
            ClusterSettings clustersettings;
            clustersettings.position = 1.5*randomSph(settings.user.plotSize);
            clustersettings.systemVelocity = {0,0,0};
            Vect3 sysVel = randomSign()*0.3*globularVel(clustersettings.position,1e33);

            // clustersettings.position = 2.5 * randomSph(settings.user.plotSize);
            clustersettings.systemVelocity = sysVel;
            clustersettings.axis = randomAxis();
            clustersettings.flattening = randomDouble(0,1);
            if (settings.user.enableCentralMass) {
                clustersettings.size = settings.user.plotSize/4;
            }
            else {
                clustersettings.size = settings.user.plotSize/10;
                clustersettings.flattening = 1;
            }
            // clustersettings.size = settings.user.plotSize/10;
            clustersettings.particleCount = split[c];
            clustersettings.rescaleMasses = false;
            clustersettings.centralMassEnabled = settings.user.enableCentralMass;
            clustersettings.centralMass = 1e33;

            if (settings.user.colourScheme == ColourOptions::Cluster) {
                clustersettings.clusterColour=colourList[c%colourList.size()];
            }

            queue.addCluster(clustersettings);
        }
    }   
}


//Simulation Class - Should be self-contained (as much as possible)
class Simulation {
private:
    AllSettings settings;

    std::vector<Particle> particles;
    // std::vector<std::vector<int>> particlePairs;

    double timestep = 0.5*10*24*60*60;
    double epsilon = 5e13;

    double totalMass = 0.0;
    double maxMass = 0.0; //Initialise maximum present mass, determined after particles created
    double minMass = 1e64; //As above for minimum mass
    double currentTime = 0.0;

    bool forcesExist = false;

    //Timing
    double lastUpdateMs = 0.0;
    double averageUpdateMs = 0.0;
    bool updateTimingExists = false;

public:
    Simulation(const AllSettings& settings, std::vector<Particle>& particles)
        : settings(settings), particles(particles) {
        // initialise_clusters();
        // initialise_particles();
        updateMassLimits();
        momentum_centre();

        // configureSoftening();
        configureTimestep();

        // create_particle_pairs();
        colour_options();

    }
private:
    void resetForces() {
        for (Particle& particle : particles) {
            particle.resetForce();
        }
    }
    double medianMass() const {
        std::vector<double> masses;
        masses.reserve(particles.size());

        for (const Particle& particle : particles)
        {
            masses.push_back(particle.getMass());
        }

        std::sort(masses.begin(), masses.end());

        return masses[masses.size() / 2];
    }
    void configureSoftening() {
        // double particleCount = static_cast<double>(particles.size());

        double approximateSpacing = settings.user.plotSize / std::cbrt(particles.size());

        // double baseSoftening = 0.05 * approximateSpacing;
        double baseSoftening = 0.3 * approximateSpacing;

        double referenceMass = medianMass();

        double minAllowedSoftening = 0.25 * baseSoftening;
        double maxAllowedSoftening = 5.0 * baseSoftening;

        for (Particle& particle : particles) {
            double massRatio = particle.getMass() / referenceMass;

            double massScale = std::cbrt(massRatio);

            massScale = clampDouble(massScale, 0.5, 10.0);

            double particleSoftening = baseSoftening * massScale;

            particleSoftening = clampDouble(particleSoftening, minAllowedSoftening, maxAllowedSoftening);

            particle.setSoftening(particleSoftening);
        }
    }
    void configureTimestep() {
        double systemRadius = 0.5 * settings.user.plotSize;

        double dynamicalTime = std::sqrt((systemRadius * systemRadius * systemRadius) / (G * totalMass));

        timestep = dynamicalTime / 3000.0;
        // timestep = dynamicalTime / 500.0;
    }
    //Colour options
    void colour_options() {
        switch (settings.user.colourScheme) {
            case (ColourOptions::White):
                break;
            case (ColourOptions::Random):
                for (Particle& particle : particles) {
                    particle.setColour(randomColour(colourList));
                }
                break;
            case (ColourOptions::Distance):
                for (Particle& particle : particles) {
                    particle.setColour(distanceColour(particle, settings.user.plotSize));
                }
                break;
            default:
                break;
        } 

        //Set central mass to yellow in single cluster mode
        if (particles.empty() == false && settings.user.enableCentralMass && settings.init.clusterCount == 1) {
            particles.back().setColour(YELLOW);
        }
    }

    void updateMassLimits() {
        totalMass = 0.0;
        minMass = 1e64;
        maxMass = 0.0;
        for (size_t i = 0; i < particles.size(); i++) {
            double particleMass = particles[i].getMass();
            totalMass += particleMass;
            if (particleMass > maxMass) {
                maxMass = particleMass;
            }
            if (particleMass < minMass) {
                minMass = particleMass;
            }
        }
    }
    //Recentre system to avoid drift
    void momentum_centre() {
        Vect3 weightedPosition{0.0, 0.0, 0.0};
        Vect3 totalMomentum{0.0, 0.0, 0.0};
        double massSum = 0.0;

        for (const Particle& particle : particles) {
            double mass = particle.getMass();

            weightedPosition += particle.getPos() * mass;
            totalMomentum += particle.getVel() * mass;
            massSum += mass;
        }

        Vect3 centreOfMass = weightedPosition / massSum;
        Vect3 centreVelocity = totalMomentum / massSum;

        for (Particle& particle : particles) {
            particle.setPos(particle.getPos() - centreOfMass);
            particle.setVel(particle.getVel() - centreVelocity);
        }
    }
    //Calculate total mass within the radius of a particle
    double enclosedMass(const Particle &particle) const {
        Vect3 ppos = particle.getPos();

        double radius = pythagoras(ppos);
        double mass = 0;
        for (size_t p = 0; p < particles.size(); p++) {
            Vect3 currentPos = particles[p].getPos();

            double currentRadius = pythagoras(currentPos);
            if (currentRadius < radius) {
                mass += particles[p].getMass();
            }
        }
        return mass;
    }
    double distanceCalc(const Particle &p1, const Particle &p2) const {
        Vect3 pos1 = p1.getPos();
        Vect3 pos2 = p2.getPos();

        Vect3 disp = pos1-pos2;

        double distance = disp.magnitude();

        return distance;
    }

    // Calculate the force between two particles
    Vect3 forceCalc(const Particle &p1, const Particle &p2) {
        double distance = distanceCalc(p1,p2);
        double mass1 = p1.getMass();
        double mass2 = p2.getMass();
        Vect3 pos1 = p1.getPos();
        Vect3 pos2 = p2.getPos();

        Vect3 force;

        double e1 = p1.getSoftening();
        double e2 = p2.getSoftening();

        double epsilonPair = sqrt(e1*e1+e2*e2);

        force = -G*mass1*mass2*(pos1-pos2)/pow((distance*distance)+(epsilonPair*epsilonPair),1.5);

        return force;
    }

    void calculateForces() {
        resetForces();

        for (size_t i = 0; i < particles.size(); i++) {
            for (size_t j = i+1; j < particles.size(); j++) {
                Vect3 force = forceCalc(particles[i],particles[j]);
                particles[i].addForce(force);
                particles[j].addForce(force, -1);
            }
        }
    }

public:
    //Update the simulation for each timestep
    void update() {

        const auto updateStart = std::chrono::steady_clock::now(); //Update timer

        if (forcesExist == false) {
            calculateForces();
            forcesExist = true;
        }

        for (size_t i = 0; i < particles.size(); i++) {
            particles[i].calcVel(timestep/2);
        }

        for (size_t i = 0; i < particles.size(); i++) {
            particles[i].calcPos(timestep);
        }

        calculateForces();

        for (size_t i = 0; i < particles.size(); i++) {
            particles[i].calcVel(timestep/2);
        }

        currentTime+=timestep;

        //Timing
        const auto updateEnd = std::chrono::steady_clock::now();

        lastUpdateMs = std::chrono::duration<double, std::milli>(updateEnd - updateStart).count();

        if (!updateTimingExists) {
            averageUpdateMs = lastUpdateMs;
            updateTimingExists = true;
        }
        else {
            averageUpdateMs = 0.95 * averageUpdateMs + 0.05 * lastUpdateMs;
        }
    }

    const std::vector<Particle>& getParticles() const {
        return particles;
    }

    double getPlotSize() const {
        return settings.user.plotSize;
    }

    double getMaxMass() const {
        return maxMass;
    }

    double getMinMass() const {
        return minMass;
    }

    double getTotalMass() const {
        return totalMass;
    }

    bool shouldDrawTrails() const {
        return settings.user.drawTrails;
    }

    double getCurrentTimeYears() const {
        return currentTime / secondsPerYear;
    }

    double getLastUpdateMs() const {
        return lastUpdateMs;
    }

    double getAverageUpdateMs() const {
        return averageUpdateMs;
    }
    double getTimestepYears() const {
        return timestep / secondsPerYear;
    }
};

struct EnergyDiagnosticsResult {
    double kineticEnergy = 0.0;
    double potentialEnergy = 0.0;
    int unboundParticles = 0;
};



class PhysicsDiagnostics {
private:
    const std::vector<Particle>& particles;

    double initialKE = 0.0;
    double initialPE = 0.0;
    double initialTotalEnergy = 0.0;

    double currentKE = 0.0;
    double currentPE = 0.0;
    double currentTotalEnergy = 0.0;

    std::deque<double> KEhistory;
    std::deque<double> PEhistory;
    std::deque<double> energyHistory;

    int maxHistory = 300;

    int currentUnboundParticles = 0;
    double currentUnboundFraction = 0.0;

    std::deque<double> unboundHistory;

public:
    PhysicsDiagnostics(const std::vector<Particle>& particles, int maxHistory = 300)
        : particles(particles), maxHistory(maxHistory) {
        initialise();
    }

    void initialise() {
        EnergyDiagnosticsResult result = calculateEnergyDiagnostics();

        initialKE = result.kineticEnergy;
        initialPE = result.potentialEnergy;
        initialTotalEnergy = initialKE + initialPE;

        currentKE = initialKE;
        currentPE = initialPE;
        currentTotalEnergy = initialTotalEnergy;

        currentUnboundParticles = result.unboundParticles;

        if (!particles.empty()) {
            currentUnboundFraction = static_cast<double>(currentUnboundParticles) / particles.size();
        }
        else {
            currentUnboundFraction = 0.0;
        }

        recordCurrentValues();
    }

    double calculateKE() const {
        double kineticEnergy = 0.0;

        for (const Particle& particle : particles) {
            double speed = particle.getVel().magnitude();
            kineticEnergy += 0.5 * particle.getMass() * speed * speed;
        }

        return kineticEnergy;
    }

    double calculatePE() const {
        double potentialEnergy = 0.0;

        for (size_t i = 0; i < particles.size(); i++) {
            for (size_t j = i + 1; j < particles.size(); j++) {
                Vect3 separation = particles[i].getPos() - particles[j].getPos();
                double r = separation.magnitude();

                double soft1 = particles[i].getSoftening();
                double soft2 = particles[j].getSoftening();
                double pairSoftening = std::sqrt(soft1 * soft1 + soft2 * soft2);

                potentialEnergy += -G * particles[i].getMass() * particles[j].getMass() / std::sqrt(r * r + pairSoftening * pairSoftening);
            }
        }

        return potentialEnergy;
    }

    void sample() {
        EnergyDiagnosticsResult result = calculateEnergyDiagnostics();

        currentKE = result.kineticEnergy;
        currentPE = result.potentialEnergy;
        currentTotalEnergy = currentKE + currentPE;

        currentUnboundParticles = result.unboundParticles;

        if (!particles.empty()) {
            currentUnboundFraction = static_cast<double>(currentUnboundParticles) / particles.size();
        }
        else {
            currentUnboundFraction = 0.0;
        }

        recordCurrentValues();
    }

private:
    Vect3 calculateCentreVelocity() const {
        Vect3 totalMomentum{0.0, 0.0, 0.0};
        double totalMass = 0.0;

        for (const Particle& particle : particles) {
            double mass = particle.getMass();

            totalMomentum += particle.getVel() * mass;
            totalMass += mass;
        }

        return totalMomentum / totalMass;
    }

    EnergyDiagnosticsResult calculateEnergyDiagnostics() const {
        EnergyDiagnosticsResult result;

        std::vector<double> specificPotential(particles.size(), 0.0);

        for (size_t i = 0; i < particles.size(); i++) {
            for (size_t j = i + 1; j < particles.size(); j++) {
                Vect3 separation = particles[i].getPos() - particles[j].getPos();
                double r = separation.magnitude();

                double soft1 = particles[i].getSoftening();
                double soft2 = particles[j].getSoftening();
                double pairSoftening = std::sqrt(soft1 * soft1 + soft2 * soft2);

                double softenedDistance = std::sqrt(r * r + pairSoftening * pairSoftening);

                double massI = particles[i].getMass();
                double massJ = particles[j].getMass();

                double pairPotentialEnergy = -G * massI * massJ / softenedDistance;

                result.potentialEnergy += pairPotentialEnergy;

                specificPotential[i] += -G * massJ / softenedDistance;
                specificPotential[j] += -G * massI / softenedDistance;
            }
        }

        Vect3 centreVelocity = calculateCentreVelocity();

        for (size_t i = 0; i < particles.size(); i++) {
            Vect3 relativeVelocity = particles[i].getVel() - centreVelocity;
            double speed = relativeVelocity.magnitude();

            double mass = particles[i].getMass();

            double specificKineticEnergy = 0.5 * speed * speed;
            double specificTotalEnergy = specificKineticEnergy + specificPotential[i];

            result.kineticEnergy += mass * specificKineticEnergy;

            if (specificTotalEnergy > 0.0) {
                result.unboundParticles++;
            }
        }

        return result;
    }

    void recordCurrentValues() {
        KEhistory.push_back(currentKE);
        PEhistory.push_back(currentPE);
        energyHistory.push_back(currentTotalEnergy);
        unboundHistory.push_back(static_cast<double>(currentUnboundParticles));

        trimHistory(KEhistory);
        trimHistory(PEhistory);
        trimHistory(energyHistory);
        trimHistory(unboundHistory);
    }

    void trimHistory(std::deque<double>& history) {
        while (history.size() > maxHistory) {
            history.pop_front();
        }
    }

public:
    int getCurrentUnboundParticles() const {
        return currentUnboundParticles;
    }

    double getCurrentUnboundFraction() const {
        return currentUnboundFraction;
    }

    const std::deque<double>& getUnboundHistory() const {
        return unboundHistory;
    }

    double getInitialKE() const {
        return initialKE;
    }

    double getInitialPE() const {
        return initialPE;
    }

    double getInitialEnergy() const {
        return initialTotalEnergy;
    }

    double getCurrentKE() const {
        return currentKE;
    }

    double getCurrentPE() const {
        return currentPE;
    }

    double getCurrentEnergy() const {
        return currentTotalEnergy;
    }

    const std::deque<double>& getKEhistory() const {
        return KEhistory;
    }

    const std::deque<double>& getPEhistory() const {
        return PEhistory;
    }

    const std::deque<double>& getEnergyHistory() const {
        return energyHistory;
    }

    double getEnergyErrorPercent() const {
        if (initialTotalEnergy == 0.0) {
            return 0.0;
        }

        return 100.0 * (currentTotalEnergy - initialTotalEnergy) / std::abs(initialTotalEnergy);
    }
};


// Drawing/Rendering Functions
Vector3 ToRenderPosition(Vect3 pos, double plotSize) {
    float renderScale = 10.0f / static_cast<float>(plotSize);

    return Vector3{static_cast<float>(pos.x * renderScale), static_cast<float>(pos.y * renderScale), static_cast<float>(pos.z * renderScale)};
}

class TrailManager {
private:
    std::vector<std::deque<Vector3>> allTrails;
    int maxTrailLength = 0;
    double plotSize = 1.0;

public:
    TrailManager(int particleCount, int maxTrailLength, double plotSize)
        : allTrails(particleCount), maxTrailLength(maxTrailLength), plotSize(plotSize) {
    }

    void record(const std::vector<Particle>& particles) {
        if (allTrails.size() != particles.size()) {
            allTrails.clear();
            allTrails.resize(particles.size());
        }

        for (size_t i = 0; i < particles.size(); i++) {
            Vector3 trailPoint = ToRenderPosition(particles[i].getPos(),plotSize);

            allTrails[i].push_back(trailPoint);

            while (static_cast<int>(allTrails[i].size()) > maxTrailLength) {
                allTrails[i].pop_front();
            }
        }
    }

    void draw(const std::vector<Particle>& particles) const {
        for (size_t p = 0; p < allTrails.size(); p++) {
            Color colour = particles[p].getColour();

            for (size_t i = 1; i < allTrails[p].size(); i++) {
                DrawLine3D(allTrails[p][i - 1], allTrails[p][i], colour);
            }
        }
    }

    void clear() {
        for (auto& trail : allTrails) {
            trail.clear();
        }
    }
};


//Enable user control of camera view (zoom and rotate)
void UpdateOrbitCamera(Camera3D& camera, float& yaw, float& pitch, float& distance) {
    Vector2 mouseDelta = GetMouseDelta();

    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT))
    {
        yaw -= mouseDelta.x * 0.01f;
        pitch += mouseDelta.y * 0.01f;
    }

    // Limit the pitch so the camera does not flip upside down.
    const float pitchLimit = 1.5f;

    if (pitch > pitchLimit) {
        pitch = pitchLimit;
    }

    if (pitch < -pitchLimit) {
        pitch = -pitchLimit;
    }

    distance -= GetMouseWheelMove() * 1.0f;

    if (distance < 0.5f) {
        distance = 0.5f;
    }

    if (distance > 100.0f) {
        distance = 100.0f;
    }

    // Move Camera Target
    float moveSpeed = 5.0f * GetFrameTime();

    Vector3 forward = {-sinf(yaw), 0.0f, -cosf(yaw)};

    Vector3 right = {cosf(yaw), 0.0f, -sinf(yaw)};

    if (IsKeyDown(KEY_UP)) {
        camera.target.x += forward.x * moveSpeed;
        camera.target.z += forward.z * moveSpeed;
    }

    if (IsKeyDown(KEY_DOWN)) {
        camera.target.x -= forward.x * moveSpeed;
        camera.target.z -= forward.z * moveSpeed;
    }

    if (IsKeyDown(KEY_RIGHT)) {
        camera.target.x += right.x * moveSpeed;
        camera.target.z += right.z * moveSpeed;
    }

    if (IsKeyDown(KEY_LEFT)) {
        camera.target.x -= right.x * moveSpeed;
        camera.target.z -= right.z * moveSpeed;
    }

    if (IsKeyPressed(KEY_O)) {
        camera.target.x = 0.0f;
        camera.target.z = 0.0f;
    }

    // Convert yaw, pitch and distance into a 3D camera position.
    camera.position.x = camera.target.x + distance * cosf(pitch) * sinf(yaw);
    camera.position.y = camera.target.y + distance * sinf(pitch);
    camera.position.z = camera.target.z + distance * cosf(pitch) * cosf(yaw);
}

//Toggle Info
bool checkHideInfo(bool hideInfo) {
    if (IsKeyPressed(KEY_F3)) {
        return !hideInfo;
    }
    else {
        return hideInfo;
    }
}

//Toggle Diagnostics
bool checkHideDiagnostics(bool hide) {
    if (IsKeyPressed(KEY_F5)) {
        return !hide;
    }
    else {
        return hide;
    }
}

//Toggle Axes
bool checkDrawAxes(bool show) {
    if (IsKeyPressed(KEY_F6)) {
        return !show;
    }
    else {
        return show;
    }
}

// Calculate animation sphere size
float GetVisualRadius(double mass, double totalMass, double minMass, double maxMass, int numParticles) {

    //Change to median!!!?

    float averageMass = (totalMass/numParticles);

    float minRadius = 0.01f;
    float maxRadius = 0.4f;

    double logMin = std::log10(minMass);
    double logMax = std::log10(maxMass);
    double logMass = std::log10(mass);

    float t = static_cast<float>((logMass - logMin) / (logMax - logMin));

    t = Clamp(t, 0.0f, 1.0f);

    t = std::pow(t, 1.4f); // Exponent = contrast better small/large mass

    // float radius = 0.1f + 0.05f * log10(mass/averageMass) / log10(maxMass/averageMass);

    // return Clamp(radius, 0.01f, 0.5f);
    return minRadius + t * (maxRadius-minRadius);
}
// Modified
float GetVisualRadius(double mass) {
    float minRadius = 0.01f;
    float maxRadius = 0.4f;

    double visualMinMass = 1e29;
    double visualMaxMass = 1e33;

    double logMin = std::log10(visualMinMass);
    double logMax = std::log10(visualMaxMass);
    double logMass = std::log10(mass);

    float t = static_cast<float>((logMass - logMin) / (logMax - logMin));

    t = Clamp(t, 0.0f, 1.0f);
    t = std::pow(t, 1.4f);

    return minRadius + t * (maxRadius - minRadius);
}

// void DrawTrail(const Particle &particle) {
//     for (size_t i = 1; i < particle.trail.size(); i++) {
//         DrawLine3D(particle.trail[i-1], particle.trail[i], particle.getColour());
//     }
// }

//Live sim controls
//Allow the user to move a particle mid-simulation with WASD 
void UpdateParticlePos(Particle &particle, UserSettings settings) {
    Vect3 ppos = particle.getPos();
    double x = ppos.x;
    double y = ppos.y;
    double z = ppos.z;
    bool changed = false;
    if (IsKeyDown(KEY_S)) {
        x = x + settings.plotSize/120;
        changed = true;
    }
    if (IsKeyDown(KEY_W)) {
        x = x - settings.plotSize/120;
        changed = true;
    }
    if (IsKeyDown(KEY_D)) {
        z = z - settings.plotSize/120;
        changed = true;
    }
    if (IsKeyDown(KEY_A)) {
        z = z + settings.plotSize/120;
        changed = true;
    }

    if (changed) {
        Vect3 newPos;
        newPos.set(x,y,z);
        particle.setPos(newPos);
    }

}

//Allow the user the change the mass of a particle with R/F
void UpdateParticleMass(Particle &particle) {
    double mass = particle.getMass();

    if (IsKeyDown(KEY_R)) {
        mass *= 1.05;
    }
    if (IsKeyDown(KEY_F)) {
        mass *= 0.95;
    }

    particle.setMass(mass);
}

//Change Sim Speed
double UpdateSimSpeed(double targetFramesPerSecond) {
    if (IsKeyPressed(KEY_EQUAL) || (IsKeyPressedRepeat(KEY_EQUAL) && IsKeyDown(KEY_EQUAL))) {
        targetFramesPerSecond += 10;
    };
    // if (IsKeyPressed(KEY_MINUS) || (IsKeyPressedRepeat(KEY_MINUS) && IsKeyDown(KEY_MINUS)) && targetFramesPerSecond > 0) {
    if (IsKeyPressed(KEY_MINUS) && targetFramesPerSecond > 0) {
        targetFramesPerSecond -= 10;
    }
    return targetFramesPerSecond;
}


//GUI Control Option Classes
class IntSlider {
private:
    std::string name;
    std::string lowLimitStr;
    std::string highLimitStr;
    std::string sliderText;
    float lowLimit;
    float highLimit;
    float sliderValue;
    // int pos;
    int* outputValue;
public:
    int height = 60;

    IntSlider(int &value, std::string lowLimitStr, std::string highLimitStr, float lowLimit, float highLimit, std::string sliderText) :
        outputValue(&value), lowLimitStr(lowLimitStr), highLimitStr(highLimitStr), lowLimit(lowLimit), highLimit(highLimit), sliderText(sliderText)  {

    }

    void draw(float yposition) {
        float sliderValue = static_cast<float>(*outputValue);

        //Draw slider - set position and limits
        GuiSlider(Rectangle{40, yposition, 300, 20}, lowLimitStr.c_str(), highLimitStr.c_str(), &sliderValue, static_cast<float>(lowLimit), static_cast<float>(highLimit));

        // //Set public value to slider value
        *outputValue = static_cast<int>(roundf(sliderValue));
        *outputValue = Clamp(*outputValue, lowLimit, highLimit);

        //Draw slider text
        DrawText(TextFormat("%s: %i", sliderText.c_str(), *outputValue), 40, yposition+30, 20, DARKGRAY);


    }
};

class LogSlider {
private:
    std::string name;
    std::string lowLimitStr;
    std::string highLimitStr;
    std::string sliderText;
    float lowLimit;
    float highLimit;
    float sliderValue;
    // int pos;
    double* outputValue;
public:
    int height = 60;

    LogSlider(double &value, std::string lowLimitStr, std::string highLimitStr, float lowLimit, float highLimit, std::string sliderText) :
        outputValue(&value), lowLimitStr(lowLimitStr), highLimitStr(highLimitStr), lowLimit(lowLimit), highLimit(highLimit), sliderText(sliderText)  {

    }

    void draw(float yposition) {
        float sliderValue = static_cast<float>(std::log10(*outputValue));

        GuiSlider(Rectangle{40, yposition, 300, 20}, lowLimitStr.c_str(), highLimitStr.c_str(), &sliderValue, std::log10(lowLimit), std::log10(highLimit));

        *outputValue = pow(10,sliderValue);

        DrawText(TextFormat("%s: 10^%.1f", sliderText.c_str(), sliderValue), 40, yposition+30, 20, DARKGRAY);

    }
};

class Tickbox {
private:
    std::string name;
    std::string tickboxText;
    bool ticked;
    // int pos;
    bool* outputValue;
public:
    int height = 40;

    Tickbox(bool &value, std::string tickboxText) :
        outputValue(&value), tickboxText(tickboxText)  {

    }

    void draw(float yposition) {
        bool ticked = *outputValue;

        GuiCheckBox(Rectangle{40, yposition, 20, 20}, TextFormat("%s", tickboxText.c_str()), &ticked);

        // //Set public value
        *outputValue = ticked;
    }
};


struct DropdownOptionSetup {
    std::string display_name;
    int option_number;
    std::string description;
};


class Dropdown {
private:
    std::string name;
    std::string dropdownText;
    bool editMode = false;
    std::string allOptionNames;
    int& selectedOption;
    std::vector<DropdownOptionSetup> options;
    float yposition;

public:
    int height = 80;

    Dropdown(int &value, std::string dropdownTitle) :
        selectedOption(value), dropdownText(dropdownTitle)  {

    }

    void addOption(DropdownOptionSetup option) {
        options.push_back(option);
        if (allOptionNames.empty()) {
            allOptionNames = option.display_name;
        }
        else {
            allOptionNames += ";" + option.display_name;
        }
    }
    void saveYpos(float ypos) {
        yposition = ypos;
    }

    void draw() {
        DrawText(TextFormat("%s:", dropdownText.c_str()), 40, yposition, 10, DARKGRAY);

        DrawText(TextFormat("%s", options[selectedOption].description.c_str()), 40, yposition+50, 10, DARKGRAY);

        if (GuiDropdownBox(Rectangle{40, yposition+15, 150, 30}, allOptionNames.c_str(), &selectedOption, editMode)) {
            editMode = !editMode;
        }
    }
};

//Info Boxes
void DrawGridUnitOverlay(double plotSize) {
    double metresPerGridUnit = plotSize / 10.0;

    DrawRectangle(10, GetScreenHeight()-80, 500, 45, Fade(BLACK, 0.6f));

    DrawText(TextFormat("1 grid unit = %.1e m / %.1e AU / %.1epc", metresPerGridUnit,metresPerGridUnit*AUpermetre,metresPerGridUnit*pcpermetre),15,GetScreenHeight()-68,20,RAYWHITE);
}

void DrawControlsOverlay() {
    DrawRectangle(GetScreenWidth() - 310, GetScreenHeight()-180, 500, 175, Fade(BLACK, 0.6f));

    DrawText("Click+Drag: Rotate View",GetScreenWidth() - 295,GetScreenHeight()-175,20,RAYWHITE);

    DrawText("Scroll: Zoom",GetScreenWidth() - 295,GetScreenHeight()-150,20,RAYWHITE);

    DrawText("Arrows: Move Camera",GetScreenWidth() - 295,GetScreenHeight()-125,20,RAYWHITE);

    DrawText("F3: Toggle Info",GetScreenWidth() - 295,GetScreenHeight()-100,20,RAYWHITE);

    DrawText("F5: Toggle Diagnostics",GetScreenWidth() - 295,GetScreenHeight()-75,20,RAYWHITE);

    DrawText("+/-: Speed",GetScreenWidth() - 295,GetScreenHeight()-50,20,RAYWHITE);
}


void DrawTimeOverlay(double elapsedYears, double timestepYears, double yearsPerSecond, int updatesPerSecond) {
    DrawRectangle(10, 10, 360, 120, Fade(BLACK, 0.6f));

    DrawText(TextFormat("Elapsed time: %.2f years", elapsedYears), 15, 22, 20, RAYWHITE);

    DrawText(TextFormat("Timestep: %.3e years", timestepYears), 15, 47, 20, RAYWHITE);

    DrawText(TextFormat("Rate: %.3e years/s", yearsPerSecond), 15, 72, 20, RAYWHITE);

    DrawText(TextFormat("Physics Updates Per Second: %i", updatesPerSecond), 15, 97, 20, RAYWHITE);
}

void DrawPerformanceOverlay(double lastUpdateMs, double averageUpdateMs) {
    DrawRectangle(GetScreenWidth() - 310, 10, 300, 85, Fade(BLACK, 0.6f));

    DrawText(TextFormat("Update: %.3f ms", lastUpdateMs), GetScreenWidth() - 295, 22, 20, RAYWHITE);

    DrawText(TextFormat("Average: %.3f ms", averageUpdateMs), GetScreenWidth() - 295, 47, 20, RAYWHITE);

    DrawText(TextFormat("Max Updates/s: %.1f", 1000/averageUpdateMs), GetScreenWidth() - 295, 72, 20, RAYWHITE);
}

void DrawRollingGraph(const std::deque<double>& history, Rectangle bounds, const char* label, Color lineColour) {
    DrawRectangleRec(bounds, Fade(BLACK, 0.45f));
    DrawRectangleLinesEx(bounds, 1.0f, Fade(RAYWHITE, 0.35f));

    DrawText(label, bounds.x + 6, bounds.y + 5, 14, RAYWHITE);

    if (history.size() < 2) {
        DrawText("Waiting for data...", bounds.x + 6, bounds.y + 28, 12, GRAY);
        return;
    }

    double minValue = history[0];
    double maxValue = history[0];

    for (double value : history) {
        minValue = std::min(minValue, value);
        maxValue = std::max(maxValue, value);
    }

    double range = maxValue - minValue;

    if (range == 0.0) {
        range = std::abs(maxValue);
    }

    if (range == 0.0) {
        range = 1.0;
    }

    float graphLeft = bounds.x + 8;
    float graphRight = bounds.x + bounds.width - 8;
    float graphTop = bounds.y + 24;
    float graphBottom = bounds.y + bounds.height - 18;

    for (size_t i = 1; i < history.size(); i++) {
        float t0 = static_cast<float>(i - 1) / static_cast<float>(history.size() - 1);
        float t1 = static_cast<float>(i) / static_cast<float>(history.size() - 1);

        float x0 = graphLeft + t0 * (graphRight - graphLeft);
        float x1 = graphLeft + t1 * (graphRight - graphLeft);

        float y0 = graphBottom - static_cast<float>((history[i - 1] - minValue) / range) * (graphBottom - graphTop);
        float y1 = graphBottom - static_cast<float>((history[i] - minValue) / range) * (graphBottom - graphTop);

        DrawLineEx(Vector2{x0, y0}, Vector2{x1, y1}, 2.0f, lineColour);
    }

    DrawText(TextFormat("min %.2e", minValue), bounds.x + 6, bounds.y + bounds.height - 15, 10, GRAY);
    DrawText(TextFormat("max %.2e", maxValue), bounds.x + bounds.width - 95, bounds.y + bounds.height - 15, 10, GRAY);
}

void DrawDiagnosticsOverlay(
    double initialKE,
    double currentKE,
    double initialPE,
    double currentPE,
    double initialEnergy,
    double currentEnergy,
    int unboundParticles,
    double unboundFraction,
    const std::deque<double>& KEhistory,
    const std::deque<double>& PEhistory,
    const std::deque<double>& energyHistory,
    const std::deque<double>& unboundHistory
) {
    float panelWidth = 360.0f;
    float panelHeight = 640.0f;

    float x = GetScreenWidth() - panelWidth - 20.0f;
    float y = (GetScreenHeight() - panelHeight) / 2.0f;

    DrawRectangle(static_cast<int>(x), static_cast<int>(y), static_cast<int>(panelWidth), static_cast<int>(panelHeight), Fade(BLACK, 0.65f));

    DrawRectangleLinesEx(Rectangle{x, y, panelWidth, panelHeight}, 1.0f, Fade(RAYWHITE, 0.5f));

    DrawText("Physics Diagnostics", x + 12, y + 10, 20, RAYWHITE);

    double energyErrorPercent = 0.0;

    if (initialEnergy != 0.0) {
        energyErrorPercent = 100.0 * (currentEnergy - initialEnergy) / std::abs(initialEnergy);
    }

    DrawText(TextFormat("KE start: %.3e J", initialKE), x + 12, y + 42, 14, RAYWHITE);
    DrawText(TextFormat("KE now:   %.3e J", currentKE), x + 12, y + 60, 14, RAYWHITE);

    DrawText(TextFormat("PE start: %.3e J", initialPE), x + 12, y + 84, 14, RAYWHITE);
    DrawText(TextFormat("PE now:   %.3e J", currentPE), x + 12, y + 102, 14, RAYWHITE);

    DrawText(TextFormat("E start:  %.3e J", initialEnergy), x + 12, y + 126, 14, RAYWHITE);
    DrawText(TextFormat("E now:    %.3e J", currentEnergy), x + 12, y + 144, 14, RAYWHITE);
    DrawText(TextFormat("Error:    %.6f %%", energyErrorPercent), x + 12, y + 166, 16, RAYWHITE);

    DrawText(TextFormat("Unbound:  %i / %.2f %%", unboundParticles, 100.0 * unboundFraction), x + 12, y + 188, 16, RAYWHITE);

    float graphX = x + 12;
    float graphW = panelWidth - 24;
    float graphH = 85;

    DrawRollingGraph(KEhistory, Rectangle{graphX, y + 220, graphW, graphH}, "Kinetic Energy", SKYBLUE);
    DrawRollingGraph(PEhistory, Rectangle{graphX, y + 315, graphW, graphH}, "Potential Energy", ORANGE);
    DrawRollingGraph(energyHistory, Rectangle{graphX, y + 410, graphW, graphH}, "Total Energy", LIME);
    DrawRollingGraph(unboundHistory, Rectangle{graphX, y + 505, graphW, graphH}, "Unbound Particles", RED);
}


//Setup GUI
void OpenSetupGUI(AllSettings& settings) {
    float WINDOW_HEIGHT = 700;
    float WINDOW_WIDTH = 500;

    InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Gravity Simulation Settings");
    SetTargetFPS(60);

    SetupScreen screen = SetupScreen::Main;

    bool startPressed = false;

    int selectedPreset = static_cast<int>(settings.init.preset);
    Dropdown presetChoice(selectedPreset, "Preset");
    presetChoice.addOption({"Default", 0, "Generate Random Clusters"});
    presetChoice.addOption({"Solar System", 1, "Model Solar System"});
    presetChoice.addOption({"Milky Way", 2, "Milky Way Model(ish)"});

    IntSlider nSlider(settings.user.particleN, "20", "1000", 20, 1000, "Total Particles");
    IntSlider trailSlider(settings.user.maxTrailLength, "1", "1000", 0.0f, 1000.0f, "Trail Length");
    Tickbox centralMassTick(settings.user.enableCentralMass, "Central Mass");
    LogSlider plotsizeSlider(settings.user.plotSize, "2e10", "2e20", 2e10, 2e20, "Plot Size (m)");
    Tickbox drawTrailsTick(settings.user.drawTrails, "Enable Trails");

    IntSlider clusterSlider(settings.init.clusterCount, "1", "10", 1, 10, "Number of Clusters");

    int selectedColourIndex = static_cast<int>(settings.user.colourScheme);
    Dropdown colourChoice(selectedColourIndex, "Colour Scheme");
    colourChoice.addOption({"White", 0, "All generated particles white"});
    colourChoice.addOption({"Random", 1, "All generated particles assigned random colours"});
    colourChoice.addOption({"Distance", 2, "Particles grouped and coloured by distance from centre"});
    colourChoice.addOption({"Cluster", 3, "Particles in a cluster will have the same colour"});

    while (!WindowShouldClose() && !startPressed) {
        BeginDrawing();
        ClearBackground(RAYWHITE);

        float nextYpos = 40;

        if (screen == SetupScreen::Main) {
            DrawText("Settings for Gravity Simulation", 40, nextYpos, 24, BLACK);
            nextYpos += 60;

            presetChoice.saveYpos(nextYpos);
            nextYpos += presetChoice.height;
            settings.init.preset = static_cast<SystemPreset>(selectedPreset);

            if (settings.init.preset != SystemPreset::SolarSystem) {
                nSlider.draw(nextYpos);
                nextYpos += nSlider.height;

                plotsizeSlider.draw(nextYpos);
                nextYpos += plotsizeSlider.height;

                centralMassTick.draw(nextYpos);
                nextYpos += centralMassTick.height;

                colourChoice.saveYpos(nextYpos);
                nextYpos += colourChoice.height;
                settings.user.colourScheme = static_cast<ColourOptions>(selectedColourIndex);

                if (GuiButton(Rectangle{40, WINDOW_HEIGHT - 100, 400, 35}, "Initial Conditions")) {
                    screen = SetupScreen::InitialConditions;
                }
            }

            drawTrailsTick.draw(nextYpos);
            nextYpos += drawTrailsTick.height;

            if (settings.user.drawTrails) {
                trailSlider.draw(nextYpos);
                nextYpos += trailSlider.height;
            }

            if (GuiButton(Rectangle{40, WINDOW_HEIGHT - 50, 140, 35}, "Start")) {
                startPressed = true;
            }

            presetChoice.draw();
            if (settings.init.preset != SystemPreset::SolarSystem) {
                colourChoice.draw();
            }
        }

        else if (screen == SetupScreen::InitialConditions) {
            DrawText("Initial Conditions", 40, nextYpos, 24, BLACK);
            nextYpos += 60;

            clusterSlider.draw(nextYpos);
            nextYpos += clusterSlider.height;

            if (GuiButton(Rectangle{40, WINDOW_HEIGHT - 50, 140, 35}, "Save")) {
                screen = SetupScreen::Main;
            }
        }

        EndDrawing();
    }

    CloseWindow();
}



int main() {

    InitialisationQueue initqueue;
    //Create settings struct
    AllSettings settings;

    // settings.init.preset = SystemPreset::Clusters;

    OpenSetupGUI(settings);
    if (settings.init.preset == SystemPreset::Clusters) {
        autoClusters(settings,initqueue);
    }
    if (settings.init.preset == SystemPreset::SolarSystem) {
        buildSolarSystem(settings,initqueue);
    }
    if (settings.init.preset == SystemPreset::MilkyWay) {
        buildMilkyWay(settings,initqueue);
    }
    

    std::vector<Particle> particles = InitialiseAll(initqueue);
    //Create Simulation
    Simulation sim(settings,particles);
    //Add diagnostics
    PhysicsDiagnostics diagnostics(sim.getParticles());

    TrailManager trailManager(sim.getParticles().size(), settings.user.maxTrailLength, settings.user.plotSize);

    int windowWidth = 1200;
    int windowHeight = 1200;


    bool hideInfo = false;
    bool hideDiagnostics = true;
    bool drawAxes = true;


    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    InitWindow(windowWidth, windowHeight, "Gravity Simulation");
    SetWindowMinSize(500,500);

    int monitor = GetCurrentMonitor();
    int monitorHeight = GetMonitorHeight(monitor);
    int monitorWidth = GetMonitorWidth(monitor);
    Vector2 monitorPos = GetMonitorPosition(monitor);

    if (monitorHeight < windowHeight || monitorWidth < windowWidth) {
        windowWidth = static_cast<int>(monitorHeight * 0.80f);
        windowHeight = static_cast<int>(monitorHeight * 0.80f);

    }
    SetWindowSize(windowWidth,windowHeight);
    SetWindowPosition(static_cast<int>(monitorPos.x + (monitorWidth - windowWidth) / 2),static_cast<int>(monitorPos.y + (monitorHeight - windowHeight) / 2));

    Camera3D camera = {};

    camera.target = { 0.0f, 0.0f, 0.0f};
    camera.up = { 0.0f, 1.0f, 0.0f};
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    float cameraYaw = 0.8f;
    float cameraPitch = 0.45f;
    float cameraDistance = 18.0f;

    SetTargetFPS(60);

    double targetUpdatesPerSecond = 60.0;
    double secondsPerUpdate = 1.0 / targetUpdatesPerSecond;

    const double maxFrameSeconds = 1.0/15;
    const int maxUpdatesPerFrame = 20;

    double accumulator = 0.0;

    double energyTimer = 0;

    while (!WindowShouldClose()) {
        double frameSeconds = GetFrameTime();
        frameSeconds = std::min(frameSeconds,maxFrameSeconds);

        accumulator += frameSeconds;

        int updatesThisFrame = 0;

        while (accumulator >= secondsPerUpdate && updatesThisFrame < maxUpdatesPerFrame) {
            sim.update();

            if (settings.user.drawTrails) {
                trailManager.record(sim.getParticles());
            }

            accumulator -= secondsPerUpdate;
            updatesThisFrame++;
        }

        if (updatesThisFrame == maxUpdatesPerFrame) {
            accumulator = 0.0;
        }

        if (!hideDiagnostics) {
            energyTimer += frameSeconds;

            if (energyTimer >= 1.0) {
                diagnostics.sample();
                energyTimer = 0.0;
            }
        }
        else {
            energyTimer = 0.0;
        }

        UpdateOrbitCamera(camera, cameraYaw, cameraPitch, cameraDistance);

        hideInfo = checkHideInfo(hideInfo);
        hideDiagnostics = checkHideDiagnostics(hideDiagnostics);
        drawAxes = checkDrawAxes(drawAxes);

        targetUpdatesPerSecond = UpdateSimSpeed(targetUpdatesPerSecond);
        secondsPerUpdate = 1.0 / targetUpdatesPerSecond;

        // UpdateParticlePos(particles[particles.size()-1], settings);

        // UpdateParticleMass(particles[particles.size()-1]);

        BeginDrawing();

        ClearBackground(BLACK);

        BeginMode3D(camera);

        
        if (drawAxes) {
            DrawGrid(20, 1.0f);

            DrawLine3D({ -10.0f, 0.0f, 0.0f }, { 10.0f, 0.0f, 0.0f }, RED);
            DrawLine3D({ 0.0f, -10.0f, 0.0f }, { 0.0f, 10.0f, 0.0f }, GREEN);
            DrawLine3D({ 0.0f, 0.0f, -10.0f }, { 0.0f, 0.0f, 10.0f }, BLUE);
        }        


        const std::vector<Particle> &particles = sim.getParticles();

        if (settings.user.drawTrails) {
            // DrawTrail(particles[i]);
            trailManager.draw(particles);
        }

        for (size_t i = 0; i < particles.size(); i++) {
            Vect3 ppos = particles[i].getPos();
            float renderScale = 10.0f / settings.user.plotSize;

            float px = static_cast<float>(ppos.x * renderScale);
            float py = static_cast<float>(ppos.y * renderScale);
            float pz = static_cast<float>(ppos.z * renderScale);

            float circleSize;

            circleSize = GetVisualRadius(particles[i].getMass());


            Vector3 pposV3 = ToVector3({px,py,pz});

            Color circleColor = particles[i].getColour();

            DrawSphere(pposV3,circleSize,circleColor);
        }

        EndMode3D();

        double yearsPerSecond = sim.getTimestepYears() * targetUpdatesPerSecond;

        if (!hideInfo) {
            DrawTimeOverlay(sim.getCurrentTimeYears(), sim.getTimestepYears(), yearsPerSecond, targetUpdatesPerSecond);
            DrawGridUnitOverlay(settings.user.plotSize);
            DrawPerformanceOverlay(sim.getLastUpdateMs(), sim.getAverageUpdateMs());
            DrawControlsOverlay();


        }
        if (!hideDiagnostics) {
            DrawDiagnosticsOverlay(
                diagnostics.getInitialKE(),
                diagnostics.getCurrentKE(),
                diagnostics.getInitialPE(),
                diagnostics.getCurrentPE(),
                diagnostics.getInitialEnergy(),
                diagnostics.getCurrentEnergy(),
                diagnostics.getCurrentUnboundParticles(),
                diagnostics.getCurrentUnboundFraction(),
                diagnostics.getKEhistory(),
                diagnostics.getPEhistory(),
                diagnostics.getEnergyHistory(),
                diagnostics.getUnboundHistory()
            );
        }

        EndDrawing();
    }

    CloseWindow();

    // OrbitParameters earthOrbit = {1.49598e11,0.0167,8.73e-7,-0.197,1.797,1.753,0};
    // double earthMass = 5.9722e24;

    // ParticleState earthState = stateFromKeplerianOrbit(earthOrbit,M0,earthMass,0);

    // std::cout << earthState.position.magnitude() << " " << earthState.velocity.magnitude();
}