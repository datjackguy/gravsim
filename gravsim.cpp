//Compile Portable -- g++ gravsim.cpp -o gravsim.exe -I include -L lib -lraylib -lgdi32 -lwinmm -static-libgcc -static-libstdc++ -static
//2 -- g++ gravsim.cpp -o gravsim.exe -I include -L lib -static -static-libgcc -static-libstdc++ -lraylib -lopengl32 -lgdi32 -lwinmm
//3 Optimised g++ -O2 gravsim.cpp -o gravsim.exe -I include -L lib -lraylib -lgdi32 -lwinmm -static-libgcc -static-libstdc++ -static
//4 g++ -O3 gravsim.cpp -o gravsim.exe -I include -L lib -lraylib -lgdi32 -lwinmm -static-libgcc -static-libstdc++ -static
//5 With Warning -- g++ -std=c++20 -Wall -Wextra -Wpedantic -O3 gravsim.cpp -o gravsim.exe -I include -L lib -lraylib -lgdi32 -lwinmm -static-libgcc -static-libstdc++ -static

//To-Do
//Colour options enum
//Presets enum - Solar System (Plus Pluto, Moon and Asteroid Belt), Jupiter system, Globular Cluster, Milky Way Scale, Milky Way / Andromeda Collision 

#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <thread>
#include <chrono>

#include "raylib.h"
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"
#include "raymath.h"

const double pi = 3.141592653589793;
const double G = 6.67430e-11;
const double epsilon = 5e13;
const double M0 = 1.989e30; //kg
const double secondsPerYear = 365.25 * 24.0 * 60.0 * 60.0;
const double AUpermetre = 1/1.496e11;
const double pcpermetre = 1/3.086e16;
// bool enableCentralMass = true;
// int selectedColourOption = 0;

//Initialise random generator
// std::mt19937 gen(3);
std::random_device rd;
std::mt19937 gen(rd());

//Generate a random double between min, max
double randomDouble(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(gen);
}

//Generate a random integer between min, max
double randomInteger(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(gen);
}

int randomSign() {
    std::uniform_int_distribution<int> dist(-1, 1);
    int sign = 0;
    while (sign == 0) {
        sign = dist(gen);
    }
    return sign;
}

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

//Compute the magnitude of a 3D vector
double pythagoras(double x, double y, double z) {
    double r = sqrt(x*x+y*y+z*z);
    return r;
}



//Overloaded for any vector length
double pythagoras(std::vector<double> vector) {
    double sumsq = 0;
    for (int i = 0; i < vector.size(); i++) {
        sumsq += vector[i]*vector[i];
    }
    return sqrt(sumsq);
}

struct Vect3 {
    double x = 0;
    double y = 0;
    double z = 0;
    double magnitude() const {
        return pythagoras(x,y,z);
    }
    void set(double x2, double y2, double z2) {
        x = x2;
        y = y2;
        z = z2;
    }
    Vect3& operator+=(const Vect3& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
    Vect3& operator-=(const Vect3& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }
    Vect3 operator+(const Vect3& other) const {
        Vect3 result = *this;
        result += other;
        return result;
    }
    Vect3 operator-(const Vect3& other) const {
        Vect3 result = *this;
        result -= other;
        return result;
    }
    Vect3 operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }
    Vect3 operator/(double scalar) const {
        return {x / scalar, y / scalar, z / scalar};
    }
};

Vect3 operator*(double scalar, const Vect3& vector) {
    return vector * scalar;
}

//Overloaded for Vect3
double pythagoras(Vect3 vector) {
    return sqrt(vector.x*vector.x+vector.y*vector.y+vector.z*vector.z);
}


//Particle Class
class Particle {
private:
    double mass;
    Vect3 pos;
    Vect3 vel;
    Vect3 force;
    Color colour = WHITE;

public:
    std::vector<Vector3> trail;
    // Constructor
    Particle(double mass, Vect3 pos, Vect3 vel)
        : mass(mass), pos(pos), vel(vel) {
    }

    // Add (or subtract) force from total force acting on particle
    void addForce(Vect3 fnew, int sign = 1) {
        force += sign * fnew;
    }

    // Calculate the new position
    void calcPos(double dt) {
        pos += vel * dt;
    }

    // Calculate the new velocity
    void calcVel(double dt) {
        vel += force * dt / mass;
    }

    // Reset force vector to 0
    void resetForce() {
        force.set(0,0,0);
    }

    //GETTERS/SETTERS
    void setColour(Color newColour) {
        colour = newColour;
    }

    Color getColour() const {
        return colour;
    }

    double getMass() const {
        return mass;
    }

    void setMass(double newMass) {
        mass = newMass;
    }

    Vect3 getForce() const {
        return force;
    }

    Vect3 getPos() const {
        return pos;
    }

    void setPos(Vect3 newPos) {
        pos = newPos;
    }

    Vect3 getVel() const {
        return vel;
    }

    void setVel(Vect3 newVel) {
        vel = newVel;
    }
};

//Kroupa 2002 IMF 
double initial_mass_function(double m) {
    if (m < 0.08) {
        return 0.0;
    }
    else if (m < 0.5) {
        return std::pow(m, -1.3);
    }
    else if (0.5 < m && m < 1) {
        return std::pow(m, -2.3);
    }
    else {
        return std::pow(m, -2.7);
    }
}

//Generate N masses according to IMF
std::vector<double> generate_masses(int N) {
    std::vector<double> masses(N);
    for (int i = 0; i < masses.size(); i++) {
        double mass_candidate = 0.0;
        while (mass_candidate == 0.0) {
            //Random Mass (M_0?)
            double x_candidate = randomDouble(0.0,20.0);
            //Random Probability Threshold?
            double y_candidate = randomDouble(0.0,1.0);
            if (y_candidate <= initial_mass_function(x_candidate)) {
                mass_candidate = x_candidate;
            }
        }
        masses[i] = 3e4*M0*mass_candidate;
    }
    return masses;
}


// Generate position of individual particle within sphere* of maxRad
Vect3 randomSph(double maxRad) {
    // double rho = randomDouble(maxRad/20,maxRad);
    // double theta = randomDouble(0.0,pi);
    // double phi = randomDouble(0.0, 2*pi);
    // double x = rho * std::sin(theta) * std::cos(phi);
    // double y = 0.4 * rho * std::sin(theta) * std::sin(phi);
    // double z = rho * std::cos(theta);
    // Vect3 vector;
    // vector.x = x;
    // vector.y = y;
    // vector.z = z;

    // return vector;

    double u = randomDouble(0.0,1.0);
    double rho = maxRad * std::cbrt(u);

    double cosTheta = randomDouble(-1.0, 1.0);
    double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    double phi = randomDouble(0.0, 2*pi);

    Vect3 vector;
    
    vector.x = rho * sinTheta * std::cos(phi);
    vector.y = 0.4 * rho * sinTheta * std::sin(phi);
    vector.z = rho * cosTheta;

    return vector;

}



//Compute the velocity of a particle at x,y,z
Vect3 globularVel(Vect3 pos, double totalMass) {
    double velScale;
    double radius = pythagoras(pos);
    double e = 1e12;
    velScale = 2* pow((G*totalMass/(radius+e)),0.5);

    double x = pos.x;
    double y = pos.y;
    double z = pos.z;

    // double velx = randomDouble(0,1)-0.5;
    // double vely = randomDouble(0,1)-0.5;
    // double velz = -(velx*x+vely*y)/z;

    // double mag = pythagoras(velx,vely,velz);

    // velx = 2.5*velx/mag*velScale;
    // vely = 0.1*2*vely/mag*velScale;
    // velz = 2.5*velz/mag*velScale;


    //Cylindrical Polar
    double r = sqrt(x*x+z*z);
    double phi = atan2(z,x);
    //r velocity ~ 0
    double dr_dt = randomDouble(-0.1,0.1);
    //y velocity ~ 0
    double dy_dt = randomDouble(-0.2,0.2);
    //phi velocity - rotation
    double omega = sqrt(G*totalMass/(radius*radius*radius));
    double velx = dr_dt*cos(phi)-r*omega*sin(phi);
    double vely = dy_dt;
    double velz = dr_dt*sin(phi)+r*omega*cos(phi);

    double mag = pythagoras(velx,vely,velz);

    // velx = velx/mag*velScale;
    // vely = vely/mag*velScale;
    // velz = velz/mag*velScale;
    

    Vect3 vector;
    vector.x = velx;
    vector.y = vely;
    vector.z = velz;

    // int sign = randomSign();
    // vector = vector * sign;

    return vector;

}




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


//Simulation Parameters
// double timestep = 0.5*10*24*60*60; //s
// double runtime = 100*365.25*24*60*60; //s
// int frame = 0;
// int frames = 10;
// std::vector<double> times;
// double currentTime = 0.0;
// std::vector<Particle> particles;
// int particleN = 300;
// double plotSize = 2e15; //m
// double totalMass = 0; //kg
// std::vector<std::vector<int>> particlePairs;
// double maxMass = 0;
// double minMass = 1e64;
// int maxTrailLength = 20;
std::vector<Color> colourList = {RED,GREEN,BLUE,ORANGE,PURPLE,GREEN,MAGENTA,SKYBLUE};
// bool randomColours = false;
// bool drawTrails = true;

struct UserSettings {
    //Directly Modifiable
    //Toggles
    // bool randomColours = false;
    // bool distanceColours = false;
    bool enableCentralMass = true;
    bool drawTrails = true;
    //Values
    int particleN = 300;
    double plotSize = 2e15;
    int maxTrailLength = 20;
    int selectedColourOption = 1;
    int selectedWindowSize = 0;
};

class Simulation {
private:
    UserSettings settings;

    std::vector<Particle> particles;
    std::vector<std::vector<int>> particlePairs;

    double timestep = 0.5*10*24*60*60;
    double totalMass = 0.0;
    double maxMass = 0.0; //Initialise maximum present mass, determined after particles created
    double minMass = 1e64; //As above for minimum mass
    double currentTime = 0.0;

    bool forcesExist = false;

public:
    Simulation(const UserSettings& settings)
        : settings(settings) {
        initialise_particles();
        momentum_centre();
        // create_particle_pairs();
        colour_options();

    }
private:
    void resetForces() {
        for (Particle& particle : particles) {
            particle.resetForce();
        }
    }
    //Colour options
    void colour_options() {
        // if (settings.randomColours) {
        //     for (Particle& particle : particles) {
        //         particle.setColour(randomColour(colourList));
        //     }
        // }

        // else if (settings.distanceColours) {
        //     for (Particle& particle : particles) {
        //         particle.setColour(distanceColour(particle, settings.plotSize));
        //     }
        // }

        switch (settings.selectedColourOption) {
            case 0:
                break;
            case 1:
                for (Particle& particle : particles) {
                    particle.setColour(randomColour(colourList));
                }
                break;
            case 2:
                for (Particle& particle : particles) {
                    particle.setColour(distanceColour(particle, settings.plotSize));
                }
                break;
            default:
                break;
        } 




        if (particles.empty() == false && settings.enableCentralMass) {
            particles.back().setColour(YELLOW);
        }
    }


    //Create all particle objects
    void initialise_particles() {
        std::vector<double> masses = generate_masses(settings.particleN);

        for (int i = 0; i < settings.particleN; i++) {
            Vect3 pos = randomSph(settings.plotSize/2);
            Vect3 vel{0,0,0};
            particles.push_back(Particle(masses[i],pos,vel));
        }
        if (settings.enableCentralMass) {
            Vect3 zero{0,0,0};
            particles.push_back(Particle(5e38, zero, zero));
        }
        // particles.push_back(Particle(6e38, Vect3{2e15, 0, 2e15}, 0.2*globularVel(Vect3{2e15, 0, 2e15}, 5e38)));
        // particles.push_back(Particle(6e38, Vect3{-2e15, 0, -2e15}, 0.22*globularVel(Vect3{-2e15, 0, -2e15}, 5e38)));
        updateMassLimits();

        for (int i = 0; i < settings.particleN; i++) {
            double encMass = enclosedMass(particles[i]);
            Vect3 currentPos = particles[i].getPos();

            Vect3 vel = globularVel(currentPos,encMass);

            particles[i].setVel(vel);
        }

    }
    void updateMassLimits() {
        totalMass = 0.0;
        minMass = 1e64;
        maxMass = 0.0;
        for (int i = 0; i < particles.size(); i++) {
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
        for (int p = 0; p < particles.size(); p++) {
            Vect3 currentPos = particles[p].getPos();

            double currentRadius = pythagoras(currentPos);
            if (currentRadius < radius) {
                mass += particles[p].getMass();
            }
        }
        return mass;
    }
    // Create list containing the indices of each pair of particles
    void create_particle_pairs() {
        particlePairs.clear();

        for (int i = 0; i < particles.size(); i++) {
            for (int j = i+1; j < particles.size(); j++) {
                std::vector<int> pair = {i,j}; 
                particlePairs.push_back(pair);
            }
        }
    }

    // Calculate the distance between two particles
    double distanceCalc(const Particle &p1, const Particle &p2) {
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

        force = -G*mass1*mass2*(pos1-pos2)/pow((distance*distance)+(epsilon*epsilon),1.5);

        return force;
    }

    void calculateForces() {
        resetForces();

        for (int i = 0; i < particles.size(); i++) {
            for (int j = i+1; j < particles.size(); j++) {
                Vect3 force = forceCalc(particles[i],particles[j]);
                particles[i].addForce(force);
                particles[j].addForce(force, -1);
            }
        }
    }

public:
    //Update the simulation for each timestep
    void update() {
        // for (int p = 0; p < particlePairs.size(); p++) {
        //     Vect3 force = forceCalc(particles[particlePairs[p][0]],particles[particlePairs[p][1]]);

        //     particles[particlePairs[p][0]].addForce(force);
        //     particles[particlePairs[p][1]].addForce(force, -1);
        // }
        if (forcesExist == false) {
            calculateForces();
            forcesExist = true;
        }

        for (int i = 0; i < particles.size(); i++) {
            particles[i].calcVel(timestep/2);
        }

        for (int i = 0; i < particles.size(); i++) {
            particles[i].calcPos(timestep);
        }

        calculateForces();

        for (int i = 0; i < particles.size(); i++) {
            particles[i].calcVel(timestep/2);
        }

        for (int i = 0; i < particles.size(); i++) {
            // particles[i].calcVel(timestep);
            // particles[i].calcPos(timestep);
            // particles[i].resetForce();

            Vect3 pos = particles[i].getPos();
            float renderScale = 10.0f / settings.plotSize;

            Vector3 trailPoint = {
                static_cast<float>(pos.x * renderScale),
                static_cast<float>(pos.y * renderScale),
                static_cast<float>(pos.z * renderScale)
            };

            particles[i].trail.push_back(trailPoint);

            if (particles[i].trail.size() > settings.maxTrailLength) {
                particles[i].trail.erase(particles[i].trail.begin());
            }
        }
        currentTime+=timestep;
    }

    const std::vector<Particle>& getParticles() const {
        return particles;
    }

    double getPlotSize() const {
        return settings.plotSize;
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
        return settings.drawTrails;
    }

    double getCurrentTimeYears() const {
        return currentTime / secondsPerYear;
    }
};












//Calculate the centre of momentum for the initial system - should produce a simulation that remains centred on origin
// void momentum_centre() {
//     std::vector<double> totalMomentum = {0,0,0};
//     std::vector<double> CoM0 = {0,0,0};
//     for (int i = 0; i < particles.size(); i++) {
//         CoM0[0] = CoM0[0] + particles[i].getPos()[0] * particles[i].getMass();
//         CoM0[1] = CoM0[1] + particles[i].getPos()[1] * particles[i].getMass();
//         CoM0[2] = CoM0[2] + particles[i].getPos()[2] * particles[i].getMass();
//         totalMomentum[0] = totalMomentum[0] + particles[i].getMass()*particles[i].getVel()[0];
//         totalMomentum[1] = totalMomentum[1] + particles[i].getMass()*particles[i].getVel()[1];
//         totalMomentum[2] = totalMomentum[2] + particles[i].getMass()*particles[i].getVel()[2];
//     }
//     CoM0[0] = CoM0[0] / totalMass;
//     CoM0[1] = CoM0[1] / totalMass;
//     CoM0[2] = CoM0[2] / totalMass;
//     for (int i = 0; i < particles.size(); i++) {
//         double newPosx = particles[i].getPos()[0] - CoM0[0];
//         double newPosy = particles[i].getPos()[1] - CoM0[1];
//         double newPosz = particles[i].getPos()[2] - CoM0[2];

//         double newVelx = particles[i].getVel()[0] - totalMomentum[0]/totalMass;
//         double newVely = particles[i].getVel()[1] - totalMomentum[1]/totalMass;
//         double newVelz = particles[i].getVel()[2] - totalMomentum[2]/totalMass;

//         particles[i].setPos(newPosx,newPosy,newPosz);
//         particles[i].setVel(newVelx,newVely,newVelz);
//     }
// }



// Drawing Functions -----
//Enable user control of camera view (zoom and rotate)
void UpdateOrbitCamera(Camera3D& camera, float& yaw, float& pitch, float& distance)
{
    Vector2 mouseDelta = GetMouseDelta();

    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT))
    {
        yaw -= mouseDelta.x * 0.01f;
        pitch += mouseDelta.y * 0.01f;
    }

    // Limit the pitch so the camera does not flip upside down.
    const float pitchLimit = 1.5f;

    if (pitch > pitchLimit)
    {
        pitch = pitchLimit;
    }

    if (pitch < -pitchLimit)
    {
        pitch = -pitchLimit;
    }

    distance -= GetMouseWheelMove() * 1.0f;

    if (distance < 3.0f)
    {
        distance = 3.0f;
    }

    if (distance > 40.0f)
    {
        distance = 40.0f;
    }

    //Move Camera Target
    if (IsKeyDown(KEY_LEFT)) {
        camera.target.x -= 0.1f;
    }
    if (IsKeyDown(KEY_RIGHT)) {
        camera.target.x += 0.1f;
    }
    if (IsKeyDown(KEY_UP)) {
        camera.target.z -= 0.1f;
    }
    if (IsKeyDown(KEY_DOWN)) {
        camera.target.z += 0.1f;
    }

    // Convert yaw, pitch and distance into a 3D camera position.
    camera.position.x = camera.target.x + distance * cosf(pitch) * sinf(yaw);
    camera.position.y = camera.target.y + distance * sinf(pitch);
    camera.position.z = camera.target.z + distance * cosf(pitch) * cosf(yaw);
}

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

    double visualMinMass = 1e33;
    double visualMaxMass = 1e40;

    double logMin = std::log10(visualMinMass);
    double logMax = std::log10(visualMaxMass);
    double logMass = std::log10(mass);

    float t = static_cast<float>((logMass - logMin) / (logMax - logMin));

    t = Clamp(t, 0.0f, 1.0f);
    t = std::pow(t, 1.4f);

    return minRadius + t * (maxRadius - minRadius);
}



void DrawTrail(const Particle &particle) {
    for (int i = 1; i < particle.trail.size(); i++) {
        DrawLine3D(particle.trail[i-1], particle.trail[i], particle.getColour());
    }
}

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
        // float sliderValue = static_cast<float>(*outputValue);

        // //Draw slider - set position and limits
        // GuiSlider(Rectangle{40, yposition, 300, 20}, lowLimitStr.c_str(), highLimitStr.c_str(), &sliderValue, static_cast<float>(lowLimit), static_cast<float>(highLimit));

        // // //Set public value to slider value
        // *outputValue = static_cast<int>(roundf(sliderValue));
        // *outputValue = Clamp(*outputValue, lowLimit, highLimit);

        // //Draw slider text
        // DrawText(TextFormat("%s: %i", sliderText.c_str(), *outputValue), 40, yposition+30, 20, DARKGRAY);

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

        //Draw slider - set position and limits
        // GuiSlider(Rectangle{40, yposition, 300, 20}, lowLimitStr.c_str(), highLimitStr.c_str(), &sliderValue, static_cast<float>(lowLimit), static_cast<float>(highLimit));
        GuiCheckBox(Rectangle{40, yposition, 20, 20}, TextFormat("%s", tickboxText.c_str()), &ticked);

        // //Set public value to slider value
        *outputValue = ticked;


    }
};

enum class ColourOptions {
    White,
    RandomColours,
    DistanceColours
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

    void draw(float yposition) {
        // bool ticked = *outputValue;

        //Draw slider - set position and limits
        // GuiSlider(Rectangle{40, yposition, 300, 20}, lowLimitStr.c_str(), highLimitStr.c_str(), &sliderValue, static_cast<float>(lowLimit), static_cast<float>(highLimit));
        // GuiCheckBox(Rectangle{40, yposition, 20, 20}, TextFormat("%s", tickboxText.c_str()), &ticked);
        DrawText(TextFormat("%s:", dropdownText.c_str()), 40, yposition, 10, DARKGRAY);



        // //Set public value to slider value
        // *outputValue = ticked;
        DrawText(TextFormat("%s", options[selectedOption].description.c_str()), 40, yposition+50, 10, DARKGRAY);


        // if (GuiDropdownBox(Rectangle{150, yposition-10, 150, 30}, allOptionNames.c_str(), &selectedOption, editMode)) {
        //     editMode = !editMode;
        // }
        if (GuiDropdownBox(Rectangle{40, yposition+15, 150, 30}, allOptionNames.c_str(), &selectedOption, editMode)) {
            editMode = !editMode;
        }
    }
};

void DrawGridUnitOverlay(double plotSize)
{
    double metresPerGridUnit = plotSize / 10.0;

    DrawRectangle(10, GetScreenHeight()-80, 500, 45, Fade(BLACK, 0.6f));

    DrawText(TextFormat("1 grid unit = %.1e m / %.1e AU / %.1epc", metresPerGridUnit,metresPerGridUnit*AUpermetre,metresPerGridUnit*pcpermetre),15,GetScreenHeight()-68,20,RAYWHITE);
}

void DrawTimeOverlay(double elapsedYears)
{
    DrawRectangle(10, 10, 250, 45, Fade(BLACK, 0.6f));

    DrawText(TextFormat("Elapsed time: %.2f years", elapsedYears),15,22,20,RAYWHITE);
}

void OpenSetupGUI(UserSettings& settings) {
    //Open Settings GUI
    float WINDOW_HEIGHT = 700;
    float WINDOW_WIDTH = 500;

    InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Gravity Simulation Settings");
    SetTargetFPS(60);

    IntSlider nSlider = IntSlider(settings.particleN, "20", "1000", 20, 1000, "Total Particles");
    IntSlider trailSlider = IntSlider(settings.maxTrailLength, "1", "100", 0.0f, 100.0f, "Trail Length");
    Tickbox centralMassTick = Tickbox(settings.enableCentralMass, "Central Mass");
    LogSlider plotsizeSlider = LogSlider(settings.plotSize, "2e14", "2e16", 2e14, 2e16, "Cluster Size (m)");
    bool startPressed = false;
    Tickbox drawTrailsTick = Tickbox(settings.drawTrails, "Enable Trails");
    // Tickbox randomColoursTick = Tickbox(settings.randomColours, "Random Colours");

    
    Dropdown colourChoice = Dropdown(settings.selectedColourOption,"Colour Scheme");
    colourChoice.addOption({"White", 0, "All generated particles white"});
    colourChoice.addOption({"Random", 1, "All generated particles assigned random colours"});
    colourChoice.addOption({"Distance", 2, "Particles grouped and coloured by distanced from centre"});


    while (!WindowShouldClose() && !startPressed) {
        BeginDrawing();

        ClearBackground(RAYWHITE);

        float nextYpos = 40;

        DrawText("Settings for Gravity Simulation", 40, nextYpos, 24, BLACK);
        nextYpos += 60;
        
        nSlider.draw(nextYpos);
        nextYpos += nSlider.height;

        plotsizeSlider.draw(nextYpos);
        nextYpos += plotsizeSlider.height;

        // GuiCheckBox(Rectangle{40, 230, 20, 20}, "Enable Trails", &drawTrails);
        drawTrailsTick.draw(nextYpos);
        nextYpos += drawTrailsTick.height;

        if (settings.drawTrails) {
            trailSlider.draw(nextYpos);
            nextYpos += trailSlider.height;
        }

        // GuiCheckBox(Rectangle{40, 320, 20, 20}, "Random Colours", &randomColours);
        // randomColoursTick.draw(nextYpos);
        // nextYpos += randomColoursTick.height;

        centralMassTick.draw(nextYpos);
        nextYpos += centralMassTick.height;

        colourChoice.draw(nextYpos);
        nextYpos += colourChoice.height;

        if (GuiButton(Rectangle{40, WINDOW_HEIGHT-50, 140, 35}, "Start")) {
            startPressed = true;
            }

        EndDrawing();
    }

    CloseWindow();
}




int main() {

    //Create settings struct
    UserSettings settings;

    OpenSetupGUI(settings);

    //Start Simulation
    Simulation sim(settings);

    // momentum_centre();


    int windowWidth = 1200;
    int windowHeight = 1200;





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

    const std::vector<Particle> &particles = sim.getParticles();
    while (!WindowShouldClose()) {
        float dt = GetFrameTime();
        float time = GetTime();

        UpdateOrbitCamera(camera, cameraYaw, cameraPitch, cameraDistance);

        // UpdateParticlePos(particles[particles.size()-1], settings);

        // UpdateParticleMass(particles[particles.size()-1]);

        BeginDrawing();

        ClearBackground(BLACK);

        BeginMode3D(camera);

        DrawGrid(20, 1.0f);


        DrawLine3D({ -10.0f, 0.0f, 0.0f }, { 10.0f, 0.0f, 0.0f }, RED);
        DrawLine3D({ 0.0f, -10.0f, 0.0f }, { 0.0f, 10.0f, 0.0f }, GREEN);
        DrawLine3D({ 0.0f, 0.0f, -10.0f }, { 0.0f, 0.0f, 10.0f }, BLUE);

        for (int i = 0; i < particles.size(); i++) {
            Vect3 ppos = particles[i].getPos();
            float renderScale = 10.0f / settings.plotSize;

            float px = static_cast<float>(ppos.x * renderScale);
            float py = static_cast<float>(ppos.y * renderScale);
            float pz = static_cast<float>(ppos.z * renderScale);

            // float circleSize = particles[i].getMass()/1.989e30; 
            float circleSize;// = std::log(particles[i].getMass())/std::log((1.989e34))*2;
            // circleSize = 0.3f * cbrt(particles[i].getMass()/totalMass);
            // circleSize = GetVisualRadius(particles[i].getMass(), sim.getTotalMass(), sim.getMinMass(), sim.getMaxMass(), particles.size());
            circleSize = GetVisualRadius(particles[i].getMass());

            // if (particles[i].getMass() > 1e37) {
            //     circleSize = 0.4f;
            // }
            Vector3 pposV3 = ToVector3({px,py,pz});

            Color circleColor = particles[i].getColour();
            // if (i == 0) {
            //     circleColor = GREEN;
            // }
            // else if (i == particles.size()-1) {
            //     circleColor = YELLOW;
            // }
            // else if (i == particles.size()-1) {
            //     circleColor = RED;
            // }
            if (settings.drawTrails) {
                DrawTrail(particles[i]);
            }
            DrawSphere(pposV3,circleSize,circleColor);
        }

        EndMode3D();

        DrawTimeOverlay(sim.getCurrentTimeYears());
        DrawGridUnitOverlay(settings.plotSize);

        EndDrawing();
        sim.update();
    }

    CloseWindow();

}