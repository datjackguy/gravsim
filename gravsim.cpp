//Compile Portable -- g++ gravsim.cpp -o gravsim.exe -I include -L lib -lraylib -lgdi32 -lwinmm -static-libgcc -static-libstdc++ -static
//2 -- g++ gravsim.cpp -o gravsim.exe -I include -L lib -static -static-libgcc -static-libstdc++ -lraylib -lopengl32 -lgdi32 -lwinmm
//3 Optimised g++ -O2 gravsim.cpp -o gravsim.exe -I include -L lib -lraylib -lgdi32 -lwinmm -static-libgcc -static-libstdc++ -static
//4 g++ -O3 gravsim.cpp -o gravsim.exe -I include -L lib -lraylib -lgdi32 -lwinmm -static-libgcc -static-libstdc++ -static
//5 With Warning -- g++ -std=c++20 -Wall -Wextra -Wpedantic -O3 gravsim.cpp -o gravsim.exe -I include -L lib -lraylib -lgdi32 -lwinmm -static-libgcc -static-libstdc++ -static

//To-Do
//Colour options enum
//Presets enum - Solar System (Plus Pluto, Moon and Asteroid Belt), Jupiter system, Globular Cluster, Milky Way Scale, Milky Way / Andromeda Collision 
//Timestep/Softening Automation

#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <thread>
#include <chrono>
#include <algorithm>

#include "raylib.h"
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"
#include "raymath.h"

const double pi = 3.141592653589793;
const double G = 6.67430e-11;
// const double epsilon = 5e13;
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

double clampDouble(double value, double low, double high)
{
    return std::max(low, std::min(value, high));
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
    Vect3 normalise() const {
        double mag = magnitude();
        return {x / mag,y / mag,z / mag};
    }
    Vect3 cross(const Vect3& other) const {
        return {y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x};
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
    bool operator==(const Vect3& other) const {
        if (x == other.x && y == other.y && z == other.z) {
            return true;
        }
        else {
            return false;
        }
    }
    double dot(const Vect3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
};

Vect3 operator*(double scalar, const Vect3& vector) {
    return vector * scalar;
}

Vect3 randomAxis() {
    double x = randomDouble(0,1);
    double y = randomDouble(0,1);
    double z = randomDouble(0,1);
    Vect3 vector{x,y,z};
    return vector.normalise();
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
    double softening = 0.0;

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

    void setSoftening(double newSoft) {
        softening = newSoft;
    }

    double getSoftening() const {
        return softening;
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
    vector.y = rho * sinTheta * std::sin(phi);
    vector.z = rho * cosTheta;

    return vector;

}

Vect3 polarToCartesian(double r, double theta, double phi) {
    Vect3 vector;
    vector.x = r * std::sin(theta)*std::cos(phi);
    vector.y = r * std::sin(theta)*std::sin(phi);
    vector.z = r * std::cos(theta);
    return vector;
}


Vect3 generateDisk(double maxRad, double flattening, Vect3 axis) {
    Vect3 n = axis.normalise();

    if (flattening < 0.0) flattening = 0.0;
    if (flattening > 1.0) flattening = 1.0;

    Vect3 temp;
    temp = {0.0, 0.0, 1.0};
    if (temp.cross(n) == Vect3{0,0,0}) {
        temp = {1.0, 0.0, 0.0};
    }

    Vect3 e1 = temp.cross(n).normalise();
    Vect3 e2 = n.cross(e1);

    double u = randomDouble(0.0, 1.0);
    double rho = maxRad * std::cbrt(u);

    double cosTheta = randomDouble(-1.0, 1.0);
    double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    double phi = randomDouble(0.0, 2.0 * pi);

    double xLocal = rho * sinTheta * std::cos(phi);
    double yLocal = rho * sinTheta * std::sin(phi);

    double zLocal = flattening * rho * cosTheta;

    Vect3 vector = e1 * xLocal + e2 * yLocal + n * zLocal;

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

Vect3 diskVelocity(Vect3 pos, double orbitMass, Vect3 axis) {
    Vect3 n = axis.normalise();

    Vect3 radial = pos - n * pos.dot(n);

    double r = radial.magnitude();

    Vect3 tangent = n.cross(radial).normalise();

    double speed = std::sqrt(G * orbitMass / r);

    return tangent * speed;
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



std::vector<Color> colourList = {RED,GREEN,BLUE,ORANGE,PURPLE,GREEN,MAGENTA,SKYBLUE};


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

struct InitSettings {
    //Dropdowns/Enum
    int preset = 1; //0 - Custom, 1 - Globular Cluster, 2 - Galaxy, 3 - Solar System etc.

    int clusterCount = 1;
    //Add cluster menu
};

struct ClusterSettings {
    double size; //Radius - determines particle count (distributed from total UserSettings.particleN)
    double distancefromorigin; //Distance of cluster from 0,0,0;
    Vect3 position; //Position of cluster centre - default determined randomly from distancefromorigin
    
    int shape = 0; //Presets - Spherical, Disk (Thickness), Lobed, Spiral etc.
    Vect3 axis; //Rotation axis of disk
    double velocityScheme = 0; //0=Perfect circular orbit, 1=Purely Random?
    double velocityScale = 0; //-1=Slower than orbit, 0=Perfect Circular, 1=Faster Than Circular
    double systemVelocity = 0; //Motion of cluster in simulation frame
    bool centralMassEnabled = true; //Central Mass in cluster
    double centralMass = 1e35;
    int IMF = 0; //0 = Uniform Random Masses, 1 = Kroupa 2002 IMF-like
    double lowMassBound = 1 * M0;
    double highMassBound = 100 * M0;
};

struct OneParticleSettings {
    double mass;
    Vect3 position;
    Vect3 velocity;

};

struct AllSettings {
    UserSettings user;
    InitSettings init;
    ClusterSettings cluster;
    OneParticleSettings particle;
};

class Cluster {
private:
    Vect3 centrepos;
    Vect3 vel;
    Vect3 axis;
    double flattening;
    int particleCount;
    double radius;

    std::vector<Particle> particles;
    bool centralMassExists = false;

public:
    Cluster(Vect3 centrepos, Vect3 vel, Vect3 axis, double flattening, int particleCount, double radius) :
        centrepos(centrepos), vel(vel), axis(axis), flattening(flattening), particleCount(particleCount), radius(radius) {
        
        initialise_particles();
    
    }
    //Calculate total mass within the radius of a particle
    double enclosedMass(const Particle& particle) const {
        Vect3 relativePos = particle.getPos() - centrepos;
        double radius = relativePos.magnitude();

        double mass = 0.0;

        for (int p = 0; p < particles.size(); p++) {
            Particle other = particles[p];
            Vect3 otherRelativePos = other.getPos() - centrepos;
            double otherRadius = otherRelativePos.magnitude();

            if (otherRadius < radius) {
                mass += other.getMass();
            }
        }

        return mass;
    }
    void AddCentralMass(double mass) {
        if (!centralMassExists) {
            particles.push_back(Particle{mass,centrepos,vel});
            centralMassExists=true;
        }

        //Recalculate Velocities
        for (int i = 0; i < particleCount; i++) {
            double encMass = enclosedMass(particles[i]);
            Vect3 currentPos = particles[i].getPos();

            Vect3 particlevel = diskVelocity(currentPos-centrepos,encMass,axis);

            particles[i].setVel(particlevel);
        }

        for (int i = 0; i < particleCount; i++) {
            particles[i].setVel(particles[i].getVel()+vel);
        }

    }
    //Create all particle objects
    void initialise_particles() {
        std::vector<double> masses = generate_masses(particleCount);

        for (int i = 0; i < particleCount; i++) {
            Vect3 particlepos = generateDisk(radius, flattening, axis);
            Vect3 particlevel{0,0,0};
            particles.push_back(Particle(masses[i],particlepos,particlevel));
        }


        for (int i = 0; i < particleCount; i++) {
            double encMass = enclosedMass(particles[i]);
            Vect3 currentPos = particles[i].getPos();

            Vect3 particlevel = diskVelocity(currentPos,encMass,axis);

            particles[i].setVel(particlevel);
        }

        for (int i = 0; i < particleCount; i++) {
            particles[i].setPos(particles[i].getPos()+centrepos);
            particles[i].setVel(particles[i].getVel()+vel);
        }

    }
    
    std::vector<Particle> getParticleList() {
        return particles;
    }
};


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

class Simulation {
private:
    AllSettings settings;

    std::vector<Particle> particles;
    std::vector<std::vector<int>> particlePairs;

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
    Simulation(const AllSettings& settings)
        : settings(settings) {
        initialise_clusters();
        // initialise_particles();
        updateMassLimits();
        momentum_centre();

        configureSoftening();
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

        double baseSoftening = 0.05 * approximateSpacing;

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

        switch (settings.user.selectedColourOption) {
            case 0:
                break;
            case 1:
                for (Particle& particle : particles) {
                    particle.setColour(randomColour(colourList));
                }
                break;
            case 2:
                for (Particle& particle : particles) {
                    particle.setColour(distanceColour(particle, settings.user.plotSize));
                }
                break;
            default:
                break;
        } 




        if (particles.empty() == false && settings.user.enableCentralMass) {
            particles.back().setColour(YELLOW);
        }
    }
    void initialise_clusters() {
        if (settings.init.clusterCount == 1) {
            Cluster cluster1(Vect3{0,0,0}, Vect3{0,0,0}, Vect3{0,1,0}, 0.5, settings.user.particleN, settings.user.plotSize/2);
            if (settings.user.enableCentralMass) {
                cluster1.AddCentralMass(5e38);
            }
            std::vector<Particle> clusterParticles = cluster1.getParticleList();
            particles.insert(particles.end(),clusterParticles.begin(),clusterParticles.end());
        }
        else {
            std::vector<int> split = GenerateClusterSizes(settings.init.clusterCount,settings.user.particleN);
            for (int c = 0; c < settings.init.clusterCount; c++) {
                Vect3 clusterCentre = randomSph(settings.user.plotSize);
                Vect3 clusterVelocity{0,0,0};
                Vect3 clusterAxis = randomAxis();
                double clusterFlattening = randomDouble(0,1);
                double clusterRadius = settings.user.plotSize/4;

                Cluster cluster(clusterCentre,clusterVelocity,clusterAxis,clusterFlattening,split[c],clusterRadius);
                if (settings.user.enableCentralMass) {
                    cluster.AddCentralMass(5e38);
                }
                std::vector<Particle> clusterParticles = cluster.getParticleList();
                particles.insert(particles.end(),clusterParticles.begin(),clusterParticles.end());
            }
        }
    }

    //Create all particle objects
    void initialise_particles() {
        std::vector<double> masses = generate_masses(settings.user.particleN);

        for (int i = 0; i < settings.user.particleN; i++) {
            Vect3 pos = randomSph(settings.user.plotSize/2);
            Vect3 vel{0,0,0};
            particles.push_back(Particle(masses[i],pos,vel));
        }
        if (settings.user.enableCentralMass) {
            Vect3 zero{0,0,0};
            particles.push_back(Particle(5e38, zero, zero));
        }
        // particles.push_back(Particle(6e38, Vect3{2e15, 0, 2e15}, 0.2*globularVel(Vect3{2e15, 0, 2e15}, 5e38)));
        // particles.push_back(Particle(6e38, Vect3{-2e15, 0, -2e15}, 0.22*globularVel(Vect3{-2e15, 0, -2e15}, 5e38)));
        updateMassLimits();

        for (int i = 0; i < settings.user.particleN; i++) {
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

        double e1 = p1.getSoftening();
        double e2 = p2.getSoftening();

        double epsilonPair = sqrt(e1*e1+e2*e2);

        force = -G*mass1*mass2*(pos1-pos2)/pow((distance*distance)+(epsilonPair*epsilonPair),1.5);

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

        const auto updateStart = std::chrono::steady_clock::now();

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
            float renderScale = 10.0f / settings.user.plotSize;

            Vector3 trailPoint = {
                static_cast<float>(pos.x * renderScale),
                static_cast<float>(pos.y * renderScale),
                static_cast<float>(pos.z * renderScale)
            };

            particles[i].trail.push_back(trailPoint);

            if (particles[i].trail.size() > settings.user.maxTrailLength) {
                particles[i].trail.erase(particles[i].trail.begin());
            }
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
};


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

    if (distance > 100.0f)
    {
        distance = 100.0f;
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

void DrawPerformanceOverlay(double lastUpdateMs, double averageUpdateMs)
{
    DrawRectangle(GetScreenWidth() - 310, 10, 300, 70, Fade(BLACK, 0.6f));

    DrawText(TextFormat("Update: %.3f ms", lastUpdateMs), GetScreenWidth() - 295, 22, 20, RAYWHITE);

    DrawText(TextFormat("Average: %.3f ms", averageUpdateMs), GetScreenWidth() - 295, 47, 20, RAYWHITE);

    DrawText(TextFormat("Updates Per Second: %.1f", 1000/averageUpdateMs), GetScreenWidth() - 295, 72, 20, RAYWHITE);
}

// void OpenInitialisationGUI(AllSettings& settings) {
//     //Open Settings GUI
//     float WINDOW_HEIGHT = 700;
//     float WINDOW_WIDTH = 500;

//     InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Initial Conditions Setup");
//     SetTargetFPS(60);

//     bool savePressed = false;
//     IntSlider clusterSlider = IntSlider(settings.init.clusterCount, "1", "10", 1, 10, "Number of Clusters");


//     while (!WindowShouldClose() && !savePressed) {
//         BeginDrawing();

//         ClearBackground(RAYWHITE);

//         float nextYpos = 40;

//         DrawText("Initial Conditions", 40, nextYpos, 24, BLACK);
//         nextYpos += 60;
        
//         clusterSlider.draw(nextYpos);
//         nextYpos += clusterSlider.height;

//         if (GuiButton(Rectangle{40, WINDOW_HEIGHT-50, 140, 35}, "Save")) {
//             savePressed = true;
//             }

//         EndDrawing();
//     }

//     CloseWindow();

// }

enum class SetupScreen {
    Main,
    InitialConditions
};

void OpenSetupGUI(AllSettings& settings) {
    float WINDOW_HEIGHT = 700;
    float WINDOW_WIDTH = 500;

    InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Gravity Simulation Settings");
    SetTargetFPS(60);

    SetupScreen screen = SetupScreen::Main;

    bool startPressed = false;

    IntSlider nSlider(settings.user.particleN, "20", "1000", 20, 1000, "Total Particles");
    IntSlider trailSlider(settings.user.maxTrailLength, "1", "100", 0.0f, 100.0f, "Trail Length");
    Tickbox centralMassTick(settings.user.enableCentralMass, "Central Mass");
    LogSlider plotsizeSlider(settings.user.plotSize, "2e14", "2e16", 2e14, 2e16, "Cluster Size (m)");
    Tickbox drawTrailsTick(settings.user.drawTrails, "Enable Trails");

    IntSlider clusterSlider(settings.init.clusterCount, "1", "10", 1, 10, "Number of Clusters");

    Dropdown colourChoice(settings.user.selectedColourOption, "Colour Scheme");
    colourChoice.addOption({"White", 0, "All generated particles white"});
    colourChoice.addOption({"Random", 1, "All generated particles assigned random colours"});
    colourChoice.addOption({"Distance", 2, "Particles grouped and coloured by distance from centre"});

    while (!WindowShouldClose() && !startPressed) {
        BeginDrawing();
        ClearBackground(RAYWHITE);

        float nextYpos = 40;

        if (screen == SetupScreen::Main) {
            DrawText("Settings for Gravity Simulation", 40, nextYpos, 24, BLACK);
            nextYpos += 60;

            nSlider.draw(nextYpos);
            nextYpos += nSlider.height;

            plotsizeSlider.draw(nextYpos);
            nextYpos += plotsizeSlider.height;

            drawTrailsTick.draw(nextYpos);
            nextYpos += drawTrailsTick.height;

            if (settings.user.drawTrails) {
                trailSlider.draw(nextYpos);
                nextYpos += trailSlider.height;
            }

            centralMassTick.draw(nextYpos);
            nextYpos += centralMassTick.height;

            colourChoice.draw(nextYpos);
            nextYpos += colourChoice.height;

            if (GuiButton(Rectangle{40, WINDOW_HEIGHT - 100, 400, 35}, "Initial Conditions")) {
                screen = SetupScreen::InitialConditions;
            }

            if (GuiButton(Rectangle{40, WINDOW_HEIGHT - 50, 140, 35}, "Start")) {
                startPressed = true;
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

    //Create settings struct
    AllSettings settings;

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
            float renderScale = 10.0f / settings.user.plotSize;

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
            if (settings.user.drawTrails) {
                DrawTrail(particles[i]);
            }
            DrawSphere(pposV3,circleSize,circleColor);
        }

        EndMode3D();

        DrawTimeOverlay(sim.getCurrentTimeYears());
        DrawGridUnitOverlay(settings.user.plotSize);
        DrawPerformanceOverlay(sim.getLastUpdateMs(), sim.getAverageUpdateMs());

        EndDrawing();
        sim.update();
    }

    CloseWindow();

}