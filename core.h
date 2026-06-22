#pragma once

#include <cmath>
#include "raylib.h"
#include <vector>


//Global Constants
constexpr double pi = 3.141592653589793;
constexpr double G = 6.67430e-11;
constexpr double M0 = 1.989e30; //kg
constexpr double secondsPerYear = 365.25 * 24.0 * 60.0 * 60.0;
constexpr double AUpermetre = 1/1.496e11;
constexpr double pcpermetre = 1/3.086e16;

//Enums
enum class ColourOptions {
    White,
    Random,
    Distance,
    Cluster,
    None
};

enum class SystemPreset {
    Clusters,
    SolarSystem,
    MilkyWay,
    JupiterMoons,
    GlobularCluster,
    DiskGalaxy,
    MWandAndromeda,
    EllipticalGalaxy,
    GalaxyCluster,
    CosmicWeb
};

enum class ClusterShape {
    Disk,
    Sphere,
};

enum class VelocityScheme {
    None,
    Circular,
    Random,
    Virialised
};

enum class MassDistribution {
    Uniform,
    KroupaIMF,
};

enum class SetupScreen {
    Main,
    InitialConditions
};

struct Vect3 {
    double x = 0;
    double y = 0;
    double z = 0;
    double magnitude() const {
        return std::sqrt(x*x+y*y+z*z);
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



inline Vect3 operator*(double scalar, const Vect3& vector) {
    return vector * scalar;
}


//Settings
struct OneParticleSettings {
    double mass;
    Vect3 position; //Set manually
    Vect3 velocity; //Set manually or automatically if autoOrbit = true;
    Vect3 axis; //Accounts for inclination and orbit direction
    double eccentricity;
    bool autoOrbit;
    Color particleColour = BLUE;
};

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
    int maxTrailLength = 200;
    // int selectedColourOption = 1;
    ColourOptions colourScheme = ColourOptions::Random;
    int selectedWindowSize = 0;
};

struct InitSettings {
    //Dropdowns/Enum
    SystemPreset preset = SystemPreset::Clusters;

    int clusterCount = 1;
    //Add cluster menu
};

struct ClusterSettings {
    double size; //Radius - determines particle count (distributed from total UserSettings.particleN)
    double distancefromorigin; //Distance of cluster from 0,0,0;
    Vect3 position; //Position of cluster centre - default determined randomly from distancefromorigin
    
    ClusterShape shape = ClusterShape::Disk;
    double flattening;
    Vect3 axis; //Rotation axis of disk
    double velocityScheme = 0; //0=Perfect circular orbit, 1=Purely Random?
    double velocityScale = 0; //-1=Slower than orbit, 0=Perfect Circular, 1=Faster Than Circular
    Vect3 systemVelocity = {0,0,0}; //Motion of cluster in simulation frame
    bool centralMassEnabled = true; //Central Mass in cluster
    double centralMass = 1e34;
    MassDistribution IMF = MassDistribution::KroupaIMF; //0 = Uniform Random Masses, 1 = Kroupa 2002 IMF-like
    double lowMassBound = 1 * M0;
    double highMassBound = 100 * M0;
    int particleCount;
    Color clusterColour = WHITE;
    double clusterGeneratedMass;
    bool rescaleMasses;
};


struct AllSettings {
    UserSettings user;
    InitSettings init;
    ClusterSettings cluster;
    OneParticleSettings particle;
};

//Holds candidate objects (particles/clusters) while user configures
class InitialisationQueue {
private:
    std::vector<OneParticleSettings> particles;
    std::vector<ClusterSettings> clusters;
public:
    std::vector<OneParticleSettings> getParticles() {
        return particles;
    }
    std::vector<ClusterSettings> getClusters() {
        return clusters;
    }
    void addParticle(OneParticleSettings settings) {
        particles.push_back(settings);
    }
    void addCluster(ClusterSettings settings) {
        clusters.push_back(settings);
    }
    void removeParticle(int index) {
        particles.erase(particles.begin() + index);
    }
    void removeCluster(int index) {
        clusters.erase(clusters.begin() + index);
    }
};