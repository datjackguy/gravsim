#include "generation.h"

#include <random>
#include <algorithm>
#include <cmath>


std::random_device rd;
std::mt19937 gen(rd());

//Compute the magnitude of a 3D vector
double pythagoras(double x, double y, double z) {
    double r = sqrt(x*x+y*y+z*z);
    return r;
}

//Overloaded for any vector length
double pythagoras(std::vector<double> vector) {
    double sumsq = 0;
    for (size_t i = 0; i < vector.size(); i++) {
        sumsq += vector[i]*vector[i];
    }
    return sqrt(sumsq);
}

//Vect3-Dependent Utilities
//Overloaded for Vect3
double pythagoras(Vect3 vector) {
    return sqrt(vector.x*vector.x+vector.y*vector.y+vector.z*vector.z);
}


double randomDouble(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(gen);
}

//Generate a random integer between min, max
int randomInteger(int min, int max) {
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

double clampDouble(double value, double low, double high) {
    return std::max(low, std::min(value, high));
}

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
    for (size_t i = 0; i < masses.size(); i++) {
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
        masses[i] = M0*mass_candidate;
    }
    return masses;
}

Vect3 randomAxis() {
    double x = randomDouble(0,1);
    double y = randomDouble(0,1);
    double z = randomDouble(0,1);
    Vect3 vector{x,y,z};
    return vector.normalise();
}

// Generate position of individual particle within sphere* of maxRad
Vect3 randomSph(double maxRad) {
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

Vect3 generateDisk(double minRad, double maxRad, double flattening, Vect3 axis) {
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

    double rho = std::sqrt(
        minRad * minRad
        + u * (maxRad * maxRad - minRad * minRad)
    );

    double cosTheta = randomDouble(-1.0, 1.0);
    double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    double phi = randomDouble(0.0, 2.0 * pi);

    double xLocal = rho * std::cos(phi);
    double yLocal = rho * std::sin(phi);

    double zLocal = flattening * rho * randomDouble(-1.0, 1.0);

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

Vect3 diskVelocity(Vect3 pos, double orbitMass, Vect3 axis, double softening) {
    Vect3 n = axis.normalise();

    Vect3 radial = pos - n * pos.dot(n);
    double orbitalRadius = radial.magnitude();

    Vect3 tangent = n.cross(radial).normalise();

    double distanceFromCentre = pos.magnitude();

    double softenedDistance2 = distanceFromCentre * distanceFromCentre + softening * softening;

    double softenedDistance = std::sqrt(softenedDistance2);

    double speed = std::sqrt(G * orbitMass * orbitalRadius * orbitalRadius / (softenedDistance2 * softenedDistance));

    return tangent * speed;
}