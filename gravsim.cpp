#include <iostream>
#include <cmath>
#include <vector>
#include <random>

const double pi = 3.141592653589793;
const double G = 6.67430e-11;
const double epsilon = 0.1;

double randomDouble(double min, double max) {
    // std::mt19937 gen(3);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(min, max);
    return dist(gen);
}

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

std::vector<double> generate_masses(int N) {
    std::vector<double> masses(N);
    // std::uniform_real_distribution<double> x_candidates(0.0, 20.0);
    // std::uniform_real_distribution<double> y_candidates(0.0, 1.0);
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
        masses[i] = mass_candidate;
    }
    return masses;
}

std::tuple<double, double, double> randomSph(double maxRad) {
    double rho = randomDouble(0.0,maxRad);
    double theta = randomDouble(0.0,pi);
    double phi = randomDouble(0.0, 2*pi);
    double x = rho * std::sin(theta) * std::cos(phi);
    double y = rho * std::sin(theta) * std::sin(phi);
    double z = rho * std::cos(theta);
    return {x, y, z};
}

std::tuple<double, double, double> globularVel() {
    return {0,0,0};
}

void mergeParticles() {
}

class Particle {
private:
    // Private member variables
    double mass;
    double posx;
    double posy;
    double posz;
    double velx;
    double vely;
    double velz;
    double pvelx;
    double pvely;
    double pvelz;

public:
    // Constructor
    Particle(double mass, double posx, double posy, double posz, double velx, double vely, double velz)
        : mass(mass), posx(posx), posy(posy), posz(posz), velx(velx), vely(vely), velz(velz), pvelx(velx), pvely(vely), pvelz(velz) {
    }

    double getMass() {
        return mass;
    }

    // Getter for pos
    std::vector<double> getPos() {
        return {posx,posy,posz};
    }

    // Getter for vel
    std::vector<double> getVel() {
        return {velx, vely, velz};
    }

    void calcPos(double dt) {
        posx = posx + velx * dt;
        posy = posy + vely * dt;
        posz = posz + velz * dt;
    }

    void calcVel(double dt, std::vector<double> netForce) {
        velx = velx + netForce[0] * dt / mass;
        vely = vely + netForce[1] * dt / mass;
        velz = velz + netForce[2] * dt / mass;
    }
};

double distanceCalc(Particle &p1, Particle &p2) {
    double distance = sqrt(pow(p1.getPos()[0]-p2.getPos()[0],2)+pow(p1.getPos()[1]-p2.getPos()[1],2)+pow(p1.getPos()[2]-p2.getPos()[2],2));
    return distance;
}

std::vector<double> forceCalc(Particle &p1, Particle &p2) {
    double distance = distanceCalc(p1,p2);
    double mass1 = p1.getMass();
    double mass2 = p2.getMass();
    std::vector<double> pos1 = p1.getPos();
    std::vector<double> pos2 = p2.getPos();
    // std::vector<double> vel1 = p1.getVel();
    // std::vector<double> vel2 = p2.getVel();

    double forcex = -G*mass1*mass2*(pos1[0]-pos2[0])/pow((distance*distance)+(epsilon*epsilon),1.5);
    double forcey = -G*mass1*mass2*(pos1[1]-pos2[1])/pow((distance*distance)+(epsilon*epsilon),1.5);
    double forcez = -G*mass1*mass2*(pos1[2]-pos2[2])/pow((distance*distance)+(epsilon*epsilon),1.5);

    return {forcex,forcey,forcez};
}

int main() {
    // Random Seed
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // Manual Seed
    // std::mt19937 gen(3);

    std::cout << "IMF Check\n";
    std::cout << initial_mass_function(0.03) << "\n";
    std::cout << initial_mass_function(0.3) << "\n";
    std::cout << initial_mass_function(0.7) << "\n";
    std::cout << initial_mass_function(1.2) << "\n";

    std::cout << "Masses Check\n";
    std::vector masses = generate_masses(10);
    for (int i = 0; i < masses.size(); i++) {
        std::cout << masses[i] << " ";
    }
    std::cout << "\n";

    std::cout << "Random Sphere Check \n";
    std::tuple<double, double, double> randomSphere = randomSph(10);
    std::cout << std::get<0>(randomSphere) << " " << std::get<1>(randomSphere) << " " << std::get<2>(randomSphere) << "\n"; 

    Particle testParticle1(10, 1, 1, 1, 0, 0.2, 0);
    Particle testParticle2(20, 3, 0, -2, 0.4, 0, 0.1);
    double particle_distance = distanceCalc(testParticle1, testParticle2);
    std::cout << "Distance Calc Check\n" << particle_distance << "\n";
    std::vector<double> f = forceCalc(testParticle1,testParticle2);
    std::cout << "Force Calc Check\n" << f[0] << " " << f[1] << " " << f[2] << "\n";

}