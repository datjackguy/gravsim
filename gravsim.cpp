#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <thread>
#include <chrono>

const double pi = 3.141592653589793;
const double G = 6.67430e-11;
const double epsilon = 5e9;

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
        masses[i] = 1.989e30*mass_candidate;
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

double pythagoras(double x, double y, double z) {
    double r = sqrt(x*x+y*y+z*z);
    return r;
}

std::tuple<double, double, double> globularVel(double x, double y, double z, double totalMass) {
    double velScale;
    double radius = pythagoras(x,y,z);
    double e = 0.5;
    velScale = 0.5 * pow((G*totalMass/(radius+e)),0.25);

    double velx = randomDouble(0,1)-0.5;
    double vely = randomDouble(0,1)-0.5;
    double velz = -(velx*x+vely*y)/z;

    double mag = pythagoras(velx,vely,velz);

    velx = velx/mag*velScale;
    vely = vely/mag*velScale;
    velz = velz/mag*velScale;

    return {velx,vely,velz};

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
    std::vector<double> force = {0,0,0};

public:
    // Constructor
    Particle(double mass, double posx, double posy, double posz, double velx, double vely, double velz)
        : mass(mass), posx(posx), posy(posy), posz(posz), velx(velx), vely(vely), velz(velz), pvelx(velx), pvely(vely), pvelz(velz) {
    }

    double getMass() {
        return mass;
    }

    void resetForce() {
        force = {0,0,0};
    }

    void addForce(std::vector<double> fnew, int sign = 1) {
        force[0]=force[0]+sign*fnew[0];
        force[1]=force[1]+sign*fnew[1];
        force[2]=force[2]+sign*fnew[2];
    }

    std::vector<double> getForce() {
        return force;
    }

    // Getter for pos
    std::vector<double> getPos() {
        return {posx,posy,posz};
    }

    void setPos(double x, double y, double z) {
        posx = x;
        posy = y;
        posz = z;
    }

    // Getter for vel
    std::vector<double> getVel() {
        return {velx, vely, velz};
    }

    void setVel(double x, double y, double z) {
        velx = x;
        vely = y;
        velz = z;
    }

    void calcPos(double dt) {
        posx = posx + velx * dt;
        posy = posy + vely * dt;
        posz = posz + velz * dt;
    }

    void calcVel(double dt) {
        velx = velx + force[0] * dt / mass;
        vely = vely + force[1] * dt / mass;
        velz = velz + force[2] * dt / mass;
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

double timestep = 10*24*60*60;
double runtime = 100*365.25*24*60*60;
int frame = 0;
int frames = 10;
std::vector<double> times;
double currentTime = 0.0;
std::vector<Particle> particles;
int particleN = 10;
double plotSize = 2e11;
double totalMass = 0;
std::vector<std::vector<int>> particlePairs;

void initialise_particles() {
    std::vector<double> masses = generate_masses(particleN);
    double totalMass = 0;
    for (int i = 0; i < particleN; i++) {
        totalMass = totalMass + masses[i];
    }
    for (int i = 0; i < particleN; i++) {
        std::tuple<double, double, double> pos = randomSph(plotSize);
        double xpos = std::get<0>(pos);
        double ypos = std::get<1>(pos);
        double zpos = std::get<2>(pos);
        
        std::tuple<double,double,double> vel = globularVel(xpos,ypos,zpos,totalMass);
        double xvel = std::get<0>(vel);
        double yvel = std::get<1>(vel);
        double zvel = std::get<2>(vel);

        particles.push_back(Particle(masses[i],xpos,ypos,zpos,xvel,yvel,zvel));
    }
}

void create_particle_pairs() {
    for (int i = 0; i < particleN; i++) {
        for (int j = i+1; j < particleN; j++) {
            std::vector<int> pair = {i,j}; 
            particlePairs.push_back(pair);
        }
    }
}

void momentum_centre() {
    std::vector<double> totalMomentum = {0,0,0};
    std::vector<double> CoM0 = {0,0,0};
    for (int i = 0; i < particleN; i++) {
        CoM0[0] = CoM0[0] + particles[i].getPos()[0] * particles[i].getMass();
        CoM0[1] = CoM0[1] + particles[i].getPos()[1] * particles[i].getMass();
        CoM0[2] = CoM0[2] + particles[i].getPos()[2] * particles[i].getMass();
        totalMomentum[0] = totalMomentum[0] + particles[i].getMass()*particles[i].getVel()[0];
        totalMomentum[0] = totalMomentum[1] + particles[i].getMass()*particles[i].getVel()[1];
        totalMomentum[0] = totalMomentum[2] + particles[i].getMass()*particles[i].getVel()[2];
    }
    CoM0[0] = CoM0[0] / totalMass;
    CoM0[1] = CoM0[1] / totalMass;
    CoM0[2] = CoM0[2] / totalMass;
    for (int i = 0; i < particleN; i++) {
        double newPosx = particles[i].getPos()[0] - CoM0[0];
        double newPosy = particles[i].getPos()[1] - CoM0[1];
        double newPosz = particles[i].getPos()[2] - CoM0[2];

        double newVelx = particles[i].getVel()[0] - totalMomentum[0]/totalMass;
        double newVely = particles[i].getVel()[1] - totalMomentum[1]/totalMass;
        double newVelz = particles[i].getVel()[2] - totalMomentum[2]/totalMass;

        particles[i].setPos(newPosx,newPosy,newPosz);
        particles[i].setVel(newVelx,newVely,newVelz);
    }
}

void update() {
    std::vector<std::vector<double>> forces;
    for (int p = 0; p < particlePairs.size(); p++) {
        std::vector<double> force = forceCalc(particles[particlePairs[p][0]],particles[particlePairs[p][1]]);
        // double forcex = force[0];
        // double forcey = force[1];
        // double forcez = force[2];
        forces.push_back(force);
        particles[particlePairs[p][0]].addForce(force);
        particles[particlePairs[p][1]].addForce(force, -1);
    }
    for (int i = 0; i < particles.size(); i++) {
        particles[i].calcPos(timestep);
        particles[i].calcVel(timestep);
        particles[i].resetForce();
    }
}

int main() {
    // Random Seed
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // Manual Seed
    // std::mt19937 gen(3);

    // std::cout << "IMF Check\n";
    // std::cout << initial_mass_function(0.03) << "\n";
    // std::cout << initial_mass_function(0.3) << "\n";
    // std::cout << initial_mass_function(0.7) << "\n";
    // std::cout << initial_mass_function(1.2) << "\n";

    // std::cout << "Masses Check\n";
    // std::vector masses = generate_masses(10);
    // for (int i = 0; i < masses.size(); i++) {
    //     std::cout << masses[i] << " ";
    // }
    // std::cout << "\n";

    // std::cout << "Random Sphere Check \n";
    // std::tuple<double, double, double> randomSphere = randomSph(10);
    // std::cout << std::get<0>(randomSphere) << " " << std::get<1>(randomSphere) << " " << std::get<2>(randomSphere) << "\n"; 

    // Particle testParticle1(10, 1, 1, 1, 0, 0.2, 0);
    // Particle testParticle2(20, 3, 0, -2, 0.4, 0, 0.1);
    // double particle_distance = distanceCalc(testParticle1, testParticle2);
    // std::cout << "Distance Calc Check\n" << particle_distance << "\n";
    // std::vector<double> f = forceCalc(testParticle1,testParticle2);
    // std::cout << "Force Calc Check\n" << f[0] << " " << f[1] << " " << f[2] << "\n";

    // std::cout << "Particle Pair Check\n";
    // create_particle_pairs();
    // for (int i = 0; i < particlePairs.size(); i++) {
    //     std::cout << "[" << particlePairs[i][0] << particlePairs[i][1] << "] ";
    // }

    initialise_particles();
    create_particle_pairs();
    while (frame < frames) {
        std::cout << particles[0].getPos()[0] << "\n";
        update();
        frame++;
        std::this_thread::sleep_for(std::chrono::seconds(1));
    }

}