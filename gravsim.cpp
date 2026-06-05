//Compile Portable -- g++ gravsim.cpp -o gravsim.exe -I include -L lib -lraylib -lgdi32 -lwinmm -static-libgcc -static-libstdc++ -static

#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <thread>
#include <chrono>

#include "raylib.h"
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

const double pi = 3.141592653589793;
const double G = 6.67430e-11;
const double epsilon = 5e13;

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
        masses[i] = 3e4*1.989e30*mass_candidate;
    }
    return masses;
}

std::tuple<double, double, double> randomSph(double maxRad) {
    // double rho = randomDouble(0.0,maxRad);
    double rho = randomDouble(maxRad/20,maxRad);
    double theta = randomDouble(0.0,pi);
    double phi = randomDouble(0.0, 2*pi);
    double x = rho * std::sin(theta) * std::cos(phi);
    double y = 0.4 * rho * std::sin(theta) * std::sin(phi);
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
    double e = 1e12;
    velScale = 2* pow((G*totalMass/(radius+e)),0.5);

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

    // std::cout << velx << " " << vely << " " << velz << "\n";

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
        // std::cout << "setPos called!";
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

double timestep = 0.5*10*24*60*60;
double runtime = 100*365.25*24*60*60;
int frame = 0;
int frames = 10;
std::vector<double> times;
double currentTime = 0.0;
std::vector<Particle> particles;
int particleN = 300;
double plotSize = 2e15;
double totalMass = 0;
std::vector<std::vector<int>> particlePairs;

double enclosedMass(Particle &particle) {
    std::vector<double> ppos = particle.getPos();
    double x = ppos[0];
    double y = ppos[1];
    double z = ppos[2];
    double radius = pythagoras(x,y,z);
    double mass = 0;
    for (int p = 0; p < particles.size(); p++) {
        std::vector<double> currentPos = particles[p].getPos();
        double cx = currentPos[0];
        double cy = currentPos[1];
        double cz = currentPos[2];
        double currentRadius = pythagoras(cx,cy,cz);
        if (currentRadius < radius) {
            mass += particles[p].getMass();
        }
    }
    return mass;
}


void initialise_particles() {
    std::vector<double> masses = generate_masses(particleN);
    // double totalMass = 0;
    for (int i = 0; i < particleN; i++) {
        totalMass = totalMass + masses[i];
    }
    // std::cout << "Total Mass: " << totalMass;
    for (int i = 0; i < particleN; i++) {
        std::tuple<double, double, double> pos = randomSph(plotSize/2);
        double xpos = std::get<0>(pos);
        double ypos = std::get<1>(pos);
        double zpos = std::get<2>(pos);
        
        std::tuple<double,double,double> vel = globularVel(xpos,ypos,zpos, totalMass);
        double xvel = std::get<0>(vel);
        double yvel = std::get<1>(vel);
        double zvel = std::get<2>(vel);

        particles.push_back(Particle(masses[i],xpos,ypos,zpos,xvel,yvel,zvel));
    }
    particles.push_back(Particle(5e38, 0, 0, 0, 0, 0, 0));
    particles.push_back(Particle(5e38, 1e16, 0, 1e16, -7e6, 5e4, -6e6));
    for (int i = 0; i < particleN; i++) {
        double encMass = enclosedMass(particles[i]);
        std::vector<double> currentPos = particles[i].getPos();
        double xpos = currentPos[0];
        double ypos = currentPos[1];
        double zpos = currentPos[2];

        std::tuple<double,double,double> vel = globularVel(xpos,ypos,zpos,encMass);
        double xvel = std::get<0>(vel);
        double yvel = std::get<1>(vel);
        double zvel = std::get<2>(vel);

        particles[i].setVel(xvel,yvel,zvel);
    }
}



void create_particle_pairs() {
    for (int i = 0; i < particles.size(); i++) {
        for (int j = i+1; j < particles.size(); j++) {
            std::vector<int> pair = {i,j}; 
            particlePairs.push_back(pair);
        }
    }
}

void momentum_centre() {
    std::vector<double> totalMomentum = {0,0,0};
    std::vector<double> CoM0 = {0,0,0};
    for (int i = 0; i < particles.size(); i++) {
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
    for (int i = 0; i < particles.size(); i++) {
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

    // Convert yaw, pitch and distance into a 3D camera position.
    camera.position.x = camera.target.x + distance * cosf(pitch) * sinf(yaw);
    camera.position.y = camera.target.y + distance * sinf(pitch);
    camera.position.z = camera.target.z + distance * cosf(pitch) * cosf(yaw);
}

void UpdateParticlePos(Particle &particle) {
    std::vector<double> ppos = particle.getPos();
    double x = ppos[0];
    double y = ppos[1];
    double z = ppos[2];
    bool changed = false;
    if (IsKeyDown(KEY_S)) {
        x = x + plotSize/120;
        changed = true;
    }
    if (IsKeyDown(KEY_W)) {
        x = x - plotSize/120;
        changed = true;
    }
    if (IsKeyDown(KEY_D)) {
        z = z - plotSize/120;
        changed = true;
    }
    if (IsKeyDown(KEY_A)) {
        z = z + plotSize/120;
        changed = true;
    }

    //...The rest
    if (changed) {
        particle.setPos(x,y,z);
    }

}

Vector3 ToVector3(const std::vector<double>& values)
{
    return Vector3{
        static_cast<float>(values[0]),
        static_cast<float>(values[1]),
        static_cast<float>(values[2])
    };
}

int main() {

    InitWindow(500, 300, "Gravity Simulation Settings");
    SetTargetFPS(60);

    bool startPressed = false;
    while (!WindowShouldClose() && !startPressed) {
        BeginDrawing();

        ClearBackground(RAYWHITE);

        DrawText("Settings for Gravity Simulation", 40, 40, 24, BLACK);

        float particleSlider = particleN;

        GuiSlider(Rectangle{40, 100, 300, 20}, "20", "400", &particleSlider, 20, 400);

        particleN = static_cast<int>(roundf(particleSlider));

        DrawText(TextFormat("Total Particles: %i", particleN), 40, 140, 20, DARKGRAY);

        if (GuiButton(Rectangle{40, 190, 140, 35}, "Start")) {
            startPressed = true;
            }

        EndDrawing();
    }

    CloseWindow;




    initialise_particles();
    create_particle_pairs();
    // momentum_centre();

    int windowWidth = 1200;
    int windowHeight = 1200;

    InitWindow(windowWidth, windowHeight, "Gravity Simulation");

    Camera3D camera = {};

    camera.target = { 0.0f, 0.0f, 0.0f};
    camera.up = { 0.0f, 1.0f, 0.0f};
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    float cameraYaw = 0.8f;
    float cameraPitch = 0.45f;
    float cameraDistance = 18.0f;

    SetTargetFPS(60);

    // std::cout << "Total Mass" << totalMass;

    while (!WindowShouldClose()) {
        float dt = GetFrameTime();
        float time = GetTime();

        UpdateOrbitCamera(camera, cameraYaw, cameraPitch, cameraDistance);

        UpdateParticlePos(particles[particles.size()-2]);

        BeginDrawing();

        ClearBackground(BLACK);

        BeginMode3D(camera);

        DrawGrid(20, 1.0f);


        DrawLine3D({ -10.0f, 0.0f, 0.0f }, { 10.0f, 0.0f, 0.0f }, RED);
        DrawLine3D({ 0.0f, -10.0f, 0.0f }, { 0.0f, 10.0f, 0.0f }, GREEN);
        DrawLine3D({ 0.0f, 0.0f, -10.0f }, { 0.0f, 0.0f, 10.0f }, BLUE);

        for (int i = 0; i < particles.size(); i++) {
            std::vector ppos = particles[i].getPos();
            float renderScale = 10.0f / plotSize;

            float px = static_cast<float>(ppos[0] * renderScale);
            float py = static_cast<float>(ppos[1] * renderScale);
            float pz = static_cast<float>(ppos[2] * renderScale);

            // float circleSize = particles[i].getMass()/1.989e30; 
            float circleSize;// = std::log(particles[i].getMass())/std::log((1.989e34))*2;
            // std::cout << circleSize << "\n";
            if (particles[i].getMass() > 1e36) {
                circleSize = 0.2f;
            }
            else {
                circleSize = 0.05f;
            }

            // DrawSphere(point.position, point.sphereRadius, point.color);
            // DrawSphere(ppos,)
            Vector3 pposV3 = ToVector3({px,py,pz});
            // pposV3.x = pposV3.x/plotSize*windowWidth+windowWidth/2;
            // pposV3.y = pposV3.y/plotSize*windowWidth+windowWidth/2;
            // pposV3.z = pposV3.z/plotSize*windowWidth+windowWidth/2;
            Color circleColor = WHITE;
            if (i == 0) {
                circleColor = GREEN;
            }
            else if (i == particles.size()-2) {
                circleColor = YELLOW;
            }
            else if (i == particles.size()-1) {
                circleColor = RED;
            }
            DrawSphere(pposV3,circleSize,circleColor);
            // DrawCircle((int)px, (int)py, circleSize, WHITE);
        }

        EndMode3D();

        EndDrawing();
        // std::cout << particles[0].getPos()[0]/plotSize*800 << "\n";
        update();


    }

    CloseWindow();
    // while (frame < frames) {
    //     std::cout << particles[0].getPos()[0] << "\n";
    //     update();
    //     frame++;
    //     std::this_thread::sleep_for(std::chrono::seconds(1));
    // }

}