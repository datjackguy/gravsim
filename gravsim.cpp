//Compile Portable -- g++ gravsim.cpp -o gravsim.exe -I include -L lib -lraylib -lgdi32 -lwinmm -static-libgcc -static-libstdc++ -static
//2 -- g++ gravsim.cpp -o gravsim.exe -I include -L lib -static -static-libgcc -static-libstdc++ -lraylib -lopengl32 -lgdi32 -lwinmm
//3 Optimised g++ -O2 gravsim.cpp -o gravsim.exe -I include -L lib -lraylib -lgdi32 -lwinmm -static-libgcc -static-libstdc++ -static

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
const double M0 = 1.989e30;
bool enableCentralMass = true;

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

struct Vect3 {
    double x = 0;
    double y = 0;
    double z = 0;
    double magnitude() {
        return pythagoras(x,y,z);
    }
    void set(double x2, double y2, double z2) {
        x = x2;
        y = y2;
        z = z2;
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
    double rho = randomDouble(maxRad/20,maxRad);
    double theta = randomDouble(0.0,pi);
    double phi = randomDouble(0.0, 2*pi);
    double x = rho * std::sin(theta) * std::cos(phi);
    double y = 0.4 * rho * std::sin(theta) * std::sin(phi);
    double z = rho * std::cos(theta);
    Vect3 vector;
    vector.x = x;
    vector.y = y;
    vector.z = z;

    return vector;
}



//Compute the velocity of a particle at x,y,z
Vect3 globularVel(double x, double y, double z, double totalMass) {
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

    Vect3 vector;
    vector.x = velx;
    vector.y = vely;
    vector.z = velz;

    return vector;

}


//Particle Class
class Particle {
private:
    double mass;
    // double posx;
    // double posy;
    // double posz;
    Vect3 pos;
    // double velx;
    // double vely;
    // double velz;
    Vect3 vel;
    double pvelx; //Previous velocity -- Currently unused
    double pvely;
    double pvelz;
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
        force.x += sign*fnew.x;
        force.y += sign*fnew.y;
        force.z += sign*fnew.z;
    }

    // Calculate the new position
    void calcPos(double dt) {
        pos.x += vel.x * dt;
        pos.y += vel.y * dt;
        pos.z += vel.z * dt;
    }

    // Calculate the new velocity
    void calcVel(double dt) {
        vel.x += force.x * dt / mass;
        vel.y += force.y * dt / mass;
        vel.z += force.z * dt / mass;
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

    double getMass() {
        return mass;
    }

    void setMass(double newMass) {
        mass = newMass;
    }

    Vect3 getForce() {
        return force;
    }

    Vect3 getPos() {
        return pos;
    }

    void setPos(Vect3 newPos) {
        pos.x = newPos.x;
        pos.y = newPos.y;
        pos.z = newPos.z;
    }

    Vect3 getVel() {
        return vel;
    }

    void setVel(Vect3 newVel) {
        vel = newVel;
    }
};

// Calculate the distance between two particles - Move to particle class?
double distanceCalc(Particle &p1, Particle &p2) {
    Vect3 pos1 = p1.getPos();
    Vect3 pos2 = p2.getPos();
    double dx = (pos1.x-pos2.x);
    double dy = (pos1.y-pos2.y);
    double dz = (pos1.z-pos2.z);

    double distance = pythagoras(dx,dy,dz);

    return distance;
}

// Calculate the force between two particles - Move to particle class?
Vect3 forceCalc(Particle &p1, Particle &p2) {
    double distance = distanceCalc(p1,p2);
    double mass1 = p1.getMass();
    double mass2 = p2.getMass();
    Vect3 pos1 = p1.getPos();
    Vect3 pos2 = p2.getPos();

    Vect3 force;
    force.x = -G*mass1*mass2*(pos1.x-pos2.x)/pow((distance*distance)+(epsilon*epsilon),1.5);
    force.y = -G*mass1*mass2*(pos1.y-pos2.y)/pow((distance*distance)+(epsilon*epsilon),1.5);
    force.z = -G*mass1*mass2*(pos1.z-pos2.z)/pow((distance*distance)+(epsilon*epsilon),1.5);

    return force;
}

//Simulation Parameters
double timestep = 0.5*10*24*60*60; //s
double runtime = 100*365.25*24*60*60; //s
int frame = 0;
int frames = 10;
std::vector<double> times;
double currentTime = 0.0;
std::vector<Particle> particles;
int particleN = 300;
double plotSize = 2e15; //m
double totalMass = 0; //kg
std::vector<std::vector<int>> particlePairs;
double maxMass = 0;
double minMass = 1e64;
int maxTrailLength = 20;
std::vector<Color> colourList = {RED,GREEN,BLUE,ORANGE,PURPLE,GREEN,MAGENTA,SKYBLUE};
bool randomColours = false;

Color randomColour() {
    int colouri = randomInteger(0,colourList.size()-1);
    return colourList[colouri];
}


//Calculate total mass within the radius of a particle
double enclosedMass(Particle &particle) {
    Vect3 ppos = particle.getPos();
    double x = ppos.x;
    double y = ppos.y;
    double z = ppos.z;
    double radius = pythagoras(x,y,z);
    double mass = 0;
    for (int p = 0; p < particles.size(); p++) {
        Vect3 currentPos = particles[p].getPos();
        double cx = currentPos.x;
        double cy = currentPos.y;
        double cz = currentPos.z;
        double currentRadius = pythagoras(cx,cy,cz);
        if (currentRadius < radius) {
            mass += particles[p].getMass();
        }
    }
    return mass;
}

//Create all particle objects
void initialise_particles() {
    std::vector<double> masses = generate_masses(particleN);
    for (int i = 0; i < particleN; i++) {
        totalMass = totalMass + masses[i];
    }

    for (int i = 0; i < particleN; i++) {
        Vect3 pos = randomSph(plotSize/2);
        double xpos = pos.x;
        double ypos = pos.y;
        double zpos = pos.z;
        
        Vect3 vel = globularVel(xpos,ypos,zpos, totalMass);
        double xvel = vel.x;
        double yvel = vel.y;
        double zvel = vel.z;

        particles.push_back(Particle(masses[i],pos,vel));
    }
    if (enableCentralMass) {
        Vect3 zero;
        zero.set(0,0,0);
        particles.push_back(Particle(5e38, zero, zero));
    }
    // particles.push_back(Particle(5e38, 1e16, 0, 1e16, -7e6, 5e4, -6e6));
    for (int i = 0; i < particleN; i++) {
        double encMass = enclosedMass(particles[i]);
        Vect3 currentPos = particles[i].getPos();
        double xpos = currentPos.x;
        double ypos = currentPos.y;
        double zpos = currentPos.z;

        Vect3 vel = globularVel(xpos,ypos,zpos,encMass);
        double xvel = vel.x;
        double yvel = vel.y;
        double zvel = vel.z;

        particles[i].setVel(vel);
    }
    for (int i = 0; i < particles.size(); i++) {
        if (particles[i].getMass() > maxMass) {
            maxMass = particles[i].getMass();
        }
        else if (particles[i].getMass() < minMass) {
            minMass = particles[i].getMass();
        }
    }
}


// Create list containing the indices of each pair of particles
void create_particle_pairs() {
    for (int i = 0; i < particles.size(); i++) {
        for (int j = i+1; j < particles.size(); j++) {
            std::vector<int> pair = {i,j}; 
            particlePairs.push_back(pair);
        }
    }
}

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

//Update the simulation for each timestep
void update() {
    for (int p = 0; p < particlePairs.size(); p++) {
        Vect3 force = forceCalc(particles[particlePairs[p][0]],particles[particlePairs[p][1]]);

        particles[particlePairs[p][0]].addForce(force);
        particles[particlePairs[p][1]].addForce(force, -1);
    }
    for (int i = 0; i < particles.size(); i++) {
        particles[i].calcPos(timestep);
        particles[i].calcVel(timestep);
        particles[i].resetForce();

        Vect3 pos = particles[i].getPos();
        float renderScale = 10.0f / plotSize;

        Vector3 trailPoint = {
            static_cast<float>(pos.x * renderScale),
            static_cast<float>(pos.y * renderScale),
            static_cast<float>(pos.z * renderScale)
        };

        particles[i].trail.push_back(trailPoint);

        if (particles[i].trail.size() > maxTrailLength) {
            particles[i].trail.erase(particles[i].trail.begin());
        }
    }
}

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

    // Convert yaw, pitch and distance into a 3D camera position.
    camera.position.x = camera.target.x + distance * cosf(pitch) * sinf(yaw);
    camera.position.y = camera.target.y + distance * sinf(pitch);
    camera.position.z = camera.target.z + distance * cosf(pitch) * cosf(yaw);
}

//Allow the user to move a particle mid-simulation with WASD
void UpdateParticlePos(Particle &particle) {
    Vect3 ppos = particle.getPos();
    double x = ppos.x;
    double y = ppos.y;
    double z = ppos.z;
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
float GetVisualRadius(double mass) {

    //Change to median!!!?

    float averageMass = (totalMass/particles.size());

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

bool drawTrails = true;


void OpenSetupGUI() {
    //Open Settings GUI
    InitWindow(500, 500, "Gravity Simulation Settings");
    SetTargetFPS(60);

    IntSlider nSlider = IntSlider(particleN, "20", "400", 20, 400, "Total Particles");
    IntSlider trailSlider = IntSlider(maxTrailLength, "1", "100", 0.0f, 100.0f, "Trail Length");
    Tickbox centralMassTick = Tickbox(enableCentralMass, "Central Mass");
    LogSlider plotsizeSlider = LogSlider(plotSize, "2e14", "2e16", 2e14, 2e16, "Cluster Size (m)");
    bool startPressed = false;
    Tickbox drawTrailsTick = Tickbox(drawTrails, "Enable Trails");
    Tickbox randomColoursTick = Tickbox(randomColours, "Random Colours");
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

        if (drawTrails) {
            trailSlider.draw(nextYpos);
            nextYpos += trailSlider.height;
        }

        // GuiCheckBox(Rectangle{40, 320, 20, 20}, "Random Colours", &randomColours);
        randomColoursTick.draw(nextYpos);
        nextYpos += randomColoursTick.height;

        centralMassTick.draw(nextYpos);
        nextYpos += centralMassTick.height;

        if (GuiButton(Rectangle{40, 400, 140, 35}, "Start")) {
            startPressed = true;
            }

        EndDrawing();
    }

    CloseWindow();
}


int main() {

    OpenSetupGUI();

    //Start Simulation
    initialise_particles();
    create_particle_pairs();
    // momentum_centre();

    if (randomColours) {
        for (int i = 0; i < particles.size(); i++) {
            Color colour = randomColour();
            particles[i].setColour(colour);
        }
    }
    particles[particles.size()-1].setColour(YELLOW);

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


    while (!WindowShouldClose()) {
        float dt = GetFrameTime();
        float time = GetTime();

        UpdateOrbitCamera(camera, cameraYaw, cameraPitch, cameraDistance);

        UpdateParticlePos(particles[particles.size()-2]);

        UpdateParticleMass(particles[particles.size()-2]);

        BeginDrawing();

        ClearBackground(BLACK);

        BeginMode3D(camera);

        DrawGrid(20, 1.0f);


        DrawLine3D({ -10.0f, 0.0f, 0.0f }, { 10.0f, 0.0f, 0.0f }, RED);
        DrawLine3D({ 0.0f, -10.0f, 0.0f }, { 0.0f, 10.0f, 0.0f }, GREEN);
        DrawLine3D({ 0.0f, 0.0f, -10.0f }, { 0.0f, 0.0f, 10.0f }, BLUE);

        for (int i = 0; i < particles.size(); i++) {
            Vect3 ppos = particles[i].getPos();
            float renderScale = 10.0f / plotSize;

            float px = static_cast<float>(ppos.x * renderScale);
            float py = static_cast<float>(ppos.y * renderScale);
            float pz = static_cast<float>(ppos.z * renderScale);

            // float circleSize = particles[i].getMass()/1.989e30; 
            float circleSize;// = std::log(particles[i].getMass())/std::log((1.989e34))*2;
            // circleSize = 0.3f * cbrt(particles[i].getMass()/totalMass);
            circleSize = GetVisualRadius(particles[i].getMass());

            // if (particles[i].getMass() > 1e37) {
            //     circleSize = 0.4f;
            // }
            Vector3 pposV3 = ToVector3({px,py,pz});

            Color circleColor = particles[i].getColour();
            if (i == 0) {
                circleColor = GREEN;
            }
            else if (i == particles.size()-1) {
                circleColor = YELLOW;
            }
            // else if (i == particles.size()-1) {
            //     circleColor = RED;
            // }
            if (drawTrails) {
                DrawTrail(particles[i]);
            }
            DrawSphere(pposV3,circleSize,circleColor);
        }

        EndMode3D();

        EndDrawing();
        update();
    }

    CloseWindow();

}