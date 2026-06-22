#pragma once

#include "core.h"

double pythagoras(double x, double y, double z);
double pythagoras(std::vector<double> vector);
double pythagoras(Vect3 vector);

double randomDouble(double min, double max);
int randomInteger(int min, int max);
int randomSign();

double clampDouble(double value, double low, double high);

double initial_mass_function(double m);
std::vector<double> generate_masses(int N);

Vect3 randomAxis();

Vect3 randomSph(double maxRad);
Vect3 generateDisk(double maxRad, double flattening, Vect3 axis);


Vect3 globularVel(Vect3 pos, double totalMass);
Vect3 diskVelocity(Vect3 pos, double orbitMass, Vect3 axis);
Vect3 diskVelocity(Vect3 pos, double orbitMass, Vect3 axis, double softening);