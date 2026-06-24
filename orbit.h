#pragma once

#include "core.h"

double wrapAngle(double angle);
double solveKepler(double meanAnomaly, double eccentricity);
ParticleState stateFromKeplerianOrbit(
    const OrbitParameters& orbit,
    double centralMass,
    double bodyMass,
    double targetTimeSeconds
);