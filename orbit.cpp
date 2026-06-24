#include "orbit.h"

double wrapAngle(double angle) {
    angle = std::fmod(angle, 2.0 * pi);

    if (angle < 0.0) {
        angle += 2.0 * pi;
    }

    return angle;
}

double solveKepler(double meanAnomaly, double eccentricity) {
    double eccentricAnomaly = meanAnomaly;

    for (int iteration = 0; iteration < 20; ++iteration) {
        double f = eccentricAnomaly - eccentricity * std::sin(eccentricAnomaly) - meanAnomaly;

        double derivative = 1.0 - eccentricity * std::cos(eccentricAnomaly);

        double correction = f / derivative;
        eccentricAnomaly -= correction;

        if (std::abs(correction) < 1e-12) {
            break;
        }
    }

    return eccentricAnomaly;
}

ParticleState stateFromKeplerianOrbit(
    const OrbitParameters& orbit,
    double centralMass,
    double bodyMass,
    double targetTimeSeconds
) {
    double mu = G * (centralMass + bodyMass);

    double meanMotion = std::sqrt(
        mu / std::pow(orbit.semiMajorAxis, 3)
    );

    double elapsedTime = targetTimeSeconds - orbit.epochSeconds;

    double meanAnomaly = wrapAngle(
        orbit.meanAnomalyAtEpoch + meanMotion * elapsedTime
    );

    double eccentricAnomaly = solveKepler(
        meanAnomaly,
        orbit.eccentricity
    );

    double denominator = 1.0 - orbit.eccentricity * std::cos(eccentricAnomaly);

    double cosTrueAnomaly = (std::cos(eccentricAnomaly) - orbit.eccentricity) / denominator;

    double sinTrueAnomaly = (std::sqrt(1.0 - orbit.eccentricity * orbit.eccentricity) * std::sin(eccentricAnomaly)) / denominator;

    double trueAnomaly = std::atan2(sinTrueAnomaly, cosTrueAnomaly);

    double p = orbit.semiMajorAxis * (1.0 - orbit.eccentricity * orbit.eccentricity);

    double radius = p / (1.0 + orbit.eccentricity * std::cos(trueAnomaly));

    double xOrbital = radius * std::cos(trueAnomaly);
    double yOrbital = radius * std::sin(trueAnomaly);

    double velocityScale = std::sqrt(mu / p);

    double vxOrbital = -velocityScale * std::sin(trueAnomaly);

    double vyOrbital = velocityScale * (orbit.eccentricity + std::cos(trueAnomaly));

    double cosNode = std::cos(orbit.ascendingNode);
    double sinNode = std::sin(orbit.ascendingNode);

    double cosInclination = std::cos(orbit.inclination);
    double sinInclination = std::sin(orbit.inclination);

    double cosPeriapsis = std::cos(orbit.argumentOfPeriapsis);
    double sinPeriapsis = std::sin(orbit.argumentOfPeriapsis);


    double r11 = cosNode * cosPeriapsis - sinNode * sinPeriapsis * cosInclination;
    double r12 = -cosNode * sinPeriapsis- sinNode * cosPeriapsis * cosInclination;

    double r21 = sinNode * cosPeriapsis + cosNode * sinPeriapsis * cosInclination;
    double r22 = -sinNode * sinPeriapsis + cosNode * cosPeriapsis * cosInclination;

    double r31 = sinPeriapsis * sinInclination;
    double r32 = cosPeriapsis * sinInclination;

    Vect3 astronomicalPosition {
        r11 * xOrbital + r12 * yOrbital,
        r21 * xOrbital + r22 * yOrbital,
        r31 * xOrbital + r32 * yOrbital
    };

    Vect3 astronomicalVelocity {
        r11 * vxOrbital + r12 * vyOrbital,
        r21 * vxOrbital + r22 * vyOrbital,
        r31 * vxOrbital + r32 * vyOrbital
    };

    Vect3 simulationPosition {
        astronomicalPosition.x,
        astronomicalPosition.z,
        -astronomicalPosition.y
    };

    Vect3 simulationVelocity {
        astronomicalVelocity.x,
        astronomicalVelocity.z,
        -astronomicalVelocity.y
    };

    return {simulationPosition, simulationVelocity};
}