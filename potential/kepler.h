#pragma once

#include "static_spherical.h"
#include <consts.h>

class Kepler : public PotentialStaticSpherical<Kepler> {

    double M;

public:
    Kepler():M(1.0) {}
    Kepler(double M):M(M) {}

    Kepler(std::istringstream& iss)
    {
        iss >> M;
        if(M < 0) {
            throw PotentialError("M must be > 0 for Kepler potential");
        }
    }

    std::string get_name() const
    {
        return "Kepler";
    }

    double pot(double r) const {
        return -G*M/r;
    }

    double ddr(double r) const {
        return G*M * pow(r, -2.0);
    }

    double d2dr2(double r) const {
        return -2*G*M * pow(r, -3.0);
    }

};
