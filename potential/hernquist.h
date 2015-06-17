#pragma once

#include "static_spherical.h"
#include <cmath>
#include <consts.h>


class Hernquist : public PotentialStaticSpherical<Hernquist>
{
    double M;
    double a;

public:

    Hernquist():M(1.0), a(1.0)
    {}

    Hernquist(double M, double a):M(M),a(a)
    {}

    Hernquist(std::istringstream& iss)
    {
        iss >> M >> a;
    }

    std::string get_name() const
    {
        return "Hernquist";
    }

    double pot(double r) const {
        return -G*M / (r + a);
    }

    double ddr(double r) const {
        return G*M * pow( r + a, -2);
    }

    double d2dr2(double r) const {
        return -2*G*M * pow(r + a, -3);
    }

};

