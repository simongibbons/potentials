#pragma once

#include "static_spherical.h"
#include <consts.h>

class Plummer : public PotentialStaticSpherical<Plummer> {

    double m, a;

public:
    Plummer():m(1.0), a(0.5) {}

    Plummer(std::istringstream& iss)
    {
        iss >> m >> a;
        if( m < 0 ) {
            throw PotentialError("m must be > 0 for Plummer Potential");
        }
        if( a < 0 ) {
            throw PotentialError("a must be > 0 for Plummer Potential");
        }
    }

    std::string get_name() const
    {
        return "Plummer";
    }

    double pot(double r) const
    {
        return -G*m / sqrt(r*r + a*a);
    }

    double ddr(double r) const
    {
        return G*m*r * pow(a*a + r*r, -3.0/2.0);
    }

    double d2dr2(double r) const
    {
        double rsoft = sqrt(a*a + r*r);
        return -G*m*( 3*r*r*pow(rsoft, -5) - pow(rsoft, -3) );
    }

};
