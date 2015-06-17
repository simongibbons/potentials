#pragma once

#include "static_spherical.h"

#include <cmath>
#include <consts.h>

class NFW : public PotentialStaticSpherical<NFW> {

    double rho0;
    double rs;
    double A;

public:

    NFW():rho0(0.0001),rs(10.0),A(4*pi*G*rho0*pow(rs,3))
    {}

    NFW(double rho0, double rs):rho0(rho0),rs(rs),A(4*pi*G*rho0*pow(rs,3))
    {}

    NFW(std::istringstream& iss)
    {
        iss >> rho0 >> rs;
        A = 4*pi*G*rho0*pow(rs,3);

        if(rho0 < 0) {
            throw PotentialError("rho0 must be > 0 for NFW Potential");
        }
        if(rs < 0 ) {
            throw PotentialError("rs must be > 0 for NFW Potential");
        }
    }

    std::string get_name() const
    {
        return "NFW";
    }

    double pot(double r) const {
        return -A*log( 1 + r/rs ) / r;
    }

    double ddr(double r) const {
        double X = r/rs;

        return A*( log(1+X) - X/(1+X) ) / pow(r, 2);
    }

    double d2dr2(double r) const {
        return A / pow(r, 3) * ( r*(3*r + 2*rs) / pow(r + rs, 2) -
                                       2*log(1 + r/rs) );
    }

};

