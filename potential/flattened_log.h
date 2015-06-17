#pragma once

#include <cmath>

#include "static_cylindrical.h"
#include <consts.h>

class FlattenedLog : public PotentialStaticCylindrical<FlattenedLog> {

    double v02;
    double q;
    double q2;

public:

    FlattenedLog():v02(pow(220.0,2)), q(0.9), q2(q*q)
    {}

    FlattenedLog(double v0, double q):v02(v0*v0), q(q), q2(q*q)
    {}

    FlattenedLog(std::istringstream& iss)
    {
        double v0, q;
        iss >> v0 >> q;
        v02 = v0*v0;
        q2 = q*q;

        if(v0 < 0) {
            throw PotentialError("v0 must be > 0 for FlattenedLog potential");
        }
    }

    std::string get_name() const
    {
        return "FlattenedLog";
    }

    double pot(double R, double z) const {
        return 0.5*v02*log( R*R + z*z / q2 );
    }

    double ddR(double R, double z) const {
        return R*v02 / ( R*R + z*z / q2 );
    }

    double ddz(double R, double z) const {
        return z*v02 / (q2 * (R*R + z*z / q2) );
    }

    double d2dR2(double R, double z) const {
        return q2*v02*(z*z - q2*R*R) / pow(q2*R*R + z*z, 2);
    }

    double d2dRdz(double R, double z) const {
        return -2*v02*R*z / (q2 * pow(R*R + z*z/q2, 2));
    }

    double d2dz2(double R, double z) const {
        return v02*(q*R - z)*(q*R + z) / pow(q2*R*R + z*z, 2);
    }

};

