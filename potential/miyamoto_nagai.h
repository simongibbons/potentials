#pragma once

#include "static_cylindrical.h"

#include <consts.h>
#include <cmath>

class MiyamotoNagai : public PotentialStaticCylindrical<MiyamotoNagai> {

    double M;
    double a;
    double b;

public:

    MiyamotoNagai():M(1.0),a(2.0),b(0.2)
    {}

    MiyamotoNagai(double M, double a, double b):M(M),a(a),b(b)
    {}

    MiyamotoNagai(std::istringstream& iss)
    {
        iss >> M >> a >> b;
        if(M < 0) {
            throw PotentialError("M must be > 0 for MiyamotoNagai potential");
        }
        if(a < 0) {
            throw PotentialError("a must be > 0 for MiyamotoNagai potential");
        }
        if(b < 0) {
            throw PotentialError("b must be > 0 for MiyamotoNagai potential");
        }
    }

    std::string get_name() const
    {
        return "MiyamotoNagai";
    }

    double pot( double R, double z ) const {
        return -G*M / sqrt( R*R + pow(a + sqrt(z*z + b*b), 2.0) );
    }


    double ddR( double R, double z ) const {
        return G*M*R / pow( R*R + pow(a + sqrt(z*z + b*b), 2.0), 1.5 );
    }

    double ddz( double R, double z ) const {
        double q = sqrt( b*b + z*z );
        return G*M*z*( a + q ) / (q*pow( R*R + pow( a + q, 2.0 ), 1.5 ) );
    }

    double d2dR2( double R, double z ) const {
        double q2 = R*R + pow(a + sqrt(b*b + z*z), 2);
        return -G*M * (3*R*R*pow(q2, -5.0/2.0) - pow(q2, -3.0/2.0) );
    }

    double d2dz2( double R, double z ) const {
        double p = sqrt(b*b + z*z);
        double q = sqrt(R*R + pow(a + p, 2));

        return -G*M * ( 3*z*z * pow(a + p, 2) * pow(p, -2) * pow(q, -5) +
                       -z*z * pow(p, -2) * pow(q, -3) +
                       z*z*(a+p) * pow(p, -3) * pow(q, -3) +
                       -(a+p) * pow(p, -1) * pow(q, -3 ) );

    }

    double d2dRdz( double R, double z ) const {
        double p = sqrt(b*b + z*z);
        double denom = p * pow( R*R + pow(a + p, 2), 5.0/2.0 );

        return -3*G*M*R*z*(a + p) / denom;
    }

};

