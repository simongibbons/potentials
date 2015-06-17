#pragma once

#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <cmath>
#include "base.h"

template <typename Derived>
class PotentialStaticCylindrical : public Potential {

public:
    double operator()(double x, double y, double z, double t = 0 ) const {
        return static_cast<const Derived* const>(this)->pot(sqrt(x*x + y*y), z);
    }

    double ddx(double x, double y, double z, double t = 0 ) const {
        double R = sqrt(x*x + y*y);
        return R > 0 ? (x/R) * static_cast<const Derived* const>(this)->ddR(R,z) : 0;
    }

    double ddy(double x, double y, double z, double t = 0 ) const {
        double R = sqrt(x*x + y*y);
        return R > 0 ? (y/R) * static_cast<const Derived* const>(this)->ddR(R,z) : 0;
    }

    double ddz(double x, double y, double z, double t = 0 ) const {
        return static_cast<const Derived* const>(this)->ddz(sqrt(x*x + y*y),z);
    }

    double d2dr2(double x, double y, double z, double t = 0 ) const {
        double R = sqrt(x*x + y*y);
        double r2 = x*x + y*y + z*z;
        return (1/r2) * ( R*R * static_cast<const Derived* const>(this)->d2dR2(R,z)+
                          R*z * static_cast<const Derived* const>(this)->d2dRdz(R,z)+
                          z*z * static_cast<const Derived* const>(this)->d2dz2(R,z) );
    }


    double d2dx2(double x, double y, double z, double t = 0 ) const {
        double R = sqrt(x*x + y*y);
        double ddR = static_cast<const Derived* const>(this)->ddR(R,z);
        double d2dR2 = static_cast<const Derived* const>(this)->d2dR2(R,z);
        if( R > 0 ) {
            return pow(R,-2) * (pow(y,2)*ddR/R + pow(x,2)*d2dR2);
        } else {
            return 0.0;
        }
    }

    double d2dxdy(double x, double y, double z, double t = 0 ) const {
        double R = sqrt(x*x + y*y);
        return x*y/pow(R,2) * ( static_cast<const Derived* const>(this)->d2dR2(R,z) -
                                static_cast<const Derived* const>(this)->ddR(R,z)/R );
    }

    double d2dxdz(double x, double y, double z, double t = 0 ) const {
        double R = sqrt(x*x + y*y);
        return (x/R) * static_cast<const Derived* const>(this)->d2dRdz(R,z);
    }

    double d2dy2(double x, double y, double z, double t = 0 ) const {
        double R = sqrt(x*x + y*y);
        double ddR = static_cast<const Derived* const>(this)->ddR(R,z);
        double d2dR2 = static_cast<const Derived* const>(this)->d2dR2(R,z);
        if( R > 0 ) {
            return pow(R,-2) * (pow(x,2)*ddR/R + pow(y,2)*d2dR2);
        } else {
            return 0.0;
        }
    }

    double d2dydz(double x, double y, double z, double t = 0 ) const {
        double R = sqrt(x*x + y*y);
        return (y/R) * static_cast<const Derived* const>(this)->d2dRdz(R,z);
    }

    double d2dz2(double x, double y, double z, double t = 0 ) const {
        return static_cast<const Derived* const>(this)->d2dz2(sqrt(x*x + y*y),z);
    }

    void acc(double x, double y, double z, double t, Vec3D& force) const {
        const double R = sqrt(x*x + y*y);
        const double ddR = static_cast<const Derived* const>(this)->ddR(R,z);

        force.x = R > 0 ? (x/R) * ddR : 0;
        force.y = R > 0 ? (y/R) * ddR : 0;
        force.z = static_cast<const Derived* const>(this)->ddz(R,z);
    } 

};

