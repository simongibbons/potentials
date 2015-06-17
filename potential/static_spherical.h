#pragma once

#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <cmath>
#include "base.h"

template <typename Derived>
class PotentialStaticSpherical : public Potential {

public:
    double operator()(double x, double y, double z, double t = 0 ) const {
        return static_cast<const Derived* const>(this)->pot(sqrt(x*x + y*y + z*z));
    }

    double ddx(double x, double y, double z, double t = 0 ) const {
        double r = sqrt(x*x + y*y + z*z);
        return (x/r) * static_cast<const Derived* const>(this)->ddr(r);
    }

    double ddy(double x, double y, double z, double t = 0 ) const {
        double r = sqrt(x*x + y*y + z*z);
        return (y/r) * static_cast<const Derived* const>(this)->ddr(r);
    }

    double ddz(double x, double y, double z, double t = 0 ) const {
        double r = sqrt(x*x + y*y + z*z);
        return (z/r) * static_cast<const Derived* const>(this)->ddr(r);
    }

    double d2dr2(double x, double y, double z, double t = 0 ) const {
        return static_cast<const Derived* const>(this)->d2dr2( sqrt( x*x + y*y + z*z ) );
    }

    double d2dx2(double x, double y, double z, double t = 0 ) const {
        double r = sqrt(x*x + y*y + z*z);
        double ddr = static_cast<const Derived* const>(this)->ddr(r);
        double d2dr2 = static_cast<const Derived* const>(this)->d2dr2(r);

        return (ddr/r)*(1 - pow(x/r, 2) ) + pow(x/r,2) * d2dr2;
    }

    double d2dxdy(double x, double y, double z, double t = 0 ) const {
        double r = sqrt(x*x + y*y + z*z);
        double ddr = static_cast<const Derived* const>(this)->ddr(r);
        double d2dr2 = static_cast<const Derived* const>(this)->d2dr2(r);

        return x*y * pow(r,-2) * ( d2dr2 - ddr / r);
    }

    double d2dxdz(double x, double y, double z, double t = 0 ) const {
        double r = sqrt(x*x + y*y + z*z);
        double ddr = static_cast<const Derived* const>(this)->ddr(r);
        double d2dr2 = static_cast<const Derived* const>(this)->d2dr2(r);

        return x*z * pow(r,-2) * ( d2dr2 - ddr / r);
    }

    double d2dy2(double x, double y, double z, double t = 0) const {
        double r = sqrt(x*x + y*y + z*z);
        double ddr = static_cast<const Derived* const>(this)->ddr(r);
        double d2dr2 = static_cast<const Derived* const>(this)->d2dr2(r);

        return (ddr/r)*(1 - pow(y/r, 2) ) + pow(y/r,2) * d2dr2;
    }

    double d2dydz(double x, double y, double z, double t = 0) const {
        double r = sqrt(x*x + y*y + z*z);
        double ddr = static_cast<const Derived* const>(this)->ddr(r);
        double d2dr2 = static_cast<const Derived* const>(this)->d2dr2(r);

        return y*z * pow(r,-2) * ( d2dr2 - ddr / r);
    }

    double d2dz2(double x, double y, double z, double t = 0) const {
        double r = sqrt(x*x + y*y + z*z);
        double ddr = static_cast<const Derived* const>(this)->ddr(r);
        double d2dr2 = static_cast<const Derived* const>(this)->d2dr2(r);

        return (ddr/r)*(1 - pow(z/r, 2) ) + pow(z/r,2) * d2dr2;
    }

    void acc(double x, double y, double z, double t, Vec3D& force) const {
        double r = sqrt(x*x + y*y + z*z);
        double ddr = static_cast<const Derived* const>(this)->ddr(r);

        force.x = (x/r) * ddr;
        force.y = (y/r) * ddr;
        force.z = (z/r) * ddr;
    }

};

