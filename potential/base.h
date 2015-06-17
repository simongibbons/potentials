#pragma once

#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <iostream>
#include <string>
#include <exception>

#include <util/derivatives.h>
#include <util/vec3d.h>

#include <armadillo>

class PotentialError : public std::exception
{
    std::string error_string;
public:
    PotentialError(const std::string& s):error_string("POTENTIAL ERROR: " + s)
    {}

    virtual const char* what() const throw()
    {
        return error_string.c_str();
    }
};

class Potential {

public:
    //virtual std::string get_name() const = 0;

    virtual double operator()(double x, double y, double z, double t ) const = 0;

    virtual double ddx(double x, double y, double z, double t ) const
    {
        return derivative(x, 0.001,
                          [=](double x){ return (*this)(x,y,z,t); });
    }

    virtual double ddy(double x, double y, double z, double t ) const
    {
        return derivative(y, 0.001,
                          [=](double y){ return (*this)(x,y,z,t); });
    }

    virtual double ddz(double x, double y, double z, double t ) const
    {
        return derivative(z, 0.001,
                          [=](double z){ return (*this)(x,y,z,t); });
    }

    virtual double d2dx2(double x, double y, double z, double t) const
    {
        return derivative(x, 0.001,
                          [=](double x){ return this->ddx(x,y,z,t); });
    }

    virtual double d2dxdy(double x, double y, double z, double t) const
    {
        return derivative(x, 0.001,
                          [=](double x){ return this->ddy(x,y,z,t); });
    }

    virtual double d2dxdz(double x, double y, double z, double t) const
    {
        return derivative(x, 0.001,
                          [=](double x){ return this->ddz(x,y,z,t); });
    }

    virtual double d2dy2(double x, double y, double z, double t) const
    {
        return derivative(y, 0.001,
                          [=](double y){ return this->ddy(x,y,z,t); });
    }

    virtual double d2dydz(double x, double y, double z, double t) const
    {
        return derivative(y, 0.001,
                          [=](double y){ return this->ddz(x,y,z,t); });
    }

    virtual double d2dz2(double x, double y, double z, double t) const
    {
        return derivative(z, 0.001,
                          [=](double z){ return this->ddz(x,y,z,t); });
    }

    virtual double d2dr2(double x, double y, double z, double t ) const = 0;

    virtual void acc(double x, double y, double z, double t, Vec3D& force) const
    {
        force.x = ddx(x, y, z, t);
        force.y = ddy(x, y, z, t);
        force.z = ddz(x, y, z, t);
    }

    // Compute the hessian of the potential at a given point
    arma::mat33 hessian(double x, double y, double z, double t) const
    {
        arma::mat33 m = arma::zeros<arma::mat>(3,3);

        m(0,0) = d2dx2(x,y,z,t);
        m(0,1) = d2dxdy(x,y,z,t);
        m(0,2) = d2dxdz(x,y,z,t);
        m(1,1) = d2dy2(x,y,z,t);
        m(1,2) = d2dydz(x,y,z,t);
        m(2,2) = d2dz2(x,y,z,t);

        m(1,0) = m(0,1);
        m(2,0) = m(0,2);
        m(2,1) = m(1,2);

        return m;
    }

    virtual ~Potential()
    {}

};

