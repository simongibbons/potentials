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


    virtual ~Potential()
    {}

};

