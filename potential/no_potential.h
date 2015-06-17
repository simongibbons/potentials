#pragma once

#include <sstream>

#include "base.h"

#pragma GCC diagnostic ignored "-Wunused-parameter"

class NoPotential : public Potential {
public:
    NoPotential()
    {}

    NoPotential(std::istringstream& iss)
    {}

    std::string get_name() const
    {
        return "NoPotential";
    }

    double operator()(double x, double y, double z, double t = 0) const {
        return 0.0;
    }

    double d2dr2(double x, double y, double z, double t = 0) const {
        return 0.0;
    }

};

