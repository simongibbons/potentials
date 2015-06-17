#pragma once

#include "static_spherical.h"
#include <consts.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

class TFPotential : public PotentialStaticSpherical<TFPotential> {

    double v02, rs, alpha;

public:
    TFPotential():v02(220.0*220.0), rs(5.0), alpha(0.6)
    {}

    TFPotential(double v0, double rs, double alpha):
        v02(v0*v0),
        rs(rs),
        alpha(alpha)
    {}

    std::string get_name() const
    {
        return "TFPotential";
    }

    TFPotential(std::istringstream& iss)
    {
        double v0;
        iss >> v0 >> rs >> alpha;
        v02 = v0*v0;

        if(v0 < 0) {
            throw PotentialError("v0 must be > 0 for TF potential");
        }
        if(rs < 0) {
            throw PotentialError("rs must be > 0 for TF potential");
        }
        if(alpha < 0 || alpha > 1) {
            throw PotentialError("0 < alpha < 1 for TF potential");
        }
    }

    double pot(double r) const {
        auto old_handler = gsl_set_error_handler_off(); 
        
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);

        double result, error;
        double params[3] = {v02, rs, alpha};

        gsl_function F;
        F.function = &pot_integrand;
        F.params = &params;

        gsl_integration_qagiu(&F, r, 1e-13, 1e-13, 5000, w, &result, &error);

        gsl_integration_workspace_free(w);

        gsl_set_error_handler(old_handler);

        return result;
    }

    double ddr(double r) const {
        return v02 * pow(1 + pow(r/rs, 2), -alpha/2) / r;
    }

    double d2dr2(double r) const {
        return -v02 * pow( 1 + pow(r/rs, 2.0), -0.5*alpha ) *
                          ((1 + alpha)*r*r + rs*rs) /
                          (r*r*(r*r + rs*rs));
    }

private:
    static double pot_integrand(double r, void *params)
    {
        double *p = static_cast<double*>(params);
        double v02 = p[0];
        double rs = p[1];
        double alpha = p[2];

        return -v02 / (r * pow(1 + pow(r/rs, 2), 0.5*alpha));
    }

};

