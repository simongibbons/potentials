%module potentials

%{
#include "potential/kepler.h"
#include "potential/nfw.h"
#include "potential/tf.h"
#include "potential/miyamoto_nagai.h"
#include "potential/hernquist.h"
#include "potential/plummer.h"
%}


%include "potential/base.h"
%include "potential/static_spherical.h"
%include "potential/static_cylindrical.h"

%template(PotentialStaticSphericalNFW) PotentialStaticSpherical<NFW>;
%include "potential/nfw.h"

%template(PotentialStaticSphericalKepler) PotentialStaticSpherical<Kepler>;
%include "potential/kepler.h"

%template(PotentialStaticSphericalTF) PotentialStaticSpherical<TFPotential>;
%include "potential/tf.h"

%template(PotentialStaticCylindricalMN) PotentialStaticCylindrical<MiyamotoNagai>;
%include "potential/miyamoto_nagai.h"

%template(PotentialStaticSphericalHernquist) PotentialStaticSpherical<Hernquist>;
%include "potential/hernquist.h"

%template(PotentialStaticSphericalPlummer) PotentialStaticSpherical<Plummer>;
%include "potential/plummer.h"


%pythoncode %{
import numpy as np
def make_hessian(potential, x, y, z, t=0.0):
    rv = np.empty((3,3))
    rv[0,0] = potential.d2dx2(x,y,z,t)
    rv[1,1] = potential.d2dy2(x,y,z,t)
    rv[2,2] = potential.d2dz2(x,y,z,t)

    rv[0,1] = potential.d2dxdy(x,y,z,t)
    rv[0,2] = potential.d2dxdz(x,y,z,t)
    rv[1,2] = potential.d2dydz(x,y,z,t)

    rv[1,0] = rv[0,1]
    rv[2,0] = rv[0,2]
    rv[2,1] = rv[1,2]

    return rv

Potential.hessian = make_hessian
%}
