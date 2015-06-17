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

