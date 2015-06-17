#pragma once

#include <unordered_map>
#include <memory>
#include <functional>

#include "base.h"
#include "composite.h"

#include "kepler.h"
#include "nfw.h"
#include "tf.h"
#include "miyamoto_nagai.h"
#include "hernquist.h"
#include "flattened_log.h"
#include "no_potential.h"
#include "plummer.h"

using pot_function_t = std::function<std::unique_ptr<Potential>(std::istringstream&)>;

template <typename T>
std::unique_ptr<Potential> potential_maker(std::istringstream& iss)
{
    return std::unique_ptr<Potential>(new T(iss));
}

const std::unordered_map<std::string, pot_function_t> potential_map = {
    {"Kepler", potential_maker<Kepler>},
    {"NFW", potential_maker<NFW>},
    {"TF", potential_maker<TFPotential>},
    {"MiyamotoNagai", potential_maker<MiyamotoNagai>},
    {"Hernquist", potential_maker<Hernquist>},
    {"FlattenedLog", potential_maker<FlattenedLog>},
    {"NoPotential", potential_maker<NoPotential>},
    {"Plummer", potential_maker<Plummer>}
};

