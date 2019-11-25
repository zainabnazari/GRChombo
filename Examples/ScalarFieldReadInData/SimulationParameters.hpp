/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "CCZ4.hpp"
#include "ScalarPotential.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // for scalar field and initial data
        pp.load("G_Newton", G_Newton, 1.0);
        pp.load("scalar_mass", potential_params.scalar_mass);
    }

    // Initial data for matter and potential
    double G_Newton;
    ScalarPotential::params_t potential_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
