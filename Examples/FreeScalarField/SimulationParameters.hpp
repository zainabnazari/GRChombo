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
#include "Potential.hpp"

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
        // Reads the problem specific params
        auto_read_params(pp);

        // Fill in the potential parameters
        potential_params.scalar_mass = scalar_mass;
    }

    void auto_read_params(GRParmParse &pp)
    {
        // Initial and SF data
        pp.load("initial_data_prefix", initial_data_prefix);
        pp.load("scalar_mass", scalar_mass);
	pp.load("G_Newton", G_Newton, 1.0);

        // for the regridding
        pp.load("regrid_threshold_phi", regrid_threshold_phi);
        pp.load("regrid_threshold_chi", regrid_threshold_chi);
    }

    // Initial data for potential
    double G_Newton;
    double scalar_mass, regrid_threshold_chi, regrid_threshold_phi;

    // location of data for the oscilloton profile in isotropic coords
    std::string initial_data_prefix;

    // Collection of parameters necessary for the problem
    Potential::params_t potential_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
