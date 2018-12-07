/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    // For the Interpolator test we don't need many parameters
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        pp.get("verbosity", verbosity);
        // Grid setup
        pp.get("L", L);
        pp.getarr("isPeriodic", isPeriodic, 0, SpaceDim);
        pp.get("num_ghosts", num_ghosts);
        pp.get("num_files", num_files);
        pp.get("start_file", start_file);
        pp.get("checkpoint_interval", checkpoint_interval);
        pp.get("chk_prefix", chk_prefix);
        pp.get("dt_multiplier", dt_multiplier);
        pp.get("tag_buffer_size", tag_buffer_size);
    }
    int verbosity;
    std::string chk_prefix;
    Real L; // Physical sidelength of the grid
    int num_ghosts, num_files, start_file, checkpoint_interval;
    std::vector<bool> isPeriodic;
    double regrid_threshold = 0.5;
    int tag_buffer_size;
    bool ignore_checkpoint_name_mismatch = false;
    double dt_multiplier; // Doesn't matter for this test
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
