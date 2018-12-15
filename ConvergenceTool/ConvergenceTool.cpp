/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to LICENSE, in Chombo's root directory.
 */
#endif

// General includes:
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sys/time.h>

#include "parstream.H" //Gives us pout()
using std::endl;
#include "GRAMR.hpp"

#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "AMRInterpolator.hpp"
#include "DefaultLevelFactory.hpp"
#include "EmptyLevel.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "UserVariables.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

int runConvergenceTool(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    std::string in_string = argv[argc - 1];
    pout() << in_string << std::endl;
    char const *in_file = argv[argc - 1];
    GRParmParse pp(0, argv + argc, NULL, in_file);
    SimulationParameters sim_params(pp);

    // Setup the initial object (from restart_file checkpoint)
    GRAMR gr_amr;
    DefaultLevelFactory<EmptyLevel> empty_level_fact(gr_amr, sim_params);
    setupAMRObject(gr_amr, empty_level_fact);

    // set up number of extraction points and variables to extract
    int num_points = 5;
    pp.query("num_points", num_points);
    int var1 = c_chi;
    int var2 = c_phi;

    // The base value is often zero - so all levels are interpolated to the
    // corner point
    double base_dx = 0.0;
    pp.query("base_dx", base_dx);

    std::unique_ptr<double[]> var1_ptr{new double[num_points]};
    std::unique_ptr<double[]> var2_ptr{new double[num_points]};
    std::unique_ptr<double[]> interp_x{new double[num_points]};
    std::unique_ptr<double[]> interp_y{new double[num_points]};
    std::unique_ptr<double[]> interp_z{new double[num_points]};

    // set up the extraction point locations - evenly spaced on grid
    for (int ir = 0; ir < num_points; ++ir)
    {
                interp_x[ir] = 0.5 * base_dx + 0.5*sim_params.L * ir / num_points;
                interp_y[ir] = 0.5 * base_dx + 0.5*sim_params.L * ir / num_points;
                interp_z[ir] = 0.5 * base_dx + 0.5*sim_params.L * ir / num_points;
    }

    // set up the query
    InterpolationQuery query(num_points);
    query.setCoords(0, interp_x.get())
        .setCoords(1, interp_y.get())
        .setCoords(2, interp_z.get())
        .addComp(var1, var1_ptr.get())
        .addComp(var2, var2_ptr.get());

    // now the interpolator object, can choose up to 4th order
    AMRInterpolator<Lagrange<2>> interpolator(gr_amr, sim_params.origin, sim_params.dx, 2);

    // now loop over chk files
    for (int ifile = sim_params.start_file; ifile < sim_params.num_files;
         ifile++)
    {
        pout() << " loop number " << ifile << endl;
        // set up the file from next checkpoint
        std::ostringstream current_file;
        current_file << std::setw(6) << std::setfill('0')
                     << ifile * sim_params.checkpoint_interval;
        std::string restart_file(sim_params.checkpoint_prefix + current_file.str() +
                                 ".3d.hdf5");
        pout() << " restart file name " << restart_file << endl;
        HDF5Handle handle(restart_file, HDF5Handle::OPEN_RDONLY);
        pout() << " setupforrestart " << endl;

        gr_amr.setupForRestart(handle);

        pout() << " close handle " << endl;

        handle.close();

        pout() << " refresh interpolator " << endl;

        // refresh the interpolator and execute the query
        interpolator.refresh();
        interpolator.interp(query);

        // only rank 0 does the write out
        int rank;
        MPI_Comm_rank(Chombo_MPI::comm, &rank);
        if (rank == 0)
        {
            // set up file names and component names
            char file_str[sizeof "my_extraction_000000.txt"];
            sprintf(file_str, "my_extraction_%06d.txt",
                    sim_params.checkpoint_interval * ifile);
            char comp_str1[20];
            sprintf(comp_str1, UserVariables::variable_names[var1]);
            char comp_str2[20];
            sprintf(comp_str2, UserVariables::variable_names[var2]);

            ofstream outfile;
            outfile.open(file_str);
            if (!outfile.is_open())
            {
                MayDay::Error(
                    "error opening output file for extraction points");
            }

            // header data
            outfile << "# file number : " << ifile << endl;
            outfile << "#" << std::setw(19) << "x";
            outfile << std::setw(20) << "y";
            outfile << std::setw(20) << "z";
            outfile << std::setw(20) << comp_str1; //<< endl;
            outfile << std::setw(20) << comp_str2 << endl;

            // data itself
            for (int ir = 0; ir < num_points; ++ir)
            {
                outfile << std::setw(20) << interp_x[ir];
                outfile << std::setw(20) << interp_y[ir];
                outfile << std::setw(20) << interp_z[ir];
                outfile << std::setw(20) << std::setprecision(9)
                        << var1_ptr[ir]; // << endl;
                outfile << std::setw(20) << std::setprecision(9)
                        << var2_ptr[ir] << endl;
            }
            outfile.close();
        }
    }

    gr_amr.conclude();

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runConvergenceTool(argc, argv);

    mainFinalize();
    return status;
}
