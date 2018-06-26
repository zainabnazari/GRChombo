/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4.hpp"

// For constraints calculation
#include "MatterConstraints.hpp"

// For tag cells
#include "OscillotonTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "OscillotonInitial.hpp"
#include "OscillotonPotential.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"
#include "SommerfeldBoundary.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, FILL_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // information about the csv file data
    int lines = 100002;
    double spacing = 0.01; // in r for the values

    std::array<double, 3> tmp = {0.0};
    std::vector<double> lapse_values, psi_values, Pi_values;

    std::string lapse_file(m_p.initial_data_prefix + "alpha001.csv");
    ifstream ifs0(lapse_file);

    std::string psi_file(m_p.initial_data_prefix + "psi001.csv");
    ifstream ifs1(psi_file);

    std::string Pi_file(m_p.initial_data_prefix + "Pi001.csv");
    ifstream ifs2(Pi_file);

    for (int i = 0; i < lines; ++i)
    {
        ifs0 >> tmp[0];
        ifs1 >> tmp[1];
        ifs2 >> tmp[2];

        lapse_values.push_back(tmp[0]);
        psi_values.push_back(tmp[1]);
        Pi_values.push_back(tmp[2]);
    }

    // First set everything to zero ... we don't want undefined values in
    // constraints etc
    m_state_new.setVal(0.);

    // Initial conditions for scalar field - Oscilloton
    BoxLoops::loop(OscillotonInitial(m_p.L, m_dx, m_p.oscilloton_center,
                                     spacing, lapse_values, psi_values,
                                     Pi_values),
                   m_state_new, m_state_new, FILL_GHOST_CELLS, disable_simd());
}

// Things to do before outputting a checkpoint file
void ScalarFieldLevel::prePlotLevel()
{
    fillAllGhosts();
    OscillotonPotential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    MatterConstraints<ScalarFieldWithPotential> constraints(scalar_field, m_dx,
                                                            m_p.G_Newton);
    BoxLoops::loop(constraints, m_state_new, m_state_new, SKIP_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Relaxation function for chi - this will eventually be done separately
    // with hdf5 as input
    if (m_time < m_p.relaxtime)
    {
        // Enforce trace free A_ij and positive chi and alpha
        BoxLoops::loop(
            make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()), a_soln,
            a_soln, FILL_GHOST_CELLS);

        // Calculate chi relaxation right hand side
        // Note this assumes conformal chi and Mom constraint trivially
        // satisfied  No evolution in other variables, which are assumed to
        // satisfy constraints per initial conditions
        OscillotonPotential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);
        ChiRelaxation<ScalarFieldWithPotential> relaxation(
            scalar_field, m_dx, m_p.relaxspeed, m_p.G_Newton);
        SetValue set_other_values_zero(0.0, Interval(c_h11, c_Mom3));
        auto compute_pack1 =
            make_compute_pack(relaxation, set_other_values_zero);
        BoxLoops::loop(compute_pack1, a_soln, a_rhs, SKIP_GHOST_CELLS);
    }
    else
    {
        // Enforce trace free A_ij and positive chi and alpha
        BoxLoops::loop(
            make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()), a_soln,
            a_soln, FILL_GHOST_CELLS);

        // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
        // We don't want undefined values floating around in the constraints so
        // zero these
        OscillotonPotential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);
        MatterCCZ4<ScalarFieldWithPotential> my_ccz4_matter(
            scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
            m_p.G_Newton);
        SetValue set_constraints_zero(0.0, Interval(c_rho, c_Mom3));
        auto compute_pack2 =
            make_compute_pack(my_ccz4_matter, set_constraints_zero);
        BoxLoops::loop(compute_pack2, a_soln, a_rhs, SKIP_GHOST_CELLS);
    }

    // Add in Sommerfeld BCs if on coarsest grid
    if (m_level == 0)
    {
        // Overwrite at boundaries with Sommerfeld Boundary conditions
        // (assumes boundary is only covered  by coarsest grid, no simd as
        // it isn't practical when you are only doing the three cells at
        // each border. As this only operates on the coarsest grid this
        // shouldn't be an issue.
        BoxLoops::loop(
            SommerfeldBoundary<ScalarFieldWithPotential>(m_p.L, m_dx, m_time),
            a_soln, a_rhs, SKIP_GHOST_CELLS, disable_simd());
    }
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, FILL_GHOST_CELLS);
}

// Specify if you want any plot files to be written, with which vars
void ScalarFieldLevel::specificWritePlotHeader(
    std::vector<int> &plot_states) const
{
    plot_states = {c_chi, c_phi, c_rho, c_Ham, c_Madm};
}

void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(OscillotonTaggingCriterion(m_dx, m_p.regrid_threshold_chi,
                                              m_p.regrid_threshold_Pi, m_level,
                                              m_p.extraction_params),
                   current_state, tagging_criterion);
}

void ScalarFieldLevel::specificPostTimeStep()
{
}
