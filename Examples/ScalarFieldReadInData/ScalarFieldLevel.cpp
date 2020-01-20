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
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "ChiRelaxation.hpp"
#include "ComputePack.hpp"
#include "MinkowskiMetric.hpp"
#include "ScalarField.hpp"
#include "ScalarPotential.hpp"
#include "SetValue.hpp"
#include "VarDataInterpolation.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // First set everything to zero, then set Minkowski
    BoxLoops::loop(make_compute_pack(SetValue(0.0), MinkowskiMetric()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Then read in the data for interpolation - need to provide the precision
    // for the coord and data values in the file
    SmallDataIO input_file_Pi("TestData/field.txt", m_dt, m_time,
                              m_restart_time, SmallDataIO::READ, true);
    std::vector<std::array<double, CH_SPACEDIM + 1>> Pi_data;
    input_file_Pi.get_data_array(Pi_data);

    SmallDataIO input_file_chi("TestData/dfield.txt", m_dt, m_time,
                               m_restart_time, SmallDataIO::READ, true);
    std::vector<std::array<double, CH_SPACEDIM + 1>> chi_data;
    input_file_chi.get_data_array(chi_data);

    // Then interpolate data
    BoxLoops::loop(make_compute_pack(
                       VarDataInterpolation(m_dx, c_Pi, m_p.center, Pi_data),
                       VarDataInterpolation(m_dx, c_chi, m_p.center, chi_data)),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
                   disable_simd());
}

// Things to do before outputting a checkpoint file
void ScalarFieldLevel::preCheckpointLevel()
{
    fillAllGhosts();
    ScalarPotential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(
                       scalar_field, m_dx, m_p.G_Newton),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    // We don't want undefined values floating around in the constraints so
    // zero these
    ScalarPotential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    MatterCCZ4<ScalarFieldWithPotential> my_ccz4_matter(
        scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
        m_p.G_Newton);
    SetValue set_constraints_zero(0.0, Interval(c_Ham, c_Mom3));
    auto compute_pack2 =
        make_compute_pack(my_ccz4_matter, set_constraints_zero);
    BoxLoops::loop(compute_pack2, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

// Specify if you want any plot files to be written, with which vars
void ScalarFieldLevel::specificWritePlotHeader(
    std::vector<int> &plot_states) const
{
    plot_states = {c_phi, c_chi};
}

void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(
        FixedGridsTaggingCriterion(m_dx, m_level, 2.0 * m_p.L, m_p.center),
        current_state, tagging_criterion, disable_simd());
}
