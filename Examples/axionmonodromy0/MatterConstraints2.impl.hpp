/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MATTERCONSTRAINTS2_HPP_)
#error "This file should only be included through MatterConstraints2.hpp"
#endif

#ifndef MATTERCONSTRAINTS2_IMPL_HPP_
#define MATTERCONSTRAINTS2_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class matter_t>
MatterConstraints2<matter_t>::MatterConstraints2(const matter_t a_matter,
                                               double dx, double G_Newton)
    : Constraints(dx, 0.0 /*No cosmological constant*/), my_matter(a_matter),
      m_G_Newton(G_Newton)
{
}

template <class matter_t>
template <class data_t>
void MatterConstraints2<matter_t>::compute(Cell<data_t> current_cell) const
{
    // Load local vars and calculate derivs
    const auto vars = current_cell.template load_vars<BSSNMatterVars>();
    const auto d1 = m_deriv.template diff1<BSSNMatterVars>(current_cell);
    const auto d2 = m_deriv.template diff2<BSSNMatterVars>(current_cell);

    // Get the non matter terms for the constraints
    const data_t chi_regularised = simd_max(1e-6, vars.chi);
    auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);
    auto ricci = CCZ4Geometry::compute_ricci(vars, d1, d2, h_UU, chris);
    auto A_UU = TensorAlgebra::raise_all(vars.A, h_UU);
    data_t tr_A2 = TensorAlgebra::compute_trace(vars.A, A_UU);

    // Ham
    data_t Ham = ricci.scalar +
              (GR_SPACEDIM - 1.) * vars.K * vars.K / GR_SPACEDIM - tr_A2;

    // Ham Tot
    data_t HamTot = abs(ricci.scalar) + (GR_SPACEDIM - 1.) * vars.K * vars.K / GR_SPACEDIM 
                  + abs(tr_A2);

    // Energy Momentum Tensor
    const auto emtensor = my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);
    Ham += -16.0 * M_PI * m_G_Newton * emtensor.rho;
    HamTot += abs(16.0 * M_PI * m_G_Newton * emtensor.rho);

    // Write the constraints into the output FArrayBox
    Ham = abs(Ham);
    current_cell.store_vars(Ham, c_absHam);
    current_cell.store_vars(HamTot, c_absHamTot);
}

#endif /* MATTERCONSTRAINTS2_IMPL_HPP_ */
