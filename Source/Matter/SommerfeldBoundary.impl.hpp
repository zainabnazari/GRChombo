/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SOMMERFELDBOUNDARY_HPP_)
#error "This file should only be included through SommerfeldBoundary.hpp"
#endif

#ifndef SOMMERFELDBOUNDARY_IMPL_HPP_
#define SOMMERFELDBOUNDARY_IMPL_HPP_

template <class matter_t>
void SommerfeldBoundary<matter_t>::compute(Cell<double> current_cell) const
{
    std::array<double, 3> center;
    center.fill(0.5 * m_L);
    const Coordinates<double> coords(current_cell, m_dx, center);
    Tensor<1, double> position = {coords.x, coords.y, coords.z};
    Tensor<1, int> boundary;
    Vars<double> boundary_rhs;

    FOR1(i)
    {
        boundary[i] = 0;
        // check if on positive boundary
        if ((m_L / 2.0 - 3.0 * m_dx) < position[i])
        {
            boundary[i] = 1;
        }

        // check if on negative boundary
        if ((-m_L / 2.0 + 3.0 * m_dx) > position[i])
        {
            boundary[i] = -1;
        }
    }

    // if on a boundary, update rhs for Sommerfeld boundary condition
    if (boundary[0] != 0 || boundary[1] != 0 || boundary[2] != 0)
    {
        // radial distance from centre
        double radius = coords.get_radius();

        // recalculate the boundaries for all vars
        boundary_rhs.enum_mapping([this, &radius, &position, &boundary,
                                   &current_cell](const int &ivar,
                                                  double &var) {
            var =
                sommerfeld_rhs(ivar, position, boundary, radius, current_cell);
        });

        // Add non zero asymptotic vars for the CCZ4 vars
        boundary_rhs.chi += 1.0 / radius;
        boundary_rhs.lapse += 1.0 / radius;
        boundary_rhs.h[0][0] += 1.0 / radius;
        boundary_rhs.h[1][1] += 1.0 / radius;
        boundary_rhs.h[2][2] += 1.0 / radius;

        // Call static functions which add asymptotic values for matter
        matter_t::add_asymptotic_vals(boundary_rhs, radius, m_time);

        // VarsTools::assign(boundary_rhs, 0.);

        // Write the rhs into the output cell
        current_cell.store_vars(boundary_rhs);

    } // if not on boundary do nothing - do not overwrite calculated rhs
}

template <class matter_t>
double SommerfeldBoundary<matter_t>::sommerfeld_rhs(
    const int ivar, const Tensor<1, double> position,
    const Tensor<1, int> boundary, const double radius,
    const Cell<double> current_cell) const
{
    // Calculate the revised rhs, assuming asymptotic values zero
    double variable = current_cell.load_vars(ivar);
    double boundary_rhs = -variable / radius;

    Tensor<1, double> d1_normal;
    Tensor<1, double> d1_positive;
    Tensor<1, double> d1_negative;

    FOR1(idir)
    {
        // now add the derivatives, depending on whether on boundary and if so
        // which for dirs on the positive boundary, use negative derivs
        if (boundary[idir] == 1)
        {
            m_deriv_onesided.diff1(d1_negative, current_cell, idir, ivar,
                                   NEGATIVE_DIRECTION);
            boundary_rhs += -d1_negative[idir] * position[idir] / radius;
        }

        // for dirs on the negative boundary, use positive derivs
        else if (boundary[idir] == -1)
        {
            m_deriv_onesided.diff1(d1_positive, current_cell, idir, ivar,
                                   POSITIVE_DIRECTION);
            boundary_rhs += -d1_positive[idir] * position[idir] / radius;
        }

        // for dirs not on the boundary, use normal derivs
        else
        {
            m_deriv.diff1(d1_normal, current_cell, idir, ivar);
            boundary_rhs += -d1_normal[idir] * position[idir] / radius;
        }
    }

    return boundary_rhs;
}

#endif /* SOMMERFELDBOUNDARY_IMPL_HPP_ */
