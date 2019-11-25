/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(VARDATAINTERPOLATION_HPP_)
#error "This file should only be included through VarDataInterpolation.hpp"
#endif

#ifndef VARDATAINTERPOLATION_IMPL_HPP_
#define VARDATAINTERPOLATION_IMPL_HPP_

// Compute function for interpolation
inline void VarDataInterpolation::compute(Cell<double> current_cell) const
{
    // Define Coordinates
    Coordinates<double> coords(current_cell, m_dx, m_center);

    // Do the read in
    double var;
    var = interpolate_variable(coords);
    var += m_offset_to_data;

    // store vars
    current_cell.store_vars(var, m_var_num);
}

// inteprolate the variable at this coord
inline double
VarDataInterpolation::interpolate_variable(Coordinates<double> coords) const
{
    // currently only works for 3d data
    CH_assert(CH_SPACEDIM == 3);

    // position and value of variable at current cell coords
    std::array<double, CH_SPACEDIM> position = {coords.x, coords.y, coords.z};
    double interpolated_var;

    // gather the surrounding points, but zero them first
    std::array<std::array<double, CH_SPACEDIM + 1>, 8> surrounding_data;
    for (int idx = 0; idx < 8; idx++)
    {
        surrounding_data[idx] = {0, 0, 0, 0};
    }
    std::array<bool, CH_SPACEDIM> exact_data_found = {false, false, false};

    // assumes a uniform spacing of data in x y and z dirs
    // so dx is same for all. This way of calculating it also assumes
    // that the data is written as sequential in at least one direction, and
    // constant in the others. It should be easy to generate data in this format
    // so there seems little point in making it more general
    double spacing = max(abs(m_var_data[0][0] - m_var_data[1][0]),
                         max(abs(m_var_data[0][1] - m_var_data[1][1]),
                             abs(m_var_data[0][2] - m_var_data[1][2])));

    // loop through the array of data
    for (int idata = 0; idata < m_var_data.size(); idata++)
    {
        std::array<double, CH_SPACEDIM + 1> data = m_var_data[idata];
        bool is_surrounding_point = true;
        int idx = 0;
        FOR1(i)
        {
            if (abs(position[i] - data[i]) >= spacing)
            {
                // is not a surrounding point so break
                is_surrounding_point = false;
                break;
            }
            else if ((position[i] - data[i] < spacing) &&
                     (position[i] - data[i] > 0))
            {
                idx += 0 * pow(2, i);
            }
            else if ((data[i] - position[i] < spacing) &&
                     (data[i] - position[i] > 0))
            {
                idx += 1 * pow(2, i);
            }
            else if (data[i] == position[i])
            {
                // assign value on low cells only
                // the use of dx = 0 below will discard the high points
                idx += 0 * pow(2, i);
                exact_data_found[i] = true;
            }
        }
        if (is_surrounding_point)
        {
            surrounding_data[idx] = data;
        }
    }

    // do the interpolation, assumes cuboid grid of points
    std::array<double, CH_SPACEDIM> x0 = {
        surrounding_data[0][0], surrounding_data[0][1], surrounding_data[0][2]};
    std::array<double, CH_SPACEDIM> dx = {0.0, 0.0, 0.0};
    FOR1(i)
    {
        if (!exact_data_found[i]) // otherwise keep as zero
        {
            dx[i] = (position[i] - x0[i]) / spacing;
        }
    }

    // first interpolation in x
    int var_idx = CH_SPACEDIM; // ie, 3, so 4th entry with zero base
    double c00 = surrounding_data[0][var_idx] * (1.0 - dx[0]) +
                 surrounding_data[1][var_idx] * dx[0];
    double c10 = surrounding_data[2][var_idx] * (1.0 - dx[0]) +
                 surrounding_data[3][var_idx] * dx[0];
    double c01 = surrounding_data[4][var_idx] * (1.0 - dx[0]) +
                 surrounding_data[5][var_idx] * dx[0];
    double c11 = surrounding_data[6][var_idx] * (1.0 - dx[0]) +
                 surrounding_data[7][var_idx] * dx[0];

    // then interpolation in y
    double c0 = c00 * (1.0 - dx[1]) + c10 * dx[1];
    double c1 = c01 * (1.0 - dx[1]) + c11 * dx[1];

    // finally interpolation in z
    interpolated_var = c0 * (1.0 - dx[2]) + c1 * dx[2];

    // return the interpolated value
    return interpolated_var;
}

#endif /* VARDATAINTERPOLATION_IMPL_HPP_ */
