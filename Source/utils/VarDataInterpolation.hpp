/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef VARDATAINTERPOLATION_HPP_
#define VARDATAINTERPOLATION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "ScalarField.hpp"
#include "SmallDataIO.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs c_NUM
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which interpolates data read in from external files
// The read in method is a bit inefficient but it only gets done once so let's
// not be too fussy. Note that it requires disable_simd() to be passed to
// BoxLoops as this makes the comparative statements more readable
class VarDataInterpolation
{
  public:
    VarDataInterpolation(
        const double a_dx, const int a_var_num,
        const std::array<double, CH_SPACEDIM> a_center,
        const std::vector<std::array<double, CH_SPACEDIM + 1>> a_var_data,
        const double a_offset_to_data = 0.0)
        : m_dx(a_dx), m_var_num(a_var_num), m_center(a_center),
          m_var_data(a_var_data), m_offset_to_data(a_offset_to_data)
    {
    }

  protected:
    const double m_dx;
    const double m_offset_to_data;
    const int m_var_num;
    const std::array<double, CH_SPACEDIM> m_center;
    const std::vector<std::array<double, CH_SPACEDIM + 1>> m_var_data;

  public:
    //! Function to compute the value of all the initial vars on the grid
    void compute(Cell<double> current_cell) const;

    // linearly interpolate var onto grid from surrounding 8 points
    double interpolate_variable(Coordinates<double> coords) const;
};

#include "VarDataInterpolation.impl.hpp"

#endif /* VARDATAINTERPOLATION_HPP_ */
