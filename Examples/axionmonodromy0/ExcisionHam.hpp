/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISIONHAM_HPP_
#define EXCISIONHAM_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Does excision for fixed BG BH solutions
//! Note that it is does not using simd so one must set disable_simd()
class ExcisionHam
{
  protected:
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The BH center

  public:
    ExcisionHam(const double a_dx,
                         const std::array<double, CH_SPACEDIM> a_center)
        : m_dx(a_dx), m_center(a_center)
    {
    }

    void compute(const Cell<double> current_cell) const
    {
        const Coordinates<double> coords(current_cell, m_dx, m_center);
        double chi = current_cell.load_vars(c_chi);
        double r = coords.get_radius();

        if ((r > 32.0) || (chi < 0.3))
        {
            current_cell.store_vars(0.0, c_absHam);
            current_cell.store_vars(0.0, c_absHamTot);
        }
        // else do nothing
    }
};

#endif /* EXCISIONHAM_HPP_ */
