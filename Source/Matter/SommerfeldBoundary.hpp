/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SOMMERFELDBOUNDARY_HPP_
#define SOMMERFELDBOUNDARY_HPP_

//#include "simd.hpp" for now not vectorised - not worth it for few points
// and unlikely to operate on whole simd object
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4.hpp"
#include "OneSidedDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS
#include "VarsTools.hpp"

//!  Calculates RHS for Sommerfeld boundaries
/*!
     The class calculates the RHS evolution for the variables at the boundaries.
     It does not assume a specific form of matter but is templated over a matter
   class matter_t. Please see the class SFMatter as an example of a matter_t.
     \sa CCZ4(), ScalarField(), VectorField()
*/

template <class matter_t> class SommerfeldBoundary
{
  public:
    // Use MatterCCZ4 as we need all vars
    template <class data_t>
    using Vars = typename MatterCCZ4<matter_t>::template Vars<data_t>;

  protected:
    const double m_L;                     //!< The Length of the grid
    const FourthOrderDerivatives m_deriv; //!< An object for calculating
                                          //!< derivatives of the vars at the
                                          //!< point
    const OneSidedDerivatives
        m_deriv_onesided; //!< An object for calculating one sided derivatives
    const double m_dx;    //!< The grid spacing
    const double m_time;  //!< The current time

  public:
    //!  Constructor of class SommerfeldBoundary
    /*!
         Inputs are the box length and grid spacing
    */
    SommerfeldBoundary(double a_L, double a_dx, double a_time = 0.0)
        : m_dx(a_dx), m_L(a_L), m_deriv(a_dx), m_time(a_time),
          m_deriv_onesided(a_dx)
    {
    }

    //!  The compute member which calculates the RHS at each point in the box
    //!  \sa matter_rhs_equation()
    void compute(Cell<double> current_cell) const;

    //! The function which calculates the new rhs term on boundaries
    double sommerfeld_rhs(const int ivar, const Tensor<1, double> position,
                          const Tensor<1, int> boundary, const double radius,
                          const Cell<double> current_cell) const;
};

#include "SommerfeldBoundary.impl.hpp"

#endif /* SOMMERFELDBOUNDARY_HPP_ */
