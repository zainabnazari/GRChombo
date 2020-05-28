/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <iostream>
#include "parstream.H"

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"

class Potential
{
  public:
    struct params_t
    {
        double overall_normalization;
	double decay_constant;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        // m^2 M^2[1-1/\sqrt(1+phi^2/M^2)]
        double pi = M_PI;
//        V_of_phi = pow(m_params.overall_normalization,2.0)(1-1/sqrt(1+pow(vars.phi,2/pow(m_params.overal_normalisation)))

        // The potential gradient at phi
        // m^2 *phi/(1+phi^2/M^2)^(3/2)
//        dVdphi = vars.phi/pow(1+pow(vars.phi,2)/(pow(m_params.overal_normalisation,2)),3/2);
    }
};

#endif /* POTENTIAL_HPP_ */
