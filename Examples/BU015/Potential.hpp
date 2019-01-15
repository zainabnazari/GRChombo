/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

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
    params_t modulus_params;

  public:
    //! The constructor
    Potential(params_t a_params) : modulus_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at phi (gs=0.15)
        double phi_min = 0.0920595855;

        V_of_phi = modulus_params.overall_normalization * (1.205759642115e-12 + 0.000827720012 * pow(phi_min + vars.phi, 4.0 * pow(3.0,-1.0)) * exp(-415.935390761567 * pow(phi_min + vars.phi, 4.0 * pow(3.0,-1.0))) - 3.104860978085e-7 * pow(phi_min + vars.phi, 4.0 * pow(3.0,-1.0)) * exp(-207.967695380784 * pow(phi_min + vars.phi, 4.0 * pow(3.0,-1.0))));
 
        // The potential gradient at phi
        dVdphi = modulus_params.overall_normalization * (0.001103626683 * pow(phi_min + vars.phi, pow(3.0,-1.0)) * exp(-415.935390761567 * pow(phi_min + vars.phi, 4.0 * pow(3.0,-1.0))) - 4.139814637446e-7 * pow(phi_min + vars.phi, pow(3.0,-1.0)) * exp(-207.967695380784 * pow(phi_min + vars.phi, 4.0 * pow(3.0,-1.0))) - 0.459037395594 * pow(phi_min + vars.phi, 5.0 * pow(3.0, -1.0)) * exp(-415.935390761567 * pow(phi_min + vars.phi, 4.0 * pow(3.0,-1.0))) + 0.000086094771 * pow(phi_min + vars.phi, 5.0 * pow(3.0, -1.0)) * exp(-207.967695380784 * pow(phi_min + vars.phi, 4.0 * pow(3.0,-1.0))) );
    }
};

#endif /* POTENTIAL_HPP_ */
