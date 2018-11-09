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
        double axion_mass;
	double decay_constant;
    };

  private:
    params_t axion_params;

  public:
    //! The constructor
    Potential(params_t a_params) : axion_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        // m^2 f^2(1 - cos(phi/f)
        V_of_phi = pow(axion_params.axion_mass * axion_params.decay_constant, 2.0) * (1.0 - cos(vars.phi * pow(axion_params.decay_constant, -1.0)));

        // The potential gradient at phi
        // m^2 f sin(phi/f)
        dVdphi = pow(axion_params.axion_mass, 2.0) * axion_params.decay_constant * sin(vars.phi * pow(axion_params.decay_constant, -1.0));
    }
};

#endif /* POTENTIAL_HPP_ */
