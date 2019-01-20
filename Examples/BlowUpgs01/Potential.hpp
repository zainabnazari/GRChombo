/*GRChombo
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
        // The potential value at phi (gs=0.1)
        double phi_min = 0.017916743164;

        V_of_phi = modulus_params.overall_normalization*(1.223107208054e-17+0.000154459879*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0))*exp(-5308.586346781779*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0)))-1.270705407945e-9*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0))*exp(-2654.293173390889*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0))));
 
        // The potential gradient at phi
        dVdphi = modulus_params.overall_normalization*(0.000205946505*pow(phi_min+vars.phi,pow(3.0,-1.0))*exp(-5308.586346781779*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0)))-1.694273877261e-9*pow(phi_min+vars.phi,pow(3.0,-1.0))*exp(-2654.293173390889*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0)))-1.093284807992*pow(phi_min+vars.phi,5.0*pow(3.0,-1.0))*exp(-5308.586346781779*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0)))+4.497099586268e-6*pow(phi_min+vars.phi,5.0*pow(3.0,-1.0))*exp(-2654.293173390889*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0))));
    }
};

#endif /* POTENTIAL_HPP_ */
