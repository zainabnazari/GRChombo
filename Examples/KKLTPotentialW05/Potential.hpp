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
//	double decay_constant;
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
        // The potential value at phi
        // W0 = e-5
        double phi_min = 1.428314650033344;
	double pi = M_PI;
        V_of_phi = modulus_params.overall_normalization*(4.868172392485633e-11*exp(-4.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
                 					   +41.887902047863*exp(-2.0*pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))*exp(-4.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
                        					*(-0.00003*exp(pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))+10.0*(3.0+pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))));

        // The potential gradient at phi
        dVdphi = modulus_params.overall_normalization*(-1.124256256505e-10*exp(-4.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
							 +41.887902047863*exp(-2.0*pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))*exp(-4.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
								*(36.275987284684*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))-0.000108828*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))*exp(pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))))
							 -96.735966092491*exp(-2.0*pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))*exp(-4.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
								*(-0.00003*exp(pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))+10.0*(3.0+pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))))
							 -303.905000414083*exp(-2.0*pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))*exp(-2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
                                                                *(-0.00003*exp(pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))+10.0*(3.0+pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))));
	}
};

#endif /* POTENTIAL_HPP_ */
