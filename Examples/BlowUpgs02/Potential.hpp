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
        // The potential value at phi (gs=0.2)
        double phi_min = 0.1992903401934;

        V_of_phi = modulus_params.overall_normalization*(3.925604834774e-10+0.002090248430*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0))*exp(-115.950986275041*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0)))-5.327012933954e-6*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0))*exp(-57.975493137520*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0))));
 
        // The potential gradient at phi
        dVdphi = modulus_params.overall_normalization*(0.002786997906*pow(phi_min+vars.phi,pow(3.0,-1.0))*exp(-115.950986275041*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0)))-7.102683911938e-6*pow(phi_min+vars.phi,pow(3.0,-1.0))*exp(-57.975493137520*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0)))-0.323155156045*pow(phi_min+vars.phi,5.0*pow(3.0,-1.0))*exp(-115.950986275041*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0)))+0.000411781602*pow(phi_min+vars.phi,5.0*pow(3.0,-1.0))*exp(-57.975493137520*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0))));
    }
};

#endif /* POTENTIAL_HPP_ */
