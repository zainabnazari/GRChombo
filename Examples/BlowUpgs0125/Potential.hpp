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
        double phi_min = 0.048388550027;

        V_of_phi = modulus_params.overall_normalization*(1.203799274656e-14+0.000414319862*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0))*exp(-1152.812791250027*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0)))-3.368178423958e-8*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0))*exp(-576.406395625014*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0))));
 
        // The potential gradient at phi
        dVdphi = modulus_params.overall_normalization*(0.000552426483*pow(phi_min+vars.phi,pow(3.0,-1.0))*exp(-1152.812791250027*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0)))-4.490904565277e-8*pow(phi_min+vars.phi,pow(3.0,-1.0))*exp(-576.406395625014*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0)))-0.636844315658*pow(phi_min+vars.phi,5.0*pow(3.0,-1.0))*exp(-1152.812791250027*pow(phi_min+vars.phi,4.0*pow(3.0,-1.0)))+0.000025885861*pow(phi_min+vars.phi,5.0*pow(3.0,-1.0))*exp(-576.406 *pow(phi_min+vars.phi,4.0*pow(3.0,-1.0))));
    }
};

#endif /* POTENTIAL_HPP_ */
