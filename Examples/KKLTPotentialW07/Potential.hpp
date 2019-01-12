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
        // V = 
        double phi_min = 1.428314650033344;
        V_of_phi = modulus_params.overall_normalization * (4.868172392485633 * 
                                  pow(10.0,-11.0) * exp(-4.0*sqrt(3.0)*(phi_min + vars.phi)) 
                 + 0.05 * exp( - 0.1*exp((2.0*(phi_min + vars.phi))/sqrt(3.0)) 
                            - (4.0*(phi_min + vars.phi))/sqrt(3.0))
                        * (- 0.0001 + exp(-0.1*exp((2.0*(phi_min + vars.phi))/sqrt(3.0))) 
                           + 0.03333333333333333*exp(- 0.1*exp((2.0*(phi_min + vars.phi))/sqrt(3.0)) 
                                                     + (2.0*(phi_min + vars.phi))/sqrt(3.0))));

        // The potential gradient at phi
        // 2 V0 a Exp[-a phi] (1-Exp[-a phi])
        dVdphi = modulus_params.overall_normalization * 
               (
                 - 1.030027101927e-8 * exp(-2.0*sqrt(3.0)*(phi_min + vars.phi)) 
                 + 0.05 * exp( - 0.1*exp((2.0*(phi_min + vars.phi))/sqrt(3.0)) 
                             - (4.0*(phi_min + vars.phi))/sqrt(3.0)) 
                      * (- 4.0/sqrt(3.0) - 0.115470053838*exp((2.0*(phi_min + vars.phi))/sqrt(3.0))) 
                      * (- 0.0001 + exp(-0.1*exp((2.0*(phi_min + vars.phi))/sqrt(3.0))) 
                         + 0.03333333333333*exp(-0.1*exp((2.0*(phi_min + vars.phi))/sqrt(3.0)) 
                                                + (2.0*(phi_min + vars.phi))/sqrt(3.0))) 
                 + 0.05 * exp( - 0.1*exp((2.0*(phi_min + vars.phi))/sqrt(3.0))
                             - (4.0*(phi_min + vars.phi))/sqrt(3.0))
                      * (- 0.115470053838*exp(- 0.1*exp((2.0*(phi_min + vars.phi))/sqrt(3.0)) 
                                              + (2.0*(phi_min + vars.phi))/sqrt(3.0)) 
                         + 0.0333333333333*exp(- 0.1*exp((2.0*(phi_min + vars.phi))/sqrt(3.0))
                                               + (2.0*(phi_min + vars.phi))/sqrt(3.0))
                                          * (2.0/sqrt(3.0) - 0.115470053838*exp((2.0*(phi_min + vars.phi))/sqrt(3.0))))
               ); 
        //pout() << "phi, V, dVdphi " << vars.phi << " " << V_of_phi << " " << dVdphi << endl;
    }
};

#endif /* POTENTIAL_HPP_ */
