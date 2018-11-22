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
        V_of_phi = modulus_params.overall_normalization * 2.9700000000000003e-9*pow(2.718281828459045,3.4641016151377544*(4.1083165570994 + vars.phi)) + 0.05*pow(2.718281828459045,-0.1*pow(2.718281828459045,1.1547005383792517*(4.1083165570994 + vars.phi)) - 2.3094010767585034*(4.1083165570994 + vars.phi))*(-0.0001 + pow(2.718281828459045,-0.1*pow(2.718281828459045,1.1547005383792517*(4.1083165570994 + vars.phi))) + 0.03333333333333333*pow(2.718281828459045,-0.1*pow(2.718281828459045,1.1547005383792517*(4.1083165570994 + vars.phi)) + 1.1547005383792517*(4.1083165570994 + vars.phi)));

        // The potential gradient at phi
        // 2 V0 a Exp[-a phi] (1-Exp[-a phi])
        dVdphi = 2 * modulus_params.overall_normalization * -1.0288381796959132e-8*pow(2.718281828459045,3.4641016151377544*(4.1083165570994 + vars.phi)) + 0.05*pow(2.718281828459045,-0.1*pow(2.718281828459045,1.1547005383792517*(4.1083165570994 + vars.phi)) - 2.3094010767585034*(4.1083165570994 + vars.phi))*(-2.3094010767585034 - 0.11547005383792518*pow(2.718281828459045,1.1547005383792517*(4.1083165570994 + vars.phi))) * (-0.0001 + pow(2.718281828459045,-0.1*pow(2.718281828459045,1.1547005383792517*(4.1083165570994 + vars.phi))) + 0.03333333333333333*pow(2.718281828459045,-0.1*pow(2.718281828459045,1.1547005383792517*(4.1083165570994 + vars.phi)) + 1.1547005383792517*(4.1083165570994 + vars.phi))) + 0.05* pow(2.718281828459045,-0.1*pow(2.718281828459045,1.1547005383792517*(4.1083165570994 + vars.phi)) - 2.3094010767585034*(4.1083165570994 + vars.phi))* (-0.11547005383792518*pow(2.718281828459045,-0.1*pow(2.718281828459045,1.1547005383792517*(4.1083165570994 + vars.phi)) +1.1547005383792517*(4.1083165570994 + vars.phi)) + 0.03333333333333333*pow(2.718281828459045,-0.1*pow(2.718281828459045,1.1547005383792517*(4.1083165570994 + vars.phi)) +1.1547005383792517*(4.1083165570994 + vars.phi))*(1.1547005383792517 -0.11547005383792518*pow(2.718281828459045,1.1547005383792517*(4.1083165570994 + vars.phi))));
    }
};

#endif /* POTENTIAL_HPP_ */
