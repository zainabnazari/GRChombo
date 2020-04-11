/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SCALARPOTENTIAL_HPP_
#define SCALARPOTENTIAL_HPP_

#include "simd.hpp"

class ScalarPotential
{
  public:
    struct params_t
    {
        double scalar_mass;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    ScalarPotential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        V_of_phi = 0.5 * pow(m_params.scalar_mass * vars.phi, 2.0);

        // The potential gradient at phi
        dVdphi = pow(m_params.scalar_mass, 2.0) * vars.phi;
    }
};

#endif /* SCALARPOTENTIAL_HPP_ */
