/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_phi,
    c_chi,
    c_lapse,
    c_Ham,
    c_VofPhi,

    NUM_VARS
};

namespace UserVariables
{
static constexpr char const *variable_names[NUM_VARS] = {
    "phi",    "chi",  "lapse",  "Ham", "VofPhi"};
}

// assign another enum to the other CCZ4 variables
enum
{
    c_blank1,
    c_blank2,
    c_blank3,
    c_blank4,
    c_blank5,

    c_h11,
    c_h12,
    c_h13,
    c_h22,
    c_h23,
    c_h33,

    c_K,

    c_A11,
    c_A12,
    c_A13,
    c_A22,
    c_A23,
    c_A33,

    c_Theta,

    c_Gamma1,
    c_Gamma2,
    c_Gamma3,

    c_shift1,
    c_shift2,
    c_shift3,

    c_B1,
    c_B2,
    c_B3,

    c_Pi,  //(minus) conjugate momentum

    c_Omega,

    c_Mom1,
    c_Mom2,
    c_Mom3,

    NUM_VARS_CCZ4
};

#endif /* USERVARIABLES_HPP */
