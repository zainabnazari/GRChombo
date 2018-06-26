/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_chi,

    c_A11, // put Aij here since it changes during relaxation
    c_A12,
    c_A13,
    c_A22,
    c_A23,
    c_A33,

    c_h11,
    c_h12,
    c_h13,
    c_h22,
    c_h23,
    c_h33,

    c_K,

    c_Theta,

    c_Gamma1,
    c_Gamma2,
    c_Gamma3,

    c_lapse,

    c_shift1,
    c_shift2,
    c_shift3,

    c_B1,
    c_B2,
    c_B3,

    c_phi, // matter field added
    c_Pi,  //(minus) conjugate momentum

    c_rho,

    c_Weyl4_Re,
    c_Weyl4_Im,
    c_Madm,

    c_Ham,

    c_Mom1,
    c_Mom2,
    c_Mom3,

    NUM_VARS
};

namespace UserVariables
{
static constexpr char const *variable_names[NUM_VARS] = {
    "chi",

    "A11",      "A12",      "A13",    "A22", "A23", "A33",

    "h11",      "h12",      "h13",    "h22", "h23", "h33",

    "K",

    "Theta",

    "Gamma1",   "Gamma2",   "Gamma3",

    "lapse",

    "shift1",   "shift2",   "shift3",

    "B1",       "B2",       "B3",

    "phi",      "Pi",

    "rho",

    "Weyl4_Re", "Weyl4_Im", "c_Madm",

    "Ham",      "Mom1",     "Mom2",   "Mom3"};
}

#endif /* USERVARIABLES_HPP */
