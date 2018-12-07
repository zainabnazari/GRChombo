/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EMPTYLEVEL_HPP_
#define EMPTYLEVEL_HPP_

#include "GRAMRLevel.hpp"

class EmptyLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<EmptyLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // initialize data
    virtual void initialData() { m_state_new.setVal(0.); }

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time)
    {
    }

    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state){};
};

#endif /* EMPTYLEVEL_HPP_ */
