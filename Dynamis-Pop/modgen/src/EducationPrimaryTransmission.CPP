#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/EducationPrimaryTransmission.mpp"
//LABEL (EducationPrimaryTransmission, EN) Primary Education Intergenerational Transmission
/* NOTE(EducationPrimaryTransmission, EN)
 
TODO adapt text

This module allows to add relative factors (odds ratios) to the modeling of primary education outcomes. Currently mother's education by sex is introduced. The relative factors are used to modify the base probabilities of school entry and school retention until graduation. Aggregated outcomes by district and sex are calibrated in order to match the base model's outcome. This calibration is performed either (1) for all years of birth, or (2) just once. In the first case, aggregate outcomes by district and sex are identical, but the model picks different children entering and graduating school, accounting for the relative differences (odds ratios) by mother's education. In the second case, the calibration is done for a selected year of birth. For all following cohorts, parameters are frozen and all trends result from composition effects due to the changing educational composition of mothers. The module only affects persons born in the simulation as mother's education is unknown for others.

**Parameters:**

- Model selection: the 3 choices are (1) Use Base Model, (2) Use the refined model calibrated for all years of birth, (3) Use the refined model calibrated once.
- Odds for starting primary education by mother's education
- Odds for graduating from primary education by mother's education
- The first birth cohort from which onwards the refined model is used, which can be any year beginning from the starting year of the simulation.

The module can be added or removed from the model without requiring modification in other modules. The only exception is the requirement to initialize mother's education.

In its core, this module consists of two functions.

- The function SetEduc1AdjustmentFactors() calculates the required adjustment factors (in the form of log odds by district and sex) to be applied in addition to the relative factors in order to get probabilities further broken down by mother's education. The calculations are performed once at the end of each relevant calendar year as at this point in time the educational composition of mothers of all born in this year is known. Adjustment factors are either calculated once or updated each year according to the user's model selection.
- The function AdjustEducOne() applies the adjustments at the Person level at the first end of the year after birth. Thereby this module modifies the three "fate" states introduced in the base model, namely the two individual progression probabilities educ_one_prob1 and educ_one_prob2, as well as random outcome educ_one_fate based on the new probabilities.

**How to change or expand the list of relative factors**

The list of population groups is declared in the classification EDUC1_GROUP. The classification can be changed or expanded whereby a person can belong only to one of the groups. (For example, when adding an additional dimension ethnicity, all combinations of ethnicity and mother's education have to be listed)
- The individual group membership is a derived state educ1_group. When changing the levels of EDUC1_GROUP, the calculation of the derived state has to be updated accordingly.
- No other code changes are required
*/


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor Sets
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//EN Babies by sex, geo , and mother's education
//actor_set Person asSimBornAge0[sex][geo][educ1_group] filter is_alive && creation_type==CT_BIRTH && integer_age==0;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Types
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//classification EDUC1_GROUP                                      //EN Population Groups
//{
//    E1G_00,                                                     //EN Mothers education low
//    E1G_01,                                                     //EN Mothers education medium
//    E1G_02                                                      //EN Mothers education high
//};

//classification EDUC1_MODEL                                      //EN Primary Education Model Selection
//{
//    E1M_BASE,                                                   //EN Use base model
//    E1M_REFINED_ALIGNALL,                                       //EN Use refined model aligned to base for all birth cohorts
//    E1M_REFINED_ALIGNONCE                                       //EN Use refined model aligned to base once
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters
//{
//    double      Educ1StartOdds[EDUC1_GROUP][SEX];               //EN Odds of starting primary school
//    double      Educ1GradOdds[EDUC1_GROUP][SEX];                //EN Odds of primary school graduation
//    EDUC1_MODEL Educ1Model;                                     //EN Model Selection
//    int         Educ1FirstCohortRefinedModel;                   //EN First birth cohort to apply refined model
//};

//parameter_group PG11_SchoolOneFateRefined                       //EN Primary School Fate - Refined Model
//{
//    Educ1StartOdds, Educ1GradOdds,
//    Educ1FirstCohortRefinedModel
//};

//parameter_group PG12_SchoolOneFate                              //EN Primary School 
//{
//    Educ1Model,
//    PG10_SchoolOneFateBase,
//    PG11_SchoolOneFateRefined,
//    PG10_ShoolOneTracking
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor states and functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Calibrator
//{
//    double alignment_educ1_medium[SEX][GEO];                    //EN Alignment low to medium
//    double alignment_educ1_high[SEX][GEO];                      //EN Alignment medium to high
//    void   SetEduc1AdjustmentFactors();                         //EN Set calibration factors at CalibratorYearEnd()
//    
//    //EN Probability adjustment
//    double AdjustedProbability(double dProb, double dLogOddEduc, double dLogOddAdjust);     
//};

//actor Person
//{
//    EDUC_ONE_LEVEL educ_mother = { EOL_LOW };                   //EN Mother's education set in Start() if born in simulation
//
//    //EN Education risk group
//    EDUC1_GROUP educ1_group = ( educ_mother == EOL_LOW ) ? E1G_00 : ( educ_mother == EOL_MEDIUM ) ?  E1G_01 : E1G_02;
//
//    void AdjustEducOne();                                       //EN Adjust education fate at year end
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Calibrator::SetEduc1AdjustmentFactors()
{
    if (  calibrator_year >= MIN(SIM_YEAR_RANGE)  && (
       (gprParam->Educ1Model == E1M_REFINED_ALIGNALL  &&  calibrator_year >= gprParam->Educ1FirstCohortRefinedModel ) ||
       (gprParam->Educ1Model == E1M_REFINED_ALIGNONCE &&  calibrator_year == gprParam->Educ1FirstCohortRefinedModel )))
    {
        for (int nSex = 0; nSex < SIZE(SEX); nSex++ )
        {
            for (int nGeo = 0; nGeo < SIZE(GEO); nGeo++ )
            {
                // calculate total population for given geo and sex
                double nTotalPop = 0.0;
                for ( int dGroup = 0; dGroup < SIZE(EDUC1_GROUP); dGroup++ )
                {
                    nTotalPop = nTotalPop + asSimBornAge0[nSex][nGeo][dGroup]->Count();
                }

                // Run Adjustment for school entry
                double dFactorEntry = 0.0;
                if ( nTotalPop > 0.0 )
                {
                    int nIterations = 0;
                    double dResultProb = 10.0;
                    double dTargetProb = gprParam->EducTrans1[CHKBND(nSex,2,"EducTrans1",0)][CHKBND(RANGE_POS(YOB_EDUC_TRANS1, calibrator_year ),63,"EducTrans1",1)][CHKBND(nGeo,13,"EducTrans1",2)];
                    double dLower = -10.0;
                    double dUpper = 10.0;
                    double dCenter = 0.0;
                    while (abs(dResultProb - dTargetProb) > 0.001 && nIterations < 1000)
                    {
                        nIterations++;
                        dCenter = (dLower + dUpper) / 2.0;
                        dResultProb = 0.0;
                        for ( int nGroup = 0; nGroup < SIZE(EDUC1_GROUP); nGroup++ )
                        {
                            dResultProb = dResultProb + ( asSimBornAge0[nSex][nGeo][nGroup]->Count() / nTotalPop ) *
                                AdjustedProbability(dTargetProb, log(gprParam->Educ1StartOdds[CHKBND(nGroup,3,"Educ1StartOdds",0)][CHKBND(nSex,2,"Educ1StartOdds",1)]), dCenter);
                        }
                        if (dTargetProb > dResultProb) dLower = dCenter;
                        else dUpper = dCenter;
                    }
                    dFactorEntry = dCenter;
                }
                // set factor
                alignment_educ1_medium[CHKBND(nSex,2,"Calibrator.alignment_educ1_medium",0)][CHKBND(nGeo,13,"Calibrator.alignment_educ1_medium",1)] = dFactorEntry;

                // Run Adjustment for school retention
                double dFactorGrad = 0.0;
                if ( nTotalPop > 0.0 )
                {
                    int nIterations = 0;
                    double dResultProb = 10.0;
                    double dTargetProb = gprParam->EducTrans2[CHKBND(nSex,2,"EducTrans2",0)][CHKBND(RANGE_POS(YOB_EDUC_TRANS2, calibrator_year ),67,"EducTrans2",1)][CHKBND(nGeo,13,"EducTrans2",2)];
                    double dLower = -10.0;
                    double dUpper = 10.0;
                    double dCenter = 0.0;
                    while (abs(dResultProb - dTargetProb) > 0.001 && nIterations < 1000)
                    {
                        nIterations++;
                        dCenter = (dLower + dUpper) / 2.0;
                        dResultProb = 0.0;
                        for ( int nGroup = 0; nGroup < SIZE(EDUC1_GROUP); nGroup++ )
                        {
                            dResultProb = dResultProb + ( asSimBornAge0[nSex][nGeo][nGroup]->Count() / nTotalPop ) *
                                AdjustedProbability(dTargetProb, log(gprParam->Educ1GradOdds[CHKBND(nGroup,3,"Educ1GradOdds",0)][CHKBND(nSex,2,"Educ1GradOdds",1)]), dCenter);
                        }
                        if (dTargetProb > dResultProb) dLower = dCenter;
                        else dUpper = dCenter;
                    }
                    dFactorGrad = dCenter;
                }
                // set factor
                alignment_educ1_high[CHKBND(nSex,2,"Calibrator.alignment_educ1_high",0)][CHKBND(nGeo,13,"Calibrator.alignment_educ1_high",1)] = dFactorGrad;
            }
        }
    }
}

double Calibrator::AdjustedProbability(double dProb, double dLogOddEduc, double dLogOddAdjust)
{
    if (dProb <= 0.0) return 0.0;
    else if (dProb >= 1.0) return 1.0;
    else
    {
        if (dProb >= 0.9999) dProb = 0.9999;
        double dValue = log(dProb / (1 - dProb)) + dLogOddEduc + dLogOddAdjust;
        if (dValue > 50) dValue = 50;
        double dExp = exp(dValue);
        return dExp / (1 + dExp);
    }
}

void Person::AdjustEducOne()
{
    if ( creation_type == CT_BIRTH &&  integer_age == 0 &&  calendar_year >= MIN(SIM_YEAR_RANGE) &&
         calendar_year >= gprParam->Educ1FirstCohortRefinedModel)
    {
        // Adjustment for selection of  refined model aligned to base for all birth cohorts
        if ( gprParam->Educ1Model == E1M_REFINED_ALIGNALL )
        {
            double dProb1 = lCalibrator ->UpdateTime()->AdjustedProbability(
                gprParam->EducTrans1[CHKBND(sex ,2,"EducTrans1",0)][CHKBND(RANGE_POS(YOB_EDUC_TRANS1,year_of_birth ),63,"EducTrans1",1)][CHKBND(geo_birth ,13,"EducTrans1",2)],
                log(gprParam->Educ1StartOdds[CHKBND(educ1_group ,3,"Educ1StartOdds",0)][CHKBND(sex ,2,"Educ1StartOdds",1)]),
                lCalibrator ->alignment_educ1_medium[sex ][geo_birth ]);
            double dProb2 = lCalibrator ->UpdateTime()->AdjustedProbability(
                gprParam->EducTrans2[CHKBND(sex ,2,"EducTrans2",0)][CHKBND(RANGE_POS(YOB_EDUC_TRANS2,year_of_birth ),67,"EducTrans2",1)][CHKBND(geo_birth ,13,"EducTrans2",2)],
                log(gprParam->Educ1GradOdds[CHKBND(educ1_group ,3,"Educ1GradOdds",0)][CHKBND(sex ,2,"Educ1GradOdds",1)]),
                lCalibrator ->alignment_educ1_high[sex ][geo_birth ]);
            EDUC_ONE_LEVEL eolFate = EOL_LOW;
            if (RandUniform(50) < dProb1) eolFate = EOL_MEDIUM;
            if (eolFate == EOL_MEDIUM && RandUniform(51) < dProb2)  eolFate = EOL_HIGH;
             poStateFunctions->Set_educ_one_fate(  eolFate );
        }
        // Adjustment for selection of  refined model aligned to base just once
        else if ( gprParam->Educ1Model == E1M_REFINED_ALIGNONCE )
        {
            double dProb1 = lCalibrator ->UpdateTime()->AdjustedProbability(
                gprParam->EducTrans1[CHKBND(sex ,2,"EducTrans1",0)][CHKBND(RANGE_POS(YOB_EDUC_TRANS1, gprParam->Educ1FirstCohortRefinedModel),63,"EducTrans1",1)][CHKBND(geo_birth ,13,"EducTrans1",2)],
                log(gprParam->Educ1StartOdds[CHKBND(educ1_group ,3,"Educ1StartOdds",0)][CHKBND(sex ,2,"Educ1StartOdds",1)]),
                lCalibrator ->alignment_educ1_medium[sex ][geo_birth ]);
            double dProb2 = lCalibrator ->UpdateTime()->AdjustedProbability(
                gprParam->EducTrans2[CHKBND(sex ,2,"EducTrans2",0)][CHKBND(RANGE_POS(YOB_EDUC_TRANS2, gprParam->Educ1FirstCohortRefinedModel),67,"EducTrans2",1)][CHKBND(geo_birth ,13,"EducTrans2",2)],
                log(gprParam->Educ1GradOdds[CHKBND(educ1_group ,3,"Educ1GradOdds",0)][CHKBND(sex ,2,"Educ1GradOdds",1)]),
                lCalibrator ->alignment_educ1_high[sex ][geo_birth ]);
            EDUC_ONE_LEVEL eolFate = EOL_LOW;
            if (RandUniform(52) < dProb1) eolFate = EOL_MEDIUM;
            if (eolFate == EOL_MEDIUM && RandUniform(53) < dProb2)  eolFate = EOL_HIGH;
             poStateFunctions->Set_educ_one_fate(  eolFate );
        }
    }
}

} // namespace 
