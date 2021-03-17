#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/MortalityBase.mpp"
//LABEL(MortalityBase, EN) The base version of mortality

/* NOTE(MortalityBase, EN)
This module contains the base version of mortality. It implements mortality based on a life table approach. People can die at any moment based on age-specific mortality rates calculated from a standard life table and a trend factor. The trend factor is found in pre-simulation and calibrates mortality to reach a target period life expectancy. The standard life table and target life expectancies are model parameters.

The MortalityBaseEvent() event calls a function HandleMortality() which is implemented in the MortalityGeneral module. HandleMortality() decides if the base mortality model is used directly - leading to the person's death - or indirectly. In the latter case, the event leads to another person's death, the person to die being selected by a more detailed model.
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters 
//{
//    double  MortalityTable[AGE_RANGE][SEX];	                        //EN Life Table (age profile of hazards)
//    double	LifeExpectancy[SIM_YEAR_RANGE][SEX];		            //EN Life Expectancy
//    model_generated double	MortalityTrend[SIM_YEAR_RANGE][SEX];	//EN Mortality trend
//};

//parameter_group PG02_MortalityBase                                   //EN Mortality base version
//{
//    MortalityTable,
//    LifeExpectancy
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor states and functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    //EN Mortality event base model
//    event timeMortalityBaseEvent, MortalityBaseEvent;    
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TIME Person::timeMortalityBaseEvent()
{
    TIME    dEventTime = TIME_INFINITE;
    double  dMortalityHazard = gprParam->MortalityTable[CHKBND(integer_age ,101,"MortalityTable",0)][CHKBND(sex ,2,"MortalityTable",1)]
          * gprParam->MortalityTrend[CHKBND(RANGE_POS(SIM_YEAR_RANGE, calendar_year ),151,"MortalityTrend",0)][CHKBND(sex ,2,"MortalityTrend",1)];    
    
    // check if a person is at risk
    if ( is_mortal && dMortalityHazard > 0.0 && 
        ( calendar_year < MIN(CHILD_MORTALITY_YEARS) || gprParam->MortalityModel == MOM_BASE ||  integer_age > MAX(CHILD_MORTALITY_AGE)))
    {
        // determine the event time
        dEventTime = WAIT(-log(RandUniform(2)) / dMortalityHazard);
    }
    return dEventTime;
}

void Person::MortalityBaseEvent() {  Finish(); }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pre-Simulation: Finding calibration factors for target life expectancy
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PreSimulation4()
{
    double	dLower, dUpper, dCenter, dResult, dTarget, dAlive, dDeaths;
    int		nIterations;
    for (int nSex = 0; nSex < SIZE(SEX); nSex++)
    {
        for (int nYear = 0; nYear < SIZE(SIM_YEAR_RANGE); nYear++)
        {
            dTarget = gprParam->LifeExpectancy[CHKBND(nYear,151,"LifeExpectancy",0)][CHKBND(nSex,2,"LifeExpectancy",1)];      // Target: life expectancy
            dResult = 0.0;                              // Search result: life expectancy
            nIterations = 10000;                        // Maximum number of iterations in search
            dLower = 0.1;                               // Lower limit of relative risk)
            dUpper = 3.0;                               // Upper limit of relative risk

            while (abs(dResult - dTarget) > 0.0001 && nIterations > 0)
            {
                nIterations--;
                dCenter = (dLower + dUpper) / 2.0;      // New calibration factor for probing

                dResult = 0.0;
                dAlive = 1.0;                           // Proportion of people still alive

                // Life expectancy calculated applying calibration factor
                for (int nAge = 0; nAge < SIZE(AGE_RANGE); nAge++)
                {
                    // proportion of deaths in year: survival = exp(-hazard)
                    dDeaths = dAlive * (1 - exp(-gprParam->MortalityTable[CHKBND(nAge,101,"MortalityTable",0)][CHKBND(nSex,2,"MortalityTable",1)] * dCenter));
                    dAlive = dAlive - dDeaths;
                    // people dying in this year are assumed to die in th emiddle of the year
                    dResult = dResult + dDeaths * 0.5 + dAlive;
                }
                // Moving the search limits for next iteration
                if (dTarget < dResult) dLower = dCenter;
                else dUpper = dCenter;
            }
            // copying the best solution into the model-generated mortality trend parameter
            gprParam->MortalityTrend[CHKBND(nYear,151,"MortalityTrend",0)][CHKBND(nSex,2,"MortalityTrend",1)] = dCenter;
        }
    }
}

} // namespace 
