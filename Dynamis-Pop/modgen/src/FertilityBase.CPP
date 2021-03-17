#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/FertilityBase.mpp"
//LABEL(FertilityBase, EN) The base version of fertility

/* NOTE(FertilityBase, EN)
This module contains the dummy base version of fertility.
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters 
//{
//    //EN Age distribution of fertility
//    double  AgeSpecificFertility[FERTILE_AGE_RANGE][SIM_YEAR_RANGE];
//
//    //EN Total Fertility Rate (TFR)
//    double  TotalFertilityRate[SIM_YEAR_RANGE];
//
//    //EN Age specific fertility rate
//    model_generated	double	AgeSpecificFertilityRate[FERTILE_AGE_RANGE][SIM_YEAR_RANGE];
//};

//parameter_group PG01_FertilityBase          //EN Fertility base version
//{
//    AgeSpecificFertility,
//    TotalFertilityRate
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor states and functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    
//    event timeFertilityBaseEvent, FertilityBaseEvent;           //EN Mortality event base model
//    void  doIncreaseParity();                                   //EN Increase parity and reset in_this_parity spell
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Person::doIncreaseParity()
{
    // reset spell to end old spell
     poStateFunctions->Set_this_parity_spell(  FALSE );
    
    // In past check if parity is not set already and then increase
    if (! is_mortal && ( creation_type == CT_START ||  creation_type == CT_POOL)
        &&  lStartValues &&  parity < lStartValues ->StartPopValue[PMC_PARITY])
    {
        if ( parity < MAX(PARITY_RANGE))  poStateFunctions->Set_parity(   parity + 1 );
    }
    // this is the regular case
    else if ( parity < MAX(PARITY_RANGE))  poStateFunctions->Set_parity(   parity + 1 );

    // reset spell to initiate new spell
     poStateFunctions->Set_this_parity_spell(  TRUE );
}

TIME Person::timeFertilityBaseEvent()
{
    TIME dEventTime = TIME_INFINITE;
    double dHazard = 0.0;
    if ( is_fertile && gprParam->FertilityModel != FEM_DETAIL)
    {
        dHazard = gprParam->AgeSpecificFertilityRate[CHKBND(RANGE_POS(FERTILE_AGE_RANGE, integer_age ),38,"AgeSpecificFertilityRate",0)][CHKBND(RANGE_POS(SIM_YEAR_RANGE, calendar_year ),151,"AgeSpecificFertilityRate",1)];
        if (dHazard > 0.0) dEventTime = WAIT(-TIME(log(RandUniform(3)) / dHazard));
    }
    return dEventTime;
}

void Person::FertilityBaseEvent()
{
    HandleFertility();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Presimulation: calculate papameter AgeSpecificFertilityRate[nAge][nYear]
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PreSimulation3()
{
    double dSum;
    for (int nYear = 0; nYear < SIZE(SIM_YEAR_RANGE); nYear++)
    {
        dSum = 0.0;
        // check if distribution parameter adds up too 1
        for (int nAge = 0; nAge < SIZE(FERTILE_AGE_RANGE); nAge++)
        {
            dSum = dSum + gprParam->AgeSpecificFertility[CHKBND(nAge,38,"AgeSpecificFertility",0)][CHKBND(nYear,151,"AgeSpecificFertility",1)];
        }
        // scale distribution to 1 and convert to fertility rates; copy to model generated parameter
        for (int nAge = 0; nAge < SIZE(FERTILE_AGE_RANGE); nAge++)
        {
            if (dSum > 0.0)
            {
                gprParam->AgeSpecificFertilityRate[CHKBND(nAge,38,"AgeSpecificFertilityRate",0)][CHKBND(nYear,151,"AgeSpecificFertilityRate",1)]
                    = gprParam->AgeSpecificFertility[CHKBND(nAge,38,"AgeSpecificFertility",0)][CHKBND(nYear,151,"AgeSpecificFertility",1)] / dSum * gprParam->TotalFertilityRate[CHKBND(nYear,151,"TotalFertilityRate",0)];
            }
            else gprParam->AgeSpecificFertilityRate[CHKBND(nAge,38,"AgeSpecificFertilityRate",0)][CHKBND(nYear,151,"AgeSpecificFertilityRate",1)] = 0.0;
        }
    }
}
} // namespace 
