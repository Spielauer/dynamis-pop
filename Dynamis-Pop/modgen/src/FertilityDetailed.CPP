#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/FertilityDetailed.mpp"
//LABEL(FertilityDetailed, EN) The detailed version of fertility

/* NOTE(FertilityDetailed, EN)

TODO update text

This module contains the dummy detailed version of fertility.

This mortality model is used in two modes:

- Unaligned: if (FertilityModel == FEM_DETAIL) this model is used as is. Fertility is triggered by the FertilityDetailedEvent. 
- Aligned: if (FertilityModel == FEM_ALIGNED) this model is used indirectely. Whenever a fertility event is produced from the base mode, the person with the shortest waiting time to giving birth according to the detailed model is identified. In order to search (and return) a link to the person with the shortest waiting time, a Calibrator function GetNextToGiveBirth(Person_ptr prPerson) is used. The parameter is the link to the calling person (the one triggering the event) in order to make personal characteristics used as alignment dimensions accessible. This calibrator function loops through all relevant persons to identify the one with the shortest waiting time. Individual waiting times are returned by the Person function getTimeToBirth().

TODO: 

- Declare appropriate alignment dimensions of asAllFertilePersonsForFertilityAlignment. The default integer age and residence status. Adapt the dimensions wherever the actorset is used.  
- Add the code to calcualte a random waiting time to next birth in getTimeToBirth(). This will probably involve adding or modifying parameters.
- The same algorithm should be used/added to calculate the random time to birth in the time function timeFertilityDetailedEvent()
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor Sets
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//EN All fertile women by age and residence status
//actor_set Person asAllFertilePersonsForFertilityAlignmentByAge[fertile_age][is_resident]
//    filter FertilityModel == FEM_ALIGNED_AGE && is_fertile;

//EN All fertile women by residence status
//actor_set Person asAllFertilePersonsForFertilityAlignment[is_resident]
//    filter FertilityModel == FEM_ALIGNED_TOTAL && is_fertile;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Types
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//classification BIRTH1_GROUP                         //EN Population Groups
//{
//    B1G_00,                                         //EN Never entered school
//    B1G_01,                                         //EN Primary dropout
//    B1G_02                                          //EN Primary graduate
//};

//range PARITY_RANGE  { 0, 15 };                      //EN Parity range 0+
//range PARITY_RANGE1 { 1, 15 };                      //EN Parity range 1+
//range PARITY_RANGE2 { 2, 15 };                      //EN Parity range 2+

//classification EVER_UNION                           //EN Ever in a union
//{
//    EU_NO,                                          //EN No     
//    EU_YES                                          //EN Yes
//};

//classification HIGHER_BIRTHS_PARA                   //EN Parameters for higher order births
//{
//    HBP_PERIOD1,                                    //EN 0-1 years after previous birth
//    HBP_PERIOD2,                                    //EN 1-3 years after previous birth
//    HBP_PERIOD3,                                    //EN 3-6 years after previous birth
//    HBP_PERIOD4,                                    //EN 6-9 years after previous birth
//    HBP_PERIOD5,                                    //EN 9-12 years after previous birth
//    HBP_PERIOD6,                                    //EN 12+ years after previous birth
//    HBP_AGE35,                                      //EN Age 35-39
//    HBP_AGE40,                                      //EN Age 40-44
//    HBP_AGE45,                                      //EN Age 45+
//    HBP_EDUC1,                                      //EN Entered primary school but dropped out
//    HBP_EDUC2                                       //EN Graduated from primary school
//};

//partition BIRTH_AGE_PART { 35, 40, 45 };            //EN Age Groups
//partition DUR_TIME_IN_PARITY { 1, 3, 6, 9, 12 };    //EN Duration episodes since last birth


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters 
//{
//    //EN First Births
//    double  FirstBirthRates[BIRTH1_GROUP][EVER_UNION][FERTILE_AGE_RANGE][BIRTH1_LOC];
//
//    //EN Birth Trends
//    double  BirthTrends[PARITY_RANGE1][SIM_YEAR_RANGE];
//
//    //EN Higher order births
//    double  HigherOrderBirthsPara[HIGHER_BIRTHS_PARA][PARITY_RANGE2];
//};

//parameter_group PG01_FertilityDetailed                          //EN Fertility detailed version
//{
//    FirstBirthRates, HigherOrderBirthsPara, BirthTrends
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor states and functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    PARITY_RANGE parity = { 0 };                                //EN Parity
//    
//    int time_in_parity = self_scheduling_split(active_spell_duration(this_parity_spell,TRUE), DUR_TIME_IN_PARITY);
//    
//    event  timeFertilityDetailedEvent, FertilityDetailedEvent;  //EN Fertility event detailed model
//    double getTimeToBirth();                                    //EN Return time to next birth
//};

//actor Calibrator
//{
//    Person_ptr GetNextToGiveBirth(Person_ptr prPerson);         //EN Identify next person to give birth
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Person::getTimeToBirth()
{
    TIME   tEventTime = TIME_INFINITE;
    double dHazard = 0.0;

    // Individual characteristics
    BIRTH1_GROUP cGroup  = ( educ_one_fate == EOL_LOW ) ? B1G_00 : ( educ_one_fate == EOL_MEDIUM ) ? B1G_01 : B1G_02;
    BIRTH1_LOC   cLocation = GEO_To_BIRTH1_LOC(geo );
    int          nBirthAgePart = SPLIT(integer_age , BIRTH_AGE_PART);    

    if (is_fertile )
    {
        // first birth
        if ( parity == 0)
        {
            dHazard = gprParam->FirstBirthRates[CHKBND(cGroup,3,"FirstBirthRates",0)][CHKBND(ever_union ,2,"FirstBirthRates",1)]
                [CHKBND(RANGE_POS(FERTILE_AGE_RANGE, integer_age ),38,"FirstBirthRates",2)][CHKBND(cLocation,5,"FirstBirthRates",3)]
                * gprParam->BirthTrends[CHKBND(RANGE_POS(PARITY_RANGE1, parity +1),15,"BirthTrends",0)][CHKBND(RANGE_POS(SIM_YEAR_RANGE, calendar_year ),151,"BirthTrends",1)];
        }
        // higher order births
        else
        {
            // Baseline hazard and trend
            dHazard = gprParam->HigherOrderBirthsPara[CHKBND(time_in_parity ,11,"HigherOrderBirthsPara",0)][CHKBND(RANGE_POS(PARITY_RANGE2, parity +1),14,"HigherOrderBirthsPara",1)]
                * gprParam->BirthTrends[CHKBND(RANGE_POS(PARITY_RANGE1, parity +1),15,"BirthTrends",0)][CHKBND(RANGE_POS(SIM_YEAR_RANGE, calendar_year ),151,"BirthTrends",1)];

            // relative risk for older age
            if (nBirthAgePart > 0) //not the baseline age group
            {
                dHazard = dHazard *
                    gprParam->HigherOrderBirthsPara[CHKBND(HBP_AGE35 - 1 + nBirthAgePart,11,"HigherOrderBirthsPara",0)][CHKBND(RANGE_POS(PARITY_RANGE2, parity +1),14,"HigherOrderBirthsPara",1)];
            }
            // relative risk for primary dropouts and graduates
            if ( educ_one_fate == EOL_MEDIUM)
            {
                dHazard = dHazard * gprParam->HigherOrderBirthsPara[CHKBND(HBP_EDUC1,11,"HigherOrderBirthsPara",0)][CHKBND(RANGE_POS(PARITY_RANGE2, parity +1),14,"HigherOrderBirthsPara",1)];
            }
            else if ( educ_one_fate == EOL_HIGH)
            {
                dHazard = dHazard * gprParam->HigherOrderBirthsPara[CHKBND(HBP_EDUC2,11,"HigherOrderBirthsPara",0)][CHKBND(RANGE_POS(PARITY_RANGE2, parity +1),14,"HigherOrderBirthsPara",1)];
            }
        }

        // schedule event
        if (dHazard > 0.0) tEventTime = WAIT(-TIME(log(RandUniform(12)) / dHazard));
    }
    return tEventTime;
}

Person_ptr Calibrator::GetNextToGiveBirth(Person_ptr prPers)
{
    long        nPopSize;
    double      dWaitToBirth = TIME_INFINITE;
    Person_ptr  ptrNextGivingBirth = NULL;
    Person_ptr  ptrThisOne = NULL;

    // Population size by alignment variant
    if (gprParam->FertilityModel == FEM_ALIGNED_AGE) nPopSize = asAllFertilePersonsForFertilityAlignmentByAge[RANGE_POS(FERTILE_AGE_RANGE, prPers->integer_age )][prPers->is_resident ]->Count();
    else nPopSize = asAllFertilePersonsForFertilityAlignment[prPers->is_resident ]->Count();

    for (long nI = 0; nI < nPopSize; nI++)
    {
        // Pointer to person from actor set at position nI by alignment variant
        if (gprParam->FertilityModel == FEM_ALIGNED_AGE) ptrThisOne = asAllFertilePersonsForFertilityAlignmentByAge[RANGE_POS(FERTILE_AGE_RANGE, prPers->integer_age )][prPers->is_resident ]->Item(nI);
        else ptrThisOne = asAllFertilePersonsForFertilityAlignment[prPers->is_resident ]->Item(nI);
        
        // Check if waiting time is shorter
        double dThisWaitTime = ptrThisOne->UpdateTime()->getTimeToBirth();
        if (dThisWaitTime <= dWaitToBirth)
        {
            dWaitToBirth = dThisWaitTime;
            ptrNextGivingBirth = ptrThisOne;
        }
    }
    return ptrNextGivingBirth;
}


TIME Person::timeFertilityDetailedEvent()
{
    TIME   tEventTime = TIME_INFINITE;
    if ( is_fertile && gprParam->FertilityModel == FEM_DETAIL)
    {
        double dHazard = 0.0;

        // Individual characteristics
        BIRTH1_GROUP cGroup = ( educ_one_fate == EOL_LOW) ? B1G_00 : ( educ_one_fate == EOL_MEDIUM) ? B1G_01 : B1G_02;
        BIRTH1_LOC   cLocation = GEO_To_BIRTH1_LOC(geo );
        int          nBirthAgePart = SPLIT(integer_age , BIRTH_AGE_PART);    //EN Age group
    
        // first birth
        if ( parity == 0)
        {
            dHazard = gprParam->FirstBirthRates[CHKBND(cGroup,3,"FirstBirthRates",0)][CHKBND(ever_union ,2,"FirstBirthRates",1)]
                [CHKBND(RANGE_POS(FERTILE_AGE_RANGE, integer_age ),38,"FirstBirthRates",2)][CHKBND(cLocation,5,"FirstBirthRates",3)]
                * gprParam->BirthTrends[CHKBND(RANGE_POS(PARITY_RANGE1, parity +1),15,"BirthTrends",0)][CHKBND(RANGE_POS(SIM_YEAR_RANGE, calendar_year ),151,"BirthTrends",1)];
        }
        // higher order births
        else
        {
            // Baseline hazard and trend
            dHazard = gprParam->HigherOrderBirthsPara[CHKBND(time_in_parity ,11,"HigherOrderBirthsPara",0)][CHKBND(RANGE_POS(PARITY_RANGE2, parity +1),14,"HigherOrderBirthsPara",1)]
                * gprParam->BirthTrends[CHKBND(RANGE_POS(PARITY_RANGE1, parity +1),15,"BirthTrends",0)][CHKBND(RANGE_POS(SIM_YEAR_RANGE, calendar_year ),151,"BirthTrends",1)];

            // relative risk for older age
            if (nBirthAgePart > 0) //not the baseline age group
            {
                dHazard = dHazard *
                    gprParam->HigherOrderBirthsPara[CHKBND(HBP_AGE35 - 1 + nBirthAgePart,11,"HigherOrderBirthsPara",0)][CHKBND(RANGE_POS(PARITY_RANGE2, parity +1),14,"HigherOrderBirthsPara",1)];
            }
            // relative risk for primary dropouts and graduates
            if ( educ_one_fate == EOL_MEDIUM)
            {
                dHazard = dHazard * gprParam->HigherOrderBirthsPara[CHKBND(HBP_EDUC1,11,"HigherOrderBirthsPara",0)][CHKBND(RANGE_POS(PARITY_RANGE2, parity +1),14,"HigherOrderBirthsPara",1)];
            }
            else if ( educ_one_fate == EOL_HIGH)
            {
                dHazard = dHazard * gprParam->HigherOrderBirthsPara[CHKBND(HBP_EDUC2,11,"HigherOrderBirthsPara",0)][CHKBND(RANGE_POS(PARITY_RANGE2, parity +1),14,"HigherOrderBirthsPara",1)];
            }
        }

        // schedule event
        if (dHazard > 0.0) tEventTime = WAIT(-TIME(log(RandUniform(13)) / dHazard));
    }
    return tEventTime;
}

void Person::FertilityDetailedEvent()
{
    HandleFertility();
}
} // namespace 
