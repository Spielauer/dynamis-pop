#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/MortalityDetailed.mpp"
//LABEL(MortalityDetailed, EN) The detailed version of mortality

/* NOTE(MortalityDetailed, EN)


*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Type
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//classification CHILD_MORTALITY_GROUP            //EN Relative risks for child mortality
//{
//    CMG_00,                                     //EN Mother age 17+ completed primary
//    CMG_01,                                     //EN Mother age 17+, primary dropout
//    CMG_02,                                     //EN Mother age 17+, neve entered primary
//    CMG_03,                                     //EN Mother age 15-16, completed primary
//    CMG_04,                                     //EN Mother age 15-16, primary dropout
//    CMG_05,                                     //EN Mother age 15-16, neve entered primary
//    CMG_06,                                     //EN Mother age <15, completed primary
//    CMG_07,                                     //EN Mother age <15, primary dropout
//    CMG_08                                      //EN Mother age <15, neve entered primary
//};

//range CHILD_MORTALITY_AGE{ 0,4 };               //EN Age

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters 
//{
//    //EN Child mortality baseline risks
//    double  ChildMortalityBaseRisk[CHILD_MORTALITY_AGE][SEX];
//
//    //EN Child mortality relative risks
//    double  ChildMortalityRelativeRisks[CHILD_MORTALITY_AGE][CHILD_MORTALITY_GROUP];
//
//    //EN Child mortality time trend
//    double ChildMortalityTrend[CHILD_MORTALITY_AGE][CHILD_MORTALITY_YEARS];
//};

//parameter_group PG02_MortalityDetailed                          //EN Mortality detailed model
//{
//    ChildMortalityBaseRisk,
//    ChildMortalityRelativeRisks,
//    ChildMortalityTrend
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor states and functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    //EN Baseline child mortality
//    double child_mortality =
//        (integer_age == 0 && sex == MALE) ? double(lCalibrator->mort_male_0) :
//        (integer_age == 1 && sex == MALE) ? double(lCalibrator->mort_male_1) :
//        (integer_age == 2 && sex == MALE) ? double(lCalibrator->mort_male_2) :
//        (integer_age == 3 && sex == MALE) ? double(lCalibrator->mort_male_3) :
//        (integer_age == 4 && sex == MALE) ? double(lCalibrator->mort_male_4) :
//        (integer_age == 0 && sex == FEMALE) ? double(lCalibrator->mort_female_0) :
//        (integer_age == 1 && sex == FEMALE) ? double(lCalibrator->mort_female_1) :
//        (integer_age == 2 && sex == FEMALE) ? double(lCalibrator->mort_female_2) :
//        (integer_age == 3 && sex == FEMALE) ? double(lCalibrator->mort_female_3) :
//        (integer_age == 4 && sex == FEMALE) ? double(lCalibrator->mort_female_4) : 0.0;
//
//    //EN Child mortality risk group
//    CHILD_MORTALITY_GROUP child_mortality_group =
//        (mother_age_at_birth >= 17 && educ_mother == EOL_HIGH) ? CMG_00 :   // 17+ completed
//        (mother_age_at_birth >= 17 && educ_mother == EOL_MEDIUM) ? CMG_01 : // 17+ dropout
//        (mother_age_at_birth >= 17 && educ_mother == EOL_LOW) ? CMG_02 :    // 17+ non
//        (mother_age_at_birth >= 15 && educ_mother == EOL_HIGH) ? CMG_03 :   // 15-16 completed
//        (mother_age_at_birth >= 15 && educ_mother == EOL_MEDIUM) ? CMG_04 : // 15-16 dropout
//        (mother_age_at_birth >= 15 && educ_mother == EOL_LOW) ? CMG_05 :    // 15-16 non
//        (educ_mother == EOL_HIGH) ? CMG_06 :                                // <15 completed
//        (educ_mother == EOL_MEDIUM) ? CMG_07 :                              // <15 dropout
//        CMG_08;                                                             // <15 non
//
//    double mother_age_at_birth = { TIME_INFINITE };                         //EN Mothers age at own birth
//    event  timeMortalityDetailedEvent, MortalityDetailedEvent;              //EN Mortality event detailed model
//};

//actor Calibrator
//{
//    double mort_male_0 = { 1.0 };                                           //EN Child mortality risk male age 0
//    double mort_male_1 = { 1.0 };                                           //EN Child mortality risk male age 1
//    double mort_male_2 = { 1.0 };                                           //EN Child mortality risk male age 2
//    double mort_male_3 = { 1.0 };                                           //EN Child mortality risk male age 3
//    double mort_male_4 = { 1.0 };                                           //EN Child mortality risk male age 4
//
//    double mort_female_0 = { 1.0 };                                         //EN Child mortality risk female age 0
//    double mort_female_1 = { 1.0 };                                         //EN Child mortality risk female age 1
//    double mort_female_2 = { 1.0 };                                         //EN Child mortality risk female age 2
//    double mort_female_3 = { 1.0 };                                         //EN Child mortality risk female age 3
//    double mort_female_4 = { 1.0 };                                         //EN Child mortality risk female age 4
//
//    logical is_calibrated_child_mortality = { FALSE };                      //EN Child mortality is calibrated
//    event timeChildMortalityCalibration, ChildMortalityCalibration;         //EN Child mortality calibration
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TIME Person::timeMortalityDetailedEvent()
{
    double dEventTime = TIME_INFINITE;
    double dHazard = 0.0;

    if (  is_mortal && gprParam->MortalityModel != MOM_BASE &&  integer_age <= MAX(CHILD_MORTALITY_AGE) 
        &&  calendar_year >= MIN(CHILD_MORTALITY_YEARS) )
    {
        // base risk * trend
        // unaligned model
        if (gprParam->MortalityModel == MOM_DETAIL)
        {
            dHazard = gprParam->ChildMortalityBaseRisk[CHKBND(integer_age ,5,"ChildMortalityBaseRisk",0)][CHKBND(sex ,2,"ChildMortalityBaseRisk",1)]
                * gprParam->ChildMortalityTrend[CHKBND(integer_age ,5,"ChildMortalityTrend",0)][CHKBND(RANGE_POS(CHILD_MORTALITY_YEARS, calendar_year ),146,"ChildMortalityTrend",1)];
        }

        // aligned with overall trend
        else if (gprParam->MortalityModel == MOM_ALIGNED_MACRO_TRENDS)
        {
            dHazard = 
                child_mortality * gprParam->MortalityTrend[CHKBND(RANGE_POS(SIM_YEAR_RANGE, calendar_year ),151,"MortalityTrend",0)][CHKBND(sex ,2,"MortalityTrend",1)]
                / gprParam->MortalityTrend[CHKBND(RANGE_POS(SIM_YEAR_RANGE, MIN(CHILD_MORTALITY_YEARS)),151,"MortalityTrend",0)][CHKBND(sex ,2,"MortalityTrend",1)];
        }

        // aligned with specific child mortality trend
        else
        {
            dHazard = 
                child_mortality * gprParam->ChildMortalityTrend[CHKBND(integer_age ,5,"ChildMortalityTrend",0)][CHKBND(RANGE_POS(CHILD_MORTALITY_YEARS, calendar_year ),146,"ChildMortalityTrend",1)];
        }

        // relativ risks
        dHazard = dHazard * gprParam->ChildMortalityRelativeRisks[CHKBND(integer_age ,5,"ChildMortalityRelativeRisks",0)][CHKBND(child_mortality_group ,9,"ChildMortalityRelativeRisks",1)];

        // if hazard is positive calculate event time
        if (dHazard > 0) dEventTime = WAIT(-log(RandUniform(9)) / dHazard);;
    }
    return dEventTime;
}

void Person::MortalityDetailedEvent() { Finish(); }


TIME Calibrator::timeChildMortalityCalibration()
{
    if (!is_calibrated_child_mortality ) return MIN(CHILD_MORTALITY_YEARS);
    else return TIME_INFINITE;
}

void Calibrator::ChildMortalityCalibration()
{
    // local matrices
    double dDeaths[SIZE(CHILD_MORTALITY_AGE)][SIZE(SEX)];      // Number expected deaths
    double dProb[SIZE(CHILD_MORTALITY_AGE)][SIZE(SEX)];        // Death probability
    double dBase[SIZE(CHILD_MORTALITY_AGE)][SIZE(SEX)];        // Baseline Hazard in simulation
    long   nPop[SIZE(CHILD_MORTALITY_AGE)][SIZE(SEX)][SIZE(CHILD_MORTALITY_GROUP)]; // Pop sizes

    // Initialize matrices
    // set population sizes nPop and expected deaths nDeaths to 0
    // sets death probabilities dProb by age and sex for the year in which the model is calibrated
    for (int nAge = 0; nAge < SIZE(CHILD_MORTALITY_AGE); nAge++)
    {
        for (int nSex = 0; nSex < SIZE(SEX); nSex++)
        {
            dProb[nAge][nSex] = 1.0 - exp(-gprParam->MortalityTable[CHKBND(nAge,101,"MortalityTable",0)][CHKBND(nSex,2,"MortalityTable",1)]
                * gprParam->MortalityTrend[CHKBND(RANGE_POS(SIM_YEAR_RANGE, MIN(CHILD_MORTALITY_YEARS)),151,"MortalityTrend",0)][CHKBND(nSex,2,"MortalityTrend",1)]);
            dDeaths[nAge][nSex] = 0.0;
            for (int nGroup = 0; nGroup < SIZE(CHILD_MORTALITY_GROUP); nGroup++)
            {
                nPop[nAge][nSex][nGroup] = 0.0;
            }
        }
    }

    // Population Count
    // calculates population sizes nPop by age, sex, risk group
    // calculates expected deaths dDeaths by age and sex
    for (long nJ = 0; nJ < asAllPerson->Count(); nJ++)
    {
        auto prPerson = asAllPerson->Item(nJ);
        if (prPerson-> integer_age <= MAX(CHILD_MORTALITY_AGE) && prPerson->creation_type ==CT_BIRTH && prPerson->is_resident )
        {
            dDeaths[prPerson->integer_age ][prPerson->sex ] = dDeaths[prPerson->integer_age ][prPerson->sex ] + dProb[prPerson->integer_age ][prPerson->sex ];
            nPop[prPerson->integer_age ][prPerson->sex ][prPerson->child_mortality_group ]++;
        }
    }

    // find calibrated baselines
    for (int nAge = 0; nAge < SIZE(CHILD_MORTALITY_AGE); nAge++)
    {
        for (int nSex = 0; nSex < SIZE(SEX); nSex++)
        {
            double  dUpper = 2.0;
            double  dLower = 0.0;
            double  dCenter = 1.0;
            double  dNumberDeaths = 0.0;
            int     nIterations = 10000;

            while ( abs(dNumberDeaths - dDeaths[nAge][nSex]) > 0.0001 && nIterations > 0 )
            {
                nIterations--;
                dBase[nAge][nSex] = ( dLower + dUpper) / 2.0;

                dNumberDeaths = 0.0;

                //Celculate numer of deaths for given dBase
                for (int nGroup = 0; nGroup < SIZE(CHILD_MORTALITY_GROUP); nGroup++)
                {
                    dNumberDeaths = dNumberDeaths + nPop[nAge][nSex][nGroup]
                        * (1-exp(-dBase[nAge][nSex] *  gprParam->ChildMortalityRelativeRisks[CHKBND(nAge,5,"ChildMortalityRelativeRisks",0)][CHKBND(nGroup,9,"ChildMortalityRelativeRisks",1)]));
                }
                // shrink search interval
                if ( dNumberDeaths > dDeaths[nAge][nSex]) dUpper = dBase[nAge][nSex];
                else dLower = dBase[nAge][nSex];
            }
        }
    }

    // set states
     poStateFunctions->Set_mort_male_0(  dBase[0][MALE] );  poStateFunctions->Set_mort_female_0(  dBase[0][FEMALE] );
     poStateFunctions->Set_mort_male_1(  dBase[1][MALE] );  poStateFunctions->Set_mort_female_1(  dBase[1][FEMALE] );
     poStateFunctions->Set_mort_male_2(  dBase[2][MALE] );  poStateFunctions->Set_mort_female_2(  dBase[2][FEMALE] );
     poStateFunctions->Set_mort_male_3(  dBase[3][MALE] );  poStateFunctions->Set_mort_female_3(  dBase[3][FEMALE] );
     poStateFunctions->Set_mort_male_4(  dBase[4][MALE] );  poStateFunctions->Set_mort_female_4(  dBase[4][FEMALE] );

     poStateFunctions->Set_is_calibrated_child_mortality(  TRUE );
}

} // namespace 
