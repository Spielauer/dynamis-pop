#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/FamilyFirstUnion.mpp"
//LABEL (FamilyFirstUnion, EN) First union formation
/* NOTE(FamilyFirstUnion, EN) 
This module implements two alternative ways of modelling women's first union formation:

- A parametric Coale McNeil model by population groups. 
- A model based on first union hazard rates by population group.

Two events handle first union events:

- For persons created form starting population file information (including population pools of immigrants), the date of first union formation in the starting population record is used to impute the union formation in the past. For efficiency reasons, the event FirstUnionImputationEvent() is handled by the temporary StartpopValues actor.
- The event FirstUnionFormationEvent() of the actor Person handles all other first union formations. It is based on one of the two alternative ways of modelling first union formation. These events happen in projected time and – in the case of immigrants from immigration pools – once entering the country. For immigrants created from scratch, all first union formation events (including those before starting the simulation or before entering the country) are handled by this process. 

**The Coale McNeil model:**

The Coale McNeil model for first union formation of women allows parameterizing union formation by three very intuitive parameters:

- The youngest age of first union formation (this is where about the first 1% have entered a first union)
- The average age of first union formation
- The proportion of the population ever entering a union

The parameterization of first union formation is by year of birth and population groups UNION1_GROUP. The classification is generic; in the current implementation, it corresponds with the primary education fate. The list can be extended or changed as long as the derived state union1_group which stores the individual group membership is adapted accordingly.

Internally the Coale McNeil model calculates period union formation rates by calendar year, age and population group in the pre-simulation phase. 

**First union formation rates:**

The Coale McNeil model typically gives a good fit for developing countries and allows easy scenario building. At the other hand, it might perform poorly for developed countries, therefore losing its attractivity as countries develop. This module's technical implementation allows the use of an alternative model that substitutes the period rates calculated in pre-simulation by rates obtained from other sources.
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Types
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//classification UNION1_CHOICE                //EN First union model choice
//{
//    U1C_CMN,                                //EN Coale McNeil
//    U1C_HAZARDS                             //EN Hazards
//};

//classification UNION1_GROUP                 //EN Population Groups
//{
//    U1G_00,                                 //EN Never entered school
//    U1G_01,                                 //EN Primary dropout
//    U1G_02                                  //EN Primary graduate
//};

//classification UNION1_PARA                  //EN First union formation parameters
//{
//    UP_MINAGE,                              //EN Age at which first 1% enter a union
//    UP_AVERAGE,                             //EN Average age at first union formation
//    UP_EVER                                 //EN Proportion ow femailes ever entering a union
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters
//{
//    //EN First union model choice
//    UNION1_CHOICE Union1Choice;
//
//    //EN Union formation Coale McNeil model
//    double Union1ParametersCMN[UNION1_GROUP][YOB_UNION][UNION1_PARA];
//
//    //EN Union formation hazard model
//    double Union1ParametersHazards[UNION1_GROUP][FEMALE_SPOUSE_AGE][ALL_YEAR_RANGE];
//
//    //EN Union formation hazards
//    model_generated double Union1FormationHazard[UNION1_GROUP][ALL_YEAR_RANGE][FEMALE_SPOUSE_AGE];
//};

//parameter_group PG06_FirstUnion                                         //EN First Union
//{
//    Union1Choice,
//    Union1ParametersCMN,
//    Union1ParametersHazards
//};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor states and functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    //EN Population group 
//    UNION1_GROUP union1_group = (educ_one_fate == EOL_LOW ) ? U1G_00 : (educ_one_fate == EOL_MEDIUM ) ? U1G_01 : U1G_02;
//
//    logical		ever_union = { FALSE };                                 //EN Ever had a partner
//    event		timeFirstUnionFormationEvent, FirstUnionFormationEvent; //EN First union formation event
//};

//actor StartpopValues
//{
//    // Setting clock for first union imputation event
//    logical is_set_first_union_clock = { FALSE };                       //EN Clock is set
//    event   timeSetFirstUnionClock, SetFirstUnionClock;                 //EN Set first union imputation clock
//
//    TIME time_first_union_imputation = { TIME_INFINITE };               //EN Time of first union imputation
//    event timeFirstUnionImputationEvent, FirstUnionImputationEvent;     //EN First union in past
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TIME StartpopValues::timeSetFirstUnionClock()
{
    if ( is_activated && !is_set_first_union_clock ) return WAIT(0);
    else return TIME_INFINITE;
}

void StartpopValues::SetFirstUnionClock()
{
    TIME dEventTime = TIME_INFINITE;
    if (lPersonStartpop -> sex == FEMALE
        && (lPersonStartpop -> creation_type == CT_START || lPersonStartpop -> creation_type == CT_POOL)
        && StartPopValue[CHKBND(PMC_UNION,18,"StartpopValues.StartPopValue",0)] > StartPopValue[CHKBND(PMC_BIRTH,18,"StartpopValues.StartPopValue",0)]
        && StartPopValue[CHKBND(PMC_UNION,18,"StartpopValues.StartPopValue",0)] < MIN(SIM_YEAR_RANGE))
    {
        dEventTime = StartPopValue[CHKBND(PMC_UNION,18,"StartpopValues.StartPopValue",0)];
        if (dEventTime == int(dEventTime)) dEventTime = dEventTime + RandUniform(30);
        //add offset
        dEventTime = dEventTime + int(lPersonStartpop ->time_of_birth ) - StartPopValue[CHKBND(PMC_BIRTH,18,"StartpopValues.StartPopValue",0)];
    }

     poStateFunctions->Set_time_first_union_imputation(  dEventTime );
     poStateFunctions->Set_is_set_first_union_clock(  TRUE );
}

TIME StartpopValues::timeFirstUnionImputationEvent()
{
    return time_first_union_imputation ;
}

void StartpopValues::FirstUnionImputationEvent()
{
    lPersonStartpop -> UpdateTime() -> poStateFunctions->Set_ever_union(  TRUE );
     poStateFunctions->Set_time_first_union_imputation(  TIME_INFINITE );
}

TIME Person::timeFirstUnionFormationEvent()
{
    if (  sex == FEMALE && ! ever_union && WITHIN(FEMALE_SPOUSE_AGE, integer_age ) && ( is_mortal ||  creation_type == CT_SCRATCH ))
    {
        double dHazard = gprParam->Union1FormationHazard[CHKBND(union1_group ,3,"Union1FormationHazard",0)][CHKBND(RANGE_POS(ALL_YEAR_RANGE, calendar_year ),252,"Union1FormationHazard",1)][CHKBND(RANGE_POS(FEMALE_SPOUSE_AGE, integer_age ),51,"Union1FormationHazard",2)];
        if (dHazard > 0.0) return WAIT(-TIME(log(RandUniform(31)) / dHazard));
        else return TIME_INFINITE;
    }
    else return TIME_INFINITE;
}

void Person::FirstUnionFormationEvent()
{
     poStateFunctions->Set_ever_union(  TRUE );
    FindSpouse();
}


/*
GetUnionFormationHazard() returns the union formation hazard for given education, age and period based
on a Coale McNeil model. This is a parametric model with three parameters
*/
double GetUnionFormationHazard(int nEduc, int nYear, int nAge);
double GetUnionFormationHazard(int nEduc, int nYear, int nAge)
{
    double	dReturnValue = 0.0;
    double	g[SIZE(FEMALE_SPOUSE_AGE)];
    double	a0 = gprParam->Union1ParametersCMN[CHKBND(nEduc,3,"Union1ParametersCMN",0)][CHKBND(RANGE_POS(YOB_UNION, nYear - nAge),87,"Union1ParametersCMN",1)][CHKBND(UP_MINAGE,3,"Union1ParametersCMN",2)];
    double	my = gprParam->Union1ParametersCMN[CHKBND(nEduc,3,"Union1ParametersCMN",0)][CHKBND(RANGE_POS(YOB_UNION, nYear - nAge),87,"Union1ParametersCMN",1)][CHKBND(UP_AVERAGE,3,"Union1ParametersCMN",2)];
    double	C = gprParam->Union1ParametersCMN[CHKBND(nEduc,3,"Union1ParametersCMN",0)][CHKBND(RANGE_POS(YOB_UNION, nYear - nAge),87,"Union1ParametersCMN",1)][CHKBND(UP_EVER,3,"Union1ParametersCMN",2)];
    double	k = (my - a0) / 11.36;
    double dThisYearCumulative = 0.0;
    if (nAge < a0) dReturnValue = 0.0;
    else if ( nAge >= a0 )
    {
        for (int nX = 0; nX < SIZE(FEMALE_SPOUSE_AGE); nX++)
        {
            g[nX] = 0.19465*exp(-0.174*(nX - 6.06) - exp(-0.288*(nX - 6.06)));
        }
        double	index = ((double)nAge - a0) / k;
        int		nIndex = int(index);		// integer part
        double	dIndex = index - nIndex;		// after comma part
        double G = 0.0;
        for (int nI = 0; nI <= nIndex; nI++)
        {
            if (nI < SIZE(FEMALE_SPOUSE_AGE)) G = G + g[nI];
        }
        if (nIndex < SIZE(FEMALE_SPOUSE_AGE) - 1) G = G + dIndex*g[nIndex + 1];
        dThisYearCumulative = G * C;
    }
    if (nAge == a0) dReturnValue = dThisYearCumulative;
    else if ( nAge > a0 )
    {
        double	index = ((double)nAge - a0 -1) / k;
        int		nIndex = int(index);		// integer part
        double	dIndex = index - nIndex;		// after comma part
        double G = 0.0;
        for (int nI = 0; nI <= nIndex; nI++)
        {
            if (nI < SIZE(FEMALE_SPOUSE_AGE)) G = G + g[nI];
        }
        if (nIndex < SIZE(FEMALE_SPOUSE_AGE) - 1) G = G + dIndex*g[nIndex + 1];
        if (1.0 - G * C > 0.0001) dReturnValue = (dThisYearCumulative - G * C) / (1 - G * C);
        else dReturnValue = 0.0;
    }
    // convert probability to hazard
    if (dReturnValue < 1.0)  dReturnValue = -log(1 - dReturnValue);
    return dReturnValue;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Pre-Simulation
////////////////////////////////////////////////////////////////////////////////////////////////////

/*
In the presimulation of the union formation module, age specific hazards for union formation
are calculated and copied into a model-generated parameter.
*/
void PreSimulation2()
{
    for (int nEduc = 0; nEduc < SIZE(UNION1_GROUP); nEduc++)
    {
        for (int nYear = MIN(ALL_YEAR_RANGE); nYear <= MAX(ALL_YEAR_RANGE); nYear++)
        {
            for (int nAge = MIN(FEMALE_SPOUSE_AGE); nAge <= MAX(FEMALE_SPOUSE_AGE); nAge++)
            {
                // Coale McNeil
                if (gprParam->Union1Choice == U1C_CMN)
                {
                    gprParam->Union1FormationHazard[CHKBND(nEduc,3,"Union1FormationHazard",0)][CHKBND(RANGE_POS(ALL_YEAR_RANGE, nYear),252,"Union1FormationHazard",1)][CHKBND(RANGE_POS(FEMALE_SPOUSE_AGE, nAge),51,"Union1FormationHazard",2)]
                        = GetUnionFormationHazard(nEduc, nYear, nAge);
                }
                // Hazard model
                else
                {
                    gprParam->Union1FormationHazard[CHKBND(nEduc,3,"Union1FormationHazard",0)][CHKBND(RANGE_POS(ALL_YEAR_RANGE, nYear),252,"Union1FormationHazard",1)][CHKBND(RANGE_POS(FEMALE_SPOUSE_AGE, nAge),51,"Union1FormationHazard",2)]
                        = gprParam->Union1ParametersHazards[CHKBND((UNION1_GROUP)nEduc,3,"Union1ParametersHazards",0)][CHKBND(RANGE_POS(FEMALE_SPOUSE_AGE, nAge),51,"Union1ParametersHazards",1)][CHKBND(RANGE_POS(ALL_YEAR_RANGE, nYear),252,"Union1ParametersHazards",2)];
                }
            }
        }
    }
}

} // namespace 
