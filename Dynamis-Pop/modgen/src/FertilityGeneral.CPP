#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/FertilityGeneral.mpp"
//LABEL(FertilityGeneral, EN) Fertility General

/* NOTE(FertilityGeneral, EN)
This module contains the general functionality of fertility
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dimensions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//classification FERTILITY_MODEL                  //EN Fertility model selection
//{
//    FEM_BASE,                                   //EN Base model
//    FEM_DETAIL,                                 //EN Detailed model
//    FEM_ALIGNED_AGE,                            //EN Detailed model aligned to base by age
//    FEM_ALIGNED_TOTAL                           //EN Detailed model aligned to base by totals
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters 
//{
//    FERTILITY_MODEL FertilityModel;             //EN Fertility model selection
//    double          SexRatio[SIM_YEAR_RANGE];   //EN Sex Ratio (males per 100 females)
//};

//parameter_group PG01_Fertility                  //EN Fertility
//{
//    FertilityModel,
//    SexRatio,
//    PG01_FertilityBase,
//    PG01_FertilityDetailed
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor states and functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person            
//{
//    //EN Person is at risk of giving birth
//    logical is_fertile = (sex == FEMALE && WITHIN(FERTILE_AGE_RANGE, integer_age) && (creation_type == CT_BIRTH
//                     || (creation_type == CT_START && calendar_year >= MIN(SIM_YEAR_RANGE))
//                     || (creation_type == CT_POOL && ever_resident)
//                     || (creation_type == CT_SCRATCH && ever_resident))) ? TRUE : FALSE;
//    
//    FERTILE_AGE_RANGE fertile_age = COERCE(FERTILE_AGE_RANGE, integer_age);     //EN Age
//    
//    void HandleFertility();             //EN Handle Fertility
//    void doGiveBirth();                 //EN Give Birth
//
//    logical this_parity_spell = { FALSE };                              //EN In this parity spell switch
//};

//actor StartpopValues
//{
//    // Setting clock for last birth imputation event
//    logical is_set_last_birth_clock = { FALSE };                        //EN Clock is set
//    event   timeSetLastBirthClock, SetLastBirthClock;                   //EN Set last birth imputation clock
//
//    TIME time_last_birth_imputation = { TIME_INFINITE };                //EN Time of last birth imputation
//    event timeLastBirthImputationEvent, LastBirthImputationEvent;       //EN First union in past
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TIME StartpopValues::timeSetLastBirthClock()
{
    if ( is_activated && !is_set_last_birth_clock ) return WAIT(0);
    else return TIME_INFINITE;
}

void StartpopValues::SetLastBirthClock()
{
    TIME dEventTime = TIME_INFINITE;
    if (lPersonStartpop -> sex == FEMALE
        && (lPersonStartpop -> creation_type == CT_START || lPersonStartpop -> creation_type == CT_POOL)
        && StartPopValue[CHKBND(PMC_LASTBIR,18,"StartpopValues.StartPopValue",0)] > StartPopValue[CHKBND(PMC_BIRTH,18,"StartpopValues.StartPopValue",0)]
        && StartPopValue[CHKBND(PMC_LASTBIR,18,"StartpopValues.StartPopValue",0)] < MIN(SIM_YEAR_RANGE))
    {
        dEventTime = StartPopValue[CHKBND(PMC_LASTBIR,18,"StartpopValues.StartPopValue",0)];
        if (dEventTime == int(dEventTime)) dEventTime = dEventTime + RandUniform(54);
        //add offset
        dEventTime = dEventTime + int(lPersonStartpop ->time_of_birth ) - StartPopValue[CHKBND(PMC_BIRTH,18,"StartpopValues.StartPopValue",0)];
    }

     poStateFunctions->Set_time_last_birth_imputation(  dEventTime );
     poStateFunctions->Set_is_set_last_birth_clock(  TRUE );
}

TIME StartpopValues::timeLastBirthImputationEvent() { return time_last_birth_imputation ; }
void StartpopValues::LastBirthImputationEvent()
{
    lPersonStartpop -> UpdateTime() -> poStateFunctions->Set_this_parity_spell(  FALSE );
    if (StartPopValue[CHKBND(PMC_PARITY,18,"StartpopValues.StartPopValue",0)] > lPersonStartpop ->parity ) lPersonStartpop -> UpdateTime() -> poStateFunctions->Set_parity(  COERCE(PARITY_RANGE,int(StartPopValue[CHKBND(PMC_PARITY,18,"StartpopValues.StartPopValue",0)])) );
    lPersonStartpop -> UpdateTime() -> poStateFunctions->Set_this_parity_spell(  TRUE );
     poStateFunctions->Set_time_last_birth_imputation(  TIME_INFINITE );
}

void Person::HandleFertility()
{
    // have a child yourself if no alignment
    if ( gprParam->FertilityModel == FEM_BASE || gprParam->FertilityModel == FEM_DETAIL ) doGiveBirth();

    // find person with shortest waiting time to birth
    else lCalibrator ->UpdateTime()->GetNextToGiveBirth(this)->UpdateTime()->doGiveBirth();
}

void Person::doGiveBirth()
{
    auto peChild = new Person;  
    
    doIncreaseParity();
        
    peChild->UpdateTime()->Start(NULL, this, TIME_INFINITE, (SEX)0);   
    doLeaveParentalHome();
}
} // namespace 
