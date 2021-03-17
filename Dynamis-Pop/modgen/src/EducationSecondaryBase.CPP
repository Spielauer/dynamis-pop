#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/EducationSecondaryBase.mpp"
//LABEL (EducationSecondaryBase, EN)
/* NOTE(EducationSecondaryBase, EN) 
 
 TODO adapt text

This is an optional module for modeling a "School Two" which can be attended after graduating from the"School One" - in the current implementation this school refers tw a secondary school. The number of grades can be set in the range EDUC_TWO_GRADE and thus is adaptable for the application context. The model follows a period approach. Each year at the end and beginning of school season it is decided if a person enters school, succeeds a grade, progresses to the next or repeats a grade, interrups studies or permanently drops out. A parameter controls the maximum number of years which can be delayed due to repetition or temporary dropout. A student interrupting education can resume studies after a year or stay out permanently. The model is driven by parameters for each calendar year and grade:

- The probability to pass a grade
- The probability to move on directly after passing a grade
- The probability to repeat a grade immediately after failing it
- The probability to resume studies after being out for a year after passing a grade
- The pribability to repeat a grade after being out for a year after failing a grade
- The maximum number of years studies can be delayed by repetition or inactive spells

Besides grade and period, the school progression parameters have two additional dimensions, one for geographical region, the other for personal characteristics. The levels of these two characteristics are model specific, they are declared in the classifications EDUC_TWO_GEO and EDUC_TWO_GROUP, with the corresponding states educ_two_geo and educ_two_group. Currently only totals are implemented.

As an optional module the module can be added or removed from the model. The only required code change in other modules is the creation of the actor SchoolTwoActor in the Simulation() function. 
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor Sets
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//EN Potential Students of School Two
//actor_set Person asPotentialSchoolTwoStudent
//    filter educ_one_grade_passed == MAX(EDUC_ONE_GRADE) && educ_two_status != ETS_OUT && educ_two_delay <= Educ2AllowedDelays;

//EN Students of School Two To Be Processed
//actor_set Person asPotentialSchoolTwoStudentToProcess filter is_alive && educ_two_to_process;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Types
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//range EDUC_TWO_GRADE { 1, 6 };          //EN Grades School Two
//range EDUC_TWO_GRADE0 { 0 , 6 };        //EN Grades School Two (incl. 0)
//range EDUC_TWO_ENTRY_DUR { 0, 2 };      //EN Possible entry points years since primary 

//classification EDUC_TWO_STATUS          //EN Current attendance status in School Two
//{
//    ETS_NEVER,                          //EN Never entered
//    ETS_WAIT,                           //EN Wait for entry (e.g. school break)
//    ETS_ATTEND,                         //EN Currently attending
//    ETS_PAUSE,                          //EN Pause
//    ETS_OUT                             //EN Left school
//};


//classification EDUC_TWO_LEVEL           //EN Secondary Education Level
//{
//    ETL_NEVER,                          //EN Never entered
//    ETL_INCOMPLETE,                     //EN Incomplete 
//    ETL_GRADUATE                        //EN Graduate
//};


//classification EDUC_TWO_GEO             //EN Regional Unit
//{
//    ETG_OO                              //EN All
//};

//classification EDUC_TWO_GROUP           //EN Population Group
//{
//    E2G_OO                              //EN All
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Prameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters
//{
//    //EN Lower Secondary Period Direct Progression Intake
//    double  Educ2DirectProgressionIntake[EDUC_TWO_GEO][EDUC_TWO_GROUP][EDUC_TWO_GRADE][SIM_YEAR_RANGE];  
//    
//    //EN Lower Secondary Period Delayed Progression Intake
//    double  Educ2DelayedProgressionIntake[EDUC_TWO_GEO][EDUC_TWO_GROUP][EDUC_TWO_GRADE][SIM_YEAR_RANGE]; 
//
//    //EN Lower Secondary Period Success
//    double  Educ2PeriodSuccess[EDUC_TWO_GEO][EDUC_TWO_GROUP][EDUC_TWO_GRADE][SIM_YEAR_RANGE];                  
//    
//    //EN Lower Secondary Period Direct Repetition Intake
//    double  Educ2DirectRepetitionIntake[EDUC_TWO_GEO][EDUC_TWO_GROUP][EDUC_TWO_GRADE][SIM_YEAR_RANGE];   
//    
//    //EN Lower Secondary Period Delayed Repetition Intake
//    double  Educ2DelayedRepetitionIntake[EDUC_TWO_GEO][EDUC_TWO_GROUP][EDUC_TWO_GRADE][SIM_YEAR_RANGE];  
//
//    //EN Maximum years of interruptions and repetition of grades
//    int     Educ2AllowedDelays;
//
//    model_generated double  StartSchoolTwoYear;                     //EN Start of school year
//    model_generated double  EndSchoolTwoYear;                       //EN End of school year 
//};

//parameter_group PG12_EducTwo                                        //EN Secondary School
//{
//    Educ2DirectProgressionIntake, Educ2DelayedProgressionIntake, Educ2PeriodSuccess,
//    Educ2DirectRepetitionIntake, Educ2DelayedRepetitionIntake, Educ2AllowedDelays
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actors
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    int                 educ_two_delay = { 0 };                     //EN Delays accumulated
//    logical             educ_two_to_process = { FALSE };            //EN Secondary to be processed flag
//    EDUC_TWO_STATUS     educ_two_status = { ETS_NEVER };            //EN Secondary education status
//
//    //EN Secondary Education Level
//    EDUC_TWO_LEVEL      educ_two_level = (educ_two_status == ETS_NEVER) ? ETL_NEVER :
//                        (educ_two_grade_passed == MAX(EDUC_TWO_GRADE)) ? ETL_INCOMPLETE : ETL_GRADUATE;
//
//
//    EDUC_TWO_GRADE0     educ_two_grade_attended = { 0 };            //EN Highest grade attended
//    EDUC_TWO_GRADE0     educ_two_grade_passed = { 0 };              //EN Highest grade passed
//    EDUC_TWO_GEO        educ_two_geo = {ETG_OO};                    //EN Geographical location
//    EDUC_TWO_GROUP      educ_two_group = {E2G_OO};                  //EN Person group
//};

//actor Clock
//{
//    TIME time_start_school_two_year = { TIME_INFINITE};             //EN Next start of school year
//    TIME time_end_school_two_year = { TIME_INFINITE};               //EN Next end of school year    
//
//    event timeStartSchoolTwoYearEvent, StartSchoolTwoYearEvent;     //EN Start school year event
//    event timeEndSchoolTwoYearEvent, EndSchoolTwoYearEvent;         //EN End school year event
//    void  SetSchoolYearTwoClock();                                  //EN Clock for scheduling school year
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation Clock Functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TIME Clock::timeStartSchoolTwoYearEvent(){return time_start_school_two_year ;}
void Clock::StartSchoolTwoYearEvent()
{
    for (long nJ = 0; nJ < asPotentialSchoolTwoStudent->Count(); nJ++)
    {
        auto prPerson = asPotentialSchoolTwoStudent->Item(nJ);
        if (prPerson-> educ_two_status == ETS_WAIT )
        {
            prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_ATTEND  );
            prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_grade_attended(  prPerson-> educ_two_grade_passed + 1 );
        }
    }
     poStateFunctions->Set_time_start_school_two_year(  TIME_INFINITE );
}

TIME Clock::timeEndSchoolTwoYearEvent(){return time_end_school_two_year ;}
void Clock::EndSchoolTwoYearEvent()
{
    // Flag the persons to be processed
    for (long nJ = 0; nJ < asPotentialSchoolTwoStudent->Count(); nJ++)
    {
        auto prPerson = asPotentialSchoolTwoStudent->Item(nJ);
        prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_to_process(  TRUE );
    }

    // Process all potential students
    for (long nJ = 0; nJ < asPotentialSchoolTwoStudentToProcess->Count(); nJ++)
    {
        auto prPerson = asPotentialSchoolTwoStudentToProcess->Item(nJ);

        int nGeo         = prPerson->educ_two_geo ;
        int nGroup       = prPerson->educ_two_group ;
        int nGradeAttend = RANGE_POS(EDUC_TWO_GRADE, prPerson->educ_two_grade_attended );
        int nYear        = RANGE_POS(SIM_YEAR_RANGE, prPerson->calendar_year );


        // Potential Student has not entered school yet 
        if (prPerson-> educ_two_status == ETS_NEVER && ( 
            ( prPerson-> educ_two_delay == 0 && RandUniform(68) < gprParam->Educ2DirectProgressionIntake[CHKBND(nGeo,1,"Educ2DirectProgressionIntake",0)][CHKBND(nGroup,1,"Educ2DirectProgressionIntake",1)][CHKBND(0,6,"Educ2DirectProgressionIntake",2)][CHKBND(nYear,151,"Educ2DirectProgressionIntake",3)] ) ||
            ( prPerson-> educ_two_delay == 1 && RandUniform(69) < gprParam->Educ2DelayedProgressionIntake[CHKBND(nGeo,1,"Educ2DelayedProgressionIntake",0)][CHKBND(nGroup,1,"Educ2DelayedProgressionIntake",1)][CHKBND(0,6,"Educ2DelayedProgressionIntake",2)][CHKBND(nYear,151,"Educ2DelayedProgressionIntake",3)] )))
        {
            prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_WAIT );
            prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_delay(  0 );           //reset at school entry
        }
        else if (prPerson-> educ_two_status == ETS_NEVER) prPerson->UpdateTime() -> poStateFunctions->PlusPlusSet_educ_two_delay();

        // Active student in last year
        else if (prPerson-> educ_two_grade_attended == MAX(EDUC_TWO_GRADE) && prPerson-> educ_two_status == ETS_ATTEND)  
        {
            // Graduate
            if ( RandUniform(70) < gprParam->Educ2PeriodSuccess[CHKBND(nGeo,1,"Educ2PeriodSuccess",0)][CHKBND(nGroup,1,"Educ2PeriodSuccess",1)][CHKBND(nGradeAttend,6,"Educ2PeriodSuccess",2)][CHKBND(nYear,151,"Educ2PeriodSuccess",3)] )  
            {
                prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_grade_passed(  prPerson->educ_two_grade_attended  );
                prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_OUT );
            }
            // Fail and repeat if still time
            else if (RandUniform(71) < gprParam->Educ2DirectRepetitionIntake[CHKBND(nGeo,1,"Educ2DirectRepetitionIntake",0)][CHKBND(nGroup,1,"Educ2DirectRepetitionIntake",1)][CHKBND(nGradeAttend,6,"Educ2DirectRepetitionIntake",2)][CHKBND(nYear,151,"Educ2DirectRepetitionIntake",3)] && prPerson-> educ_two_delay < gprParam->Educ2AllowedDelays) 
            {
                prPerson->UpdateTime() -> poStateFunctions->PlusPlusSet_educ_two_delay();
                prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_WAIT );
            }
            // Fail and wait
            else if ( prPerson-> educ_two_delay + 1 < gprParam->Educ2AllowedDelays ) 
            {
                prPerson->UpdateTime() -> poStateFunctions->PlusPlusSet_educ_two_delay();
                prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_PAUSE );
            }
            // Fail and finally out
            else 
            {
                prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_OUT );
            }
        }

        // Active student before last year
        else if (prPerson-> educ_two_grade_attended < MAX(EDUC_TWO_GRADE) && prPerson-> educ_two_status == ETS_ATTEND)  
        {
            // Decide grade success
            if ( RandUniform(72) < gprParam->Educ2PeriodSuccess[CHKBND(nGeo,1,"Educ2PeriodSuccess",0)][CHKBND(nGroup,1,"Educ2PeriodSuccess",1)][CHKBND(nGradeAttend,6,"Educ2PeriodSuccess",2)][CHKBND(nYear,151,"Educ2PeriodSuccess",3)] )  
            {
                prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_grade_passed(  prPerson->educ_two_grade_attended  );
            }
            // Success
            if ( prPerson-> educ_two_grade_passed == prPerson-> educ_two_grade_attended )
            {
                if (RandUniform(73) < gprParam->Educ2DirectProgressionIntake[CHKBND(nGeo,1,"Educ2DirectProgressionIntake",0)][CHKBND(nGroup,1,"Educ2DirectProgressionIntake",1)][CHKBND(nGradeAttend+1,6,"Educ2DirectProgressionIntake",2)][CHKBND(nYear,151,"Educ2DirectProgressionIntake",3)]) 
                {
                    prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_WAIT );
                }
                else if (prPerson-> educ_two_delay < gprParam->Educ2AllowedDelays)
                {
                    prPerson->UpdateTime() -> poStateFunctions->PlusPlusSet_educ_two_delay();
                    prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_PAUSE );
                }
                else prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_OUT );
            }
            // no success
            else
                if (RandUniform(74) < gprParam->Educ2DirectRepetitionIntake[CHKBND(nGeo,1,"Educ2DirectRepetitionIntake",0)][CHKBND(nGroup,1,"Educ2DirectRepetitionIntake",1)][CHKBND(nGradeAttend,6,"Educ2DirectRepetitionIntake",2)][CHKBND(nYear,151,"Educ2DirectRepetitionIntake",3)] 
                    && prPerson-> educ_two_delay < gprParam->Educ2AllowedDelays) 
                {
                    prPerson->UpdateTime() -> poStateFunctions->PlusPlusSet_educ_two_delay();
                    prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_WAIT );
                }
                else if (prPerson-> educ_two_delay + 1 < gprParam->Educ2AllowedDelays)
                {
                    prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_delay(  prPerson-> educ_two_delay + 2 );
                    prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_PAUSE );
                }
                else prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_OUT );
        }
        // Previousely inactive
        else
        {
            // Progression
            if ( prPerson-> educ_two_grade_attended == prPerson-> 
                educ_two_grade_passed && RandUniform(75) < gprParam->Educ2DelayedProgressionIntake[CHKBND(nGeo,1,"Educ2DelayedProgressionIntake",0)][CHKBND(nGroup,1,"Educ2DelayedProgressionIntake",1)][CHKBND(nGradeAttend+1,6,"Educ2DelayedProgressionIntake",2)][CHKBND(nYear,151,"Educ2DelayedProgressionIntake",3)] ) 
            {
                prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_WAIT );
            }
            // Repetition
            else if (RandUniform(76) < gprParam->Educ2DelayedRepetitionIntake[CHKBND(nGeo,1,"Educ2DelayedRepetitionIntake",0)][CHKBND(nGroup,1,"Educ2DelayedRepetitionIntake",1)][CHKBND(nGradeAttend,6,"Educ2DelayedRepetitionIntake",2)][CHKBND(nYear,151,"Educ2DelayedRepetitionIntake",3)] ) 
            {
                prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_WAIT );
            }
            else prPerson-> UpdateTime() -> poStateFunctions->Set_educ_two_status(  ETS_OUT );
        }
    }

    // Unflag all persons
    while (asPotentialSchoolTwoStudentToProcess->Count() > 0)
    {
        asPotentialSchoolTwoStudentToProcess->Item(0)-> UpdateTime() -> poStateFunctions->Set_educ_two_to_process(  FALSE );
    }
     poStateFunctions->Set_time_end_school_two_year(  TIME_INFINITE );
}

void Clock::SetSchoolYearTwoClock()
{
     poStateFunctions->Set_time_start_school_two_year(  WAIT(gprParam->StartSchoolTwoYear) );
     poStateFunctions->Set_time_end_school_two_year(  WAIT(gprParam->EndSchoolTwoYear) );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pre-Simulation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PreSimulation1()
{
    gprParam->StartSchoolTwoYear = gprParam->StartSchoolOneYear + 0.000114155;  // one hour later than school one
    gprParam->EndSchoolTwoYear = gprParam->EndSchoolOneYear + 0.000114155;      // one hour later than school one
};



} // namespace 
