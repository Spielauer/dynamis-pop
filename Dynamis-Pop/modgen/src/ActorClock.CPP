#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/ActorClock.mpp"
//LABEL (ActorClock, EN) Definition and initialization of Clock actor

/* NOTE(ActorClock, EN)
This module contains the basic information associated with a Clock actor. The clock actor is used to schedule clock events which affect all persons of the simulation. This prevents clogging the event queue with thousands of individual-level events scheduled at the same time. Currently, the following clocks are implemented:

- Year end: calls a corresponding YearEnd() function of all persons in the simulation. This should be the last function called in a calendar year. At this time, the end of the year is reached, but the calendar year is still not incremented.
- Year start: calls a corresponding YearStart() function of all persons in the simulation. This is the first function called in a calendar year. At this time, the calendar year is set to the new year already.
- Midyear: calls a corresponding MidYear() function of all persons in the simulation.

The Clock actor is the first actor created in the simulation. All other actors are linked to the clock when created.
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor sets
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor_set Clock asClock;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor states and functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Clock //EN Actor for synchronising population-level activities
//{
//    void    Start();                                        //EN Starts the Clock
//    void    Finish();                                       //EN Finishes the Clock
//
//    int     clock_year = { 1900 };                          //EN Calendar Year
//
//    TIME    next_clock_year_end = { 1900 };                 //EN Time of next calendar clock year end
//    event   timeClockYearEndEvent, ClockYearEndEvent;       //EN Clock year end event
//
//    TIME    next_clock_year_start = { 1900 };               //EN Time of next calendar clock year start
//    event   timeClockYearStartEvent, ClockYearStartEvent;   //EN Clock year start event
//
//    TIME    next_midyear_clock_event = { TIME_INFINITE };   //EN Time of next midyear clock event
//    event   timeClockMidyearEvent, ClockMidyearEvent;       //EN Midyear clock event
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Clock Creation and termination
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Clock::Start( )
{ if (__finished) {StartAfterFinish(); return; }; StartSpecial1();
    // Initialize all attributes (OpenM++).
    initialize_attributes();

    // Setting the actor weight (Modgen only)
    Set_actor_weight(gprParam->ScalingFactor); Set_actor_subsample_weight(gprParam->ScalingFactor);

     poStateFunctions->Set_time(  (TIME) MIN(ALL_YEAR_RANGE), true );
     poStateFunctions->Set_clock_year(  MIN(ALL_YEAR_RANGE) );
     poStateFunctions->Set_age(  0, true );
     poStateFunctions->Set_next_clock_year_end(  WAIT(1) );
     poStateFunctions->Set_next_clock_year_start(  TIME_INFINITE );

    // Have the entity enter the simulation (OpenM++).
    enter_simulation();
StartSpecial2();}

void Clock::Finish()
{
    // Have the entity exit the simulation (OpenM++).
    exit_simulation();
if (__finished) {	DoubleFinish(case_seed);	return; } else __finished = TRUE;FinishSpecial();}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation of clock events
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TIME Clock::timeClockMidyearEvent() { return next_midyear_clock_event ; }
void Clock::ClockMidyearEvent()
{
    // Call things happening midyear
    UpdatePartnershipStatus();

    // Call midyear event for all persons
    int nPerson = asAllPerson->Count();
    for ( int nJ = 0; nJ < nPerson; nJ++ ) 
    {
		auto paPerson = asAllPerson->Item( nJ );
        paPerson->UpdateTime()->MidYear();
    }
    // reset clock; next clock event set in ClockYearStartEvent()
     poStateFunctions->Set_next_midyear_clock_event(  TIME_INFINITE );
}


TIME Clock::timeClockYearStartEvent() { return next_clock_year_start ; }
void Clock::ClockYearStartEvent()
{
    poStateFunctions->PlusPlusSet_clock_year();

    // Age all Persons to the current time.
    int nPerson = asAllPerson->Count();
    for ( int nJ = 0; nJ < nPerson; nJ++ ) 
    {
		auto paPerson = asAllPerson->Item( nJ );
        paPerson-> UpdateTime() -> poStateFunctions->Set_calendar_year(  (ALL_YEAR_RANGE)(int)clock_year  );
        paPerson->UpdateTime()->YearStart();
    }

    // Schedule the next year change.
     poStateFunctions->Set_next_clock_year_end(  WAIT(1.0) );
     poStateFunctions->Set_next_clock_year_start(  TIME_INFINITE );
     poStateFunctions->Set_next_midyear_clock_event(  WAIT(0.5) );

    // Calls
    SetSchoolYearOneClock();
    SetSchoolYearTwoClock();
}

TIME Clock::timeClockYearEndEvent() { return next_clock_year_end ; }
void Clock::ClockYearEndEvent()
{
    // Age all Persons to the current time.
    int nPerson = asAllPerson->Count();
    for ( int nJ = 0; nJ < nPerson; nJ++ ) 
    {
		auto paPerson = asAllPerson->Item( nJ );
        paPerson->UpdateTime()->YearEnd();
    }

    // Call End of year of Calibrator
    lClockToCalibrator ->UpdateTime()->CalibratorYearEnd();

    // Schedule the next year change.
     poStateFunctions->Set_next_clock_year_end(  TIME_INFINITE );
     poStateFunctions->Set_next_clock_year_start(  WAIT(0.0) );
}



} // namespace 
