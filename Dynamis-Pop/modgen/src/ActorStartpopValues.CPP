#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/ActorStartpopValues.mpp"
//LABEL (ActorStartpopValues, EN) Definition and initialization of the StartpopValues actor

/*NOTE(ActorStartpopValues, EN)
This module introduces a StartpopValues actor. The StartpopValues actor is temporary and contains information from the starting population record of persons created from the starting population file. It is deleted once the simulated time is entered or � for tose from population pools - once the country is entered. 

The actor can also contain events and functions only relevant for those created from the starting population file and only used before entering the simulation time (or until first immigration). For example, it is used by the ImmigrationGeneral module for scheduling residential moves which occurred before entering the simulated time. 
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Links
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//EN Link between Person and StartpopValues
//link StartpopValues.lPersonStartpop Person.lStartValues;	

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor states and functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor StartpopValues //EN Actor storing starting population values
//{
//    logical     is_activated = { FALSE };                       //EN Is activated (synchronized when Person set alive) 
//    double      StartPopValue[PERSON_MICRODATA_COLUMNS];        //EN Startpop values
//
//    void        Start(TIME dTime);                              //EN Starts the Clock
//    void        Finish();                                       //EN Finishes the Clock
//
//    logical     is_removeable = (lPersonStartpop->is_mortal);   //EN Can be removed 
//    event       timeSuicideEvent, SuicideEvent;                 //EN Actor self-destruction when not required anymore
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void StartpopValues::Start(TIME dTime)
{ if (__finished) {StartAfterFinish(); return; }; StartSpecial1();
    // Initialize all attributes (OpenM++).
    initialize_attributes();

     poStateFunctions->Set_time(  dTime, true );
     poStateFunctions->Set_age(  0, true );

    // Have the entity enter the simulation (OpenM++).
    enter_simulation();
StartSpecial2();}
 

void StartpopValues::Finish()
{
    // Have the entity exit the simulation (OpenM++).
    exit_simulation();
if (__finished) {	DoubleFinish(case_seed);	return; } else __finished = TRUE;FinishSpecial();}

TIME StartpopValues::timeSuicideEvent() { if (is_removeable ) return WAIT(0.0); else return TIME_INFINITE; }
void StartpopValues::SuicideEvent() { Finish(); }


} // namespace 
