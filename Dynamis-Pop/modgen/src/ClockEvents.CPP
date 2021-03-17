#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/ClockEvents.mpp"
//LABEL(ClockEvents, EN) Clock Events

/* NOTE(ClockEvents, EN)
This module contains periodic clock events like the end of a calendar year, the start of a new calendar year, and birthdays. 
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor states and functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    AGE_RANGE       integer_age = { 0 };                        //EN Age
//    TIME            time_next_birthday = { TIME_INFINITE };     //EN Time of next birthday
//    event           timeBirthdayEvent, BirthdayEvent;           //EN Birthday Event
//
//    void            YearEnd();                                  //EN Year End Function
//    void            YearStart();                                //EN Year Start Function
//    void            MidYear();                                  //EN Mid year
//
//    ALL_YEAR_RANGE  calendar_year;                              //EN Calendar year
//};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simulation function
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Person::MidYear()
{
    // Code to be performed at each midyear
    MidYearIMPLEMENT_HOOK();
}

void Person::YearEnd()
{
    // Code to be performed at each end of the year can be placed here
    AdjustEducOne();
    
    YearEndIMPLEMENT_HOOK();
}

void Person::YearStart()
{
    // Code to be performed at each start of the year can be placed here
    YearStartIMPLEMENT_HOOK();
}


TIME Person::timeBirthdayEvent()
{
    if (  integer_age == 0 ) return  time_of_birth + 1.0;
    else return time_next_birthday ;
}

void Person::BirthdayEvent()
{
    // Increment integer age
    if ( integer_age < MAX(AGE_RANGE)) poStateFunctions->PlusPlusSet_integer_age();

    // Code to be performed at each birthday can be entered here
    BirthdayEventIMPLEMENT_HOOK();

    // Set clock for next birthday
     poStateFunctions->Set_time_next_birthday(  WAIT(1) );
}

} // namespace 
