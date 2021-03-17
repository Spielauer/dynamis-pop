#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/MigrationGeneral.mpp"
//LABEL(MigrationGeneral, EN) Migration General

/* NOTE(MigrationGeneral, EN)
This module implements general migration stuff
*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor sets
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor_set Person asWantToMove[geo][geo_want_to_move][sex][age_mig] filter is_alive;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Types
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//partition   AGE_MIG_PART { 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };	    //EN Age Groups
//range       AGE_MIG { 0, 12 };                                                  //EN Age index

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters
//{
//    logical MigrationTryKeepingFamiliesTogether;
//};

//parameter_group PG05_Migration                                                  //EN Migration
//{
//    ModelMigration,
//    MigrationTryKeepingFamiliesTogether,
//    MigrationProbability,
//    MigrationDestination
//};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor states and functions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    AGE_MIG age_mig = COERCE(AGE_MIG,self_scheduling_split(age,AGE_MIG_PART));  //EN Age index
//
//    TIME time_first_move  = { TIME_INFINITE };                                  //EN Time of first residential move
//    TIME time_last_move  = { TIME_INFINITE };                                   //EN Time of last residential move
//
//    logical first_move_done = { FALSE };
//
//    TIME time_of_first_immigration = { TIME_INFINITE };                         //EN Time of first immigration
//    void doResidentialMove(GEO cDestination);                                   //EN Immigration
//
//    GEO geo_want_to_move = { GEO_00 };                                          //EN Place a person wants to move
//
//    GEO_NAT     geo_nat = aggregate( geo, GEO_NAT );                            //EN Geo of residence 
//
//    REGION_INT  region_int = aggregate( geo, REGION_INT);                       //EN Region of residence
//    REGION_NAT  region_nat = aggregate( region_int, REGION_NAT);	            //EN Region of residence
//    REGION_INT  region_birth = aggregate( geo_birth, REGION_INT);               //EN Region of birth
//
//};

//actor StartpopValues
//{
//    // Setting clock for (past) residential moves
//    logical is_set_residential_move_clock = { FALSE };                          //EN Clock is set
//    event   timeSetResidentialMoveClock, SetResidentialMoveClock;               //EN 
//
//    // Clock events for (past) residential moves 
//    TIME    scheduled_time_first_move = { TIME_INFINITE };                      //EN Time of first move
//    TIME    scheduled_time_last_move = { TIME_INFINITE };                       //EN Time of last move
//    event   timeFirstMoveEvent, FirstMoveEvent;                                 //EN First residential move event
//    event   timeLastMoveEvent, LastMoveEvent;                                   //EN Last residential move event
//};

//actor Calibrator
//{
//    integer migration_counter = {0};                                            //EN Number of internal migration moves
//    integer immigration_counter = {0};                                          //EN Number of immigrants
//    integer emigration_counter = {0};                                           //EN Number of emigrants
//};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Person::doResidentialMove(GEO cDestination)
{
    // Update geography
     poStateFunctions->Set_geo_prev(  geo  );
     poStateFunctions->Set_geo(  cDestination );

    // Update timestamps
     poStateFunctions->Set_time_last_move(  time  );
    if (  time_first_move == TIME_INFINITE )  time_first_move =   time  ;

    // Update ever resident status
    if (! ever_resident &&  geo < SIZE(GEO_NAT))  poStateFunctions->Set_ever_resident(  TRUE );    

    // Increase counters Calibrator level
    if ( geo_prev < SIZE(GEO_NAT) &&  geo < SIZE(GEO_NAT)) lCalibrator ->UpdateTime() -> migration_counter++;
    else if ( geo_prev >= SIZE(GEO_NAT) &&  geo < SIZE(GEO_NAT)) lCalibrator ->UpdateTime() -> immigration_counter++;
    else if ( geo_prev < SIZE(GEO_NAT) &&  geo >= SIZE(GEO_NAT)) lCalibrator ->UpdateTime() -> emigration_counter++;

    // Tell family to follow
    if (gprParam->MigrationTryKeepingFamiliesTogether)
    {
        if (lSpouse ) lSpouse -> UpdateTime() -> poStateFunctions->Set_geo_want_to_move(  geo  );
        if (lHHFather ) lHHFather -> UpdateTime() -> poStateFunctions->Set_geo_want_to_move(  geo  );
        if (lHHMother ) lHHMother -> UpdateTime() -> poStateFunctions->Set_geo_want_to_move(  geo  );
        if ( children_in_household > 0)
        {
            int nIndex;
            auto prChild = ( sex == FEMALE) ? mlHHMotherChildren ->GetNext(0, &nIndex) : mlHHFatherChildren ->GetNext(0, &nIndex);
            while (prChild)
            {
                prChild-> UpdateTime() -> poStateFunctions->Set_geo_want_to_move(  geo  );
                // set pointer to next child
                prChild = ( sex == FEMALE) ? mlHHMotherChildren ->GetNext(nIndex + 1, &nIndex) : mlHHFatherChildren ->GetNext(nIndex + 1, &nIndex);
            }
        }
    }
}

TIME StartpopValues::timeSetResidentialMoveClock()
{
    if ( is_activated && !is_set_residential_move_clock ) return WAIT(0);
    else return TIME_INFINITE;
}

void StartpopValues::SetResidentialMoveClock()
{
    // set scheduled_time_first_move if geo at birth different from previous geo
    if (lPersonStartpop -> creation_type == CT_START && StartPopValue[CHKBND(PMC_GEOBIR,18,"StartpopValues.StartPopValue",0)] != StartPopValue[CHKBND(PMC_GEOPRE,18,"StartpopValues.StartPopValue",0)])
    {
        TIME dTimeFirst = StartPopValue[CHKBND(PMC_MOVEFIRST,18,"StartpopValues.StartPopValue",0)];
        if (dTimeFirst > lPersonStartpop -> time_of_birth && dTimeFirst < MIN(SIM_YEAR_RANGE))  poStateFunctions->Set_scheduled_time_first_move(  dTimeFirst );
        else  poStateFunctions->Set_scheduled_time_first_move(  (lPersonStartpop -> time_of_birth + (int)lPersonStartpop -> time_of_birth + 1) / 2.0 ); 
    }
    else  poStateFunctions->Set_scheduled_time_first_move(  TIME_INFINITE );

    // set scheduled_time_last_move if previous geo different from geo
    if (lPersonStartpop -> creation_type == CT_START && StartPopValue[CHKBND(PMC_GEO,18,"StartpopValues.StartPopValue",0)] != StartPopValue[CHKBND(PMC_GEOPRE,18,"StartpopValues.StartPopValue",0)])
    {
        TIME dTimeLast = StartPopValue[CHKBND(PMC_MOVELAST,18,"StartpopValues.StartPopValue",0)];
        if (dTimeLast > lPersonStartpop -> time_of_birth && dTimeLast < MIN(SIM_YEAR_RANGE) && ( scheduled_time_first_move == TIME_INFINITE ||  scheduled_time_first_move < dTimeLast))  poStateFunctions->Set_scheduled_time_last_move(  dTimeLast );
        else
        {
            if ( scheduled_time_first_move == TIME_INFINITE)  poStateFunctions->Set_scheduled_time_last_move(  (lPersonStartpop -> time_of_birth + MIN(SIM_YEAR_RANGE)) / 2.0 );
            else  poStateFunctions->Set_scheduled_time_last_move(  ( scheduled_time_first_move + MIN(SIM_YEAR_RANGE)) / 2.0 ); 
        }
    }
    else  poStateFunctions->Set_scheduled_time_last_move(  TIME_INFINITE );

    // do not call again
     poStateFunctions->Set_is_set_residential_move_clock(  TRUE );
}

TIME StartpopValues::timeFirstMoveEvent() 
{
    //return TIME_INFINITE;
    return scheduled_time_first_move ; 
}

void StartpopValues::FirstMoveEvent()
{
    lPersonStartpop ->UpdateTime()->doResidentialMove((GEO)(int)StartPopValue[CHKBND(PMC_GEOPRE,18,"StartpopValues.StartPopValue",0)]);
     poStateFunctions->Set_scheduled_time_first_move(  TIME_INFINITE );
}


TIME StartpopValues::timeLastMoveEvent() 
{
    //return TIME_INFINITE;
    return scheduled_time_last_move ;
}
void StartpopValues::LastMoveEvent() 
{ 
    lPersonStartpop ->UpdateTime()->doResidentialMove((GEO)(int)StartPopValue[CHKBND(PMC_GEO,18,"StartpopValues.StartPopValue",0)]);
     poStateFunctions->Set_scheduled_time_last_move(  TIME_INFINITE );
}


} // namespace 
