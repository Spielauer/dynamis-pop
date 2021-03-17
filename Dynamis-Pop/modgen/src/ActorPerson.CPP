#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/ActorPerson.mpp"
//LABEL(ActorPerson, EN) Core functionality of the Person actor

/* NOTE(ActorPerson, EN)
This module introduces the main actor of the model: Person. It contains the core functionality and general states of the actor Person. The essential functions of this module are the Start() and Finish() functions. Start() is called when an actor is created and initialises all states. Finish() terminates the Person at death and cleans up the memory space. 

Persons are created in four different ways, each corresponding to a creation_type:

- Persons from the starting population (CT_START) are created based on a corresponding starting population record. These persons can be both residents or living abroad. They are created on their date of birth, i.e., in the past. If the date of birth in the starting population file is an integer value, a random moment within this year is assigned. 
- Persons from an immigration pool (CT_POOL) are future immigrants created from starting population records marked as belonging to a specific pool, such as e.g. international students or asylum seekers. The same pools are used for sampling new immigrants each year, whereby time-values (like time of birth) are adapted to the respective year of immigration. For example, a person currently 30 years old in the starting population file will also be 30 years old when immigrating in 2050, the year of birth being shifted accordingly. Like all persons, they are created on their date of birth. They enter the actual simulation at the moment of first immigration; their lives between birth and immigration have to be modelled separately. 
- Persons born in the simulation (CT_BIRTH) are created at a birth event during the simulation. The life-courses of these people are entirely synthetic. 
- Immigrants created from scratch (CT_SCRATCH) are immigrants not corresponding to starting population records but created ‘from scratch’ their lives before immigration having to be modelled separately (or characteristics being sampled from some simulated host population). Typically only a few personal characteristics are known from parameters, e.g. their time and place of birth. 

The Start() function initialises the following central states:

- creation_type
- time
- time_of_birth (same as time)
- sex
- ever_resident
- time_of_first_immigration (for immigrants from immigration pools or immigrants created from scratch)
- geo_of_birth
- geo (same as geo_of_birth) 
- geo_prev (same as geo_of_birth) 
- calendar_year
- ethnicity (except for immigrants from scratch who get it assigned at immigration)

Once leaving the Start() function, the time of birth and the current point in time in the simulation are known. This allows flagging an actor to be alive and setting all family links at the right point in time. For that reason, the event SetAliveEvent()is scheduled with zero waiting time at the end of the Start() function. 

For those created from the starting population, the SetAliveEvent() function also creates a temporary actor available via the link lStartPop. This actor stores the starting population information. 

- lStartValues->StartPopValue[nJ]  

The fields of the starting population file record are available until the start of the simulated time. For those created from pools, they are available until first immigration.
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor sets
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor_set Person asAllPerson filter is_alive;   //EN All Person actors

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Types
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//classification SEX                              //EN Sex
//{
//    FEMALE,                                     //EN Female
//    MALE                                        //EN Male
//};

//classification CREATION_TYPE                    //EN Person creation type
//{
//    CT_START,                                   //EN Person from starting population
//    CT_POOL,                                    //EN Person from an immigration pool
//    CT_BIRTH,                                   //EN Person born in simulation
//    CT_SCRATCH                                  //EN Immigrant created from scratch
//};
                                                     
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor states and functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person                                            
//{
//    TIME            time_of_birth = { 0.0 };                //EN Time of birth
//    SEX             sex = { FEMALE };                       //EN Sex
//
//    CREATION_TYPE   creation_type = { CT_START };           //EN Creation type
//
//    GEO             geo_birth = { GEO_00 };                 //EN Place of birth
//    GEO             geo  = { GEO_00 };                      //EN Current place of residence
//    GEO             geo_prev  = { GEO_00 };                 //EN Previous place of residence
//
//    logical         is_alive = { FALSE };                   //EN is alive
//    logical         ever_resident = { FALSE };              //EN Ever a resident
//    logical         is_ready_for_birthtables = { FALSE };   //EN All characteristics set at birth
//
//    logical         is_resident = (geo < SIZE(GEO_NAT));    //EN Is resident
//    
//    logical         ready_to_set_alive = { FALSE };         //EN Ready to set alive
//    event           timeSetAliveEvent, SetAliveEvent;       //EN Set alive
//    
//    Person_ptr      ptr_creator;                            //EN Pointer to creator (family head)
//    Observation_ptr ptr_observation;                        //EN Pointer to observation
//
//    //EN Year of birth
//    ALL_YEAR_RANGE year_of_birth = COERCE(ALL_YEAR_RANGE,int(time_of_birth));    
//
//    //EN Start - the first function called after an actor is created
//    void Start(Observation_ptr peObs, Person *peCreator, int nYearOfImmigration, SEX nImmiSex);
//    
//    //EN Finish - the last function called before an actor is destroyed
//    void Finish();    
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Person::Start(Observation_ptr peObs, Person *peCreator, int nYearOfImmigration, SEX nImmiSex)
{ if (__finished) {StartAfterFinish(); return; }; StartSpecial1();
    // Initialize all attributes (OpenM++).
    initialize_attributes();

    // Setting the actor weight (Modgen only)
    Set_actor_weight(gprParam->ScalingFactor); Set_actor_subsample_weight(gprParam->ScalingFactor);

    TIME dTime;
    ptr_observation = peObs;

    if (peCreator) ptr_creator = peCreator;
    else ptr_creator = NULL;

    // Identify creation type
    if (peObs && peObs-> pop_pool == PP_NON)  poStateFunctions->Set_creation_type(  CT_START );     // starting population
    else if (peObs && peObs-> pop_pool != PP_NON)  poStateFunctions->Set_creation_type(  CT_POOL ); // immigrant from pool
    else if (!peObs && peCreator)  poStateFunctions->Set_creation_type(  CT_BIRTH );               // birth
    else  poStateFunctions->Set_creation_type(  CT_SCRATCH );                                      // immigrant modeled from scratch

    // Persons from start only
    if ( creation_type == CT_START)
    {
        dTime = peObs->pmc[PMC_BIRTH];
        if (int(dTime) == dTime) dTime = dTime + RandUniform(11);

        // cannot be born before creator
        if (ptr_creator && dTime < ptr_creator->time_of_birth ) dTime = ptr_creator-> time_of_birth + RandUniform(15)/10000.0;
         poStateFunctions->Set_ever_resident(  (peObs->pmc[PMC_GEO] < SIZE(GEO_NAT)
            || peObs->pmc[PMC_GEOBIR] < SIZE(GEO_NAT)
            || peObs->pmc[PMC_GEOPRE] < SIZE(GEO_NAT)) );
    }

    // Persons from pools only
    if ( creation_type == CT_POOL)
    {
        // time (now at birth)
        dTime = peObs->pmc[PMC_BIRTH] + 1 + RANGE_POS(SIM_YEAR_RANGE, nYearOfImmigration);
        if (int(dTime) == dTime) dTime = dTime + RandUniform(17);
        
        // time at immigration
        double  dTimeImmi = (ptr_creator) ? TIME_INFINITE : nYearOfImmigration + RandUniform(19);
        if (!ptr_creator && dTimeImmi < dTime) dTimeImmi = dTime + RandUniform(23) / 1000;

        // cannot be born before creator (creator is oldest)
        if (ptr_creator && dTime < ptr_creator->time_of_birth ) dTime = ptr_creator-> time_of_birth + RandUniform(21)/10000.0;
        // cannot be born after creator immigrates
        if (ptr_creator && dTime > ptr_creator->time_of_first_immigration ) dTime = ptr_creator-> time_of_first_immigration - RandUniform(22)/10000.0;
        
         poStateFunctions->Set_time_of_first_immigration(  dTimeImmi );
         poStateFunctions->Set_ever_resident(  FALSE );
    }

    // Persons from starting population including pools -> create start values actor
    if ( creation_type == CT_START ||  creation_type == CT_POOL)
    {
        // Initialise states
               poStateFunctions->Set_sex(  (SEX)(int)peObs->pmc[PMC_SEX] );
         poStateFunctions->Set_geo_birth(  (GEO)(int)peObs->pmc[PMC_GEOBIR] );
         poStateFunctions->Set_ethnicity(  (ETHNICITY)peObs->pmc[PMC_ETHNO] );
    }

    // If immigrant from scratch 
    if ( creation_type == CT_SCRATCH)
    {
         poStateFunctions->Set_sex(  nImmiSex );
         poStateFunctions->Set_time_of_first_immigration(  nYearOfImmigration + RandUniform(7) );
         poStateFunctions->Set_ever_resident(  FALSE );
        SetGeobirthTimeofbirthCtScratch();
        dTime = time_of_birth ;
    }

    // If birth in simulation
    if ( creation_type == CT_BIRTH)
    {
        dTime               = peCreator->UpdateTime() -> time ;
                   poStateFunctions->Set_geo_birth(  peCreator->geo  );
               poStateFunctions->Set_ever_resident(  peCreator->is_resident  );
        if (RandUniform(20) < 100.0 / (100.0 + gprParam->SexRatio[CHKBND(RANGE_POS(SIM_YEAR_RANGE,(int)dTime),151,"SexRatio",0)]))  poStateFunctions->Set_sex(  FEMALE ); 
        else             poStateFunctions->Set_sex(  MALE );
                   poStateFunctions->Set_ethnicity(  GetInheritedEthnicity(peCreator->ethnicity ) );
                 poStateFunctions->Set_educ_mother(  peCreator->educ_one_fate  );
         poStateFunctions->Set_mother_age_at_birth(  peCreator->UpdateTime() -> age  );
                      poStateFunctions->Set_parity(  0 );
    }

    // All
                     poStateFunctions->Set_age(  0, true );
                    poStateFunctions->Set_time(  dTime, true );
                     poStateFunctions->Set_geo(  geo_birth  );
                poStateFunctions->Set_geo_prev(  geo  );
           poStateFunctions->Set_calendar_year(  int(time ) );
           poStateFunctions->Set_time_of_birth(  time  );
        poStateFunctions->Set_geo_want_to_move(  geo  );
      poStateFunctions->Set_ready_to_set_alive(  TRUE );

    // Have the entity enter the simulation (OpenM++).
    enter_simulation();
StartSpecial2();}

void Person::Finish()
{
    // last things still to do
    doMaintainLinksAtDeath();
    CalculateHCIVariables();

    // ready for end of life tables
     poStateFunctions->Set_is_alive(  FALSE );

    if (lStartValues ) lStartValues ->UpdateTime()->Finish();

    // Have the entity exit the simulation (OpenM++).
    exit_simulation();
if (__finished) {	DoubleFinish(case_seed);	return; } else __finished = TRUE;FinishSpecial();}

TIME Person::timeSetAliveEvent() 
{
    if (ready_to_set_alive ) return WAIT(0.0);
    else return TIME_INFINITE;
}


void Person::SetAliveEvent()
{
     Set_lCalibrator(  asCalibrator->Item(0) );
     poStateFunctions->Set_is_alive(  TRUE );
     poStateFunctions->Set_ready_to_set_alive(  FALSE );

    // For those created from the starting population create an actor for temproarily storing starting population values
    // and link the person to other family members. (The oldest actor is created first and available via the pointer ptr_creator.
    if ( creation_type == CT_START ||  creation_type == CT_POOL)
    {
        // Create a temporary actor StartPopValues
        auto prStartpopValues = new StartpopValues(); prStartpopValues->UpdateTime()->Start(time );
         Set_lStartValues(  prStartpopValues );
        for (int nJ = 0; nJ < SIZE(PERSON_MICRODATA_COLUMNS); nJ++)
        {
            lStartValues ->StartPopValue[nJ] = ptr_observation->pmc[nJ];
        }
        lStartValues -> UpdateTime() -> poStateFunctions->Set_is_activated(  TRUE );
        
        // Link actor to other family members
        if (ptr_creator) doLinkToFamilyAtStart();
    }
    if ( creation_type == CT_BIRTH)
    {
        doLinkToFamilyAtStart();
    }
    if ( creation_type == CT_SCRATCH)
    {
        FindImmigrantMother();
    }

    // Calls for setting states at birth
    SetEduc1BaseFate();
    SetEducOneEntryAgeDroputGrade(); 
    DecideImmunizationStatusResidents();
    DecideImmunizationStatusImmigrants();
    DecideStuntingFate();

    // all set
     poStateFunctions->Set_is_ready_for_birthtables(  TRUE );
}


} // namespace 
