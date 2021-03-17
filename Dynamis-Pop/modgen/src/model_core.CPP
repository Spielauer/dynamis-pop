#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/model_core.mpp"
//LABEL(model_core, EN) Core simulation functions

/* NOTE(model_core, EN)
This module contains the core simulation function and model settings. The Simulation() function is creating all actors known at the start of the simulation:

- A single central Clock actor handling peridic clock events affecting all other actors
- A Calibrator actor, used for alignment and calibration routines which require knowledge on current population level outcomes
- Temporary Observation actors corresponding to th erecords of the starting population file
- Person actors representing the simulated population. Persons created already at the start of teh simulation are persons from the starting population and immigrants created from scratch. Persons from the starting population can contain both residents and non-residenst (e.g. temporary work migrants) as well as immigrants created from immigration pools. In the latter case, calendar-time sensitive starting population characteristics (like time of birth) are to be shifted according the year of immigration. 

*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Settings
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//model_type time_based, just_in_time;        // The model type
//time_type double;                           // The data type used to represent time
//languages 
//{
//    EN // English
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simulation function
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Simulation()
{
    extern void LogSimulationStart(); 
    extern void SimulateEvents(); 

    // Write simulation start message to log
    LogSimulationStart();

    ////////////////////////////////////////////////////////////////////
    // Start of the model-specific part of the simulation engine      //
    ////////////////////////////////////////////////////////////////////
    
    // Create the single Clock actor
    auto paClock = new Clock(); paClock->UpdateTime()->Start();

    // Create the single Clock actor
    auto paCalibrator = new Calibrator(); paCalibrator->UpdateTime()->Start();

    // Create observations
    input_csv in_csv; in_csv.open(gprParam->MicroDataInputFile); in_csv.read_header();
    for (long nJ = 0; nJ < gprParam->MicroDataInputFileSize; nJ++)
    {
        in_csv.read_record(nJ);
        auto paObservation = new Observation();
        paObservation->UpdateTime()->Start(in_csv);
    }
    in_csv.close();

    // Identify oldest of each family
    for (long nJ = 0; nJ < gprParam->MicroDataInputFileSize; nJ++)
    {
        // Nobody identified as obs_oldest yet
        if (asObservationByFamOldest[asObservations->Item(nJ)->fam_id ]->Count() == 0)
        {
            Observation_ptr ptrThisObservation = asObservationByFam[asObservations->Item(nJ)->fam_id ]->Item(0);
            Observation_ptr ptrOldestObservation = ptrThisObservation;
            double dEarliestBirthdate = ptrOldestObservation->obs_birth ;
            for (int nI = 0; nI < asObservationByFam[asObservations->Item(nJ)->fam_id ]->Count(); nI++)
            {
                ptrThisObservation = asObservationByFam[asObservations->Item(nJ)->fam_id ]->Item(nI);
                if (ptrThisObservation-> obs_birth <= dEarliestBirthdate)
                {
                    dEarliestBirthdate = ptrThisObservation->obs_birth ;
                    ptrOldestObservation = ptrThisObservation;
                }
            }
            ptrOldestObservation-> UpdateTime() -> poStateFunctions->Set_obs_oldest(  TRUE );
        }
    }

    // Create the population not belonging to an immigration pool
    /////////////////////////////////////////////////////////////

    // Set observation weights of observations not belonging to an immigration pool
    for (long nJ = 0; nJ < asObservationHeads[PP_NON]->Count(); nJ++)
    {
        // Set integer weights obs_weight of heads..
        double  dWeight = asObservationHeads[PP_NON]->Item(nJ)->pmc[PMC_WEIGHT] / gprParam->ScalingFactor;
        int nWeight = (int)dWeight; if (RandUniform(1) < dWeight - nWeight) nWeight++;
        asObservationHeads[PP_NON]->Item(nJ)-> UpdateTime() -> poStateFunctions->Set_obs_weight(  nWeight );
        // ... and now the same for their family members
        for (int nX = 0; nX < asObservationNonHeads[asObservationHeads[PP_NON]->Item(nJ)->fam_id ][PP_NON]->Count(); nX++)
        {
            asObservationNonHeads[asObservationHeads[PP_NON]->Item(nJ)->fam_id ][PP_NON]->Item(nX)-> UpdateTime() -> poStateFunctions->Set_obs_weight(  nWeight );
        }
    }

    // Create the starting population excl. immigrant population pools
    while (asSimulatedObservationHeads[PP_NON]->Count() > 0)
    {
        // First create the head
        auto paObservation = asSimulatedObservationHeads[PP_NON]->Item(0);
        auto paPerson = new Person();
        paPerson->UpdateTime()->Start(paObservation, NULL, 0, MALE);
        // And now all other members of this family
        for (int nJ = 0; nJ < asObservationNonHeads[paObservation->fam_id ][PP_NON]->Count(); nJ++)
        {
            auto paOtherPerson = new Person();
            paOtherPerson->UpdateTime()->Start(asObservationNonHeads[paObservation->fam_id ][PP_NON]->Item(nJ), paPerson, 0, MALE);
        }
        paObservation->UpdateTime() -> poStateFunctions->MinusMinusSet_obs_weight();
    }

    // Delete all observation actors not belonging to an immigration pool
    while (asObservationAll[PP_NON]->Count() > 0) asObservationAll[PP_NON]->Item(0)->UpdateTime()->Finish();

    // Create the potential immigrant population from immigration pools
    ////////////////////////////////////////////////////////////////////

    if (gprParam->ModelImmigrationFromPools)
    {
        for (int nPool = 1; nPool < SIZE(POP_POOL); nPool++)
        {
    		for (int nYear = MIN(SIM_YEAR_RANGE); nYear < int(SIMULATION_END()); nYear++)
            {
                double dTargetPop = (gprParam->ImmiPoolSize[CHKBND(nPool - 1,3,"ImmiPoolSize",0)][CHKBND(RANGE_POS(SIM_YEAR_RANGE, nYear),151,"ImmiPoolSize",1)]) / gprParam->ScalingFactor;
                double dSumWeights = 0.0;
                for (long nJ = 0; nJ < asObservationHeads[nPool]->Count(); nJ++)
                {
                    dSumWeights = dSumWeights + asObservationHeads[nPool]->Item(nJ)->pmc[PMC_WEIGHT] *
                        (1 + asObservationNonHeads[asObservationHeads[nPool]->Item(nJ)->fam_id ][nPool]->Count());
                }
                for (long nJ = 0; nJ < asObservationHeads[nPool]->Count(); nJ++)
                {
                    double dWeight = asObservationHeads[nPool]->Item(nJ)->pmc[PMC_WEIGHT] * dTargetPop / dSumWeights;
                    int nWeight = int(dWeight);
                    if (RandUniform(6) < (dWeight - nWeight))
                    {
                        nWeight++;
                    };
                    asObservationHeads[nPool]->Item(nJ)-> UpdateTime() -> poStateFunctions->Set_obs_weight(  nWeight );
                }
                while (asSimulatedObservationHeads[nPool]->Count() > 0)
                {
                    // First create the head
                    auto paObservation = asSimulatedObservationHeads[nPool]->Item(0);
                    auto paPerson = new Person();
                    paPerson->UpdateTime()->Start(paObservation, NULL, nYear, MALE);

                    // And now all other members of this family
                    for (int nJ = 0; nJ < asObservationNonHeads[paObservation->fam_id ][nPool]->Count(); nJ++)
                    {
                        auto paOtherPerson = new Person();
                        paOtherPerson->UpdateTime()->Start(asObservationNonHeads[paObservation->fam_id ][nPool]->Item(nJ), 
                            paPerson, nYear, MALE);
                    }
                    paObservation->UpdateTime() -> poStateFunctions->MinusMinusSet_obs_weight();
                }
            }
        }
    }

    // Delete all observation actors left
    while (asObservations->Count() > 0) asObservations->Item(0)->UpdateTime()->Finish();


    // Create Immigrants from scratch
    ////////////////////////////////////////////////////////////////////

	if (gprParam->ModelImmigrationFromScratch)
	{
        for (int nYear = MIN(SIM_YEAR_RANGE); nYear < int(SIMULATION_END()); nYear++)
        {
            for (int nSex = 0; nSex < SIZE(SEX); nSex++)
            {
                int nCohortSize = int(gprParam->NumberImmigrantsFromScratch[CHKBND(RANGE_POS(SIM_YEAR_RANGE, nYear),151,"NumberImmigrantsFromScratch",0)][CHKBND((SEX)nSex,2,"NumberImmigrantsFromScratch",1)]
                    / gprParam->ScalingFactor);
                double dCohortSize = gprParam->NumberImmigrantsFromScratch[CHKBND(RANGE_POS(SIM_YEAR_RANGE, nYear),151,"NumberImmigrantsFromScratch",0)][CHKBND((SEX)nSex,2,"NumberImmigrantsFromScratch",1)]
                    / gprParam->ScalingFactor;
                if (RandUniform(4) < dCohortSize - nCohortSize) nCohortSize++;

                for (int nIndex = 0; nIndex < nCohortSize; nIndex++)
                {
                    auto prPerson = new Person();
                    prPerson->UpdateTime()->Start(NULL, NULL, nYear, (SEX)nSex);
                }
            }
        }
	}

    ////////////////////////////////////////////////////////////////////
    // End of the model-specific part of the simulation engine        //
    ////////////////////////////////////////////////////////////////////

    // Simulate events until there are no more.
    SimulateEvents();
}

} // namespace 
