/**
 * @file om_developer.cpp
 * Developer-supplied C++ functions
 */

#include "omc/omPch.h"
#include "omc/omSimulation.h"
using namespace openm;

#line 1 "../code/ActorCalibrator.mpp"





















	
              




















void Calibrator::CalibratorYearEnd()
{
    
    SetEduc1AdjustmentFactors();
}

void Calibrator::Start( )
{
    
    initialize_attributes();

    
    Set_actor_weight(ScalingFactor); Set_actor_subsample_weight(ScalingFactor);

    lCalibratorToClock = asClock->Item(0);
    time = lCalibratorToClock->time;

    
    enter_simulation();
}

void Calibrator::Finish()
{
    
    exit_simulation();
}

#line 1 "../code/ActorClock.mpp"











































void Clock::Start( )
{
    
    initialize_attributes();

    
    Set_actor_weight(ScalingFactor); Set_actor_subsample_weight(ScalingFactor);

    time = (TIME) MIN(ALL_YEAR_RANGE);
    clock_year = MIN(ALL_YEAR_RANGE);
    age = 0;
    next_clock_year_end = WAIT(1);
    next_clock_year_start = TIME_INFINITE;

    
    enter_simulation();
}

void Clock::Finish()
{
    
    exit_simulation();
}





TIME Clock::timeClockMidyearEvent() { return next_midyear_clock_event; }
void Clock::ClockMidyearEvent()
{
    
    UpdatePartnershipStatus();

    
    int nPerson = asAllPerson->Count();
    for ( int nJ = 0; nJ < nPerson; nJ++ ) 
    {
		auto paPerson = asAllPerson->Item( nJ );
        paPerson->MidYear();
    }
    
    next_midyear_clock_event = TIME_INFINITE;
}


TIME Clock::timeClockYearStartEvent() { return next_clock_year_start; }
void Clock::ClockYearStartEvent()
{
    clock_year++;

    
    int nPerson = asAllPerson->Count();
    for ( int nJ = 0; nJ < nPerson; nJ++ ) 
    {
		auto paPerson = asAllPerson->Item( nJ );
        paPerson->calendar_year = (ALL_YEAR_RANGE)(int)clock_year;
        paPerson->YearStart();
    }

    
    next_clock_year_end = WAIT(1.0);
    next_clock_year_start = TIME_INFINITE;
    next_midyear_clock_event = WAIT(0.5);

    
    SetSchoolYearOneClock();
    SetSchoolYearTwoClock();
}

TIME Clock::timeClockYearEndEvent() { return next_clock_year_end; }
void Clock::ClockYearEndEvent()
{
    
    int nPerson = asAllPerson->Count();
    for ( int nJ = 0; nJ < nPerson; nJ++ ) 
    {
		auto paPerson = asAllPerson->Item( nJ );
        paPerson->YearEnd();
    }

    
    lClockToCalibrator->CalibratorYearEnd();

    
    next_clock_year_end = TIME_INFINITE;
    next_clock_year_start = WAIT(0.0);
}



#line 1 "../code/ActorObservation.mpp"

















































                         










































void Observation::Start(const input_csv& in_csv)
{
    
    initialize_attributes();

    for (int nJ = 0; nJ < SIZE(PERSON_MICRODATA_COLUMNS); nJ++)
    {
        pmc[nJ] = in_csv[nJ];
    }
    fam_id      = (int)pmc[PMC_FAMID];
    pop_pool    = (POP_POOL)(int)pmc[PMC_POOL];
    obs_birth   = pmc[PMC_BIRTH];

    time = MIN(ALL_YEAR_RANGE); 

    
    enter_simulation();
};

void Observation::Finish()
{
    
    exit_simulation();
};





void om_PreSimulation_0()
{
    
    input_csv inCsv;
    inCsv.open(MicroDataInputFile);

    
    long    lRecordCount = 0;
    double  dResidentPopSize = 0.0;
    while (inCsv.read_record(lRecordCount))
    {
        lRecordCount++;
        
        if ((int)inCsv[PMC_POOL]==PP_NON && (int)inCsv[PMC_GEO]< SIZE(GEO_NAT))
        {
            dResidentPopSize = dResidentPopSize + inCsv[PMC_WEIGHT];
        }
    }
    ScalingFactor = dResidentPopSize / StartPopSampleSize;
    MicroDataInputFileSize = lRecordCount - 1;

    
    inCsv.close();
}


#line 1 "../code/ActorPerson.mpp"







































   


















                                                     









































void Person::Start(Observation_ptr peObs, Person *peCreator, int nYearOfImmigration, SEX nImmiSex)
{
    
    initialize_attributes();

    
    Set_actor_weight(ScalingFactor); Set_actor_subsample_weight(ScalingFactor);

    TIME dTime;
    ptr_observation = peObs;

    if (peCreator) ptr_creator = peCreator;
    else ptr_creator = NULL;

    
    if (peObs && peObs->pop_pool == PP_NON) creation_type = CT_START;     
    else if (peObs && peObs->pop_pool != PP_NON) creation_type = CT_POOL; 
    else if (!peObs && peCreator) creation_type = CT_BIRTH;               
    else creation_type = CT_SCRATCH;                                      

    
    if (creation_type == CT_START)
    {
        dTime = peObs->pmc[PMC_BIRTH];
        if (int(dTime) == dTime) dTime = dTime + RandUniform(11);

        
        if (ptr_creator && dTime < ptr_creator->time_of_birth) dTime = ptr_creator->time_of_birth + RandUniform(15)/10000.0;
        ever_resident = (peObs->pmc[PMC_GEO] < SIZE(GEO_NAT)
            || peObs->pmc[PMC_GEOBIR] < SIZE(GEO_NAT)
            || peObs->pmc[PMC_GEOPRE] < SIZE(GEO_NAT));
    }

    
    if (creation_type == CT_POOL)
    {
        
        dTime = peObs->pmc[PMC_BIRTH] + 1 + RANGE_POS(SIM_YEAR_RANGE, nYearOfImmigration);
        if (int(dTime) == dTime) dTime = dTime + RandUniform(17);
        
        
        double  dTimeImmi = (ptr_creator) ? TIME_INFINITE : nYearOfImmigration + RandUniform(19);
        if (!ptr_creator && dTimeImmi < dTime) dTimeImmi = dTime + RandUniform(23) / 1000;

        
        if (ptr_creator && dTime < ptr_creator->time_of_birth) dTime = ptr_creator->time_of_birth + RandUniform(21)/10000.0;
        
        if (ptr_creator && dTime > ptr_creator->time_of_first_immigration) dTime = ptr_creator->time_of_first_immigration - RandUniform(22)/10000.0;
        
        time_of_first_immigration = dTimeImmi;
        ever_resident = FALSE;
    }

    
    if (creation_type == CT_START || creation_type == CT_POOL)
    {
        
        sex       = (SEX)(int)peObs->pmc[PMC_SEX];
        geo_birth = (GEO)(int)peObs->pmc[PMC_GEOBIR];
        ethnicity = (ETHNICITY)peObs->pmc[PMC_ETHNO];
    }

    
    if (creation_type == CT_SCRATCH)
    {
        sex = nImmiSex;
        time_of_first_immigration = nYearOfImmigration + RandUniform(7);
        ever_resident = FALSE;
        SetGeobirthTimeofbirthCtScratch();
        dTime = time_of_birth;
    }

    
    if (creation_type == CT_BIRTH)
    {
        dTime               = peCreator->time;
        geo_birth           = peCreator->geo;
        ever_resident       = peCreator->is_resident;
        if (RandUniform(20) < 100.0 / (100.0 + SexRatio[RANGE_POS(SIM_YEAR_RANGE,(int)dTime)])) sex = FEMALE; 
        else sex            = MALE;
        ethnicity           = GetInheritedEthnicity(peCreator->ethnicity);
        educ_mother         = peCreator->educ_one_fate;
        mother_age_at_birth = peCreator->age;
        parity              = 0;
    }

    
    age                 = 0;
    time                = dTime;
    geo                 = geo_birth;
    geo_prev            = geo;
    calendar_year       = int(time);
    time_of_birth       = time;
    geo_want_to_move    = geo;
    ready_to_set_alive  = TRUE;

    
    enter_simulation();
}

void Person::Finish()
{
    
    doMaintainLinksAtDeath();
    CalculateHCIVariables();

    
    is_alive = FALSE;

    if (lStartValues) lStartValues->Finish();

    
    exit_simulation();
}

TIME Person::timeSetAliveEvent() 
{
    if (ready_to_set_alive) return WAIT(0.0);
    else return TIME_INFINITE;
}


void Person::SetAliveEvent()
{
    lCalibrator = asCalibrator->Item(0);
    is_alive = TRUE;
    ready_to_set_alive = FALSE;

    
    
    if (creation_type == CT_START || creation_type == CT_POOL)
    {
        
        auto prStartpopValues = new StartpopValues(); prStartpopValues->Start(time);
        lStartValues = prStartpopValues;
        for (int nJ = 0; nJ < SIZE(PERSON_MICRODATA_COLUMNS); nJ++)
        {
            lStartValues->StartPopValue[nJ] = ptr_observation->pmc[nJ];
        }
        lStartValues->is_activated = TRUE;
        
        
        if (ptr_creator) doLinkToFamilyAtStart();
    }
    if (creation_type == CT_BIRTH)
    {
        doLinkToFamilyAtStart();
    }
    if (creation_type == CT_SCRATCH)
    {
        FindImmigrantMother();
    }

    
    SetEduc1BaseFate();
    SetEducOneEntryAgeDroputGrade(); 
    DecideImmunizationStatusResidents();
    DecideImmunizationStatusImmigrants();
    DecideStuntingFate();

    
    is_ready_for_birthtables = TRUE;
}


#line 1 "../code/ActorStartpopValues.mpp"













	





















void StartpopValues::Start(TIME dTime)
{
    
    initialize_attributes();

    time = dTime;
    age = 0;

    
    enter_simulation();
}
 

void StartpopValues::Finish()
{
    
    exit_simulation();
}

TIME StartpopValues::timeSuicideEvent() { if (is_removeable) return WAIT(0.0); else return TIME_INFINITE; }
void StartpopValues::SuicideEvent() { Finish(); }


#line 1 "../code/ChildVaccination.mpp"






























































void Person::DecideImmunizationStatusResidents()
{
    if (in_projected_time && is_resident)
    {
        int nIndex = SPLIT(time - mother_age_at_birth, IMMU_YOB_PART); 
        int nYear = RANGE_POS(SIM_YEAR_RANGE, calendar_year);

        
        double dOddsCare = PreNatalCareOdds[nYear][PP_CONSTANT];                                      
        if (educ_mother == EOL_MEDIUM) dOddsCare = dOddsCare * PreNatalCareOdds[nYear][PP_EDUCMO_1];  
        if (educ_mother == EOL_HIGH)   dOddsCare = dOddsCare * PreNatalCareOdds[nYear][PP_EDUCMO_2];  
        if (region_nat > REGN_00)      dOddsCare = dOddsCare * PreNatalCareOdds[nYear][region_nat+2]; 
        if (nIndex > 0)                dOddsCare = dOddsCare
            * PreNatalCareOdds[nYear][nIndex + SIZE(REGION_NAT) + 1];
        if (mother_age_at_birth < 18)  dOddsCare = dOddsCare * PreNatalCareOdds[nYear][PP_YOUNGMO];

        double dProbCare = dOddsCare / (dOddsCare + 1);
        if (RandUniform(79) < dProbCare) got_prenat_care = GPC_YES;

        
        double dOdds = ChildVaccinationOdds[got_prenat_care][nYear][IP_CONSTANT];                                  
        if (sex==MALE)                 dOdds = dOdds * ChildVaccinationOdds[got_prenat_care][nYear][IP_MALE];      
        if (educ_mother == EOL_MEDIUM) dOdds = dOdds * ChildVaccinationOdds[got_prenat_care][nYear][IP_EDUCMO_1];  
        if (educ_mother == EOL_HIGH)   dOdds = dOdds * ChildVaccinationOdds[got_prenat_care][nYear][IP_EDUCMO_2];  
        if (region_nat > REGN_00)      dOdds = dOdds * ChildVaccinationOdds[got_prenat_care][nYear][region_nat+3]; 
        if (ethnicity_short > ES_00)   dOdds = dOdds
            * ChildVaccinationOdds[got_prenat_care][nYear][ethnicity_short + SIZE(REGION_NAT) + 2];
        if (nIndex > 0)                dOdds = dOdds
            * ChildVaccinationOdds[got_prenat_care][nYear][nIndex + SIZE(REGION_NAT) + SIZE(ETHNICITY_SHORT) + 1];
        if (mother_age_at_birth < 18)  dOdds = dOdds * ChildVaccinationOdds[got_prenat_care][nYear][IP_YOUNGMO];

        double dProb = dOdds / (dOdds + 1);
        if (RandUniform(80) < dProb) is_immunized = TRUE;
    }
}

void Person::DecideImmunizationStatusImmigrants()
{
    if (in_projected_time && !is_resident)
    {
        if (asResidentBabies->Count() > 0)
        {
            auto prDonor = asResidentBabies->GetRandom(RandUniform(77));
            got_prenat_care = prDonor->got_prenat_care;
            is_immunized = prDonor->is_immunized;
        }
        else if (RandUniform(78) < 0.5) is_immunized = TRUE; 
    }
}


#line 1 "../code/ClockEvents.mpp"




























void Person::MidYear()
{
    
    IMPLEMENT_HOOK();
}

void Person::YearEnd()
{
    
    AdjustEducOne();
    
    IMPLEMENT_HOOK();
}

void Person::YearStart()
{
    
    IMPLEMENT_HOOK();
}


TIME Person::timeBirthdayEvent()
{
    if ( integer_age == 0 ) return time_of_birth + 1.0;
    else return time_next_birthday;
}

void Person::BirthdayEvent()
{
    
    if (integer_age < MAX(AGE_RANGE)) integer_age++;

    
    IMPLEMENT_HOOK();

    
    time_next_birthday = WAIT(1);
}

#line 1 "../code/EducationPreSchool.mpp"














































TIME Person::timeSetPreschoolYearsEvent()
{
    if (!preschool_is_decided && creation_type == CT_BIRTH && is_resident && educ_one_grade_attended == 1) return WAIT(0);
    else return TIME_INFINITE;
}

void Person::SetPreschoolYearsEvent()
{   
    if (educ_one_fate != EOL_LOW)
    {
        if (RandUniform(10) < PreSchoolAttendance[sex][region_nat][RANGE_POS(SIM_YEAR_RANGE, calendar_year)][PLP_ANY])
        {
            years_preschool++;
            if (RandUniform(59) < PreSchoolAttendance[sex][region_nat][RANGE_POS(SIM_YEAR_RANGE, calendar_year)][PLS_TWO]) years_preschool++;
        }
    }
    preschool_is_decided = TRUE;
}

#line 1 "../code/EducationPrimaryBase.mpp"









































































void Person::SetEduc1BaseFate()
{
    
    if ( creation_type == CT_BIRTH )
    {
        EDUC_ONE_LEVEL eolFate = EOL_LOW;
        if ( RandUniform(41) < EducTrans1[sex][RANGE_POS(YOB_EDUC_TRANS1,year_of_birth)][geo_birth] ) eolFate = EOL_MEDIUM;
        if ( eolFate == EOL_MEDIUM 
            && RandUniform(42) < EducTrans2[sex][RANGE_POS(YOB_EDUC_TRANS2,year_of_birth)][geo_birth] )  eolFate = EOL_HIGH;
        educ_one_fate = eolFate;
    }

    
    else if ( creation_type == CT_START || creation_type == CT_POOL )
    {
        
        EDUC_ONE_LEVEL eolStart = EOL_LOW;
        if ( int(lStartValues->StartPopValue[PMC_EDUC]) == 1 ) eolStart = EOL_MEDIUM;
        else if ( int(lStartValues->StartPopValue[PMC_EDUC]) == 2 ) eolStart = EOL_HIGH;

        
        int nYOB = int(lStartValues->StartPopValue[PMC_BIRTH]);

        
        if ( nYOB < MIN(YOB_EDUC_TRANS2)) educ_one_fate = eolStart;
        
        
        else if ( nYOB < MIN(YOB_EDUC_TRANS1))
        {
            EDUC_ONE_LEVEL eolFate = EOL_LOW;
            if ( eolStart != EOL_LOW ) eolFate = EOL_MEDIUM;
            if ( eolFate == EOL_MEDIUM 
                && RandUniform(43) < EducTrans2[sex][RANGE_POS(YOB_EDUC_TRANS2,nYOB)][geo_birth] )  eolFate = EOL_HIGH;
            educ_one_fate = eolFate;
        }
        
        else
        {
            EDUC_ONE_LEVEL eolFate = EOL_LOW;
            if ( RandUniform(44) < EducTrans1[sex][RANGE_POS(YOB_EDUC_TRANS1,nYOB)][geo_birth] ) eolFate = EOL_MEDIUM;
            if ( eolFate == EOL_MEDIUM 
                && RandUniform(45) < EducTrans2[sex][RANGE_POS(YOB_EDUC_TRANS2,nYOB)][geo_birth] )  eolFate = EOL_HIGH;
            educ_one_fate = eolFate;
        }
    }

    
    else if (creation_type==CT_SCRATCH)
    {
        
        if (year_of_birth >= MIN(YOB_EDUC_TRANS1) || (asResidentsAge0SexGeo[sex][geo_birth]->Count() == 0 && asResidentsAge0Sex[sex]->Count() == 0))
        {
            EDUC_ONE_LEVEL eolFate = EOL_LOW;
            if (RandUniform(46) < EducTrans1[sex][RANGE_POS(YOB_EDUC_TRANS1, year_of_birth)][geo_birth]) eolFate = EOL_MEDIUM;
            if (eolFate == EOL_MEDIUM && RandUniform(47) < EducTrans2[sex][RANGE_POS(YOB_EDUC_TRANS2, year_of_birth)][geo_birth])  eolFate = EOL_HIGH;
            educ_one_fate = eolFate;
        }
        
        else if (asResidentsAge0SexGeo[sex][geo_birth]->Count() > 0) 
        {
            educ_one_fate = asResidentsAge0SexGeo[sex][geo_birth]->GetRandom(RandUniform(48))->educ_one_fate;
        }
        
        else educ_one_fate = asResidentsAge0Sex[sex]->GetRandom(RandUniform(49))->educ_one_fate; 
    }
}

#line 1 "../code/EducationPrimaryPlanning.mpp"















































#line 1 "../code/EducationPrimaryTracking.mpp"

































          
        
      





























































































void Person::SetEducOneEntryAgeDroputGrade()
{
    int nAge;
    Lookup_EducOneEntryAge(RandUniform(62), (int)educ_one_geo, (int)educ_one_group,
        RANGE_POS(YOB_EDUC_TRANS2,year_of_birth), &nAge);
    educ_one_entry_age = MIN(EDUC_ONE_ENTRY_AGE) + nAge;

    int nGrade;
    Lookup_EducOneDropoutGrade(RandUniform(63), (int)educ_one_geo, (int)educ_one_group,
        RANGE_POS(YOB_EDUC_TRANS2,year_of_birth), &nGrade);
    educ_one_grade_fate = MIN(EDUC_ONE_GRADE) + nGrade;
}





TIME Clock::timeStartSchoolOneYearEvent(){return time_start_school_one_year;}
void Clock::StartSchoolOneYearEvent()
{
    for (long nJ = 0; nJ < asPotentialSchoolOneStudent->Count(); nJ++)
    {
        auto prPerson = asPotentialSchoolOneStudent->Item(nJ);
        if (prPerson->educ_one_status == EOS_WAIT )
        {
            prPerson->educ_one_status = EOS_ATTEND ;
            prPerson->educ_one_grade_attended = prPerson->educ_one_grade_passed + 1;
        }
    }
    time_start_school_one_year = TIME_INFINITE;
}

TIME Clock::timeEndSchoolOneYearEvent(){return time_end_school_one_year;}
void Clock::EndSchoolOneYearEvent()
{
    
    for (long nJ = 0; nJ < asPotentialSchoolOneStudent->Count(); nJ++)
    {
        auto prPerson = asPotentialSchoolOneStudent->Item(nJ);
        prPerson->educ_one_to_process = TRUE;
    }

    
    for (long nJ = 0; nJ < asPotentialSchoolOneStudentToProcess->Count(); nJ++)
    {
        auto prPerson = asPotentialSchoolOneStudentToProcess->Item(nJ);
        bool bIsDropout = (prPerson->educ_one_fate == EOL_MEDIUM );

        
        if (prPerson->integer_age == prPerson->educ_one_entry_age || prPerson->educ_one_status == EOS_PAUSE )
        {
            prPerson->educ_one_status = EOS_WAIT;
        }

        
        else if (prPerson->educ_one_grade_attended == MAX(EDUC_ONE_GRADE))
        {
            
            if ( bIsDropout )  prPerson->educ_one_status = EOS_OUT;
            
            else if (RandUniform(64) > SchoolOneRepetitionRate[prPerson->educ_one_geo]
                [prPerson->educ_one_group][RANGE_POS(SIM_YEAR_RANGE, prPerson->calendar_year)])
            {
                prPerson->educ_one_grade_passed = prPerson->educ_one_grade_attended;
                prPerson->educ_one_status = EOS_OUT;
            }
            
            else if (RandUniform(65) < SchoolOneInterruptionRate[prPerson->educ_one_geo]
                [prPerson->educ_one_group][RANGE_POS(SIM_YEAR_RANGE, prPerson->calendar_year)])
            {
                prPerson->educ_one_status = EOS_PAUSE;
            }
            else prPerson->educ_one_status = EOS_WAIT;
        }

        
        else if (prPerson->integer_age > prPerson->educ_one_entry_age)
        {
            
            if (RandUniform(66) > SchoolOneRepetitionRate[prPerson->educ_one_geo]
                [prPerson->educ_one_group][RANGE_POS(SIM_YEAR_RANGE, prPerson->calendar_year)])
            {
                prPerson->educ_one_grade_passed = prPerson->educ_one_grade_attended;
            }

            
            if (bIsDropout && prPerson->educ_one_grade_fate == prPerson->educ_one_grade_attended)
            {
                prPerson->educ_one_status = EOS_OUT;
            }
            
            else if (RandUniform(67) > SchoolOneInterruptionRate[prPerson->educ_one_geo]
                [prPerson->educ_one_group][RANGE_POS(SIM_YEAR_RANGE, prPerson->calendar_year)])
            {
                prPerson->educ_one_status = EOS_WAIT;
            }
            else prPerson->educ_one_status = EOS_PAUSE;
        }
    }

    
    while (asPotentialSchoolOneStudentToProcess->Count() > 0)
    {
        asPotentialSchoolOneStudentToProcess->Item(0)->educ_one_to_process = FALSE;
    }
    time_end_school_one_year = TIME_INFINITE;
}





void Clock::SetSchoolYearOneClock()
{
    time_start_school_one_year = WAIT(StartSchoolOneYear);
    time_end_school_one_year = WAIT(EndSchoolOneYear);
}

#line 1 "../code/EducationPrimaryTransmission.mpp"












































































































void Calibrator::SetEduc1AdjustmentFactors()
{
    if ( calibrator_year >= MIN(SIM_YEAR_RANGE)  && (
       (Educ1Model == E1M_REFINED_ALIGNALL  && calibrator_year >= Educ1FirstCohortRefinedModel ) ||
       (Educ1Model == E1M_REFINED_ALIGNONCE && calibrator_year == Educ1FirstCohortRefinedModel )))
    {
        for (int nSex = 0; nSex < SIZE(SEX); nSex++ )
        {
            for (int nGeo = 0; nGeo < SIZE(GEO); nGeo++ )
            {
                
                double nTotalPop = 0.0;
                for ( int dGroup = 0; dGroup < SIZE(EDUC1_GROUP); dGroup++ )
                {
                    nTotalPop = nTotalPop + asSimBornAge0[nSex][nGeo][dGroup]->Count();
                }

                
                double dFactorEntry = 0.0;
                if ( nTotalPop > 0.0 )
                {
                    int nIterations = 0;
                    double dResultProb = 10.0;
                    double dTargetProb = EducTrans1[nSex][RANGE_POS(YOB_EDUC_TRANS1, calibrator_year)][nGeo];
                    double dLower = -10.0;
                    double dUpper = 10.0;
                    double dCenter = 0.0;
                    while (abs(dResultProb - dTargetProb) > 0.001 && nIterations < 1000)
                    {
                        nIterations++;
                        dCenter = (dLower + dUpper) / 2.0;
                        dResultProb = 0.0;
                        for ( int nGroup = 0; nGroup < SIZE(EDUC1_GROUP); nGroup++ )
                        {
                            dResultProb = dResultProb + ( asSimBornAge0[nSex][nGeo][nGroup]->Count() / nTotalPop ) *
                                AdjustedProbability(dTargetProb, log(Educ1StartOdds[nGroup][nSex]), dCenter);
                        }
                        if (dTargetProb > dResultProb) dLower = dCenter;
                        else dUpper = dCenter;
                    }
                    dFactorEntry = dCenter;
                }
                
                alignment_educ1_medium[nSex][nGeo] = dFactorEntry;

                
                double dFactorGrad = 0.0;
                if ( nTotalPop > 0.0 )
                {
                    int nIterations = 0;
                    double dResultProb = 10.0;
                    double dTargetProb = EducTrans2[nSex][RANGE_POS(YOB_EDUC_TRANS2, calibrator_year)][nGeo];
                    double dLower = -10.0;
                    double dUpper = 10.0;
                    double dCenter = 0.0;
                    while (abs(dResultProb - dTargetProb) > 0.001 && nIterations < 1000)
                    {
                        nIterations++;
                        dCenter = (dLower + dUpper) / 2.0;
                        dResultProb = 0.0;
                        for ( int nGroup = 0; nGroup < SIZE(EDUC1_GROUP); nGroup++ )
                        {
                            dResultProb = dResultProb + ( asSimBornAge0[nSex][nGeo][nGroup]->Count() / nTotalPop ) *
                                AdjustedProbability(dTargetProb, log(Educ1GradOdds[nGroup][nSex]), dCenter);
                        }
                        if (dTargetProb > dResultProb) dLower = dCenter;
                        else dUpper = dCenter;
                    }
                    dFactorGrad = dCenter;
                }
                
                alignment_educ1_high[nSex][nGeo] = dFactorGrad;
            }
        }
    }
}

double Calibrator::AdjustedProbability(double dProb, double dLogOddEduc, double dLogOddAdjust)
{
    if (dProb <= 0.0) return 0.0;
    else if (dProb >= 1.0) return 1.0;
    else
    {
        if (dProb >= 0.9999) dProb = 0.9999;
        double dValue = log(dProb / (1 - dProb)) + dLogOddEduc + dLogOddAdjust;
        if (dValue > 50) dValue = 50;
        double dExp = exp(dValue);
        return dExp / (1 + dExp);
    }
}

void Person::AdjustEducOne()
{
    if (creation_type == CT_BIRTH && integer_age == 0 && calendar_year >= MIN(SIM_YEAR_RANGE) &&
        calendar_year >= Educ1FirstCohortRefinedModel)
    {
        
        if ( Educ1Model == E1M_REFINED_ALIGNALL )
        {
            double dProb1 = lCalibrator->AdjustedProbability(
                EducTrans1[sex][RANGE_POS(YOB_EDUC_TRANS1,year_of_birth)][geo_birth],
                log(Educ1StartOdds[educ1_group][sex]),
                lCalibrator->alignment_educ1_medium[sex][geo_birth]);
            double dProb2 = lCalibrator->AdjustedProbability(
                EducTrans2[sex][RANGE_POS(YOB_EDUC_TRANS2,year_of_birth)][geo_birth],
                log(Educ1GradOdds[educ1_group][sex]),
                lCalibrator->alignment_educ1_high[sex][geo_birth]);
            EDUC_ONE_LEVEL eolFate = EOL_LOW;
            if (RandUniform(50) < dProb1) eolFate = EOL_MEDIUM;
            if (eolFate == EOL_MEDIUM && RandUniform(51) < dProb2)  eolFate = EOL_HIGH;
            educ_one_fate = eolFate;
        }
        
        else if ( Educ1Model == E1M_REFINED_ALIGNONCE )
        {
            double dProb1 = lCalibrator->AdjustedProbability(
                EducTrans1[sex][RANGE_POS(YOB_EDUC_TRANS1, Educ1FirstCohortRefinedModel)][geo_birth],
                log(Educ1StartOdds[educ1_group][sex]),
                lCalibrator->alignment_educ1_medium[sex][geo_birth]);
            double dProb2 = lCalibrator->AdjustedProbability(
                EducTrans2[sex][RANGE_POS(YOB_EDUC_TRANS2, Educ1FirstCohortRefinedModel)][geo_birth],
                log(Educ1GradOdds[educ1_group][sex]),
                lCalibrator->alignment_educ1_high[sex][geo_birth]);
            EDUC_ONE_LEVEL eolFate = EOL_LOW;
            if (RandUniform(52) < dProb1) eolFate = EOL_MEDIUM;
            if (eolFate == EOL_MEDIUM && RandUniform(53) < dProb2)  eolFate = EOL_HIGH;
            educ_one_fate = eolFate;
        }
    }
}

#line 1 "../code/EducationSecondaryBase.mpp"


































          
        
      


































































































TIME Clock::timeStartSchoolTwoYearEvent(){return time_start_school_two_year;}
void Clock::StartSchoolTwoYearEvent()
{
    for (long nJ = 0; nJ < asPotentialSchoolTwoStudent->Count(); nJ++)
    {
        auto prPerson = asPotentialSchoolTwoStudent->Item(nJ);
        if (prPerson->educ_two_status == ETS_WAIT )
        {
            prPerson->educ_two_status = ETS_ATTEND ;
            prPerson->educ_two_grade_attended = prPerson->educ_two_grade_passed + 1;
        }
    }
    time_start_school_two_year = TIME_INFINITE;
}

TIME Clock::timeEndSchoolTwoYearEvent(){return time_end_school_two_year;}
void Clock::EndSchoolTwoYearEvent()
{
    
    for (long nJ = 0; nJ < asPotentialSchoolTwoStudent->Count(); nJ++)
    {
        auto prPerson = asPotentialSchoolTwoStudent->Item(nJ);
        prPerson->educ_two_to_process = TRUE;
    }

    
    for (long nJ = 0; nJ < asPotentialSchoolTwoStudentToProcess->Count(); nJ++)
    {
        auto prPerson = asPotentialSchoolTwoStudentToProcess->Item(nJ);

        int nGeo         = prPerson->educ_two_geo;
        int nGroup       = prPerson->educ_two_group;
        int nGradeAttend = RANGE_POS(EDUC_TWO_GRADE, prPerson->educ_two_grade_attended);
        int nYear        = RANGE_POS(SIM_YEAR_RANGE, prPerson->calendar_year);


        
        if (prPerson->educ_two_status == ETS_NEVER && ( 
            ( prPerson->educ_two_delay == 0 && RandUniform(68) < Educ2DirectProgressionIntake[nGeo][nGroup][0][nYear] ) ||
            ( prPerson->educ_two_delay == 1 && RandUniform(69) < Educ2DelayedProgressionIntake[nGeo][nGroup][0][nYear] )))
        {
            prPerson->educ_two_status = ETS_WAIT;
            prPerson->educ_two_delay = 0;           
        }
        else if (prPerson->educ_two_status == ETS_NEVER) prPerson->educ_two_delay++;

        
        else if (prPerson->educ_two_grade_attended == MAX(EDUC_TWO_GRADE) && prPerson->educ_two_status == ETS_ATTEND)  
        {
            
            if ( RandUniform(70) < Educ2PeriodSuccess[nGeo][nGroup][nGradeAttend][nYear] )  
            {
                prPerson->educ_two_grade_passed = prPerson->educ_two_grade_attended;
                prPerson->educ_two_status = ETS_OUT;
            }
            
            else if (RandUniform(71) < Educ2DirectRepetitionIntake[nGeo][nGroup][nGradeAttend][nYear] && prPerson->educ_two_delay < Educ2AllowedDelays) 
            {
                prPerson->educ_two_delay++;
                prPerson->educ_two_status = ETS_WAIT;
            }
            
            else if ( prPerson->educ_two_delay + 1 < Educ2AllowedDelays ) 
            {
                prPerson->educ_two_delay++;
                prPerson->educ_two_status = ETS_PAUSE;
            }
            
            else 
            {
                prPerson->educ_two_status = ETS_OUT;
            }
        }

        
        else if (prPerson->educ_two_grade_attended < MAX(EDUC_TWO_GRADE) && prPerson->educ_two_status == ETS_ATTEND)  
        {
            
            if ( RandUniform(72) < Educ2PeriodSuccess[nGeo][nGroup][nGradeAttend][nYear] )  
            {
                prPerson->educ_two_grade_passed = prPerson->educ_two_grade_attended;
            }
            
            if ( prPerson->educ_two_grade_passed == prPerson->educ_two_grade_attended )
            {
                if (RandUniform(73) < Educ2DirectProgressionIntake[nGeo][nGroup][nGradeAttend+1][nYear]) 
                {
                    prPerson->educ_two_status = ETS_WAIT;
                }
                else if (prPerson->educ_two_delay < Educ2AllowedDelays)
                {
                    prPerson->educ_two_delay++;
                    prPerson->educ_two_status = ETS_PAUSE;
                }
                else prPerson->educ_two_status = ETS_OUT;
            }
            
            else
                if (RandUniform(74) < Educ2DirectRepetitionIntake[nGeo][nGroup][nGradeAttend][nYear] 
                    && prPerson->educ_two_delay < Educ2AllowedDelays) 
                {
                    prPerson->educ_two_delay++;
                    prPerson->educ_two_status = ETS_WAIT;
                }
                else if (prPerson->educ_two_delay + 1 < Educ2AllowedDelays)
                {
                    prPerson->educ_two_delay = prPerson->educ_two_delay + 2;
                    prPerson->educ_two_status = ETS_PAUSE;
                }
                else prPerson->educ_two_status = ETS_OUT;
        }
        
        else
        {
            
            if ( prPerson->educ_two_grade_attended == prPerson->educ_two_grade_passed 
                && RandUniform(75) < Educ2DelayedProgressionIntake[nGeo][nGroup][nGradeAttend+1][nYear] ) 
            {
                prPerson->educ_two_status = ETS_WAIT;
            }
            
            else if (RandUniform(76) < Educ2DelayedRepetitionIntake[nGeo][nGroup][nGradeAttend][nYear] ) 
            {
                prPerson->educ_two_status = ETS_WAIT;
            }
            else prPerson->educ_two_status = ETS_OUT;
        }
    }

    
    while (asPotentialSchoolTwoStudentToProcess->Count() > 0)
    {
        asPotentialSchoolTwoStudentToProcess->Item(0)->educ_two_to_process = FALSE;
    }
    time_end_school_two_year = TIME_INFINITE;
}

void Clock::SetSchoolYearTwoClock()
{
    time_start_school_two_year = WAIT(StartSchoolTwoYear);
    time_end_school_two_year = WAIT(EndSchoolTwoYear);
}





void om_PreSimulation_1()
{
    StartSchoolTwoYear = StartSchoolOneYear + 0.000114155;  
    EndSchoolTwoYear = EndSchoolOneYear + 0.000114155;      
};



#line 1 "../code/EmigrationBase.mpp"






                                                     










































TIME Person::timeEmigrationEvent()
{
    if (ModelEmigration && is_resident && is_mortal && EmigrationRatesDistrict[sex][age_mig][(GEO_NAT)geo] > 0.0)
    {
        return WAIT(-log(RandUniform(38)) / EmigrationRatesDistrict[sex][age_mig][(GEO_NAT)geo]);
    }
    else return TIME_INFINITE;
}

void Person::EmigrationEvent()
{
    
    int nDestination; Lookup_EmigrationDestination(RandUniform(56), &nDestination);
    
    
    lCalibrator->GetNextToEmigrate(this,(GEO)(SIZE(GEO_NAT) + nDestination))->doResidentialMove((GEO)(SIZE(GEO_NAT) + nDestination));
}


Person_ptr Calibrator::GetNextToEmigrate(Person_ptr ptrPerson, GEO cGeoTo)
{
    if (MigrationTryKeepingFamiliesTogether)
    {
        if (asWantToMove[ptrPerson->geo][(GEO)cGeoTo][ptrPerson->sex][ptrPerson->age_mig]->Count() > 0)
        {
            return asWantToMove[ptrPerson->geo][(GEO)cGeoTo][ptrPerson->sex][ptrPerson->age_mig]->GetRandom(RandUniform(27));
        }
        else return ptrPerson;
    }
    else return ptrPerson;
}


#line 1 "../code/Ethnicity.mpp"


















































ETHNICITY Person::GetInheritedEthnicity(ETHNICITY eMothersEthnicity)
{
    int nSex = (int)sex;
    int nEthnicity;
    int nMothersEthnicity = (int)eMothersEthnicity;
    Lookup_EthnicTransmission(RandUniform(28), nSex, eMothersEthnicity, &nEthnicity);
    return (ETHNICITY)nEthnicity;
}

ETHNICITY Person::GetImmigrantsScratchEthnicity(GEO_NAT toGeo)
{
    int nSex = (int)sex;
    int nDistrict = (int)toGeo;
    int nEthnicity;
    Lookup_EthnicityImmigrantsScratch(RandUniform(29), nSex, nDistrict, &nEthnicity);
    return (ETHNICITY)nEthnicity;
}


#line 1 "../code/FamilyFemalePartnershipStatus.mpp"






















 










        
    




























































































void Clock::UpdatePartnershipStatus()
{
    long nTarget;
    if (clock_year >= MIN(SIM_YEAR_RANGE))
    {
        
        for (long nJ = 0; nJ < asAllPerson->Count(); nJ++) asAllPerson->Item(nJ)->is_blocked_from_marriage=FALSE;

        
        for (int nGroup = 0; nGroup < SIZE(UNION1_GROUP); nGroup++)
        {
            for (int nChildAge = 0; nChildAge < SIZE(CHILD_AGEGR); nChildAge++)
            {
                for (int nMothAge = 0; nMothAge < SIZE(MOTH_AGEGR); nMothAge++)
                {
                    long nGroupSize = asWomenWithChildren[nGroup][nChildAge][nMothAge][FALSE]->Count()
                                    + asWomenWithChildren[nGroup][nChildAge][nMothAge][TRUE]->Count();

                    nTarget = round(InUnionProbWithChildren[nGroup][nChildAge][nMothAge] * nGroupSize);

                    if (nTarget > asWomenWithChildren[nGroup][nChildAge][nMothAge][TRUE]->Count())
                    {
                        
                        long nSize = asWomenWithChildren[nGroup][nChildAge][nMothAge][FALSE]->Count();
                        long nIndex = 0;
                        for (long nJ = 0; nJ < nSize; nJ++)
                        {
                            auto nPers = asWomenWithChildren[nGroup][nChildAge][nMothAge][FALSE]->Item(nIndex);
                            if (nPers->ever_union) nIndex++;
                            else nPers->is_blocked_from_marriage = TRUE;
                        }

                        
                        while (nTarget > asWomenWithChildren[nGroup][nChildAge][nMothAge][TRUE]->Count() &&
                            asWomenWithChildren[nGroup][nChildAge][nMothAge][FALSE]->Count() > 0) 
                        {
                            auto prFam = asWomenWithChildren[nGroup][nChildAge][nMothAge][FALSE]->GetRandom(RandUniform(55));
                            if (!prFam->FindSpouse()) prFam->is_blocked_from_marriage = TRUE;
                        }
                    }
                    else if (nTarget < asWomenWithChildren[nGroup][nChildAge][nMothAge][TRUE]->Count())
                    {
                        while (nTarget < asWomenWithChildren[nGroup][nChildAge][nMothAge][TRUE]->Count() &&
                            asWomenWithChildren[nGroup][nChildAge][nMothAge][TRUE]->Count() > 0)
                        {
                            auto prFam = asWomenWithChildren[nGroup][nChildAge][nMothAge][TRUE]->GetRandom(RandUniform(57));
                            prFam->doDissolveUnion();
                        }
                    }
                }
            }
        }
        
        for (int nGroup = 0; nGroup < SIZE(UNION1_GROUP); nGroup++)
        {
            for (int nAge = 0; nAge < SIZE(FEMALE_SPOUSE_AGE); nAge++)
            {
                nTarget = round(InUnionProbNoChildren[nAge][nGroup] 
                    * (asWomenNoChildren[nGroup][nAge][FALSE]->Count() + asWomenNoChildren[nGroup][nAge][TRUE]->Count()));

                if (nTarget > asWomenNoChildren[nGroup][nAge][TRUE]->Count())
                {
                    
                    long nSize = asWomenNoChildren[nGroup][nAge][FALSE]->Count();
                    long nIndex = 0;
                    for (long nJ = 0; nJ < nSize; nJ++)
                    {
                        auto nPers = asWomenNoChildren[nGroup][nAge][FALSE]->Item(nIndex);
                        if (nPers->ever_union) nIndex++;
                        else nPers->is_blocked_from_marriage = TRUE;
                    }

                    
                    while (nTarget > asWomenNoChildren[nGroup][nAge][TRUE]->Count() &&
                        asWomenNoChildren[nGroup][nAge][FALSE]->Count() > 0) 
                    {
                        auto prFam = asWomenNoChildren[nGroup][nAge][FALSE]->GetRandom(RandUniform(58));
                        if (!prFam->FindSpouse()) prFam->is_blocked_from_marriage = TRUE;
                    }
                }
                else if (nTarget < asWomenNoChildren[nGroup][nAge][TRUE]->Count())
                {
                    while (nTarget < asWomenNoChildren[nGroup][nAge][TRUE]->Count() &&
                        asWomenNoChildren[nGroup][nAge][TRUE]->Count() > 0)
                    {
                        auto prFam = asWomenNoChildren[nGroup][nAge][TRUE]->GetRandom(RandUniform(39));
                        prFam->doDissolveUnion();
                    }
                }
            }
        }
    }
}





#line 1 "../code/FamilyFirstUnion.mpp"











































































































TIME StartpopValues::timeSetFirstUnionClock()
{
    if (is_activated && !is_set_first_union_clock) return WAIT(0);
    else return TIME_INFINITE;
}

void StartpopValues::SetFirstUnionClock()
{
    TIME dEventTime = TIME_INFINITE;
    if (lPersonStartpop->sex == FEMALE
        && (lPersonStartpop->creation_type == CT_START || lPersonStartpop->creation_type == CT_POOL)
        && StartPopValue[PMC_UNION] > StartPopValue[PMC_BIRTH]
        && StartPopValue[PMC_UNION] < MIN(SIM_YEAR_RANGE))
    {
        dEventTime = StartPopValue[PMC_UNION];
        if (dEventTime == int(dEventTime)) dEventTime = dEventTime + RandUniform(30);
        
        dEventTime = dEventTime + int(lPersonStartpop->time_of_birth) - StartPopValue[PMC_BIRTH];
    }

    time_first_union_imputation = dEventTime;
    is_set_first_union_clock = TRUE;
}

TIME StartpopValues::timeFirstUnionImputationEvent()
{
    return time_first_union_imputation;
}

void StartpopValues::FirstUnionImputationEvent()
{
    lPersonStartpop->ever_union = TRUE;
    time_first_union_imputation = TIME_INFINITE;
}

TIME Person::timeFirstUnionFormationEvent()
{
    if ( sex == FEMALE && !ever_union && WITHIN(FEMALE_SPOUSE_AGE, integer_age) && (is_mortal || creation_type == CT_SCRATCH ))
    {
        double dHazard = Union1FormationHazard[union1_group][RANGE_POS(ALL_YEAR_RANGE, calendar_year)][RANGE_POS(FEMALE_SPOUSE_AGE, integer_age)];
        if (dHazard > 0.0) return WAIT(-TIME(log(RandUniform(31)) / dHazard));
        else return TIME_INFINITE;
    }
    else return TIME_INFINITE;
}

void Person::FirstUnionFormationEvent()
{
    ever_union = TRUE;
    FindSpouse();
}






double GetUnionFormationHazard(int nEduc, int nYear, int nAge);
double GetUnionFormationHazard(int nEduc, int nYear, int nAge)
{
    double	dReturnValue = 0.0;
    double	g[SIZE(FEMALE_SPOUSE_AGE)];
    double	a0 = Union1ParametersCMN[nEduc][RANGE_POS(YOB_UNION, nYear - nAge)][UP_MINAGE];
    double	my = Union1ParametersCMN[nEduc][RANGE_POS(YOB_UNION, nYear - nAge)][UP_AVERAGE];
    double	C = Union1ParametersCMN[nEduc][RANGE_POS(YOB_UNION, nYear - nAge)][UP_EVER];
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
        int		nIndex = int(index);		
        double	dIndex = index - nIndex;		
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
        int		nIndex = int(index);		
        double	dIndex = index - nIndex;		
        double G = 0.0;
        for (int nI = 0; nI <= nIndex; nI++)
        {
            if (nI < SIZE(FEMALE_SPOUSE_AGE)) G = G + g[nI];
        }
        if (nIndex < SIZE(FEMALE_SPOUSE_AGE) - 1) G = G + dIndex*g[nIndex + 1];
        if (1.0 - G * C > 0.0001) dReturnValue = (dThisYearCumulative - G * C) / (1 - G * C);
        else dReturnValue = 0.0;
    }
    
    if (dReturnValue < 1.0)  dReturnValue = -log(1 - dReturnValue);
    return dReturnValue;
}










void om_PreSimulation_2()
{
    for (int nEduc = 0; nEduc < SIZE(UNION1_GROUP); nEduc++)
    {
        for (int nYear = MIN(ALL_YEAR_RANGE); nYear <= MAX(ALL_YEAR_RANGE); nYear++)
        {
            for (int nAge = MIN(FEMALE_SPOUSE_AGE); nAge <= MAX(FEMALE_SPOUSE_AGE); nAge++)
            {
                
                if (Union1Choice == U1C_CMN)
                {
                    Union1FormationHazard[nEduc][RANGE_POS(ALL_YEAR_RANGE, nYear)][RANGE_POS(FEMALE_SPOUSE_AGE, nAge)]
                        = GetUnionFormationHazard(nEduc, nYear, nAge);
                }
                
                else
                {
                    Union1FormationHazard[nEduc][RANGE_POS(ALL_YEAR_RANGE, nYear)][RANGE_POS(FEMALE_SPOUSE_AGE, nAge)]
                        = Union1ParametersHazards[(UNION1_GROUP)nEduc][RANGE_POS(FEMALE_SPOUSE_AGE, nAge)][RANGE_POS(ALL_YEAR_RANGE, nYear)];
                }
            }
        }
    }
}

#line 1 "../code/FamilyGeneral.mpp"




















    
    

      
      

                                    





























































void Person::doStartUnion(Person_ptr ptrPartner)
{
    
    lSpouse = ptrPartner;
    lSpouse->ever_union = TRUE;

    
    doLeaveParentalHome(); lSpouse->doLeaveParentalHome();
    
    
    if (lSpouse->children_in_household)
    {
        int nIndex;
        auto prChild = lSpouse->mlHHFatherChildren->GetNext( 0, &nIndex );
        while ( prChild != NULL )
        {
            prChild->lHHMother = this;
            prChild = lSpouse->mlHHFatherChildren->GetNext(nIndex + 1, &nIndex);
        }
    }

    
    if (children_in_household)
    {
        int nIndex;
        auto prChild = mlHHMotherChildren->GetNext( 0, &nIndex );
        while ( prChild != NULL )
        {
            prChild->lHHFather = lSpouse;
            
            if (!prChild->lBioFather && prChild->integer_age == 0)
            {
                prChild->lBioFather = lSpouse;
            }
            prChild = mlHHMotherChildren->GetNext(nIndex + 1, &nIndex);
        }
    }
}

void Person::doDissolveUnion()
{
    
    if (children_in_household)
    {
        int nIndex;
        bool bStayWithMother = (RandUniform(32) < ProbStayWithMother);

        
        auto prChild = mlHHMotherChildren->GetNext(0, &nIndex);
        while (prChild)
        {
            if (bStayWithMother && ((prChild->lBioMother && prChild->lBioMother == (Person_ptr)this)
                || !prChild->lBioFather || (prChild->lBioFather && prChild->lBioFather != lSpouse)))
            {
                prChild->lHHFather = NULL;
            }
            prChild = mlHHMotherChildren->GetNext(nIndex + 1, &nIndex);
        }

        
        if (lSpouse->children_in_household)
        {
            auto prChild = lSpouse->mlHHFatherChildren->GetNext(0, &nIndex);
            while (prChild)
            {
                prChild->lHHMother = NULL;
                prChild = lSpouse->mlHHFatherChildren->GetNext(nIndex + 1, &nIndex);
            }
        }
    }

    
    
    
    if (!lSpouse->ever_resident && lSpouse->creation_type == CT_SCRATCH && lSpouse->time_of_first_immigration < time
        && int(lSpouse->time_of_first_immigration) == calendar_year)  
    {
        lSpouse->time_of_first_immigration = time; 
    }

    
    lSpouse = NULL;
}

void Person::doLeaveParentalHome()
{
    
    
    
    if (!ever_resident && creation_type == CT_SCRATCH && time_of_first_immigration < time
        && int(time_of_first_immigration) == calendar_year)  
    {
        time_of_first_immigration = time; 
    }

    
    lHHFather = NULL;  
    lHHMother = NULL;  
}

void Person::doMaintainLinksAtDeath()
{
    if (!lSpouse && children_in_household > 0)
    {
        
        int nIndex;
        Person_ptr prGuardian;
        auto prChild = (sex==FEMALE) ? mlHHMotherChildren->GetNext( 0, &nIndex ) : mlHHFatherChildren->GetNext( 0, &nIndex );
        while ( prChild )
        {
            
            if (sex == FEMALE && prChild->lBioFather) prGuardian = prChild->lBioFather;      
            else if (sex == MALE && prChild->lBioMother) prGuardian = prChild->lBioMother;   
            else if (lBioMother) prGuardian = lBioMother;                                    
            else if (lBioFather) prGuardian = lBioFather;                                    
            else prGuardian = NULL;                                                          

            
            if (prGuardian && prGuardian->sex == MALE)
            {
                prChild->lHHFather = prGuardian;
                if (prGuardian->lSpouse) prChild->lHHMother = prGuardian->lSpouse;
            }
            else if (prGuardian && prGuardian->sex == FEMALE)
            {
                prChild->lHHMother = prGuardian;
                if (prGuardian->lSpouse) prChild->lHHFather = prGuardian->lSpouse;
            }

            
            prChild = (sex==FEMALE) ? mlHHMotherChildren->GetNext( nIndex+1, &nIndex ) : mlHHFatherChildren->GetNext( nIndex+1, &nIndex );
        }
    }
    lSpouse = NULL; 
}

void Person::doLinkToFamilyAtStart()
{
    if (creation_type == CT_START || creation_type == CT_POOL)
    {
        if (lStartValues->StartPopValue[PMC_ROLE] == FR_HEAD || lStartValues->StartPopValue[PMC_ROLE] == FR_SPOUSE)
        {
            lSpouse = ptr_creator;
        }
        else if (lStartValues->StartPopValue[PMC_ROLE] == FR_CHILD && ptr_creator->sex == MALE)
        {
            lHHFather = ptr_creator;
            lBioFather = ptr_creator;
        }
        else if (lStartValues->StartPopValue[PMC_ROLE] == FR_CHILD && ptr_creator->sex == FEMALE)
        {
            lHHMother = ptr_creator;
            lBioMother = ptr_creator;
            lBioMother->doIncreaseParity();
            mother_age_at_birth = lBioMother->age;
        }
    }
    else if (creation_type == CT_BIRTH)
    {
        lHHMother = ptr_creator;
        lBioMother = ptr_creator;
        if (lHHMother->lSpouse)
        {
            lHHFather = lHHMother->lSpouse;
            lBioFather = lHHMother->lSpouse;
        }
    }
}
#line 1 "../code/FamilyLeavingHome.mpp"

































TIME Person::timeLeavingHomeEvent()
{
    if (family_role == FR_CHILD) return time_of_birth + AgeLeavingHome;
    else return TIME_INFINITE;
}

void Person::LeavingHomeEvent()
{
    doLeaveParentalHome();
}
#line 1 "../code/FamilyPartnerMatching.mpp"













 



 



   





 



 



   



























































logical Person::FindSpouse()
{
    bool        bFoundSpouse = FALSE;
    double      dExpectedPartners[SIZE(MALE_SPOUSE_AGE)];
    double      dObservedPartners[SIZE(MALE_SPOUSE_AGE)];
    double      dSumExpectedPartners = 0;
    double      dSumObservedPartners = 0;
    double      dGap = 0.0;
    double      dLargestGap = 0.0;
    int         nAgePartner;
    int         nGroupPartner;
    Person_ptr  ptrSpouse = NULL;
    
    
    for (long nI = 0; nI < SIZE(MALE_SPOUSE_AGE); nI++)
    {
        dExpectedPartners[nI] = PartnerAgeDistribution[RANGE_POS(FEMALE_SPOUSE_AGE, integer_age)][nI];
        
        if (is_mortal) dObservedPartners[nI] = asFemaleInUnionByAgeAndPartnerAge[RANGE_POS(FEMALE_SPOUSE_AGE, integer_age)][nI]->Count();
        else dObservedPartners[nI] = asFemaleInUnionByAgeAndPartnerAgeImmiScratch[RANGE_POS(FEMALE_SPOUSE_AGE, integer_age)][nI][RANGE_POS(SIM_YEAR_RANGE,year_of_first_immigration)]->Count();

        dSumExpectedPartners = dSumExpectedPartners + dExpectedPartners[nI];
        dSumObservedPartners = dSumObservedPartners + dObservedPartners[nI];
    }
    for (long nI = 0; nI < SIZE(MALE_SPOUSE_AGE); nI++)
    {
        if (dSumObservedPartners == 0.0) dSumObservedPartners = 1.0;
        dExpectedPartners[nI] = 1.001 * dExpectedPartners[nI] / dSumExpectedPartners;
        dObservedPartners[nI] = dObservedPartners[nI] / dSumObservedPartners;
        dGap = dExpectedPartners[nI] - dObservedPartners[nI];
        if (dExpectedPartners[nI] > 0.0 && dGap > dLargestGap 
            && ((is_mortal && asAvailableMale[nI][geo]->Count() > 0) 
                || (!is_mortal && asAvailableMaleImmiScratch[nI][geo][RANGE_POS(SIM_YEAR_RANGE,year_of_first_immigration)]->Count() > 0)))
        {
            dLargestGap = dGap;
            bFoundSpouse = TRUE;
            nAgePartner = MIN(MALE_SPOUSE_AGE) + nI;
        }
    }
    
    if (bFoundSpouse)
    {
        Lookup_PartnerCharacteristicDistribution(RandUniform(33), spouse_group, &nGroupPartner);

        
        if (is_mortal && asAvailableMaleByType[RANGE_POS(MALE_SPOUSE_AGE,nAgePartner)][nGroupPartner][geo]->Count() > 0 )
        {
            ptrSpouse = asAvailableMaleByType[RANGE_POS(MALE_SPOUSE_AGE,nAgePartner)][nGroupPartner][geo]->GetRandom(RandUniform(36));
        }
        else if (!is_mortal && asAvailableMaleByTypeImmiScratch[RANGE_POS(MALE_SPOUSE_AGE,nAgePartner)][nGroupPartner][geo][RANGE_POS(SIM_YEAR_RANGE,year_of_first_immigration)]->Count() > 0 )
        {
            ptrSpouse = asAvailableMaleByTypeImmiScratch[RANGE_POS(MALE_SPOUSE_AGE,nAgePartner)][nGroupPartner][geo][RANGE_POS(SIM_YEAR_RANGE,year_of_first_immigration)]->GetRandom(RandUniform(34));
        }
        
        
        else if (is_mortal) ptrSpouse = asAvailableMale[RANGE_POS(MALE_SPOUSE_AGE,nAgePartner)][geo]->GetRandom(RandUniform(37));
        else  ptrSpouse = asAvailableMaleImmiScratch[RANGE_POS(MALE_SPOUSE_AGE,nAgePartner)][geo][RANGE_POS(SIM_YEAR_RANGE,year_of_first_immigration)]->GetRandom(RandUniform(35));
        
        doStartUnion(ptrSpouse);
    }
    return bFoundSpouse;
}


#line 1 "../code/FertilityBase.mpp"











































void Person::doIncreaseParity()
{
    
    this_parity_spell = FALSE;
    
    
    if (!is_mortal && (creation_type == CT_START || creation_type == CT_POOL)
        && lStartValues && parity < lStartValues->StartPopValue[PMC_PARITY])
    {
        if (parity < MAX(PARITY_RANGE)) parity = parity + 1;
    }
    
    else if (parity < MAX(PARITY_RANGE)) parity = parity + 1;

    
    this_parity_spell = TRUE;
}

TIME Person::timeFertilityBaseEvent()
{
    TIME dEventTime = TIME_INFINITE;
    double dHazard = 0.0;
    if (is_fertile && FertilityModel != FEM_DETAIL)
    {
        dHazard = AgeSpecificFertilityRate[RANGE_POS(FERTILE_AGE_RANGE, integer_age)][RANGE_POS(SIM_YEAR_RANGE, calendar_year)];
        if (dHazard > 0.0) dEventTime = WAIT(-TIME(log(RandUniform(3)) / dHazard));
    }
    return dEventTime;
}

void Person::FertilityBaseEvent()
{
    HandleFertility();
}





void om_PreSimulation_3()
{
    double dSum;
    for (int nYear = 0; nYear < SIZE(SIM_YEAR_RANGE); nYear++)
    {
        dSum = 0.0;
        
        for (int nAge = 0; nAge < SIZE(FERTILE_AGE_RANGE); nAge++)
        {
            dSum = dSum + AgeSpecificFertility[nAge][nYear];
        }
        
        for (int nAge = 0; nAge < SIZE(FERTILE_AGE_RANGE); nAge++)
        {
            if (dSum > 0.0)
            {
                AgeSpecificFertilityRate[nAge][nYear]
                    = AgeSpecificFertility[nAge][nYear] / dSum * TotalFertilityRate[nYear];
            }
            else AgeSpecificFertilityRate[nAge][nYear] = 0.0;
        }
    }
}
#line 1 "../code/FertilityDetailed.mpp"












































                      
                      
                      






















            
    














































double Person::getTimeToBirth()
{
    TIME   tEventTime = TIME_INFINITE;
    double dHazard = 0.0;

    
    BIRTH1_GROUP cGroup  = (educ_one_fate == EOL_LOW ) ? B1G_00 : (educ_one_fate == EOL_MEDIUM ) ? B1G_01 : B1G_02;
    BIRTH1_LOC   cLocation = GEO_To_BIRTH1_LOC(geo);
    int          nBirthAgePart = SPLIT(integer_age, BIRTH_AGE_PART);    

    if (is_fertile)
    {
        
        if (parity == 0)
        {
            dHazard = FirstBirthRates[cGroup][ever_union]
                [RANGE_POS(FERTILE_AGE_RANGE, integer_age)][cLocation]
                * BirthTrends[RANGE_POS(PARITY_RANGE1, parity+1)][RANGE_POS(SIM_YEAR_RANGE, calendar_year)];
        }
        
        else
        {
            
            dHazard = HigherOrderBirthsPara[time_in_parity][RANGE_POS(PARITY_RANGE2, parity+1)]
                * BirthTrends[RANGE_POS(PARITY_RANGE1, parity+1)][RANGE_POS(SIM_YEAR_RANGE, calendar_year)];

            
            if (nBirthAgePart > 0) 
            {
                dHazard = dHazard *
                    HigherOrderBirthsPara[HBP_AGE35 - 1 + nBirthAgePart][RANGE_POS(PARITY_RANGE2, parity+1)];
            }
            
            if (educ_one_fate == EOL_MEDIUM)
            {
                dHazard = dHazard * HigherOrderBirthsPara[HBP_EDUC1][RANGE_POS(PARITY_RANGE2, parity+1)];
            }
            else if (educ_one_fate == EOL_HIGH)
            {
                dHazard = dHazard * HigherOrderBirthsPara[HBP_EDUC2][RANGE_POS(PARITY_RANGE2, parity+1)];
            }
        }

        
        if (dHazard > 0.0) tEventTime = WAIT(-TIME(log(RandUniform(12)) / dHazard));
    }
    return tEventTime;
}

Person_ptr Calibrator::GetNextToGiveBirth(Person_ptr prPers)
{
    long        nPopSize;
    double      dWaitToBirth = TIME_INFINITE;
    Person_ptr  ptrNextGivingBirth = NULL;
    Person_ptr  ptrThisOne = NULL;

    
    if (FertilityModel == FEM_ALIGNED_AGE) nPopSize = asAllFertilePersonsForFertilityAlignmentByAge[RANGE_POS(FERTILE_AGE_RANGE, prPers->integer_age)][prPers->is_resident]->Count();
    else nPopSize = asAllFertilePersonsForFertilityAlignment[prPers->is_resident]->Count();

    for (long nI = 0; nI < nPopSize; nI++)
    {
        
        if (FertilityModel == FEM_ALIGNED_AGE) ptrThisOne = asAllFertilePersonsForFertilityAlignmentByAge[RANGE_POS(FERTILE_AGE_RANGE, prPers->integer_age)][prPers->is_resident]->Item(nI);
        else ptrThisOne = asAllFertilePersonsForFertilityAlignment[prPers->is_resident]->Item(nI);
        
        
        double dThisWaitTime = ptrThisOne->getTimeToBirth();
        if (dThisWaitTime <= dWaitToBirth)
        {
            dWaitToBirth = dThisWaitTime;
            ptrNextGivingBirth = ptrThisOne;
        }
    }
    return ptrNextGivingBirth;
}


TIME Person::timeFertilityDetailedEvent()
{
    TIME   tEventTime = TIME_INFINITE;
    if (is_fertile && FertilityModel == FEM_DETAIL)
    {
        double dHazard = 0.0;

        
        BIRTH1_GROUP cGroup = (educ_one_fate == EOL_LOW) ? B1G_00 : (educ_one_fate == EOL_MEDIUM) ? B1G_01 : B1G_02;
        BIRTH1_LOC   cLocation = GEO_To_BIRTH1_LOC(geo);
        int          nBirthAgePart = SPLIT(integer_age, BIRTH_AGE_PART);    
    
        
        if (parity == 0)
        {
            dHazard = FirstBirthRates[cGroup][ever_union]
                [RANGE_POS(FERTILE_AGE_RANGE, integer_age)][cLocation]
                * BirthTrends[RANGE_POS(PARITY_RANGE1, parity+1)][RANGE_POS(SIM_YEAR_RANGE, calendar_year)];
        }
        
        else
        {
            
            dHazard = HigherOrderBirthsPara[time_in_parity][RANGE_POS(PARITY_RANGE2, parity+1)]
                * BirthTrends[RANGE_POS(PARITY_RANGE1, parity+1)][RANGE_POS(SIM_YEAR_RANGE, calendar_year)];

            
            if (nBirthAgePart > 0) 
            {
                dHazard = dHazard *
                    HigherOrderBirthsPara[HBP_AGE35 - 1 + nBirthAgePart][RANGE_POS(PARITY_RANGE2, parity+1)];
            }
            
            if (educ_one_fate == EOL_MEDIUM)
            {
                dHazard = dHazard * HigherOrderBirthsPara[HBP_EDUC1][RANGE_POS(PARITY_RANGE2, parity+1)];
            }
            else if (educ_one_fate == EOL_HIGH)
            {
                dHazard = dHazard * HigherOrderBirthsPara[HBP_EDUC2][RANGE_POS(PARITY_RANGE2, parity+1)];
            }
        }

        
        if (dHazard > 0.0) tEventTime = WAIT(-TIME(log(RandUniform(13)) / dHazard));
    }
    return tEventTime;
}

void Person::FertilityDetailedEvent()
{
    HandleFertility();
}
#line 1 "../code/FertilityGeneral.mpp"






































































TIME StartpopValues::timeSetLastBirthClock()
{
    if (is_activated && !is_set_last_birth_clock) return WAIT(0);
    else return TIME_INFINITE;
}

void StartpopValues::SetLastBirthClock()
{
    TIME dEventTime = TIME_INFINITE;
    if (lPersonStartpop->sex == FEMALE
        && (lPersonStartpop->creation_type == CT_START || lPersonStartpop->creation_type == CT_POOL)
        && StartPopValue[PMC_LASTBIR] > StartPopValue[PMC_BIRTH]
        && StartPopValue[PMC_LASTBIR] < MIN(SIM_YEAR_RANGE))
    {
        dEventTime = StartPopValue[PMC_LASTBIR];
        if (dEventTime == int(dEventTime)) dEventTime = dEventTime + RandUniform(54);
        
        dEventTime = dEventTime + int(lPersonStartpop->time_of_birth) - StartPopValue[PMC_BIRTH];
    }

    time_last_birth_imputation = dEventTime;
    is_set_last_birth_clock = TRUE;
}

TIME StartpopValues::timeLastBirthImputationEvent() { return time_last_birth_imputation; }
void StartpopValues::LastBirthImputationEvent()
{
    lPersonStartpop->this_parity_spell = FALSE;
    if (StartPopValue[PMC_PARITY] > lPersonStartpop->parity) lPersonStartpop->parity = COERCE(PARITY_RANGE,int(StartPopValue[PMC_PARITY]));
    lPersonStartpop->this_parity_spell = TRUE;
    time_last_birth_imputation = TIME_INFINITE;
}

void Person::HandleFertility()
{
    
    if ( FertilityModel == FEM_BASE || FertilityModel == FEM_DETAIL ) doGiveBirth();

    
    else lCalibrator->GetNextToGiveBirth(this)->doGiveBirth();
}

void Person::doGiveBirth()
{
    auto peChild = new Person;  
    
    doIncreaseParity();
        
    peChild->Start(NULL, this, TIME_INFINITE, (SEX)0);   
    doLeaveParentalHome();
}
#line 1 "../code/HumanCapitalIndex.mpp"
















































































































void Person::CalculateHCIVariables()
{
    if (in_hci_sample)
    {
        
        double dQ = SchoolQuality[region_birth][SQP_AV]
            + RandNormal(82) * SchoolQuality[region_birth][SQP_SD];
        if (dQ > 0.0 && age >= 5) quality_of_schooling = dQ;

        
        
        
        if (integer_age >= 60) adult_survival = 1; else adult_survival = 0.0;

        
        ind_hci = survived_early_years
            * exp(HCICoefficients[HCI_EDUC] * (years_of_schooling * quality_of_schooling - 14.0))
            * exp((HCICoefficients[HCI_ASR] * (adult_survival - 1.0) + HCICoefficients[HCI_STUNT] * (!is_stunted - 1.0)) / 2.0);


    }
}



#line 1 "../code/ImmigrationBackMigration.mpp"



































TIME Person::timeBackMigrationEvent()
{
    if (!is_resident && ModelBackmigration && ever_resident && calendar_year >= MIN(SIM_YEAR_RANGE) && BackMigrationHazard > 0)
    {
        return WAIT(-TIME(log(RandUniform(16)) / BackMigrationHazard));
    }
    else return TIME_INFINITE;
};

void Person::BackMigrationEvent()
{
    doResidentialMove(geo_prev);
}
#line 1 "../code/ImmigrationFromPool.mpp"






































TIME Person::timeFirstImmigrationFromPoolEvent()
{
    
    if (!ever_resident && creation_type == CT_POOL ) return time_of_first_immigration;
    else return TIME_INFINITE;
}

void Person::FirstImmigrationFromPoolEvent()
{
    
    int nGeo;

    Lookup_ImmiPoolDestination(RandUniform(18), (int)lStartValues->StartPopValue[PMC_POOL] - 1, RANGE_POS(SIM_YEAR_RANGE,calendar_year), &nGeo);

    
    doResidentialMove((GEO)nGeo);

    
    if (lSpouse) lSpouse->doResidentialMove((GEO)nGeo);
    if (children_in_household > 0)
    {
        int nIndex;
        auto prChild = (sex == FEMALE) ? mlHHMotherChildren->GetNext(0, &nIndex) : mlHHFatherChildren->GetNext(0, &nIndex);
        while (prChild)
        {
            prChild->doResidentialMove((GEO)nGeo);
            
            prChild = (sex == FEMALE) ? mlHHMotherChildren->GetNext(nIndex + 1, &nIndex) : mlHHFatherChildren->GetNext(nIndex + 1, &nIndex);
        }
    }
}

#line 1 "../code/ImmigrationFromScratch.mpp"






































































void Person::ImputeCharacteristicsAtFirstImmigrationScratch(GEO_NAT cGeoDestination)
{
    ethnicity = GetImmigrantsScratchEthnicity(cGeoDestination);
}

void Person::SetGeobirthTimeofbirthCtScratch()
{
    int nAgeAtImmigration;
    Lookup_AgeImmigrantsScratch(RandUniform(14), sex, &nAgeAtImmigration);
    time_of_birth = time_of_first_immigration - nAgeAtImmigration - RandUniform(5);
    
    
    geo_birth = (GEO)int(SIZE(GEO_NAT) + int(RandUniform(8) * (SIZE(GEO) - SIZE(GEO_NAT))));
}

TIME Person::timeFirstImmigrationFromScratchEvent()
{
    if (!ever_resident && creation_type == CT_SCRATCH && family_role == FR_HEAD)
    {
        return time_of_first_immigration;
    }
    else return TIME_INFINITE;
}





















void Person::FirstImmigrationFromScratchEvent()
{
    
    int nGeo = 0; Lookup_ImmiScratchDestination(RandUniform(24),sex, age_mig, &nGeo);

    
    ImputeCharacteristicsAtFirstImmigrationScratch((GEO_NAT)nGeo);
    doResidentialMove((GEO)nGeo);

    
    if (lSpouse)
    {
        lSpouse->ImputeCharacteristicsAtFirstImmigrationScratch((GEO_NAT)nGeo);
        lSpouse->doResidentialMove((GEO)nGeo);
    }
    if (children_in_household > 0)
    {
        int nIndex;
        auto prChild = (sex == FEMALE) ? mlHHMotherChildren->GetNext(0, &nIndex) : mlHHFatherChildren->GetNext(0, &nIndex);
        while (prChild)
        {
            prChild->ImputeCharacteristicsAtFirstImmigrationScratch((GEO_NAT)nGeo);
            prChild->doResidentialMove((GEO)nGeo);
            
            prChild = (sex == FEMALE) ? mlHHMotherChildren->GetNext(nIndex + 1, &nIndex) : mlHHFatherChildren->GetNext(nIndex + 1, &nIndex);
        }
    }
}

void Person::FindImmigrantMother()
{
    if (time_of_first_immigration - time_of_birth < AgeImmiSearchMother)
    {
        int nAge;
        int nCount = 0;
        bool bFound = FALSE;
        while (!bFound && nCount < 100)
        {
            Lookup_AgeOfImmigrantMother(RandUniform(26), &nAge);
            nAge = nAge + MIN(FERTILE_AGE_RANGE);
            if (asPotentialImmigrantMothers[nAge][RANGE_POS(SIM_YEAR_RANGE, year_of_first_immigration)][geo_birth]->Count() > 0) bFound = TRUE;
            nCount++;
        }
        if (bFound)
        {
            auto prMother = asPotentialImmigrantMothers[nAge][RANGE_POS(SIM_YEAR_RANGE, year_of_first_immigration)][geo_birth]->Item(RandUniform(40));
            lBioMother = prMother;
            lHHMother = prMother;
            lHHMother->doIncreaseParity();
            mother_age_at_birth = lBioMother->age;
        }
    }
}


#line 1 "../code/MigrationBase.mpp"









































TIME Person::timeMigrationEvent()
{
    
    if (ModelMigration && is_mortal && is_resident && calendar_year != (int)time_last_move)
    {
        
        double dMoveProb = MigrationProbability[sex][age_mig][(GEO_NAT)geo];

        if (dMoveProb <= 0.0) return TIME_INFINITE;     
        else if (dMoveProb >= 1.0) return WAIT(0);      
        else                                            
        {
            
            
            return WAIT(-log(RandUniform(60)) / -log(1 - dMoveProb));
        }
    }
    return TIME_INFINITE;
}

void Person::MigrationEvent()
{
    int nDestination;

    
    Lookup_MigrationDestination(RandUniform(61), sex, (GEO_NAT)geo, int(age_mig), &nDestination);

    
    lCalibrator->GetNextToMigrate(this,(GEO_NAT)nDestination)->doResidentialMove((GEO)nDestination);
}

Person_ptr Calibrator::GetNextToMigrate(Person_ptr ptrPerson, GEO_NAT cGeoTo)
{
    if (MigrationTryKeepingFamiliesTogether)
    {
        if (asWantToMove[ptrPerson->geo][(GEO)cGeoTo][ptrPerson->sex][ptrPerson->age_mig]->Count() > 0)
        {
            return asWantToMove[ptrPerson->geo][(GEO)cGeoTo][ptrPerson->sex][ptrPerson->age_mig]->GetRandom(RandUniform(25));
        }
        else return ptrPerson;
    }
    else return ptrPerson;
}




#line 1 "../code/MigrationGeneral.mpp"
















	    
                                                  




































































void Person::doResidentialMove(GEO cDestination)
{
    
    geo_prev = geo;
    geo = cDestination;

    
    time_last_move = time;
    if ( time_first_move == TIME_INFINITE ) time_first_move = time;

    
    if (!ever_resident && geo < SIZE(GEO_NAT)) ever_resident = TRUE;    

    
    if (geo_prev < SIZE(GEO_NAT) && geo < SIZE(GEO_NAT)) lCalibrator->migration_counter++;
    else if (geo_prev >= SIZE(GEO_NAT) && geo < SIZE(GEO_NAT)) lCalibrator->immigration_counter++;
    else if (geo_prev < SIZE(GEO_NAT) && geo >= SIZE(GEO_NAT)) lCalibrator->emigration_counter++;

    
    if (MigrationTryKeepingFamiliesTogether)
    {
        if (lSpouse) lSpouse->geo_want_to_move = geo;
        if (lHHFather) lHHFather->geo_want_to_move = geo;
        if (lHHMother) lHHMother->geo_want_to_move = geo;
        if (children_in_household > 0)
        {
            int nIndex;
            auto prChild = (sex == FEMALE) ? mlHHMotherChildren->GetNext(0, &nIndex) : mlHHFatherChildren->GetNext(0, &nIndex);
            while (prChild)
            {
                prChild->geo_want_to_move = geo;
                
                prChild = (sex == FEMALE) ? mlHHMotherChildren->GetNext(nIndex + 1, &nIndex) : mlHHFatherChildren->GetNext(nIndex + 1, &nIndex);
            }
        }
    }
}

TIME StartpopValues::timeSetResidentialMoveClock()
{
    if (is_activated && !is_set_residential_move_clock) return WAIT(0);
    else return TIME_INFINITE;
}

void StartpopValues::SetResidentialMoveClock()
{
    
    if (lPersonStartpop->creation_type == CT_START && StartPopValue[PMC_GEOBIR] != StartPopValue[PMC_GEOPRE])
    {
        TIME dTimeFirst = StartPopValue[PMC_MOVEFIRST];
        if (dTimeFirst > lPersonStartpop->time_of_birth && dTimeFirst < MIN(SIM_YEAR_RANGE)) scheduled_time_first_move = dTimeFirst;
        else scheduled_time_first_move = (lPersonStartpop->time_of_birth + (int)lPersonStartpop->time_of_birth + 1) / 2.0; 
    }
    else scheduled_time_first_move = TIME_INFINITE;

    
    if (lPersonStartpop->creation_type == CT_START && StartPopValue[PMC_GEO] != StartPopValue[PMC_GEOPRE])
    {
        TIME dTimeLast = StartPopValue[PMC_MOVELAST];
        if (dTimeLast > lPersonStartpop->time_of_birth && dTimeLast < MIN(SIM_YEAR_RANGE) && (scheduled_time_first_move == TIME_INFINITE || scheduled_time_first_move < dTimeLast)) scheduled_time_last_move = dTimeLast;
        else
        {
            if (scheduled_time_first_move == TIME_INFINITE) scheduled_time_last_move = (lPersonStartpop->time_of_birth + MIN(SIM_YEAR_RANGE)) / 2.0;
            else scheduled_time_last_move = (scheduled_time_first_move + MIN(SIM_YEAR_RANGE)) / 2.0; 
        }
    }
    else scheduled_time_last_move = TIME_INFINITE;

    
    is_set_residential_move_clock = TRUE;
}

TIME StartpopValues::timeFirstMoveEvent() 
{
    
    return scheduled_time_first_move; 
}

void StartpopValues::FirstMoveEvent()
{
    lPersonStartpop->doResidentialMove((GEO)(int)StartPopValue[PMC_GEOPRE]);
    scheduled_time_first_move = TIME_INFINITE;
}


TIME StartpopValues::timeLastMoveEvent() 
{
    
    return scheduled_time_last_move;
}
void StartpopValues::LastMoveEvent() 
{ 
    lPersonStartpop->doResidentialMove((GEO)(int)StartPopValue[PMC_GEO]);
    scheduled_time_last_move = TIME_INFINITE;
}


#line 1 "../code/MortalityBase.mpp"







































TIME Person::timeMortalityBaseEvent()
{
    TIME    dEventTime = TIME_INFINITE;
    double  dMortalityHazard = MortalityTable[integer_age][sex]
          * MortalityTrend[RANGE_POS(SIM_YEAR_RANGE, calendar_year)][sex];    
    
    
    if (is_mortal && dMortalityHazard > 0.0 && 
        (calendar_year < MIN(CHILD_MORTALITY_YEARS) || MortalityModel == MOM_BASE || integer_age > MAX(CHILD_MORTALITY_AGE)))
    {
        
        dEventTime = WAIT(-log(RandUniform(2)) / dMortalityHazard);
    }
    return dEventTime;
}

void Person::MortalityBaseEvent() {  Finish(); }





void om_PreSimulation_4()
{
    double	dLower, dUpper, dCenter, dResult, dTarget, dAlive, dDeaths;
    int		nIterations;
    for (int nSex = 0; nSex < SIZE(SEX); nSex++)
    {
        for (int nYear = 0; nYear < SIZE(SIM_YEAR_RANGE); nYear++)
        {
            dTarget = LifeExpectancy[nYear][nSex];      
            dResult = 0.0;                              
            nIterations = 10000;                        
            dLower = 0.1;                               
            dUpper = 3.0;                               

            while (abs(dResult - dTarget) > 0.0001 && nIterations > 0)
            {
                nIterations--;
                dCenter = (dLower + dUpper) / 2.0;      

                dResult = 0.0;
                dAlive = 1.0;                           

                
                for (int nAge = 0; nAge < SIZE(AGE_RANGE); nAge++)
                {
                    
                    dDeaths = dAlive * (1 - exp(-MortalityTable[nAge][nSex] * dCenter));
                    dAlive = dAlive - dDeaths;
                    
                    dResult = dResult + dDeaths * 0.5 + dAlive;
                }
                
                if (dTarget < dResult) dLower = dCenter;
                else dUpper = dCenter;
            }
            
            MortalityTrend[nYear][nSex] = dCenter;
        }
    }
}

#line 1 "../code/MortalityDetailed.mpp"
























               


















































































TIME Person::timeMortalityDetailedEvent()
{
    double dEventTime = TIME_INFINITE;
    double dHazard = 0.0;

    if ( is_mortal && MortalityModel != MOM_BASE && integer_age <= MAX(CHILD_MORTALITY_AGE) 
        && calendar_year >= MIN(CHILD_MORTALITY_YEARS) )
    {
        
        
        if (MortalityModel == MOM_DETAIL)
        {
            dHazard = ChildMortalityBaseRisk[integer_age][sex]
                * ChildMortalityTrend[integer_age][RANGE_POS(CHILD_MORTALITY_YEARS, calendar_year)];
        }

        
        else if (MortalityModel == MOM_ALIGNED_MACRO_TRENDS)
        {
            dHazard = child_mortality
                * MortalityTrend[RANGE_POS(SIM_YEAR_RANGE, calendar_year)][sex]
                / MortalityTrend[RANGE_POS(SIM_YEAR_RANGE, MIN(CHILD_MORTALITY_YEARS))][sex];
        }

        
        else
        {
            dHazard = child_mortality
                * ChildMortalityTrend[integer_age][RANGE_POS(CHILD_MORTALITY_YEARS, calendar_year)];
        }

        
        dHazard = dHazard * ChildMortalityRelativeRisks[integer_age][child_mortality_group];

        
        if (dHazard > 0) dEventTime = WAIT(-log(RandUniform(9)) / dHazard);;
    }
    return dEventTime;
}

void Person::MortalityDetailedEvent() { Finish(); }


TIME Calibrator::timeChildMortalityCalibration()
{
    if (!is_calibrated_child_mortality) return MIN(CHILD_MORTALITY_YEARS);
    else return TIME_INFINITE;
}

void Calibrator::ChildMortalityCalibration()
{
    
    double dDeaths[SIZE(CHILD_MORTALITY_AGE)][SIZE(SEX)];      
    double dProb[SIZE(CHILD_MORTALITY_AGE)][SIZE(SEX)];        
    double dBase[SIZE(CHILD_MORTALITY_AGE)][SIZE(SEX)];        
    long   nPop[SIZE(CHILD_MORTALITY_AGE)][SIZE(SEX)][SIZE(CHILD_MORTALITY_GROUP)]; 

    
    
    
    for (int nAge = 0; nAge < SIZE(CHILD_MORTALITY_AGE); nAge++)
    {
        for (int nSex = 0; nSex < SIZE(SEX); nSex++)
        {
            dProb[nAge][nSex] = 1.0 - exp(-MortalityTable[nAge][nSex]
                * MortalityTrend[RANGE_POS(SIM_YEAR_RANGE, MIN(CHILD_MORTALITY_YEARS))][nSex]);
            dDeaths[nAge][nSex] = 0.0;
            for (int nGroup = 0; nGroup < SIZE(CHILD_MORTALITY_GROUP); nGroup++)
            {
                nPop[nAge][nSex][nGroup] = 0.0;
            }
        }
    }

    
    
    
    for (long nJ = 0; nJ < asAllPerson->Count(); nJ++)
    {
        auto prPerson = asAllPerson->Item(nJ);
        if (prPerson->integer_age <= MAX(CHILD_MORTALITY_AGE) && prPerson->creation_type==CT_BIRTH && prPerson->is_resident)
        {
            dDeaths[prPerson->integer_age][prPerson->sex] = dDeaths[prPerson->integer_age][prPerson->sex] + dProb[prPerson->integer_age][prPerson->sex];
            nPop[prPerson->integer_age][prPerson->sex][prPerson->child_mortality_group]++;
        }
    }

    
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

                
                for (int nGroup = 0; nGroup < SIZE(CHILD_MORTALITY_GROUP); nGroup++)
                {
                    dNumberDeaths = dNumberDeaths + nPop[nAge][nSex][nGroup]
                        * (1-exp(-dBase[nAge][nSex] *  ChildMortalityRelativeRisks[nAge][nGroup]));
                }
                
                if ( dNumberDeaths > dDeaths[nAge][nSex]) dUpper = dBase[nAge][nSex];
                else dLower = dBase[nAge][nSex];
            }
        }
    }

    
    mort_male_0 = dBase[0][MALE]; mort_female_0 = dBase[0][FEMALE];
    mort_male_1 = dBase[1][MALE]; mort_female_1 = dBase[1][FEMALE];
    mort_male_2 = dBase[2][MALE]; mort_female_2 = dBase[2][FEMALE];
    mort_male_3 = dBase[3][MALE]; mort_female_3 = dBase[3][FEMALE];
    mort_male_4 = dBase[4][MALE]; mort_female_4 = dBase[4][FEMALE];

    is_calibrated_child_mortality = TRUE;
}

#line 1 "../code/MortalityGeneral.mpp"


























































TIME Person::timeDeathAtMaxLifespanEvent() { return time_of_birth + MAX(AGE_RANGE) + 0.99999; }
void Person::DeathAtMaxLifespanEvent() { Finish(); }




#line 1 "../code/Stunting.mpp"





































void Person::DecideStuntingFate()
{
    if (in_projected_time && creation_type == CT_BIRTH && is_resident)
    {
        if (RandUniform(81) < ProportionStunting[sex][region_nat][educ_mother])
        {
            is_stunted = TRUE;
        }
    }
};

#line 1 "../code/TablesEducation.mpp"









	






































































































































































































































#line 1 "../code/TablesFamily.mpp"






















#line 1 "../code/TablesHumanCapital.mpp"






































































#line 1 "../code/TablesPopulation.mpp"










                                       





































































#line 1 "../code/TablesStunting.mpp"




































#line 1 "../code/TablesVaccination.mpp"























































#line 1 "../code/TrackingList.mpp"













#line 1 "../code/_ContextABC.mpp"


















































                 
                 
                      
               
               
          


          













































































































































































































































                 
                  

                         
                 















#line 1 "../code/_ContextNPL.mpp"














































































































































































































































































































































































































































































































































































































































#line 1 "../code/model_core.mpp"
















        
                           









void Simulation()
{
    extern void LogSimulationStart(); 
    extern void SimulateEvents(); 

    
    LogSimulationStart();

    
    
    
    
    
    auto paClock = new Clock(); paClock->Start();

    
    auto paCalibrator = new Calibrator(); paCalibrator->Start();

    
    input_csv in_csv; in_csv.open(MicroDataInputFile); in_csv.read_header();
    for (long nJ = 0; nJ < MicroDataInputFileSize; nJ++)
    {
        in_csv.read_record(nJ);
        auto paObservation = new Observation();
        paObservation->Start(in_csv);
    }
    in_csv.close();

    
    for (long nJ = 0; nJ < MicroDataInputFileSize; nJ++)
    {
        
        if (asObservationByFamOldest[asObservations->Item(nJ)->fam_id]->Count() == 0)
        {
            Observation_ptr ptrThisObservation = asObservationByFam[asObservations->Item(nJ)->fam_id]->Item(0);
            Observation_ptr ptrOldestObservation = ptrThisObservation;
            double dEarliestBirthdate = ptrOldestObservation->obs_birth;
            for (int nI = 0; nI < asObservationByFam[asObservations->Item(nJ)->fam_id]->Count(); nI++)
            {
                ptrThisObservation = asObservationByFam[asObservations->Item(nJ)->fam_id]->Item(nI);
                if (ptrThisObservation->obs_birth <= dEarliestBirthdate)
                {
                    dEarliestBirthdate = ptrThisObservation->obs_birth;
                    ptrOldestObservation = ptrThisObservation;
                }
            }
            ptrOldestObservation->obs_oldest = TRUE;
        }
    }

    
    

    
    for (long nJ = 0; nJ < asObservationHeads[PP_NON]->Count(); nJ++)
    {
        
        double  dWeight = asObservationHeads[PP_NON]->Item(nJ)->pmc[PMC_WEIGHT] / ScalingFactor;
        int nWeight = (int)dWeight; if (RandUniform(1) < dWeight - nWeight) nWeight++;
        asObservationHeads[PP_NON]->Item(nJ)->obs_weight = nWeight;
        
        for (int nX = 0; nX < asObservationNonHeads[asObservationHeads[PP_NON]->Item(nJ)->fam_id][PP_NON]->Count(); nX++)
        {
            asObservationNonHeads[asObservationHeads[PP_NON]->Item(nJ)->fam_id][PP_NON]->Item(nX)->obs_weight = nWeight;
        }
    }

    
    while (asSimulatedObservationHeads[PP_NON]->Count() > 0)
    {
        
        auto paObservation = asSimulatedObservationHeads[PP_NON]->Item(0);
        auto paPerson = new Person();
        paPerson->Start(paObservation, NULL, 0, MALE);
        
        for (int nJ = 0; nJ < asObservationNonHeads[paObservation->fam_id][PP_NON]->Count(); nJ++)
        {
            auto paOtherPerson = new Person();
            paOtherPerson->Start(asObservationNonHeads[paObservation->fam_id][PP_NON]->Item(nJ), paPerson, 0, MALE);
        }
        paObservation->obs_weight--;
    }

    
    while (asObservationAll[PP_NON]->Count() > 0) asObservationAll[PP_NON]->Item(0)->Finish();

    
    

    if (ModelImmigrationFromPools)
    {
        for (int nPool = 1; nPool < SIZE(POP_POOL); nPool++)
        {
    		for (int nYear = MIN(SIM_YEAR_RANGE); nYear < int(SIMULATION_END()); nYear++)
            {
                double dTargetPop = (ImmiPoolSize[nPool - 1][RANGE_POS(SIM_YEAR_RANGE, nYear)]) / ScalingFactor;
                double dSumWeights = 0.0;
                for (long nJ = 0; nJ < asObservationHeads[nPool]->Count(); nJ++)
                {
                    dSumWeights = dSumWeights + asObservationHeads[nPool]->Item(nJ)->pmc[PMC_WEIGHT] *
                        (1 + asObservationNonHeads[asObservationHeads[nPool]->Item(nJ)->fam_id][nPool]->Count());
                }
                for (long nJ = 0; nJ < asObservationHeads[nPool]->Count(); nJ++)
                {
                    double dWeight = asObservationHeads[nPool]->Item(nJ)->pmc[PMC_WEIGHT] * dTargetPop / dSumWeights;
                    int nWeight = int(dWeight);
                    if (RandUniform(6) < (dWeight - nWeight))
                    {
                        nWeight++;
                    };
                    asObservationHeads[nPool]->Item(nJ)->obs_weight = nWeight;
                }
                while (asSimulatedObservationHeads[nPool]->Count() > 0)
                {
                    
                    auto paObservation = asSimulatedObservationHeads[nPool]->Item(0);
                    auto paPerson = new Person();
                    paPerson->Start(paObservation, NULL, nYear, MALE);

                    
                    for (int nJ = 0; nJ < asObservationNonHeads[paObservation->fam_id][nPool]->Count(); nJ++)
                    {
                        auto paOtherPerson = new Person();
                        paOtherPerson->Start(asObservationNonHeads[paObservation->fam_id][nPool]->Item(nJ), 
                            paPerson, nYear, MALE);
                    }
                    paObservation->obs_weight--;
                }
            }
        }
    }

    
    while (asObservations->Count() > 0) asObservations->Item(0)->Finish();


    
    

	if (ModelImmigrationFromScratch)
	{
        for (int nYear = MIN(SIM_YEAR_RANGE); nYear < int(SIMULATION_END()); nYear++)
        {
            for (int nSex = 0; nSex < SIZE(SEX); nSex++)
            {
                int nCohortSize = int(NumberImmigrantsFromScratch[RANGE_POS(SIM_YEAR_RANGE, nYear)][(SEX)nSex]
                    / ScalingFactor);
                double dCohortSize = NumberImmigrantsFromScratch[RANGE_POS(SIM_YEAR_RANGE, nYear)][(SEX)nSex]
                    / ScalingFactor;
                if (RandUniform(4) < dCohortSize - nCohortSize) nCohortSize++;

                for (int nIndex = 0; nIndex < nCohortSize; nIndex++)
                {
                    auto prPerson = new Person();
                    prPerson->Start(NULL, NULL, nYear, (SEX)nSex);
                }
            }
        }
	}

    
    
    

    
    SimulateEvents();
}

#line 1 "../code/model_info.mpp"




















#line 1 "../code/ompp_framework.ompp"










#include "omc/omSimulation.h" // For IDE
#if 0 // Hide from IDE






#endif





namespace fmk 
{
    
    const int model_streams = 100;
}

double population_scaling_factor()
{
    return ScalingFactor;
}

#line 1 "C:/Users/Martin/Dropbox/openM++/use/common.ompp"










namespace fmk {

    


    thread_local openm::IModel * i_model = nullptr;

    


    int simulation_members = 0;

    


    thread_local int simulation_member = 0;

    


    thread_local long member_entity_counter = 0;

    



    const long lcg_modulus = 2147483647;

    


    const long model_stream_seed_generator = 376740290;

    







    thread_local long long combined_seed = 0;

    



















    thread_local long master_seed = 0;

    




    long set_population_value = 0;

    





    thread_local bool do_exit_simulation_all = false;

    


    const int progress_percent_default = 1;
}










long long get_combined_seed()
{
    return fmk::combined_seed;
}







void initialize_model_streams()
{
    using namespace fmk;

    long stream_seed = fmk::master_seed;
    for (int model_stream = 0; model_stream < model_streams; model_stream++) {
        
        
        initialize_stream(model_stream, simulation_member, stream_seed);
		
		
        long long product = model_stream_seed_generator;
        product *= stream_seed;
        stream_seed = product % lcg_modulus;
    }
}




void handle_backwards_time(
    double the_event_time,
    double the_current_time, 
    int the_event, 
    int the_entity)
{
    
    
    extern const char * event_id_to_name(int event_id); 
    std::stringstream ss;
    ss  << std::setprecision(std::numeric_limits<long double>::digits10 + 1) 
        << LT("error : Event time ") << std::showpoint << the_event_time
        << LT(" is earlier than current time ") << the_current_time
        << LT(" in event ") << event_id_to_name(the_event)
        << LT(" in entity_id ") << the_entity
        << LT(" in simulation member ") << get_simulation_member()
        << LT(" with combined seed ") << get_combined_seed()
        ;
    ModelExit(ss.str().c_str());
}




void handle_streams_exceeded(
    int strm, 
    int model_streams)
{
    
    std::stringstream ss;
    ss  << LT("error : stream number ") << strm
        << LT(" exceeds the maximum number of streams ") << model_streams << LT(".")
        << LT(" Increase the number of streams in ompp_framework.ompp.")
        ;
    ModelExit(ss.str().c_str());
}






int get_simulation_members()
{
    return fmk::simulation_members;
}






int get_simulation_member()
{
    return fmk::simulation_member;
}










int get_next_entity_id()
{
    fmk::member_entity_counter++;
    return fmk::member_entity_counter * fmk::simulation_members + fmk::simulation_member;
}







void ModelExit(const char* msg)
{
    theLog->logMsg(msg);
    throw openm::SimulationException(msg);
}





void report_simulation_progress(int member, int percent)
{
    theLog->logFormatted("member=%d Simulation progress=%d%%", member, percent);
    report_simulation_progress_beat(percent);
}





void report_simulation_progress_beat(int percent, double value)
{
    fmk::i_model->updateProgress(percent, value);
}









void ProgressMsg(const char* msg)
{
    
}









void TimeReport(double tm)
{
    
}






void SetMaxTime(double max_value)
{
    fixed_precision<Time::type>::set_max(max_value);
}






void StartEventTrace()
{
    if (BaseEvent::trace_event_enabled) {
        BaseEvent::trace_event_on = true;
    }
}






void StopEventTrace()
{
    if (BaseEvent::trace_event_enabled) {
        BaseEvent::trace_event_on = false;
    }
}






int GetThreads()
{
    return 1;
}






int GetThreadNumber()
{
    return 1;
}






void SetPopulation(long lPopulation)
{
    fmk::set_population_value = lPopulation;
}






long GetPopulation()
{
    return fmk::set_population_value;
}




void signal_exit_simulation_all()
{
    fmk::do_exit_simulation_all = true;
}












double PieceLinearLookup(double x, const double *ax, const double *ay, int n)
{
    double y = 0.0;

    if (x <= ax[0]) {
        
        y = ay[0];
    }
    else {
        
        bool found = false;
        for (int j = 1; j < n; ++j) {
            if (!(ax[j] > ax[j - 1])) {
                
                ModelExit("error : non-increasing x in PieceLinearLookup");
                
            }
            if (x < ax[j]) {
                found = true;
                
                y = ay[j - 1] + (ay[j] - ay[j - 1]) / (ax[j] - ax[j - 1]) * (x - ax[j - 1]);
                break;
            }
        }
        if (!found) {
            
            y = ay[n - 1];
        }
    }
    return y;
}










double PieceLinearLookup(double x, const double *axy, int n)
{
    assert(0 == n % 2);
    
    n = n / 2;
    double y = 0.0;

    if (x <= axy[0]) {
        
        y = axy[1];
    }
    else {
        
        bool found = false;
        for (int j = 1; j < n; ++j) {
            int jx = 2 * j;    
            int jy = jx + 1;   
            int jx_p = jx - 2; 
            int jy_p = jy - 2; 
            if (!(axy[jx] > axy[jx_p])) {
                
                ModelExit("error : non-increasing x in PieceLinearLookup");
                
            }
            if (x < axy[jx]) {
                found = true;
                
                y = axy[jy_p] + (axy[jy] - axy[jy_p]) / (axy[jx] - axy[jx_p]) * (x - axy[jx_p]);
                break;
            }
        }
        if (!found) {
            
            int jy = 2 * n - 1;
            y = axy[jy];
        }
    }
    return y;
}

#line 1 "C:/Users/Martin/Dropbox/openM++/use/common_modgen.ompp"










namespace fmk {

} 
#line 1 "C:/Users/Martin/Dropbox/openM++/use/random/random_lcg41.ompp"

















namespace fmk {

    


    const int max_stream_generators = 41;

    


    const long model_stream_generators[max_stream_generators] = {
        16807,
        1826645050,
        519701294,
        1912518406,
        87921397,
        755482893,
        673205363,
        727452832,
        630360016,
        1142281875,
        219667202,
        200558872,
        1185331463,
        573186566,
        396907481,
        1106264918,
        1605529283,
        1902548864,
        1444095898,
        1600915560,
        1987505485,
        1323051066,
        1715488211,
        1289290241,
        967740346,
        1644645313,
        2142074246,
        1397488348,
        97473033,
        1210640156,
        990191797,
        640039787,
        1141672104,
        2081478048,
        1236995837,
        1985494258,
        84845685,
        184528125,
        1303680654,
        61496220,
        1096609123,
    };

    
    thread_local long stream_generator = 0;

    
    thread_local long stream_seeds[model_streams];

    
    thread_local bool other_normal_valid[model_streams] = { false };

    
    thread_local double other_normal[model_streams];

} 


void new_streams()
{
	
}


void delete_streams()
{
	
}


void initialize_stream(int model_stream, int member, long seed)
{
	
	
	fmk::stream_generator = fmk::model_stream_generators[member % fmk::max_stream_generators];

    fmk::stream_seeds[model_stream] = seed;
    fmk::other_normal_valid[model_stream] = false;
}










random_state serialize_random_state()
{
    random_state rs;
    rs.reserve(2 * fmk::model_streams);
    const size_t bufsize = 50;
    char wrk[bufsize];
    for (int j = 0; j < fmk::model_streams; ++j) {
        snprintf(wrk, bufsize, "%ld", fmk::stream_seeds[j]);
        rs.push_back(wrk);
        if (fmk::other_normal_valid[j]) {
            snprintf(wrk, bufsize, "%a", fmk::other_normal[j]);
        }
        else {
            
            wrk[0] = '\0';
        }
        rs.push_back(wrk);
    }
    return rs;
}






void deserialize_random_state(const random_state & rs)
{
    for (int j = 0, k = 0; j < fmk::model_streams; ++j) {
        fmk::stream_seeds[j] = atol(rs[k++].c_str());
        auto wrk = rs[k++];
        fmk::other_normal_valid[j] = (wrk != "");
        if (fmk::other_normal_valid[j]) {
            fmk::other_normal[j] = atof(wrk.c_str());
        }
        else {
            fmk::other_normal[j] = 0.0; 
        }
    }
}

double RandUniform(int strm)
{
    using namespace fmk;

    assert(strm < model_streams); 
    if (strm > model_streams) {
        
        handle_streams_exceeded(strm, model_streams);
        
    }

    long seed = stream_seeds[strm];
    long long product = stream_generator;
    product *= seed;
    seed = product % lcg_modulus;
    stream_seeds[strm] = seed;
    return (double)seed / (double)lcg_modulus;
}

double RandNormal(int strm)
{
    using namespace fmk;

    assert(strm < model_streams); 
    if (strm > model_streams) {
        
        handle_streams_exceeded(strm, model_streams);
        
    }

    if (other_normal_valid[strm]) {
        other_normal_valid[strm] = false;
        return other_normal[strm];
    }
    else {
        double r2 = 1;
        double x = 0;
        double y = 0;
        while (r2 >= 1) {
            x = 2.0 * RandUniform(strm) - 1.0;
            y = 2.0 * RandUniform(strm) - 1.0;
            r2 = x * x + y * y;
        }
        double scale = sqrt(-2.0 * log(r2) / r2);
        double n1 = scale * x;
        double n2 = scale * y;
        other_normal[strm] = n2;
        other_normal_valid[strm] = true;
        return n1;
    }
}

double RandLogistic(int strm)
{
    using namespace fmk;

    assert(strm < model_streams); 
    if (strm > model_streams) {
        
        handle_streams_exceeded(strm, model_streams);
        
    }

    double p = RandUniform(strm);
    double odds_ratio = p / (1.0 - p);
    double x = log(odds_ratio);
    return x;
}


#line 1 "C:/Users/Martin/Dropbox/openM++/use/time_based/time_based_common.ompp"
















#line 1 "C:/Users/Martin/Dropbox/openM++/use/time_based/time_based_core.ompp"
















namespace fmk {

    


    Time simulation_end;

    


    thread_local std::chrono::system_clock::time_point clock_time_start_events;

} 





void report_simulation_progress(int member, int percent, Time i_time)
{
    theLog->logFormatted("member=%d Simulation progress=%d%% time=%g", member, percent, i_time);
    report_simulation_progress_beat(percent, (double)i_time);
}




void before_presimulation(int mem_id, int mem_count)
{
    
    fmk::simulation_members = mem_count;

    
    fmk::simulation_member = mem_id;

    
    
    
    
    fmk::master_seed = SimulationSeed;

    
    
    new_streams();

    
    initialize_model_streams(); 
}




void after_presimulation()
{
    fmk::master_seed = 0;

    
    
    delete_streams();
}

void LogSimulationStart()
{
    
    BaseEvent::set_global_time(-time_infinite);
    theLog->logFormatted("Member=%d Create starting population", fmk::simulation_member);
}

void SimulateEvents()
{
    using namespace fmk;

    theLog->logFormatted("Member=%d Simulate events", fmk::simulation_member);

    
    
    
    fmk::clock_time_start_events = std::chrono::system_clock::now();

    bool is_first_event = true;
    Time time_first_event(0); 

    
    
    
    
    
    
    bool is_step_progress = i_model->runOptions()->progressStep > 1.0;

    if (i_model->runOptions()->progressStep != 0.0  && i_model->runOptions()->progressStep < FLT_MIN) {
        is_step_progress = false;
        theLog->logFormatted("Warning: incorrect value of progress step reporting: %g", i_model->runOptions()->progressStep);
    }
    double step_progress = is_step_progress ? i_model->runOptions()->progressStep : 0.0;

    
    
    bool is_percent_progress = i_model->runOptions()->progressPercent > 0;

    if (i_model->runOptions()->progressPercent < 0) {
        theLog->logFormatted("Warning: incorrect value of progress percent reporting: %d", i_model->runOptions()->progressPercent);
    }
    int percent_progress = is_percent_progress ? i_model->runOptions()->progressPercent : progress_percent_default;

    is_percent_progress |= !is_step_progress;   

    
    double next_step_progress = step_progress;
    int next_percent_progress = percent_progress;
    bool is_100_percent_done = false;
    int64_t next_progress_beat = 0;
    int64_t next_ms_progress_beat = getMilliseconds() + OM_STATE_BEAT_TIME;
    
    report_simulation_progress(simulation_member, 0);   
    
    
    while ( true ) {
        
        Time the_time = BaseEvent::time_next_event();

        if (is_first_event) {
            
            time_first_event = the_time;
            is_first_event = false;
        }

        if (the_time > simulation_end) {
            
            BaseAgent::age_all_agents(simulation_end);
            break;
        }

        
        if (!BaseEvent::do_next_event()) {
            
            break;
        }

        
        if (simulation_end - time_first_event > 0) {

            bool is_do_percent_progress = false;
            bool is_do_step_progress = false;

            double td = the_time - time_first_event;
            int percent_done = (int)(100.0 * (td / (simulation_end - time_first_event)));

            if (is_percent_progress) {
                is_do_percent_progress = percent_done >= next_percent_progress;
                if (is_do_percent_progress) {
                    next_percent_progress = (percent_done / percent_progress) * percent_progress + percent_progress;
                }
            }
            if (!is_do_percent_progress && is_step_progress) {
                is_do_step_progress = td >= next_step_progress;
                if (is_do_step_progress) {
                    next_step_progress = floor(td / step_progress) * step_progress + step_progress;
                }
            }
            if (is_do_percent_progress || is_do_step_progress) {
                is_100_percent_done = percent_done >= 100;
                report_simulation_progress(simulation_member, percent_done, the_time);
                next_progress_beat = 0;
                next_ms_progress_beat = getMilliseconds() + OM_STATE_BEAT_TIME;
            }
            else {
                if (++next_progress_beat > 1000) {
                    next_progress_beat = 0;
                    int64_t ms = getMilliseconds();
                    if (ms > next_ms_progress_beat) {
                        report_simulation_progress_beat(percent_done, (double)the_time);
                        next_ms_progress_beat = ms + OM_STATE_BEAT_TIME;
                    }
                }
            }
        }
    }

    
    BaseAgent::exit_simulation_all();
    BaseEvent::clean_all();
    BaseAgent::free_all_zombies();

    
    if (!is_100_percent_done) {
        report_simulation_progress(simulation_member, 100, simulation_end);
    }
}







void RunSimulation(int mem_id, int mem_count, IModel * const i_model)
{
    using namespace fmk;

    auto clock_time_start = std::chrono::system_clock::now();

    
    fmk::i_model = i_model;

    
    member_entity_counter = 0;

    
    simulation_members = mem_count;

    
    simulation_member = mem_id;

    simulation_end = SimulationEnd; 

    
    
    
    
    
    fmk::master_seed = SimulationSeed;

    
    fmk::combined_seed = fmk::master_seed + simulation_member * ((long long)lcg_modulus + 1);

	
	
	new_streams();

    initialize_model_streams(); 

    
    BaseEvent::set_global_time(-time_infinite);
     
    
    Simulation();

	
	
	delete_streams();

    {
        
        auto clock_time_end = std::chrono::system_clock::now();
        std::chrono::duration<double> startpop_seconds = fmk::clock_time_start_events - clock_time_start;
        double clock_time_startpop = (double)startpop_seconds.count();
        std::chrono::duration<double> events_seconds = clock_time_end - fmk::clock_time_start_events;
        double clock_time_events = (double)events_seconds.count();

        auto event_count = BaseEvent::global_event_counter;
        double events_per_entity = (double)event_count / (double)fmk::member_entity_counter;
        theLog->logFormatted(
            "member=%d Simulation summary: entities=%ld, events/entity=%.1f, elapsed(startpop)=%.6fs, elapsed(events)=%.6fs",
            fmk::simulation_member,
            fmk::member_entity_counter,
            events_per_entity,
            clock_time_startpop,
            clock_time_events
        );
    }

}


#line 1 "C:/Users/Martin/Dropbox/openM++/use/time_based/time_based_modgen.ompp"








double GetCaseSeed()
{
    
    
    
    return 0;
}








int GetReplicate()
{
    return fmk::simulation_member;
}








int GetUserTableReplicate()
{
    return fmk::simulation_member;
}








int GetReplicates()
{
    return fmk::simulation_members;
}

double SIMULATION_END()
{
    return SimulationEnd;
}

void Set_actor_weight(double weight)
{
    
    
    
}

void Set_actor_subsample_weight(double weight)
{
    
    
    
}
#line 1 "../parameters/Default/Framework.odat"






#line 1 "../parameters/Default/PersonCore.dat"































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































#line 1 "../parameters/Default/scenario_info.odat"





