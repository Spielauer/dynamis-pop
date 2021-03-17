#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/MigrationBase.mpp"
//LABEL(MigrationBase, EN) Internal Migration Base Version

/* NOTE(MigrationBase, EN)
This module implements internal migration. It is based on age and sex specific transition matrices, which is a typical approach in cohort-component models. As the transitions are typically obtained from Census data recording the place of residence now and one year ago, only one transition per year is allowed. The module has three parameters, one to switch migration on/off, one containing the probabilities to move away (by origin, age, sex), and one to sample the destination (by origin, age and sex).
*/


////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters
//{
//    //EN Internal migration on/off
//    logical ModelMigration;
//
//    //EN Migration probability
//    double  MigrationProbability[SEX][AGE_MIG][GEO_NAT];
//
//    //EN Migration Destination
//    cumrate MigrationDestination[SEX][GEO_NAT][AGE_MIG][GEO_NAT];
//};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor block
////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    event           timeMigrationEvent, MigrationEvent;         //EN Migration Event
//};

//actor Calibrator
//{
//   Person_ptr GetNextToMigrate(Person_ptr prPerson, GEO_NAT cGeoTo);   
//};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////////////////////////

TIME Person::timeMigrationEvent()
{
    // Check if a person is at risk for moving
    if (gprParam->ModelMigration &&  is_mortal &&  is_resident &&  calendar_year != (int)time_last_move )
    {
        // get the probability to move
        double dMoveProb = gprParam->MigrationProbability[CHKBND(sex ,2,"MigrationProbability",0)][CHKBND(age_mig ,13,"MigrationProbability",1)][CHKBND((GEO_NAT)geo ,12,"MigrationProbability",2)];

        if (dMoveProb <= 0.0) return TIME_INFINITE;     // Not at risk
        else if (dMoveProb >= 1.0) return WAIT(0);      // 100% move immediately
        else                                            // positive probability below 100%
        {
            // calculate a random waiting time based on the given probability converted to
            // a hazard rate = -log(1-probability)
            return WAIT(-log(RandUniform(60)) / -log(1 - dMoveProb));
        }
    }
    return TIME_INFINITE;
}

void Person::MigrationEvent()
{
    int nDestination;

    // Sample the destination
    Lookup_MigrationDestination(RandUniform(61), sex , (GEO_NAT)geo , int(age_mig ), &nDestination);

    //move
    lCalibrator ->UpdateTime()->GetNextToMigrate(this,(GEO_NAT)nDestination)->UpdateTime()->doResidentialMove((GEO)nDestination);
}

Person_ptr Calibrator::GetNextToMigrate(Person_ptr ptrPerson, GEO_NAT cGeoTo)
{
    if (gprParam->MigrationTryKeepingFamiliesTogether)
    {
        if (asWantToMove[ptrPerson->geo ][(GEO)cGeoTo][ptrPerson->sex ][ptrPerson->age_mig ]->Count() > 0)
        {
            return asWantToMove[ptrPerson->geo ][(GEO)cGeoTo][ptrPerson->sex ][ptrPerson->age_mig ]->GetRandom(RandUniform(25));
        }
        else return ptrPerson;
    }
    else return ptrPerson;
}




} // namespace 
