#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/EmigrationBase.mpp"
//LABEL(EmigrationBase, EN) Emigration General

/* NOTE(EmigrationBase, EN)
This module implements emigration
*/

                                                     
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// parameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters
//{
//    //EN Switch emigration on/off
//    logical	ModelEmigration;
//
//    //EN Emigration rates on district level
//    double  EmigrationRatesDistrict[SEX][AGE_MIG][GEO_NAT];
//
//    //EN Destination
//    cumrate EmigrationDestination[GEO_ABR];
//};

//parameter_group PG03_Emigration                             //EN Emigration
//{
//    ModelEmigration,
//    EmigrationRatesDistrict, 
//    EmigrationDestination
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    event   timeEmigrationEvent, EmigrationEvent;           //EN Emigration event
//};

//actor Calibrator
//{
//    //EN Idenitify and return next person to emigrate
//    Person_ptr GetNextToEmigrate(Person_ptr prPerson, GEO toGeo);   
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TIME Person::timeEmigrationEvent()
{
    if (gprParam->ModelEmigration &&  is_resident &&  is_mortal && gprParam->EmigrationRatesDistrict[CHKBND(sex ,2,"EmigrationRatesDistrict",0)][CHKBND(age_mig ,13,"EmigrationRatesDistrict",1)][CHKBND((GEO_NAT)geo ,12,"EmigrationRatesDistrict",2)] > 0.0)
    {
        return WAIT(-log(RandUniform(38)) / gprParam->EmigrationRatesDistrict[CHKBND(sex ,2,"EmigrationRatesDistrict",0)][CHKBND(age_mig ,13,"EmigrationRatesDistrict",1)][CHKBND((GEO_NAT)geo ,12,"EmigrationRatesDistrict",2)]);
    }
    else return TIME_INFINITE;
}

void Person::EmigrationEvent()
{
    // Sample the destination
    int nDestination; Lookup_EmigrationDestination(RandUniform(56), &nDestination);
    
    // move
    lCalibrator ->UpdateTime()->GetNextToEmigrate(this,(GEO)(SIZE(GEO_NAT) + nDestination))->UpdateTime()->doResidentialMove((GEO)(SIZE(GEO_NAT) + nDestination));
}


Person_ptr Calibrator::GetNextToEmigrate(Person_ptr ptrPerson, GEO cGeoTo)
{
    if (gprParam->MigrationTryKeepingFamiliesTogether)
    {
        if (asWantToMove[ptrPerson->geo ][(GEO)cGeoTo][ptrPerson->sex ][ptrPerson->age_mig ]->Count() > 0)
        {
            return asWantToMove[ptrPerson->geo ][(GEO)cGeoTo][ptrPerson->sex ][ptrPerson->age_mig ]->GetRandom(RandUniform(27));
        }
        else return ptrPerson;
    }
    else return ptrPerson;
}


} // namespace 
