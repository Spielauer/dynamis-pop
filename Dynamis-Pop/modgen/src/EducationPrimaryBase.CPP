#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/EducationPrimaryBase.mpp"
//LABEL (EducationPrimaryBase, EN) Primary Education Base Version
/* NOTE(EducationPrimaryBase, EN)

TODO update text

This module sets the education fate of primary education. The education fate is decided at birth. There are three different outcomes: 

- low: never entered primary school
- medium: entered but did not finish primary school
- high: finished primary school

Probabilities to progress from low to medium and from medium to high are given by two parameters by year of birth, district of birth, and sex. The way how the education fate is decided depends on the year of birth and the person type.

- Persons from the starting population: the education outcome is taken from the starting population file if the person was born before the years of birth covered by the parameters. This means, up to two age cut-offs, the information of the starting population is ignored at all, or only used to decide the first progression (low to medium, but not medium to high which is modeled).
- Persons born in the simulation: the education fate is decided based on the parameters.
- Immigrants from scratch: for years of birth covered by parameters the education fate is modeled. For immigrants born before, the education state is sampled from foreign born residents born in the previous 12 months.
- Immigrants from pools: the education fate is modeld like for persons from the starting population, whereby the year of birth is shifted to reflect the original record in the population file. 

The function SetEduc1BaseFate() which decides the fate  is hooked to the SetAlive() Event of the PersonCore.mpp module, i.e. is called directly at birth after the actor enters the simulation and all actor sets are available. 
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor Sets
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//EN Residents of age 0 by sex and place of birth
//actor_set Person asResidentsAge0SexGeo[sex][geo_birth] filter is_alive && is_resident && integer_age==0;

//EN Residents of age 0 by sex
//actor_set Person asResidentsAge0Sex[sex] filter is_alive && integer_age==0 && is_resident;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Types
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//classification EDUC_ONE_LEVEL                                   //EN Primary Education level
//{
//    EOL_LOW,                                                    //EN Never entered primary school
//    EOL_MEDIUM,                                                 //EN Entered primary school
//    EOL_HIGH                                                    //EN Graduated from primary school
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters
//{
//    double  EducTrans1[SEX][YOB_EDUC_TRANS1][GEO];              //EN Probability to start primary school
//    double  EducTrans2[SEX][YOB_EDUC_TRANS2][GEO];              //EN Probability to graduate from primary
//};

//parameter_group PG10_SchoolOneFateBase                          //EN Primary School Fate - Base Model
//{
//    EducTrans1,
//    EducTrans2
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor States and Events
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    EDUC_ONE_LEVEL  educ_one_fate = { EOL_LOW };        //EN Primary Education Fate
//
//    void SetEduc1BaseFate();                            //EN Setting the primary school fate
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Person::SetEduc1BaseFate()
{
    // Person born in simulation
    if (  creation_type == CT_BIRTH )
    {
        EDUC_ONE_LEVEL eolFate = EOL_LOW;
        if ( RandUniform(41) < gprParam->EducTrans1[CHKBND(sex ,2,"EducTrans1",0)][CHKBND(RANGE_POS(YOB_EDUC_TRANS1,year_of_birth ),63,"EducTrans1",1)][CHKBND(geo_birth ,13,"EducTrans1",2)] ) eolFate = EOL_MEDIUM;
        if ( eolFate == EOL_MEDIUM 
            && RandUniform(42) < gprParam->EducTrans2[CHKBND(sex ,2,"EducTrans2",0)][CHKBND(RANGE_POS(YOB_EDUC_TRANS2,year_of_birth ),67,"EducTrans2",1)][CHKBND(geo_birth ,13,"EducTrans2",2)] )  eolFate = EOL_HIGH;
         poStateFunctions->Set_educ_one_fate(  eolFate );
    }

    // Person from starting population incl. immigration pools
    else if (  creation_type == CT_START ||  creation_type == CT_POOL )
    {
        // Look up edcuation in starting population file record
        EDUC_ONE_LEVEL eolStart = EOL_LOW;
        if ( int(lStartValues ->StartPopValue[PMC_EDUC]) == 1 ) eolStart = EOL_MEDIUM;
        else if ( int(lStartValues ->StartPopValue[PMC_EDUC]) == 2 ) eolStart = EOL_HIGH;

        // Look up relevant year of birth from starting population file (important for those in pool)
        int nYOB = int(lStartValues ->StartPopValue[PMC_BIRTH]);

        // Case 1: take the education from the starting population
        if ( nYOB < MIN(YOB_EDUC_TRANS2))  poStateFunctions->Set_educ_one_fate(  eolStart );
        
        // Case 2: take first progression from the starting population
        else if ( nYOB < MIN(YOB_EDUC_TRANS1))
        {
            EDUC_ONE_LEVEL eolFate = EOL_LOW;
            if ( eolStart != EOL_LOW ) eolFate = EOL_MEDIUM;
            if ( eolFate == EOL_MEDIUM 
                && RandUniform(43) < gprParam->EducTrans2[CHKBND(sex ,2,"EducTrans2",0)][CHKBND(RANGE_POS(YOB_EDUC_TRANS2,nYOB),67,"EducTrans2",1)][CHKBND(geo_birth ,13,"EducTrans2",2)] )  eolFate = EOL_HIGH;
             poStateFunctions->Set_educ_one_fate(  eolFate );
        }
        // Case 3: model from parameters ignoring starting population values
        else
        {
            EDUC_ONE_LEVEL eolFate = EOL_LOW;
            if ( RandUniform(44) < gprParam->EducTrans1[CHKBND(sex ,2,"EducTrans1",0)][CHKBND(RANGE_POS(YOB_EDUC_TRANS1,nYOB),63,"EducTrans1",1)][CHKBND(geo_birth ,13,"EducTrans1",2)] ) eolFate = EOL_MEDIUM;
            if ( eolFate == EOL_MEDIUM 
                && RandUniform(45) < gprParam->EducTrans2[CHKBND(sex ,2,"EducTrans2",0)][CHKBND(RANGE_POS(YOB_EDUC_TRANS2,nYOB),67,"EducTrans2",1)][CHKBND(geo_birth ,13,"EducTrans2",2)] )  eolFate = EOL_HIGH;
             poStateFunctions->Set_educ_one_fate(  eolFate );
        }
    }

    // Person is an immigrant
    else if (creation_type ==CT_SCRATCH)
    {
        // Parameters available for this year of birth (or nobody to clone from, the use last availabel probabilities from parameter)
        if ( year_of_birth >= MIN(YOB_EDUC_TRANS1) || (asResidentsAge0SexGeo[sex ][geo_birth ]->Count() == 0 && asResidentsAge0Sex[sex ]->Count() == 0))
        {
            EDUC_ONE_LEVEL eolFate = EOL_LOW;
            if (RandUniform(46) < gprParam->EducTrans1[CHKBND(sex ,2,"EducTrans1",0)][CHKBND(RANGE_POS(YOB_EDUC_TRANS1, year_of_birth ),63,"EducTrans1",1)][CHKBND(geo_birth ,13,"EducTrans1",2)]) eolFate = EOL_MEDIUM;
            if (eolFate == EOL_MEDIUM && RandUniform(47) < gprParam->EducTrans2[CHKBND(sex ,2,"EducTrans2",0)][CHKBND(RANGE_POS(YOB_EDUC_TRANS2, year_of_birth ),67,"EducTrans2",1)][CHKBND(geo_birth ,13,"EducTrans2",2)])  eolFate = EOL_HIGH;
             poStateFunctions->Set_educ_one_fate(  eolFate );
        }
        // Education to be sampled from resident immigrants
        else if (asResidentsAge0SexGeo[sex ][geo_birth ]->Count() > 0) //there are 0 year old residents born in same foreign geo
        {
             poStateFunctions->Set_educ_one_fate(  asResidentsAge0SexGeo[sex ][geo_birth ]->GetRandom(RandUniform(48))->educ_one_fate  );
        }
        // else sample from resident population
        else  poStateFunctions->Set_educ_one_fate(  asResidentsAge0Sex[sex ]->GetRandom(RandUniform(49))->educ_one_fate  ); 
    }
}

} // namespace 
