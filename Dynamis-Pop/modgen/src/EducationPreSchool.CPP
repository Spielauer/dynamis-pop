#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/EducationPreSchool.mpp"
//LABEL (EducationPreSchool, EN) Pre-School Education

/* NOTE(EducationPreSchool, EN)
The pre-school module is a simple module implementing up to 2 years of pre-school experience. This module is an add-on developed when introducing the calculations of the Human Capital Index (HCI). It is assumed, that all children attending pre-school enter primary education. Pre-school experience is decided at age 4. 
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dimensions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//classification PRE_SCHOOL_PARA      //EN Preschool
//{
//    PLP_ANY,                        //EN Proportion of primary entrants ever attending pre-school
//    PLS_TWO                         //EN Proportion of pre-school attendants attending 2 years
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters
//{
//    double PreSchoolAttendance[SEX][REGION_NAT][SIM_YEAR_RANGE][PRE_SCHOOL_PARA];   //EN Pre-School Attendance
//};

//parameter_group PG15_PreSchool                                                      //EN Pre-School
//{
//    PreSchoolAttendance
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor States and Functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    int years_preschool = { 0 };                                                    //EN Years in pre-school
//
//    logical preschool_is_decided = { FALSE };                                       //EN Pre-school fate is decided
//    event timeSetPreschoolYearsEvent, SetPreschoolYearsEvent;                       //EN Decide pre-school
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TIME Person::timeSetPreschoolYearsEvent()
{
    if (! preschool_is_decided &&  creation_type == CT_BIRTH &&  is_resident &&  educ_one_grade_attended == 1) return WAIT(0);
    else return TIME_INFINITE;
}

void Person::SetPreschoolYearsEvent()
{   
    if ( educ_one_fate != EOL_LOW)
    {
        if (RandUniform(10) < gprParam->PreSchoolAttendance[CHKBND(sex ,2,"PreSchoolAttendance",0)][CHKBND(region_nat ,5,"PreSchoolAttendance",1)][CHKBND(RANGE_POS(SIM_YEAR_RANGE, calendar_year ),151,"PreSchoolAttendance",2)][CHKBND(PLP_ANY,2,"PreSchoolAttendance",3)])
        {
            poStateFunctions->PlusPlusSet_years_preschool();
            if (RandUniform(59) < gprParam->PreSchoolAttendance[CHKBND(sex ,2,"PreSchoolAttendance",0)][CHKBND(region_nat ,5,"PreSchoolAttendance",1)][CHKBND(RANGE_POS(SIM_YEAR_RANGE, calendar_year ),151,"PreSchoolAttendance",2)][CHKBND(PLS_TWO,2,"PreSchoolAttendance",3)]) poStateFunctions->PlusPlusSet_years_preschool();
        }
    }
     poStateFunctions->Set_preschool_is_decided(  TRUE );
}

} // namespace 
