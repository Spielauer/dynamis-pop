#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/EducationPrimaryPlanning.mpp"
//LABEL (EducationPrimaryPlanning, EN)
/*NOTE (EducationPrimaryPlanning, EN) 
This module is an add-on to the primary education module used for planning of necessary school infrastructure investments. Based on parameters of current and future target teacher to student ratios and classroom sizes, it calculates the number of teachers and classrooms required for each calendar year. The module is optional and can be removed without damage to the model.
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Types
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//classification SCHOOL_INFRA             //EN School infrastructure measures
//{
//    SI_TEACHERS,                        //EN Pupil to teacher ratio
//    SI_ROOMS                            //EN Classroom sizes
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters
//{
//    //EN Pupil to teacher and classroom ratios
//    double  Educ1Infrastructure[SCHOOL_INFRA][GEO_NAT][SIM_YEAR_RANGE];
//};

//parameter_group PG15_Educ1Ressources    //EN Ressource planning
//{
//    Educ1Infrastructure
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor Person
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    //EN Required teachers
//    double  educ_one_required_teachers = (is_educ_one_student && is_resident) ? 1.0
//    / Educ1Infrastructure[SI_TEACHERS][geo_nat][RANGE_POS(SIM_YEAR_RANGE,calendar_year)] : 0.0;
//
//    //EN Required classrooms
//    double  educ_one_required_rooms = (is_educ_one_student && is_resident) ? 1.0
//    / Educ1Infrastructure[SI_ROOMS][geo_nat][RANGE_POS(SIM_YEAR_RANGE, calendar_year)] : 0.0;
//
//};


} // namespace 
