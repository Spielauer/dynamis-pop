#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/Stunting.mpp"
//LABEL (Stunting, EN) Stunting
/* NOTE(Stunting, EN)
Stunting is the impaired growth and development that children experience from poor nutrition, repeated infection, and inadequate psychosocial stimulation. Children are defined as stunted if their height-for-age is more than two standard deviations below the WHO Child Growth Standards median. Stunting in early life has adverse functional consequences on the child. Some of those consequences include poor cognition and educational performance, low adult wages, lost productivity and, when accompanied by excessive weight gain later in childhood, an increased risk of nutrition-related chronic diseases in adult life. (WHO)

This module implements stunting as an individual level 'fate' decided at birth based on sex, region, and mother's education. Stunting is modeled for babies born within the country in simulated time.

In the current application, stunting has no time dimension thus rates are assumed to stay unchanged for given sex, region, and mother's education. All changes in aggregate stunting rates thus result from composition effects. Accordingly, the module has only one parameter: the proportion of children affected by stunting by sex, region, and mother's education.
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters
//{
//    double ProportionStunting[SEX][REGION_NAT][EDUC_ONE_LEVEL];         //EN Stunting by mother's education, region and sex
//};

//parameter_group PG_Stunting                                             //EN Stunting
//{
//    ProportionStunting
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Person actor
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    logical is_stunted = { FALSE };                                     //EN Person with impared growth (stunting) age<5
//    void    DecideStuntingFate();                                       //EN Decide stunting fate at birth
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Person::DecideStuntingFate()
{
    if ( in_projected_time &&  creation_type == CT_BIRTH && is_resident )
    {
        if (RandUniform(81) < gprParam->ProportionStunting[CHKBND(sex ,2,"ProportionStunting",0)][CHKBND(region_nat ,5,"ProportionStunting",1)][CHKBND(educ_mother ,3,"ProportionStunting",2)])
        {
             poStateFunctions->Set_is_stunted(  TRUE );
        }
    }
};

} // namespace 
