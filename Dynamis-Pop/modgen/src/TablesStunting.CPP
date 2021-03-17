#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/TablesStunting.mpp"
//LABEL (TablesStunting, EN) Tables Stunting
/* NOTE(TablesStunting, EN) 
Tables about stunting of national born 
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Stunting tables
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//table_group TG15_Stunting                   //EN Stunting
//{
//    tabStuntingSexRegYob,
//    tabStuntingSexRegMotherYob
//};

//table Person tabStuntingSexRegYob            //EN Stunting by region
//[in_projected_time && is_resident && creation_type == CT_BIRTH && trigger_entrances(is_ready_for_birthtables,TRUE)]
//{
//    {
//        value_in(is_stunted) / unit         //EN Rate stunted decimals=3
//    }
//    *region_nat +
//    *sim_year
//};

//table Person tabStuntingSexRegMotherYob      //EN Stunting detailed 
//[in_projected_time && is_resident && creation_type == CT_BIRTH && trigger_entrances(is_ready_for_birthtables,TRUE)]
//{
//    sex + *
//    educ_mother + *
//    {
//        value_in(is_stunted) / unit         //EN Rate stunted decimals=3
//    }
//    *region_nat +
//    *sim_year
//};
} // namespace 
