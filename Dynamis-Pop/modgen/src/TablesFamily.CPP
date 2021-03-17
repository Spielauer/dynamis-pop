#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/TablesFamily.mpp"
//LABEL (TablesFamily, EN) Family Tables
/* NOTE(TablesFamily, EN) 
Family-related table output
*/

//table_group TG02_Family                                         //EN Family
//{
//    tabHavingSpouse
//
//};

//table Person tabHavingSpouse                                    //EN Partnerships                    
//[WITHIN(SIM_YEAR_RANGE,calendar_year) && is_resident]
//{
//    sex + *                                                     //EN Sex
//    {   
//        duration(has_spouse,TRUE) / duration(),                 //EN Currently married decimals=4
//        duration(ever_union,TRUE) / duration()                  //EN Ever married decimals=4
//    }
//    * integer_age                                               //EN Age
//    * sim_year                                                  //EN Year
//};
} // namespace 
