#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/TablesVaccination.mpp"
//LABEL (TablesVaccination, EN) Tables Child Immunization
/* NOTE(TablesVaccination, EN) 
A collection of tables on child vaccination (immunization)
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Child vaccination tables
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//table_group TG12_Vaccination                                //EN Child Vaccination
//{
//    TabPrenatCare,
//    TabImmunizationChildren,
//    TabChildVaccination
//};

//table Person TabPrenatCare                                  //EN Proportion of prenatal care
//[in_projected_time && integer_age == 0 && creation_type == CT_BIRTH && is_resident]
//{
//    region_nat+*                                            //EN Region
//    {
//        
//        duration(got_prenat_care,GPC_YES) / duration()      //EN Proportion prenatal care decimals=3
//    }
//    * educ_mother+                                          //EN Mother's education
//    * sim_year                                              //EN Year
//};

//table Person TabImmunizationChildren                        //EN Proportion of chidren vaccinated detailed
//[in_projected_time && integer_age == 0 && creation_type == CT_BIRTH && is_resident]
//{
//    sex+ *                                                  //EN Sex
//    ethnicity_short+*                                       //EN Ethnicity
//    got_prenat_care+*                                       //EN Mother received prenatal care    
//    region_nat+*                                            //EN Region
//    {
//        
//        duration(is_immunized,TRUE) / duration()            //EN Proportion immunized decimals=3
//    }
//    * educ_mother+                                          //EN Mother's education
//    * sim_year                                              //EN Year
//};

//table Person TabChildVaccination                            //EN Proportion of chidren vaccinated by region
//[in_projected_time && creation_type == CT_BIRTH && is_resident && integer_age == 0]
//{
//    {
//        duration(is_immunized,TRUE) / duration() 	        //EN Vaccination rate decimals=3
//    }
//    * sim_year                                              //EN Year
//    * region_nat+                                           //EN Region
//};


} // namespace 
