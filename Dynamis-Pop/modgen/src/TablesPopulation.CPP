#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/TablesPopulation.mpp"
//LABEL(TablesPopulation, EN) Population Tables Output

/* NOTE(TablesPopulation, EN)
This module contains table output focusing on general population dnamiscs.
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Popultion tables
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//partition tab_age15 { 15,30,60 };                                       //EN Age group
//classification SELECTED_YEAR10                                          //EN Selected years
//{
//    SY00,                                                               //EN Start Year
//    SY10,                                                               //EN Start Year + 10
//    SY20,                                                               //EN Start Year + 20
//    SY30,                                                               //EN Start Year + 30
//    SY40                                                                //EN Start Year + 40
//};

//actor Person
//{
//    ALL_YEAR_RANGE all_year = COERCE(ALL_YEAR_RANGE, calendar_year);    //EN Year
//    SIM_YEAR_RANGE sim_year = COERCE(SIM_YEAR_RANGE, calendar_year);    //EN Year
//
//    //EN Selected year
//    SELECTED_YEAR10 selected_year10 = 
//        (calendar_year == MIN(SIM_YEAR_RANGE)) ? SY00 :
//        (calendar_year == MIN(SIM_YEAR_RANGE) + 10) ? SY10 :
//        (calendar_year == MIN(SIM_YEAR_RANGE) + 20) ? SY20 :
//        (calendar_year == MIN(SIM_YEAR_RANGE) + 30) ? SY30 : SY40;
//
//    //EN In selected year
//    logical in_selected_year = (selected_year10 != SY40 || calendar_year == MIN(SIM_YEAR_RANGE) + 40);
//};

//table_group TG00_Population                                     //EN Population
//{
//    tabPopulationYearPlace,
//    tabBirthsYearPlace,
//    tabMigrationOriginDestination
//};

//table Person tabPopulationYearPlace                             //EN Population by year and location
//[ WITHIN(SIM_YEAR_RANGE,calendar_year) ]
//{ 
//    {
//        duration()                                              //EN Population     
//    }
//    * geo+                                                      //EN Place
//    * sim_year                                                  //EN Year
//};

//table Person tabBirthsYearPlace                                 //EN Births by year and location
//[ WITHIN(SIM_YEAR_RANGE,calendar_year) ]
//{ 
//    {
//        transitions(is_alive,FALSE,TRUE)                        //EN Births
//    }
//    * geo+                                                      //EN Place
//    * sim_year                                                  //EN Year
//};

//table Person tabMigrationOriginDestination                      //EN Migration Origin Destination selected years
//[ WITHIN(SIM_YEAR_RANGE,calendar_year) 
//  && trigger_changes(geo) && in_selected_year ]
//{ 
//    selected_year10 *                                           //EN Year
//    {
//        value_in(is_alive)                                      //EN Moves
//    }
//    * geo_prev+                                                 //EN From
//    * geo                                                       //EN To
//};






} // namespace 
