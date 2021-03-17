#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/_ContextABC.mpp"
//LABEL (_ContextABC, EN)
/*NOTE (_ContextABC, EN) 
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Starting population 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//classification PERSON_MICRODATA_COLUMNS             //EN Starting population variables
//{
//    PMC_FAMID,                                      //EN Family ID
//    PMC_WEIGHT,                                     //EN Weight
//    PMC_ROLE,                                       //EN Role in family
//    PMC_BIRTH,                                      //EN Time of birth
//    PMC_SEX,                                        //EN Sex
//    PMC_EDUC,                                       //EN Primary education status
//    PMC_SECGRAD,                                    //EN Highest secondary grade attended
//    PMC_SECATTEND,                                  //EN Currently attending secondary education
//    PMC_UNION,                                      //EN Time of first union formation
//    PMC_PARITY,                                     //EN Parity
//    PMC_LASTBIR,                                    //EN Time of last birth
//    PMC_ETHNO,                                      //EN Ethnicity
//    PMC_GEO,                                        //EN Place of residence
//    PMC_GEOBIR,                                     //EN Place of residence at birth
//    PMC_GEOPRE,                                     //EN Pervious place of residence
//    PMC_MOVEFIRST,                                  //EN Time of first move
//    PMC_MOVELAST,                                   //EN Time of last move
//    PMC_POOL                                        //EN Immigration Pool
//};

//actor Person
//{
//    //EN In selected year
//    logical in_sel_years = calendar_year == 2000 || calendar_year == 2010 
//        || calendar_year == 2020 || calendar_year == 2030 || calendar_year == 2040 
//        || calendar_year == 2050 || calendar_year == 2060;
//
//    //EN Year
//    TAB_SEL_YEARS sel_years = (calendar_year == 2000) ? TSY_2000 : 
//        (calendar_year == 2010) ? TSY_2010 :
//        (calendar_year == 2020) ? TSY_2020 :
//        (calendar_year == 2030) ? TSY_2030 :
//        (calendar_year == 2040) ? TSY_2040 :   
//        (calendar_year == 2050) ? TSY_2050 : TSY_2060;
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// New Dimensions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//range ALL_YEAR_RANGE{ 1899, 2150 };                 //EN All calendar years
//range SIM_YEAR_RANGE{ 2000, 2150 };                 //EN Simulated years
//range YOB_UNION{ 1951, 2037 };                      //EN Year of birth
//range YOB_EDUC_TRANS1 { 1988, 2050 };               //EN Year of birth
//range YOB_EDUC_TRANS2 { 1984, 2050 };               //EN Year of birth
//range CHILD_MORTALITY_YEARS{ 2005, 2150 };          //EN Year range for which child mortality is modeled


//partition IMMU_YOB_PART{ 1971,1976,1981 };          //EN Birth Cohort

//classification VACCINATION_PARA                     //EN Immunization
//{
//    IP_CONSTANT,    //EN Constant
//    IP_MALE,        //EN Male 
//    IP_EDUCMO_1,    //EN Mother education 'medium' 
//    IP_EDUCMO_2,    //EN Mother education 'high'
//    IP_REGION_1,    //EN West (Far West is reference)
//    IP_REGION_2,    //EN Central
//    IP_REGION_3,    //EN Capital
//    IP_REGION_4,    //EN East
//    IP_ETNO_1,      //EN Ethnicity B (Ethnicity A is reference)
//    IP_ETNO_2,      //EN Ethnicity C
//    IP_ETNO_3,      //EN Ethnicity D
//    IP_ETNO_4,      //EN Ethnicity E
//    IP_ETNO_5,      //EN Ethnicity F
//    IP_ETNO_6,      //EN Ethnicity G
//    IP_COHORT_1,    //EN Birth cohorts of mother 1970-1974 (<1970 is reference)
//    IP_COHORT_2,    //EN Birth cohorts of mother  1975-1979
//    IP_COHORT_3,    //EN Birth cohorts of mother  1980+
//    IP_YOUNGMO      //EN Young Mother (age < 18)
//};

//classification PRENATCARE_PARA //EN Prenatal care parameters
//{
//    PP_CONSTANT,    //EN Constant
//    PP_EDUCMO_1,    //EN Mother education 'medium' 
//    PP_EDUCMO_2,    //EN Mother education 'high'
//    PP_REGION_1,    //EN West (Far West is reference)
//    PP_REGION_2,    //EN Central
//    PP_REGION_3,    //EN Capital
//    PP_REGION_4,    //EN East
//    PP_COHORT_1,    //EN Birth cohorts of mother 1970-1974 (<1970 is reference)
//    PP_COHORT_2,    //EN Birth cohorts of mother  1975-1979
//    PP_COHORT_3,    //EN Birth cohorts of mother  1980+
//    PP_YOUNGMO      //EN Young Mother (age < 18)
//};

//classification ETHNICITY                            //EN Ethnic Group
//{
//    ETHNO_00,       //EN Ethnicity 0
//    ETHNO_01,       //EN Ethnicity 1
//    ETHNO_02,       //EN Ethnicity 2
//    ETHNO_03,       //EN Ethnicity 3
//    ETHNO_04,       //EN Ethnicity 4
//    ETHNO_05,       //EN Ethnicity 5
//    ETHNO_06,       //EN Ethnicity 6
//    ETHNO_07,       //EN Ethnicity 7
//    ETHNO_08,       //EN Ethnicity 8
//    ETHNO_09        //EN Ethnicity 9
//};

//classification ETHNICITY_SHORT                      //EN Ethnic Group Aggregated
//{
//    ES_00,          //EN Ethnicity A
//    ES_01,          //EN Ethnicity B
//    ES_02,          //EN Ethnicity C
//    ES_03,          //EN Ethnicity D
//    ES_04,          //EN Ethnicity E
//    ES_05,          //EN Ethnicity F
//    ES_06           //EN Ethnicity G
//};

//EN Aggregation Ethnicity to simplified Ethnicity
//aggregation ETHNICITY_SHORT, ETHNICITY
//{
//    ES_00, ETHNO_00,         // Ethnicity A
//    ES_00, ETHNO_03,         // Ethnicity A
//    ES_01, ETHNO_04,         // Ethnicity B
//    ES_02, ETHNO_02,         // Ethnicity C
//    ES_03, ETHNO_05,         // Ethnicity D
//    ES_04, ETHNO_01,         // Ethnicity E
//    ES_04, ETHNO_06,         // Ethnicity E
//    ES_05, ETHNO_07,         // Ethnicity F
//    ES_06, ETHNO_08,         // Ethnicity G
//    ES_06, ETHNO_09          // Ethnicity G
//};

//classification GEO            //EN Smallest Geographical Units
//{
//    GEO_00,                   //EN Far-West-Mountains
//    GEO_01,                   //EN Far-West-Hills
//    GEO_02,                   //EN Far-West-Lowlands
//    GEO_03,                   //EN West-Mountains
//    GEO_04,                   //EN West-Hills
//    GEO_05,                   //EN West-Lowlands
//    GEO_06,                   //EN Central-Mountains
//    GEO_07,                   //EN Central-Hills
//    GEO_08,                   //EN Central-Lowlands
//    GEO_09,                   //EN East-Mountains
//    GEO_10,                   //EN Capital
//    GEO_11,                   //EN East-Lowlands
//    GEO_ABROAD                //EN Abroad
//};

//classification GEO_ABR        //EN Foreign Geo
//{
//    GEOA_00                   //EN Abroad  
//};

//classification GEO_NAT        //EN Districts
//{
//    GEON_00,                   //EN Far-West-Mountains
//    GEON_01,                   //EN Far-West-Hills
//    GEON_02,                   //EN Far-West-Lowlands
//    GEON_03,                   //EN West-Mountains
//    GEON_04,                   //EN West-Hills
//    GEON_05,                   //EN West-Lowlands
//    GEON_06,                   //EN Central-Mountains
//    GEON_07,                   //EN Central-Hills
//    GEON_08,                   //EN Central-Lowlands
//    GEON_09,                   //EN East-Mountains
//    GEON_10,                   //EN Capital
//    GEON_11                    //EN East-Lowlnds
//};

//classification REGION_NAT       //EN Region
//{
//    REGN_00,                    //EN Far-West
//    REGN_01,                    //EN West
//    REGN_02,                    //EN Central
//    REGN_03,                    //EN Capital
//    REGN_04                     //EN East
//};

//classification REGION_INT       //EN Region
//{
//    REGI_00,                    //EN Far-West
//    REGI_01,                    //EN West
//    REGI_02,                    //EN Central
//    REGI_03,                    //EN Capital
//    REGI_04,                    //EN East
//    REGI_ABROAD                 //EN Abroad
//};


//EN Aggregation Region International to National (use with filter is_resident only )
//aggregation REGION_NAT, REGION_INT
//{
//    REGN_00, REGI_00,
//    REGN_01, REGI_01,
//    REGN_02, REGI_02,
//    REGN_03, REGI_03,
//    REGN_04, REGI_04,
//    REGN_04, REGI_ABROAD
//};

//EN Aggregation District to Region
//aggregation REGION_INT, GEO
//{
//    REGI_00, GEO_00,
//    REGI_00, GEO_01,
//    REGI_00, GEO_02,
//    REGI_01, GEO_03,
//    REGI_01, GEO_04,
//    REGI_01, GEO_05,
//    REGI_02, GEO_06,
//    REGI_02, GEO_07,
//    REGI_02, GEO_08,
//    REGI_04, GEO_09,
//    REGI_03, GEO_10,
//    REGI_04, GEO_11,
//    REGI_ABROAD, GEO_ABROAD
//};

//classification BIRTH1_LOC                           //EN Geographical Location
//{
//    B1L_00,                                         //EN Far-West
//    B1L_01,                                         //EN West
//    B1L_02,                                         //EN Central
//    B1L_03,                                         //EN Capital
//    B1L_04                                          //EN East
//};

//EN Aggregation Location (used in) First Birth to Districts
//aggregation BIRTH1_LOC, GEO
//{
//    B1L_00, GEO_00,
//    B1L_00, GEO_01,
//    B1L_00, GEO_02,
//    B1L_01, GEO_03,
//    B1L_01, GEO_04,
//    B1L_01, GEO_05,
//    B1L_02, GEO_06,
//    B1L_02, GEO_07,
//    B1L_02, GEO_08,
//    B1L_04, GEO_09,
//    B1L_03, GEO_10,
//    B1L_04, GEO_11,
//    B1L_04, GEO_ABROAD
//};

//classification POP_POOL                             //EN Population Pool
//{
//    PP_NON,                                         //EN Not a pool
//    PP_STUDENTS,                                    //EN Students
//    PP_WORKER,                                      //EN Worker
//    PP_REFUGEES                                     //EN Refugees
//};

//classification IMMI_POOL                            //EN Immigration Pool
//{
//    IP_STUDENTS,                                    //EN Students
//    IP_WORKER,                                      //EN Worker
//    IP_REFUGEES                                     //EN Refugees
//};

//classification  TAB_SEL_YEARS           //EN Selected Years
//{
//    TSY_2000,                           //EN 2000    
//    TSY_2010,                           //EN 2010
//    TSY_2020,                           //EN 2020
//    TSY_2030,                           //EN 2030
//    TSY_2040,                           //EN 2040
//    TSY_2050,                           //EN 2050
//    TSY_2060                            //EN 2060
//};

//EN Aggregation District International to National (use with filter is_resident only)
//aggregation GEO_NAT, GEO
//{
//    GEON_00, GEO_00,
//    GEON_01, GEO_01,
//    GEON_02, GEO_02,
//    GEON_03, GEO_03,
//    GEON_04, GEO_04,
//    GEON_05, GEO_05,
//    GEON_06, GEO_06,
//    GEON_07, GEO_07,
//    GEON_08, GEO_08,
//    GEON_09, GEO_09,
//    GEON_10, GEO_10,
//    GEON_11, GEO_11,
//    GEON_11, GEO_ABROAD
//};


//range FEMALE_SPOUSE_AGE { 10, 60 };                 //EN Age
//range MALE_SPOUSE_AGE { 15, 100 };                  //EN Age 

//range AGE_RANGE { 0, 100 };                         //EN Age range
//range FERTILE_AGE_RANGE { 12, 49 };                 //EN Fertile age range















} // namespace 
