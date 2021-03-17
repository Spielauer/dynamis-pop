namespace mm {


// one-word synonyms for unsigned integral types
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

// TIME definition
typedef double TIME;
typedef TIME TIME_t;

// real definition
typedef double real;

// counter definition
typedef int counter;

// integer definition
typedef int integer;

// index definition
typedef int index;

// logical definition
typedef int logical;

// time constants
CONSTANT( TIME, ANY_AGE, (TIME) -1 );
CONSTANT( TIME, TIME_UNDEF, (TIME) -1 );
CONSTANT( TIME, TIME_INFINITE, (TIME) SHRT_MAX );

extern Thread EventQueue *gpoEventQueue;

void SetCaseWeight( double dCaseWeight, double dCaseSubsampleWeight = UNUSED_VALUE );

// classifications
enum BIRTH1_GROUP {
B1G_00,
B1G_01,
B1G_02 };
typedef BIRTH1_GROUP BIRTH1_GROUP_t;
CONSTANT( int, rc_size_BIRTH1_GROUP, 3 );
enum BIRTH1_LOC {
B1L_00,
B1L_01,
B1L_02,
B1L_03,
B1L_04 };
typedef BIRTH1_LOC BIRTH1_LOC_t;
CONSTANT( int, rc_size_BIRTH1_LOC, 5 );
enum CHILD_AGEGR {
CA00,
CA01,
CA03,
CA06,
CA09,
CA12,
CA15 };
typedef CHILD_AGEGR CHILD_AGEGR_t;
CONSTANT( int, rc_size_CHILD_AGEGR, 7 );
enum CHILD_MORTALITY_GROUP {
CMG_00,
CMG_01,
CMG_02,
CMG_03,
CMG_04,
CMG_05,
CMG_06,
CMG_07,
CMG_08 };
typedef CHILD_MORTALITY_GROUP CHILD_MORTALITY_GROUP_t;
CONSTANT( int, rc_size_CHILD_MORTALITY_GROUP, 9 );
enum CREATION_TYPE {
CT_START,
CT_POOL,
CT_BIRTH,
CT_SCRATCH };
typedef CREATION_TYPE CREATION_TYPE_t;
CONSTANT( int, rc_size_CREATION_TYPE, 4 );
enum EDUC1_GROUP {
E1G_00,
E1G_01,
E1G_02 };
typedef EDUC1_GROUP EDUC1_GROUP_t;
CONSTANT( int, rc_size_EDUC1_GROUP, 3 );
enum EDUC1_MODEL {
E1M_BASE,
E1M_REFINED_ALIGNALL,
E1M_REFINED_ALIGNONCE };
typedef EDUC1_MODEL EDUC1_MODEL_t;
CONSTANT( int, rc_size_EDUC1_MODEL, 3 );
enum EDUC_ONE_GEO {
EOG_OO };
typedef EDUC_ONE_GEO EDUC_ONE_GEO_t;
CONSTANT( int, rc_size_EDUC_ONE_GEO, 1 );
enum EDUC_ONE_GROUP {
E1G_OO };
typedef EDUC_ONE_GROUP EDUC_ONE_GROUP_t;
CONSTANT( int, rc_size_EDUC_ONE_GROUP, 1 );
enum EDUC_ONE_LEVEL {
EOL_LOW,
EOL_MEDIUM,
EOL_HIGH };
typedef EDUC_ONE_LEVEL EDUC_ONE_LEVEL_t;
CONSTANT( int, rc_size_EDUC_ONE_LEVEL, 3 );
enum EDUC_ONE_STATUS {
EOS_NEVER,
EOS_WAIT,
EOS_ATTEND,
EOS_PAUSE,
EOS_OUT };
typedef EDUC_ONE_STATUS EDUC_ONE_STATUS_t;
CONSTANT( int, rc_size_EDUC_ONE_STATUS, 5 );
enum EDUC_TWO_GEO {
ETG_OO };
typedef EDUC_TWO_GEO EDUC_TWO_GEO_t;
CONSTANT( int, rc_size_EDUC_TWO_GEO, 1 );
enum EDUC_TWO_GROUP {
E2G_OO };
typedef EDUC_TWO_GROUP EDUC_TWO_GROUP_t;
CONSTANT( int, rc_size_EDUC_TWO_GROUP, 1 );
enum EDUC_TWO_LEVEL {
ETL_NEVER,
ETL_INCOMPLETE,
ETL_GRADUATE };
typedef EDUC_TWO_LEVEL EDUC_TWO_LEVEL_t;
CONSTANT( int, rc_size_EDUC_TWO_LEVEL, 3 );
enum EDUC_TWO_STATUS {
ETS_NEVER,
ETS_WAIT,
ETS_ATTEND,
ETS_PAUSE,
ETS_OUT };
typedef EDUC_TWO_STATUS EDUC_TWO_STATUS_t;
CONSTANT( int, rc_size_EDUC_TWO_STATUS, 5 );
enum ETHNICITY {
ETHNO_00,
ETHNO_01,
ETHNO_02,
ETHNO_03,
ETHNO_04,
ETHNO_05,
ETHNO_06,
ETHNO_07,
ETHNO_08,
ETHNO_09 };
typedef ETHNICITY ETHNICITY_t;
CONSTANT( int, rc_size_ETHNICITY, 10 );
enum ETHNICITY_SHORT {
ES_00,
ES_01,
ES_02,
ES_03,
ES_04,
ES_05,
ES_06 };
typedef ETHNICITY_SHORT ETHNICITY_SHORT_t;
CONSTANT( int, rc_size_ETHNICITY_SHORT, 7 );
enum EVER_UNION {
EU_NO,
EU_YES };
typedef EVER_UNION EVER_UNION_t;
CONSTANT( int, rc_size_EVER_UNION, 2 );
enum FAM_ROLE {
FR_HEAD,
FR_SPOUSE,
FR_CHILD };
typedef FAM_ROLE FAM_ROLE_t;
CONSTANT( int, rc_size_FAM_ROLE, 3 );
enum FERTILITY_MODEL {
FEM_BASE,
FEM_DETAIL,
FEM_ALIGNED_AGE,
FEM_ALIGNED_TOTAL };
typedef FERTILITY_MODEL FERTILITY_MODEL_t;
CONSTANT( int, rc_size_FERTILITY_MODEL, 4 );
enum GEO {
GEO_00,
GEO_01,
GEO_02,
GEO_03,
GEO_04,
GEO_05,
GEO_06,
GEO_07,
GEO_08,
GEO_09,
GEO_10,
GEO_11,
GEO_ABROAD };
typedef GEO GEO_t;
CONSTANT( int, rc_size_GEO, 13 );
enum GEO_ABR {
GEOA_00 };
typedef GEO_ABR GEO_ABR_t;
CONSTANT( int, rc_size_GEO_ABR, 1 );
enum GEO_NAT {
GEON_00,
GEON_01,
GEON_02,
GEON_03,
GEON_04,
GEON_05,
GEON_06,
GEON_07,
GEON_08,
GEON_09,
GEON_10,
GEON_11 };
typedef GEO_NAT GEO_NAT_t;
CONSTANT( int, rc_size_GEO_NAT, 12 );
enum GOT_PRENAT_CARE {
GPC_NO,
GPC_YES };
typedef GOT_PRENAT_CARE GOT_PRENAT_CARE_t;
CONSTANT( int, rc_size_GOT_PRENAT_CARE, 2 );
enum HCI_COEF {
HCI_EDUC,
HCI_ASR,
HCI_STUNT };
typedef HCI_COEF HCI_COEF_t;
CONSTANT( int, rc_size_HCI_COEF, 3 );
enum HIGHER_BIRTHS_PARA {
HBP_PERIOD1,
HBP_PERIOD2,
HBP_PERIOD3,
HBP_PERIOD4,
HBP_PERIOD5,
HBP_PERIOD6,
HBP_AGE35,
HBP_AGE40,
HBP_AGE45,
HBP_EDUC1,
HBP_EDUC2 };
typedef HIGHER_BIRTHS_PARA HIGHER_BIRTHS_PARA_t;
CONSTANT( int, rc_size_HIGHER_BIRTHS_PARA, 11 );
enum IMMI_POOL {
IP_STUDENTS,
IP_WORKER,
IP_REFUGEES };
typedef IMMI_POOL IMMI_POOL_t;
CONSTANT( int, rc_size_IMMI_POOL, 3 );
enum MORTALITY_MODEL {
MOM_BASE,
MOM_DETAIL,
MOM_ALIGNED_MACRO_TRENDS,
MOM_ALIGNED_MICRO_TRENDS };
typedef MORTALITY_MODEL MORTALITY_MODEL_t;
CONSTANT( int, rc_size_MORTALITY_MODEL, 4 );
enum MOTH_AGEGR {
CMA20,
CMA25,
CMA30,
CMA35,
CMA40,
CMA40P };
typedef MOTH_AGEGR MOTH_AGEGR_t;
CONSTANT( int, rc_size_MOTH_AGEGR, 6 );
enum PERSON_MICRODATA_COLUMNS {
PMC_FAMID,
PMC_WEIGHT,
PMC_ROLE,
PMC_BIRTH,
PMC_SEX,
PMC_EDUC,
PMC_SECGRAD,
PMC_SECATTEND,
PMC_UNION,
PMC_PARITY,
PMC_LASTBIR,
PMC_ETHNO,
PMC_GEO,
PMC_GEOBIR,
PMC_GEOPRE,
PMC_MOVEFIRST,
PMC_MOVELAST,
PMC_POOL };
typedef PERSON_MICRODATA_COLUMNS PERSON_MICRODATA_COLUMNS_t;
CONSTANT( int, rc_size_PERSON_MICRODATA_COLUMNS, 18 );
enum POP_POOL {
PP_NON,
PP_STUDENTS,
PP_WORKER,
PP_REFUGEES };
typedef POP_POOL POP_POOL_t;
CONSTANT( int, rc_size_POP_POOL, 4 );
enum PRE_SCHOOL_PARA {
PLP_ANY,
PLS_TWO };
typedef PRE_SCHOOL_PARA PRE_SCHOOL_PARA_t;
CONSTANT( int, rc_size_PRE_SCHOOL_PARA, 2 );
enum PRENATCARE_PARA {
PP_CONSTANT,
PP_EDUCMO_1,
PP_EDUCMO_2,
PP_REGION_1,
PP_REGION_2,
PP_REGION_3,
PP_REGION_4,
PP_COHORT_1,
PP_COHORT_2,
PP_COHORT_3,
PP_YOUNGMO };
typedef PRENATCARE_PARA PRENATCARE_PARA_t;
CONSTANT( int, rc_size_PRENATCARE_PARA, 11 );
enum REGION_INT {
REGI_00,
REGI_01,
REGI_02,
REGI_03,
REGI_04,
REGI_ABROAD };
typedef REGION_INT REGION_INT_t;
CONSTANT( int, rc_size_REGION_INT, 6 );
enum REGION_NAT {
REGN_00,
REGN_01,
REGN_02,
REGN_03,
REGN_04 };
typedef REGION_NAT REGION_NAT_t;
CONSTANT( int, rc_size_REGION_NAT, 5 );
enum SCHOOL_INFRA {
SI_TEACHERS,
SI_ROOMS };
typedef SCHOOL_INFRA SCHOOL_INFRA_t;
CONSTANT( int, rc_size_SCHOOL_INFRA, 2 );
enum SCHOOL_QUALITY_PARA {
SQP_AV,
SQP_SD };
typedef SCHOOL_QUALITY_PARA SCHOOL_QUALITY_PARA_t;
CONSTANT( int, rc_size_SCHOOL_QUALITY_PARA, 2 );
enum SELECTED_YEAR10 {
SY00,
SY10,
SY20,
SY30,
SY40 };
typedef SELECTED_YEAR10 SELECTED_YEAR10_t;
CONSTANT( int, rc_size_SELECTED_YEAR10, 5 );
enum SEX {
FEMALE,
MALE };
typedef SEX SEX_t;
CONSTANT( int, rc_size_SEX, 2 );
enum SPOUSE_GROUP {
SMG_00,
SMG_01,
SMG_02 };
typedef SPOUSE_GROUP SPOUSE_GROUP_t;
CONSTANT( int, rc_size_SPOUSE_GROUP, 3 );
enum TAB_SEL_YEARS {
TSY_2000,
TSY_2010,
TSY_2020,
TSY_2030,
TSY_2040,
TSY_2050,
TSY_2060 };
typedef TAB_SEL_YEARS TAB_SEL_YEARS_t;
CONSTANT( int, rc_size_TAB_SEL_YEARS, 7 );
enum UNION1_CHOICE {
U1C_CMN,
U1C_HAZARDS };
typedef UNION1_CHOICE UNION1_CHOICE_t;
CONSTANT( int, rc_size_UNION1_CHOICE, 2 );
enum UNION1_GROUP {
U1G_00,
U1G_01,
U1G_02 };
typedef UNION1_GROUP UNION1_GROUP_t;
CONSTANT( int, rc_size_UNION1_GROUP, 3 );
enum UNION1_PARA {
UP_MINAGE,
UP_AVERAGE,
UP_EVER };
typedef UNION1_PARA UNION1_PARA_t;
CONSTANT( int, rc_size_UNION1_PARA, 3 );
enum VACCINATION_PARA {
IP_CONSTANT,
IP_MALE,
IP_EDUCMO_1,
IP_EDUCMO_2,
IP_REGION_1,
IP_REGION_2,
IP_REGION_3,
IP_REGION_4,
IP_ETNO_1,
IP_ETNO_2,
IP_ETNO_3,
IP_ETNO_4,
IP_ETNO_5,
IP_ETNO_6,
IP_COHORT_1,
IP_COHORT_2,
IP_COHORT_3,
IP_YOUNGMO };
typedef VACCINATION_PARA VACCINATION_PARA_t;
CONSTANT( int, rc_size_VACCINATION_PARA, 18 );

// ranges
typedef int AGE_MIG;
CONSTANT( int, r_min_AGE_MIG, 0 );
CONSTANT( int, r_max_AGE_MIG, 12 );
CONSTANT( int, rc_size_AGE_MIG, 13 );
typedef AGE_MIG AGE_MIG_t;
typedef int AGE_RANGE;
CONSTANT( int, r_min_AGE_RANGE, 0 );
CONSTANT( int, r_max_AGE_RANGE, 100 );
CONSTANT( int, rc_size_AGE_RANGE, 101 );
typedef AGE_RANGE AGE_RANGE_t;
typedef int ALL_YEAR_RANGE;
CONSTANT( int, r_min_ALL_YEAR_RANGE, 1899 );
CONSTANT( int, r_max_ALL_YEAR_RANGE, 2150 );
CONSTANT( int, rc_size_ALL_YEAR_RANGE, 252 );
typedef ALL_YEAR_RANGE ALL_YEAR_RANGE_t;
typedef int CHILD_MORTALITY_AGE;
CONSTANT( int, r_min_CHILD_MORTALITY_AGE, 0 );
CONSTANT( int, r_max_CHILD_MORTALITY_AGE, 4 );
CONSTANT( int, rc_size_CHILD_MORTALITY_AGE, 5 );
typedef CHILD_MORTALITY_AGE CHILD_MORTALITY_AGE_t;
typedef int CHILD_MORTALITY_YEARS;
CONSTANT( int, r_min_CHILD_MORTALITY_YEARS, 2005 );
CONSTANT( int, r_max_CHILD_MORTALITY_YEARS, 2150 );
CONSTANT( int, rc_size_CHILD_MORTALITY_YEARS, 146 );
typedef CHILD_MORTALITY_YEARS CHILD_MORTALITY_YEARS_t;
typedef int EDUC_ONE_ENTRY_AGE;
CONSTANT( int, r_min_EDUC_ONE_ENTRY_AGE, 5 );
CONSTANT( int, r_max_EDUC_ONE_ENTRY_AGE, 8 );
CONSTANT( int, rc_size_EDUC_ONE_ENTRY_AGE, 4 );
typedef EDUC_ONE_ENTRY_AGE EDUC_ONE_ENTRY_AGE_t;
typedef int EDUC_ONE_GRADE;
CONSTANT( int, r_min_EDUC_ONE_GRADE, 1 );
CONSTANT( int, r_max_EDUC_ONE_GRADE, 6 );
CONSTANT( int, rc_size_EDUC_ONE_GRADE, 6 );
typedef EDUC_ONE_GRADE EDUC_ONE_GRADE_t;
typedef int EDUC_ONE_GRADE0;
CONSTANT( int, r_min_EDUC_ONE_GRADE0, 0 );
CONSTANT( int, r_max_EDUC_ONE_GRADE0, 6 );
CONSTANT( int, rc_size_EDUC_ONE_GRADE0, 7 );
typedef EDUC_ONE_GRADE0 EDUC_ONE_GRADE0_t;
typedef int EDUC_TWO_ENTRY_DUR;
CONSTANT( int, r_min_EDUC_TWO_ENTRY_DUR, 0 );
CONSTANT( int, r_max_EDUC_TWO_ENTRY_DUR, 2 );
CONSTANT( int, rc_size_EDUC_TWO_ENTRY_DUR, 3 );
typedef EDUC_TWO_ENTRY_DUR EDUC_TWO_ENTRY_DUR_t;
typedef int EDUC_TWO_GRADE;
CONSTANT( int, r_min_EDUC_TWO_GRADE, 1 );
CONSTANT( int, r_max_EDUC_TWO_GRADE, 6 );
CONSTANT( int, rc_size_EDUC_TWO_GRADE, 6 );
typedef EDUC_TWO_GRADE EDUC_TWO_GRADE_t;
typedef int EDUC_TWO_GRADE0;
CONSTANT( int, r_min_EDUC_TWO_GRADE0, 0 );
CONSTANT( int, r_max_EDUC_TWO_GRADE0, 6 );
CONSTANT( int, rc_size_EDUC_TWO_GRADE0, 7 );
typedef EDUC_TWO_GRADE0 EDUC_TWO_GRADE0_t;
typedef int FAM_ID;
CONSTANT( int, r_min_FAM_ID, 0 );
CONSTANT( int, r_max_FAM_ID, 2000000 );
CONSTANT( int, rc_size_FAM_ID, 2000001 );
typedef FAM_ID FAM_ID_t;
typedef int FEMALE_SPOUSE_AGE;
CONSTANT( int, r_min_FEMALE_SPOUSE_AGE, 10 );
CONSTANT( int, r_max_FEMALE_SPOUSE_AGE, 60 );
CONSTANT( int, rc_size_FEMALE_SPOUSE_AGE, 51 );
typedef FEMALE_SPOUSE_AGE FEMALE_SPOUSE_AGE_t;
typedef int FERTILE_AGE_RANGE;
CONSTANT( int, r_min_FERTILE_AGE_RANGE, 12 );
CONSTANT( int, r_max_FERTILE_AGE_RANGE, 49 );
CONSTANT( int, rc_size_FERTILE_AGE_RANGE, 38 );
typedef FERTILE_AGE_RANGE FERTILE_AGE_RANGE_t;
typedef int MALE_SPOUSE_AGE;
CONSTANT( int, r_min_MALE_SPOUSE_AGE, 15 );
CONSTANT( int, r_max_MALE_SPOUSE_AGE, 100 );
CONSTANT( int, rc_size_MALE_SPOUSE_AGE, 86 );
typedef MALE_SPOUSE_AGE MALE_SPOUSE_AGE_t;
typedef int PARITY_RANGE;
CONSTANT( int, r_min_PARITY_RANGE, 0 );
CONSTANT( int, r_max_PARITY_RANGE, 15 );
CONSTANT( int, rc_size_PARITY_RANGE, 16 );
typedef PARITY_RANGE PARITY_RANGE_t;
typedef int PARITY_RANGE1;
CONSTANT( int, r_min_PARITY_RANGE1, 1 );
CONSTANT( int, r_max_PARITY_RANGE1, 15 );
CONSTANT( int, rc_size_PARITY_RANGE1, 15 );
typedef PARITY_RANGE1 PARITY_RANGE1_t;
typedef int PARITY_RANGE2;
CONSTANT( int, r_min_PARITY_RANGE2, 2 );
CONSTANT( int, r_max_PARITY_RANGE2, 15 );
CONSTANT( int, rc_size_PARITY_RANGE2, 14 );
typedef PARITY_RANGE2 PARITY_RANGE2_t;
typedef int SIM_YEAR_RANGE;
CONSTANT( int, r_min_SIM_YEAR_RANGE, 2000 );
CONSTANT( int, r_max_SIM_YEAR_RANGE, 2150 );
CONSTANT( int, rc_size_SIM_YEAR_RANGE, 151 );
typedef SIM_YEAR_RANGE SIM_YEAR_RANGE_t;
typedef int YOB_EDUC_TRANS1;
CONSTANT( int, r_min_YOB_EDUC_TRANS1, 1988 );
CONSTANT( int, r_max_YOB_EDUC_TRANS1, 2050 );
CONSTANT( int, rc_size_YOB_EDUC_TRANS1, 63 );
typedef YOB_EDUC_TRANS1 YOB_EDUC_TRANS1_t;
typedef int YOB_EDUC_TRANS2;
CONSTANT( int, r_min_YOB_EDUC_TRANS2, 1984 );
CONSTANT( int, r_max_YOB_EDUC_TRANS2, 2050 );
CONSTANT( int, rc_size_YOB_EDUC_TRANS2, 67 );
typedef YOB_EDUC_TRANS2 YOB_EDUC_TRANS2_t;
typedef int YOB_UNION;
CONSTANT( int, r_min_YOB_UNION, 1951 );
CONSTANT( int, r_max_YOB_UNION, 2037 );
CONSTANT( int, rc_size_YOB_UNION, 87 );
typedef YOB_UNION YOB_UNION_t;

// partitions
typedef int AGE_MIG_PART;
CONSTANT( int, rc_size_AGE_MIG_PART, 13 );
extern real p_points_AGE_MIG_PART[];
typedef AGE_MIG_PART AGE_MIG_PART_t;
typedef int BIRTH_AGE_PART;
CONSTANT( int, rc_size_BIRTH_AGE_PART, 4 );
extern real p_points_BIRTH_AGE_PART[];
typedef BIRTH_AGE_PART BIRTH_AGE_PART_t;
typedef int CHILD_AGEGR_PART;
CONSTANT( int, rc_size_CHILD_AGEGR_PART, 7 );
extern real p_points_CHILD_AGEGR_PART[];
typedef CHILD_AGEGR_PART CHILD_AGEGR_PART_t;
typedef int DUR_TIME_IN_PARITY;
CONSTANT( int, rc_size_DUR_TIME_IN_PARITY, 6 );
extern real p_points_DUR_TIME_IN_PARITY[];
typedef DUR_TIME_IN_PARITY DUR_TIME_IN_PARITY_t;
typedef int IMMU_YOB_PART;
CONSTANT( int, rc_size_IMMU_YOB_PART, 4 );
extern real p_points_IMMU_YOB_PART[];
typedef IMMU_YOB_PART IMMU_YOB_PART_t;
typedef int MOTH_AGEGR_PART;
CONSTANT( int, rc_size_MOTH_AGEGR_PART, 6 );
extern real p_points_MOTH_AGEGR_PART[];
typedef MOTH_AGEGR_PART MOTH_AGEGR_PART_t;
typedef int tab_age15;
CONSTANT( int, rc_size_tab_age15, 4 );
extern real p_points_tab_age15[];
typedef tab_age15 tab_age15_t;
typedef int TAB_AGEC5;
CONSTANT( int, rc_size_TAB_AGEC5, 21 );
extern real p_points_TAB_AGEC5[];
typedef TAB_AGEC5 TAB_AGEC5_t;

class Parameters {
	public:
		Parameters();

		long	AgeImmigrantsScratch
[rc_size_SEX][CUM_RATES_ROWS_CNT]
[rc_size_AGE_RANGE];
		double	AgeImmiSearchMother;
		double	AgeLeavingHome;
		long	AgeOfImmigrantMother[CUM_RATES_ROWS_CNT]
[rc_size_FERTILE_AGE_RANGE];
		double	AgeSpecificFertility
[rc_size_FERTILE_AGE_RANGE]
[rc_size_SIM_YEAR_RANGE];
		double	AgeSpecificFertilityRate
[rc_size_FERTILE_AGE_RANGE]
[rc_size_SIM_YEAR_RANGE];
		double	BackMigrationHazard;
		double	BirthTrends
[rc_size_PARITY_RANGE1]
[rc_size_SIM_YEAR_RANGE];
		double	ChildMortalityBaseRisk
[rc_size_CHILD_MORTALITY_AGE]
[rc_size_SEX];
		double	ChildMortalityRelativeRisks
[rc_size_CHILD_MORTALITY_AGE]
[rc_size_CHILD_MORTALITY_GROUP];
		double	ChildMortalityTrend
[rc_size_CHILD_MORTALITY_AGE]
[rc_size_CHILD_MORTALITY_YEARS];
		double	ChildVaccinationOdds
[rc_size_GOT_PRENAT_CARE]
[rc_size_SIM_YEAR_RANGE]
[rc_size_VACCINATION_PARA];
		int	Educ1FirstCohortRefinedModel;
		double	Educ1GradOdds
[rc_size_EDUC1_GROUP]
[rc_size_SEX];
		double	Educ1Infrastructure
[rc_size_SCHOOL_INFRA]
[rc_size_GEO_NAT]
[rc_size_SIM_YEAR_RANGE];
		EDUC1_MODEL	Educ1Model;
		double	Educ1StartOdds
[rc_size_EDUC1_GROUP]
[rc_size_SEX];
		int	Educ2AllowedDelays;
		double	Educ2DelayedProgressionIntake
[rc_size_EDUC_TWO_GEO]
[rc_size_EDUC_TWO_GROUP]
[rc_size_EDUC_TWO_GRADE]
[rc_size_SIM_YEAR_RANGE];
		double	Educ2DelayedRepetitionIntake
[rc_size_EDUC_TWO_GEO]
[rc_size_EDUC_TWO_GROUP]
[rc_size_EDUC_TWO_GRADE]
[rc_size_SIM_YEAR_RANGE];
		double	Educ2DirectProgressionIntake
[rc_size_EDUC_TWO_GEO]
[rc_size_EDUC_TWO_GROUP]
[rc_size_EDUC_TWO_GRADE]
[rc_size_SIM_YEAR_RANGE];
		double	Educ2DirectRepetitionIntake
[rc_size_EDUC_TWO_GEO]
[rc_size_EDUC_TWO_GROUP]
[rc_size_EDUC_TWO_GRADE]
[rc_size_SIM_YEAR_RANGE];
		double	Educ2PeriodSuccess
[rc_size_EDUC_TWO_GEO]
[rc_size_EDUC_TWO_GROUP]
[rc_size_EDUC_TWO_GRADE]
[rc_size_SIM_YEAR_RANGE];
		long	EducOneDropoutGrade
[rc_size_EDUC_ONE_GEO]
[rc_size_EDUC_ONE_GROUP]
[rc_size_YOB_EDUC_TRANS2][CUM_RATES_ROWS_CNT]
[rc_size_EDUC_ONE_GRADE];
		long	EducOneEntryAge
[rc_size_EDUC_ONE_GEO]
[rc_size_EDUC_ONE_GROUP]
[rc_size_YOB_EDUC_TRANS2][CUM_RATES_ROWS_CNT]
[rc_size_EDUC_ONE_ENTRY_AGE];
		double	EducTrans1
[rc_size_SEX]
[rc_size_YOB_EDUC_TRANS1]
[rc_size_GEO];
		double	EducTrans2
[rc_size_SEX]
[rc_size_YOB_EDUC_TRANS2]
[rc_size_GEO];
		long	EmigrationDestination[CUM_RATES_ROWS_CNT]
[rc_size_GEO_ABR];
		double	EmigrationRatesDistrict
[rc_size_SEX]
[rc_size_AGE_MIG]
[rc_size_GEO_NAT];
		double	EndSchoolOneYear;
		double	EndSchoolTwoYear;
		long	EthnicityImmigrantsScratch
[rc_size_SEX]
[rc_size_GEO_NAT][CUM_RATES_ROWS_CNT]
[rc_size_ETHNICITY];
		long	EthnicTransmission
[rc_size_SEX]
[rc_size_ETHNICITY][CUM_RATES_ROWS_CNT]
[rc_size_ETHNICITY];
		FERTILITY_MODEL	FertilityModel;
		double	FirstBirthRates
[rc_size_BIRTH1_GROUP]
[rc_size_EVER_UNION]
[rc_size_FERTILE_AGE_RANGE]
[rc_size_BIRTH1_LOC];
		double	HCICoefficients
[rc_size_HCI_COEF];
		double	HigherOrderBirthsPara
[rc_size_HIGHER_BIRTHS_PARA]
[rc_size_PARITY_RANGE2];
		long	ImmiPoolDestination
[rc_size_IMMI_POOL]
[rc_size_SIM_YEAR_RANGE][CUM_RATES_ROWS_CNT]
[rc_size_GEO_NAT];
		double	ImmiPoolSize
[rc_size_IMMI_POOL]
[rc_size_SIM_YEAR_RANGE];
		long	ImmiScratchDestination
[rc_size_SEX]
[rc_size_AGE_MIG][CUM_RATES_ROWS_CNT]
[rc_size_GEO_NAT];
		double	InUnionProbNoChildren
[rc_size_FEMALE_SPOUSE_AGE]
[rc_size_UNION1_GROUP];
		double	InUnionProbWithChildren
[rc_size_UNION1_GROUP]
[rc_size_CHILD_AGEGR]
[rc_size_MOTH_AGEGR];
		double	LifeExpectancy
[rc_size_SIM_YEAR_RANGE]
[rc_size_SEX];
		CString	MicroDataInputFile;
		long	MicroDataInputFileSize;
		long	MigrationDestination
[rc_size_SEX]
[rc_size_GEO_NAT]
[rc_size_AGE_MIG][CUM_RATES_ROWS_CNT]
[rc_size_GEO_NAT];
		double	MigrationProbability
[rc_size_SEX]
[rc_size_AGE_MIG]
[rc_size_GEO_NAT];
		logical	MigrationTryKeepingFamiliesTogether;
		logical	ModelBackmigration;
		logical	ModelEmigration;
		logical	ModelImmigrationFromPools;
		logical	ModelImmigrationFromScratch;
		logical	ModelMigration;
		MORTALITY_MODEL	MortalityModel;
		double	MortalityTable
[rc_size_AGE_RANGE]
[rc_size_SEX];
		double	MortalityTrend
[rc_size_SIM_YEAR_RANGE]
[rc_size_SEX];
		double	NumberImmigrantsFromScratch
[rc_size_SIM_YEAR_RANGE]
[rc_size_SEX];
		double	PartnerAgeDistribution
[rc_size_FEMALE_SPOUSE_AGE]
[rc_size_MALE_SPOUSE_AGE];
		long	PartnerCharacteristicDistribution
[rc_size_SPOUSE_GROUP][CUM_RATES_ROWS_CNT]
[rc_size_SPOUSE_GROUP];
		double	PreNatalCareOdds
[rc_size_SIM_YEAR_RANGE]
[rc_size_PRENATCARE_PARA];
		double	PreSchoolAttendance
[rc_size_SEX]
[rc_size_REGION_NAT]
[rc_size_SIM_YEAR_RANGE]
[rc_size_PRE_SCHOOL_PARA];
		double	ProbStayWithMother;
		double	ProportionStunting
[rc_size_SEX]
[rc_size_REGION_NAT]
[rc_size_EDUC_ONE_LEVEL];
		double	ScalingFactor;
		double	SchoolOneInterruptionRate
[rc_size_EDUC_ONE_GEO]
[rc_size_EDUC_ONE_GROUP]
[rc_size_SIM_YEAR_RANGE];
		double	SchoolOneRepetitionRate
[rc_size_EDUC_ONE_GEO]
[rc_size_EDUC_ONE_GROUP]
[rc_size_SIM_YEAR_RANGE];
		double	SchoolQuality
[rc_size_REGION_NAT]
[rc_size_SCHOOL_QUALITY_PARA];
		double	SexRatio
[rc_size_SIM_YEAR_RANGE];
		double	StartPopSampleSize;
		double	StartSchoolOneYear;
		double	StartSchoolTwoYear;
		double	TotalFertilityRate
[rc_size_SIM_YEAR_RANGE];
		UNION1_CHOICE	Union1Choice;
		double	Union1FormationHazard
[rc_size_UNION1_GROUP]
[rc_size_ALL_YEAR_RANGE]
[rc_size_FEMALE_SPOUSE_AGE];
		double	Union1ParametersCMN
[rc_size_UNION1_GROUP]
[rc_size_YOB_UNION]
[rc_size_UNION1_PARA];
		double	Union1ParametersHazards
[rc_size_UNION1_GROUP]
[rc_size_FEMALE_SPOUSE_AGE]
[rc_size_ALL_YEAR_RANGE];
};

extern Parameters	*gprParam;

bool Lookup_AgeImmigrantsScratch( long lNumber, int nInd0, int *pnDest0 );
bool Lookup_AgeImmigrantsScratch( double dNumber, int nInd0, int *pnDest0 );
int EmptyDistributions_AgeImmigrantsScratch( void );
bool Lookup_AgeOfImmigrantMother( long lNumber, int *pnDest0 );
bool Lookup_AgeOfImmigrantMother( double dNumber, int *pnDest0 );
int EmptyDistributions_AgeOfImmigrantMother( void );
bool Lookup_EducOneDropoutGrade( long lNumber, int nInd0, int nInd1, int nInd2, int *pnDest0 );
bool Lookup_EducOneDropoutGrade( double dNumber, int nInd0, int nInd1, int nInd2, int *pnDest0 );
int EmptyDistributions_EducOneDropoutGrade( void );
bool Lookup_EducOneEntryAge( long lNumber, int nInd0, int nInd1, int nInd2, int *pnDest0 );
bool Lookup_EducOneEntryAge( double dNumber, int nInd0, int nInd1, int nInd2, int *pnDest0 );
int EmptyDistributions_EducOneEntryAge( void );
bool Lookup_EmigrationDestination( long lNumber, int *pnDest0 );
bool Lookup_EmigrationDestination( double dNumber, int *pnDest0 );
int EmptyDistributions_EmigrationDestination( void );
bool Lookup_EthnicityImmigrantsScratch( long lNumber, int nInd0, int nInd1, int *pnDest0 );
bool Lookup_EthnicityImmigrantsScratch( double dNumber, int nInd0, int nInd1, int *pnDest0 );
int EmptyDistributions_EthnicityImmigrantsScratch( void );
bool Lookup_EthnicTransmission( long lNumber, int nInd0, int nInd1, int *pnDest0 );
bool Lookup_EthnicTransmission( double dNumber, int nInd0, int nInd1, int *pnDest0 );
int EmptyDistributions_EthnicTransmission( void );
bool Lookup_ImmiPoolDestination( long lNumber, int nInd0, int nInd1, int *pnDest0 );
bool Lookup_ImmiPoolDestination( double dNumber, int nInd0, int nInd1, int *pnDest0 );
int EmptyDistributions_ImmiPoolDestination( void );
bool Lookup_ImmiScratchDestination( long lNumber, int nInd0, int nInd1, int *pnDest0 );
bool Lookup_ImmiScratchDestination( double dNumber, int nInd0, int nInd1, int *pnDest0 );
int EmptyDistributions_ImmiScratchDestination( void );
bool Lookup_MigrationDestination( long lNumber, int nInd0, int nInd1, int nInd2, int *pnDest0 );
bool Lookup_MigrationDestination( double dNumber, int nInd0, int nInd1, int nInd2, int *pnDest0 );
int EmptyDistributions_MigrationDestination( void );
bool Lookup_PartnerCharacteristicDistribution( long lNumber, int nInd0, int *pnDest0 );
bool Lookup_PartnerCharacteristicDistribution( double dNumber, int nInd0, int *pnDest0 );
int EmptyDistributions_PartnerCharacteristicDistribution( void );
ETHNICITY_SHORT ETHNICITY_To_ETHNICITY_SHORT( ETHNICITY eLevel );
REGION_INT GEO_To_REGION_INT( GEO eLevel );
BIRTH1_LOC GEO_To_BIRTH1_LOC( GEO eLevel );
GEO_NAT GEO_To_GEO_NAT( GEO eLevel );
REGION_NAT REGION_INT_To_REGION_NAT( REGION_INT eLevel );
bool IntIs_BIRTH1_GROUP( int nNumber );
BIRTH1_GROUP IntTo_BIRTH1_GROUP( int nNumber );
int IntFrom_BIRTH1_GROUP( BIRTH1_GROUP nLevel );
bool IntIs_BIRTH1_LOC( int nNumber );
BIRTH1_LOC IntTo_BIRTH1_LOC( int nNumber );
int IntFrom_BIRTH1_LOC( BIRTH1_LOC nLevel );
bool IntIs_CHILD_AGEGR( int nNumber );
CHILD_AGEGR IntTo_CHILD_AGEGR( int nNumber );
int IntFrom_CHILD_AGEGR( CHILD_AGEGR nLevel );
bool IntIs_CHILD_MORTALITY_GROUP( int nNumber );
CHILD_MORTALITY_GROUP IntTo_CHILD_MORTALITY_GROUP( int nNumber );
int IntFrom_CHILD_MORTALITY_GROUP( CHILD_MORTALITY_GROUP nLevel );
bool IntIs_EDUC1_GROUP( int nNumber );
EDUC1_GROUP IntTo_EDUC1_GROUP( int nNumber );
int IntFrom_EDUC1_GROUP( EDUC1_GROUP nLevel );
bool IntIs_ETHNICITY( int nNumber );
ETHNICITY IntTo_ETHNICITY( int nNumber );
int IntFrom_ETHNICITY( ETHNICITY nLevel );
bool IntIs_ETHNICITY_SHORT( int nNumber );
ETHNICITY_SHORT IntTo_ETHNICITY_SHORT( int nNumber );
int IntFrom_ETHNICITY_SHORT( ETHNICITY_SHORT nLevel );
bool IntIs_GEO( int nNumber );
GEO IntTo_GEO( int nNumber );
int IntFrom_GEO( GEO nLevel );
bool IntIs_GEO_ABR( int nNumber );
GEO_ABR IntTo_GEO_ABR( int nNumber );
int IntFrom_GEO_ABR( GEO_ABR nLevel );
bool IntIs_GEO_NAT( int nNumber );
GEO_NAT IntTo_GEO_NAT( int nNumber );
int IntFrom_GEO_NAT( GEO_NAT nLevel );
bool IntIs_HIGHER_BIRTHS_PARA( int nNumber );
HIGHER_BIRTHS_PARA IntTo_HIGHER_BIRTHS_PARA( int nNumber );
int IntFrom_HIGHER_BIRTHS_PARA( HIGHER_BIRTHS_PARA nLevel );
bool IntIs_MOTH_AGEGR( int nNumber );
MOTH_AGEGR IntTo_MOTH_AGEGR( int nNumber );
int IntFrom_MOTH_AGEGR( MOTH_AGEGR nLevel );
bool IntIs_PRENATCARE_PARA( int nNumber );
PRENATCARE_PARA IntTo_PRENATCARE_PARA( int nNumber );
int IntFrom_PRENATCARE_PARA( PRENATCARE_PARA nLevel );
bool IntIs_REGION_INT( int nNumber );
REGION_INT IntTo_REGION_INT( int nNumber );
int IntFrom_REGION_INT( REGION_INT nLevel );
bool IntIs_REGION_NAT( int nNumber );
REGION_NAT IntTo_REGION_NAT( int nNumber );
int IntFrom_REGION_NAT( REGION_NAT nLevel );
bool IntIs_SELECTED_YEAR10( int nNumber );
SELECTED_YEAR10 IntTo_SELECTED_YEAR10( int nNumber );
int IntFrom_SELECTED_YEAR10( SELECTED_YEAR10 nLevel );
bool IntIs_SPOUSE_GROUP( int nNumber );
SPOUSE_GROUP IntTo_SPOUSE_GROUP( int nNumber );
int IntFrom_SPOUSE_GROUP( SPOUSE_GROUP nLevel );
bool IntIs_TAB_SEL_YEARS( int nNumber );
TAB_SEL_YEARS IntTo_TAB_SEL_YEARS( int nNumber );
int IntFrom_TAB_SEL_YEARS( TAB_SEL_YEARS nLevel );
bool IntIs_UNION1_GROUP( int nNumber );
UNION1_GROUP IntTo_UNION1_GROUP( int nNumber );
int IntFrom_UNION1_GROUP( UNION1_GROUP nLevel );
bool IntIs_VACCINATION_PARA( int nNumber );
VACCINATION_PARA IntTo_VACCINATION_PARA( int nNumber );
int IntFrom_VACCINATION_PARA( VACCINATION_PARA nLevel );

class Calibrator : public ActorClass {
	public:
		TIME	__time;
		double	__events;
		class CalibratorStateFunctions	*poStateFunctions;
		class CalibratorDerivedStates	*poDerivedStates;
		void *operator new( size_t count );

		// links to other actors
		class Clock	*lCalibratorToClock;
		void	Set_lCalibratorToClock( Clock *value );
		class Calibrator_mlCalibratorToPerson	*mlCalibratorToPerson;

		// user-supplied data members
		double	alignment_educ1_high[rc_size_SEX][rc_size_GEO];
		double	alignment_educ1_medium[rc_size_SEX][rc_size_GEO];

		// user-supplied and required functions
		Calibrator() { };
		double AdjustedProbability ( double dProb , double dLogOddEduc , double dLogOddAdjust ) ; 
		void CalibratorYearEnd ( ) ; 
		void Finish ( ) ; 
		Person_ptr GetNextToEmigrate ( Person_ptr prPerson , GEO toGeo ) ; 
		Person_ptr GetNextToGiveBirth ( Person_ptr prPerson ) ; 
		Person_ptr GetNextToMigrate ( Person_ptr prPerson , GEO_NAT cGeoTo ) ; 
		void SetEduc1AdjustmentFactors ( ) ; 
		void Start ( ) ; 

		// special functions generated by the pre-compiler
		Calibrator( BOOL bInit ) { InitActor(); };
		~Calibrator() { DeleteActor(); };
		void	InitActor();
		void	DeleteActor();
		void	InitializeStates();
		void	UpdateLinkedStates();
		void	StartSpecial1();
		void	StartSpecial2();
		void	FinishSpecial();
		void	WaitUntil( double target_time );
		Calibrator	*UpdateTime( );
		void	WaitUntilThisActor(double target_time);
		void	WaitSpecial( TIME wait_time );

		// actor set nodes
		ActorSetNode	*poasCalibrator;

		// event functions and variables
		void __fastcall	ChildMortalityCalibration();
		TIME __fastcall	timeChildMortalityCalibration();
		Event	*poChildMortalityCalibration;

		// event processing functions
		double	EventTime( int nEventNum, int *pnEventInfo );
		void	Implement( int nEventNum, int nEventInfo );

		// state variables and their special functions
		long	Get_actor_id() { return actor_id; };
		double	Set_actor_subsample_weight( double value );
		double	Set_actor_weight( double value );

		// begin MembersDump
		//   group=<unknown token type>
		//   group=double
		double	actor_subsample_weight;
		double	actor_weight;
		TIME	age;
		double	case_seed;
		double	mort_female_0;
		double	mort_female_1;
		double	mort_female_2;
		double	mort_female_3;
		double	mort_female_4;
		double	mort_male_0;
		double	mort_male_1;
		double	mort_male_2;
		double	mort_male_3;
		double	mort_male_4;
		TIME	time;
		//   group=float
		//   group=long
		long	actor_id;
		long	case_id;
		//   group=ulong
		//   group=int
		integer	emigration_counter;
		integer	immigration_counter;
		logical	is_calibrated_child_mortality;
		integer	migration_counter;
		logical	actorset_filter_asCalibrator;
		ALL_YEAR_RANGE	calibrator_year;
		int	lCalibratorToClock_clock_year;
		logical	__finished;
		//   group=uint
		//   group=short
		//   group=ushort
		//   group=char
		//   group=uchar
		// end MembersDump

};

class CalibratorStateFunctions {
	public:
		CalibratorStateFunctions() { };
		class Calibrator	*poParent;
		class CalibratorDerivedStates	*poDerivedStates;

		// state special functions
		double	Set_actor_subsample_weight( double value ){ return poParent->actor_subsample_weight = value; };
		double	Set_actor_weight( double value ){ return poParent->actor_weight = value; };
		TIME	Set_age( TIME value, bool bStartFunction = false ){ return poParent->age = value; };
		logical	Set_is_calibrated_child_mortality( logical value );
		double	Set_mort_female_0( double value );
		double	Set_mort_female_1( double value );
		double	Set_mort_female_2( double value );
		double	Set_mort_female_3( double value );
		double	Set_mort_female_4( double value );
		double	Set_mort_male_0( double value );
		double	Set_mort_male_1( double value );
		double	Set_mort_male_2( double value );
		double	Set_mort_male_3( double value );
		double	Set_mort_male_4( double value );
		TIME	Set_time( TIME value, bool bStartFunction = false ){ return poParent->time = value; };
		logical	Set_actorset_filter_asCalibrator( logical value );
		ALL_YEAR_RANGE	Set_calibrator_year( ALL_YEAR_RANGE value );
		int	Set_lCalibratorToClock_clock_year( int value );
};

class CalibratorDerivedStates {
	public:
		CalibratorDerivedStates() { };
		class Calibrator	*poParent;
		class CalibratorStateFunctions	*poStateFunctions;

		// derived state variables and their special functions

		// begin MembersDump
		//   group=<unknown token type>
		//   group=double
		//   group=float
		//   group=long
		//   group=ulong
		//   group=int
		//   group=uint
		//   group=short
		//   group=ushort
		//   group=char
		//   group=uchar
		// end MembersDump

};

class Clock : public ActorClass {
	public:
		TIME	__time;
		double	__events;
		class ClockStateFunctions	*poStateFunctions;
		class ClockDerivedStates	*poDerivedStates;
		void *operator new( size_t count );

		// links to other actors
		class Calibrator	*lClockToCalibrator;
		void	Set_lClockToCalibrator( Calibrator *value );

		// user-supplied data members

		// user-supplied and required functions
		Clock() { };
		void Finish ( ) ; 
		void SetSchoolYearOneClock ( ) ; 
		void SetSchoolYearTwoClock ( ) ; 
		void Start ( ) ; 
		void UpdatePartnershipStatus ( ) ; 

		// special functions generated by the pre-compiler
		Clock( BOOL bInit ) { InitActor(); };
		~Clock() { DeleteActor(); };
		void	InitActor();
		void	DeleteActor();
		void	InitializeStates();
		void	UpdateLinkedStates();
		void	StartSpecial1();
		void	StartSpecial2();
		void	FinishSpecial();
		void	WaitUntil( double target_time );
		Clock	*UpdateTime( );
		void	WaitUntilThisActor(double target_time);
		void	WaitSpecial( TIME wait_time );

		// actor set nodes
		ActorSetNode	*poasClock;

		// event functions and variables
		void __fastcall	ClockMidyearEvent();
		TIME __fastcall	timeClockMidyearEvent();
		Event	*poClockMidyearEvent;
		void __fastcall	ClockYearEndEvent();
		TIME __fastcall	timeClockYearEndEvent();
		Event	*poClockYearEndEvent;
		void __fastcall	ClockYearStartEvent();
		TIME __fastcall	timeClockYearStartEvent();
		Event	*poClockYearStartEvent;
		void __fastcall	EndSchoolOneYearEvent();
		TIME __fastcall	timeEndSchoolOneYearEvent();
		Event	*poEndSchoolOneYearEvent;
		void __fastcall	EndSchoolTwoYearEvent();
		TIME __fastcall	timeEndSchoolTwoYearEvent();
		Event	*poEndSchoolTwoYearEvent;
		void __fastcall	StartSchoolOneYearEvent();
		TIME __fastcall	timeStartSchoolOneYearEvent();
		Event	*poStartSchoolOneYearEvent;
		void __fastcall	StartSchoolTwoYearEvent();
		TIME __fastcall	timeStartSchoolTwoYearEvent();
		Event	*poStartSchoolTwoYearEvent;

		// event processing functions
		double	EventTime( int nEventNum, int *pnEventInfo );
		void	Implement( int nEventNum, int nEventInfo );

		// state variables and their special functions
		long	Get_actor_id() { return actor_id; };
		double	Set_actor_subsample_weight( double value );
		double	Set_actor_weight( double value );

		// begin MembersDump
		//   group=<unknown token type>
		//   group=double
		double	actor_subsample_weight;
		double	actor_weight;
		TIME	age;
		double	case_seed;
		TIME	next_clock_year_end;
		TIME	next_clock_year_start;
		TIME	next_midyear_clock_event;
		TIME	time;
		TIME	time_end_school_one_year;
		TIME	time_end_school_two_year;
		TIME	time_start_school_one_year;
		TIME	time_start_school_two_year;
		//   group=float
		//   group=long
		long	actor_id;
		long	case_id;
		//   group=ulong
		//   group=int
		int	clock_year;
		logical	actorset_filter_asClock;
		logical	__finished;
		//   group=uint
		//   group=short
		//   group=ushort
		//   group=char
		//   group=uchar
		// end MembersDump

};

class ClockStateFunctions {
	public:
		ClockStateFunctions() { };
		class Clock	*poParent;
		class ClockDerivedStates	*poDerivedStates;

		// state special functions
		double	Set_actor_subsample_weight( double value ){ return poParent->actor_subsample_weight = value; };
		double	Set_actor_weight( double value ){ return poParent->actor_weight = value; };
		TIME	Set_age( TIME value, bool bStartFunction = false ){ return poParent->age = value; };
		int	Set_clock_year( int value );
		int	PlusPlusSet_clock_year() { int old_value; old_value = poParent->clock_year;
			Set_clock_year((int) (poParent->clock_year + 1)); return old_value; };
		TIME	Set_next_clock_year_end( TIME value );
		TIME	Set_next_clock_year_start( TIME value );
		TIME	Set_next_midyear_clock_event( TIME value );
		TIME	Set_time( TIME value, bool bStartFunction = false ){ return poParent->time = value; };
		TIME	Set_time_end_school_one_year( TIME value );
		TIME	Set_time_end_school_two_year( TIME value );
		TIME	Set_time_start_school_one_year( TIME value );
		TIME	Set_time_start_school_two_year( TIME value );
		logical	Set_actorset_filter_asClock( logical value );
};

class ClockDerivedStates {
	public:
		ClockDerivedStates() { };
		class Clock	*poParent;
		class ClockStateFunctions	*poStateFunctions;

		// derived state variables and their special functions

		// begin MembersDump
		//   group=<unknown token type>
		//   group=double
		//   group=float
		//   group=long
		//   group=ulong
		//   group=int
		//   group=uint
		//   group=short
		//   group=ushort
		//   group=char
		//   group=uchar
		// end MembersDump

};

class Observation : public ActorClass {
	public:
		TIME	__time;
		double	__events;
		class ObservationStateFunctions	*poStateFunctions;
		class ObservationDerivedStates	*poDerivedStates;
		void *operator new( size_t count );

		// links to other actors

		// user-supplied data members
		double	pmc[rc_size_PERSON_MICRODATA_COLUMNS];

		// user-supplied and required functions
		Observation() { };
		void Finish ( ) ; 
		void Start ( const input_csv & input ) ; 

		// special functions generated by the pre-compiler
		Observation( BOOL bInit ) { InitActor(); };
		~Observation() { DeleteActor(); };
		void	InitActor();
		void	DeleteActor();
		void	InitializeStates();
		void	UpdateLinkedStates();
		void	StartSpecial1();
		void	StartSpecial2();
		void	FinishSpecial();
		void	WaitUntil( double target_time );
		Observation	*UpdateTime( );
		void	WaitUntilThisActor(double target_time);
		void	WaitSpecial( TIME wait_time );

		// actor set nodes
		ActorSetNode	*poasObservationAll;
		ActorSetNode	*poasObservationByFam;
		ActorSetNode	*poasObservationByFamOldest;
		ActorSetNode	*poasObservationHeads;
		ActorSetNode	*poasObservationNonHeads;
		ActorSetNode	*poasObservations;
		ActorSetNode	*poasSimulatedObservationHeads;

		// event functions and variables

		// event processing functions
		double	EventTime( int nEventNum, int *pnEventInfo );
		void	Implement( int nEventNum, int nEventInfo );

		// state variables and their special functions
		long	Get_actor_id() { return actor_id; };
		double	Set_actor_subsample_weight( double value );
		double	Set_actor_weight( double value );

		// begin MembersDump
		//   group=<unknown token type>
		//   group=double
		double	actor_subsample_weight;
		double	actor_weight;
		TIME	age;
		double	case_seed;
		TIME	obs_birth;
		TIME	time;
		//   group=float
		//   group=long
		long	actor_id;
		long	case_id;
		//   group=ulong
		//   group=int
		FAM_ID	fam_id;
		logical	obs_oldest;
		integer	obs_weight;
		POP_POOL	pop_pool;
		logical	actorset_filter_asObservationAll;
		logical	actorset_filter_asObservationByFam;
		logical	actorset_filter_asObservationByFamOldest;
		logical	actorset_filter_asObservationHeads;
		logical	actorset_filter_asObservationNonHeads;
		logical	actorset_filter_asObservations;
		logical	actorset_filter_asSimulatedObservationHeads;
		logical	__finished;
		//   group=uint
		//   group=short
		//   group=ushort
		//   group=char
		//   group=uchar
		// end MembersDump

};

class ObservationStateFunctions {
	public:
		ObservationStateFunctions() { };
		class Observation	*poParent;
		class ObservationDerivedStates	*poDerivedStates;

		// state special functions
		double	Set_actor_subsample_weight( double value ){ return poParent->actor_subsample_weight = value; };
		double	Set_actor_weight( double value ){ return poParent->actor_weight = value; };
		TIME	Set_age( TIME value, bool bStartFunction = false ){ return poParent->age = value; };
		FAM_ID	Set_fam_id( FAM_ID value );
		logical	Set_obs_oldest( logical value );
		integer	Set_obs_weight( integer value );
		integer	MinusMinusSet_obs_weight() { integer old_value; old_value = poParent->obs_weight;
			Set_obs_weight((integer) (poParent->obs_weight - 1)); return old_value; };
		POP_POOL	Set_pop_pool( POP_POOL value );
		TIME	Set_time( TIME value, bool bStartFunction = false ){ return poParent->time = value; };
		logical	Set_actorset_filter_asObservationAll( logical value );
		logical	Set_actorset_filter_asObservationByFam( logical value );
		logical	Set_actorset_filter_asObservationByFamOldest( logical value );
		logical	Set_actorset_filter_asObservationHeads( logical value );
		logical	Set_actorset_filter_asObservationNonHeads( logical value );
		logical	Set_actorset_filter_asObservations( logical value );
		logical	Set_actorset_filter_asSimulatedObservationHeads( logical value );
};

class ObservationDerivedStates {
	public:
		ObservationDerivedStates() { };
		class Observation	*poParent;
		class ObservationStateFunctions	*poStateFunctions;

		// derived state variables and their special functions

		// begin MembersDump
		//   group=<unknown token type>
		//   group=double
		//   group=float
		//   group=long
		//   group=ulong
		//   group=int
		//   group=uint
		//   group=short
		//   group=ushort
		//   group=char
		//   group=uchar
		// end MembersDump

};

class PersonTable {
	public:
		PersonTable() { };

		// begin MembersDump
		//   group=<unknown token type>
		//   group=double
		double	tab0events;
		real	tab0open0;
		real	tab0close0;
		real	tab0open1;
		real	tab0close1;
		real	tab0open2;
		real	tab0close2;
		real	tab0open3;
		real	tab0close3;
		double	tab1events;
		real	tab1open0;
		real	tab1close0;
		real	tab1open1;
		real	tab1close1;
		real	tab1open2;
		real	tab1close2;
		real	tab1open3;
		real	tab1close3;
		double	tab2events;
		real	tab2open0;
		real	tab2close0;
		real	tab2open1;
		real	tab2close1;
		real	tab2open2;
		real	tab2close2;
		double	tab3events;
		real	tab3open0;
		real	tab3close0;
		double	tab4events;
		real	tab4open0;
		real	tab4close0;
		real	tab4open1;
		real	tab4close1;
		real	tab4open2;
		real	tab4close2;
		double	tab5events;
		real	tab5open0;
		real	tab5close0;
		real	tab5open1;
		real	tab5close1;
		real	tab5open2;
		real	tab5close2;
		real	tab5open3;
		real	tab5close3;
		real	tab5open4;
		real	tab5close4;
		double	tab6events;
		real	tab6open0;
		real	tab6close0;
		double	tab7events;
		real	tab7open0;
		real	tab7close0;
		real	tab7open1;
		real	tab7close1;
		real	tab7open2;
		real	tab7close2;
		real	tab7open3;
		real	tab7close3;
		double	tab8events;
		real	tab8open0;
		real	tab8close0;
		real	tab8open1;
		real	tab8close1;
		real	tab8open2;
		real	tab8close2;
		real	tab8open3;
		real	tab8close3;
		double	tab9events;
		real	tab9open0;
		real	tab9close0;
		double	tab10events;
		real	tab10open0;
		real	tab10close0;
		double	tab11events;
		real	tab11open0;
		real	tab11close0;
		real	tab11open1;
		real	tab11close1;
		double	tab12events;
		real	tab12open0;
		real	tab12close0;
		double	tab13events;
		real	tab13open0;
		real	tab13close0;
		real	tab13open1;
		real	tab13close1;
		real	tab13open2;
		real	tab13close2;
		real	tab13open3;
		real	tab13close3;
		double	tab14events;
		real	tab14open0;
		real	tab14close0;
		real	tab14open1;
		real	tab14close1;
		real	tab14open2;
		real	tab14close2;
		double	tab15events;
		real	tab15open0;
		real	tab15close0;
		real	tab15open1;
		real	tab15close1;
		real	tab15open2;
		real	tab15close2;
		real	tab15open3;
		real	tab15close3;
		real	tab15open4;
		real	tab15close4;
		real	tab15open5;
		real	tab15close5;
		real	tab15open6;
		real	tab15close6;
		real	tab15open7;
		real	tab15close7;
		real	tab15open8;
		real	tab15close8;
		double	tab16events;
		real	tab16open0;
		real	tab16close0;
		real	tab16open1;
		real	tab16close1;
		real	tab16open2;
		real	tab16close2;
		real	tab16open3;
		real	tab16close3;
		real	tab16open4;
		real	tab16close4;
		real	tab16open5;
		real	tab16close5;
		real	tab16open6;
		real	tab16close6;
		real	tab16open7;
		real	tab16close7;
		real	tab16open8;
		real	tab16close8;
		double	tab17events;
		real	tab17open0;
		real	tab17close0;
		double	tab18events;
		real	tab18open0;
		real	tab18close0;
		double	tab19events;
		real	tab19open0;
		real	tab19close0;
		double	tab20events;
		real	tab20open0;
		real	tab20close0;
		real	tab20open1;
		real	tab20close1;
		double	tab21events;
		real	tab21open0;
		real	tab21close0;
		real	tab21open1;
		real	tab21close1;
		double	tab22events;
		real	tab22open0;
		real	tab22close0;
		real	tab22open1;
		real	tab22close1;
		double	tab23events;
		real	tab23open0;
		real	tab23close0;
		real	tab23open1;
		real	tab23close1;
		double	tab24events;
		real	tab24open0;
		real	tab24close0;
		real	tab24open1;
		real	tab24close1;
		//   group=float
		//   group=long
		//   group=ulong
		//   group=int
		logical	tab0dim_filter;
		index	tab0cell;
		index	tab0total_cell;
		logical	tab1dim_filter;
		index	tab1cell;
		index	tab1total_cell;
		logical	tab2dim_filter;
		index	tab2cell;
		index	tab2total_cell;
		logical	tab3dim_filter;
		index	tab3cell;
		index	tab3total_cell;
		logical	tab4dim_filter;
		index	tab4cell;
		index	tab4total_cell;
		logical	tab5dim_filter;
		index	tab5cell;
		index	tab5total_cell;
		logical	tab6dim_filter;
		index	tab6cell;
		index	tab6total_cell;
		logical	tab7dim_filter;
		index	tab7cell;
		index	tab7total_cell;
		logical	tab8dim_filter;
		index	tab8cell;
		index	tab8total_cell;
		logical	tab9dim_filter;
		index	tab9cell;
		index	tab9total_cell;
		logical	tab10dim_filter;
		index	tab10cell;
		index	tab10total_cell;
		logical	tab11dim_filter;
		index	tab11cell;
		index	tab11total_cell;
		logical	tab12dim_filter;
		index	tab12cell;
		index	tab12total_cell;
		logical	tab13dim_filter;
		index	tab13cell;
		index	tab13total_cell;
		logical	tab14dim_filter;
		index	tab14cell;
		index	tab14total_cell;
		logical	tab15dim_filter;
		index	tab15cell;
		index	tab15total_cell;
		logical	tab16dim_filter;
		index	tab16cell;
		index	tab16total_cell;
		logical	tab17dim_filter;
		index	tab17cell;
		index	tab17total_cell;
		logical	tab18dim_filter;
		index	tab18cell;
		index	tab18total_cell;
		logical	tab19dim_filter;
		index	tab19cell;
		index	tab19total_cell;
		logical	tab20dim_filter;
		index	tab20cell;
		index	tab20total_cell;
		logical	tab21dim_filter;
		index	tab21cell;
		index	tab21total_cell;
		logical	tab22dim_filter;
		index	tab22cell;
		index	tab22total_cell;
		logical	tab23dim_filter;
		index	tab23cell;
		index	tab23total_cell;
		logical	tab24dim_filter;
		index	tab24cell;
		index	tab24total_cell;
		//   group=uint
		//   group=short
		//   group=ushort
		//   group=char
		//   group=uchar
		// end MembersDump

};

class Person : public ActorClass {
	public:
		TIME	__time;
		double	__events;
		class PersonStateFunctions	*poStateFunctions;
		class PersonDerivedStates	*poDerivedStates;
		class PersonTable	tableinfo;
		void *operator new( size_t count );

		// links to other actors
		class Person	*lBioFather;
		void	Set_lBioFather( Person *value );
		class Person	*lBioMother;
		void	Set_lBioMother( Person *value );
		class Calibrator	*lCalibrator;
		void	Set_lCalibrator( Calibrator *value );
		class Person	*lHHFather;
		void	Set_lHHFather( Person *value );
		class Person	*lHHMother;
		void	Set_lHHMother( Person *value );
		class Person	*lSpouse;
		void	Set_lSpouse( Person *value );
		class StartpopValues	*lStartValues;
		void	Set_lStartValues( StartpopValues *value );
		class Person_mlBioFatherChildren	*mlBioFatherChildren;
		class Person_mlBioMotherChildren	*mlBioMotherChildren;
		class Person_mlHHFatherChildren	*mlHHFatherChildren;
		class Person_mlHHMotherChildren	*mlHHMotherChildren;

		// user-supplied data members
		Person_ptr	ptr_creator;
		Observation_ptr	ptr_observation;

		// user-supplied and required functions
		Person() { };
		void AdjustEducOne ( ) ; 
		void CalculateHCIVariables ( ) ; 
		void DecideImmunizationStatusImmigrants ( ) ; 
		void DecideImmunizationStatusResidents ( ) ; 
		void DecideStuntingFate ( ) ; 
		void doDissolveUnion ( ) ; 
		void doGiveBirth ( ) ; 
		void doIncreaseParity ( ) ; 
		void doLeaveParentalHome ( ) ; 
		void doLinkToFamilyAtStart ( ) ; 
		void doMaintainLinksAtDeath ( ) ; 
		void doResidentialMove ( GEO cDestination ) ; 
		void doStartUnion ( Person_ptr ptrPartner ) ; 
		void FindImmigrantMother ( ) ; 
		logical FindSpouse ( ) ; 
		void Finish ( ) ; 
		ETHNICITY GetImmigrantsScratchEthnicity ( GEO_NAT toGeo ) ; 
		ETHNICITY GetInheritedEthnicity ( ETHNICITY eMothersEthnicity ) ; 
		double getTimeToBirth ( ) ; 
		void HandleFertility ( ) ; 
		void HandleMortality ( ) ; 
		void ImputeCharacteristicsAtFirstImmigrationScratch ( GEO_NAT cGeoDestination ) ; 
		void MidYear ( ) ; 
		void SetEduc1BaseFate ( ) ; 
		void SetEducOneEntryAgeDroputGrade ( ) ; 
		void SetGeobirthTimeofbirthCtScratch ( ) ; 
		void Start ( Observation_ptr peObs , Person * peCreator , int nYearOfImmigration , SEX nImmiSex ) ; 
		void YearEnd ( ) ; 
		void YearStart ( ) ; 
		void	MidYearIMPLEMENT_HOOK();
		void	YearEndIMPLEMENT_HOOK();
		void	YearStartIMPLEMENT_HOOK();
		void	BirthdayEventIMPLEMENT_HOOK();

		// special functions generated by the pre-compiler
		Person( BOOL bInit ) { InitActor(); };
		~Person() { DeleteActor(); };
		void	InitActor();
		void	DeleteActor();
		void	InitializeStates();
		void	UpdateLinkedStates();
		void	StartSpecial1();
		void	StartSpecial2();
		void	FinishSpecial();
		void	WaitUntil( double target_time );
		Person	*UpdateTime( );
		void	WaitUntilThisActor(double target_time);
		void	WaitSpecial( TIME wait_time );

		// actor set nodes
		ActorSetNode	*poasAllFertilePersonsForFertilityAlignment;
		ActorSetNode	*poasAllFertilePersonsForFertilityAlignmentByAge;
		ActorSetNode	*poasAllPerson;
		ActorSetNode	*poasAvailableMale;
		ActorSetNode	*poasAvailableMaleByType;
		ActorSetNode	*poasAvailableMaleByTypeImmiScratch;
		ActorSetNode	*poasAvailableMaleImmiScratch;
		ActorSetNode	*poasFemaleInUnionByAgeAndPartnerAge;
		ActorSetNode	*poasFemaleInUnionByAgeAndPartnerAgeImmiScratch;
		ActorSetNode	*poasPotentialImmigrantMothers;
		ActorSetNode	*poasPotentialSchoolOneStudent;
		ActorSetNode	*poasPotentialSchoolOneStudentToProcess;
		ActorSetNode	*poasPotentialSchoolTwoStudent;
		ActorSetNode	*poasPotentialSchoolTwoStudentToProcess;
		ActorSetNode	*poasResidentBabies;
		ActorSetNode	*poasResidentsAge0Sex;
		ActorSetNode	*poasResidentsAge0SexGeo;
		ActorSetNode	*poasSimBornAge0;
		ActorSetNode	*poasWantToMove;
		ActorSetNode	*poasWomenNoChildren;
		ActorSetNode	*poasWomenWithChildren;

		// event functions and variables
		void __fastcall	ScheduleEvent();
		TIME __fastcall	timeScheduleEvent();
		Event	*poScheduleEvent;
		TIME schedule_event_time;
		void __fastcall	BackMigrationEvent();
		TIME __fastcall	timeBackMigrationEvent();
		Event	*poBackMigrationEvent;
		void __fastcall	BirthdayEvent();
		TIME __fastcall	timeBirthdayEvent();
		Event	*poBirthdayEvent;
		void __fastcall	DeathAtMaxLifespanEvent();
		TIME __fastcall	timeDeathAtMaxLifespanEvent();
		Event	*poDeathAtMaxLifespanEvent;
		void __fastcall	EmigrationEvent();
		TIME __fastcall	timeEmigrationEvent();
		Event	*poEmigrationEvent;
		void __fastcall	FertilityBaseEvent();
		TIME __fastcall	timeFertilityBaseEvent();
		Event	*poFertilityBaseEvent;
		void __fastcall	FertilityDetailedEvent();
		TIME __fastcall	timeFertilityDetailedEvent();
		Event	*poFertilityDetailedEvent;
		void __fastcall	FirstImmigrationFromPoolEvent();
		TIME __fastcall	timeFirstImmigrationFromPoolEvent();
		Event	*poFirstImmigrationFromPoolEvent;
		void __fastcall	FirstImmigrationFromScratchEvent();
		TIME __fastcall	timeFirstImmigrationFromScratchEvent();
		Event	*poFirstImmigrationFromScratchEvent;
		void __fastcall	FirstUnionFormationEvent();
		TIME __fastcall	timeFirstUnionFormationEvent();
		Event	*poFirstUnionFormationEvent;
		void __fastcall	LeavingHomeEvent();
		TIME __fastcall	timeLeavingHomeEvent();
		Event	*poLeavingHomeEvent;
		void __fastcall	MigrationEvent();
		TIME __fastcall	timeMigrationEvent();
		Event	*poMigrationEvent;
		void __fastcall	MortalityBaseEvent();
		TIME __fastcall	timeMortalityBaseEvent();
		Event	*poMortalityBaseEvent;
		void __fastcall	MortalityDetailedEvent();
		TIME __fastcall	timeMortalityDetailedEvent();
		Event	*poMortalityDetailedEvent;
		void __fastcall	SetAliveEvent();
		TIME __fastcall	timeSetAliveEvent();
		Event	*poSetAliveEvent;
		void __fastcall	SetPreschoolYearsEvent();
		TIME __fastcall	timeSetPreschoolYearsEvent();
		Event	*poSetPreschoolYearsEvent;

		// event processing functions
		double	EventTime( int nEventNum, int *pnEventInfo );
		void	Implement( int nEventNum, int nEventInfo );

		// state variables and their special functions
		long	Get_actor_id() { return actor_id; };
		double	Set_actor_subsample_weight( double value );
		double	Set_actor_weight( double value );

		// begin MembersDump
		//   group=<unknown token type>
		//   group=double
		double	actor_subsample_weight;
		double	actor_weight;
		double	adult_survival;
		TIME	age;
		double	case_seed;
		double	ind_hci;
		double	mother_age_at_birth;
		double	quality_adjusted_schooling;
		double	quality_of_schooling;
		TIME	time;
		TIME	time_first_move;
		TIME	time_last_move;
		TIME	time_next_birthday;
		TIME	time_of_birth;
		TIME	time_of_first_immigration;
		double	age_last_birth;
		double	age_youngest_child;
		double	child_mortality;
		double	educ_one_required_rooms;
		double	educ_one_required_rooms_value_out;
		double	educ_one_required_teachers;
		double	educ_one_required_teachers_value_out;
		double	lCalibrator_mort_female_0;
		double	lCalibrator_mort_female_1;
		double	lCalibrator_mort_female_2;
		double	lCalibrator_mort_female_3;
		double	lCalibrator_mort_female_4;
		double	lCalibrator_mort_male_0;
		double	lCalibrator_mort_male_1;
		double	lCalibrator_mort_male_2;
		double	lCalibrator_mort_male_3;
		double	lCalibrator_mort_male_4;
		//   group=float
		//   group=long
		long	actor_id;
		long	case_id;
		//   group=ulong
		//   group=int
		ALL_YEAR_RANGE	calendar_year;
		CREATION_TYPE	creation_type;
		EDUC_ONE_LEVEL	educ_mother;
		EDUC_ONE_ENTRY_AGE	educ_one_entry_age;
		EDUC_ONE_LEVEL	educ_one_fate;
		EDUC_ONE_GEO	educ_one_geo;
		EDUC_ONE_GRADE0	educ_one_grade_attended;
		EDUC_ONE_GRADE	educ_one_grade_fate;
		EDUC_ONE_GRADE0	educ_one_grade_passed;
		EDUC_ONE_GROUP	educ_one_group;
		EDUC_ONE_STATUS	educ_one_status;
		logical	educ_one_to_process;
		int	educ_two_delay;
		EDUC_TWO_GEO	educ_two_geo;
		EDUC_TWO_GRADE0	educ_two_grade_attended;
		EDUC_TWO_GRADE0	educ_two_grade_passed;
		EDUC_TWO_GROUP	educ_two_group;
		EDUC_TWO_STATUS	educ_two_status;
		logical	educ_two_to_process;
		ETHNICITY	ethnicity;
		logical	ever_resident;
		logical	ever_union;
		logical	first_move_done;
		GEO	geo;
		GEO	geo_birth;
		GEO	geo_prev;
		GEO	geo_want_to_move;
		GOT_PRENAT_CARE	got_prenat_care;
		AGE_RANGE	integer_age;
		logical	is_alive;
		logical	is_blocked_from_marriage;
		logical	is_immunized;
		logical	is_ready_for_birthtables;
		logical	is_stunted;
		PARITY_RANGE	parity;
		logical	preschool_is_decided;
		logical	ready_to_set_alive;
		SEX	sex;
		counter	table_unit_TabEducFateByGroup;
		counter	table_unit_TabEducFateDistrYob;
		counter	table_unit_TabPrimarySchoolPlanning;
		counter	table_unit_TabSchool2TrackTab;
		counter	table_unit_tabStuntingSexRegMotherYob;
		counter	table_unit_tabStuntingSexRegYob;
		logical	this_parity_spell;
		int	years_preschool;
		logical	actorset_filter_asAllFertilePersonsForFertilityAlignment;
		logical	actorset_filter_asAllFertilePersonsForFertilityAlignmentByAge;
		logical	actorset_filter_asAllPerson;
		logical	actorset_filter_asAvailableMale;
		logical	actorset_filter_asAvailableMaleByType;
		logical	actorset_filter_asAvailableMaleByTypeImmiScratch;
		logical	actorset_filter_asAvailableMaleImmiScratch;
		logical	actorset_filter_asFemaleInUnionByAgeAndPartnerAge;
		logical	actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch;
		logical	actorset_filter_asPotentialImmigrantMothers;
		logical	actorset_filter_asPotentialSchoolOneStudent;
		logical	actorset_filter_asPotentialSchoolOneStudentToProcess;
		logical	actorset_filter_asPotentialSchoolTwoStudent;
		logical	actorset_filter_asPotentialSchoolTwoStudentToProcess;
		logical	actorset_filter_asResidentBabies;
		logical	actorset_filter_asResidentsAge0Sex;
		logical	actorset_filter_asResidentsAge0SexGeo;
		logical	actorset_filter_asSimBornAge0;
		logical	actorset_filter_asWantToMove;
		logical	actorset_filter_asWomenNoChildren;
		logical	actorset_filter_asWomenWithChildren;
		AGE_MIG	age_mig;
		ALL_YEAR_RANGE	all_year;
		CHILD_AGEGR	child_agegr;
		int	child_agegr_part;
		CHILD_MORTALITY_GROUP	child_mortality_group;
		EDUC1_GROUP	educ1_group;
		EDUC_ONE_GRADE	educ_one_grade;
		EDUC_ONE_LEVEL	educ_one_level;
		EDUC_TWO_LEVEL	educ_two_level;
		ETHNICITY_SHORT	ethnicity_short;
		FAM_ROLE	family_role;
		FEMALE_SPOUSE_AGE	female_spouse_age;
		FERTILE_AGE_RANGE	fertile_age;
		GEO_NAT	geo_nat;
		logical	has_passed_primary;
		logical	has_spouse;
		logical	in_hci_sample;
		logical	in_projected_time;
		logical	in_sel_years;
		logical	in_selected_year;
		logical	is_aded0;
		logical	is_aded1;
		logical	is_aded2;
		logical	is_adult;
		logical	is_alive_value_in;
		logical	is_child;
		logical	is_educ_one_student;
		logical	is_fertile;
		logical	is_mortal;
		logical	is_old;
		logical	is_resident;
		logical	is_stunted_value_in;
		logical	lives_with_dependent_child;
		MALE_SPOUSE_AGE	male_age;
		MALE_SPOUSE_AGE	male_spouse_age;
		MOTH_AGEGR	moth_agegr;
		int	moth_agegr_part;
		REGION_INT	region_birth;
		REGION_INT	region_int;
		REGION_NAT	region_nat;
		TAB_SEL_YEARS	sel_years;
		SELECTED_YEAR10	selected_year10;
		SIM_YEAR_RANGE	sim_year;
		SPOUSE_GROUP	spouse_group;
		logical	survived_early_years;
		logical	tab_out_of_primary;
		SIM_YEAR_RANGE	tab_sim_yob;
		logical	table_filter_PopPyramidByEduc;
		logical	table_filter_tabBirthsYearPlace;
		logical	table_filter_TabChildVaccination;
		logical	table_filter_TabEduc15ByDistrict;
		logical	table_filter_TabEduc15ByDistrictBirth;
		logical	table_filter_tabEducationFateGeobirYob;
		logical	table_filter_TabEducFateByGroup;
		logical	table_filter_TabEducFateDistrYob;
		logical	table_filter_tabHavingSpouse;
		logical	table_filter_tabHCI;
		logical	table_filter_TabHCIDistrict;
		logical	table_filter_TabImmunizationChildren;
		logical	table_filter_tabMigrationOriginDestination;
		logical	table_filter_TabPopProvAgeEducSex;
		logical	table_filter_tabPopulationYearPlace;
		logical	table_filter_TabPrenatCare;
		logical	table_filter_tabPreSchool;
		logical	table_filter_TabPrimarySchoolPlanning;
		logical	table_filter_TabPrimSchoolEntries;
		logical	table_filter_TabPrimSchoolGraduations;
		logical	table_filter_TabPrimSchoolOutOfSchool9to11;
		logical	table_filter_TabSchool2AttainmentsTab;
		logical	table_filter_TabSchool2TrackTab;
		logical	table_filter_tabStuntingSexRegMotherYob;
		logical	table_filter_tabStuntingSexRegYob;
		int	time_in_parity;
		UNION1_GROUP	union1_group;
		ALL_YEAR_RANGE	year_of_birth;
		SIM_YEAR_RANGE	year_of_first_immigration;
		int	years_of_schooling;
		counter	count_mlHHFatherChildren;
		counter	count_mlHHMotherChildren;
		EDUC_ONE_LEVEL	lSpouse_educ_one_fate;
		AGE_RANGE	lSpouse_integer_age;
		AGE_RANGE	min_over_mlHHMotherChildren_integer_age;
		integer	sum_over_mlHHFatherChildren_is_alive;
		integer	sum_over_mlHHMotherChildren_is_alive;
		logical	__finished;
		//   group=uint
		//   group=short
		short	children_in_household;
		//   group=ushort
		//   group=char
		//   group=uchar
		// end MembersDump

};

class PersonStateFunctions {
	public:
		PersonStateFunctions() { };
		class Person	*poParent;
		class PersonDerivedStates	*poDerivedStates;

		// state special functions
		double	Set_actor_subsample_weight( double value ){ return poParent->actor_subsample_weight = value; };
		double	Set_actor_weight( double value ){ return poParent->actor_weight = value; };
		TIME	Set_age( TIME value, bool bStartFunction = false );
		ALL_YEAR_RANGE	Set_calendar_year( ALL_YEAR_RANGE value );
		CREATION_TYPE	Set_creation_type( CREATION_TYPE value );
		EDUC_ONE_LEVEL	Set_educ_mother( EDUC_ONE_LEVEL value );
		EDUC_ONE_ENTRY_AGE	Set_educ_one_entry_age( EDUC_ONE_ENTRY_AGE value );
		EDUC_ONE_LEVEL	Set_educ_one_fate( EDUC_ONE_LEVEL value );
		EDUC_ONE_GRADE0	Set_educ_one_grade_attended( EDUC_ONE_GRADE0 value );
		EDUC_ONE_GRADE	Set_educ_one_grade_fate( EDUC_ONE_GRADE value );
		EDUC_ONE_GRADE0	Set_educ_one_grade_passed( EDUC_ONE_GRADE0 value );
		EDUC_ONE_STATUS	Set_educ_one_status( EDUC_ONE_STATUS value );
		logical	Set_educ_one_to_process( logical value );
		int	Set_educ_two_delay( int value );
		int	PlusPlusSet_educ_two_delay() { int old_value; old_value = poParent->educ_two_delay;
			Set_educ_two_delay((int) (poParent->educ_two_delay + 1)); return old_value; };
		EDUC_TWO_GRADE0	Set_educ_two_grade_attended( EDUC_TWO_GRADE0 value );
		EDUC_TWO_GRADE0	Set_educ_two_grade_passed( EDUC_TWO_GRADE0 value );
		EDUC_TWO_STATUS	Set_educ_two_status( EDUC_TWO_STATUS value );
		logical	Set_educ_two_to_process( logical value );
		ETHNICITY	Set_ethnicity( ETHNICITY value );
		logical	Set_ever_resident( logical value );
		logical	Set_ever_union( logical value );
		GEO	Set_geo( GEO value );
		GEO	Set_geo_birth( GEO value );
		GEO	Set_geo_prev( GEO value );
		GEO	Set_geo_want_to_move( GEO value );
		GOT_PRENAT_CARE	Set_got_prenat_care( GOT_PRENAT_CARE value );
		AGE_RANGE	Set_integer_age( AGE_RANGE value );
		AGE_RANGE	PlusPlusSet_integer_age() { AGE_RANGE old_value; old_value = poParent->integer_age;
			Set_integer_age((AGE_RANGE) (poParent->integer_age + 1)); return old_value; };
		logical	Set_is_alive( logical value );
		logical	Set_is_blocked_from_marriage( logical value );
		logical	Set_is_immunized( logical value );
		logical	Set_is_ready_for_birthtables( logical value );
		logical	Set_is_stunted( logical value );
		double	Set_mother_age_at_birth( double value );
		PARITY_RANGE	Set_parity( PARITY_RANGE value );
		logical	Set_preschool_is_decided( logical value );
		logical	Set_ready_to_set_alive( logical value );
		SEX	Set_sex( SEX value );
		counter	Set_table_unit_TabEducFateByGroup( counter value );
		counter	Set_table_unit_TabEducFateDistrYob( counter value );
		counter	Set_table_unit_TabPrimarySchoolPlanning( counter value );
		counter	Set_table_unit_TabSchool2TrackTab( counter value );
		counter	Set_table_unit_tabStuntingSexRegMotherYob( counter value );
		counter	Set_table_unit_tabStuntingSexRegYob( counter value );
		logical	Set_this_parity_spell( logical value );
		TIME	Set_time( TIME value, bool bStartFunction = false ){ return poParent->time = value; };
		TIME	Set_time_last_move( TIME value );
		TIME	Set_time_next_birthday( TIME value );
		TIME	Set_time_of_birth( TIME value );
		TIME	Set_time_of_first_immigration( TIME value );
		int	Set_years_preschool( int value );
		int	PlusPlusSet_years_preschool() { int old_value; old_value = poParent->years_preschool;
			Set_years_preschool((int) (poParent->years_preschool + 1)); return old_value; };
		logical	Set_actorset_filter_asAllFertilePersonsForFertilityAlignment( logical value );
		logical	Set_actorset_filter_asAllFertilePersonsForFertilityAlignmentByAge( logical value );
		logical	Set_actorset_filter_asAllPerson( logical value );
		logical	Set_actorset_filter_asAvailableMale( logical value );
		logical	Set_actorset_filter_asAvailableMaleByType( logical value );
		logical	Set_actorset_filter_asAvailableMaleByTypeImmiScratch( logical value );
		logical	Set_actorset_filter_asAvailableMaleImmiScratch( logical value );
		logical	Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAge( logical value );
		logical	Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch( logical value );
		logical	Set_actorset_filter_asPotentialImmigrantMothers( logical value );
		logical	Set_actorset_filter_asPotentialSchoolOneStudent( logical value );
		logical	Set_actorset_filter_asPotentialSchoolOneStudentToProcess( logical value );
		logical	Set_actorset_filter_asPotentialSchoolTwoStudent( logical value );
		logical	Set_actorset_filter_asPotentialSchoolTwoStudentToProcess( logical value );
		logical	Set_actorset_filter_asResidentBabies( logical value );
		logical	Set_actorset_filter_asResidentsAge0Sex( logical value );
		logical	Set_actorset_filter_asResidentsAge0SexGeo( logical value );
		logical	Set_actorset_filter_asSimBornAge0( logical value );
		logical	Set_actorset_filter_asWantToMove( logical value );
		logical	Set_actorset_filter_asWomenNoChildren( logical value );
		logical	Set_actorset_filter_asWomenWithChildren( logical value );
		double	Set_age_last_birth( double value );
		AGE_MIG	Set_age_mig( AGE_MIG value );
		double	Set_age_youngest_child( double value );
		ALL_YEAR_RANGE	Set_all_year( ALL_YEAR_RANGE value );
		CHILD_AGEGR	Set_child_agegr( CHILD_AGEGR value );
		int	Set_child_agegr_part( int value );
		double	Set_child_mortality( double value );
		CHILD_MORTALITY_GROUP	Set_child_mortality_group( CHILD_MORTALITY_GROUP value );
		EDUC1_GROUP	Set_educ1_group( EDUC1_GROUP value );
		EDUC_ONE_GRADE	Set_educ_one_grade( EDUC_ONE_GRADE value );
		EDUC_ONE_LEVEL	Set_educ_one_level( EDUC_ONE_LEVEL value );
		double	Set_educ_one_required_rooms( double value );
		double	Set_educ_one_required_rooms_value_out( double value );
		double	Set_educ_one_required_teachers( double value );
		double	Set_educ_one_required_teachers_value_out( double value );
		EDUC_TWO_LEVEL	Set_educ_two_level( EDUC_TWO_LEVEL value );
		ETHNICITY_SHORT	Set_ethnicity_short( ETHNICITY_SHORT value );
		FAM_ROLE	Set_family_role( FAM_ROLE value );
		FEMALE_SPOUSE_AGE	Set_female_spouse_age( FEMALE_SPOUSE_AGE value );
		FERTILE_AGE_RANGE	Set_fertile_age( FERTILE_AGE_RANGE value );
		GEO_NAT	Set_geo_nat( GEO_NAT value );
		logical	Set_has_passed_primary( logical value );
		logical	Set_has_spouse( logical value );
		logical	Set_in_hci_sample( logical value );
		logical	Set_in_projected_time( logical value );
		logical	Set_in_sel_years( logical value );
		logical	Set_in_selected_year( logical value );
		logical	Set_is_aded0( logical value );
		logical	Set_is_aded1( logical value );
		logical	Set_is_aded2( logical value );
		logical	Set_is_alive_value_in( logical value );
		logical	Set_is_child( logical value );
		logical	Set_is_educ_one_student( logical value );
		logical	Set_is_fertile( logical value );
		logical	Set_is_mortal( logical value );
		logical	Set_is_old( logical value );
		logical	Set_is_resident( logical value );
		logical	Set_is_stunted_value_in( logical value );
		logical	Set_lives_with_dependent_child( logical value );
		MALE_SPOUSE_AGE	Set_male_age( MALE_SPOUSE_AGE value );
		MALE_SPOUSE_AGE	Set_male_spouse_age( MALE_SPOUSE_AGE value );
		MOTH_AGEGR	Set_moth_agegr( MOTH_AGEGR value );
		int	Set_moth_agegr_part( int value );
		REGION_INT	Set_region_birth( REGION_INT value );
		REGION_INT	Set_region_int( REGION_INT value );
		REGION_NAT	Set_region_nat( REGION_NAT value );
		TAB_SEL_YEARS	Set_sel_years( TAB_SEL_YEARS value );
		SELECTED_YEAR10	Set_selected_year10( SELECTED_YEAR10 value );
		SIM_YEAR_RANGE	Set_sim_year( SIM_YEAR_RANGE value );
		SPOUSE_GROUP	Set_spouse_group( SPOUSE_GROUP value );
		logical	Set_survived_early_years( logical value );
		logical	Set_tab_out_of_primary( logical value );
		SIM_YEAR_RANGE	Set_tab_sim_yob( SIM_YEAR_RANGE value );
		logical	Set_table_filter_PopPyramidByEduc( logical value );
		logical	Set_table_filter_tabBirthsYearPlace( logical value );
		logical	Set_table_filter_TabChildVaccination( logical value );
		logical	Set_table_filter_TabEduc15ByDistrict( logical value );
		logical	Set_table_filter_TabEduc15ByDistrictBirth( logical value );
		logical	Set_table_filter_tabEducationFateGeobirYob( logical value );
		logical	Set_table_filter_TabEducFateByGroup( logical value );
		logical	Set_table_filter_TabEducFateDistrYob( logical value );
		logical	Set_table_filter_tabHavingSpouse( logical value );
		logical	Set_table_filter_tabHCI( logical value );
		logical	Set_table_filter_TabHCIDistrict( logical value );
		logical	Set_table_filter_TabImmunizationChildren( logical value );
		logical	Set_table_filter_tabMigrationOriginDestination( logical value );
		logical	Set_table_filter_TabPopProvAgeEducSex( logical value );
		logical	Set_table_filter_tabPopulationYearPlace( logical value );
		logical	Set_table_filter_TabPrenatCare( logical value );
		logical	Set_table_filter_tabPreSchool( logical value );
		logical	Set_table_filter_TabPrimarySchoolPlanning( logical value );
		logical	Set_table_filter_TabPrimSchoolEntries( logical value );
		logical	Set_table_filter_TabPrimSchoolGraduations( logical value );
		logical	Set_table_filter_TabPrimSchoolOutOfSchool9to11( logical value );
		logical	Set_table_filter_TabSchool2AttainmentsTab( logical value );
		logical	Set_table_filter_TabSchool2TrackTab( logical value );
		logical	Set_table_filter_tabStuntingSexRegMotherYob( logical value );
		logical	Set_table_filter_tabStuntingSexRegYob( logical value );
		int	Set_time_in_parity( int value );
		UNION1_GROUP	Set_union1_group( UNION1_GROUP value );
		ALL_YEAR_RANGE	Set_year_of_birth( ALL_YEAR_RANGE value );
		SIM_YEAR_RANGE	Set_year_of_first_immigration( SIM_YEAR_RANGE value );
		counter	Set_count_mlHHFatherChildren( counter value );
		counter	Set_count_mlHHMotherChildren( counter value );
		double	Set_lCalibrator_mort_female_0( double value );
		double	Set_lCalibrator_mort_female_1( double value );
		double	Set_lCalibrator_mort_female_2( double value );
		double	Set_lCalibrator_mort_female_3( double value );
		double	Set_lCalibrator_mort_female_4( double value );
		double	Set_lCalibrator_mort_male_0( double value );
		double	Set_lCalibrator_mort_male_1( double value );
		double	Set_lCalibrator_mort_male_2( double value );
		double	Set_lCalibrator_mort_male_3( double value );
		double	Set_lCalibrator_mort_male_4( double value );
		EDUC_ONE_LEVEL	Set_lSpouse_educ_one_fate( EDUC_ONE_LEVEL value );
		AGE_RANGE	Set_lSpouse_integer_age( AGE_RANGE value );
		AGE_RANGE	Set_min_over_mlHHMotherChildren_integer_age( AGE_RANGE value );
		integer	Set_sum_over_mlHHFatherChildren_is_alive( integer value );
		integer	Set_sum_over_mlHHMotherChildren_is_alive( integer value );
};

class PersonDerivedStates {
	public:
		PersonDerivedStates() { };
		class Person	*poParent;
		class PersonStateFunctions	*poStateFunctions;

		// derived state variables and their special functions
		TIME	Get_active_spell_duration_this_parity_spell_TRUE() { return active_spell_duration_this_parity_spell_TRUE + ( ( poParent->this_parity_spell == TRUE ) ? ( poParent->__time - active_spell_duration_this_parity_spell_TRUE_time ) : 0 ); };
		ETHNICITY_SHORT	Set_aggregate_ethnicity_ETHNICITY_SHORT( ETHNICITY_SHORT value );
		GEO_NAT	Set_aggregate_geo_birth_GEO_NAT( GEO_NAT value );
		REGION_INT	Set_aggregate_geo_birth_REGION_INT( REGION_INT value );
		GEO_NAT	Set_aggregate_geo_GEO_NAT( GEO_NAT value );
		REGION_INT	Set_aggregate_geo_REGION_INT( REGION_INT value );
		REGION_NAT	Set_aggregate_region_birth_REGION_NAT( REGION_NAT value );
		REGION_NAT	Set_aggregate_region_int_REGION_NAT( REGION_NAT value );
		TIME	Set_duration_( TIME value );
		TIME	Get_duration_() { return duration_ + poParent->__time - duration__time; };
		TIME	Set_duration_educ_one_fate_EOL_HIGH( TIME value );
		TIME	Get_duration_educ_one_fate_EOL_HIGH() { return duration_educ_one_fate_EOL_HIGH + ( ( poParent->educ_one_fate == EOL_HIGH ) ? ( poParent->__time - duration_educ_one_fate_EOL_HIGH_time ) : 0 ); };
		TIME	Set_duration_educ_one_fate_EOL_LOW( TIME value );
		TIME	Get_duration_educ_one_fate_EOL_LOW() { return duration_educ_one_fate_EOL_LOW + ( ( poParent->educ_one_fate == EOL_LOW ) ? ( poParent->__time - duration_educ_one_fate_EOL_LOW_time ) : 0 ); };
		TIME	Set_duration_educ_one_fate_EOL_MEDIUM( TIME value );
		TIME	Get_duration_educ_one_fate_EOL_MEDIUM() { return duration_educ_one_fate_EOL_MEDIUM + ( ( poParent->educ_one_fate == EOL_MEDIUM ) ? ( poParent->__time - duration_educ_one_fate_EOL_MEDIUM_time ) : 0 ); };
		TIME	Set_duration_educ_one_level_EOL_HIGH( TIME value );
		TIME	Get_duration_educ_one_level_EOL_HIGH() { return duration_educ_one_level_EOL_HIGH + ( ( poParent->educ_one_level == EOL_HIGH ) ? ( poParent->__time - duration_educ_one_level_EOL_HIGH_time ) : 0 ); };
		TIME	Set_duration_educ_one_level_EOL_LOW( TIME value );
		TIME	Get_duration_educ_one_level_EOL_LOW() { return duration_educ_one_level_EOL_LOW + ( ( poParent->educ_one_level == EOL_LOW ) ? ( poParent->__time - duration_educ_one_level_EOL_LOW_time ) : 0 ); };
		TIME	Set_duration_educ_one_level_EOL_MEDIUM( TIME value );
		TIME	Get_duration_educ_one_level_EOL_MEDIUM() { return duration_educ_one_level_EOL_MEDIUM + ( ( poParent->educ_one_level == EOL_MEDIUM ) ? ( poParent->__time - duration_educ_one_level_EOL_MEDIUM_time ) : 0 ); };
		TIME	Set_duration_educ_two_level_ETL_GRADUATE( TIME value );
		TIME	Get_duration_educ_two_level_ETL_GRADUATE() { return duration_educ_two_level_ETL_GRADUATE + ( ( poParent->educ_two_level == ETL_GRADUATE ) ? ( poParent->__time - duration_educ_two_level_ETL_GRADUATE_time ) : 0 ); };
		TIME	Set_duration_educ_two_level_ETL_INCOMPLETE( TIME value );
		TIME	Get_duration_educ_two_level_ETL_INCOMPLETE() { return duration_educ_two_level_ETL_INCOMPLETE + ( ( poParent->educ_two_level == ETL_INCOMPLETE ) ? ( poParent->__time - duration_educ_two_level_ETL_INCOMPLETE_time ) : 0 ); };
		TIME	Set_duration_educ_two_level_ETL_NEVER( TIME value );
		TIME	Get_duration_educ_two_level_ETL_NEVER() { return duration_educ_two_level_ETL_NEVER + ( ( poParent->educ_two_level == ETL_NEVER ) ? ( poParent->__time - duration_educ_two_level_ETL_NEVER_time ) : 0 ); };
		TIME	Set_duration_ever_union_TRUE( TIME value );
		TIME	Get_duration_ever_union_TRUE() { return duration_ever_union_TRUE + ( ( poParent->ever_union == TRUE ) ? ( poParent->__time - duration_ever_union_TRUE_time ) : 0 ); };
		TIME	Set_duration_got_prenat_care_GPC_YES( TIME value );
		TIME	Get_duration_got_prenat_care_GPC_YES() { return duration_got_prenat_care_GPC_YES + ( ( poParent->got_prenat_care == GPC_YES ) ? ( poParent->__time - duration_got_prenat_care_GPC_YES_time ) : 0 ); };
		TIME	Set_duration_has_spouse_TRUE( TIME value );
		TIME	Get_duration_has_spouse_TRUE() { return duration_has_spouse_TRUE + ( ( poParent->has_spouse == TRUE ) ? ( poParent->__time - duration_has_spouse_TRUE_time ) : 0 ); };
		TIME	Set_duration_is_aded0_TRUE( TIME value );
		TIME	Get_duration_is_aded0_TRUE() { return duration_is_aded0_TRUE + ( ( poParent->is_aded0 == TRUE ) ? ( poParent->__time - duration_is_aded0_TRUE_time ) : 0 ); };
		TIME	Set_duration_is_aded1_TRUE( TIME value );
		TIME	Get_duration_is_aded1_TRUE() { return duration_is_aded1_TRUE + ( ( poParent->is_aded1 == TRUE ) ? ( poParent->__time - duration_is_aded1_TRUE_time ) : 0 ); };
		TIME	Set_duration_is_aded2_TRUE( TIME value );
		TIME	Get_duration_is_aded2_TRUE() { return duration_is_aded2_TRUE + ( ( poParent->is_aded2 == TRUE ) ? ( poParent->__time - duration_is_aded2_TRUE_time ) : 0 ); };
		TIME	Set_duration_is_child_TRUE( TIME value );
		TIME	Get_duration_is_child_TRUE() { return duration_is_child_TRUE + ( ( poParent->is_child == TRUE ) ? ( poParent->__time - duration_is_child_TRUE_time ) : 0 ); };
		TIME	Set_duration_is_immunized_TRUE( TIME value );
		TIME	Get_duration_is_immunized_TRUE() { return duration_is_immunized_TRUE + ( ( poParent->is_immunized == TRUE ) ? ( poParent->__time - duration_is_immunized_TRUE_time ) : 0 ); };
		TIME	Set_duration_is_old_TRUE( TIME value );
		TIME	Get_duration_is_old_TRUE() { return duration_is_old_TRUE + ( ( poParent->is_old == TRUE ) ? ( poParent->__time - duration_is_old_TRUE_time ) : 0 ); };
		TIME	Set_duration_tab_out_of_primary_TRUE( TIME value );
		TIME	Get_duration_tab_out_of_primary_TRUE() { return duration_tab_out_of_primary_TRUE + ( ( poParent->tab_out_of_primary == TRUE ) ? ( poParent->__time - duration_tab_out_of_primary_TRUE_time ) : 0 ); };
		TIME	Set_duration_years_preschool_0( TIME value );
		TIME	Get_duration_years_preschool_0() { return duration_years_preschool_0 + ( ( poParent->years_preschool == 0 ) ? ( poParent->__time - duration_years_preschool_0_time ) : 0 ); };
		TIME	Set_duration_years_preschool_1( TIME value );
		TIME	Get_duration_years_preschool_1() { return duration_years_preschool_1 + ( ( poParent->years_preschool == 1 ) ? ( poParent->__time - duration_years_preschool_1_time ) : 0 ); };
		TIME	Set_duration_years_preschool_2( TIME value );
		TIME	Get_duration_years_preschool_2() { return duration_years_preschool_2 + ( ( poParent->years_preschool == 2 ) ? ( poParent->__time - duration_years_preschool_2_time ) : 0 ); };
		counter	Set_entrances_educ_one_grade_attended_1( counter value );
		counter	Set_entrances_has_passed_primary_TRUE( counter value );
		MOTH_AGEGR_PART	Set_split_age_last_birth_MOTH_AGEGR_PART( MOTH_AGEGR_PART value );
		CHILD_AGEGR_PART	Set_split_age_youngest_child_CHILD_AGEGR_PART( CHILD_AGEGR_PART value );
		TAB_AGEC5	Set_split_integer_age_TAB_AGEC5( TAB_AGEC5 value );
		DUR_TIME_IN_PARITY	Set_sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY( DUR_TIME_IN_PARITY value );
		AGE_MIG_PART	Set_sssplit_age_AGE_MIG_PART( AGE_MIG_PART value );
		counter	Set_transitions_integer_age_14_15( counter value );
		counter	Set_transitions_is_alive_FALSE_TRUE( counter value );
		counter	Set_transitions_is_alive_TRUE_FALSE( counter value );
		counter	Set_transitions_survived_early_years_FALSE_TRUE( counter value );
		logical	Set_trigger_changes_calendar_year( logical value );
		logical	Set_trigger_changes_geo( logical value );
		logical	Set_trigger_entrances_is_ready_for_birthtables_TRUE( logical value );
		real	Set_value_at_transitions_is_alive_TRUE_FALSE_adult_survival( real value );
		real	Set_value_at_transitions_is_alive_TRUE_FALSE_ind_hci( real value );
		integer	Set_value_at_transitions_is_alive_TRUE_FALSE_is_stunted( integer value );
		real	Set_value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling( real value );
		integer	Set_value_at_transitions_is_alive_TRUE_FALSE_survived_early_years( integer value );
		integer	Set_value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling( integer value );

		// begin MembersDump
		//   group=<unknown token type>
		//   group=double
		TIME	active_spell_duration_this_parity_spell_TRUE;
		TIME	active_spell_duration_this_parity_spell_TRUE_time;
		TIME	completed_spell_duration_this_parity_spell_TRUE;
		TIME	duration_;
		TIME	duration__time;
		TIME	duration_educ_one_fate_EOL_HIGH;
		TIME	duration_educ_one_fate_EOL_HIGH_time;
		TIME	duration_educ_one_fate_EOL_LOW;
		TIME	duration_educ_one_fate_EOL_LOW_time;
		TIME	duration_educ_one_fate_EOL_MEDIUM;
		TIME	duration_educ_one_fate_EOL_MEDIUM_time;
		TIME	duration_educ_one_level_EOL_HIGH;
		TIME	duration_educ_one_level_EOL_HIGH_time;
		TIME	duration_educ_one_level_EOL_LOW;
		TIME	duration_educ_one_level_EOL_LOW_time;
		TIME	duration_educ_one_level_EOL_MEDIUM;
		TIME	duration_educ_one_level_EOL_MEDIUM_time;
		TIME	duration_educ_two_level_ETL_GRADUATE;
		TIME	duration_educ_two_level_ETL_GRADUATE_time;
		TIME	duration_educ_two_level_ETL_INCOMPLETE;
		TIME	duration_educ_two_level_ETL_INCOMPLETE_time;
		TIME	duration_educ_two_level_ETL_NEVER;
		TIME	duration_educ_two_level_ETL_NEVER_time;
		TIME	duration_ever_union_TRUE;
		TIME	duration_ever_union_TRUE_time;
		TIME	duration_got_prenat_care_GPC_YES;
		TIME	duration_got_prenat_care_GPC_YES_time;
		TIME	duration_has_spouse_TRUE;
		TIME	duration_has_spouse_TRUE_time;
		TIME	duration_is_aded0_TRUE;
		TIME	duration_is_aded0_TRUE_time;
		TIME	duration_is_aded1_TRUE;
		TIME	duration_is_aded1_TRUE_time;
		TIME	duration_is_aded2_TRUE;
		TIME	duration_is_aded2_TRUE_time;
		TIME	duration_is_child_TRUE;
		TIME	duration_is_child_TRUE_time;
		TIME	duration_is_immunized_TRUE;
		TIME	duration_is_immunized_TRUE_time;
		TIME	duration_is_old_TRUE;
		TIME	duration_is_old_TRUE_time;
		TIME	duration_tab_out_of_primary_TRUE;
		TIME	duration_tab_out_of_primary_TRUE_time;
		TIME	duration_years_preschool_0;
		TIME	duration_years_preschool_0_time;
		TIME	duration_years_preschool_1;
		TIME	duration_years_preschool_1_time;
		TIME	duration_years_preschool_2;
		TIME	duration_years_preschool_2_time;
		TIME	sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time;
		TIME	sssplit_age_AGE_MIG_PART_time;
		TIME	trigger_changes_calendar_year_time;
		TIME	trigger_changes_geo_time;
		TIME	trigger_entrances_is_ready_for_birthtables_TRUE_time;
		real	value_at_transitions_is_alive_TRUE_FALSE_adult_survival;
		real	value_at_transitions_is_alive_TRUE_FALSE_ind_hci;
		real	value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling;
		//   group=float
		//   group=long
		//   group=ulong
		//   group=int
		ETHNICITY_SHORT	aggregate_ethnicity_ETHNICITY_SHORT;
		GEO_NAT	aggregate_geo_birth_GEO_NAT;
		REGION_INT	aggregate_geo_birth_REGION_INT;
		GEO_NAT	aggregate_geo_GEO_NAT;
		REGION_INT	aggregate_geo_REGION_INT;
		REGION_NAT	aggregate_region_birth_REGION_NAT;
		REGION_NAT	aggregate_region_int_REGION_NAT;
		counter	entrances_educ_one_grade_attended_1;
		counter	entrances_has_passed_primary_TRUE;
		MOTH_AGEGR_PART	split_age_last_birth_MOTH_AGEGR_PART;
		CHILD_AGEGR_PART	split_age_youngest_child_CHILD_AGEGR_PART;
		TAB_AGEC5	split_integer_age_TAB_AGEC5;
		DUR_TIME_IN_PARITY	sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY;
		AGE_MIG_PART	sssplit_age_AGE_MIG_PART;
		counter	transitions_integer_age_14_15;
		counter	transitions_is_alive_FALSE_TRUE;
		counter	transitions_is_alive_TRUE_FALSE;
		counter	transitions_survived_early_years_FALSE_TRUE;
		logical	trigger_changes_calendar_year;
		logical	trigger_changes_geo;
		logical	trigger_entrances_is_ready_for_birthtables_TRUE;
		integer	value_at_transitions_is_alive_TRUE_FALSE_is_stunted;
		integer	value_at_transitions_is_alive_TRUE_FALSE_survived_early_years;
		integer	value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling;
		//   group=uint
		//   group=short
		//   group=ushort
		//   group=char
		//   group=uchar
		// end MembersDump

};

class StartpopValues : public ActorClass {
	public:
		TIME	__time;
		double	__events;
		class StartpopValuesStateFunctions	*poStateFunctions;
		class StartpopValuesDerivedStates	*poDerivedStates;
		void *operator new( size_t count );

		// links to other actors
		class Person	*lPersonStartpop;
		void	Set_lPersonStartpop( Person *value );

		// user-supplied data members
		double	StartPopValue[rc_size_PERSON_MICRODATA_COLUMNS];

		// user-supplied and required functions
		StartpopValues() { };
		void Finish ( ) ; 
		void Start ( TIME dTime ) ; 

		// special functions generated by the pre-compiler
		StartpopValues( BOOL bInit ) { InitActor(); };
		~StartpopValues() { DeleteActor(); };
		void	InitActor();
		void	DeleteActor();
		void	InitializeStates();
		void	UpdateLinkedStates();
		void	StartSpecial1();
		void	StartSpecial2();
		void	FinishSpecial();
		void	WaitUntil( double target_time );
		StartpopValues	*UpdateTime( );
		void	WaitUntilThisActor(double target_time);
		void	WaitSpecial( TIME wait_time );

		// event functions and variables
		void __fastcall	FirstMoveEvent();
		TIME __fastcall	timeFirstMoveEvent();
		Event	*poFirstMoveEvent;
		void __fastcall	FirstUnionImputationEvent();
		TIME __fastcall	timeFirstUnionImputationEvent();
		Event	*poFirstUnionImputationEvent;
		void __fastcall	LastBirthImputationEvent();
		TIME __fastcall	timeLastBirthImputationEvent();
		Event	*poLastBirthImputationEvent;
		void __fastcall	LastMoveEvent();
		TIME __fastcall	timeLastMoveEvent();
		Event	*poLastMoveEvent;
		void __fastcall	SetFirstUnionClock();
		TIME __fastcall	timeSetFirstUnionClock();
		Event	*poSetFirstUnionClock;
		void __fastcall	SetLastBirthClock();
		TIME __fastcall	timeSetLastBirthClock();
		Event	*poSetLastBirthClock;
		void __fastcall	SetResidentialMoveClock();
		TIME __fastcall	timeSetResidentialMoveClock();
		Event	*poSetResidentialMoveClock;
		void __fastcall	SuicideEvent();
		TIME __fastcall	timeSuicideEvent();
		Event	*poSuicideEvent;

		// event processing functions
		double	EventTime( int nEventNum, int *pnEventInfo );
		void	Implement( int nEventNum, int nEventInfo );

		// state variables and their special functions
		long	Get_actor_id() { return actor_id; };
		double	Set_actor_subsample_weight( double value );
		double	Set_actor_weight( double value );

		// begin MembersDump
		//   group=<unknown token type>
		//   group=double
		double	actor_subsample_weight;
		double	actor_weight;
		TIME	age;
		double	case_seed;
		TIME	scheduled_time_first_move;
		TIME	scheduled_time_last_move;
		TIME	time;
		TIME	time_first_union_imputation;
		TIME	time_last_birth_imputation;
		//   group=float
		//   group=long
		long	actor_id;
		long	case_id;
		//   group=ulong
		//   group=int
		logical	is_activated;
		logical	is_set_first_union_clock;
		logical	is_set_last_birth_clock;
		logical	is_set_residential_move_clock;
		logical	is_removeable;
		logical	lPersonStartpop_is_mortal;
		logical	__finished;
		//   group=uint
		//   group=short
		//   group=ushort
		//   group=char
		//   group=uchar
		// end MembersDump

};

class StartpopValuesStateFunctions {
	public:
		StartpopValuesStateFunctions() { };
		class StartpopValues	*poParent;
		class StartpopValuesDerivedStates	*poDerivedStates;

		// state special functions
		double	Set_actor_subsample_weight( double value ){ return poParent->actor_subsample_weight = value; };
		double	Set_actor_weight( double value ){ return poParent->actor_weight = value; };
		TIME	Set_age( TIME value, bool bStartFunction = false ){ return poParent->age = value; };
		logical	Set_is_activated( logical value );
		logical	Set_is_set_first_union_clock( logical value );
		logical	Set_is_set_last_birth_clock( logical value );
		logical	Set_is_set_residential_move_clock( logical value );
		TIME	Set_scheduled_time_first_move( TIME value );
		TIME	Set_scheduled_time_last_move( TIME value );
		TIME	Set_time( TIME value, bool bStartFunction = false ){ return poParent->time = value; };
		TIME	Set_time_first_union_imputation( TIME value );
		TIME	Set_time_last_birth_imputation( TIME value );
		logical	Set_is_removeable( logical value );
		logical	Set_lPersonStartpop_is_mortal( logical value );
};

class StartpopValuesDerivedStates {
	public:
		StartpopValuesDerivedStates() { };
		class StartpopValues	*poParent;
		class StartpopValuesStateFunctions	*poStateFunctions;

		// derived state variables and their special functions

		// begin MembersDump
		//   group=<unknown token type>
		//   group=double
		//   group=float
		//   group=long
		//   group=ulong
		//   group=int
		//   group=uint
		//   group=short
		//   group=ushort
		//   group=char
		//   group=uchar
		// end MembersDump

};

// global actor declarations
extern Thread Calibrator	*gpoEmptyCalibrator;
typedef CTypedPtrArray<CObArray, Calibrator*> CalibratorArray;
void	DeleteAllCalibratorActors();
extern Thread Clock	*gpoEmptyClock;
typedef CTypedPtrArray<CObArray, Clock*> ClockArray;
void	DeleteAllClockActors();
extern Thread Observation	*gpoEmptyObservation;
typedef CTypedPtrArray<CObArray, Observation*> ObservationArray;
void	DeleteAllObservationActors();
extern Thread Person	*gpoEmptyPerson;
typedef CTypedPtrArray<CObArray, Person*> PersonArray;
void	DeleteAllPersonActors();
extern Thread StartpopValues	*gpoEmptyStartpopValues;
typedef CTypedPtrArray<CObArray, StartpopValues*> StartpopValuesArray;
void	DeleteAllStartpopValuesActors();

class Calibrator_mlCalibratorToPerson : public CObject {
	public:
		Calibrator_mlCalibratorToPerson( Calibrator *pObject ) { prParentCalibrator = pObject; };
		void	Add( Person *prPerson );
		void	Remove( Person *prPerson );
		void	RemoveAll();
		void	Finish();
		void	FinishAll();
		Person	*GetNext( int nInitPos, int *pnPos );
		void	UpdateReferences();
		Calibrator	*prParentCalibrator;
		PersonArray	arData;
};

class Person_mlBioFatherChildren : public CObject {
	public:
		Person_mlBioFatherChildren( Person *pObject ) { prParentPerson = pObject; };
		void	Add( Person *prPerson );
		void	Remove( Person *prPerson );
		void	RemoveAll();
		void	Finish();
		void	FinishAll();
		Person	*GetNext( int nInitPos, int *pnPos );
		void	UpdateReferences();
		Person	*prParentPerson;
		PersonArray	arData;
};

class Person_mlBioMotherChildren : public CObject {
	public:
		Person_mlBioMotherChildren( Person *pObject ) { prParentPerson = pObject; };
		void	Add( Person *prPerson );
		void	Remove( Person *prPerson );
		void	RemoveAll();
		void	Finish();
		void	FinishAll();
		Person	*GetNext( int nInitPos, int *pnPos );
		void	UpdateReferences();
		Person	*prParentPerson;
		PersonArray	arData;
};

class Person_mlHHFatherChildren : public CObject {
	public:
		Person_mlHHFatherChildren( Person *pObject ) { prParentPerson = pObject; };
		void	Add( Person *prPerson );
		void	Remove( Person *prPerson );
		void	RemoveAll();
		void	Finish();
		void	FinishAll();
		Person	*GetNext( int nInitPos, int *pnPos );
		counter	Count();
		integer	SumOver_is_alive();
		void	UpdateReferences();
		Person	*prParentPerson;
		PersonArray	arData;
};

class Person_mlHHMotherChildren : public CObject {
	public:
		Person_mlHHMotherChildren( Person *pObject ) { prParentPerson = pObject; };
		void	Add( Person *prPerson );
		void	Remove( Person *prPerson );
		void	RemoveAll();
		void	Finish();
		void	FinishAll();
		Person	*GetNext( int nInitPos, int *pnPos );
		counter	Count();
		AGE_RANGE	MinOver_integer_age( Person *prSpecialPerson, AGE_RANGE value );
		integer	SumOver_is_alive();
		void	UpdateReferences();
		Person	*prParentPerson;
		PersonArray	arData;
};

class asAllFertilePersonsForFertilityAlignmentActorSetTree : public ActorSetTree	{
	public:
		asAllFertilePersonsForFertilityAlignmentActorSetTree() {m_nId = 0;};
		~asAllFertilePersonsForFertilityAlignmentActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asAllFertilePersonsForFertilityAlignmentActorSetTree	*asAllFertilePersonsForFertilityAlignment[2];

class asAllFertilePersonsForFertilityAlignmentByAgeActorSetTree : public ActorSetTree	{
	public:
		asAllFertilePersonsForFertilityAlignmentByAgeActorSetTree() {m_nId = 1;};
		~asAllFertilePersonsForFertilityAlignmentByAgeActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asAllFertilePersonsForFertilityAlignmentByAgeActorSetTree	*asAllFertilePersonsForFertilityAlignmentByAge[rc_size_FERTILE_AGE_RANGE][2];

class asAllPersonActorSetTree : public ActorSetTree	{
	public:
		asAllPersonActorSetTree() {m_nId = 2;};
		~asAllPersonActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asAllPersonActorSetTree	*asAllPerson;

class asAvailableMaleActorSetTree : public ActorSetTree	{
	public:
		asAvailableMaleActorSetTree() {m_nId = 3;};
		~asAvailableMaleActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asAvailableMaleActorSetTree	*asAvailableMale[rc_size_MALE_SPOUSE_AGE][rc_size_GEO];

class asAvailableMaleByTypeActorSetTree : public ActorSetTree	{
	public:
		asAvailableMaleByTypeActorSetTree() {m_nId = 4;};
		~asAvailableMaleByTypeActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asAvailableMaleByTypeActorSetTree	*asAvailableMaleByType[rc_size_MALE_SPOUSE_AGE][rc_size_SPOUSE_GROUP][rc_size_GEO];

class asAvailableMaleByTypeImmiScratchActorSetTree : public ActorSetTree	{
	public:
		asAvailableMaleByTypeImmiScratchActorSetTree() {m_nId = 5;};
		~asAvailableMaleByTypeImmiScratchActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asAvailableMaleByTypeImmiScratchActorSetTree	*asAvailableMaleByTypeImmiScratch[rc_size_MALE_SPOUSE_AGE][rc_size_SPOUSE_GROUP][rc_size_GEO][rc_size_SIM_YEAR_RANGE];

class asAvailableMaleImmiScratchActorSetTree : public ActorSetTree	{
	public:
		asAvailableMaleImmiScratchActorSetTree() {m_nId = 6;};
		~asAvailableMaleImmiScratchActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asAvailableMaleImmiScratchActorSetTree	*asAvailableMaleImmiScratch[rc_size_MALE_SPOUSE_AGE][rc_size_GEO][rc_size_SIM_YEAR_RANGE];

class asCalibratorActorSetTree : public ActorSetTree	{
	public:
		asCalibratorActorSetTree() {m_nId = 7;};
		~asCalibratorActorSetTree() {};
		Calibrator	* Item(int nIndex) {return (Calibrator *) ActorSetTree::Item(nIndex);};
		Calibrator	* GetRandom(double dRandom) {return (Calibrator *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asCalibratorActorSetTree	*asCalibrator;

class asClockActorSetTree : public ActorSetTree	{
	public:
		asClockActorSetTree() {m_nId = 8;};
		~asClockActorSetTree() {};
		Clock	* Item(int nIndex) {return (Clock *) ActorSetTree::Item(nIndex);};
		Clock	* GetRandom(double dRandom) {return (Clock *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asClockActorSetTree	*asClock;

class asFemaleInUnionByAgeAndPartnerAgeActorSetTree : public ActorSetTree	{
	public:
		asFemaleInUnionByAgeAndPartnerAgeActorSetTree() {m_nId = 9;};
		~asFemaleInUnionByAgeAndPartnerAgeActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asFemaleInUnionByAgeAndPartnerAgeActorSetTree	*asFemaleInUnionByAgeAndPartnerAge[rc_size_FEMALE_SPOUSE_AGE][rc_size_MALE_SPOUSE_AGE];

class asFemaleInUnionByAgeAndPartnerAgeImmiScratchActorSetTree : public ActorSetTree	{
	public:
		asFemaleInUnionByAgeAndPartnerAgeImmiScratchActorSetTree() {m_nId = 10;};
		~asFemaleInUnionByAgeAndPartnerAgeImmiScratchActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asFemaleInUnionByAgeAndPartnerAgeImmiScratchActorSetTree	*asFemaleInUnionByAgeAndPartnerAgeImmiScratch[rc_size_FEMALE_SPOUSE_AGE][rc_size_MALE_SPOUSE_AGE][rc_size_SIM_YEAR_RANGE];

class asObservationAllActorSetTree : public ActorSetTree	{
	public:
		asObservationAllActorSetTree() {m_nId = 11;};
		~asObservationAllActorSetTree() {};
		Observation	* Item(int nIndex) {return (Observation *) ActorSetTree::Item(nIndex);};
		Observation	* GetRandom(double dRandom) {return (Observation *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asObservationAllActorSetTree	*asObservationAll[rc_size_POP_POOL];

class asObservationByFamActorSetTree : public ActorSetTree	{
	public:
		asObservationByFamActorSetTree() {m_nId = 12;};
		~asObservationByFamActorSetTree() {};
		Observation	* Item(int nIndex) {return (Observation *) ActorSetTree::Item(nIndex);};
		Observation	* GetRandom(double dRandom) {return (Observation *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asObservationByFamActorSetTree	*asObservationByFam[rc_size_FAM_ID];

class asObservationByFamOldestActorSetTree : public ActorSetTree	{
	public:
		asObservationByFamOldestActorSetTree() {m_nId = 13;};
		~asObservationByFamOldestActorSetTree() {};
		Observation	* Item(int nIndex) {return (Observation *) ActorSetTree::Item(nIndex);};
		Observation	* GetRandom(double dRandom) {return (Observation *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asObservationByFamOldestActorSetTree	*asObservationByFamOldest[rc_size_FAM_ID];

class asObservationHeadsActorSetTree : public ActorSetTree	{
	public:
		asObservationHeadsActorSetTree() {m_nId = 14;};
		~asObservationHeadsActorSetTree() {};
		Observation	* Item(int nIndex) {return (Observation *) ActorSetTree::Item(nIndex);};
		Observation	* GetRandom(double dRandom) {return (Observation *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asObservationHeadsActorSetTree	*asObservationHeads[rc_size_POP_POOL];

class asObservationNonHeadsActorSetTree : public ActorSetTree	{
	public:
		asObservationNonHeadsActorSetTree() {m_nId = 15;};
		~asObservationNonHeadsActorSetTree() {};
		Observation	* Item(int nIndex) {return (Observation *) ActorSetTree::Item(nIndex);};
		Observation	* GetRandom(double dRandom) {return (Observation *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asObservationNonHeadsActorSetTree	*asObservationNonHeads[rc_size_FAM_ID][rc_size_POP_POOL];

class asObservationsActorSetTree : public ActorSetTree	{
	public:
		asObservationsActorSetTree() {m_nId = 16;};
		~asObservationsActorSetTree() {};
		Observation	* Item(int nIndex) {return (Observation *) ActorSetTree::Item(nIndex);};
		Observation	* GetRandom(double dRandom) {return (Observation *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asObservationsActorSetTree	*asObservations;

class asPotentialImmigrantMothersActorSetTree : public ActorSetTree	{
	public:
		asPotentialImmigrantMothersActorSetTree() {m_nId = 17;};
		~asPotentialImmigrantMothersActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asPotentialImmigrantMothersActorSetTree	*asPotentialImmigrantMothers[rc_size_AGE_RANGE][rc_size_SIM_YEAR_RANGE][rc_size_GEO];

class asPotentialSchoolOneStudentActorSetTree : public ActorSetTree	{
	public:
		asPotentialSchoolOneStudentActorSetTree() {m_nId = 18;};
		~asPotentialSchoolOneStudentActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asPotentialSchoolOneStudentActorSetTree	*asPotentialSchoolOneStudent;

class asPotentialSchoolOneStudentToProcessActorSetTree : public ActorSetTree	{
	public:
		asPotentialSchoolOneStudentToProcessActorSetTree() {m_nId = 19;};
		~asPotentialSchoolOneStudentToProcessActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asPotentialSchoolOneStudentToProcessActorSetTree	*asPotentialSchoolOneStudentToProcess;

class asPotentialSchoolTwoStudentActorSetTree : public ActorSetTree	{
	public:
		asPotentialSchoolTwoStudentActorSetTree() {m_nId = 20;};
		~asPotentialSchoolTwoStudentActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asPotentialSchoolTwoStudentActorSetTree	*asPotentialSchoolTwoStudent;

class asPotentialSchoolTwoStudentToProcessActorSetTree : public ActorSetTree	{
	public:
		asPotentialSchoolTwoStudentToProcessActorSetTree() {m_nId = 21;};
		~asPotentialSchoolTwoStudentToProcessActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asPotentialSchoolTwoStudentToProcessActorSetTree	*asPotentialSchoolTwoStudentToProcess;

class asResidentBabiesActorSetTree : public ActorSetTree	{
	public:
		asResidentBabiesActorSetTree() {m_nId = 22;};
		~asResidentBabiesActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asResidentBabiesActorSetTree	*asResidentBabies;

class asResidentsAge0SexActorSetTree : public ActorSetTree	{
	public:
		asResidentsAge0SexActorSetTree() {m_nId = 23;};
		~asResidentsAge0SexActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asResidentsAge0SexActorSetTree	*asResidentsAge0Sex[rc_size_SEX];

class asResidentsAge0SexGeoActorSetTree : public ActorSetTree	{
	public:
		asResidentsAge0SexGeoActorSetTree() {m_nId = 24;};
		~asResidentsAge0SexGeoActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asResidentsAge0SexGeoActorSetTree	*asResidentsAge0SexGeo[rc_size_SEX][rc_size_GEO];

class asSimBornAge0ActorSetTree : public ActorSetTree	{
	public:
		asSimBornAge0ActorSetTree() {m_nId = 25;};
		~asSimBornAge0ActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asSimBornAge0ActorSetTree	*asSimBornAge0[rc_size_SEX][rc_size_GEO][rc_size_EDUC1_GROUP];

class asSimulatedObservationHeadsActorSetTree : public ActorSetTree	{
	public:
		asSimulatedObservationHeadsActorSetTree() {m_nId = 26;};
		~asSimulatedObservationHeadsActorSetTree() {};
		Observation	* Item(int nIndex) {return (Observation *) ActorSetTree::Item(nIndex);};
		Observation	* GetRandom(double dRandom) {return (Observation *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asSimulatedObservationHeadsActorSetTree	*asSimulatedObservationHeads[rc_size_POP_POOL];

class asWantToMoveActorSetTree : public ActorSetTree	{
	public:
		asWantToMoveActorSetTree() {m_nId = 27;};
		~asWantToMoveActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asWantToMoveActorSetTree	*asWantToMove[rc_size_GEO][rc_size_GEO][rc_size_SEX][rc_size_AGE_MIG];

class asWomenNoChildrenActorSetTree : public ActorSetTree	{
	public:
		asWomenNoChildrenActorSetTree() {m_nId = 28;};
		~asWomenNoChildrenActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asWomenNoChildrenActorSetTree	*asWomenNoChildren[rc_size_UNION1_GROUP][rc_size_FEMALE_SPOUSE_AGE][2];

class asWomenWithChildrenActorSetTree : public ActorSetTree	{
	public:
		asWomenWithChildrenActorSetTree() {m_nId = 29;};
		~asWomenWithChildrenActorSetTree() {};
		Person	* Item(int nIndex) {return (Person *) ActorSetTree::Item(nIndex);};
		Person	* GetRandom(double dRandom) {return (Person *) ActorSetTree::GetRandom(dRandom);};
		int	ActorCompare(ActorClass *poActor1, ActorClass *poActor2);

	private:
		ActorSetNode *RemoveAll() {return ActorSetTree::RemoveAll();}
};
extern Thread asWomenWithChildrenActorSetTree	*asWomenWithChildren[rc_size_UNION1_GROUP][rc_size_CHILD_AGEGR][rc_size_MOTH_AGEGR][2];

void PreSimulation0();
void PreSimulation1();
void PreSimulation2();
void PreSimulation3();
void PreSimulation4();

} // namespace 
