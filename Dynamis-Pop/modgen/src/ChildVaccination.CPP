#include "model.h"
static int HaveAGoodDay;
namespace mm {
#line 1 "../code/ChildVaccination.mpp"
//LABEL (ChildVaccination, EN)
/* NOTE(ChildVaccination, EN)
This module implements child vaccination. A child is assumed to be immunized if it received the following 8 doses of vaccines during the first year of life:

- 1 dose of Bacillus of Calmette and Guérin (BCG) – a vaccine against tuberculosis
- 3 doses of DPT – a mixture of three vaccines to immunize against diphtheria, pertussis and tetanus
- 3 doses of polio;
- 1 dose of measles.

Immunization is decided at birth and depends on a set of individual and mother's characteristics. A key predictor for immunization is if the mother has received prenatal care. Therefore, immunization is decided in two steps based on logistic regression models. Based on a first model we decide if a mother has received prenatal care. Immunization is then decided by models estimated separately for children of mothers who have received prenatal care and those whose mothers have not received prenatal care. The model has two parameters - PreNatalCareOdds and ChildVaccinationOdds - which contain sets of odds of receiving prenatal care respectively receiving all required vaccines. The lists of parameters are declared in classifications -- PRENATCARE_PARA,  VACCINATION_PARA -- which are part of the _CountryContext module, as they contains country specific categories like region and ethnicity. Immunization is decided according these models for all children born as residents during the simulation. Children born abroad during the simulation sample the immunization status of resident 0 year old peers living in the region assigned to them as region of entry.
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor Sets
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//EN Resident Babies born in simulation
//actor_set Person asResidentBabies filter is_alive && integer_age == 0 && is_resident && creation_type == CT_BIRTH;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Types
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//classification GOT_PRENAT_CARE                                                      //EN Received prenatal care
//{
//    GPC_NO,                                                                         //EN No prenatal care was received
//    GPC_YES                                                                         //EN Some prenatal care was received
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters
//{
//    double PreNatalCareOdds[SIM_YEAR_RANGE][PRENATCARE_PARA];                       //EN Pre-natal Care (odds)
//    double ChildVaccinationOdds[GOT_PRENAT_CARE][SIM_YEAR_RANGE][VACCINATION_PARA]; //EN Child Vaccination (odds)
//};

//parameter_group PG_ChildVaccination                                                 //EN Child Vaccination (Immunization)
//{
//    PreNatalCareOdds, ChildVaccinationOdds
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Actor Declarations
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//actor Person
//{
//    
//    ETHNICITY_SHORT ethnicity_short = aggregate( ethnicity, ETHNICITY_SHORT);       //EN Ethnicity
//    logical is_immunized = { FALSE };                                               //EN Received all vaccines 
//    GOT_PRENAT_CARE got_prenat_care = { GPC_NO };                                   //EN Mother received prenatal care
//    void DecideImmunizationStatusResidents();                                       //EN Decide immunization residents
//    void DecideImmunizationStatusImmigrants();                                      //EN Decide immunization foreign born
//};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Person::DecideImmunizationStatusResidents()
{
    if ( in_projected_time && is_resident )
    {
        int nIndex = SPLIT( time - mother_age_at_birth , IMMU_YOB_PART); // mothers yob index
        int nYear = RANGE_POS(SIM_YEAR_RANGE, calendar_year );

        // Prenatal care
        double dOddsCare = gprParam->PreNatalCareOdds[CHKBND(nYear,151,"PreNatalCareOdds",0)][CHKBND(PP_CONSTANT,11,"PreNatalCareOdds",1)];                                      // 0
        if ( educ_mother == EOL_MEDIUM) dOddsCare = dOddsCare * gprParam->PreNatalCareOdds[CHKBND(nYear,151,"PreNatalCareOdds",0)][CHKBND(PP_EDUCMO_1,11,"PreNatalCareOdds",1)];  // 1
        if ( educ_mother == EOL_HIGH)   dOddsCare = dOddsCare * gprParam->PreNatalCareOdds[CHKBND(nYear,151,"PreNatalCareOdds",0)][CHKBND(PP_EDUCMO_2,11,"PreNatalCareOdds",1)];  // 2
        if ( region_nat > REGN_00)      dOddsCare = dOddsCare * gprParam->PreNatalCareOdds[CHKBND(nYear,151,"PreNatalCareOdds",0)][CHKBND(region_nat +2,11,"PreNatalCareOdds",1)]; // 3+
        if (nIndex > 0)                dOddsCare = dOddsCare
            * gprParam->PreNatalCareOdds[CHKBND(nYear,151,"PreNatalCareOdds",0)][CHKBND(nIndex + SIZE(REGION_NAT) + 1,11,"PreNatalCareOdds",1)];
        if ( mother_age_at_birth < 18)  dOddsCare = dOddsCare * gprParam->PreNatalCareOdds[CHKBND(nYear,151,"PreNatalCareOdds",0)][CHKBND(PP_YOUNGMO,11,"PreNatalCareOdds",1)];

        double dProbCare = dOddsCare / (dOddsCare + 1);
        if (RandUniform(79) < dProbCare)  poStateFunctions->Set_got_prenat_care(  GPC_YES );

        // Immunization
        double dOdds = gprParam->ChildVaccinationOdds[CHKBND(got_prenat_care ,2,"ChildVaccinationOdds",0)][CHKBND(nYear,151,"ChildVaccinationOdds",1)][CHKBND(IP_CONSTANT,18,"ChildVaccinationOdds",2)];                                  // 0
        if (sex ==MALE)                 dOdds = dOdds * gprParam->ChildVaccinationOdds[CHKBND(got_prenat_care ,2,"ChildVaccinationOdds",0)][CHKBND(nYear,151,"ChildVaccinationOdds",1)][CHKBND(IP_MALE,18,"ChildVaccinationOdds",2)];      // 1
        if ( educ_mother == EOL_MEDIUM) dOdds = dOdds * gprParam->ChildVaccinationOdds[CHKBND(got_prenat_care ,2,"ChildVaccinationOdds",0)][CHKBND(nYear,151,"ChildVaccinationOdds",1)][CHKBND(IP_EDUCMO_1,18,"ChildVaccinationOdds",2)];  // 2
        if ( educ_mother == EOL_HIGH)   dOdds = dOdds * gprParam->ChildVaccinationOdds[CHKBND(got_prenat_care ,2,"ChildVaccinationOdds",0)][CHKBND(nYear,151,"ChildVaccinationOdds",1)][CHKBND(IP_EDUCMO_2,18,"ChildVaccinationOdds",2)];  // 3
        if ( region_nat > REGN_00)      dOdds = dOdds * gprParam->ChildVaccinationOdds[CHKBND(got_prenat_care ,2,"ChildVaccinationOdds",0)][CHKBND(nYear,151,"ChildVaccinationOdds",1)][CHKBND(region_nat +3,18,"ChildVaccinationOdds",2)]; // 4+
        if ( ethnicity_short > ES_00)   dOdds = dOdds
            * gprParam->ChildVaccinationOdds[CHKBND(got_prenat_care ,2,"ChildVaccinationOdds",0)][CHKBND(nYear,151,"ChildVaccinationOdds",1)][CHKBND( ethnicity_short + SIZE(REGION_NAT) + 2,18,"ChildVaccinationOdds",2)];
        if (nIndex > 0)                dOdds = dOdds
            * gprParam->ChildVaccinationOdds[CHKBND(got_prenat_care ,2,"ChildVaccinationOdds",0)][CHKBND(nYear,151,"ChildVaccinationOdds",1)][CHKBND(nIndex + SIZE(REGION_NAT) + SIZE(ETHNICITY_SHORT) + 1,18,"ChildVaccinationOdds",2)];
        if ( mother_age_at_birth < 18)  dOdds = dOdds * gprParam->ChildVaccinationOdds[CHKBND(got_prenat_care ,2,"ChildVaccinationOdds",0)][CHKBND(nYear,151,"ChildVaccinationOdds",1)][CHKBND(IP_YOUNGMO,18,"ChildVaccinationOdds",2)];

        double dProb = dOdds / (dOdds + 1);
        if (RandUniform(80) < dProb)  poStateFunctions->Set_is_immunized(  TRUE );
    }
}

void Person::DecideImmunizationStatusImmigrants()
{
    if ( in_projected_time && !is_resident )
    {
        if (asResidentBabies->Count() > 0)
        {
            auto prDonor = asResidentBabies->GetRandom(RandUniform(77));
             poStateFunctions->Set_got_prenat_care(  prDonor->got_prenat_care  );
             poStateFunctions->Set_is_immunized(  prDonor->is_immunized  );
        }
        else if (RandUniform(78) < 0.5)  poStateFunctions->Set_is_immunized(  TRUE ); //this should not happen
    }
}


} // namespace 
