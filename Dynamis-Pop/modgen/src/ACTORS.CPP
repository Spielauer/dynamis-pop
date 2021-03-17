#include "model.h"

namespace mm {

#if defined(MG_CHKSUM)
	bool ::gbCheckSumActivated = true;
#else
	bool ::gbCheckSumActivated = false;
#endif

#pragma warning(disable : 4723)

ModelApp gApp( Simulation, InitParameters, DeleteParameters,
	ProcessModelGeneratedParameters, InitActors, DeleteActors,
	InitTables, PreSimulation, PostSimulation, UserTables,
	VerifyParameters, ValidateParameters, InitActorSets, DeleteActorSets );

void AFXAPI AfxAbort()
{
	TRACE0("AfxAbort called\n");
	DebugBreak();
#ifdef _WINDOWS
	AfxWinTerm();
#endif
	abort();
}

Thread EventQueue *gpoEventQueue;

Thread double ::gdThreadCurrentTime = 0;

Thread double ::gdThreadEvents = 0;

Thread double ::gdEventsForTabulation = 0;

Parameters	*gprParam;

BOOL VerifyParameters( CString *pszError)
{
	BOOL	bResult = TRUE;

	(*pszError).Empty();
	if ( !VerifyParameter( pszError, _T("AgeImmigrantsScratch") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("AgeImmiSearchMother") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("AgeLeavingHome") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("AgeOfImmigrantMother") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("AgeSpecificFertility") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("BackMigrationHazard") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("BirthTrends") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("ChildMortalityBaseRisk") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("ChildMortalityRelativeRisks") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("ChildMortalityTrend") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("ChildVaccinationOdds") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("Educ1FirstCohortRefinedModel") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("Educ1GradOdds") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("Educ1Infrastructure") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("Educ1Model") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("Educ1StartOdds") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("Educ2AllowedDelays") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("Educ2DelayedProgressionIntake") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("Educ2DelayedRepetitionIntake") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("Educ2DirectProgressionIntake") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("Educ2DirectRepetitionIntake") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("Educ2PeriodSuccess") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("EducOneDropoutGrade") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("EducOneEntryAge") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("EducTrans1") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("EducTrans2") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("EmigrationDestination") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("EmigrationRatesDistrict") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("EndSchoolOneYear") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("EthnicityImmigrantsScratch") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("EthnicTransmission") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("FertilityModel") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("FirstBirthRates") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("HCICoefficients") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("HigherOrderBirthsPara") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("ImmiPoolDestination") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("ImmiPoolSize") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("ImmiScratchDestination") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("InUnionProbNoChildren") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("InUnionProbWithChildren") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("LifeExpectancy") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("MicroDataInputFile") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("MigrationDestination") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("MigrationProbability") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("MigrationTryKeepingFamiliesTogether") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("ModelBackmigration") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("ModelEmigration") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("ModelImmigrationFromPools") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("ModelImmigrationFromScratch") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("ModelMigration") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("MortalityModel") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("MortalityTable") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("NumberImmigrantsFromScratch") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("PartnerAgeDistribution") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("PartnerCharacteristicDistribution") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("PreNatalCareOdds") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("PreSchoolAttendance") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("ProbStayWithMother") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("ProportionStunting") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("SchoolOneInterruptionRate") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("SchoolOneRepetitionRate") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("SchoolQuality") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("SexRatio") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("StartPopSampleSize") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("StartSchoolOneYear") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("TotalFertilityRate") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("Union1Choice") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("Union1ParametersCMN") ) ) bResult = FALSE;
	if ( !VerifyParameter( pszError, _T("Union1ParametersHazards") ) ) bResult = FALSE;
	return bResult;
}

void InitParameters()
{
	gprParam = new Parameters();
}

void DeleteParameters()
{
	delete gprParam;
}

void SetCaseWeight( double dCaseWeight, double dCaseSubsampleWeight )
{
	gpoEventQueue->SetCaseWeight( dCaseWeight, dCaseSubsampleWeight ); 
}

// partitions
real p_points_AGE_MIG_PART[] = { 5, 10, 15, 20, 25, 30, 35, 40, 45, 50,
	 55, 60 };
real p_points_BIRTH_AGE_PART[] = { 35, 40, 45 };
real p_points_CHILD_AGEGR_PART[] = { 1, 3, 6, 9, 12, 15 };
real p_points_DUR_TIME_IN_PARITY[] = { 1, 3, 6, 9, 12 };
real p_points_IMMU_YOB_PART[] = { 1971, 1976, 1981 };
real p_points_MOTH_AGEGR_PART[] = { 20, 25, 30, 35, 40 };
real p_points_tab_age15[] = { 15, 30, 60 };
real p_points_TAB_AGEC5[] = { 1, 5, 10, 15, 20, 25, 30, 35, 40, 45,
	 50, 55, 60, 65, 70, 75, 80, 85, 90, 95 };

Parameters::Parameters()
{
	Initializer( &AgeImmigrantsScratch, _T("AgeImmigrantsScratch") );
	Initializer( &AgeImmiSearchMother, _T("AgeImmiSearchMother") );
	Initializer( &AgeLeavingHome, _T("AgeLeavingHome") );
	Initializer( &AgeOfImmigrantMother, _T("AgeOfImmigrantMother") );
	Initializer( &AgeSpecificFertility, _T("AgeSpecificFertility") );
	memset( &AgeSpecificFertilityRate, '\0', 5738 * sizeof( double ) );
	Initializer( &BackMigrationHazard, _T("BackMigrationHazard") );
	Initializer( &BirthTrends, _T("BirthTrends") );
	Initializer( &ChildMortalityBaseRisk, _T("ChildMortalityBaseRisk") );
	Initializer( &ChildMortalityRelativeRisks, _T("ChildMortalityRelativeRisks") );
	Initializer( &ChildMortalityTrend, _T("ChildMortalityTrend") );
	Initializer( &ChildVaccinationOdds, _T("ChildVaccinationOdds") );
	Initializer( &Educ1FirstCohortRefinedModel, _T("Educ1FirstCohortRefinedModel") );
	Initializer( &Educ1GradOdds, _T("Educ1GradOdds") );
	Initializer( &Educ1Infrastructure, _T("Educ1Infrastructure") );
	Initializer( &Educ1Model, _T("Educ1Model") );
	Initializer( &Educ1StartOdds, _T("Educ1StartOdds") );
	Initializer( &Educ2AllowedDelays, _T("Educ2AllowedDelays") );
	Initializer( &Educ2DelayedProgressionIntake, _T("Educ2DelayedProgressionIntake") );
	Initializer( &Educ2DelayedRepetitionIntake, _T("Educ2DelayedRepetitionIntake") );
	Initializer( &Educ2DirectProgressionIntake, _T("Educ2DirectProgressionIntake") );
	Initializer( &Educ2DirectRepetitionIntake, _T("Educ2DirectRepetitionIntake") );
	Initializer( &Educ2PeriodSuccess, _T("Educ2PeriodSuccess") );
	Initializer( &EducOneDropoutGrade, _T("EducOneDropoutGrade") );
	Initializer( &EducOneEntryAge, _T("EducOneEntryAge") );
	Initializer( &EducTrans1, _T("EducTrans1") );
	Initializer( &EducTrans2, _T("EducTrans2") );
	Initializer( &EmigrationDestination, _T("EmigrationDestination") );
	Initializer( &EmigrationRatesDistrict, _T("EmigrationRatesDistrict") );
	Initializer( &EndSchoolOneYear, _T("EndSchoolOneYear") );
	memset( &EndSchoolTwoYear, '\0', 1 * sizeof( double ) );
	Initializer( &EthnicityImmigrantsScratch, _T("EthnicityImmigrantsScratch") );
	Initializer( &EthnicTransmission, _T("EthnicTransmission") );
	Initializer( &FertilityModel, _T("FertilityModel") );
	Initializer( &FirstBirthRates, _T("FirstBirthRates") );
	Initializer( &HCICoefficients, _T("HCICoefficients") );
	Initializer( &HigherOrderBirthsPara, _T("HigherOrderBirthsPara") );
	Initializer( &ImmiPoolDestination, _T("ImmiPoolDestination") );
	Initializer( &ImmiPoolSize, _T("ImmiPoolSize") );
	Initializer( &ImmiScratchDestination, _T("ImmiScratchDestination") );
	Initializer( &InUnionProbNoChildren, _T("InUnionProbNoChildren") );
	Initializer( &InUnionProbWithChildren, _T("InUnionProbWithChildren") );
	Initializer( &LifeExpectancy, _T("LifeExpectancy") );
	Initializer( &MicroDataInputFile, _T("MicroDataInputFile") );
	memset( &MicroDataInputFileSize, '\0', 1 * sizeof( long ) );
	Initializer( &MigrationDestination, _T("MigrationDestination") );
	Initializer( &MigrationProbability, _T("MigrationProbability") );
	Initializer( &MigrationTryKeepingFamiliesTogether, _T("MigrationTryKeepingFamiliesTogether") );
	Initializer( &ModelBackmigration, _T("ModelBackmigration") );
	Initializer( &ModelEmigration, _T("ModelEmigration") );
	Initializer( &ModelImmigrationFromPools, _T("ModelImmigrationFromPools") );
	Initializer( &ModelImmigrationFromScratch, _T("ModelImmigrationFromScratch") );
	Initializer( &ModelMigration, _T("ModelMigration") );
	Initializer( &MortalityModel, _T("MortalityModel") );
	Initializer( &MortalityTable, _T("MortalityTable") );
	memset( &MortalityTrend, '\0', 302 * sizeof( double ) );
	Initializer( &NumberImmigrantsFromScratch, _T("NumberImmigrantsFromScratch") );
	Initializer( &PartnerAgeDistribution, _T("PartnerAgeDistribution") );
	Initializer( &PartnerCharacteristicDistribution, _T("PartnerCharacteristicDistribution") );
	Initializer( &PreNatalCareOdds, _T("PreNatalCareOdds") );
	Initializer( &PreSchoolAttendance, _T("PreSchoolAttendance") );
	Initializer( &ProbStayWithMother, _T("ProbStayWithMother") );
	Initializer( &ProportionStunting, _T("ProportionStunting") );
	memset( &ScalingFactor, '\0', 1 * sizeof( double ) );
	Initializer( &SchoolOneInterruptionRate, _T("SchoolOneInterruptionRate") );
	Initializer( &SchoolOneRepetitionRate, _T("SchoolOneRepetitionRate") );
	Initializer( &SchoolQuality, _T("SchoolQuality") );
	Initializer( &SexRatio, _T("SexRatio") );
	Initializer( &StartPopSampleSize, _T("StartPopSampleSize") );
	Initializer( &StartSchoolOneYear, _T("StartSchoolOneYear") );
	memset( &StartSchoolTwoYear, '\0', 1 * sizeof( double ) );
	Initializer( &TotalFertilityRate, _T("TotalFertilityRate") );
	Initializer( &Union1Choice, _T("Union1Choice") );
	memset( &Union1FormationHazard, '\0', 38556 * sizeof( double ) );
	Initializer( &Union1ParametersCMN, _T("Union1ParametersCMN") );
	Initializer( &Union1ParametersHazards, _T("Union1ParametersHazards") );
}

void ProcessModelGeneratedParameters()
{
	ProcessParameter( &gprParam->AgeSpecificFertilityRate, _T("AgeSpecificFertilityRate"), NULL );
	ProcessParameter( &gprParam->EndSchoolTwoYear, _T("EndSchoolTwoYear"), NULL );
	ProcessParameter( &gprParam->MicroDataInputFileSize, _T("MicroDataInputFileSize"), NULL );
	ProcessParameter( &gprParam->MortalityTrend, _T("MortalityTrend"), NULL );
	ProcessParameter( &gprParam->ScalingFactor, _T("ScalingFactor"), NULL );
	ProcessParameter( &gprParam->StartSchoolTwoYear, _T("StartSchoolTwoYear"), NULL );
	ProcessParameter( &gprParam->Union1FormationHazard, _T("Union1FormationHazard"), NULL );
}

bool Lookup_AgeImmigrantsScratch( long lNumber, int nInd0, int *pnDest0 )
{
	return CumRateSearch( lNumber, (long *) &gprParam->AgeImmigrantsScratch[nInd0], 1, 101, pnDest0 );
}

bool Lookup_AgeImmigrantsScratch( double dNumber, int nInd0, int *pnDest0 )
{

	long lNumber;
	lNumber = (long) ( dNumber * RG_M_DOUBLE);

	return CumRateSearch( lNumber, (long *) &gprParam->AgeImmigrantsScratch[nInd0], 1, 101, pnDest0 );
}

int EmptyDistributions_AgeImmigrantsScratch( void )
{
	return ((ParameterSymbol *)gppiModelInfo->aoParameters[0])->m_nEmptyDistributions;
}

bool Lookup_AgeOfImmigrantMother( long lNumber, int *pnDest0 )
{
	return CumRateSearch( lNumber, (long *) &gprParam->AgeOfImmigrantMother, 1, 38, pnDest0 );
}

bool Lookup_AgeOfImmigrantMother( double dNumber, int *pnDest0 )
{

	long lNumber;
	lNumber = (long) ( dNumber * RG_M_DOUBLE);

	return CumRateSearch( lNumber, (long *) &gprParam->AgeOfImmigrantMother, 1, 38, pnDest0 );
}

int EmptyDistributions_AgeOfImmigrantMother( void )
{
	return ((ParameterSymbol *)gppiModelInfo->aoParameters[3])->m_nEmptyDistributions;
}

bool Lookup_EducOneDropoutGrade( long lNumber, int nInd0, int nInd1, int nInd2, int *pnDest0 )
{
	return CumRateSearch( lNumber, (long *) &gprParam->EducOneDropoutGrade[nInd0][nInd1][nInd2], 1, 6, pnDest0 );
}

bool Lookup_EducOneDropoutGrade( double dNumber, int nInd0, int nInd1, int nInd2, int *pnDest0 )
{

	long lNumber;
	lNumber = (long) ( dNumber * RG_M_DOUBLE);

	return CumRateSearch( lNumber, (long *) &gprParam->EducOneDropoutGrade[nInd0][nInd1][nInd2], 1, 6, pnDest0 );
}

int EmptyDistributions_EducOneDropoutGrade( void )
{
	return ((ParameterSymbol *)gppiModelInfo->aoParameters[23])->m_nEmptyDistributions;
}

bool Lookup_EducOneEntryAge( long lNumber, int nInd0, int nInd1, int nInd2, int *pnDest0 )
{
	return CumRateSearch( lNumber, (long *) &gprParam->EducOneEntryAge[nInd0][nInd1][nInd2], 1, 4, pnDest0 );
}

bool Lookup_EducOneEntryAge( double dNumber, int nInd0, int nInd1, int nInd2, int *pnDest0 )
{

	long lNumber;
	lNumber = (long) ( dNumber * RG_M_DOUBLE);

	return CumRateSearch( lNumber, (long *) &gprParam->EducOneEntryAge[nInd0][nInd1][nInd2], 1, 4, pnDest0 );
}

int EmptyDistributions_EducOneEntryAge( void )
{
	return ((ParameterSymbol *)gppiModelInfo->aoParameters[24])->m_nEmptyDistributions;
}

bool Lookup_EmigrationDestination( long lNumber, int *pnDest0 )
{
	return CumRateSearch( lNumber, (long *) &gprParam->EmigrationDestination, 1, 1, pnDest0 );
}

bool Lookup_EmigrationDestination( double dNumber, int *pnDest0 )
{

	long lNumber;
	lNumber = (long) ( dNumber * RG_M_DOUBLE);

	return CumRateSearch( lNumber, (long *) &gprParam->EmigrationDestination, 1, 1, pnDest0 );
}

int EmptyDistributions_EmigrationDestination( void )
{
	return ((ParameterSymbol *)gppiModelInfo->aoParameters[27])->m_nEmptyDistributions;
}

bool Lookup_EthnicityImmigrantsScratch( long lNumber, int nInd0, int nInd1, int *pnDest0 )
{
	return CumRateSearch( lNumber, (long *) &gprParam->EthnicityImmigrantsScratch[nInd0][nInd1], 1, 10, pnDest0 );
}

bool Lookup_EthnicityImmigrantsScratch( double dNumber, int nInd0, int nInd1, int *pnDest0 )
{

	long lNumber;
	lNumber = (long) ( dNumber * RG_M_DOUBLE);

	return CumRateSearch( lNumber, (long *) &gprParam->EthnicityImmigrantsScratch[nInd0][nInd1], 1, 10, pnDest0 );
}

int EmptyDistributions_EthnicityImmigrantsScratch( void )
{
	return ((ParameterSymbol *)gppiModelInfo->aoParameters[31])->m_nEmptyDistributions;
}

bool Lookup_EthnicTransmission( long lNumber, int nInd0, int nInd1, int *pnDest0 )
{
	return CumRateSearch( lNumber, (long *) &gprParam->EthnicTransmission[nInd0][nInd1], 1, 10, pnDest0 );
}

bool Lookup_EthnicTransmission( double dNumber, int nInd0, int nInd1, int *pnDest0 )
{

	long lNumber;
	lNumber = (long) ( dNumber * RG_M_DOUBLE);

	return CumRateSearch( lNumber, (long *) &gprParam->EthnicTransmission[nInd0][nInd1], 1, 10, pnDest0 );
}

int EmptyDistributions_EthnicTransmission( void )
{
	return ((ParameterSymbol *)gppiModelInfo->aoParameters[32])->m_nEmptyDistributions;
}

bool Lookup_ImmiPoolDestination( long lNumber, int nInd0, int nInd1, int *pnDest0 )
{
	return CumRateSearch( lNumber, (long *) &gprParam->ImmiPoolDestination[nInd0][nInd1], 1, 12, pnDest0 );
}

bool Lookup_ImmiPoolDestination( double dNumber, int nInd0, int nInd1, int *pnDest0 )
{

	long lNumber;
	lNumber = (long) ( dNumber * RG_M_DOUBLE);

	return CumRateSearch( lNumber, (long *) &gprParam->ImmiPoolDestination[nInd0][nInd1], 1, 12, pnDest0 );
}

int EmptyDistributions_ImmiPoolDestination( void )
{
	return ((ParameterSymbol *)gppiModelInfo->aoParameters[37])->m_nEmptyDistributions;
}

bool Lookup_ImmiScratchDestination( long lNumber, int nInd0, int nInd1, int *pnDest0 )
{
	return CumRateSearch( lNumber, (long *) &gprParam->ImmiScratchDestination[nInd0][nInd1], 1, 12, pnDest0 );
}

bool Lookup_ImmiScratchDestination( double dNumber, int nInd0, int nInd1, int *pnDest0 )
{

	long lNumber;
	lNumber = (long) ( dNumber * RG_M_DOUBLE);

	return CumRateSearch( lNumber, (long *) &gprParam->ImmiScratchDestination[nInd0][nInd1], 1, 12, pnDest0 );
}

int EmptyDistributions_ImmiScratchDestination( void )
{
	return ((ParameterSymbol *)gppiModelInfo->aoParameters[39])->m_nEmptyDistributions;
}

bool Lookup_MigrationDestination( long lNumber, int nInd0, int nInd1, int nInd2, int *pnDest0 )
{
	return CumRateSearch( lNumber, (long *) &gprParam->MigrationDestination[nInd0][nInd1][nInd2], 1, 12, pnDest0 );
}

bool Lookup_MigrationDestination( double dNumber, int nInd0, int nInd1, int nInd2, int *pnDest0 )
{

	long lNumber;
	lNumber = (long) ( dNumber * RG_M_DOUBLE);

	return CumRateSearch( lNumber, (long *) &gprParam->MigrationDestination[nInd0][nInd1][nInd2], 1, 12, pnDest0 );
}

int EmptyDistributions_MigrationDestination( void )
{
	return ((ParameterSymbol *)gppiModelInfo->aoParameters[45])->m_nEmptyDistributions;
}

bool Lookup_PartnerCharacteristicDistribution( long lNumber, int nInd0, int *pnDest0 )
{
	return CumRateSearch( lNumber, (long *) &gprParam->PartnerCharacteristicDistribution[nInd0], 1, 3, pnDest0 );
}

bool Lookup_PartnerCharacteristicDistribution( double dNumber, int nInd0, int *pnDest0 )
{

	long lNumber;
	lNumber = (long) ( dNumber * RG_M_DOUBLE);

	return CumRateSearch( lNumber, (long *) &gprParam->PartnerCharacteristicDistribution[nInd0], 1, 3, pnDest0 );
}

int EmptyDistributions_PartnerCharacteristicDistribution( void )
{
	return ((ParameterSymbol *)gppiModelInfo->aoParameters[58])->m_nEmptyDistributions;
}

ETHNICITY_SHORT ETHNICITY_To_ETHNICITY_SHORT( ETHNICITY eLevel )
{
	int nResult;

	switch (eLevel) {
		case 0: nResult = 0; break;
		case 8: nResult = 6; break;
		case 5: nResult = 3; break;
		case 2: nResult = 2; break;
		case 7: nResult = 5; break;
		case 4: nResult = 1; break;
		case 1: nResult = 4; break;
		case 9: nResult = 6; break;
		case 6: nResult = 4; break;
		case 3: nResult = 0; break;
	};
	return ETHNICITY_SHORT( nResult );
}

REGION_INT GEO_To_REGION_INT( GEO eLevel )
{
	int nResult;

	switch (eLevel) {
		case 0: nResult = 0; break;
		case 11: nResult = 4; break;
		case 8: nResult = 2; break;
		case 5: nResult = 1; break;
		case 2: nResult = 0; break;
		case 10: nResult = 3; break;
		case 7: nResult = 2; break;
		case 4: nResult = 1; break;
		case 1: nResult = 0; break;
		case 12: nResult = 5; break;
		case 9: nResult = 4; break;
		case 6: nResult = 2; break;
		case 3: nResult = 1; break;
	};
	return REGION_INT( nResult );
}

BIRTH1_LOC GEO_To_BIRTH1_LOC( GEO eLevel )
{
	int nResult;

	switch (eLevel) {
		case 0: nResult = 0; break;
		case 11: nResult = 4; break;
		case 8: nResult = 2; break;
		case 5: nResult = 1; break;
		case 2: nResult = 0; break;
		case 10: nResult = 3; break;
		case 7: nResult = 2; break;
		case 4: nResult = 1; break;
		case 1: nResult = 0; break;
		case 12: nResult = 4; break;
		case 9: nResult = 4; break;
		case 6: nResult = 2; break;
		case 3: nResult = 1; break;
	};
	return BIRTH1_LOC( nResult );
}

GEO_NAT GEO_To_GEO_NAT( GEO eLevel )
{
	int nResult;

	switch (eLevel) {
		case 0: nResult = 0; break;
		case 11: nResult = 11; break;
		case 8: nResult = 8; break;
		case 5: nResult = 5; break;
		case 2: nResult = 2; break;
		case 10: nResult = 10; break;
		case 7: nResult = 7; break;
		case 4: nResult = 4; break;
		case 1: nResult = 1; break;
		case 12: nResult = 11; break;
		case 9: nResult = 9; break;
		case 6: nResult = 6; break;
		case 3: nResult = 3; break;
	};
	return GEO_NAT( nResult );
}

REGION_NAT REGION_INT_To_REGION_NAT( REGION_INT eLevel )
{
	int nResult;

	switch (eLevel) {
		case 0: nResult = 0; break;
		case 5: nResult = 4; break;
		case 2: nResult = 2; break;
		case 4: nResult = 4; break;
		case 1: nResult = 1; break;
		case 3: nResult = 3; break;
	};
	return REGION_NAT( nResult );
}

bool IntIs_BIRTH1_GROUP( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
		case 1: return true;
		case 2: return true;
	};
	return false;
}

BIRTH1_GROUP IntTo_BIRTH1_GROUP( int nNumber )
{
	switch (nNumber) {
		case 0: return B1G_00;
		case 1: return B1G_01;
		case 2: return B1G_02;
	};
	return B1G_02;
}

int IntFrom_BIRTH1_GROUP( BIRTH1_GROUP nLevel )
{
	switch (nLevel) {
		case B1G_00: return 0;
		case B1G_01: return 1;
		case B1G_02: return 2;
	};
	return -1;
}

bool IntIs_BIRTH1_LOC( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
		case 1: return true;
		case 2: return true;
		case 3: return true;
		case 4: return true;
	};
	return false;
}

BIRTH1_LOC IntTo_BIRTH1_LOC( int nNumber )
{
	switch (nNumber) {
		case 0: return B1L_00;
		case 1: return B1L_01;
		case 2: return B1L_02;
		case 3: return B1L_03;
		case 4: return B1L_04;
	};
	return B1L_04;
}

int IntFrom_BIRTH1_LOC( BIRTH1_LOC nLevel )
{
	switch (nLevel) {
		case B1L_00: return 0;
		case B1L_01: return 1;
		case B1L_02: return 2;
		case B1L_03: return 3;
		case B1L_04: return 4;
	};
	return -1;
}

bool IntIs_CHILD_AGEGR( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
		case 1: return true;
		case 3: return true;
		case 6: return true;
		case 9: return true;
		case 12: return true;
		case 15: return true;
	};
	return false;
}

CHILD_AGEGR IntTo_CHILD_AGEGR( int nNumber )
{
	switch (nNumber) {
		case 0: return CA00;
		case 1: return CA01;
		case 3: return CA03;
		case 6: return CA06;
		case 9: return CA09;
		case 12: return CA12;
		case 15: return CA15;
	};
	return CA15;
}

int IntFrom_CHILD_AGEGR( CHILD_AGEGR nLevel )
{
	switch (nLevel) {
		case CA00: return 0;
		case CA01: return 1;
		case CA03: return 3;
		case CA06: return 6;
		case CA09: return 9;
		case CA12: return 12;
		case CA15: return 15;
	};
	return -1;
}

bool IntIs_CHILD_MORTALITY_GROUP( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
		case 1: return true;
		case 2: return true;
		case 3: return true;
		case 4: return true;
		case 5: return true;
		case 6: return true;
		case 7: return true;
		case 8: return true;
	};
	return false;
}

CHILD_MORTALITY_GROUP IntTo_CHILD_MORTALITY_GROUP( int nNumber )
{
	switch (nNumber) {
		case 0: return CMG_00;
		case 1: return CMG_01;
		case 2: return CMG_02;
		case 3: return CMG_03;
		case 4: return CMG_04;
		case 5: return CMG_05;
		case 6: return CMG_06;
		case 7: return CMG_07;
		case 8: return CMG_08;
	};
	return CMG_08;
}

int IntFrom_CHILD_MORTALITY_GROUP( CHILD_MORTALITY_GROUP nLevel )
{
	switch (nLevel) {
		case CMG_00: return 0;
		case CMG_01: return 1;
		case CMG_02: return 2;
		case CMG_03: return 3;
		case CMG_04: return 4;
		case CMG_05: return 5;
		case CMG_06: return 6;
		case CMG_07: return 7;
		case CMG_08: return 8;
	};
	return -1;
}

bool IntIs_EDUC1_GROUP( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
		case 1: return true;
		case 2: return true;
	};
	return false;
}

EDUC1_GROUP IntTo_EDUC1_GROUP( int nNumber )
{
	switch (nNumber) {
		case 0: return E1G_00;
		case 1: return E1G_01;
		case 2: return E1G_02;
	};
	return E1G_02;
}

int IntFrom_EDUC1_GROUP( EDUC1_GROUP nLevel )
{
	switch (nLevel) {
		case E1G_00: return 0;
		case E1G_01: return 1;
		case E1G_02: return 2;
	};
	return -1;
}

bool IntIs_ETHNICITY( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
		case 1: return true;
		case 2: return true;
		case 3: return true;
		case 4: return true;
		case 5: return true;
		case 6: return true;
		case 7: return true;
		case 8: return true;
		case 9: return true;
	};
	return false;
}

ETHNICITY IntTo_ETHNICITY( int nNumber )
{
	switch (nNumber) {
		case 0: return ETHNO_00;
		case 1: return ETHNO_01;
		case 2: return ETHNO_02;
		case 3: return ETHNO_03;
		case 4: return ETHNO_04;
		case 5: return ETHNO_05;
		case 6: return ETHNO_06;
		case 7: return ETHNO_07;
		case 8: return ETHNO_08;
		case 9: return ETHNO_09;
	};
	return ETHNO_09;
}

int IntFrom_ETHNICITY( ETHNICITY nLevel )
{
	switch (nLevel) {
		case ETHNO_00: return 0;
		case ETHNO_01: return 1;
		case ETHNO_02: return 2;
		case ETHNO_03: return 3;
		case ETHNO_04: return 4;
		case ETHNO_05: return 5;
		case ETHNO_06: return 6;
		case ETHNO_07: return 7;
		case ETHNO_08: return 8;
		case ETHNO_09: return 9;
	};
	return -1;
}

bool IntIs_ETHNICITY_SHORT( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
		case 1: return true;
		case 2: return true;
		case 3: return true;
		case 4: return true;
		case 5: return true;
		case 6: return true;
	};
	return false;
}

ETHNICITY_SHORT IntTo_ETHNICITY_SHORT( int nNumber )
{
	switch (nNumber) {
		case 0: return ES_00;
		case 1: return ES_01;
		case 2: return ES_02;
		case 3: return ES_03;
		case 4: return ES_04;
		case 5: return ES_05;
		case 6: return ES_06;
	};
	return ES_06;
}

int IntFrom_ETHNICITY_SHORT( ETHNICITY_SHORT nLevel )
{
	switch (nLevel) {
		case ES_00: return 0;
		case ES_01: return 1;
		case ES_02: return 2;
		case ES_03: return 3;
		case ES_04: return 4;
		case ES_05: return 5;
		case ES_06: return 6;
	};
	return -1;
}

bool IntIs_GEO( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
		case 1: return true;
		case 2: return true;
		case 3: return true;
		case 4: return true;
		case 5: return true;
		case 6: return true;
		case 7: return true;
		case 8: return true;
		case 9: return true;
		case 10: return true;
		case 11: return true;
	};
	return false;
}

GEO IntTo_GEO( int nNumber )
{
	switch (nNumber) {
		case 0: return GEO_00;
		case 1: return GEO_01;
		case 2: return GEO_02;
		case 3: return GEO_03;
		case 4: return GEO_04;
		case 5: return GEO_05;
		case 6: return GEO_06;
		case 7: return GEO_07;
		case 8: return GEO_08;
		case 9: return GEO_09;
		case 10: return GEO_10;
		case 11: return GEO_11;
	};
	return GEO_ABROAD;
}

int IntFrom_GEO( GEO nLevel )
{
	switch (nLevel) {
		case GEO_00: return 0;
		case GEO_01: return 1;
		case GEO_02: return 2;
		case GEO_03: return 3;
		case GEO_04: return 4;
		case GEO_05: return 5;
		case GEO_06: return 6;
		case GEO_07: return 7;
		case GEO_08: return 8;
		case GEO_09: return 9;
		case GEO_10: return 10;
		case GEO_11: return 11;
	};
	return -1;
}

bool IntIs_GEO_ABR( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
	};
	return false;
}

GEO_ABR IntTo_GEO_ABR( int nNumber )
{
	switch (nNumber) {
		case 0: return GEOA_00;
	};
	return GEOA_00;
}

int IntFrom_GEO_ABR( GEO_ABR nLevel )
{
	switch (nLevel) {
		case GEOA_00: return 0;
	};
	return -1;
}

bool IntIs_GEO_NAT( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
		case 1: return true;
		case 2: return true;
		case 3: return true;
		case 4: return true;
		case 5: return true;
		case 6: return true;
		case 7: return true;
		case 8: return true;
		case 9: return true;
		case 10: return true;
		case 11: return true;
	};
	return false;
}

GEO_NAT IntTo_GEO_NAT( int nNumber )
{
	switch (nNumber) {
		case 0: return GEON_00;
		case 1: return GEON_01;
		case 2: return GEON_02;
		case 3: return GEON_03;
		case 4: return GEON_04;
		case 5: return GEON_05;
		case 6: return GEON_06;
		case 7: return GEON_07;
		case 8: return GEON_08;
		case 9: return GEON_09;
		case 10: return GEON_10;
		case 11: return GEON_11;
	};
	return GEON_11;
}

int IntFrom_GEO_NAT( GEO_NAT nLevel )
{
	switch (nLevel) {
		case GEON_00: return 0;
		case GEON_01: return 1;
		case GEON_02: return 2;
		case GEON_03: return 3;
		case GEON_04: return 4;
		case GEON_05: return 5;
		case GEON_06: return 6;
		case GEON_07: return 7;
		case GEON_08: return 8;
		case GEON_09: return 9;
		case GEON_10: return 10;
		case GEON_11: return 11;
	};
	return -1;
}

bool IntIs_HIGHER_BIRTHS_PARA( int nNumber )
{
	switch (nNumber) {
		case 1: return true;
		case 2: return true;
		case 3: return true;
		case 4: return true;
		case 5: return true;
		case 6: return true;
		case 35: return true;
		case 40: return true;
		case 45: return true;
	};
	return false;
}

HIGHER_BIRTHS_PARA IntTo_HIGHER_BIRTHS_PARA( int nNumber )
{
	switch (nNumber) {
		case 1: return HBP_PERIOD1;
		case 2: return HBP_PERIOD2;
		case 3: return HBP_PERIOD3;
		case 4: return HBP_PERIOD4;
		case 5: return HBP_PERIOD5;
		case 6: return HBP_PERIOD6;
		case 35: return HBP_AGE35;
		case 40: return HBP_AGE40;
		case 45: return HBP_AGE45;
	};
	return HBP_EDUC2;
}

int IntFrom_HIGHER_BIRTHS_PARA( HIGHER_BIRTHS_PARA nLevel )
{
	switch (nLevel) {
		case HBP_PERIOD1: return 1;
		case HBP_PERIOD2: return 2;
		case HBP_PERIOD3: return 3;
		case HBP_PERIOD4: return 4;
		case HBP_PERIOD5: return 5;
		case HBP_PERIOD6: return 6;
		case HBP_AGE35: return 35;
		case HBP_AGE40: return 40;
		case HBP_AGE45: return 45;
		case HBP_EDUC1: return 1;
		case HBP_EDUC2: return 2;
	};
	return -1;
}

bool IntIs_MOTH_AGEGR( int nNumber )
{
	switch (nNumber) {
		case 20: return true;
		case 25: return true;
		case 30: return true;
		case 35: return true;
		case 40: return true;
	};
	return false;
}

MOTH_AGEGR IntTo_MOTH_AGEGR( int nNumber )
{
	switch (nNumber) {
		case 20: return CMA20;
		case 25: return CMA25;
		case 30: return CMA30;
		case 35: return CMA35;
		case 40: return CMA40;
	};
	return CMA40P;
}

int IntFrom_MOTH_AGEGR( MOTH_AGEGR nLevel )
{
	switch (nLevel) {
		case CMA20: return 20;
		case CMA25: return 25;
		case CMA30: return 30;
		case CMA35: return 35;
		case CMA40: return 40;
	};
	return -1;
}

bool IntIs_PRENATCARE_PARA( int nNumber )
{
	switch (nNumber) {
		case 1: return true;
		case 2: return true;
		case 3: return true;
		case 4: return true;
	};
	return false;
}

PRENATCARE_PARA IntTo_PRENATCARE_PARA( int nNumber )
{
	switch (nNumber) {
		case 1: return PP_EDUCMO_1;
		case 2: return PP_EDUCMO_2;
		case 3: return PP_REGION_3;
		case 4: return PP_REGION_4;
	};
	return PP_YOUNGMO;
}

int IntFrom_PRENATCARE_PARA( PRENATCARE_PARA nLevel )
{
	switch (nLevel) {
		case PP_EDUCMO_1: return 1;
		case PP_EDUCMO_2: return 2;
		case PP_REGION_1: return 1;
		case PP_REGION_2: return 2;
		case PP_REGION_3: return 3;
		case PP_REGION_4: return 4;
		case PP_COHORT_1: return 1;
		case PP_COHORT_2: return 2;
		case PP_COHORT_3: return 3;
	};
	return -1;
}

bool IntIs_REGION_INT( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
		case 1: return true;
		case 2: return true;
		case 3: return true;
		case 4: return true;
	};
	return false;
}

REGION_INT IntTo_REGION_INT( int nNumber )
{
	switch (nNumber) {
		case 0: return REGI_00;
		case 1: return REGI_01;
		case 2: return REGI_02;
		case 3: return REGI_03;
		case 4: return REGI_04;
	};
	return REGI_ABROAD;
}

int IntFrom_REGION_INT( REGION_INT nLevel )
{
	switch (nLevel) {
		case REGI_00: return 0;
		case REGI_01: return 1;
		case REGI_02: return 2;
		case REGI_03: return 3;
		case REGI_04: return 4;
	};
	return -1;
}

bool IntIs_REGION_NAT( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
		case 1: return true;
		case 2: return true;
		case 3: return true;
		case 4: return true;
	};
	return false;
}

REGION_NAT IntTo_REGION_NAT( int nNumber )
{
	switch (nNumber) {
		case 0: return REGN_00;
		case 1: return REGN_01;
		case 2: return REGN_02;
		case 3: return REGN_03;
		case 4: return REGN_04;
	};
	return REGN_04;
}

int IntFrom_REGION_NAT( REGION_NAT nLevel )
{
	switch (nLevel) {
		case REGN_00: return 0;
		case REGN_01: return 1;
		case REGN_02: return 2;
		case REGN_03: return 3;
		case REGN_04: return 4;
	};
	return -1;
}

bool IntIs_SELECTED_YEAR10( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
		case 10: return true;
		case 20: return true;
		case 30: return true;
		case 40: return true;
	};
	return false;
}

SELECTED_YEAR10 IntTo_SELECTED_YEAR10( int nNumber )
{
	switch (nNumber) {
		case 0: return SY00;
		case 10: return SY10;
		case 20: return SY20;
		case 30: return SY30;
		case 40: return SY40;
	};
	return SY40;
}

int IntFrom_SELECTED_YEAR10( SELECTED_YEAR10 nLevel )
{
	switch (nLevel) {
		case SY00: return 0;
		case SY10: return 10;
		case SY20: return 20;
		case SY30: return 30;
		case SY40: return 40;
	};
	return -1;
}

bool IntIs_SPOUSE_GROUP( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
		case 1: return true;
		case 2: return true;
	};
	return false;
}

SPOUSE_GROUP IntTo_SPOUSE_GROUP( int nNumber )
{
	switch (nNumber) {
		case 0: return SMG_00;
		case 1: return SMG_01;
		case 2: return SMG_02;
	};
	return SMG_02;
}

int IntFrom_SPOUSE_GROUP( SPOUSE_GROUP nLevel )
{
	switch (nLevel) {
		case SMG_00: return 0;
		case SMG_01: return 1;
		case SMG_02: return 2;
	};
	return -1;
}

bool IntIs_TAB_SEL_YEARS( int nNumber )
{
	switch (nNumber) {
		case 2000: return true;
		case 2010: return true;
		case 2020: return true;
		case 2030: return true;
		case 2040: return true;
		case 2050: return true;
		case 2060: return true;
	};
	return false;
}

TAB_SEL_YEARS IntTo_TAB_SEL_YEARS( int nNumber )
{
	switch (nNumber) {
		case 2000: return TSY_2000;
		case 2010: return TSY_2010;
		case 2020: return TSY_2020;
		case 2030: return TSY_2030;
		case 2040: return TSY_2040;
		case 2050: return TSY_2050;
		case 2060: return TSY_2060;
	};
	return TSY_2060;
}

int IntFrom_TAB_SEL_YEARS( TAB_SEL_YEARS nLevel )
{
	switch (nLevel) {
		case TSY_2000: return 2000;
		case TSY_2010: return 2010;
		case TSY_2020: return 2020;
		case TSY_2030: return 2030;
		case TSY_2040: return 2040;
		case TSY_2050: return 2050;
		case TSY_2060: return 2060;
	};
	return -1;
}

bool IntIs_UNION1_GROUP( int nNumber )
{
	switch (nNumber) {
		case 0: return true;
		case 1: return true;
		case 2: return true;
	};
	return false;
}

UNION1_GROUP IntTo_UNION1_GROUP( int nNumber )
{
	switch (nNumber) {
		case 0: return U1G_00;
		case 1: return U1G_01;
		case 2: return U1G_02;
	};
	return U1G_02;
}

int IntFrom_UNION1_GROUP( UNION1_GROUP nLevel )
{
	switch (nLevel) {
		case U1G_00: return 0;
		case U1G_01: return 1;
		case U1G_02: return 2;
	};
	return -1;
}

bool IntIs_VACCINATION_PARA( int nNumber )
{
	switch (nNumber) {
		case 1: return true;
		case 2: return true;
		case 3: return true;
		case 4: return true;
		case 5: return true;
		case 6: return true;
	};
	return false;
}

VACCINATION_PARA IntTo_VACCINATION_PARA( int nNumber )
{
	switch (nNumber) {
		case 1: return IP_EDUCMO_1;
		case 2: return IP_EDUCMO_2;
		case 3: return IP_REGION_3;
		case 4: return IP_REGION_4;
		case 5: return IP_ETNO_5;
		case 6: return IP_ETNO_6;
	};
	return IP_YOUNGMO;
}

int IntFrom_VACCINATION_PARA( VACCINATION_PARA nLevel )
{
	switch (nLevel) {
		case IP_EDUCMO_1: return 1;
		case IP_EDUCMO_2: return 2;
		case IP_REGION_1: return 1;
		case IP_REGION_2: return 2;
		case IP_REGION_3: return 3;
		case IP_REGION_4: return 4;
		case IP_ETNO_1: return 1;
		case IP_ETNO_2: return 2;
		case IP_ETNO_3: return 3;
		case IP_ETNO_4: return 4;
		case IP_ETNO_5: return 5;
		case IP_ETNO_6: return 6;
		case IP_COHORT_1: return 1;
		case IP_COHORT_2: return 2;
		case IP_COHORT_3: return 3;
	};
	return -1;
}

void InitActors()
{
	gpoEventQueue = new EventQueue();
	gpoEmptyCalibrator = NULL;
	gpoEmptyClock = NULL;
	gpoEmptyObservation = NULL;
	gpoEmptyPerson = NULL;
	gpoEmptyStartpopValues = NULL;
}

void DeleteActors()
{
	delete gpoEventQueue;
	gctrlSection.Lock();
	gctrlSection.Unlock();



	Calibrator *poCalibrator = gpoEmptyCalibrator;
	while ( poCalibrator != NULL ) {
		Calibrator *poCalibratorTemp = (Calibrator *) poCalibrator->poNext;
		::delete poCalibrator;
		poCalibrator = poCalibratorTemp;
	}


	Clock *poClock = gpoEmptyClock;
	while ( poClock != NULL ) {
		Clock *poClockTemp = (Clock *) poClock->poNext;
		::delete poClock;
		poClock = poClockTemp;
	}


	Observation *poObservation = gpoEmptyObservation;
	while ( poObservation != NULL ) {
		Observation *poObservationTemp = (Observation *) poObservation->poNext;
		::delete poObservation;
		poObservation = poObservationTemp;
	}


	Person *poPerson = gpoEmptyPerson;
	while ( poPerson != NULL ) {
		Person *poPersonTemp = (Person *) poPerson->poNext;
		::delete poPerson;
		poPerson = poPersonTemp;
	}


	StartpopValues *poStartpopValues = gpoEmptyStartpopValues;
	while ( poStartpopValues != NULL ) {
		StartpopValues *poStartpopValuesTemp = (StartpopValues *) poStartpopValues->poNext;
		::delete poStartpopValues;
		poStartpopValues = poStartpopValuesTemp;
	}
}

// sizeof actors
int ::ganSizeOfActors[5] = {
	sizeof(Calibrator) + sizeof(CalibratorDerivedStates),
	sizeof(Clock) + sizeof(ClockDerivedStates),
	sizeof(Observation) + sizeof(ObservationDerivedStates),
	sizeof(Person) + sizeof(PersonDerivedStates) + sizeof(PersonTable),
	sizeof(StartpopValues) + sizeof(StartpopValuesDerivedStates),
};

void *Calibrator::operator new( size_t count )
{
	Calibrator *poActor;
	if ( gpoEmptyCalibrator != NULL ) {
		poActor = gpoEmptyCalibrator;
		gpoEmptyCalibrator = (Calibrator *) gpoEmptyCalibrator->poNext;
	}
	else {
		poActor = ::new Calibrator( TRUE ) ;
	}
	poActor->__finished = false;
	return poActor;
}

void Calibrator::Set_lCalibratorToClock( Clock *value )
{
	Clock *poOldLink;

	if ( lCalibratorToClock != value ) {
		if (lCalibratorToClock != NULL)
		{
			lCalibratorToClock->UpdateTime();
		}
		if (value != NULL)
		{
			value->UpdateTime();
		}
		poOldLink = lCalibratorToClock;
		lCalibratorToClock = value;
		UpdateLinkedStates();
		if ( poOldLink != NULL && poOldLink->lClockToCalibrator == this) {
			poOldLink->Set_lClockToCalibrator( NULL );
		}
		if ( value != NULL && value->lClockToCalibrator != this) {
			value->Set_lClockToCalibrator( this );
		}
	}
}

void Calibrator::InitActor()
{
	poDerivedStates = new CalibratorDerivedStates(); 
	poDerivedStates->poParent = this;
	poStateFunctions = new CalibratorStateFunctions();
	poStateFunctions->poParent = this;
	poStateFunctions->poDerivedStates = poDerivedStates;
	poDerivedStates->poStateFunctions = poStateFunctions;
	bUpdates = FALSE;
	nActorNumber = 0;
	InitializeStates();
	lCalibratorToClock = NULL;
	mlCalibratorToPerson = new Calibrator_mlCalibratorToPerson( this );
}

void Calibrator::DeleteActor()
{
	delete poDerivedStates;
	delete poStateFunctions;
	delete mlCalibratorToPerson;
}

#pragma optimize( "", off )

void Calibrator::InitializeStates()
{
	// initialization of states
	actor_id = (long) 0;
	actor_subsample_weight = (double) 0;
	actor_weight = (double) 0;
	age = (TIME) 0;
	case_id = (long) 0;
	case_seed = (double) 0;
#line 80 "../code/MigrationGeneral.mpp"
	emigration_counter = (integer) 0;
#line 1483 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 79 "../code/MigrationGeneral.mpp"
	immigration_counter = (integer) 0;
#line 1486 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 99 "../code/MortalityDetailed.mpp"
	is_calibrated_child_mortality = FALSE != 0;
#line 1489 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 78 "../code/MigrationGeneral.mpp"
	migration_counter = (integer) 0;
#line 1492 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 93 "../code/MortalityDetailed.mpp"
	mort_female_0 = (double) 1.0;
#line 1495 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 94 "../code/MortalityDetailed.mpp"
	mort_female_1 = (double) 1.0;
#line 1498 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 95 "../code/MortalityDetailed.mpp"
	mort_female_2 = (double) 1.0;
#line 1501 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 96 "../code/MortalityDetailed.mpp"
	mort_female_3 = (double) 1.0;
#line 1504 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 97 "../code/MortalityDetailed.mpp"
	mort_female_4 = (double) 1.0;
#line 1507 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 87 "../code/MortalityDetailed.mpp"
	mort_male_0 = (double) 1.0;
#line 1510 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 88 "../code/MortalityDetailed.mpp"
	mort_male_1 = (double) 1.0;
#line 1513 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 89 "../code/MortalityDetailed.mpp"
	mort_male_2 = (double) 1.0;
#line 1516 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 90 "../code/MortalityDetailed.mpp"
	mort_male_3 = (double) 1.0;
#line 1519 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 91 "../code/MortalityDetailed.mpp"
	mort_male_4 = (double) 1.0;
#line 1522 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	time = (TIME) gdThreadCurrentTime;
	actorset_filter_asCalibrator = FALSE;
	calibrator_year = r_min_ALL_YEAR_RANGE;
	lCalibratorToClock_clock_year = (int) 0;
	__finished = FALSE;
	__time = (TIME) gdThreadCurrentTime;
	__events = 0;
}

#pragma optimize( "", on )

void Calibrator::UpdateLinkedStates()
{
	if ( lCalibratorToClock != NULL ) {
		poStateFunctions->Set_lCalibratorToClock_clock_year( lCalibratorToClock->clock_year );
	}
	else {
		poStateFunctions->Set_lCalibratorToClock_clock_year( int( 0 ) );
	}
}

void Calibrator::StartSpecial1()
{
	nSubSample = GetSubSample();
	InitializeStates();
	case_seed = GetCaseSeed();
	case_id = GetCaseID();
	lActorId = actor_id = GetObjectID();
	GetCaseWeight( &actor_weight, &actor_subsample_weight );
	lCalibratorToClock = NULL;
	mlCalibratorToPerson->RemoveAll();

	// initialization of actor set nodes
	poasCalibrator = NULL;

	// initialization of events
	poChildMortalityCalibration = gpoEventQueue->NewEvent( this, 2, 0 );
	//Initalize derived states
#line 34 "../code/ActorCalibrator.mpp"
	poStateFunctions->Set_calibrator_year( ALL_YEAR_RANGE( COERCE(ALL_YEAR_RANGE,lCalibratorToClock_clock_year) ) );
#line 1563 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	calibrator_year = calibrator_year < r_min_ALL_YEAR_RANGE ? r_min_ALL_YEAR_RANGE : calibrator_year > r_max_ALL_YEAR_RANGE ? r_max_ALL_YEAR_RANGE : calibrator_year;
	CountActor(0);
}

#pragma optimize( "", off )

void Calibrator::StartSpecial2()
{
	// initialization of expression and history states
	poStateFunctions->Set_time(TIME( CoarsenMantissa( time ) ), true);
	poStateFunctions->Set_age(TIME( CoarsenMantissa( age ) ), true);
	gdThreadCurrentTime = time;
	__time = time;
#line 15 "../code/ActorCalibrator.mpp"
	poStateFunctions->Set_actorset_filter_asCalibrator( logical( 1 ) );
#line 1579 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	UpdateLinkedStates();
	poPrev = NULL;
	poNext = gpoEventQueue->poHeadActor;
	if ( gpoEventQueue->poHeadActor != NULL ) {
		gpoEventQueue->poHeadActor->poPrev = this;
	}
	gpoEventQueue->poHeadActor = this;
	bUpdates = TRUE;
}

void Calibrator::FinishSpecial()
{
	CountFinishedActor(0);
	bUpdates = FALSE;
	gpoEventQueue->DestroyEvent( poChildMortalityCalibration );
	Set_lCalibratorToClock( NULL );
	mlCalibratorToPerson->Finish();
	if (poasCalibrator != NULL)
	{	asCalibrator->Remove(poasCalibrator);
		poasCalibrator = NULL;
	}

	// Remove actor from master list of active actors
	if ( poPrev != NULL )
	{
		poPrev->poNext = poNext;
	}
	else
	{
		if (gpoEventQueue->poHeadActor != NULL && gpoEventQueue->poHeadActor == this)
		{
			gpoEventQueue->poHeadActor = poNext;
		}
	}
	if ( poNext != NULL )
	{
		poNext->poPrev = poPrev;
	}
	poPrev = NULL;

	// Add actor to head of empty list of actors of this type
	poNext = gpoEmptyCalibrator;
	gpoEmptyCalibrator = this;
}

#pragma optimize( "", on )

void Calibrator::WaitUntil( double target_time )
{
	//don't increase gdEventsForTabulation if special events because they are all treated as one event for tabulation purposes
	if (!gpoEventQueue->m_bSpecialEvents)
	{
		gdEventsForTabulation++;
	}
	__events = gdEventsForTabulation;
	if ( time < (TIME) target_time )
	{
		__time = (TIME) target_time;
		WaitSpecial( (TIME) target_time - time );
		poStateFunctions->Set_time( (TIME) target_time );
	}
}

Calibrator *Calibrator::UpdateTime(  )
{
	__events = gdEventsForTabulation;
	if ( time < (TIME) gdThreadCurrentTime ) {
		__time = (TIME) gdThreadCurrentTime;
		WaitSpecial( (TIME) gdThreadCurrentTime - time );
		poStateFunctions->Set_time( (TIME) gdThreadCurrentTime );
	}
	return this;
}

void Calibrator::WaitUntilThisActor( double target_time )
{
	WaitUntil(target_time);
}

void Calibrator::WaitSpecial( TIME wait_time )
{

	// increment the age
		poStateFunctions->Set_age( age + wait_time );
}

double Calibrator::EventTime( int nEventNum, int *pnEventInfo )
{
	double	event_time = double( TIME_INFINITE );

	gdRuntimeTime = (double) time;
	gnRuntimeActorId = actor_id;
	switch( nEventNum ) {
		case 2:
			gszRuntimeEventName = _T("Calibrator.ChildMortalityCalibration(time)");
			event_time = CoarsenMantissa( (TIME) timeChildMortalityCalibration() );
			TRCEVNT(case_seed, _T("Calibrator"), actor_id, _T("timeChildMortalityCalibration"), event_time);
			break;
		case -1:
		default:;
	}
	return (TIME) event_time;
}

void Calibrator::Implement( int nEventNum, int nEventInfo )
{
	gdThreadEvents += 1;
	gdRuntimeTime = (double) time;
	gnRuntimeActorId = actor_id;
	switch( nEventNum )
	{
		case 2:
			CHKSUM(time, 0);
			TRCEVNT(case_seed, _T("Calibrator"), actor_id, _T("Calibrator.ChildMortalityCalibration"), time);
			gszRuntimeEventName = _T("Calibrator.ChildMortalityCalibration(implement)");
			ChildMortalityCalibration();
			break;
		case -1:
		default:;
	}
}

double Calibrator::Set_actor_subsample_weight( double value )
{
	return poStateFunctions->Set_actor_subsample_weight( value);
}

double Calibrator::Set_actor_weight( double value )
{
	return poStateFunctions->Set_actor_weight( value);
}

logical CalibratorStateFunctions::Set_is_calibrated_child_mortality( logical value )
{

	value = value != 0;
	if ( poParent->is_calibrated_child_mortality != value )
	{


		logical	OldValue;
		OldValue = poParent->is_calibrated_child_mortality;
		poParent->is_calibrated_child_mortality = value;

		if ( poParent->poChildMortalityCalibration->cRecalc == 0 ) {
			poParent->poChildMortalityCalibration->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poChildMortalityCalibration );
		};
	}
	return poParent->is_calibrated_child_mortality = value;
}

double CalibratorStateFunctions::Set_mort_female_0( double value )
{

	if ( poParent->mort_female_0 != value )
	{


		double	OldValue;
		OldValue = poParent->mort_female_0;
		poParent->mort_female_0 = value;

		{
			int	nIndex;
			int	nSize;
			Person	*prPerson;

			nSize = poParent->mlCalibratorToPerson->arData.GetSize();
			for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
				prPerson = poParent->mlCalibratorToPerson->arData[nIndex];
				if ( prPerson != NULL ) {
					prPerson->UpdateTime();
					prPerson->poStateFunctions->Set_lCalibrator_mort_female_0( poParent->mort_female_0 );
				}
			}
		}
	}
	return poParent->mort_female_0 = value;
}

double CalibratorStateFunctions::Set_mort_female_1( double value )
{

	if ( poParent->mort_female_1 != value )
	{


		double	OldValue;
		OldValue = poParent->mort_female_1;
		poParent->mort_female_1 = value;

		{
			int	nIndex;
			int	nSize;
			Person	*prPerson;

			nSize = poParent->mlCalibratorToPerson->arData.GetSize();
			for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
				prPerson = poParent->mlCalibratorToPerson->arData[nIndex];
				if ( prPerson != NULL ) {
					prPerson->UpdateTime();
					prPerson->poStateFunctions->Set_lCalibrator_mort_female_1( poParent->mort_female_1 );
				}
			}
		}
	}
	return poParent->mort_female_1 = value;
}

double CalibratorStateFunctions::Set_mort_female_2( double value )
{

	if ( poParent->mort_female_2 != value )
	{


		double	OldValue;
		OldValue = poParent->mort_female_2;
		poParent->mort_female_2 = value;

		{
			int	nIndex;
			int	nSize;
			Person	*prPerson;

			nSize = poParent->mlCalibratorToPerson->arData.GetSize();
			for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
				prPerson = poParent->mlCalibratorToPerson->arData[nIndex];
				if ( prPerson != NULL ) {
					prPerson->UpdateTime();
					prPerson->poStateFunctions->Set_lCalibrator_mort_female_2( poParent->mort_female_2 );
				}
			}
		}
	}
	return poParent->mort_female_2 = value;
}

double CalibratorStateFunctions::Set_mort_female_3( double value )
{

	if ( poParent->mort_female_3 != value )
	{


		double	OldValue;
		OldValue = poParent->mort_female_3;
		poParent->mort_female_3 = value;

		{
			int	nIndex;
			int	nSize;
			Person	*prPerson;

			nSize = poParent->mlCalibratorToPerson->arData.GetSize();
			for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
				prPerson = poParent->mlCalibratorToPerson->arData[nIndex];
				if ( prPerson != NULL ) {
					prPerson->UpdateTime();
					prPerson->poStateFunctions->Set_lCalibrator_mort_female_3( poParent->mort_female_3 );
				}
			}
		}
	}
	return poParent->mort_female_3 = value;
}

double CalibratorStateFunctions::Set_mort_female_4( double value )
{

	if ( poParent->mort_female_4 != value )
	{


		double	OldValue;
		OldValue = poParent->mort_female_4;
		poParent->mort_female_4 = value;

		{
			int	nIndex;
			int	nSize;
			Person	*prPerson;

			nSize = poParent->mlCalibratorToPerson->arData.GetSize();
			for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
				prPerson = poParent->mlCalibratorToPerson->arData[nIndex];
				if ( prPerson != NULL ) {
					prPerson->UpdateTime();
					prPerson->poStateFunctions->Set_lCalibrator_mort_female_4( poParent->mort_female_4 );
				}
			}
		}
	}
	return poParent->mort_female_4 = value;
}

double CalibratorStateFunctions::Set_mort_male_0( double value )
{

	if ( poParent->mort_male_0 != value )
	{


		double	OldValue;
		OldValue = poParent->mort_male_0;
		poParent->mort_male_0 = value;

		{
			int	nIndex;
			int	nSize;
			Person	*prPerson;

			nSize = poParent->mlCalibratorToPerson->arData.GetSize();
			for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
				prPerson = poParent->mlCalibratorToPerson->arData[nIndex];
				if ( prPerson != NULL ) {
					prPerson->UpdateTime();
					prPerson->poStateFunctions->Set_lCalibrator_mort_male_0( poParent->mort_male_0 );
				}
			}
		}
	}
	return poParent->mort_male_0 = value;
}

double CalibratorStateFunctions::Set_mort_male_1( double value )
{

	if ( poParent->mort_male_1 != value )
	{


		double	OldValue;
		OldValue = poParent->mort_male_1;
		poParent->mort_male_1 = value;

		{
			int	nIndex;
			int	nSize;
			Person	*prPerson;

			nSize = poParent->mlCalibratorToPerson->arData.GetSize();
			for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
				prPerson = poParent->mlCalibratorToPerson->arData[nIndex];
				if ( prPerson != NULL ) {
					prPerson->UpdateTime();
					prPerson->poStateFunctions->Set_lCalibrator_mort_male_1( poParent->mort_male_1 );
				}
			}
		}
	}
	return poParent->mort_male_1 = value;
}

double CalibratorStateFunctions::Set_mort_male_2( double value )
{

	if ( poParent->mort_male_2 != value )
	{


		double	OldValue;
		OldValue = poParent->mort_male_2;
		poParent->mort_male_2 = value;

		{
			int	nIndex;
			int	nSize;
			Person	*prPerson;

			nSize = poParent->mlCalibratorToPerson->arData.GetSize();
			for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
				prPerson = poParent->mlCalibratorToPerson->arData[nIndex];
				if ( prPerson != NULL ) {
					prPerson->UpdateTime();
					prPerson->poStateFunctions->Set_lCalibrator_mort_male_2( poParent->mort_male_2 );
				}
			}
		}
	}
	return poParent->mort_male_2 = value;
}

double CalibratorStateFunctions::Set_mort_male_3( double value )
{

	if ( poParent->mort_male_3 != value )
	{


		double	OldValue;
		OldValue = poParent->mort_male_3;
		poParent->mort_male_3 = value;

		{
			int	nIndex;
			int	nSize;
			Person	*prPerson;

			nSize = poParent->mlCalibratorToPerson->arData.GetSize();
			for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
				prPerson = poParent->mlCalibratorToPerson->arData[nIndex];
				if ( prPerson != NULL ) {
					prPerson->UpdateTime();
					prPerson->poStateFunctions->Set_lCalibrator_mort_male_3( poParent->mort_male_3 );
				}
			}
		}
	}
	return poParent->mort_male_3 = value;
}

double CalibratorStateFunctions::Set_mort_male_4( double value )
{

	if ( poParent->mort_male_4 != value )
	{


		double	OldValue;
		OldValue = poParent->mort_male_4;
		poParent->mort_male_4 = value;

		{
			int	nIndex;
			int	nSize;
			Person	*prPerson;

			nSize = poParent->mlCalibratorToPerson->arData.GetSize();
			for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
				prPerson = poParent->mlCalibratorToPerson->arData[nIndex];
				if ( prPerson != NULL ) {
					prPerson->UpdateTime();
					prPerson->poStateFunctions->Set_lCalibrator_mort_male_4( poParent->mort_male_4 );
				}
			}
		}
	}
	return poParent->mort_male_4 = value;
}

logical CalibratorStateFunctions::Set_actorset_filter_asCalibrator( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asCalibrator != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asCalibrator;
		poParent->actorset_filter_asCalibrator = value;

		// actorset filter
		if (poParent->actorset_filter_asCalibrator)
		{
			poParent->poasCalibrator = asCalibrator->NewNode(poParent);
		}
		else
		{
			asCalibrator->Remove(poParent->poasCalibrator);
			poParent->poasCalibrator = NULL;
		}
	}
	return poParent->actorset_filter_asCalibrator = value;
}

ALL_YEAR_RANGE CalibratorStateFunctions::Set_calibrator_year( ALL_YEAR_RANGE value )
{

	CHKLMT(value, calibrator_year, r_min_ALL_YEAR_RANGE, r_max_ALL_YEAR_RANGE);
	if ( value < r_min_ALL_YEAR_RANGE ) {
		value = r_min_ALL_YEAR_RANGE;
	}
	else if ( value > r_max_ALL_YEAR_RANGE ) {
		value = r_max_ALL_YEAR_RANGE;
	}
	if ( poParent->calibrator_year != value )
	{


		ALL_YEAR_RANGE	OldValue;
		OldValue = poParent->calibrator_year;
		poParent->calibrator_year = value;

	}
	return poParent->calibrator_year = value;
}

int CalibratorStateFunctions::Set_lCalibratorToClock_clock_year( int value )
{

	if ( poParent->lCalibratorToClock_clock_year != value )
	{


		int	OldValue;
		OldValue = poParent->lCalibratorToClock_clock_year;
		poParent->lCalibratorToClock_clock_year = value;

#line 34 "../code/ActorCalibrator.mpp"
		Set_calibrator_year( COERCE(ALL_YEAR_RANGE,poParent->lCalibratorToClock_clock_year) );
#line 2083 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lCalibratorToClock_clock_year = value;
}

void *Clock::operator new( size_t count )
{
	Clock *poActor;
	if ( gpoEmptyClock != NULL ) {
		poActor = gpoEmptyClock;
		gpoEmptyClock = (Clock *) gpoEmptyClock->poNext;
	}
	else {
		poActor = ::new Clock( TRUE ) ;
	}
	poActor->__finished = false;
	return poActor;
}

void Clock::Set_lClockToCalibrator( Calibrator *value )
{
	Calibrator *poOldLink;

	if ( lClockToCalibrator != value ) {
		if (lClockToCalibrator != NULL)
		{
			lClockToCalibrator->UpdateTime();
		}
		if (value != NULL)
		{
			value->UpdateTime();
		}
		poOldLink = lClockToCalibrator;
		lClockToCalibrator = value;
		UpdateLinkedStates();
		if ( poOldLink != NULL && poOldLink->lCalibratorToClock == this) {
			poOldLink->Set_lCalibratorToClock( NULL );
		}
		if ( value != NULL && value->lCalibratorToClock != this) {
			value->Set_lCalibratorToClock( this );
		}
	}
}

void Clock::InitActor()
{
	poDerivedStates = new ClockDerivedStates(); 
	poDerivedStates->poParent = this;
	poStateFunctions = new ClockStateFunctions();
	poStateFunctions->poParent = this;
	poStateFunctions->poDerivedStates = poDerivedStates;
	poDerivedStates->poStateFunctions = poStateFunctions;
	bUpdates = FALSE;
	nActorNumber = 1;
	InitializeStates();
	lClockToCalibrator = NULL;
}

void Clock::DeleteActor()
{
	delete poDerivedStates;
	delete poStateFunctions;
}

#pragma optimize( "", off )

void Clock::InitializeStates()
{
	// initialization of states
	actor_id = (long) 0;
	actor_subsample_weight = (double) 0;
	actor_weight = (double) 0;
	age = (TIME) 0;
	case_id = (long) 0;
	case_seed = (double) 0;
#line 28 "../code/ActorClock.mpp"
	clock_year = (int) 1900;
#line 2160 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 30 "../code/ActorClock.mpp"
	next_clock_year_end = (TIME) 1900;
#line 2163 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 33 "../code/ActorClock.mpp"
	next_clock_year_start = (TIME) 1900;
#line 2166 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 36 "../code/ActorClock.mpp"
	next_midyear_clock_event = (TIME) TIME_INFINITE;
#line 2169 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	time = (TIME) gdThreadCurrentTime;
#line 118 "../code/EducationPrimaryTracking.mpp"
	time_end_school_one_year = (TIME) TIME_INFINITE;
#line 2173 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 125 "../code/EducationSecondaryBase.mpp"
	time_end_school_two_year = (TIME) TIME_INFINITE;
#line 2176 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 117 "../code/EducationPrimaryTracking.mpp"
	time_start_school_one_year = (TIME) TIME_INFINITE;
#line 2179 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 124 "../code/EducationSecondaryBase.mpp"
	time_start_school_two_year = (TIME) TIME_INFINITE;
#line 2182 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	actorset_filter_asClock = FALSE;
	__finished = FALSE;
	__time = (TIME) gdThreadCurrentTime;
	__events = 0;
}

#pragma optimize( "", on )

void Clock::UpdateLinkedStates()
{
}

void Clock::StartSpecial1()
{
	nSubSample = GetSubSample();
	InitializeStates();
	case_seed = GetCaseSeed();
	case_id = GetCaseID();
	lActorId = actor_id = GetObjectID();
	GetCaseWeight( &actor_weight, &actor_subsample_weight );
	lClockToCalibrator = NULL;

	// initialization of actor set nodes
	poasClock = NULL;

	// initialization of events
	poClockMidyearEvent = gpoEventQueue->NewEvent( this, 3, 0 );
	poClockYearEndEvent = gpoEventQueue->NewEvent( this, 4, 0 );
	poClockYearStartEvent = gpoEventQueue->NewEvent( this, 5, 0 );
	poEndSchoolOneYearEvent = gpoEventQueue->NewEvent( this, 8, 0 );
	poEndSchoolTwoYearEvent = gpoEventQueue->NewEvent( this, 9, 0 );
	poStartSchoolOneYearEvent = gpoEventQueue->NewEvent( this, 29, 0 );
	poStartSchoolTwoYearEvent = gpoEventQueue->NewEvent( this, 30, 0 );
	//Initalize derived states
	CountActor(1);
}

#pragma optimize( "", off )

void Clock::StartSpecial2()
{
	// initialization of expression and history states
	poStateFunctions->Set_time(TIME( CoarsenMantissa( time ) ), true);
	poStateFunctions->Set_age(TIME( CoarsenMantissa( age ) ), true);
	gdThreadCurrentTime = time;
	__time = time;
#line 17 "../code/ActorClock.mpp"
	poStateFunctions->Set_actorset_filter_asClock( logical( 1 ) );
#line 2231 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	UpdateLinkedStates();
	poPrev = NULL;
	poNext = gpoEventQueue->poHeadActor;
	if ( gpoEventQueue->poHeadActor != NULL ) {
		gpoEventQueue->poHeadActor->poPrev = this;
	}
	gpoEventQueue->poHeadActor = this;
	bUpdates = TRUE;
}

void Clock::FinishSpecial()
{
	CountFinishedActor(1);
	bUpdates = FALSE;
	gpoEventQueue->DestroyEvent( poClockMidyearEvent );
	gpoEventQueue->DestroyEvent( poClockYearEndEvent );
	gpoEventQueue->DestroyEvent( poClockYearStartEvent );
	gpoEventQueue->DestroyEvent( poEndSchoolOneYearEvent );
	gpoEventQueue->DestroyEvent( poEndSchoolTwoYearEvent );
	gpoEventQueue->DestroyEvent( poStartSchoolOneYearEvent );
	gpoEventQueue->DestroyEvent( poStartSchoolTwoYearEvent );
	Set_lClockToCalibrator( NULL );
	if (poasClock != NULL)
	{	asClock->Remove(poasClock);
		poasClock = NULL;
	}

	// Remove actor from master list of active actors
	if ( poPrev != NULL )
	{
		poPrev->poNext = poNext;
	}
	else
	{
		if (gpoEventQueue->poHeadActor != NULL && gpoEventQueue->poHeadActor == this)
		{
			gpoEventQueue->poHeadActor = poNext;
		}
	}
	if ( poNext != NULL )
	{
		poNext->poPrev = poPrev;
	}
	poPrev = NULL;

	// Add actor to head of empty list of actors of this type
	poNext = gpoEmptyClock;
	gpoEmptyClock = this;
}

#pragma optimize( "", on )

void Clock::WaitUntil( double target_time )
{
	//don't increase gdEventsForTabulation if special events because they are all treated as one event for tabulation purposes
	if (!gpoEventQueue->m_bSpecialEvents)
	{
		gdEventsForTabulation++;
	}
	__events = gdEventsForTabulation;
	if ( time < (TIME) target_time )
	{
		__time = (TIME) target_time;
		WaitSpecial( (TIME) target_time - time );
		poStateFunctions->Set_time( (TIME) target_time );
	}
}

Clock *Clock::UpdateTime(  )
{
	__events = gdEventsForTabulation;
	if ( time < (TIME) gdThreadCurrentTime ) {
		__time = (TIME) gdThreadCurrentTime;
		WaitSpecial( (TIME) gdThreadCurrentTime - time );
		poStateFunctions->Set_time( (TIME) gdThreadCurrentTime );
	}
	return this;
}

void Clock::WaitUntilThisActor( double target_time )
{
	WaitUntil(target_time);
}

void Clock::WaitSpecial( TIME wait_time )
{

	// increment the age
		poStateFunctions->Set_age( age + wait_time );
}

double Clock::EventTime( int nEventNum, int *pnEventInfo )
{
	double	event_time = double( TIME_INFINITE );

	gdRuntimeTime = (double) time;
	gnRuntimeActorId = actor_id;
	switch( nEventNum ) {
		case 3:
			gszRuntimeEventName = _T("Clock.ClockMidyearEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeClockMidyearEvent() );
			TRCEVNT(case_seed, _T("Clock"), actor_id, _T("timeClockMidyearEvent"), event_time);
			break;
		case 4:
			gszRuntimeEventName = _T("Clock.ClockYearEndEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeClockYearEndEvent() );
			TRCEVNT(case_seed, _T("Clock"), actor_id, _T("timeClockYearEndEvent"), event_time);
			break;
		case 5:
			gszRuntimeEventName = _T("Clock.ClockYearStartEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeClockYearStartEvent() );
			TRCEVNT(case_seed, _T("Clock"), actor_id, _T("timeClockYearStartEvent"), event_time);
			break;
		case 8:
			gszRuntimeEventName = _T("Clock.EndSchoolOneYearEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeEndSchoolOneYearEvent() );
			TRCEVNT(case_seed, _T("Clock"), actor_id, _T("timeEndSchoolOneYearEvent"), event_time);
			break;
		case 9:
			gszRuntimeEventName = _T("Clock.EndSchoolTwoYearEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeEndSchoolTwoYearEvent() );
			TRCEVNT(case_seed, _T("Clock"), actor_id, _T("timeEndSchoolTwoYearEvent"), event_time);
			break;
		case 29:
			gszRuntimeEventName = _T("Clock.StartSchoolOneYearEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeStartSchoolOneYearEvent() );
			TRCEVNT(case_seed, _T("Clock"), actor_id, _T("timeStartSchoolOneYearEvent"), event_time);
			break;
		case 30:
			gszRuntimeEventName = _T("Clock.StartSchoolTwoYearEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeStartSchoolTwoYearEvent() );
			TRCEVNT(case_seed, _T("Clock"), actor_id, _T("timeStartSchoolTwoYearEvent"), event_time);
			break;
		case -1:
		default:;
	}
	return (TIME) event_time;
}

void Clock::Implement( int nEventNum, int nEventInfo )
{
	gdThreadEvents += 1;
	gdRuntimeTime = (double) time;
	gnRuntimeActorId = actor_id;
	switch( nEventNum )
	{
		case 3:
			CHKSUM(time, 0);
			TRCEVNT(case_seed, _T("Clock"), actor_id, _T("Clock.ClockMidyearEvent"), time);
			gszRuntimeEventName = _T("Clock.ClockMidyearEvent(implement)");
			ClockMidyearEvent();
			break;
		case 4:
			CHKSUM(time, 1);
			TRCEVNT(case_seed, _T("Clock"), actor_id, _T("Clock.ClockYearEndEvent"), time);
			gszRuntimeEventName = _T("Clock.ClockYearEndEvent(implement)");
			ClockYearEndEvent();
			break;
		case 5:
			CHKSUM(time, 2);
			TRCEVNT(case_seed, _T("Clock"), actor_id, _T("Clock.ClockYearStartEvent"), time);
			gszRuntimeEventName = _T("Clock.ClockYearStartEvent(implement)");
			ClockYearStartEvent();
			break;
		case 8:
			CHKSUM(time, 3);
			TRCEVNT(case_seed, _T("Clock"), actor_id, _T("Clock.EndSchoolOneYearEvent"), time);
			gszRuntimeEventName = _T("Clock.EndSchoolOneYearEvent(implement)");
			EndSchoolOneYearEvent();
			break;
		case 9:
			CHKSUM(time, 4);
			TRCEVNT(case_seed, _T("Clock"), actor_id, _T("Clock.EndSchoolTwoYearEvent"), time);
			gszRuntimeEventName = _T("Clock.EndSchoolTwoYearEvent(implement)");
			EndSchoolTwoYearEvent();
			break;
		case 29:
			CHKSUM(time, 5);
			TRCEVNT(case_seed, _T("Clock"), actor_id, _T("Clock.StartSchoolOneYearEvent"), time);
			gszRuntimeEventName = _T("Clock.StartSchoolOneYearEvent(implement)");
			StartSchoolOneYearEvent();
			break;
		case 30:
			CHKSUM(time, 6);
			TRCEVNT(case_seed, _T("Clock"), actor_id, _T("Clock.StartSchoolTwoYearEvent"), time);
			gszRuntimeEventName = _T("Clock.StartSchoolTwoYearEvent(implement)");
			StartSchoolTwoYearEvent();
			break;
		case -1:
		default:;
	}
}

double Clock::Set_actor_subsample_weight( double value )
{
	return poStateFunctions->Set_actor_subsample_weight( value);
}

double Clock::Set_actor_weight( double value )
{
	return poStateFunctions->Set_actor_weight( value);
}

int ClockStateFunctions::Set_clock_year( int value )
{

	if ( poParent->clock_year != value )
	{


		int	OldValue;
		OldValue = poParent->clock_year;
		poParent->clock_year = value;

		if ( poParent->lClockToCalibrator != NULL ) {
			if ( poParent->lClockToCalibrator->bUpdates ) {
				poParent->lClockToCalibrator->UpdateTime();
				poParent->lClockToCalibrator->poStateFunctions->Set_lCalibratorToClock_clock_year( poParent->clock_year );
			}
		}
	}
	return poParent->clock_year = value;
}

TIME ClockStateFunctions::Set_next_clock_year_end( TIME value )
{

	if ( poParent->next_clock_year_end != value )
	{


		TIME	OldValue;
		OldValue = poParent->next_clock_year_end;
		poParent->next_clock_year_end = value;

		if ( poParent->poClockYearEndEvent->cRecalc == 0 ) {
			poParent->poClockYearEndEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poClockYearEndEvent );
		};
	}
	return poParent->next_clock_year_end = value;
}

TIME ClockStateFunctions::Set_next_clock_year_start( TIME value )
{

	if ( poParent->next_clock_year_start != value )
	{


		TIME	OldValue;
		OldValue = poParent->next_clock_year_start;
		poParent->next_clock_year_start = value;

		if ( poParent->poClockYearStartEvent->cRecalc == 0 ) {
			poParent->poClockYearStartEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poClockYearStartEvent );
		};
	}
	return poParent->next_clock_year_start = value;
}

TIME ClockStateFunctions::Set_next_midyear_clock_event( TIME value )
{

	if ( poParent->next_midyear_clock_event != value )
	{


		TIME	OldValue;
		OldValue = poParent->next_midyear_clock_event;
		poParent->next_midyear_clock_event = value;

		if ( poParent->poClockMidyearEvent->cRecalc == 0 ) {
			poParent->poClockMidyearEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poClockMidyearEvent );
		};
	}
	return poParent->next_midyear_clock_event = value;
}

TIME ClockStateFunctions::Set_time_end_school_one_year( TIME value )
{

	if ( poParent->time_end_school_one_year != value )
	{


		TIME	OldValue;
		OldValue = poParent->time_end_school_one_year;
		poParent->time_end_school_one_year = value;

		if ( poParent->poEndSchoolOneYearEvent->cRecalc == 0 ) {
			poParent->poEndSchoolOneYearEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poEndSchoolOneYearEvent );
		};
	}
	return poParent->time_end_school_one_year = value;
}

TIME ClockStateFunctions::Set_time_end_school_two_year( TIME value )
{

	if ( poParent->time_end_school_two_year != value )
	{


		TIME	OldValue;
		OldValue = poParent->time_end_school_two_year;
		poParent->time_end_school_two_year = value;

		if ( poParent->poEndSchoolTwoYearEvent->cRecalc == 0 ) {
			poParent->poEndSchoolTwoYearEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poEndSchoolTwoYearEvent );
		};
	}
	return poParent->time_end_school_two_year = value;
}

TIME ClockStateFunctions::Set_time_start_school_one_year( TIME value )
{

	if ( poParent->time_start_school_one_year != value )
	{


		TIME	OldValue;
		OldValue = poParent->time_start_school_one_year;
		poParent->time_start_school_one_year = value;

		if ( poParent->poStartSchoolOneYearEvent->cRecalc == 0 ) {
			poParent->poStartSchoolOneYearEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poStartSchoolOneYearEvent );
		};
	}
	return poParent->time_start_school_one_year = value;
}

TIME ClockStateFunctions::Set_time_start_school_two_year( TIME value )
{

	if ( poParent->time_start_school_two_year != value )
	{


		TIME	OldValue;
		OldValue = poParent->time_start_school_two_year;
		poParent->time_start_school_two_year = value;

		if ( poParent->poStartSchoolTwoYearEvent->cRecalc == 0 ) {
			poParent->poStartSchoolTwoYearEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poStartSchoolTwoYearEvent );
		};
	}
	return poParent->time_start_school_two_year = value;
}

logical ClockStateFunctions::Set_actorset_filter_asClock( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asClock != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asClock;
		poParent->actorset_filter_asClock = value;

		// actorset filter
		if (poParent->actorset_filter_asClock)
		{
			poParent->poasClock = asClock->NewNode(poParent);
		}
		else
		{
			asClock->Remove(poParent->poasClock);
			poParent->poasClock = NULL;
		}
	}
	return poParent->actorset_filter_asClock = value;
}

void *Observation::operator new( size_t count )
{
	Observation *poActor;
	if ( gpoEmptyObservation != NULL ) {
		poActor = gpoEmptyObservation;
		gpoEmptyObservation = (Observation *) gpoEmptyObservation->poNext;
	}
	else {
		poActor = ::new Observation( TRUE ) ;
	}
	poActor->__finished = false;
	return poActor;
}

void Observation::InitActor()
{
	poDerivedStates = new ObservationDerivedStates(); 
	poDerivedStates->poParent = this;
	poStateFunctions = new ObservationStateFunctions();
	poStateFunctions->poParent = this;
	poStateFunctions->poDerivedStates = poDerivedStates;
	poDerivedStates->poStateFunctions = poStateFunctions;
	bUpdates = FALSE;
	nActorNumber = 2;
	InitializeStates();
}

void Observation::DeleteActor()
{
	delete poDerivedStates;
	delete poStateFunctions;
}

#pragma optimize( "", off )

void Observation::InitializeStates()
{
	// initialization of states
	actor_id = (long) 0;
	actor_subsample_weight = (double) 0;
	actor_weight = (double) 0;
	age = (TIME) 0;
	case_id = (long) 0;
	case_seed = (double) 0;
#line 80 "../code/ActorObservation.mpp"
	CHKLMT(0, fam_id, r_min_FAM_ID, r_max_FAM_ID);
	fam_id = 0 < r_min_FAM_ID ? r_min_FAM_ID : 0 > r_max_FAM_ID ? r_max_FAM_ID : (FAM_ID) 0;
#line 2662 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 82 "../code/ActorObservation.mpp"
	obs_birth = (TIME) TIME_INFINITE;
#line 2665 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 83 "../code/ActorObservation.mpp"
	obs_oldest = FALSE != 0;
#line 2668 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 78 "../code/ActorObservation.mpp"
	obs_weight = (integer) 0;
#line 2671 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 81 "../code/ActorObservation.mpp"
	pop_pool = PP_NON;
#line 2674 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	time = (TIME) gdThreadCurrentTime;
	actorset_filter_asObservationAll = FALSE;
	actorset_filter_asObservationByFam = FALSE;
	actorset_filter_asObservationByFamOldest = FALSE;
	actorset_filter_asObservationHeads = FALSE;
	actorset_filter_asObservationNonHeads = FALSE;
	actorset_filter_asObservations = FALSE;
	actorset_filter_asSimulatedObservationHeads = FALSE;
	__finished = FALSE;
	__time = (TIME) gdThreadCurrentTime;
	__events = 0;
}

#pragma optimize( "", on )

void Observation::UpdateLinkedStates()
{
}

void Observation::StartSpecial1()
{
	nSubSample = GetSubSample();
	InitializeStates();
	case_seed = GetCaseSeed();
	case_id = GetCaseID();
	lActorId = actor_id = GetObjectID();
	GetCaseWeight( &actor_weight, &actor_subsample_weight );

	// initialization of actor set nodes
	poasObservationAll = NULL;
	poasObservationByFam = NULL;
	poasObservationByFamOldest = NULL;
	poasObservationHeads = NULL;
	poasObservationNonHeads = NULL;
	poasObservations = NULL;
	poasSimulatedObservationHeads = NULL;

	// initialization of events
	//Initalize derived states
	CountActor(2);
}

#pragma optimize( "", off )

void Observation::StartSpecial2()
{
	// initialization of expression and history states
	poStateFunctions->Set_time(TIME( CoarsenMantissa( time ) ), true);
	poStateFunctions->Set_age(TIME( CoarsenMantissa( age ) ), true);
	gdThreadCurrentTime = time;
	__time = time;
#line 34 "../code/ActorObservation.mpp"
	poStateFunctions->Set_actorset_filter_asObservationAll( logical( 1 ) );
#line 2728 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 40 "../code/ActorObservation.mpp"
	poStateFunctions->Set_actorset_filter_asObservationByFam( logical( 1 ) );
#line 2731 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 43 "../code/ActorObservation.mpp"
	poStateFunctions->Set_actorset_filter_asObservationByFamOldest( logical( obs_oldest ) );
#line 2734 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 24 "../code/ActorObservation.mpp"
	poStateFunctions->Set_actorset_filter_asObservationHeads( logical( obs_oldest ) );
#line 2737 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 31 "../code/ActorObservation.mpp"
	poStateFunctions->Set_actorset_filter_asObservationNonHeads( logical( !( obs_oldest ) ) );
#line 2740 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 37 "../code/ActorObservation.mpp"
	poStateFunctions->Set_actorset_filter_asObservations( logical( 1 ) );
#line 2743 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/ActorObservation.mpp"
	poStateFunctions->Set_actorset_filter_asSimulatedObservationHeads( logical( obs_oldest&& obs_weight> 0 ) );
#line 2746 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	UpdateLinkedStates();
	poPrev = NULL;
	poNext = gpoEventQueue->poHeadActor;
	if ( gpoEventQueue->poHeadActor != NULL ) {
		gpoEventQueue->poHeadActor->poPrev = this;
	}
	gpoEventQueue->poHeadActor = this;
	bUpdates = TRUE;
}

void Observation::FinishSpecial()
{
	CountFinishedActor(2);
	bUpdates = FALSE;
	if (poasObservationAll != NULL)
	{	asObservationAll[pop_pool]->Remove(poasObservationAll);
		poasObservationAll = NULL;
	}
	if (poasObservationByFam != NULL)
	{	asObservationByFam[RANGE_POS(FAM_ID,fam_id)]->Remove(poasObservationByFam);
		poasObservationByFam = NULL;
	}
	if (poasObservationByFamOldest != NULL)
	{	asObservationByFamOldest[RANGE_POS(FAM_ID,fam_id)]->Remove(poasObservationByFamOldest);
		poasObservationByFamOldest = NULL;
	}
	if (poasObservationHeads != NULL)
	{	asObservationHeads[pop_pool]->Remove(poasObservationHeads);
		poasObservationHeads = NULL;
	}
	if (poasObservationNonHeads != NULL)
	{	asObservationNonHeads[RANGE_POS(FAM_ID,fam_id)][pop_pool]->Remove(poasObservationNonHeads);
		poasObservationNonHeads = NULL;
	}
	if (poasObservations != NULL)
	{	asObservations->Remove(poasObservations);
		poasObservations = NULL;
	}
	if (poasSimulatedObservationHeads != NULL)
	{	asSimulatedObservationHeads[pop_pool]->Remove(poasSimulatedObservationHeads);
		poasSimulatedObservationHeads = NULL;
	}

	// Remove actor from master list of active actors
	if ( poPrev != NULL )
	{
		poPrev->poNext = poNext;
	}
	else
	{
		if (gpoEventQueue->poHeadActor != NULL && gpoEventQueue->poHeadActor == this)
		{
			gpoEventQueue->poHeadActor = poNext;
		}
	}
	if ( poNext != NULL )
	{
		poNext->poPrev = poPrev;
	}
	poPrev = NULL;

	// Add actor to head of empty list of actors of this type
	poNext = gpoEmptyObservation;
	gpoEmptyObservation = this;
}

#pragma optimize( "", on )

void Observation::WaitUntil( double target_time )
{
	//don't increase gdEventsForTabulation if special events because they are all treated as one event for tabulation purposes
	if (!gpoEventQueue->m_bSpecialEvents)
	{
		gdEventsForTabulation++;
	}
	__events = gdEventsForTabulation;
	if ( time < (TIME) target_time )
	{
		__time = (TIME) target_time;
		WaitSpecial( (TIME) target_time - time );
		poStateFunctions->Set_time( (TIME) target_time );
	}
}

Observation *Observation::UpdateTime(  )
{
	__events = gdEventsForTabulation;
	if ( time < (TIME) gdThreadCurrentTime ) {
		__time = (TIME) gdThreadCurrentTime;
		WaitSpecial( (TIME) gdThreadCurrentTime - time );
		poStateFunctions->Set_time( (TIME) gdThreadCurrentTime );
	}
	return this;
}

void Observation::WaitUntilThisActor( double target_time )
{
	WaitUntil(target_time);
}

void Observation::WaitSpecial( TIME wait_time )
{

	// increment the age
		poStateFunctions->Set_age( age + wait_time );
}

double Observation::EventTime( int nEventNum, int *pnEventInfo )
{
	double	event_time = double( TIME_INFINITE );

	gdRuntimeTime = (double) time;
	gnRuntimeActorId = actor_id;
	switch( nEventNum ) {
		case -1:
		default:;
	}
	return (TIME) event_time;
}

void Observation::Implement( int nEventNum, int nEventInfo )
{
	gdThreadEvents += 1;
	gdRuntimeTime = (double) time;
	gnRuntimeActorId = actor_id;
	switch( nEventNum )
	{
		case -1:
		default:;
	}
}

double Observation::Set_actor_subsample_weight( double value )
{
	return poStateFunctions->Set_actor_subsample_weight( value);
}

double Observation::Set_actor_weight( double value )
{
	return poStateFunctions->Set_actor_weight( value);
}

FAM_ID ObservationStateFunctions::Set_fam_id( FAM_ID value )
{

	CHKLMT(value, fam_id, r_min_FAM_ID, r_max_FAM_ID);
	if ( value < r_min_FAM_ID ) {
		value = r_min_FAM_ID;
	}
	else if ( value > r_max_FAM_ID ) {
		value = r_max_FAM_ID;
	}
	if ( poParent->fam_id != value )
	{


		FAM_ID	OldValue;
		OldValue = poParent->fam_id;
		poParent->fam_id = value;

		// actorset dimension or order state
		if (poParent->actorset_filter_asObservationByFam)
		{
			asObservationByFam[RANGE_POS(FAM_ID,OldValue)]->Remove(poParent->poasObservationByFam);
			poParent->poasObservationByFam = asObservationByFam[RANGE_POS(FAM_ID, poParent->fam_id)]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asObservationByFamOldest)
		{
			asObservationByFamOldest[RANGE_POS(FAM_ID,OldValue)]->Remove(poParent->poasObservationByFamOldest);
			poParent->poasObservationByFamOldest = asObservationByFamOldest[RANGE_POS(FAM_ID, poParent->fam_id)]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asObservationNonHeads)
		{
			asObservationNonHeads[RANGE_POS(FAM_ID,OldValue)][poParent->pop_pool]->Remove(poParent->poasObservationNonHeads);
			poParent->poasObservationNonHeads = asObservationNonHeads[RANGE_POS(FAM_ID, poParent->fam_id)][poParent->pop_pool]->NewNode(poParent);
		}
	}
	return poParent->fam_id = value;
}

logical ObservationStateFunctions::Set_obs_oldest( logical value )
{

	value = value != 0;
	if ( poParent->obs_oldest != value )
	{
#line 43 "../code/ActorObservation.mpp"
		Set_actorset_filter_asObservationByFamOldest( value );
#line 2935 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 24 "../code/ActorObservation.mpp"
		Set_actorset_filter_asObservationHeads( value );
#line 2938 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 31 "../code/ActorObservation.mpp"
		Set_actorset_filter_asObservationNonHeads( !( value ) );
#line 2941 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/ActorObservation.mpp"
		Set_actorset_filter_asSimulatedObservationHeads( value&& poParent->obs_weight> 0 );
#line 2944 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		logical	OldValue;
		OldValue = poParent->obs_oldest;
		poParent->obs_oldest = value;

	}
	return poParent->obs_oldest = value;
}

integer ObservationStateFunctions::Set_obs_weight( integer value )
{

	if ( poParent->obs_weight != value )
	{
#line 28 "../code/ActorObservation.mpp"
		Set_actorset_filter_asSimulatedObservationHeads( poParent->obs_oldest&& value> 0 );
#line 2962 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		integer	OldValue;
		OldValue = poParent->obs_weight;
		poParent->obs_weight = value;

	}
	return poParent->obs_weight = value;
}

POP_POOL ObservationStateFunctions::Set_pop_pool( POP_POOL value )
{

	CHKLMT(value, pop_pool, 0, rc_size_POP_POOL - 1);
	if ( value < 0 ) {
		value = (POP_POOL) 0;
	}
	else if ( value >= rc_size_POP_POOL ) {
		value = (POP_POOL) ( rc_size_POP_POOL - 1 );
	}
	if ( poParent->pop_pool != value )
	{


		POP_POOL	OldValue;
		OldValue = poParent->pop_pool;
		poParent->pop_pool = value;

		// actorset dimension or order state
		if (poParent->actorset_filter_asObservationAll)
		{
			asObservationAll[OldValue]->Remove(poParent->poasObservationAll);
			poParent->poasObservationAll = asObservationAll[poParent->pop_pool]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asObservationHeads)
		{
			asObservationHeads[OldValue]->Remove(poParent->poasObservationHeads);
			poParent->poasObservationHeads = asObservationHeads[poParent->pop_pool]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asObservationNonHeads)
		{
			asObservationNonHeads[RANGE_POS(FAM_ID, poParent->fam_id)][OldValue]->Remove(poParent->poasObservationNonHeads);
			poParent->poasObservationNonHeads = asObservationNonHeads[RANGE_POS(FAM_ID, poParent->fam_id)][poParent->pop_pool]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asSimulatedObservationHeads)
		{
			asSimulatedObservationHeads[OldValue]->Remove(poParent->poasSimulatedObservationHeads);
			poParent->poasSimulatedObservationHeads = asSimulatedObservationHeads[poParent->pop_pool]->NewNode(poParent);
		}
	}
	return poParent->pop_pool = value;
}

logical ObservationStateFunctions::Set_actorset_filter_asObservationAll( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asObservationAll != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asObservationAll;
		poParent->actorset_filter_asObservationAll = value;

		// actorset filter
		if (poParent->actorset_filter_asObservationAll)
		{
			poParent->poasObservationAll = asObservationAll[poParent->pop_pool]->NewNode(poParent);
		}
		else
		{
			asObservationAll[poParent->pop_pool]->Remove(poParent->poasObservationAll);
			poParent->poasObservationAll = NULL;
		}
	}
	return poParent->actorset_filter_asObservationAll = value;
}

logical ObservationStateFunctions::Set_actorset_filter_asObservationByFam( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asObservationByFam != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asObservationByFam;
		poParent->actorset_filter_asObservationByFam = value;

		// actorset filter
		if (poParent->actorset_filter_asObservationByFam)
		{
			poParent->poasObservationByFam = asObservationByFam[RANGE_POS(FAM_ID, poParent->fam_id)]->NewNode(poParent);
		}
		else
		{
			asObservationByFam[RANGE_POS(FAM_ID, poParent->fam_id)]->Remove(poParent->poasObservationByFam);
			poParent->poasObservationByFam = NULL;
		}
	}
	return poParent->actorset_filter_asObservationByFam = value;
}

logical ObservationStateFunctions::Set_actorset_filter_asObservationByFamOldest( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asObservationByFamOldest != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asObservationByFamOldest;
		poParent->actorset_filter_asObservationByFamOldest = value;

		// actorset filter
		if (poParent->actorset_filter_asObservationByFamOldest)
		{
			poParent->poasObservationByFamOldest = asObservationByFamOldest[RANGE_POS(FAM_ID, poParent->fam_id)]->NewNode(poParent);
		}
		else
		{
			asObservationByFamOldest[RANGE_POS(FAM_ID, poParent->fam_id)]->Remove(poParent->poasObservationByFamOldest);
			poParent->poasObservationByFamOldest = NULL;
		}
	}
	return poParent->actorset_filter_asObservationByFamOldest = value;
}

logical ObservationStateFunctions::Set_actorset_filter_asObservationHeads( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asObservationHeads != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asObservationHeads;
		poParent->actorset_filter_asObservationHeads = value;

		// actorset filter
		if (poParent->actorset_filter_asObservationHeads)
		{
			poParent->poasObservationHeads = asObservationHeads[poParent->pop_pool]->NewNode(poParent);
		}
		else
		{
			asObservationHeads[poParent->pop_pool]->Remove(poParent->poasObservationHeads);
			poParent->poasObservationHeads = NULL;
		}
	}
	return poParent->actorset_filter_asObservationHeads = value;
}

logical ObservationStateFunctions::Set_actorset_filter_asObservationNonHeads( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asObservationNonHeads != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asObservationNonHeads;
		poParent->actorset_filter_asObservationNonHeads = value;

		// actorset filter
		if (poParent->actorset_filter_asObservationNonHeads)
		{
			poParent->poasObservationNonHeads = asObservationNonHeads[RANGE_POS(FAM_ID, poParent->fam_id)][poParent->pop_pool]->NewNode(poParent);
		}
		else
		{
			asObservationNonHeads[RANGE_POS(FAM_ID, poParent->fam_id)][poParent->pop_pool]->Remove(poParent->poasObservationNonHeads);
			poParent->poasObservationNonHeads = NULL;
		}
	}
	return poParent->actorset_filter_asObservationNonHeads = value;
}

logical ObservationStateFunctions::Set_actorset_filter_asObservations( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asObservations != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asObservations;
		poParent->actorset_filter_asObservations = value;

		// actorset filter
		if (poParent->actorset_filter_asObservations)
		{
			poParent->poasObservations = asObservations->NewNode(poParent);
		}
		else
		{
			asObservations->Remove(poParent->poasObservations);
			poParent->poasObservations = NULL;
		}
	}
	return poParent->actorset_filter_asObservations = value;
}

logical ObservationStateFunctions::Set_actorset_filter_asSimulatedObservationHeads( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asSimulatedObservationHeads != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asSimulatedObservationHeads;
		poParent->actorset_filter_asSimulatedObservationHeads = value;

		// actorset filter
		if (poParent->actorset_filter_asSimulatedObservationHeads)
		{
			poParent->poasSimulatedObservationHeads = asSimulatedObservationHeads[poParent->pop_pool]->NewNode(poParent);
		}
		else
		{
			asSimulatedObservationHeads[poParent->pop_pool]->Remove(poParent->poasSimulatedObservationHeads);
			poParent->poasSimulatedObservationHeads = NULL;
		}
	}
	return poParent->actorset_filter_asSimulatedObservationHeads = value;
}

void *Person::operator new( size_t count )
{
	Person *poActor;
	if ( gpoEmptyPerson != NULL ) {
		poActor = gpoEmptyPerson;
		gpoEmptyPerson = (Person *) gpoEmptyPerson->poNext;
	}
	else {
		poActor = ::new Person( TRUE ) ;
	}
	poActor->__finished = false;
	return poActor;
}

void Person::Set_lBioFather( Person *value )
{
	Person *poOldLink;

	if ( lBioFather != value ) {
		if (lBioFather != NULL)
		{
			lBioFather->UpdateTime();
		}
		if (value != NULL)
		{
			value->UpdateTime();
		}
		poOldLink = lBioFather;
		lBioFather = value;
		UpdateLinkedStates();
		if ( poOldLink != NULL ) {
			poOldLink->mlBioFatherChildren->Remove( this );
		}
		if ( value != NULL ) {
			value->mlBioFatherChildren->Add( this );
		}
	}
}

void Person::Set_lBioMother( Person *value )
{
	Person *poOldLink;

	if ( lBioMother != value ) {
		if (lBioMother != NULL)
		{
			lBioMother->UpdateTime();
		}
		if (value != NULL)
		{
			value->UpdateTime();
		}
		poOldLink = lBioMother;
		lBioMother = value;
		UpdateLinkedStates();
		if ( poOldLink != NULL ) {
			poOldLink->mlBioMotherChildren->Remove( this );
		}
		if ( value != NULL ) {
			value->mlBioMotherChildren->Add( this );
		}
	}
}

void Person::Set_lCalibrator( Calibrator *value )
{
	Calibrator *poOldLink;

	if ( lCalibrator != value ) {
		if (lCalibrator != NULL)
		{
			lCalibrator->UpdateTime();
		}
		if (value != NULL)
		{
			value->UpdateTime();
		}
		poOldLink = lCalibrator;
		lCalibrator = value;
		UpdateLinkedStates();
		if ( poOldLink != NULL ) {
			poOldLink->mlCalibratorToPerson->Remove( this );
		}
		if ( value != NULL ) {
			value->mlCalibratorToPerson->Add( this );
		}
	}
}

void Person::Set_lHHFather( Person *value )
{
	Person *poOldLink;

	if ( lHHFather != value ) {
		if (lHHFather != NULL)
		{
			lHHFather->UpdateTime();
		}
		if (value != NULL)
		{
			value->UpdateTime();
		}
		poStateFunctions->Set_family_role( ( lSpouse&& sex== MALE )? FR_SPOUSE: ( value|| lHHMother )? FR_CHILD: FR_HEAD );
		poOldLink = lHHFather;
		lHHFather = value;
		UpdateLinkedStates();
		if ( poOldLink != NULL ) {
			poOldLink->mlHHFatherChildren->Remove( this );
		}
		if ( value != NULL ) {
			value->mlHHFatherChildren->Add( this );
		}
	}
}

void Person::Set_lHHMother( Person *value )
{
	Person *poOldLink;

	if ( lHHMother != value ) {
		if (lHHMother != NULL)
		{
			lHHMother->UpdateTime();
		}
		if (value != NULL)
		{
			value->UpdateTime();
		}
		poStateFunctions->Set_family_role( ( lSpouse&& sex== MALE )? FR_SPOUSE: ( lHHFather|| value )? FR_CHILD: FR_HEAD );
		poOldLink = lHHMother;
		lHHMother = value;
		UpdateLinkedStates();
		if ( poOldLink != NULL ) {
			poOldLink->mlHHMotherChildren->Remove( this );
		}
		if ( value != NULL ) {
			value->mlHHMotherChildren->Add( this );
		}
	}
}

void Person::Set_lSpouse( Person *value )
{
	Person *poOldLink;

	if ( lSpouse != value ) {
		if (lSpouse != NULL)
		{
			lSpouse->UpdateTime();
		}
		if (value != NULL)
		{
			value->UpdateTime();
		}
		poStateFunctions->Set_family_role( ( value&& sex== MALE )? FR_SPOUSE: ( lHHFather|| lHHMother )? FR_CHILD: FR_HEAD );
		poStateFunctions->Set_has_spouse( ( value )? TRUE: FALSE );
		poOldLink = lSpouse;
		lSpouse = value;
		UpdateLinkedStates();
		if ( poOldLink != NULL && poOldLink->lSpouse == this) {
			poOldLink->Set_lSpouse( NULL );
		}
		if ( value != NULL && value->lSpouse != this) {
			value->Set_lSpouse( this );
		}
	}
}

void Person::Set_lStartValues( StartpopValues *value )
{
	StartpopValues *poOldLink;

	if ( lStartValues != value ) {
		if (lStartValues != NULL)
		{
			lStartValues->UpdateTime();
		}
		if (value != NULL)
		{
			value->UpdateTime();
		}
		poOldLink = lStartValues;
		lStartValues = value;
		UpdateLinkedStates();
		if ( poOldLink != NULL && poOldLink->lPersonStartpop == this) {
			poOldLink->Set_lPersonStartpop( NULL );
		}
		if ( value != NULL && value->lPersonStartpop != this) {
			value->Set_lPersonStartpop( this );
		}
	}
}

void Person::MidYearIMPLEMENT_HOOK()
{
}

void Person::YearEndIMPLEMENT_HOOK()
{
}

void Person::YearStartIMPLEMENT_HOOK()
{
}

void Person::BirthdayEventIMPLEMENT_HOOK()
{
}

void Person::InitActor()
{
	poDerivedStates = new PersonDerivedStates(); 
	poDerivedStates->poParent = this;
	poStateFunctions = new PersonStateFunctions();
	poStateFunctions->poParent = this;
	poStateFunctions->poDerivedStates = poDerivedStates;
	poDerivedStates->poStateFunctions = poStateFunctions;
	bUpdates = FALSE;
	nActorNumber = 3;
	InitializeStates();
	lBioFather = NULL;
	lBioMother = NULL;
	lCalibrator = NULL;
	lHHFather = NULL;
	lHHMother = NULL;
	lSpouse = NULL;
	lStartValues = NULL;
	mlBioFatherChildren = new Person_mlBioFatherChildren( this );
	mlBioMotherChildren = new Person_mlBioMotherChildren( this );
	mlHHFatherChildren = new Person_mlHHFatherChildren( this );
	mlHHMotherChildren = new Person_mlHHMotherChildren( this );
}

void Person::DeleteActor()
{
	delete poDerivedStates;
	delete poStateFunctions;
	delete mlBioFatherChildren;
	delete mlBioMotherChildren;
	delete mlHHFatherChildren;
	delete mlHHMotherChildren;
}

#pragma optimize( "", off )

void Person::InitializeStates()
{
	// initialization of states
	actor_id = (long) 0;
	actor_subsample_weight = (double) 0;
	actor_weight = (double) 0;
#line 100 "../code/HumanCapitalIndex.mpp"
	adult_survival = (double) 0.0;
#line 3451 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	age = (TIME) 0;
	calendar_year = r_min_ALL_YEAR_RANGE;
	case_id = (long) 0;
	case_seed = (double) 0;
#line 69 "../code/ActorPerson.mpp"
	creation_type = CT_START;
#line 3458 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 97 "../code/EducationPrimaryTransmission.mpp"
	educ_mother = EOL_LOW;
#line 3461 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 97 "../code/EducationPrimaryTracking.mpp"
	CHKLMT(6, educ_one_entry_age, r_min_EDUC_ONE_ENTRY_AGE, r_max_EDUC_ONE_ENTRY_AGE);
	educ_one_entry_age = 6 < r_min_EDUC_ONE_ENTRY_AGE ? r_min_EDUC_ONE_ENTRY_AGE : 6 > r_max_EDUC_ONE_ENTRY_AGE ? r_max_EDUC_ONE_ENTRY_AGE : (EDUC_ONE_ENTRY_AGE) 6;
#line 3465 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 65 "../code/EducationPrimaryBase.mpp"
	educ_one_fate = EOL_LOW;
#line 3468 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 99 "../code/EducationPrimaryTracking.mpp"
	educ_one_geo = EOG_OO;
#line 3471 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 94 "../code/EducationPrimaryTracking.mpp"
	CHKLMT(0, educ_one_grade_attended, r_min_EDUC_ONE_GRADE0, r_max_EDUC_ONE_GRADE0);
	educ_one_grade_attended = 0 < r_min_EDUC_ONE_GRADE0 ? r_min_EDUC_ONE_GRADE0 : 0 > r_max_EDUC_ONE_GRADE0 ? r_max_EDUC_ONE_GRADE0 : (EDUC_ONE_GRADE0) 0;
#line 3475 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 98 "../code/EducationPrimaryTracking.mpp"
	CHKLMT(1, educ_one_grade_fate, r_min_EDUC_ONE_GRADE, r_max_EDUC_ONE_GRADE);
	educ_one_grade_fate = 1 < r_min_EDUC_ONE_GRADE ? r_min_EDUC_ONE_GRADE : 1 > r_max_EDUC_ONE_GRADE ? r_max_EDUC_ONE_GRADE : (EDUC_ONE_GRADE) 1;
#line 3479 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 95 "../code/EducationPrimaryTracking.mpp"
	CHKLMT(0, educ_one_grade_passed, r_min_EDUC_ONE_GRADE0, r_max_EDUC_ONE_GRADE0);
	educ_one_grade_passed = 0 < r_min_EDUC_ONE_GRADE0 ? r_min_EDUC_ONE_GRADE0 : 0 > r_max_EDUC_ONE_GRADE0 ? r_max_EDUC_ONE_GRADE0 : (EDUC_ONE_GRADE0) 0;
#line 3483 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 100 "../code/EducationPrimaryTracking.mpp"
	educ_one_group = E1G_OO;
#line 3486 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 96 "../code/EducationPrimaryTracking.mpp"
	educ_one_status = EOS_NEVER;
#line 3489 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 93 "../code/EducationPrimaryTracking.mpp"
	educ_one_to_process = FALSE != 0;
#line 3492 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 107 "../code/EducationSecondaryBase.mpp"
	educ_two_delay = (int) 0;
#line 3495 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 118 "../code/EducationSecondaryBase.mpp"
	educ_two_geo = ETG_OO;
#line 3498 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 116 "../code/EducationSecondaryBase.mpp"
	CHKLMT(0, educ_two_grade_attended, r_min_EDUC_TWO_GRADE0, r_max_EDUC_TWO_GRADE0);
	educ_two_grade_attended = 0 < r_min_EDUC_TWO_GRADE0 ? r_min_EDUC_TWO_GRADE0 : 0 > r_max_EDUC_TWO_GRADE0 ? r_max_EDUC_TWO_GRADE0 : (EDUC_TWO_GRADE0) 0;
#line 3502 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 117 "../code/EducationSecondaryBase.mpp"
	CHKLMT(0, educ_two_grade_passed, r_min_EDUC_TWO_GRADE0, r_max_EDUC_TWO_GRADE0);
	educ_two_grade_passed = 0 < r_min_EDUC_TWO_GRADE0 ? r_min_EDUC_TWO_GRADE0 : 0 > r_max_EDUC_TWO_GRADE0 ? r_max_EDUC_TWO_GRADE0 : (EDUC_TWO_GRADE0) 0;
#line 3506 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 119 "../code/EducationSecondaryBase.mpp"
	educ_two_group = E2G_OO;
#line 3509 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 109 "../code/EducationSecondaryBase.mpp"
	educ_two_status = ETS_NEVER;
#line 3512 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 108 "../code/EducationSecondaryBase.mpp"
	educ_two_to_process = FALSE != 0;
#line 3515 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 42 "../code/Ethnicity.mpp"
	ethnicity = ETHNO_00;
#line 3518 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 76 "../code/ActorPerson.mpp"
	ever_resident = FALSE != 0;
#line 3521 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 90 "../code/FamilyFirstUnion.mpp"
	ever_union = FALSE != 0;
#line 3524 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 48 "../code/MigrationGeneral.mpp"
	first_move_done = FALSE != 0;
#line 3527 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 72 "../code/ActorPerson.mpp"
	geo = GEO_00;
#line 3530 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 71 "../code/ActorPerson.mpp"
	geo_birth = GEO_00;
#line 3533 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 73 "../code/ActorPerson.mpp"
	geo_prev = GEO_00;
#line 3536 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 53 "../code/MigrationGeneral.mpp"
	geo_want_to_move = GEO_00;
#line 3539 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 54 "../code/ChildVaccination.mpp"
	got_prenat_care = GPC_NO;
#line 3542 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 106 "../code/HumanCapitalIndex.mpp"
	ind_hci = (double) 0.0;
#line 3545 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 13 "../code/ClockEvents.mpp"
	CHKLMT(0, integer_age, r_min_AGE_RANGE, r_max_AGE_RANGE);
	integer_age = 0 < r_min_AGE_RANGE ? r_min_AGE_RANGE : 0 > r_max_AGE_RANGE ? r_max_AGE_RANGE : (AGE_RANGE) 0;
#line 3549 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 75 "../code/ActorPerson.mpp"
	is_alive = FALSE != 0;
#line 3552 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 116 "../code/FamilyFemalePartnershipStatus.mpp"
	is_blocked_from_marriage = FALSE != 0;
#line 3555 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 53 "../code/ChildVaccination.mpp"
	is_immunized = FALSE != 0;
#line 3558 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 77 "../code/ActorPerson.mpp"
	is_ready_for_birthtables = FALSE != 0;
#line 3561 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 30 "../code/Stunting.mpp"
	is_stunted = FALSE != 0;
#line 3564 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 81 "../code/MortalityDetailed.mpp"
	mother_age_at_birth = (double) TIME_INFINITE;
#line 3567 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 101 "../code/FertilityDetailed.mpp"
	CHKLMT(0, parity, r_min_PARITY_RANGE, r_max_PARITY_RANGE);
	parity = 0 < r_min_PARITY_RANGE ? r_min_PARITY_RANGE : 0 > r_max_PARITY_RANGE ? r_max_PARITY_RANGE : (PARITY_RANGE) 0;
#line 3571 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 39 "../code/EducationPreSchool.mpp"
	preschool_is_decided = FALSE != 0;
#line 3574 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 91 "../code/HumanCapitalIndex.mpp"
	quality_adjusted_schooling = (double) 0.0;
#line 3577 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 88 "../code/HumanCapitalIndex.mpp"
	quality_of_schooling = (double) 0.0;
#line 3580 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 81 "../code/ActorPerson.mpp"
	ready_to_set_alive = FALSE != 0;
#line 3583 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 67 "../code/ActorPerson.mpp"
	sex = FEMALE;
#line 3586 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	table_unit_TabEducFateByGroup = (counter) 0;
	table_unit_TabEducFateDistrYob = (counter) 0;
	table_unit_TabPrimarySchoolPlanning = (counter) 0;
	table_unit_TabSchool2TrackTab = (counter) 0;
	table_unit_tabStuntingSexRegMotherYob = (counter) 0;
	table_unit_tabStuntingSexRegYob = (counter) 0;
#line 54 "../code/FertilityGeneral.mpp"
	this_parity_spell = FALSE != 0;
#line 3595 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	time = (TIME) gdThreadCurrentTime;
#line 45 "../code/MigrationGeneral.mpp"
	time_first_move = (TIME) TIME_INFINITE;
#line 3599 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 46 "../code/MigrationGeneral.mpp"
	time_last_move = (TIME) TIME_INFINITE;
#line 3602 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 14 "../code/ClockEvents.mpp"
	time_next_birthday = (TIME) TIME_INFINITE;
#line 3605 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 66 "../code/ActorPerson.mpp"
	time_of_birth = (TIME) 0.0;
#line 3608 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 50 "../code/MigrationGeneral.mpp"
	time_of_first_immigration = (TIME) TIME_INFINITE;
#line 3611 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 37 "../code/EducationPreSchool.mpp"
	years_preschool = (int) 0;
#line 3614 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	poDerivedStates->active_spell_duration_this_parity_spell_TRUE = (TIME) 0;
	poDerivedStates->active_spell_duration_this_parity_spell_TRUE_time = (TIME) gdThreadCurrentTime;
	actorset_filter_asAllFertilePersonsForFertilityAlignment = FALSE;
	actorset_filter_asAllFertilePersonsForFertilityAlignmentByAge = FALSE;
	actorset_filter_asAllPerson = FALSE;
	actorset_filter_asAvailableMale = FALSE;
	actorset_filter_asAvailableMaleByType = FALSE;
	actorset_filter_asAvailableMaleByTypeImmiScratch = FALSE;
	actorset_filter_asAvailableMaleImmiScratch = FALSE;
	actorset_filter_asFemaleInUnionByAgeAndPartnerAge = FALSE;
	actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch = FALSE;
	actorset_filter_asPotentialImmigrantMothers = FALSE;
	actorset_filter_asPotentialSchoolOneStudent = FALSE;
	actorset_filter_asPotentialSchoolOneStudentToProcess = FALSE;
	actorset_filter_asPotentialSchoolTwoStudent = FALSE;
	actorset_filter_asPotentialSchoolTwoStudentToProcess = FALSE;
	actorset_filter_asResidentBabies = FALSE;
	actorset_filter_asResidentsAge0Sex = FALSE;
	actorset_filter_asResidentsAge0SexGeo = FALSE;
	actorset_filter_asSimBornAge0 = FALSE;
	actorset_filter_asWantToMove = FALSE;
	actorset_filter_asWomenNoChildren = FALSE;
	actorset_filter_asWomenWithChildren = FALSE;
	age_last_birth = (double) 0;
	age_mig = r_min_AGE_MIG;
	age_youngest_child = (double) 0;
	poDerivedStates->aggregate_ethnicity_ETHNICITY_SHORT = (ETHNICITY_SHORT) 0;
	poDerivedStates->aggregate_geo_birth_GEO_NAT = (GEO_NAT) 0;
	poDerivedStates->aggregate_geo_birth_REGION_INT = (REGION_INT) 0;
	poDerivedStates->aggregate_geo_GEO_NAT = (GEO_NAT) 0;
	poDerivedStates->aggregate_geo_REGION_INT = (REGION_INT) 0;
	poDerivedStates->aggregate_region_birth_REGION_NAT = (REGION_NAT) 0;
	poDerivedStates->aggregate_region_int_REGION_NAT = (REGION_NAT) 0;
	all_year = r_min_ALL_YEAR_RANGE;
	child_agegr = (CHILD_AGEGR) 0;
	child_agegr_part = (int) 0;
	child_mortality = (double) 0;
	child_mortality_group = (CHILD_MORTALITY_GROUP) 0;
	children_in_household = (short) 0;
	poDerivedStates->completed_spell_duration_this_parity_spell_TRUE = (TIME) 0;
	poDerivedStates->duration_ = (TIME) 0;
	poDerivedStates->duration__time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_educ_one_fate_EOL_HIGH = (TIME) 0;
	poDerivedStates->duration_educ_one_fate_EOL_HIGH_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_educ_one_fate_EOL_LOW = (TIME) 0;
	poDerivedStates->duration_educ_one_fate_EOL_LOW_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_educ_one_fate_EOL_MEDIUM = (TIME) 0;
	poDerivedStates->duration_educ_one_fate_EOL_MEDIUM_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_educ_one_level_EOL_HIGH = (TIME) 0;
	poDerivedStates->duration_educ_one_level_EOL_HIGH_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_educ_one_level_EOL_LOW = (TIME) 0;
	poDerivedStates->duration_educ_one_level_EOL_LOW_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_educ_one_level_EOL_MEDIUM = (TIME) 0;
	poDerivedStates->duration_educ_one_level_EOL_MEDIUM_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_educ_two_level_ETL_GRADUATE = (TIME) 0;
	poDerivedStates->duration_educ_two_level_ETL_GRADUATE_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_educ_two_level_ETL_INCOMPLETE = (TIME) 0;
	poDerivedStates->duration_educ_two_level_ETL_INCOMPLETE_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_educ_two_level_ETL_NEVER = (TIME) 0;
	poDerivedStates->duration_educ_two_level_ETL_NEVER_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_ever_union_TRUE = (TIME) 0;
	poDerivedStates->duration_ever_union_TRUE_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_got_prenat_care_GPC_YES = (TIME) 0;
	poDerivedStates->duration_got_prenat_care_GPC_YES_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_has_spouse_TRUE = (TIME) 0;
	poDerivedStates->duration_has_spouse_TRUE_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_is_aded0_TRUE = (TIME) 0;
	poDerivedStates->duration_is_aded0_TRUE_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_is_aded1_TRUE = (TIME) 0;
	poDerivedStates->duration_is_aded1_TRUE_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_is_aded2_TRUE = (TIME) 0;
	poDerivedStates->duration_is_aded2_TRUE_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_is_child_TRUE = (TIME) 0;
	poDerivedStates->duration_is_child_TRUE_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_is_immunized_TRUE = (TIME) 0;
	poDerivedStates->duration_is_immunized_TRUE_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_is_old_TRUE = (TIME) 0;
	poDerivedStates->duration_is_old_TRUE_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_tab_out_of_primary_TRUE = (TIME) 0;
	poDerivedStates->duration_tab_out_of_primary_TRUE_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_years_preschool_0 = (TIME) 0;
	poDerivedStates->duration_years_preschool_0_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_years_preschool_1 = (TIME) 0;
	poDerivedStates->duration_years_preschool_1_time = (TIME) gdThreadCurrentTime;
	poDerivedStates->duration_years_preschool_2 = (TIME) 0;
	poDerivedStates->duration_years_preschool_2_time = (TIME) gdThreadCurrentTime;
	educ1_group = (EDUC1_GROUP) 0;
	educ_one_grade = r_min_EDUC_ONE_GRADE;
	educ_one_level = (EDUC_ONE_LEVEL) 0;
	educ_one_required_rooms = (double) 0;
	educ_one_required_rooms_value_out = (double) 0;
	educ_one_required_teachers = (double) 0;
	educ_one_required_teachers_value_out = (double) 0;
	educ_two_level = (EDUC_TWO_LEVEL) 0;
	poDerivedStates->entrances_educ_one_grade_attended_1 = (counter) 0;
	poDerivedStates->entrances_has_passed_primary_TRUE = (counter) 0;
	ethnicity_short = (ETHNICITY_SHORT) 0;
	family_role = (FAM_ROLE) 0;
	female_spouse_age = r_min_FEMALE_SPOUSE_AGE;
	fertile_age = r_min_FERTILE_AGE_RANGE;
	geo_nat = (GEO_NAT) 0;
	has_passed_primary = FALSE;
	has_spouse = FALSE;
	in_hci_sample = FALSE;
	in_projected_time = FALSE;
	in_sel_years = FALSE;
	in_selected_year = FALSE;
	is_aded0 = FALSE;
	is_aded1 = FALSE;
	is_aded2 = FALSE;
	is_adult = FALSE;
	is_alive_value_in = FALSE;
	is_child = FALSE;
	is_educ_one_student = FALSE;
	is_fertile = FALSE;
	is_mortal = FALSE;
	is_old = FALSE;
	is_resident = FALSE;
	is_stunted_value_in = FALSE;
	lives_with_dependent_child = FALSE;
	male_age = r_min_MALE_SPOUSE_AGE;
	male_spouse_age = r_min_MALE_SPOUSE_AGE;
	moth_agegr = (MOTH_AGEGR) 0;
	moth_agegr_part = (int) 0;
	region_birth = (REGION_INT) 0;
	region_int = (REGION_INT) 0;
	region_nat = (REGION_NAT) 0;
	sel_years = (TAB_SEL_YEARS) 0;
	selected_year10 = (SELECTED_YEAR10) 0;
	sim_year = r_min_SIM_YEAR_RANGE;
	poDerivedStates->split_age_last_birth_MOTH_AGEGR_PART = (MOTH_AGEGR_PART) 0;
	poDerivedStates->split_age_youngest_child_CHILD_AGEGR_PART = (CHILD_AGEGR_PART) 0;
	poDerivedStates->split_integer_age_TAB_AGEC5 = (TAB_AGEC5) 0;
	spouse_group = (SPOUSE_GROUP) 0;
	poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY = (DUR_TIME_IN_PARITY) 0;
	poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time = (TIME) TIME_INFINITE;
	poDerivedStates->sssplit_age_AGE_MIG_PART = (AGE_MIG_PART) 0;
	poDerivedStates->sssplit_age_AGE_MIG_PART_time = (TIME) TIME_INFINITE;
	survived_early_years = FALSE;
	tab_out_of_primary = FALSE;
	tab_sim_yob = r_min_SIM_YEAR_RANGE;
	table_filter_PopPyramidByEduc = FALSE;
	table_filter_tabBirthsYearPlace = FALSE;
	table_filter_TabChildVaccination = FALSE;
	table_filter_TabEduc15ByDistrict = FALSE;
	table_filter_TabEduc15ByDistrictBirth = FALSE;
	table_filter_tabEducationFateGeobirYob = FALSE;
	table_filter_TabEducFateByGroup = FALSE;
	table_filter_TabEducFateDistrYob = FALSE;
	table_filter_tabHavingSpouse = FALSE;
	table_filter_tabHCI = FALSE;
	table_filter_TabHCIDistrict = FALSE;
	table_filter_TabImmunizationChildren = FALSE;
	table_filter_tabMigrationOriginDestination = FALSE;
	table_filter_TabPopProvAgeEducSex = FALSE;
	table_filter_tabPopulationYearPlace = FALSE;
	table_filter_TabPrenatCare = FALSE;
	table_filter_tabPreSchool = FALSE;
	table_filter_TabPrimarySchoolPlanning = FALSE;
	table_filter_TabPrimSchoolEntries = FALSE;
	table_filter_TabPrimSchoolGraduations = FALSE;
	table_filter_TabPrimSchoolOutOfSchool9to11 = FALSE;
	table_filter_TabSchool2AttainmentsTab = FALSE;
	table_filter_TabSchool2TrackTab = FALSE;
	table_filter_tabStuntingSexRegMotherYob = FALSE;
	table_filter_tabStuntingSexRegYob = FALSE;
	time_in_parity = (int) 0;
	poDerivedStates->transitions_integer_age_14_15 = (counter) 0;
	poDerivedStates->transitions_is_alive_FALSE_TRUE = (counter) 0;
	poDerivedStates->transitions_is_alive_TRUE_FALSE = (counter) 0;
	poDerivedStates->transitions_survived_early_years_FALSE_TRUE = (counter) 0;
	poDerivedStates->trigger_changes_calendar_year = (logical) 0;
	poDerivedStates->trigger_changes_calendar_year_time = (TIME) TIME_INFINITE;
	poDerivedStates->trigger_changes_geo = (logical) 0;
	poDerivedStates->trigger_changes_geo_time = (TIME) TIME_INFINITE;
	poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE = (logical) 0;
	poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE_time = (TIME) TIME_INFINITE;
	union1_group = (UNION1_GROUP) 0;
	poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_adult_survival = (real) 0;
	poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_ind_hci = (real) 0;
	poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_is_stunted = (integer) 0;
	poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling = (real) 0;
	poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_survived_early_years = (integer) 0;
	poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling = (integer) 0;
	year_of_birth = r_min_ALL_YEAR_RANGE;
	year_of_first_immigration = r_min_SIM_YEAR_RANGE;
	years_of_schooling = (int) 0;
	count_mlHHFatherChildren = (counter) 0;
	count_mlHHMotherChildren = (counter) 0;
	lCalibrator_mort_female_0 = (double) 0;
	lCalibrator_mort_female_1 = (double) 0;
	lCalibrator_mort_female_2 = (double) 0;
	lCalibrator_mort_female_3 = (double) 0;
	lCalibrator_mort_female_4 = (double) 0;
	lCalibrator_mort_male_0 = (double) 0;
	lCalibrator_mort_male_1 = (double) 0;
	lCalibrator_mort_male_2 = (double) 0;
	lCalibrator_mort_male_3 = (double) 0;
	lCalibrator_mort_male_4 = (double) 0;
	lSpouse_educ_one_fate = (EDUC_ONE_LEVEL) 0;
	lSpouse_integer_age = r_min_AGE_RANGE;
	min_over_mlHHMotherChildren_integer_age = r_min_AGE_RANGE;
	sum_over_mlHHFatherChildren_is_alive = (integer) 0;
	sum_over_mlHHMotherChildren_is_alive = (integer) 0;
	__finished = FALSE;
	__time = (TIME) gdThreadCurrentTime;
	__events = 0;
}

#pragma optimize( "", on )

void Person::UpdateLinkedStates()
{
	if ( lCalibrator != NULL ) {
		poStateFunctions->Set_lCalibrator_mort_female_0( lCalibrator->mort_female_0 );
		poStateFunctions->Set_lCalibrator_mort_female_1( lCalibrator->mort_female_1 );
		poStateFunctions->Set_lCalibrator_mort_female_2( lCalibrator->mort_female_2 );
		poStateFunctions->Set_lCalibrator_mort_female_3( lCalibrator->mort_female_3 );
		poStateFunctions->Set_lCalibrator_mort_female_4( lCalibrator->mort_female_4 );
		poStateFunctions->Set_lCalibrator_mort_male_0( lCalibrator->mort_male_0 );
		poStateFunctions->Set_lCalibrator_mort_male_1( lCalibrator->mort_male_1 );
		poStateFunctions->Set_lCalibrator_mort_male_2( lCalibrator->mort_male_2 );
		poStateFunctions->Set_lCalibrator_mort_male_3( lCalibrator->mort_male_3 );
		poStateFunctions->Set_lCalibrator_mort_male_4( lCalibrator->mort_male_4 );
	}
	else {
		poStateFunctions->Set_lCalibrator_mort_female_0( double( 0 ) );
		poStateFunctions->Set_lCalibrator_mort_female_1( double( 0 ) );
		poStateFunctions->Set_lCalibrator_mort_female_2( double( 0 ) );
		poStateFunctions->Set_lCalibrator_mort_female_3( double( 0 ) );
		poStateFunctions->Set_lCalibrator_mort_female_4( double( 0 ) );
		poStateFunctions->Set_lCalibrator_mort_male_0( double( 0 ) );
		poStateFunctions->Set_lCalibrator_mort_male_1( double( 0 ) );
		poStateFunctions->Set_lCalibrator_mort_male_2( double( 0 ) );
		poStateFunctions->Set_lCalibrator_mort_male_3( double( 0 ) );
		poStateFunctions->Set_lCalibrator_mort_male_4( double( 0 ) );
	}
	if ( lSpouse != NULL ) {
		poStateFunctions->Set_lSpouse_educ_one_fate( lSpouse->educ_one_fate );
		poStateFunctions->Set_lSpouse_integer_age( lSpouse->integer_age );
	}
	else {
		poStateFunctions->Set_lSpouse_educ_one_fate( EDUC_ONE_LEVEL( 0 ) );
		poStateFunctions->Set_lSpouse_integer_age( AGE_RANGE( 0 ) );
	}
	if ( mlHHFatherChildren != NULL ) {
		poStateFunctions->Set_count_mlHHFatherChildren( mlHHFatherChildren->Count() );
		poStateFunctions->Set_sum_over_mlHHFatherChildren_is_alive( mlHHFatherChildren->SumOver_is_alive() );
	}
	else {
		poStateFunctions->Set_count_mlHHFatherChildren( counter( 0 ) );
		poStateFunctions->Set_sum_over_mlHHFatherChildren_is_alive( integer( 0 ) );
	}
	if ( mlHHMotherChildren != NULL ) {
		poStateFunctions->Set_count_mlHHMotherChildren( mlHHMotherChildren->Count() );
		poStateFunctions->Set_min_over_mlHHMotherChildren_integer_age( mlHHMotherChildren->MinOver_integer_age( NULL, 0 ) );
		poStateFunctions->Set_sum_over_mlHHMotherChildren_is_alive( mlHHMotherChildren->SumOver_is_alive() );
	}
	else {
		poStateFunctions->Set_count_mlHHMotherChildren( counter( 0 ) );
		poStateFunctions->Set_min_over_mlHHMotherChildren_integer_age( AGE_RANGE( 0 ) );
		poStateFunctions->Set_sum_over_mlHHMotherChildren_is_alive( integer( 0 ) );
	}
}

void Person::StartSpecial1()
{
	nSubSample = GetSubSample();
	InitializeStates();
	case_seed = GetCaseSeed();
	case_id = GetCaseID();
	lActorId = actor_id = GetObjectID();
	GetCaseWeight( &actor_weight, &actor_subsample_weight );
	lBioFather = NULL;
	lBioMother = NULL;
	lCalibrator = NULL;
	lHHFather = NULL;
	lHHMother = NULL;
	lSpouse = NULL;
	lStartValues = NULL;
	mlBioFatherChildren->RemoveAll();
	mlBioMotherChildren->RemoveAll();
	mlHHFatherChildren->RemoveAll();
	mlHHMotherChildren->RemoveAll();

	// initialization of actor set nodes
	poasAllFertilePersonsForFertilityAlignment = NULL;
	poasAllFertilePersonsForFertilityAlignmentByAge = NULL;
	poasAllPerson = NULL;
	poasAvailableMale = NULL;
	poasAvailableMaleByType = NULL;
	poasAvailableMaleByTypeImmiScratch = NULL;
	poasAvailableMaleImmiScratch = NULL;
	poasFemaleInUnionByAgeAndPartnerAge = NULL;
	poasFemaleInUnionByAgeAndPartnerAgeImmiScratch = NULL;
	poasPotentialImmigrantMothers = NULL;
	poasPotentialSchoolOneStudent = NULL;
	poasPotentialSchoolOneStudentToProcess = NULL;
	poasPotentialSchoolTwoStudent = NULL;
	poasPotentialSchoolTwoStudentToProcess = NULL;
	poasResidentBabies = NULL;
	poasResidentsAge0Sex = NULL;
	poasResidentsAge0SexGeo = NULL;
	poasSimBornAge0 = NULL;
	poasWantToMove = NULL;
	poasWomenNoChildren = NULL;
	poasWomenWithChildren = NULL;

	// initialization of events
	poBackMigrationEvent = gpoEventQueue->NewEvent( this, 0, 0 );
	poBirthdayEvent = gpoEventQueue->NewEvent( this, 1, 0 );
	poDeathAtMaxLifespanEvent = gpoEventQueue->NewEvent( this, 6, 0 );
	poEmigrationEvent = gpoEventQueue->NewEvent( this, 7, 0 );
	poFertilityBaseEvent = gpoEventQueue->NewEvent( this, 10, 0 );
	poFertilityDetailedEvent = gpoEventQueue->NewEvent( this, 11, 0 );
	poFirstImmigrationFromPoolEvent = gpoEventQueue->NewEvent( this, 12, 0 );
	poFirstImmigrationFromScratchEvent = gpoEventQueue->NewEvent( this, 13, 0 );
	poFirstUnionFormationEvent = gpoEventQueue->NewEvent( this, 15, 0 );
	poLeavingHomeEvent = gpoEventQueue->NewEvent( this, 19, 0 );
	poMigrationEvent = gpoEventQueue->NewEvent( this, 20, 0 );
	poMortalityBaseEvent = gpoEventQueue->NewEvent( this, 21, 0 );
	poMortalityDetailedEvent = gpoEventQueue->NewEvent( this, 22, 0 );
	poSetAliveEvent = gpoEventQueue->NewEvent( this, 24, 0 );
	poSetPreschoolYearsEvent = gpoEventQueue->NewEvent( this, 27, 0 );
	poScheduleEvent = gpoEventQueue->NewEvent( this, 23, 255 );
	schedule_event_time = TIME_INFINITE;
	//Initalize derived states
#line 96 "../code/FamilyFemalePartnershipStatus.mpp"
	poStateFunctions->Set_age_last_birth( double( ( lives_with_dependent_child )? integer_age- age_youngest_child: TIME_INFINITE ) );
#line 3944 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 43 "../code/MigrationGeneral.mpp"
	poStateFunctions->Set_age_mig( AGE_MIG( COERCE(AGE_MIG,poDerivedStates->sssplit_age_AGE_MIG_PART) ) );
#line 3947 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	age_mig = age_mig < r_min_AGE_MIG ? r_min_AGE_MIG : age_mig > r_max_AGE_MIG ? r_max_AGE_MIG : age_mig;
#line 85 "../code/FamilyFemalePartnershipStatus.mpp"
	poStateFunctions->Set_age_youngest_child( double( ( sex== FEMALE&& count_mlHHMotherChildren> 0 )? double(min_over_mlHHMotherChildren_integer_age): TIME_INFINITE ) );
#line 3951 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	poDerivedStates->Set_aggregate_ethnicity_ETHNICITY_SHORT( ETHNICITY_To_ETHNICITY_SHORT( ethnicity ) );
	poDerivedStates->Set_aggregate_geo_birth_GEO_NAT( GEO_To_GEO_NAT( geo_birth ) );
	poDerivedStates->Set_aggregate_geo_birth_REGION_INT( GEO_To_REGION_INT( geo_birth ) );
	poDerivedStates->Set_aggregate_geo_GEO_NAT( GEO_To_GEO_NAT( geo ) );
	poDerivedStates->Set_aggregate_geo_REGION_INT( GEO_To_REGION_INT( geo ) );
	poDerivedStates->Set_aggregate_region_birth_REGION_NAT( REGION_INT_To_REGION_NAT( region_birth ) );
	poDerivedStates->Set_aggregate_region_int_REGION_NAT( REGION_INT_To_REGION_NAT( region_int ) );
#line 23 "../code/TablesPopulation.mpp"
	poStateFunctions->Set_all_year( ALL_YEAR_RANGE( COERCE(ALL_YEAR_RANGE,calendar_year) ) );
#line 3961 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	all_year = all_year < r_min_ALL_YEAR_RANGE ? r_min_ALL_YEAR_RANGE : all_year > r_max_ALL_YEAR_RANGE ? r_max_ALL_YEAR_RANGE : all_year;
#line 114 "../code/FamilyFemalePartnershipStatus.mpp"
	poStateFunctions->Set_child_agegr( CHILD_AGEGR( ( child_agegr_part== 0 )? CA00: ( child_agegr_part== 1 )? CA01: ( child_agegr_part== 2 )? CA03: ( child_agegr_part== 3 )? CA06: ( child_agegr_part== 4 )? CA09: ( child_agegr_part== 5 )? CA12: CA15 ) );
#line 3965 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 88 "../code/FamilyFemalePartnershipStatus.mpp"
	poStateFunctions->Set_child_agegr_part( int( poDerivedStates->split_age_youngest_child_CHILD_AGEGR_PART ) );
#line 3968 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 67 "../code/MortalityDetailed.mpp"
	poStateFunctions->Set_child_mortality( double( ( integer_age== 0&& sex== MALE )? double(lCalibrator_mort_male_0): ( integer_age== 1&& sex== MALE )? double(lCalibrator_mort_male_1): ( integer_age== 2&& sex== MALE )? double(lCalibrator_mort_male_2): ( integer_age== 3&& sex== MALE )? double(lCalibrator_mort_male_3): ( integer_age== 4&& sex== MALE )? double(lCalibrator_mort_male_4): ( integer_age== 0&& sex== FEMALE )? double(lCalibrator_mort_female_0): ( integer_age== 1&& sex== FEMALE )? double(lCalibrator_mort_female_1): ( integer_age== 2&& sex== FEMALE )? double(lCalibrator_mort_female_2): ( integer_age== 3&& sex== FEMALE )? double(lCalibrator_mort_female_3): ( integer_age== 4&& sex== FEMALE )? double(lCalibrator_mort_female_4): 0.0 ) );
#line 3971 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 79 "../code/MortalityDetailed.mpp"
	poStateFunctions->Set_child_mortality_group( CHILD_MORTALITY_GROUP( ( mother_age_at_birth>= 17&& educ_mother== EOL_HIGH )? CMG_00: ( mother_age_at_birth>= 17&& educ_mother== EOL_MEDIUM )? CMG_01: ( mother_age_at_birth>= 17&& educ_mother== EOL_LOW )? CMG_02: ( mother_age_at_birth>= 15&& educ_mother== EOL_HIGH )? CMG_03: ( mother_age_at_birth>= 15&& educ_mother== EOL_MEDIUM )? CMG_04: ( mother_age_at_birth>= 15&& educ_mother== EOL_LOW )? CMG_05: ( educ_mother== EOL_HIGH )? CMG_06: ( educ_mother== EOL_MEDIUM )? CMG_07: CMG_08 ) );
#line 3974 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 75 "../code/FamilyGeneral.mpp"
	children_in_household = short( ( sex== FEMALE )? sum_over_mlHHMotherChildren_is_alive: sum_over_mlHHFatherChildren_is_alive );
#line 3977 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 100 "../code/EducationPrimaryTransmission.mpp"
	poStateFunctions->Set_educ1_group( EDUC1_GROUP( ( educ_mother== EOL_LOW )? E1G_00: ( educ_mother== EOL_MEDIUM )? E1G_01: E1G_02 ) );
#line 3980 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 107 "../code/EducationPrimaryTracking.mpp"
	poStateFunctions->Set_educ_one_grade( EDUC_ONE_GRADE( COERCE(EDUC_ONE_GRADE,educ_one_grade_attended) ) );
#line 3983 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	educ_one_grade = educ_one_grade < r_min_EDUC_ONE_GRADE ? r_min_EDUC_ONE_GRADE : educ_one_grade > r_max_EDUC_ONE_GRADE ? r_max_EDUC_ONE_GRADE : educ_one_grade;
#line 112 "../code/EducationPrimaryTracking.mpp"
	poStateFunctions->Set_educ_one_level( EDUC_ONE_LEVEL( ( educ_one_grade_passed== MAX(EDUC_ONE_GRADE) )? EOL_HIGH: ( educ_one_grade_attended== 0 )? EOL_LOW: EOL_MEDIUM ) );
#line 3987 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 43 "../code/EducationPrimaryPlanning.mpp"
	poStateFunctions->Set_educ_one_required_rooms( double( ( is_educ_one_student&& is_resident )? 1.0/ gprParam->Educ1Infrastructure[CHKBND(SI_ROOMS, 2, "Educ1Infrastructure", 0)][CHKBND(geo_nat, 12, "Educ1Infrastructure", 1)][CHKBND(RANGE_POS(SIM_YEAR_RANGE,calendar_year), 151, "Educ1Infrastructure", 2)]: 0.0 ) );
#line 3990 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 135 "../code/TablesEducation.mpp"
	poStateFunctions->Set_educ_one_required_rooms_value_out( double( educ_one_required_rooms ) );
#line 3993 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 39 "../code/EducationPrimaryPlanning.mpp"
	poStateFunctions->Set_educ_one_required_teachers( double( ( is_educ_one_student&& is_resident )? 1.0/ gprParam->Educ1Infrastructure[CHKBND(SI_TEACHERS, 2, "Educ1Infrastructure", 0)][CHKBND(geo_nat, 12, "Educ1Infrastructure", 1)][CHKBND(RANGE_POS(SIM_YEAR_RANGE,calendar_year), 151, "Educ1Infrastructure", 2)]: 0.0 ) );
#line 3996 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 135 "../code/TablesEducation.mpp"
	poStateFunctions->Set_educ_one_required_teachers_value_out( double( educ_one_required_teachers ) );
#line 3999 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 113 "../code/EducationSecondaryBase.mpp"
	poStateFunctions->Set_educ_two_level( EDUC_TWO_LEVEL( ( educ_two_status== ETS_NEVER )? ETL_NEVER: ( educ_two_grade_passed== MAX(EDUC_TWO_GRADE) )? ETL_INCOMPLETE: ETL_GRADUATE ) );
#line 4002 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 52 "../code/ChildVaccination.mpp"
	poStateFunctions->Set_ethnicity_short( ETHNICITY_SHORT( poDerivedStates->aggregate_ethnicity_ETHNICITY_SHORT ) );
#line 4005 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 71 "../code/FamilyGeneral.mpp"
	poStateFunctions->Set_family_role( FAM_ROLE( ( lSpouse&& sex== MALE )? FR_SPOUSE: ( lHHFather|| lHHMother )? FR_CHILD: FR_HEAD ) );
#line 4008 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 81 "../code/FamilyPartnerMatching.mpp"
	poStateFunctions->Set_female_spouse_age( FEMALE_SPOUSE_AGE( ( sex== FEMALE )? COERCE(FEMALE_SPOUSE_AGE,integer_age): MAX(FEMALE_SPOUSE_AGE) ) );
#line 4011 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	female_spouse_age = female_spouse_age < r_min_FEMALE_SPOUSE_AGE ? r_min_FEMALE_SPOUSE_AGE : female_spouse_age > r_max_FEMALE_SPOUSE_AGE ? r_max_FEMALE_SPOUSE_AGE : female_spouse_age;
#line 49 "../code/FertilityGeneral.mpp"
	poStateFunctions->Set_fertile_age( FERTILE_AGE_RANGE( COERCE(FERTILE_AGE_RANGE,integer_age) ) );
#line 4015 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	fertile_age = fertile_age < r_min_FERTILE_AGE_RANGE ? r_min_FERTILE_AGE_RANGE : fertile_age > r_max_FERTILE_AGE_RANGE ? r_max_FERTILE_AGE_RANGE : fertile_age;
#line 55 "../code/MigrationGeneral.mpp"
	poStateFunctions->Set_geo_nat( GEO_NAT( poDerivedStates->aggregate_geo_GEO_NAT ) );
#line 4019 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 25 "../code/TablesEducation.mpp"
	poStateFunctions->Set_has_passed_primary( logical( ( educ_one_grade_passed== MAX(EDUC_ONE_GRADE) ) ) );
#line 4022 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 83 "../code/FamilyPartnerMatching.mpp"
	poStateFunctions->Set_has_spouse( logical( ( lSpouse )? TRUE: FALSE ) );
#line 4025 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 82 "../code/HumanCapitalIndex.mpp"
	poStateFunctions->Set_in_hci_sample( logical( ( creation_type== CT_BIRTH&& geo_birth< SIZE(GEO_NAT) )? TRUE: FALSE ) );
#line 4028 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 18 "../code/TablesEducation.mpp"
	poStateFunctions->Set_in_projected_time( logical( ( calendar_year>= MIN(SIM_YEAR_RANGE) ) ) );
#line 4031 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 36 "../code/_ContextABC.mpp"
	poStateFunctions->Set_in_sel_years( logical( calendar_year== 2000|| calendar_year== 2010|| calendar_year== 2020|| calendar_year== 2030|| calendar_year== 2040|| calendar_year== 2050|| calendar_year== 2060 ) );
#line 4034 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 34 "../code/TablesPopulation.mpp"
	poStateFunctions->Set_in_selected_year( logical( ( selected_year10!= SY40|| calendar_year== MIN(SIM_YEAR_RANGE)+ 40 ) ) );
#line 4037 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 22 "../code/TablesEducation.mpp"
	poStateFunctions->Set_is_aded0( logical( ( !( is_child )&& !( is_old )&& educ_one_fate== EOL_LOW ) ) );
#line 4040 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 23 "../code/TablesEducation.mpp"
	poStateFunctions->Set_is_aded1( logical( ( !( is_child )&& !( is_old )&& educ_one_fate== EOL_MEDIUM ) ) );
#line 4043 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 24 "../code/TablesEducation.mpp"
	poStateFunctions->Set_is_aded2( logical( ( !( is_child )&& !( is_old )&& educ_one_fate== EOL_HIGH ) ) );
#line 4046 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 103 "../code/HumanCapitalIndex.mpp"
	is_adult = logical( ( integer_age>= 15 )? TRUE: FALSE );
#line 4049 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 73 "../code/TablesPopulation.mpp"
	poStateFunctions->Set_is_alive_value_in( logical( is_alive ) );
#line 4052 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 20 "../code/TablesEducation.mpp"
	poStateFunctions->Set_is_child( logical( ( integer_age< 15 ) ) );
#line 4055 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 104 "../code/EducationPrimaryTracking.mpp"
	poStateFunctions->Set_is_educ_one_student( logical( ( educ_one_status== EOS_ATTEND|| ( educ_one_status== EOS_WAIT&& educ_one_grade_attended> 0 ) )? TRUE: FALSE ) );
#line 4058 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 47 "../code/FertilityGeneral.mpp"
	poStateFunctions->Set_is_fertile( logical( ( sex== FEMALE&& WITHIN(FERTILE_AGE_RANGE,integer_age)&& ( creation_type== CT_BIRTH|| ( creation_type== CT_START&& calendar_year>= MIN(SIM_YEAR_RANGE) )|| ( creation_type== CT_POOL&& ever_resident )|| ( creation_type== CT_SCRATCH&& ever_resident ) ) )? TRUE: FALSE ) );
#line 4061 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 45 "../code/MortalityGeneral.mpp"
	poStateFunctions->Set_is_mortal( logical( ( creation_type== CT_BIRTH|| ( creation_type== CT_START&& calendar_year>= MIN(SIM_YEAR_RANGE) )|| ( creation_type== CT_POOL&& ever_resident )|| ( creation_type== CT_SCRATCH&& ever_resident ) )? TRUE: FALSE ) );
#line 4064 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 21 "../code/TablesEducation.mpp"
	poStateFunctions->Set_is_old( logical( ( integer_age>= 60 ) ) );
#line 4067 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 79 "../code/ActorPerson.mpp"
	poStateFunctions->Set_is_resident( logical( ( geo< SIZE(GEO_NAT) ) ) );
#line 4070 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 23 "../code/TablesStunting.mpp"
	poStateFunctions->Set_is_stunted_value_in( logical( is_stunted ) );
#line 4073 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 92 "../code/FamilyFemalePartnershipStatus.mpp"
	poStateFunctions->Set_lives_with_dependent_child( logical( ( ( sex== FEMALE&& count_mlHHMotherChildren> 0 )|| ( sex== MALE&& count_mlHHFatherChildren> 0 ) )? TRUE: FALSE ) );
#line 4076 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 78 "../code/FamilyPartnerMatching.mpp"
	poStateFunctions->Set_male_age( MALE_SPOUSE_AGE( ( sex== MALE )? COERCE(MALE_SPOUSE_AGE,integer_age): MAX(MALE_SPOUSE_AGE) ) );
#line 4079 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	male_age = male_age < r_min_MALE_SPOUSE_AGE ? r_min_MALE_SPOUSE_AGE : male_age > r_max_MALE_SPOUSE_AGE ? r_max_MALE_SPOUSE_AGE : male_age;
#line 75 "../code/FamilyPartnerMatching.mpp"
	poStateFunctions->Set_male_spouse_age( MALE_SPOUSE_AGE( ( has_spouse&& sex== FEMALE )? COERCE(MALE_SPOUSE_AGE,lSpouse_integer_age): MAX(MALE_SPOUSE_AGE) ) );
#line 4083 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	male_spouse_age = male_spouse_age < r_min_MALE_SPOUSE_AGE ? r_min_MALE_SPOUSE_AGE : male_spouse_age > r_max_MALE_SPOUSE_AGE ? r_max_MALE_SPOUSE_AGE : male_spouse_age;
#line 106 "../code/FamilyFemalePartnershipStatus.mpp"
	poStateFunctions->Set_moth_agegr( MOTH_AGEGR( ( moth_agegr_part== 0 )? CMA20: ( moth_agegr_part== 1 )? CMA25: ( moth_agegr_part== 2 )? CMA30: ( moth_agegr_part== 3 )? CMA35: ( moth_agegr_part== 4 )? CMA40: CMA40P ) );
#line 4087 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 99 "../code/FamilyFemalePartnershipStatus.mpp"
	poStateFunctions->Set_moth_agegr_part( int( poDerivedStates->split_age_last_birth_MOTH_AGEGR_PART ) );
#line 4090 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 59 "../code/MigrationGeneral.mpp"
	poStateFunctions->Set_region_birth( REGION_INT( poDerivedStates->aggregate_geo_birth_REGION_INT ) );
#line 4093 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 57 "../code/MigrationGeneral.mpp"
	poStateFunctions->Set_region_int( REGION_INT( poDerivedStates->aggregate_geo_REGION_INT ) );
#line 4096 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 58 "../code/MigrationGeneral.mpp"
	poStateFunctions->Set_region_nat( REGION_NAT( poDerivedStates->aggregate_region_int_REGION_NAT ) );
#line 4099 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 44 "../code/_ContextABC.mpp"
	poStateFunctions->Set_sel_years( TAB_SEL_YEARS( ( calendar_year== 2000 )? TSY_2000: ( calendar_year== 2010 )? TSY_2010: ( calendar_year== 2020 )? TSY_2020: ( calendar_year== 2030 )? TSY_2030: ( calendar_year== 2040 )? TSY_2040: ( calendar_year== 2050 )? TSY_2050: TSY_2060 ) );
#line 4102 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 31 "../code/TablesPopulation.mpp"
	poStateFunctions->Set_selected_year10( SELECTED_YEAR10( ( calendar_year== MIN(SIM_YEAR_RANGE) )? SY00: ( calendar_year== MIN(SIM_YEAR_RANGE)+ 10 )? SY10: ( calendar_year== MIN(SIM_YEAR_RANGE)+ 20 )? SY20: ( calendar_year== MIN(SIM_YEAR_RANGE)+ 30 )? SY30: SY40 ) );
#line 4105 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 24 "../code/TablesPopulation.mpp"
	poStateFunctions->Set_sim_year( SIM_YEAR_RANGE( COERCE(SIM_YEAR_RANGE,calendar_year) ) );
#line 4108 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	sim_year = sim_year < r_min_SIM_YEAR_RANGE ? r_min_SIM_YEAR_RANGE : sim_year > r_max_SIM_YEAR_RANGE ? r_max_SIM_YEAR_RANGE : sim_year;
	poDerivedStates->Set_split_age_last_birth_MOTH_AGEGR_PART( Split( (double) age_last_birth, poDerivedStates->split_age_last_birth_MOTH_AGEGR_PART, 5, p_points_MOTH_AGEGR_PART ) );
	poDerivedStates->Set_split_age_youngest_child_CHILD_AGEGR_PART( Split( (double) age_youngest_child, poDerivedStates->split_age_youngest_child_CHILD_AGEGR_PART, 6, p_points_CHILD_AGEGR_PART ) );
	poDerivedStates->Set_split_integer_age_TAB_AGEC5( Split( (double) integer_age, poDerivedStates->split_integer_age_TAB_AGEC5, 20, p_points_TAB_AGEC5 ) );
#line 87 "../code/FamilyPartnerMatching.mpp"
	poStateFunctions->Set_spouse_group( SPOUSE_GROUP( ( has_spouse&& lSpouse_educ_one_fate== EOL_HIGH )? SMG_02: ( has_spouse&& lSpouse_educ_one_fate== EOL_MEDIUM )? SMG_01: SMG_00 ) );
#line 4115 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 97 "../code/HumanCapitalIndex.mpp"
	poStateFunctions->Set_survived_early_years( logical( ( integer_age>= 5 )? TRUE: FALSE ) );
#line 4118 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/TablesEducation.mpp"
	poStateFunctions->Set_tab_out_of_primary( logical( ( educ_one_status!= EOS_WAIT&& educ_one_status!= EOS_ATTEND&& !( has_passed_primary ) ) ) );
#line 4121 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 19 "../code/TablesEducation.mpp"
	poStateFunctions->Set_tab_sim_yob( SIM_YEAR_RANGE( COERCE(SIM_YEAR_RANGE,year_of_birth) ) );
#line 4124 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	tab_sim_yob = tab_sim_yob < r_min_SIM_YEAR_RANGE ? r_min_SIM_YEAR_RANGE : tab_sim_yob > r_max_SIM_YEAR_RANGE ? r_max_SIM_YEAR_RANGE : tab_sim_yob;
#line 103 "../code/FertilityDetailed.mpp"
	poStateFunctions->Set_time_in_parity( int( poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY ) );
#line 4128 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 88 "../code/FamilyFirstUnion.mpp"
	poStateFunctions->Set_union1_group( UNION1_GROUP( ( educ_one_fate== EOL_LOW )? U1G_00: ( educ_one_fate== EOL_MEDIUM )? U1G_01: U1G_02 ) );
#line 4131 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 88 "../code/ActorPerson.mpp"
	poStateFunctions->Set_year_of_birth( ALL_YEAR_RANGE( COERCE(ALL_YEAR_RANGE,int(time_of_birth)) ) );
#line 4134 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	year_of_birth = year_of_birth < r_min_ALL_YEAR_RANGE ? r_min_ALL_YEAR_RANGE : year_of_birth > r_max_ALL_YEAR_RANGE ? r_max_ALL_YEAR_RANGE : year_of_birth;
#line 52 "../code/ImmigrationFromScratch.mpp"
	poStateFunctions->Set_year_of_first_immigration( SIM_YEAR_RANGE( COERCE(SIM_YEAR_RANGE,int(time_of_first_immigration)) ) );
#line 4138 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	year_of_first_immigration = year_of_first_immigration < r_min_SIM_YEAR_RANGE ? r_min_SIM_YEAR_RANGE : year_of_first_immigration > r_max_SIM_YEAR_RANGE ? r_max_SIM_YEAR_RANGE : year_of_first_immigration;
#line 85 "../code/HumanCapitalIndex.mpp"
	years_of_schooling = int( years_preschool+ educ_one_grade_passed+ educ_two_grade_passed );
#line 4142 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	CountActor(3);
}

#pragma optimize( "", off )

void Person::StartSpecial2()
{
	// initialization of expression and history states
	poStateFunctions->Set_time(TIME( CoarsenMantissa( time ) ), true);
	poStateFunctions->Set_age(TIME( CoarsenMantissa( age ) ), true);
	gdThreadCurrentTime = time;
	__time = time;
	poDerivedStates->active_spell_duration_this_parity_spell_TRUE = (TIME) 0;
	poDerivedStates->active_spell_duration_this_parity_spell_TRUE_time = time;
#line 31 "../code/FertilityDetailed.mpp"
	poStateFunctions->Set_actorset_filter_asAllFertilePersonsForFertilityAlignment( logical( gprParam->FertilityModel== FEM_ALIGNED_TOTAL&& is_fertile ) );
#line 4159 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 27 "../code/FertilityDetailed.mpp"
	poStateFunctions->Set_actorset_filter_asAllFertilePersonsForFertilityAlignmentByAge( logical( gprParam->FertilityModel== FEM_ALIGNED_AGE&& is_fertile ) );
#line 4162 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 40 "../code/ActorPerson.mpp"
	poStateFunctions->Set_actorset_filter_asAllPerson( logical( is_alive ) );
#line 4165 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 14 "../code/FamilyPartnerMatching.mpp"
	poStateFunctions->Set_actorset_filter_asAvailableMale( logical( is_alive&& sex== MALE&& !( has_spouse )&& WITHIN(MALE_SPOUSE_AGE,integer_age)&& is_mortal ) );
#line 4168 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 18 "../code/FamilyPartnerMatching.mpp"
	poStateFunctions->Set_actorset_filter_asAvailableMaleByType( logical( is_alive&& sex== MALE&& !( has_spouse )&& WITHIN(MALE_SPOUSE_AGE,integer_age)&& is_mortal ) );
#line 4171 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 32 "../code/FamilyPartnerMatching.mpp"
	poStateFunctions->Set_actorset_filter_asAvailableMaleByTypeImmiScratch( logical( is_alive&& creation_type== CT_SCRATCH&& sex== MALE&& !( has_spouse )&& WITHIN(MALE_SPOUSE_AGE,integer_age)&& !( is_mortal ) ) );
#line 4174 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/FamilyPartnerMatching.mpp"
	poStateFunctions->Set_actorset_filter_asAvailableMaleImmiScratch( logical( is_alive&& creation_type== CT_SCRATCH&& sex== MALE&& !( has_spouse )&& WITHIN(MALE_SPOUSE_AGE,integer_age)&& !( is_mortal ) ) );
#line 4177 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 22 "../code/FamilyPartnerMatching.mpp"
	poStateFunctions->Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAge( logical( is_alive&& sex== FEMALE&& has_spouse&& WITHIN(FEMALE_SPOUSE_AGE,integer_age)&& is_mortal ) );
#line 4180 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 36 "../code/FamilyPartnerMatching.mpp"
	poStateFunctions->Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch( logical( is_alive&& creation_type== CT_SCRATCH&& sex== FEMALE&& has_spouse&& WITHIN(FEMALE_SPOUSE_AGE,integer_age)&& !( is_mortal ) ) );
#line 4183 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 10 "../code/ImmigrationFromScratch.mpp"
	poStateFunctions->Set_actorset_filter_asPotentialImmigrantMothers( logical( is_alive&& sex== FEMALE&& creation_type== CT_SCRATCH&& !( ever_resident )&& family_role!= FR_CHILD&& ever_union ) );
#line 4186 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 25 "../code/EducationPrimaryTracking.mpp"
	poStateFunctions->Set_actorset_filter_asPotentialSchoolOneStudent( logical( is_alive&& educ_one_status!= EOS_OUT&& educ_one_fate!= EOL_LOW ) );
#line 4189 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/EducationPrimaryTracking.mpp"
	poStateFunctions->Set_actorset_filter_asPotentialSchoolOneStudentToProcess( logical( is_alive&& educ_one_to_process ) );
#line 4192 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 26 "../code/EducationSecondaryBase.mpp"
	poStateFunctions->Set_actorset_filter_asPotentialSchoolTwoStudent( logical( educ_one_grade_passed== MAX(EDUC_ONE_GRADE)&& educ_two_status!= ETS_OUT&& educ_two_delay<= gprParam->Educ2AllowedDelays ) );
#line 4195 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 29 "../code/EducationSecondaryBase.mpp"
	poStateFunctions->Set_actorset_filter_asPotentialSchoolTwoStudentToProcess( logical( is_alive&& educ_two_to_process ) );
#line 4198 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 18 "../code/ChildVaccination.mpp"
	poStateFunctions->Set_actorset_filter_asResidentBabies( logical( is_alive&& integer_age== 0&& is_resident&& creation_type== CT_BIRTH ) );
#line 4201 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 30 "../code/EducationPrimaryBase.mpp"
	poStateFunctions->Set_actorset_filter_asResidentsAge0Sex( logical( is_alive&& integer_age== 0&& is_resident ) );
#line 4204 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 27 "../code/EducationPrimaryBase.mpp"
	poStateFunctions->Set_actorset_filter_asResidentsAge0SexGeo( logical( is_alive&& is_resident&& integer_age== 0 ) );
#line 4207 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 35 "../code/EducationPrimaryTransmission.mpp"
	poStateFunctions->Set_actorset_filter_asSimBornAge0( logical( is_alive&& creation_type== CT_BIRTH&& integer_age== 0 ) );
#line 4210 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 11 "../code/MigrationGeneral.mpp"
	poStateFunctions->Set_actorset_filter_asWantToMove( logical( is_alive ) );
#line 4213 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/FamilyFemalePartnershipStatus.mpp"
	poStateFunctions->Set_actorset_filter_asWomenNoChildren( logical( is_alive&& sex== FEMALE&& WITHIN(FEMALE_SPOUSE_AGE,integer_age)&& !( lives_with_dependent_child )&& is_mortal&& !( is_blocked_from_marriage ) ) );
#line 4216 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 23 "../code/FamilyFemalePartnershipStatus.mpp"
	poStateFunctions->Set_actorset_filter_asWomenWithChildren( logical( is_alive&& sex== FEMALE&& WITHIN(FEMALE_SPOUSE_AGE,integer_age)&& lives_with_dependent_child&& is_mortal&& !( is_blocked_from_marriage ) ) );
#line 4219 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	poDerivedStates->completed_spell_duration_this_parity_spell_TRUE = (TIME) 0;
	poDerivedStates->Set_duration_( (TIME) 0 );
	poDerivedStates->duration__time = time;
	poDerivedStates->Set_duration_educ_one_fate_EOL_HIGH( (TIME) 0 );
	poDerivedStates->duration_educ_one_fate_EOL_HIGH_time = time;
	poDerivedStates->Set_duration_educ_one_fate_EOL_LOW( (TIME) 0 );
	poDerivedStates->duration_educ_one_fate_EOL_LOW_time = time;
	poDerivedStates->Set_duration_educ_one_fate_EOL_MEDIUM( (TIME) 0 );
	poDerivedStates->duration_educ_one_fate_EOL_MEDIUM_time = time;
	poDerivedStates->Set_duration_educ_one_level_EOL_HIGH( (TIME) 0 );
	poDerivedStates->duration_educ_one_level_EOL_HIGH_time = time;
	poDerivedStates->Set_duration_educ_one_level_EOL_LOW( (TIME) 0 );
	poDerivedStates->duration_educ_one_level_EOL_LOW_time = time;
	poDerivedStates->Set_duration_educ_one_level_EOL_MEDIUM( (TIME) 0 );
	poDerivedStates->duration_educ_one_level_EOL_MEDIUM_time = time;
	poDerivedStates->Set_duration_educ_two_level_ETL_GRADUATE( (TIME) 0 );
	poDerivedStates->duration_educ_two_level_ETL_GRADUATE_time = time;
	poDerivedStates->Set_duration_educ_two_level_ETL_INCOMPLETE( (TIME) 0 );
	poDerivedStates->duration_educ_two_level_ETL_INCOMPLETE_time = time;
	poDerivedStates->Set_duration_educ_two_level_ETL_NEVER( (TIME) 0 );
	poDerivedStates->duration_educ_two_level_ETL_NEVER_time = time;
	poDerivedStates->Set_duration_ever_union_TRUE( (TIME) 0 );
	poDerivedStates->duration_ever_union_TRUE_time = time;
	poDerivedStates->Set_duration_got_prenat_care_GPC_YES( (TIME) 0 );
	poDerivedStates->duration_got_prenat_care_GPC_YES_time = time;
	poDerivedStates->Set_duration_has_spouse_TRUE( (TIME) 0 );
	poDerivedStates->duration_has_spouse_TRUE_time = time;
	poDerivedStates->Set_duration_is_aded0_TRUE( (TIME) 0 );
	poDerivedStates->duration_is_aded0_TRUE_time = time;
	poDerivedStates->Set_duration_is_aded1_TRUE( (TIME) 0 );
	poDerivedStates->duration_is_aded1_TRUE_time = time;
	poDerivedStates->Set_duration_is_aded2_TRUE( (TIME) 0 );
	poDerivedStates->duration_is_aded2_TRUE_time = time;
	poDerivedStates->Set_duration_is_child_TRUE( (TIME) 0 );
	poDerivedStates->duration_is_child_TRUE_time = time;
	poDerivedStates->Set_duration_is_immunized_TRUE( (TIME) 0 );
	poDerivedStates->duration_is_immunized_TRUE_time = time;
	poDerivedStates->Set_duration_is_old_TRUE( (TIME) 0 );
	poDerivedStates->duration_is_old_TRUE_time = time;
	poDerivedStates->Set_duration_tab_out_of_primary_TRUE( (TIME) 0 );
	poDerivedStates->duration_tab_out_of_primary_TRUE_time = time;
	poDerivedStates->Set_duration_years_preschool_0( (TIME) 0 );
	poDerivedStates->duration_years_preschool_0_time = time;
	poDerivedStates->Set_duration_years_preschool_1( (TIME) 0 );
	poDerivedStates->duration_years_preschool_1_time = time;
	poDerivedStates->Set_duration_years_preschool_2( (TIME) 0 );
	poDerivedStates->duration_years_preschool_2_time = time;
	poDerivedStates->Set_entrances_educ_one_grade_attended_1( (counter) 0 );
	poDerivedStates->Set_entrances_has_passed_primary_TRUE( (counter) 0 );
	poDerivedStates->Set_sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY( 0 );
	if ( this_parity_spell == TRUE )
	{
		poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time = WAIT( POINTS(DUR_TIME_IN_PARITY)[0] );
		if (schedule_event_time > poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time)
		{
			schedule_event_time = poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time;
		}
	}
	else
	{
		poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time = TIME_INFINITE;
	}
#line 101 "../code/TablesEducation.mpp"
	poStateFunctions->Set_table_filter_PopPyramidByEduc( logical( in_sel_years&& is_resident ) );
#line 4284 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 55 "../code/TablesPopulation.mpp"
	poStateFunctions->Set_table_filter_tabBirthsYearPlace( logical( WITHIN(SIM_YEAR_RANGE,calendar_year) ) );
#line 4287 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 46 "../code/TablesVaccination.mpp"
	poStateFunctions->Set_table_filter_TabChildVaccination( logical( in_projected_time&& creation_type== CT_BIRTH&& is_resident&& integer_age== 0 ) );
#line 4290 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 165 "../code/TablesEducation.mpp"
	poStateFunctions->Set_table_filter_TabEduc15ByDistrict( logical( integer_age== 15&& in_projected_time&& is_resident ) );
#line 4293 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 178 "../code/TablesEducation.mpp"
	poStateFunctions->Set_table_filter_TabEduc15ByDistrictBirth( logical( integer_age== 15&& in_projected_time&& is_resident ) );
#line 4296 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 88 "../code/TablesEducation.mpp"
	poStateFunctions->Set_table_filter_tabEducationFateGeobirYob( logical( integer_age== 1 ) );
#line 4299 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 115 "../code/TablesEducation.mpp"
	poStateFunctions->Set_table_filter_TabEducFateByGroup( logical( in_projected_time&& poDerivedStates->trigger_changes_calendar_year&& integer_age== 0&& creation_type== CT_BIRTH ) );
#line 4302 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 154 "../code/TablesEducation.mpp"
	poStateFunctions->Set_table_filter_TabEducFateDistrYob( logical( poDerivedStates->trigger_changes_calendar_year&& integer_age== 0&& is_resident ) );
#line 4305 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 13 "../code/TablesFamily.mpp"
	poStateFunctions->Set_table_filter_tabHavingSpouse( logical( WITHIN(SIM_YEAR_RANGE,calendar_year)&& is_resident ) );
#line 4308 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 17 "../code/TablesHumanCapital.mpp"
	poStateFunctions->Set_table_filter_tabHCI( logical( in_projected_time&& in_hci_sample&& year_of_birth>= MIN(SIM_YEAR_RANGE) ) );
#line 4311 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 45 "../code/TablesHumanCapital.mpp"
	poStateFunctions->Set_table_filter_TabHCIDistrict( logical( in_projected_time&& in_hci_sample&& year_of_birth>= MIN(SIM_YEAR_RANGE) ) );
#line 4314 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 31 "../code/TablesVaccination.mpp"
	poStateFunctions->Set_table_filter_TabImmunizationChildren( logical( in_projected_time&& integer_age== 0&& creation_type== CT_BIRTH&& is_resident ) );
#line 4317 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 66 "../code/TablesPopulation.mpp"
	poStateFunctions->Set_table_filter_tabMigrationOriginDestination( logical( WITHIN(SIM_YEAR_RANGE,calendar_year)&& poDerivedStates->trigger_changes_geo&& in_selected_year ) );
#line 4320 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 139 "../code/TablesEducation.mpp"
	poStateFunctions->Set_table_filter_TabPopProvAgeEducSex( logical( in_projected_time&& is_resident ) );
#line 4323 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 45 "../code/TablesPopulation.mpp"
	poStateFunctions->Set_table_filter_tabPopulationYearPlace( logical( WITHIN(SIM_YEAR_RANGE,calendar_year) ) );
#line 4326 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 19 "../code/TablesVaccination.mpp"
	poStateFunctions->Set_table_filter_TabPrenatCare( logical( in_projected_time&& integer_age== 0&& creation_type== CT_BIRTH&& is_resident ) );
#line 4329 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 76 "../code/TablesEducation.mpp"
	poStateFunctions->Set_table_filter_tabPreSchool( logical( WITHIN(SIM_YEAR_RANGE,year_of_birth)&& integer_age== 10 ) );
#line 4332 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 126 "../code/TablesEducation.mpp"
	poStateFunctions->Set_table_filter_TabPrimarySchoolPlanning( logical( in_projected_time&& poDerivedStates->trigger_changes_calendar_year&& is_educ_one_student&& is_resident ) );
#line 4335 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 191 "../code/TablesEducation.mpp"
	poStateFunctions->Set_table_filter_TabPrimSchoolEntries( logical( in_projected_time&& is_resident ) );
#line 4338 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 201 "../code/TablesEducation.mpp"
	poStateFunctions->Set_table_filter_TabPrimSchoolGraduations( logical( in_projected_time&& is_resident ) );
#line 4341 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 211 "../code/TablesEducation.mpp"
	poStateFunctions->Set_table_filter_TabPrimSchoolOutOfSchool9to11( logical( in_projected_time&& is_resident&& integer_age>= 9&& integer_age<= 11 ) );
#line 4344 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 231 "../code/TablesEducation.mpp"
	poStateFunctions->Set_table_filter_TabSchool2AttainmentsTab( logical( in_projected_time&& is_resident&& integer_age>= 20&& integer_age< 25 ) );
#line 4347 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 221 "../code/TablesEducation.mpp"
	poStateFunctions->Set_table_filter_TabSchool2TrackTab( logical( poDerivedStates->trigger_changes_calendar_year&& in_projected_time&& is_resident&& educ_two_status== ETS_ATTEND ) );
#line 4350 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 27 "../code/TablesStunting.mpp"
	poStateFunctions->Set_table_filter_tabStuntingSexRegMotherYob( logical( in_projected_time&& is_resident&& creation_type== CT_BIRTH&& poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE ) );
#line 4353 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 17 "../code/TablesStunting.mpp"
	poStateFunctions->Set_table_filter_tabStuntingSexRegYob( logical( in_projected_time&& is_resident&& creation_type== CT_BIRTH&& poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE ) );
#line 4356 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	poDerivedStates->Set_transitions_integer_age_14_15( (counter) 0 );
	poDerivedStates->Set_transitions_is_alive_FALSE_TRUE( (counter) 0 );
	poDerivedStates->Set_transitions_is_alive_TRUE_FALSE( (counter) 0 );
	poDerivedStates->Set_transitions_survived_early_years_FALSE_TRUE( (counter) 0 );
	poDerivedStates->Set_trigger_changes_calendar_year( FALSE );
	poDerivedStates->trigger_changes_calendar_year_time = TIME_INFINITE;
	poDerivedStates->Set_trigger_changes_geo( FALSE );
	poDerivedStates->trigger_changes_geo_time = TIME_INFINITE;
	poDerivedStates->Set_trigger_entrances_is_ready_for_birthtables_TRUE( FALSE );
	poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE_time = TIME_INFINITE;
	poDerivedStates->Set_value_at_transitions_is_alive_TRUE_FALSE_adult_survival( (real) 0 );
	poDerivedStates->Set_value_at_transitions_is_alive_TRUE_FALSE_ind_hci( (real) 0 );
	poDerivedStates->Set_value_at_transitions_is_alive_TRUE_FALSE_is_stunted( (integer) 0 );
	poDerivedStates->Set_value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling( (real) 0 );
	poDerivedStates->Set_value_at_transitions_is_alive_TRUE_FALSE_survived_early_years( (integer) 0 );
	poDerivedStates->Set_value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling( (integer) 0 );
	UpdateLinkedStates();

	// initialization of cross-tabulation updates
	tableinfo.tab0events = -1;
	tableinfo.tab0open0 = (real) ( poDerivedStates->duration_ );
	tableinfo.tab0open1 = (real) ( poDerivedStates->duration_years_preschool_0 );
	tableinfo.tab0open2 = (real) ( poDerivedStates->duration_years_preschool_1 );
	tableinfo.tab0open3 = (real) ( poDerivedStates->duration_years_preschool_2 );
	tableinfo.tab1events = -1;
	tableinfo.tab1open0 = (real) ( poDerivedStates->duration_ );
	tableinfo.tab1open1 = (real) ( poDerivedStates->duration_educ_one_fate_EOL_HIGH );
	tableinfo.tab1open2 = (real) ( poDerivedStates->duration_educ_one_fate_EOL_LOW );
	tableinfo.tab1open3 = (real) ( poDerivedStates->duration_educ_one_fate_EOL_MEDIUM );
	tableinfo.tab2events = -1;
	tableinfo.tab2open0 = (real) ( poDerivedStates->duration_educ_one_level_EOL_HIGH );
	tableinfo.tab2open1 = (real) ( poDerivedStates->duration_educ_one_level_EOL_LOW );
	tableinfo.tab2open2 = (real) ( poDerivedStates->duration_educ_one_level_EOL_MEDIUM );
	tableinfo.tab3events = -1;
	tableinfo.tab3open0 = (real) ( table_unit_TabEducFateByGroup );
	tableinfo.tab4events = -1;
	tableinfo.tab4open0 = (real) ( educ_one_required_rooms_value_out );
	tableinfo.tab4open1 = (real) ( educ_one_required_teachers_value_out );
	tableinfo.tab4open2 = (real) ( table_unit_TabPrimarySchoolPlanning );
	tableinfo.tab5events = -1;
	tableinfo.tab5open0 = (real) ( poDerivedStates->duration_is_aded0_TRUE );
	tableinfo.tab5open1 = (real) ( poDerivedStates->duration_is_aded1_TRUE );
	tableinfo.tab5open2 = (real) ( poDerivedStates->duration_is_aded2_TRUE );
	tableinfo.tab5open3 = (real) ( poDerivedStates->duration_is_child_TRUE );
	tableinfo.tab5open4 = (real) ( poDerivedStates->duration_is_old_TRUE );
	tableinfo.tab6events = -1;
	tableinfo.tab6open0 = (real) ( table_unit_TabEducFateDistrYob );
	tableinfo.tab7events = -1;
	tableinfo.tab7open0 = (real) ( poDerivedStates->duration_ );
	tableinfo.tab7open1 = (real) ( poDerivedStates->duration_educ_one_fate_EOL_HIGH );
	tableinfo.tab7open2 = (real) ( poDerivedStates->duration_educ_one_fate_EOL_LOW );
	tableinfo.tab7open3 = (real) ( poDerivedStates->duration_educ_one_fate_EOL_MEDIUM );
	tableinfo.tab8events = -1;
	tableinfo.tab8open0 = (real) ( poDerivedStates->duration_ );
	tableinfo.tab8open1 = (real) ( poDerivedStates->duration_educ_one_fate_EOL_HIGH );
	tableinfo.tab8open2 = (real) ( poDerivedStates->duration_educ_one_fate_EOL_LOW );
	tableinfo.tab8open3 = (real) ( poDerivedStates->duration_educ_one_fate_EOL_MEDIUM );
	tableinfo.tab9events = -1;
	tableinfo.tab9open0 = (real) ( poDerivedStates->entrances_educ_one_grade_attended_1 );
	tableinfo.tab10events = -1;
	tableinfo.tab10open0 = (real) ( poDerivedStates->entrances_has_passed_primary_TRUE );
	tableinfo.tab11events = -1;
	tableinfo.tab11open0 = (real) ( poDerivedStates->duration_ );
	tableinfo.tab11open1 = (real) ( poDerivedStates->duration_tab_out_of_primary_TRUE );
	tableinfo.tab12events = -1;
	tableinfo.tab12open0 = (real) ( table_unit_TabSchool2TrackTab );
	tableinfo.tab13events = -1;
	tableinfo.tab13open0 = (real) ( poDerivedStates->duration_ );
	tableinfo.tab13open1 = (real) ( poDerivedStates->duration_educ_two_level_ETL_GRADUATE );
	tableinfo.tab13open2 = (real) ( poDerivedStates->duration_educ_two_level_ETL_INCOMPLETE );
	tableinfo.tab13open3 = (real) ( poDerivedStates->duration_educ_two_level_ETL_NEVER );
	tableinfo.tab14events = -1;
	tableinfo.tab14open0 = (real) ( poDerivedStates->duration_ );
	tableinfo.tab14open1 = (real) ( poDerivedStates->duration_ever_union_TRUE );
	tableinfo.tab14open2 = (real) ( poDerivedStates->duration_has_spouse_TRUE );
	tableinfo.tab15events = -1;
	tableinfo.tab15open0 = (real) ( poDerivedStates->transitions_integer_age_14_15 );
	tableinfo.tab15open1 = (real) ( poDerivedStates->transitions_is_alive_TRUE_FALSE );
	tableinfo.tab15open2 = (real) ( poDerivedStates->transitions_survived_early_years_FALSE_TRUE );
	tableinfo.tab15open3 = (real) ( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_adult_survival );
	tableinfo.tab15open4 = (real) ( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_ind_hci );
	tableinfo.tab15open5 = (real) ( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_is_stunted );
	tableinfo.tab15open6 = (real) ( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling );
	tableinfo.tab15open7 = (real) ( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_survived_early_years );
	tableinfo.tab15open8 = (real) ( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling );
	tableinfo.tab16events = -1;
	tableinfo.tab16open0 = (real) ( poDerivedStates->transitions_integer_age_14_15 );
	tableinfo.tab16open1 = (real) ( poDerivedStates->transitions_is_alive_TRUE_FALSE );
	tableinfo.tab16open2 = (real) ( poDerivedStates->transitions_survived_early_years_FALSE_TRUE );
	tableinfo.tab16open3 = (real) ( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_adult_survival );
	tableinfo.tab16open4 = (real) ( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_ind_hci );
	tableinfo.tab16open5 = (real) ( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_is_stunted );
	tableinfo.tab16open6 = (real) ( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling );
	tableinfo.tab16open7 = (real) ( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_survived_early_years );
	tableinfo.tab16open8 = (real) ( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling );
	tableinfo.tab17events = -1;
	tableinfo.tab17open0 = (real) ( poDerivedStates->duration_ );
	tableinfo.tab18events = -1;
	tableinfo.tab18open0 = (real) ( poDerivedStates->transitions_is_alive_FALSE_TRUE );
	tableinfo.tab19events = -1;
	tableinfo.tab19open0 = (real) ( is_alive_value_in );
	tableinfo.tab20events = -1;
	tableinfo.tab20open0 = (real) ( is_stunted_value_in );
	tableinfo.tab20open1 = (real) ( table_unit_tabStuntingSexRegYob );
	tableinfo.tab21events = -1;
	tableinfo.tab21open0 = (real) ( is_stunted_value_in );
	tableinfo.tab21open1 = (real) ( table_unit_tabStuntingSexRegMotherYob );
	tableinfo.tab22events = -1;
	tableinfo.tab22open0 = (real) ( poDerivedStates->duration_ );
	tableinfo.tab22open1 = (real) ( poDerivedStates->duration_got_prenat_care_GPC_YES );
	tableinfo.tab23events = -1;
	tableinfo.tab23open0 = (real) ( poDerivedStates->duration_ );
	tableinfo.tab23open1 = (real) ( poDerivedStates->duration_is_immunized_TRUE );
	tableinfo.tab24events = -1;
	tableinfo.tab24open0 = (real) ( poDerivedStates->duration_ );
	tableinfo.tab24open1 = (real) ( poDerivedStates->duration_is_immunized_TRUE );
	poPrev = NULL;
	poNext = gpoEventQueue->poHeadActor;
	if ( gpoEventQueue->poHeadActor != NULL ) {
		gpoEventQueue->poHeadActor->poPrev = this;
	}
	gpoEventQueue->poHeadActor = this;
	bUpdates = TRUE;
}

void Person::FinishSpecial()
{
	CountFinishedActor(3);
	bUpdates = FALSE;
	gpoEventQueue->DestroyEvent( poBackMigrationEvent );
	gpoEventQueue->DestroyEvent( poBirthdayEvent );
	gpoEventQueue->DestroyEvent( poDeathAtMaxLifespanEvent );
	gpoEventQueue->DestroyEvent( poEmigrationEvent );
	gpoEventQueue->DestroyEvent( poFertilityBaseEvent );
	gpoEventQueue->DestroyEvent( poFertilityDetailedEvent );
	gpoEventQueue->DestroyEvent( poFirstImmigrationFromPoolEvent );
	gpoEventQueue->DestroyEvent( poFirstImmigrationFromScratchEvent );
	gpoEventQueue->DestroyEvent( poFirstUnionFormationEvent );
	gpoEventQueue->DestroyEvent( poLeavingHomeEvent );
	gpoEventQueue->DestroyEvent( poMigrationEvent );
	gpoEventQueue->DestroyEvent( poMortalityBaseEvent );
	gpoEventQueue->DestroyEvent( poMortalityDetailedEvent );
	gpoEventQueue->DestroyEvent( poSetAliveEvent );
	gpoEventQueue->DestroyEvent( poSetPreschoolYearsEvent );
	gpoEventQueue->DestroyEvent(poScheduleEvent);

	// close all updates
	if ( gplTableIndexes[0] >= 0 && ( tableinfo.tab0dim_filter || table_filter_tabPreSchool ) ) {
		if ( tableinfo.tab0events != -1 ) {
			tableinfo.tab0open0 = (real) CloseUpdate( 0, 0, tableinfo.tab0open0, tableinfo.tab0close0, tableinfo.tab0dim_filter, tableinfo.tab0cell, tableinfo.tab0total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab0open1 = (real) CloseUpdate( 0, 1, tableinfo.tab0open1, tableinfo.tab0close1, tableinfo.tab0dim_filter, tableinfo.tab0cell, tableinfo.tab0total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab0open2 = (real) CloseUpdate( 0, 2, tableinfo.tab0open2, tableinfo.tab0close2, tableinfo.tab0dim_filter, tableinfo.tab0cell, tableinfo.tab0total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab0open3 = (real) CloseUpdate( 0, 3, tableinfo.tab0open3, tableinfo.tab0close3, tableinfo.tab0dim_filter, tableinfo.tab0cell, tableinfo.tab0total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab0cell = (index) TableCellIndex( 0, 2,
			geo,
			tab_sim_yob - r_min_SIM_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab0total_cell = (index) tableinfo.tab0cell;
		}
		else {
			tableinfo.tab0total_cell = (index) TableTotalCellIndex( 0, 2,
			geo,
			tab_sim_yob - r_min_SIM_YEAR_RANGE );
		}
		CloseUpdate( 0, 0, tableinfo.tab0open0, poDerivedStates->Get_duration_(), table_filter_tabPreSchool, tableinfo.tab0cell, tableinfo.tab0total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 0, 1, tableinfo.tab0open1, poDerivedStates->Get_duration_years_preschool_0(), table_filter_tabPreSchool, tableinfo.tab0cell, tableinfo.tab0total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 0, 2, tableinfo.tab0open2, poDerivedStates->Get_duration_years_preschool_1(), table_filter_tabPreSchool, tableinfo.tab0cell, tableinfo.tab0total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 0, 3, tableinfo.tab0open3, poDerivedStates->Get_duration_years_preschool_2(), table_filter_tabPreSchool, tableinfo.tab0cell, tableinfo.tab0total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[1] >= 0 && ( tableinfo.tab1dim_filter || table_filter_tabEducationFateGeobirYob ) ) {
		if ( tableinfo.tab1events != -1 ) {
			tableinfo.tab1open0 = (real) CloseUpdate( 1, 0, tableinfo.tab1open0, tableinfo.tab1close0, tableinfo.tab1dim_filter, tableinfo.tab1cell, tableinfo.tab1total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab1open1 = (real) CloseUpdate( 1, 1, tableinfo.tab1open1, tableinfo.tab1close1, tableinfo.tab1dim_filter, tableinfo.tab1cell, tableinfo.tab1total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab1open2 = (real) CloseUpdate( 1, 2, tableinfo.tab1open2, tableinfo.tab1close2, tableinfo.tab1dim_filter, tableinfo.tab1cell, tableinfo.tab1total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab1open3 = (real) CloseUpdate( 1, 3, tableinfo.tab1open3, tableinfo.tab1close3, tableinfo.tab1dim_filter, tableinfo.tab1cell, tableinfo.tab1total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab1cell = (index) TableCellIndex( 1, 3,
			sex,
			geo_birth,
			year_of_birth - r_min_ALL_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab1total_cell = (index) tableinfo.tab1cell;
		}
		else {
			tableinfo.tab1total_cell = (index) TableTotalCellIndex( 1, 3,
			sex,
			geo_birth,
			year_of_birth - r_min_ALL_YEAR_RANGE );
		}
		CloseUpdate( 1, 0, tableinfo.tab1open0, poDerivedStates->Get_duration_(), table_filter_tabEducationFateGeobirYob, tableinfo.tab1cell, tableinfo.tab1total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 1, 1, tableinfo.tab1open1, poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH(), table_filter_tabEducationFateGeobirYob, tableinfo.tab1cell, tableinfo.tab1total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 1, 2, tableinfo.tab1open2, poDerivedStates->Get_duration_educ_one_fate_EOL_LOW(), table_filter_tabEducationFateGeobirYob, tableinfo.tab1cell, tableinfo.tab1total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 1, 3, tableinfo.tab1open3, poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM(), table_filter_tabEducationFateGeobirYob, tableinfo.tab1cell, tableinfo.tab1total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[2] >= 0 && ( tableinfo.tab2dim_filter || table_filter_PopPyramidByEduc ) ) {
		if ( tableinfo.tab2events != -1 ) {
			tableinfo.tab2open0 = (real) CloseUpdate( 2, 0, tableinfo.tab2open0, tableinfo.tab2close0, tableinfo.tab2dim_filter, tableinfo.tab2cell, tableinfo.tab2total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab2open1 = (real) CloseUpdate( 2, 1, tableinfo.tab2open1, tableinfo.tab2close1, tableinfo.tab2dim_filter, tableinfo.tab2cell, tableinfo.tab2total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab2open2 = (real) CloseUpdate( 2, 2, tableinfo.tab2open2, tableinfo.tab2close2, tableinfo.tab2dim_filter, tableinfo.tab2cell, tableinfo.tab2total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab2cell = (index) TableCellIndex( 2, 4,
			sel_years,
			sex,
			geo_nat,
			poDerivedStates->split_integer_age_TAB_AGEC5 );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab2total_cell = (index) tableinfo.tab2cell;
		}
		else {
			tableinfo.tab2total_cell = (index) TableTotalCellIndex( 2, 4,
			sel_years,
			sex,
			geo_nat,
			poDerivedStates->split_integer_age_TAB_AGEC5 );
		}
		CloseUpdate( 2, 0, tableinfo.tab2open0, poDerivedStates->Get_duration_educ_one_level_EOL_HIGH(), table_filter_PopPyramidByEduc, tableinfo.tab2cell, tableinfo.tab2total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 2, 1, tableinfo.tab2open1, poDerivedStates->Get_duration_educ_one_level_EOL_LOW(), table_filter_PopPyramidByEduc, tableinfo.tab2cell, tableinfo.tab2total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 2, 2, tableinfo.tab2open2, poDerivedStates->Get_duration_educ_one_level_EOL_MEDIUM(), table_filter_PopPyramidByEduc, tableinfo.tab2cell, tableinfo.tab2total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[3] >= 0 && ( tableinfo.tab3dim_filter || table_filter_TabEducFateByGroup ) ) {
		if ( tableinfo.tab3events != -1 ) {
			poStateFunctions->Set_table_unit_TabEducFateByGroup(table_unit_TabEducFateByGroup+1);
			tableinfo.tab3close0 = (real) table_unit_TabEducFateByGroup;
			tableinfo.tab3open0 = (real) CloseUpdate( 3, 0, tableinfo.tab3open0, tableinfo.tab3close0, tableinfo.tab3dim_filter, tableinfo.tab3cell, tableinfo.tab3total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab3cell = (index) TableCellIndex( 3, 3,
			educ1_group,
			tab_sim_yob - r_min_SIM_YEAR_RANGE,
			educ_one_fate );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab3total_cell = (index) tableinfo.tab3cell;
		}
		else {
			tableinfo.tab3total_cell = (index) TableTotalCellIndex( 3, 3,
			educ1_group,
			tab_sim_yob - r_min_SIM_YEAR_RANGE,
			educ_one_fate );
		}
		poStateFunctions->Set_table_unit_TabEducFateByGroup(table_unit_TabEducFateByGroup+1);
		CloseUpdate( 3, 0, tableinfo.tab3open0, table_unit_TabEducFateByGroup, table_filter_TabEducFateByGroup, tableinfo.tab3cell, tableinfo.tab3total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[4] >= 0 && ( tableinfo.tab4dim_filter || table_filter_TabPrimarySchoolPlanning ) ) {
		if ( tableinfo.tab4events != -1 ) {
			tableinfo.tab4open0 = (real) CloseUpdate( 4, 0, tableinfo.tab4open0, tableinfo.tab4close0, tableinfo.tab4dim_filter, tableinfo.tab4cell, tableinfo.tab4total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab4open1 = (real) CloseUpdate( 4, 1, tableinfo.tab4open1, tableinfo.tab4close1, tableinfo.tab4dim_filter, tableinfo.tab4cell, tableinfo.tab4total_cell, actor_weight, actor_subsample_weight );
			poStateFunctions->Set_table_unit_TabPrimarySchoolPlanning(table_unit_TabPrimarySchoolPlanning+1);
			tableinfo.tab4close2 = (real) table_unit_TabPrimarySchoolPlanning;
			tableinfo.tab4open2 = (real) CloseUpdate( 4, 2, tableinfo.tab4open2, tableinfo.tab4close2, tableinfo.tab4dim_filter, tableinfo.tab4cell, tableinfo.tab4total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab4cell = (index) TableCellIndex( 4, 3,
			geo_nat,
			sim_year - r_min_SIM_YEAR_RANGE,
			educ_one_grade - r_min_EDUC_ONE_GRADE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab4total_cell = (index) tableinfo.tab4cell;
		}
		else {
			tableinfo.tab4total_cell = (index) TableTotalCellIndex( 4, 3,
			geo_nat,
			sim_year - r_min_SIM_YEAR_RANGE,
			educ_one_grade - r_min_EDUC_ONE_GRADE );
		}
		CloseUpdate( 4, 0, tableinfo.tab4open0, educ_one_required_rooms_value_out, table_filter_TabPrimarySchoolPlanning, tableinfo.tab4cell, tableinfo.tab4total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 4, 1, tableinfo.tab4open1, educ_one_required_teachers_value_out, table_filter_TabPrimarySchoolPlanning, tableinfo.tab4cell, tableinfo.tab4total_cell, actor_weight, actor_subsample_weight );
		poStateFunctions->Set_table_unit_TabPrimarySchoolPlanning(table_unit_TabPrimarySchoolPlanning+1);
		CloseUpdate( 4, 2, tableinfo.tab4open2, table_unit_TabPrimarySchoolPlanning, table_filter_TabPrimarySchoolPlanning, tableinfo.tab4cell, tableinfo.tab4total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[5] >= 0 && ( tableinfo.tab5dim_filter || table_filter_TabPopProvAgeEducSex ) ) {
		if ( tableinfo.tab5events != -1 ) {
			tableinfo.tab5open0 = (real) CloseUpdate( 5, 0, tableinfo.tab5open0, tableinfo.tab5close0, tableinfo.tab5dim_filter, tableinfo.tab5cell, tableinfo.tab5total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab5open1 = (real) CloseUpdate( 5, 1, tableinfo.tab5open1, tableinfo.tab5close1, tableinfo.tab5dim_filter, tableinfo.tab5cell, tableinfo.tab5total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab5open2 = (real) CloseUpdate( 5, 2, tableinfo.tab5open2, tableinfo.tab5close2, tableinfo.tab5dim_filter, tableinfo.tab5cell, tableinfo.tab5total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab5open3 = (real) CloseUpdate( 5, 3, tableinfo.tab5open3, tableinfo.tab5close3, tableinfo.tab5dim_filter, tableinfo.tab5cell, tableinfo.tab5total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab5open4 = (real) CloseUpdate( 5, 4, tableinfo.tab5open4, tableinfo.tab5close4, tableinfo.tab5dim_filter, tableinfo.tab5cell, tableinfo.tab5total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab5cell = (index) TableCellIndex( 5, 3,
			geo_nat,
			sex,
			sim_year - r_min_SIM_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab5total_cell = (index) tableinfo.tab5cell;
		}
		else {
			tableinfo.tab5total_cell = (index) TableTotalCellIndex( 5, 3,
			geo_nat,
			sex,
			sim_year - r_min_SIM_YEAR_RANGE );
		}
		CloseUpdate( 5, 0, tableinfo.tab5open0, poDerivedStates->Get_duration_is_aded0_TRUE(), table_filter_TabPopProvAgeEducSex, tableinfo.tab5cell, tableinfo.tab5total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 5, 1, tableinfo.tab5open1, poDerivedStates->Get_duration_is_aded1_TRUE(), table_filter_TabPopProvAgeEducSex, tableinfo.tab5cell, tableinfo.tab5total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 5, 2, tableinfo.tab5open2, poDerivedStates->Get_duration_is_aded2_TRUE(), table_filter_TabPopProvAgeEducSex, tableinfo.tab5cell, tableinfo.tab5total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 5, 3, tableinfo.tab5open3, poDerivedStates->Get_duration_is_child_TRUE(), table_filter_TabPopProvAgeEducSex, tableinfo.tab5cell, tableinfo.tab5total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 5, 4, tableinfo.tab5open4, poDerivedStates->Get_duration_is_old_TRUE(), table_filter_TabPopProvAgeEducSex, tableinfo.tab5cell, tableinfo.tab5total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[6] >= 0 && ( tableinfo.tab6dim_filter || table_filter_TabEducFateDistrYob ) ) {
		if ( tableinfo.tab6events != -1 ) {
			poStateFunctions->Set_table_unit_TabEducFateDistrYob(table_unit_TabEducFateDistrYob+1);
			tableinfo.tab6close0 = (real) table_unit_TabEducFateDistrYob;
			tableinfo.tab6open0 = (real) CloseUpdate( 6, 0, tableinfo.tab6open0, tableinfo.tab6close0, tableinfo.tab6dim_filter, tableinfo.tab6cell, tableinfo.tab6total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab6cell = (index) TableCellIndex( 6, 3,
			poDerivedStates->aggregate_region_birth_REGION_NAT,
			year_of_birth - r_min_ALL_YEAR_RANGE,
			educ_one_fate );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab6total_cell = (index) tableinfo.tab6cell;
		}
		else {
			tableinfo.tab6total_cell = (index) TableTotalCellIndex( 6, 3,
			poDerivedStates->aggregate_region_birth_REGION_NAT,
			year_of_birth - r_min_ALL_YEAR_RANGE,
			educ_one_fate );
		}
		poStateFunctions->Set_table_unit_TabEducFateDistrYob(table_unit_TabEducFateDistrYob+1);
		CloseUpdate( 6, 0, tableinfo.tab6open0, table_unit_TabEducFateDistrYob, table_filter_TabEducFateDistrYob, tableinfo.tab6cell, tableinfo.tab6total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[7] >= 0 && ( tableinfo.tab7dim_filter || table_filter_TabEduc15ByDistrict ) ) {
		if ( tableinfo.tab7events != -1 ) {
			tableinfo.tab7open0 = (real) CloseUpdate( 7, 0, tableinfo.tab7open0, tableinfo.tab7close0, tableinfo.tab7dim_filter, tableinfo.tab7cell, tableinfo.tab7total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab7open1 = (real) CloseUpdate( 7, 1, tableinfo.tab7open1, tableinfo.tab7close1, tableinfo.tab7dim_filter, tableinfo.tab7cell, tableinfo.tab7total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab7open2 = (real) CloseUpdate( 7, 2, tableinfo.tab7open2, tableinfo.tab7close2, tableinfo.tab7dim_filter, tableinfo.tab7cell, tableinfo.tab7total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab7open3 = (real) CloseUpdate( 7, 3, tableinfo.tab7open3, tableinfo.tab7close3, tableinfo.tab7dim_filter, tableinfo.tab7cell, tableinfo.tab7total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab7cell = (index) TableCellIndex( 7, 3,
			geo_nat,
			sex,
			sim_year - r_min_SIM_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab7total_cell = (index) tableinfo.tab7cell;
		}
		else {
			tableinfo.tab7total_cell = (index) TableTotalCellIndex( 7, 3,
			geo_nat,
			sex,
			sim_year - r_min_SIM_YEAR_RANGE );
		}
		CloseUpdate( 7, 0, tableinfo.tab7open0, poDerivedStates->Get_duration_(), table_filter_TabEduc15ByDistrict, tableinfo.tab7cell, tableinfo.tab7total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 7, 1, tableinfo.tab7open1, poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH(), table_filter_TabEduc15ByDistrict, tableinfo.tab7cell, tableinfo.tab7total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 7, 2, tableinfo.tab7open2, poDerivedStates->Get_duration_educ_one_fate_EOL_LOW(), table_filter_TabEduc15ByDistrict, tableinfo.tab7cell, tableinfo.tab7total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 7, 3, tableinfo.tab7open3, poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM(), table_filter_TabEduc15ByDistrict, tableinfo.tab7cell, tableinfo.tab7total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[8] >= 0 && ( tableinfo.tab8dim_filter || table_filter_TabEduc15ByDistrictBirth ) ) {
		if ( tableinfo.tab8events != -1 ) {
			tableinfo.tab8open0 = (real) CloseUpdate( 8, 0, tableinfo.tab8open0, tableinfo.tab8close0, tableinfo.tab8dim_filter, tableinfo.tab8cell, tableinfo.tab8total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab8open1 = (real) CloseUpdate( 8, 1, tableinfo.tab8open1, tableinfo.tab8close1, tableinfo.tab8dim_filter, tableinfo.tab8cell, tableinfo.tab8total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab8open2 = (real) CloseUpdate( 8, 2, tableinfo.tab8open2, tableinfo.tab8close2, tableinfo.tab8dim_filter, tableinfo.tab8cell, tableinfo.tab8total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab8open3 = (real) CloseUpdate( 8, 3, tableinfo.tab8open3, tableinfo.tab8close3, tableinfo.tab8dim_filter, tableinfo.tab8cell, tableinfo.tab8total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab8cell = (index) TableCellIndex( 8, 3,
			geo_birth,
			sex,
			sim_year - r_min_SIM_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab8total_cell = (index) tableinfo.tab8cell;
		}
		else {
			tableinfo.tab8total_cell = (index) TableTotalCellIndex( 8, 3,
			geo_birth,
			sex,
			sim_year - r_min_SIM_YEAR_RANGE );
		}
		CloseUpdate( 8, 0, tableinfo.tab8open0, poDerivedStates->Get_duration_(), table_filter_TabEduc15ByDistrictBirth, tableinfo.tab8cell, tableinfo.tab8total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 8, 1, tableinfo.tab8open1, poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH(), table_filter_TabEduc15ByDistrictBirth, tableinfo.tab8cell, tableinfo.tab8total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 8, 2, tableinfo.tab8open2, poDerivedStates->Get_duration_educ_one_fate_EOL_LOW(), table_filter_TabEduc15ByDistrictBirth, tableinfo.tab8cell, tableinfo.tab8total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 8, 3, tableinfo.tab8open3, poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM(), table_filter_TabEduc15ByDistrictBirth, tableinfo.tab8cell, tableinfo.tab8total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[9] >= 0 && ( tableinfo.tab9dim_filter || table_filter_TabPrimSchoolEntries ) ) {
		if ( tableinfo.tab9events != -1 ) {
			tableinfo.tab9open0 = (real) CloseUpdate( 9, 0, tableinfo.tab9open0, tableinfo.tab9close0, tableinfo.tab9dim_filter, tableinfo.tab9cell, tableinfo.tab9total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab9cell = (index) TableCellIndex( 9, 2,
			sim_year - r_min_SIM_YEAR_RANGE,
			geo_nat );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab9total_cell = (index) tableinfo.tab9cell;
		}
		else {
			tableinfo.tab9total_cell = (index) TableTotalCellIndex( 9, 2,
			sim_year - r_min_SIM_YEAR_RANGE,
			geo_nat );
		}
		CloseUpdate( 9, 0, tableinfo.tab9open0, poDerivedStates->entrances_educ_one_grade_attended_1, table_filter_TabPrimSchoolEntries, tableinfo.tab9cell, tableinfo.tab9total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[10] >= 0 && ( tableinfo.tab10dim_filter || table_filter_TabPrimSchoolGraduations ) ) {
		if ( tableinfo.tab10events != -1 ) {
			tableinfo.tab10open0 = (real) CloseUpdate( 10, 0, tableinfo.tab10open0, tableinfo.tab10close0, tableinfo.tab10dim_filter, tableinfo.tab10cell, tableinfo.tab10total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab10cell = (index) TableCellIndex( 10, 2,
			sim_year - r_min_SIM_YEAR_RANGE,
			geo_nat );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab10total_cell = (index) tableinfo.tab10cell;
		}
		else {
			tableinfo.tab10total_cell = (index) TableTotalCellIndex( 10, 2,
			sim_year - r_min_SIM_YEAR_RANGE,
			geo_nat );
		}
		CloseUpdate( 10, 0, tableinfo.tab10open0, poDerivedStates->entrances_has_passed_primary_TRUE, table_filter_TabPrimSchoolGraduations, tableinfo.tab10cell, tableinfo.tab10total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[11] >= 0 && ( tableinfo.tab11dim_filter || table_filter_TabPrimSchoolOutOfSchool9to11 ) ) {
		if ( tableinfo.tab11events != -1 ) {
			tableinfo.tab11open0 = (real) CloseUpdate( 11, 0, tableinfo.tab11open0, tableinfo.tab11close0, tableinfo.tab11dim_filter, tableinfo.tab11cell, tableinfo.tab11total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab11open1 = (real) CloseUpdate( 11, 1, tableinfo.tab11open1, tableinfo.tab11close1, tableinfo.tab11dim_filter, tableinfo.tab11cell, tableinfo.tab11total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab11cell = (index) TableCellIndex( 11, 2,
			sim_year - r_min_SIM_YEAR_RANGE,
			geo_nat );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab11total_cell = (index) tableinfo.tab11cell;
		}
		else {
			tableinfo.tab11total_cell = (index) TableTotalCellIndex( 11, 2,
			sim_year - r_min_SIM_YEAR_RANGE,
			geo_nat );
		}
		CloseUpdate( 11, 0, tableinfo.tab11open0, poDerivedStates->Get_duration_(), table_filter_TabPrimSchoolOutOfSchool9to11, tableinfo.tab11cell, tableinfo.tab11total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 11, 1, tableinfo.tab11open1, poDerivedStates->Get_duration_tab_out_of_primary_TRUE(), table_filter_TabPrimSchoolOutOfSchool9to11, tableinfo.tab11cell, tableinfo.tab11total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[12] >= 0 && ( tableinfo.tab12dim_filter || table_filter_TabSchool2TrackTab ) ) {
		if ( tableinfo.tab12events != -1 ) {
			poStateFunctions->Set_table_unit_TabSchool2TrackTab(table_unit_TabSchool2TrackTab+1);
			tableinfo.tab12close0 = (real) table_unit_TabSchool2TrackTab;
			tableinfo.tab12open0 = (real) CloseUpdate( 12, 0, tableinfo.tab12open0, tableinfo.tab12close0, tableinfo.tab12dim_filter, tableinfo.tab12cell, tableinfo.tab12total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab12cell = (index) TableCellIndex( 12, 2,
			sim_year - r_min_SIM_YEAR_RANGE,
			educ_two_grade_attended - r_min_EDUC_TWO_GRADE0 );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab12total_cell = (index) tableinfo.tab12cell;
		}
		else {
			tableinfo.tab12total_cell = (index) TableTotalCellIndex( 12, 2,
			sim_year - r_min_SIM_YEAR_RANGE,
			educ_two_grade_attended - r_min_EDUC_TWO_GRADE0 );
		}
		poStateFunctions->Set_table_unit_TabSchool2TrackTab(table_unit_TabSchool2TrackTab+1);
		CloseUpdate( 12, 0, tableinfo.tab12open0, table_unit_TabSchool2TrackTab, table_filter_TabSchool2TrackTab, tableinfo.tab12cell, tableinfo.tab12total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[13] >= 0 && ( tableinfo.tab13dim_filter || table_filter_TabSchool2AttainmentsTab ) ) {
		if ( tableinfo.tab13events != -1 ) {
			tableinfo.tab13open0 = (real) CloseUpdate( 13, 0, tableinfo.tab13open0, tableinfo.tab13close0, tableinfo.tab13dim_filter, tableinfo.tab13cell, tableinfo.tab13total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab13open1 = (real) CloseUpdate( 13, 1, tableinfo.tab13open1, tableinfo.tab13close1, tableinfo.tab13dim_filter, tableinfo.tab13cell, tableinfo.tab13total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab13open2 = (real) CloseUpdate( 13, 2, tableinfo.tab13open2, tableinfo.tab13close2, tableinfo.tab13dim_filter, tableinfo.tab13cell, tableinfo.tab13total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab13open3 = (real) CloseUpdate( 13, 3, tableinfo.tab13open3, tableinfo.tab13close3, tableinfo.tab13dim_filter, tableinfo.tab13cell, tableinfo.tab13total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab13cell = (index) TableCellIndex( 13, 2,
			sex,
			sim_year - r_min_SIM_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab13total_cell = (index) tableinfo.tab13cell;
		}
		else {
			tableinfo.tab13total_cell = (index) TableTotalCellIndex( 13, 2,
			sex,
			sim_year - r_min_SIM_YEAR_RANGE );
		}
		CloseUpdate( 13, 0, tableinfo.tab13open0, poDerivedStates->Get_duration_(), table_filter_TabSchool2AttainmentsTab, tableinfo.tab13cell, tableinfo.tab13total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 13, 1, tableinfo.tab13open1, poDerivedStates->Get_duration_educ_two_level_ETL_GRADUATE(), table_filter_TabSchool2AttainmentsTab, tableinfo.tab13cell, tableinfo.tab13total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 13, 2, tableinfo.tab13open2, poDerivedStates->Get_duration_educ_two_level_ETL_INCOMPLETE(), table_filter_TabSchool2AttainmentsTab, tableinfo.tab13cell, tableinfo.tab13total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 13, 3, tableinfo.tab13open3, poDerivedStates->Get_duration_educ_two_level_ETL_NEVER(), table_filter_TabSchool2AttainmentsTab, tableinfo.tab13cell, tableinfo.tab13total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[14] >= 0 && ( tableinfo.tab14dim_filter || table_filter_tabHavingSpouse ) ) {
		if ( tableinfo.tab14events != -1 ) {
			tableinfo.tab14open0 = (real) CloseUpdate( 14, 0, tableinfo.tab14open0, tableinfo.tab14close0, tableinfo.tab14dim_filter, tableinfo.tab14cell, tableinfo.tab14total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab14open1 = (real) CloseUpdate( 14, 1, tableinfo.tab14open1, tableinfo.tab14close1, tableinfo.tab14dim_filter, tableinfo.tab14cell, tableinfo.tab14total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab14open2 = (real) CloseUpdate( 14, 2, tableinfo.tab14open2, tableinfo.tab14close2, tableinfo.tab14dim_filter, tableinfo.tab14cell, tableinfo.tab14total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab14cell = (index) TableCellIndex( 14, 3,
			sex,
			integer_age - r_min_AGE_RANGE,
			sim_year - r_min_SIM_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab14total_cell = (index) tableinfo.tab14cell;
		}
		else {
			tableinfo.tab14total_cell = (index) TableTotalCellIndex( 14, 3,
			sex,
			integer_age - r_min_AGE_RANGE,
			sim_year - r_min_SIM_YEAR_RANGE );
		}
		CloseUpdate( 14, 0, tableinfo.tab14open0, poDerivedStates->Get_duration_(), table_filter_tabHavingSpouse, tableinfo.tab14cell, tableinfo.tab14total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 14, 1, tableinfo.tab14open1, poDerivedStates->Get_duration_ever_union_TRUE(), table_filter_tabHavingSpouse, tableinfo.tab14cell, tableinfo.tab14total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 14, 2, tableinfo.tab14open2, poDerivedStates->Get_duration_has_spouse_TRUE(), table_filter_tabHavingSpouse, tableinfo.tab14cell, tableinfo.tab14total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[15] >= 0 && ( tableinfo.tab15dim_filter || table_filter_tabHCI ) ) {
		if ( tableinfo.tab15events != -1 ) {
			tableinfo.tab15open0 = (real) CloseUpdate( 15, 0, tableinfo.tab15open0, tableinfo.tab15close0, tableinfo.tab15dim_filter, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab15open1 = (real) CloseUpdate( 15, 1, tableinfo.tab15open1, tableinfo.tab15close1, tableinfo.tab15dim_filter, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab15open2 = (real) CloseUpdate( 15, 2, tableinfo.tab15open2, tableinfo.tab15close2, tableinfo.tab15dim_filter, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab15open3 = (real) CloseUpdate( 15, 3, tableinfo.tab15open3, tableinfo.tab15close3, tableinfo.tab15dim_filter, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab15open4 = (real) CloseUpdate( 15, 4, tableinfo.tab15open4, tableinfo.tab15close4, tableinfo.tab15dim_filter, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab15open5 = (real) CloseUpdate( 15, 5, tableinfo.tab15open5, tableinfo.tab15close5, tableinfo.tab15dim_filter, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab15open6 = (real) CloseUpdate( 15, 6, tableinfo.tab15open6, tableinfo.tab15close6, tableinfo.tab15dim_filter, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab15open7 = (real) CloseUpdate( 15, 7, tableinfo.tab15open7, tableinfo.tab15close7, tableinfo.tab15dim_filter, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab15open8 = (real) CloseUpdate( 15, 8, tableinfo.tab15open8, tableinfo.tab15close8, tableinfo.tab15dim_filter, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab15cell = (index) TableCellIndex( 15, 3,
			sex,
			poDerivedStates->aggregate_region_birth_REGION_NAT,
			tab_sim_yob - r_min_SIM_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab15total_cell = (index) tableinfo.tab15cell;
		}
		else {
			tableinfo.tab15total_cell = (index) TableTotalCellIndex( 15, 3,
			sex,
			poDerivedStates->aggregate_region_birth_REGION_NAT,
			tab_sim_yob - r_min_SIM_YEAR_RANGE );
		}
		CloseUpdate( 15, 0, tableinfo.tab15open0, poDerivedStates->transitions_integer_age_14_15, table_filter_tabHCI, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 15, 1, tableinfo.tab15open1, poDerivedStates->transitions_is_alive_TRUE_FALSE, table_filter_tabHCI, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 15, 2, tableinfo.tab15open2, poDerivedStates->transitions_survived_early_years_FALSE_TRUE, table_filter_tabHCI, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 15, 3, tableinfo.tab15open3, poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_adult_survival, table_filter_tabHCI, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 15, 4, tableinfo.tab15open4, poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_ind_hci, table_filter_tabHCI, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 15, 5, tableinfo.tab15open5, poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_is_stunted, table_filter_tabHCI, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 15, 6, tableinfo.tab15open6, poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling, table_filter_tabHCI, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 15, 7, tableinfo.tab15open7, poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_survived_early_years, table_filter_tabHCI, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 15, 8, tableinfo.tab15open8, poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling, table_filter_tabHCI, tableinfo.tab15cell, tableinfo.tab15total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[16] >= 0 && ( tableinfo.tab16dim_filter || table_filter_TabHCIDistrict ) ) {
		if ( tableinfo.tab16events != -1 ) {
			tableinfo.tab16open0 = (real) CloseUpdate( 16, 0, tableinfo.tab16open0, tableinfo.tab16close0, tableinfo.tab16dim_filter, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab16open1 = (real) CloseUpdate( 16, 1, tableinfo.tab16open1, tableinfo.tab16close1, tableinfo.tab16dim_filter, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab16open2 = (real) CloseUpdate( 16, 2, tableinfo.tab16open2, tableinfo.tab16close2, tableinfo.tab16dim_filter, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab16open3 = (real) CloseUpdate( 16, 3, tableinfo.tab16open3, tableinfo.tab16close3, tableinfo.tab16dim_filter, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab16open4 = (real) CloseUpdate( 16, 4, tableinfo.tab16open4, tableinfo.tab16close4, tableinfo.tab16dim_filter, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab16open5 = (real) CloseUpdate( 16, 5, tableinfo.tab16open5, tableinfo.tab16close5, tableinfo.tab16dim_filter, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab16open6 = (real) CloseUpdate( 16, 6, tableinfo.tab16open6, tableinfo.tab16close6, tableinfo.tab16dim_filter, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab16open7 = (real) CloseUpdate( 16, 7, tableinfo.tab16open7, tableinfo.tab16close7, tableinfo.tab16dim_filter, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab16open8 = (real) CloseUpdate( 16, 8, tableinfo.tab16open8, tableinfo.tab16close8, tableinfo.tab16dim_filter, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab16cell = (index) TableCellIndex( 16, 3,
			sex,
			poDerivedStates->aggregate_geo_birth_GEO_NAT,
			tab_sim_yob - r_min_SIM_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab16total_cell = (index) tableinfo.tab16cell;
		}
		else {
			tableinfo.tab16total_cell = (index) TableTotalCellIndex( 16, 3,
			sex,
			poDerivedStates->aggregate_geo_birth_GEO_NAT,
			tab_sim_yob - r_min_SIM_YEAR_RANGE );
		}
		CloseUpdate( 16, 0, tableinfo.tab16open0, poDerivedStates->transitions_integer_age_14_15, table_filter_TabHCIDistrict, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 16, 1, tableinfo.tab16open1, poDerivedStates->transitions_is_alive_TRUE_FALSE, table_filter_TabHCIDistrict, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 16, 2, tableinfo.tab16open2, poDerivedStates->transitions_survived_early_years_FALSE_TRUE, table_filter_TabHCIDistrict, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 16, 3, tableinfo.tab16open3, poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_adult_survival, table_filter_TabHCIDistrict, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 16, 4, tableinfo.tab16open4, poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_ind_hci, table_filter_TabHCIDistrict, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 16, 5, tableinfo.tab16open5, poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_is_stunted, table_filter_TabHCIDistrict, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 16, 6, tableinfo.tab16open6, poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling, table_filter_TabHCIDistrict, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 16, 7, tableinfo.tab16open7, poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_survived_early_years, table_filter_TabHCIDistrict, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 16, 8, tableinfo.tab16open8, poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling, table_filter_TabHCIDistrict, tableinfo.tab16cell, tableinfo.tab16total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[17] >= 0 && ( tableinfo.tab17dim_filter || table_filter_tabPopulationYearPlace ) ) {
		if ( tableinfo.tab17events != -1 ) {
			tableinfo.tab17open0 = (real) CloseUpdate( 17, 0, tableinfo.tab17open0, tableinfo.tab17close0, tableinfo.tab17dim_filter, tableinfo.tab17cell, tableinfo.tab17total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab17cell = (index) TableCellIndex( 17, 2,
			geo,
			sim_year - r_min_SIM_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab17total_cell = (index) tableinfo.tab17cell;
		}
		else {
			tableinfo.tab17total_cell = (index) TableTotalCellIndex( 17, 2,
			geo,
			sim_year - r_min_SIM_YEAR_RANGE );
		}
		CloseUpdate( 17, 0, tableinfo.tab17open0, poDerivedStates->Get_duration_(), table_filter_tabPopulationYearPlace, tableinfo.tab17cell, tableinfo.tab17total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[18] >= 0 && ( tableinfo.tab18dim_filter || table_filter_tabBirthsYearPlace ) ) {
		if ( tableinfo.tab18events != -1 ) {
			tableinfo.tab18open0 = (real) CloseUpdate( 18, 0, tableinfo.tab18open0, tableinfo.tab18close0, tableinfo.tab18dim_filter, tableinfo.tab18cell, tableinfo.tab18total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab18cell = (index) TableCellIndex( 18, 2,
			geo,
			sim_year - r_min_SIM_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab18total_cell = (index) tableinfo.tab18cell;
		}
		else {
			tableinfo.tab18total_cell = (index) TableTotalCellIndex( 18, 2,
			geo,
			sim_year - r_min_SIM_YEAR_RANGE );
		}
		CloseUpdate( 18, 0, tableinfo.tab18open0, poDerivedStates->transitions_is_alive_FALSE_TRUE, table_filter_tabBirthsYearPlace, tableinfo.tab18cell, tableinfo.tab18total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[19] >= 0 && ( tableinfo.tab19dim_filter || table_filter_tabMigrationOriginDestination ) ) {
		if ( tableinfo.tab19events != -1 ) {
			tableinfo.tab19open0 = (real) CloseUpdate( 19, 0, tableinfo.tab19open0, tableinfo.tab19close0, tableinfo.tab19dim_filter, tableinfo.tab19cell, tableinfo.tab19total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab19cell = (index) TableCellIndex( 19, 3,
			selected_year10,
			geo_prev,
			geo );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab19total_cell = (index) tableinfo.tab19cell;
		}
		else {
			tableinfo.tab19total_cell = (index) TableTotalCellIndex( 19, 3,
			selected_year10,
			geo_prev,
			geo );
		}
		CloseUpdate( 19, 0, tableinfo.tab19open0, is_alive_value_in, table_filter_tabMigrationOriginDestination, tableinfo.tab19cell, tableinfo.tab19total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[20] >= 0 && ( tableinfo.tab20dim_filter || table_filter_tabStuntingSexRegYob ) ) {
		if ( tableinfo.tab20events != -1 ) {
			tableinfo.tab20open0 = (real) CloseUpdate( 20, 0, tableinfo.tab20open0, tableinfo.tab20close0, tableinfo.tab20dim_filter, tableinfo.tab20cell, tableinfo.tab20total_cell, actor_weight, actor_subsample_weight );
			poStateFunctions->Set_table_unit_tabStuntingSexRegYob(table_unit_tabStuntingSexRegYob+1);
			tableinfo.tab20close1 = (real) table_unit_tabStuntingSexRegYob;
			tableinfo.tab20open1 = (real) CloseUpdate( 20, 1, tableinfo.tab20open1, tableinfo.tab20close1, tableinfo.tab20dim_filter, tableinfo.tab20cell, tableinfo.tab20total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab20cell = (index) TableCellIndex( 20, 2,
			region_nat,
			sim_year - r_min_SIM_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab20total_cell = (index) tableinfo.tab20cell;
		}
		else {
			tableinfo.tab20total_cell = (index) TableTotalCellIndex( 20, 2,
			region_nat,
			sim_year - r_min_SIM_YEAR_RANGE );
		}
		CloseUpdate( 20, 0, tableinfo.tab20open0, is_stunted_value_in, table_filter_tabStuntingSexRegYob, tableinfo.tab20cell, tableinfo.tab20total_cell, actor_weight, actor_subsample_weight );
		poStateFunctions->Set_table_unit_tabStuntingSexRegYob(table_unit_tabStuntingSexRegYob+1);
		CloseUpdate( 20, 1, tableinfo.tab20open1, table_unit_tabStuntingSexRegYob, table_filter_tabStuntingSexRegYob, tableinfo.tab20cell, tableinfo.tab20total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[21] >= 0 && ( tableinfo.tab21dim_filter || table_filter_tabStuntingSexRegMotherYob ) ) {
		if ( tableinfo.tab21events != -1 ) {
			tableinfo.tab21open0 = (real) CloseUpdate( 21, 0, tableinfo.tab21open0, tableinfo.tab21close0, tableinfo.tab21dim_filter, tableinfo.tab21cell, tableinfo.tab21total_cell, actor_weight, actor_subsample_weight );
			poStateFunctions->Set_table_unit_tabStuntingSexRegMotherYob(table_unit_tabStuntingSexRegMotherYob+1);
			tableinfo.tab21close1 = (real) table_unit_tabStuntingSexRegMotherYob;
			tableinfo.tab21open1 = (real) CloseUpdate( 21, 1, tableinfo.tab21open1, tableinfo.tab21close1, tableinfo.tab21dim_filter, tableinfo.tab21cell, tableinfo.tab21total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab21cell = (index) TableCellIndex( 21, 4,
			sex,
			educ_mother,
			region_nat,
			sim_year - r_min_SIM_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab21total_cell = (index) tableinfo.tab21cell;
		}
		else {
			tableinfo.tab21total_cell = (index) TableTotalCellIndex( 21, 4,
			sex,
			educ_mother,
			region_nat,
			sim_year - r_min_SIM_YEAR_RANGE );
		}
		CloseUpdate( 21, 0, tableinfo.tab21open0, is_stunted_value_in, table_filter_tabStuntingSexRegMotherYob, tableinfo.tab21cell, tableinfo.tab21total_cell, actor_weight, actor_subsample_weight );
		poStateFunctions->Set_table_unit_tabStuntingSexRegMotherYob(table_unit_tabStuntingSexRegMotherYob+1);
		CloseUpdate( 21, 1, tableinfo.tab21open1, table_unit_tabStuntingSexRegMotherYob, table_filter_tabStuntingSexRegMotherYob, tableinfo.tab21cell, tableinfo.tab21total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[22] >= 0 && ( tableinfo.tab22dim_filter || table_filter_TabPrenatCare ) ) {
		if ( tableinfo.tab22events != -1 ) {
			tableinfo.tab22open0 = (real) CloseUpdate( 22, 0, tableinfo.tab22open0, tableinfo.tab22close0, tableinfo.tab22dim_filter, tableinfo.tab22cell, tableinfo.tab22total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab22open1 = (real) CloseUpdate( 22, 1, tableinfo.tab22open1, tableinfo.tab22close1, tableinfo.tab22dim_filter, tableinfo.tab22cell, tableinfo.tab22total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab22cell = (index) TableCellIndex( 22, 3,
			region_nat,
			educ_mother,
			sim_year - r_min_SIM_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab22total_cell = (index) tableinfo.tab22cell;
		}
		else {
			tableinfo.tab22total_cell = (index) TableTotalCellIndex( 22, 3,
			region_nat,
			educ_mother,
			sim_year - r_min_SIM_YEAR_RANGE );
		}
		CloseUpdate( 22, 0, tableinfo.tab22open0, poDerivedStates->Get_duration_(), table_filter_TabPrenatCare, tableinfo.tab22cell, tableinfo.tab22total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 22, 1, tableinfo.tab22open1, poDerivedStates->Get_duration_got_prenat_care_GPC_YES(), table_filter_TabPrenatCare, tableinfo.tab22cell, tableinfo.tab22total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[23] >= 0 && ( tableinfo.tab23dim_filter || table_filter_TabImmunizationChildren ) ) {
		if ( tableinfo.tab23events != -1 ) {
			tableinfo.tab23open0 = (real) CloseUpdate( 23, 0, tableinfo.tab23open0, tableinfo.tab23close0, tableinfo.tab23dim_filter, tableinfo.tab23cell, tableinfo.tab23total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab23open1 = (real) CloseUpdate( 23, 1, tableinfo.tab23open1, tableinfo.tab23close1, tableinfo.tab23dim_filter, tableinfo.tab23cell, tableinfo.tab23total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab23cell = (index) TableCellIndex( 23, 6,
			sex,
			ethnicity_short,
			got_prenat_care,
			region_nat,
			educ_mother,
			sim_year - r_min_SIM_YEAR_RANGE );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab23total_cell = (index) tableinfo.tab23cell;
		}
		else {
			tableinfo.tab23total_cell = (index) TableTotalCellIndex( 23, 6,
			sex,
			ethnicity_short,
			got_prenat_care,
			region_nat,
			educ_mother,
			sim_year - r_min_SIM_YEAR_RANGE );
		}
		CloseUpdate( 23, 0, tableinfo.tab23open0, poDerivedStates->Get_duration_(), table_filter_TabImmunizationChildren, tableinfo.tab23cell, tableinfo.tab23total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 23, 1, tableinfo.tab23open1, poDerivedStates->Get_duration_is_immunized_TRUE(), table_filter_TabImmunizationChildren, tableinfo.tab23cell, tableinfo.tab23total_cell, actor_weight, actor_subsample_weight );
	}
	if ( gplTableIndexes[24] >= 0 && ( tableinfo.tab24dim_filter || table_filter_TabChildVaccination ) ) {
		if ( tableinfo.tab24events != -1 ) {
			tableinfo.tab24open0 = (real) CloseUpdate( 24, 0, tableinfo.tab24open0, tableinfo.tab24close0, tableinfo.tab24dim_filter, tableinfo.tab24cell, tableinfo.tab24total_cell, actor_weight, actor_subsample_weight );
			tableinfo.tab24open1 = (real) CloseUpdate( 24, 1, tableinfo.tab24open1, tableinfo.tab24close1, tableinfo.tab24dim_filter, tableinfo.tab24cell, tableinfo.tab24total_cell, actor_weight, actor_subsample_weight );
		}
		tableinfo.tab24cell = (index) TableCellIndex( 24, 2,
			sim_year - r_min_SIM_YEAR_RANGE,
			region_nat );
		if ( gpModelApp->m_nSubSamples == 1 ) {
			tableinfo.tab24total_cell = (index) tableinfo.tab24cell;
		}
		else {
			tableinfo.tab24total_cell = (index) TableTotalCellIndex( 24, 2,
			sim_year - r_min_SIM_YEAR_RANGE,
			region_nat );
		}
		CloseUpdate( 24, 0, tableinfo.tab24open0, poDerivedStates->Get_duration_(), table_filter_TabChildVaccination, tableinfo.tab24cell, tableinfo.tab24total_cell, actor_weight, actor_subsample_weight );
		CloseUpdate( 24, 1, tableinfo.tab24open1, poDerivedStates->Get_duration_is_immunized_TRUE(), table_filter_TabChildVaccination, tableinfo.tab24cell, tableinfo.tab24total_cell, actor_weight, actor_subsample_weight );
	}
	Set_lBioFather( NULL );
	Set_lBioMother( NULL );
	Set_lCalibrator( NULL );
	Set_lHHFather( NULL );
	Set_lHHMother( NULL );
	Set_lSpouse( NULL );
	Set_lStartValues( NULL );
	mlBioFatherChildren->Finish();
	mlBioMotherChildren->Finish();
	mlHHFatherChildren->Finish();
	mlHHMotherChildren->Finish();
	if (poasAllFertilePersonsForFertilityAlignment != NULL)
	{	asAllFertilePersonsForFertilityAlignment[is_resident]->Remove(poasAllFertilePersonsForFertilityAlignment);
		poasAllFertilePersonsForFertilityAlignment = NULL;
	}
	if (poasAllFertilePersonsForFertilityAlignmentByAge != NULL)
	{	asAllFertilePersonsForFertilityAlignmentByAge[RANGE_POS(FERTILE_AGE_RANGE,fertile_age)][is_resident]->Remove(poasAllFertilePersonsForFertilityAlignmentByAge);
		poasAllFertilePersonsForFertilityAlignmentByAge = NULL;
	}
	if (poasAllPerson != NULL)
	{	asAllPerson->Remove(poasAllPerson);
		poasAllPerson = NULL;
	}
	if (poasAvailableMale != NULL)
	{	asAvailableMale[RANGE_POS(MALE_SPOUSE_AGE,male_age)][geo]->Remove(poasAvailableMale);
		poasAvailableMale = NULL;
	}
	if (poasAvailableMaleByType != NULL)
	{	asAvailableMaleByType[RANGE_POS(MALE_SPOUSE_AGE,male_age)][spouse_group][geo]->Remove(poasAvailableMaleByType);
		poasAvailableMaleByType = NULL;
	}
	if (poasAvailableMaleByTypeImmiScratch != NULL)
	{	asAvailableMaleByTypeImmiScratch[RANGE_POS(MALE_SPOUSE_AGE,male_age)][spouse_group][geo][RANGE_POS(SIM_YEAR_RANGE,year_of_first_immigration)]->Remove(poasAvailableMaleByTypeImmiScratch);
		poasAvailableMaleByTypeImmiScratch = NULL;
	}
	if (poasAvailableMaleImmiScratch != NULL)
	{	asAvailableMaleImmiScratch[RANGE_POS(MALE_SPOUSE_AGE,male_age)][geo][RANGE_POS(SIM_YEAR_RANGE,year_of_first_immigration)]->Remove(poasAvailableMaleImmiScratch);
		poasAvailableMaleImmiScratch = NULL;
	}
	if (poasFemaleInUnionByAgeAndPartnerAge != NULL)
	{	asFemaleInUnionByAgeAndPartnerAge[RANGE_POS(FEMALE_SPOUSE_AGE,female_spouse_age)][RANGE_POS(MALE_SPOUSE_AGE,male_spouse_age)]->Remove(poasFemaleInUnionByAgeAndPartnerAge);
		poasFemaleInUnionByAgeAndPartnerAge = NULL;
	}
	if (poasFemaleInUnionByAgeAndPartnerAgeImmiScratch != NULL)
	{	asFemaleInUnionByAgeAndPartnerAgeImmiScratch[RANGE_POS(FEMALE_SPOUSE_AGE,female_spouse_age)][RANGE_POS(MALE_SPOUSE_AGE,male_spouse_age)][RANGE_POS(SIM_YEAR_RANGE,year_of_first_immigration)]->Remove(poasFemaleInUnionByAgeAndPartnerAgeImmiScratch);
		poasFemaleInUnionByAgeAndPartnerAgeImmiScratch = NULL;
	}
	if (poasPotentialImmigrantMothers != NULL)
	{	asPotentialImmigrantMothers[RANGE_POS(AGE_RANGE,integer_age)][RANGE_POS(SIM_YEAR_RANGE,year_of_first_immigration)][geo]->Remove(poasPotentialImmigrantMothers);
		poasPotentialImmigrantMothers = NULL;
	}
	if (poasPotentialSchoolOneStudent != NULL)
	{	asPotentialSchoolOneStudent->Remove(poasPotentialSchoolOneStudent);
		poasPotentialSchoolOneStudent = NULL;
	}
	if (poasPotentialSchoolOneStudentToProcess != NULL)
	{	asPotentialSchoolOneStudentToProcess->Remove(poasPotentialSchoolOneStudentToProcess);
		poasPotentialSchoolOneStudentToProcess = NULL;
	}
	if (poasPotentialSchoolTwoStudent != NULL)
	{	asPotentialSchoolTwoStudent->Remove(poasPotentialSchoolTwoStudent);
		poasPotentialSchoolTwoStudent = NULL;
	}
	if (poasPotentialSchoolTwoStudentToProcess != NULL)
	{	asPotentialSchoolTwoStudentToProcess->Remove(poasPotentialSchoolTwoStudentToProcess);
		poasPotentialSchoolTwoStudentToProcess = NULL;
	}
	if (poasResidentBabies != NULL)
	{	asResidentBabies->Remove(poasResidentBabies);
		poasResidentBabies = NULL;
	}
	if (poasResidentsAge0Sex != NULL)
	{	asResidentsAge0Sex[sex]->Remove(poasResidentsAge0Sex);
		poasResidentsAge0Sex = NULL;
	}
	if (poasResidentsAge0SexGeo != NULL)
	{	asResidentsAge0SexGeo[sex][geo_birth]->Remove(poasResidentsAge0SexGeo);
		poasResidentsAge0SexGeo = NULL;
	}
	if (poasSimBornAge0 != NULL)
	{	asSimBornAge0[sex][geo][educ1_group]->Remove(poasSimBornAge0);
		poasSimBornAge0 = NULL;
	}
	if (poasWantToMove != NULL)
	{	asWantToMove[geo][geo_want_to_move][sex][RANGE_POS(AGE_MIG,age_mig)]->Remove(poasWantToMove);
		poasWantToMove = NULL;
	}
	if (poasWomenNoChildren != NULL)
	{	asWomenNoChildren[union1_group][RANGE_POS(FEMALE_SPOUSE_AGE,female_spouse_age)][has_spouse]->Remove(poasWomenNoChildren);
		poasWomenNoChildren = NULL;
	}
	if (poasWomenWithChildren != NULL)
	{	asWomenWithChildren[union1_group][child_agegr][moth_agegr][has_spouse]->Remove(poasWomenWithChildren);
		poasWomenWithChildren = NULL;
	}

	// Remove actor from master list of active actors
	if ( poPrev != NULL )
	{
		poPrev->poNext = poNext;
	}
	else
	{
		if (gpoEventQueue->poHeadActor != NULL && gpoEventQueue->poHeadActor == this)
		{
			gpoEventQueue->poHeadActor = poNext;
		}
	}
	if ( poNext != NULL )
	{
		poNext->poPrev = poPrev;
	}
	poPrev = NULL;

	// Add actor to head of empty list of actors of this type
	poNext = gpoEmptyPerson;
	gpoEmptyPerson = this;
}

#pragma optimize( "", on )

void Person::WaitUntil( double target_time )
{
	//don't increase gdEventsForTabulation if special events because they are all treated as one event for tabulation purposes
	if (!gpoEventQueue->m_bSpecialEvents)
	{
		gdEventsForTabulation++;
	}
	__events = gdEventsForTabulation;
	if ( time < (TIME) target_time )
	{
		__time = (TIME) target_time;
		WaitSpecial( (TIME) target_time - time );
		poStateFunctions->Set_time( (TIME) target_time );
	}
}

Person *Person::UpdateTime(  )
{
	__events = gdEventsForTabulation;
	if ( time < (TIME) gdThreadCurrentTime ) {
		__time = (TIME) gdThreadCurrentTime;
		WaitSpecial( (TIME) gdThreadCurrentTime - time );
		poStateFunctions->Set_time( (TIME) gdThreadCurrentTime );
	}
	return this;
}

void Person::WaitUntilThisActor( double target_time )
{
	WaitUntil(target_time);
}

void Person::WaitSpecial( TIME wait_time )
{

	// increment the age
		poStateFunctions->Set_age( age + wait_time );
}

TIME Person::timeScheduleEvent()
{
	if (schedule_event_time == TIME_INFINITE)
	{
		TIME event_time = TIME_INFINITE;
		event_time = min( poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time, poDerivedStates->sssplit_age_AGE_MIG_PART_time );
		event_time = min( event_time, poDerivedStates->trigger_changes_calendar_year_time );
		event_time = min( event_time, poDerivedStates->trigger_changes_geo_time );
		event_time = min( event_time, poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE_time );
		schedule_event_time = event_time;
	}
	return schedule_event_time;
}

void Person::ScheduleEvent()
{

	if (poDerivedStates->trigger_changes_calendar_year_time == schedule_event_time)
	{
		//CHKSUM(time, 160);
		TRCEVNT(case_seed, _T("Person"), actor_id, _T("scheduled - 2"), time);
		gszRuntimeEventName = _T("Person.scheduled - 2");
		poDerivedStates->trigger_changes_calendar_year_time = TIME_INFINITE;
		poDerivedStates->Set_trigger_changes_calendar_year( FALSE );

	}

	if (poDerivedStates->trigger_changes_geo_time == schedule_event_time)
	{
		//CHKSUM(time, 161);
		TRCEVNT(case_seed, _T("Person"), actor_id, _T("scheduled - 3"), time);
		gszRuntimeEventName = _T("Person.scheduled - 3");
		poDerivedStates->trigger_changes_geo_time = TIME_INFINITE;
		poDerivedStates->Set_trigger_changes_geo( FALSE );

	}

	if (poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE_time == schedule_event_time)
	{
		//CHKSUM(time, 162);
		TRCEVNT(case_seed, _T("Person"), actor_id, _T("scheduled - 4"), time);
		gszRuntimeEventName = _T("Person.scheduled - 4");
		poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE_time = TIME_INFINITE;
		poDerivedStates->Set_trigger_entrances_is_ready_for_birthtables_TRUE( FALSE );

	}

	if (poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time == schedule_event_time)
	{
		//CHKSUM(time, 125);
		TRCEVNT(case_seed, _T("Person"), actor_id, _T("scheduled - 0"), time);
		gszRuntimeEventName = _T("Person.scheduled - 0");
		poDerivedStates->Set_sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY( poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY + 1 );
		if ( poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY < 5 ) {
		poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time = WAIT(TIME(POINTS(DUR_TIME_IN_PARITY)[poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY] - POINTS(DUR_TIME_IN_PARITY)[poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY - 1]));
		}
		else
		{		poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time = TIME_INFINITE;
		}
	}

	if (poDerivedStates->sssplit_age_AGE_MIG_PART_time == schedule_event_time)
	{
		//CHKSUM(time, 126);
		TRCEVNT(case_seed, _T("Person"), actor_id, _T("scheduled - 1"), time);
		gszRuntimeEventName = _T("Person.scheduled - 1");
		poDerivedStates->Set_sssplit_age_AGE_MIG_PART( poDerivedStates->sssplit_age_AGE_MIG_PART + 1 );
		if ( poDerivedStates->sssplit_age_AGE_MIG_PART < 12 ) {
		poDerivedStates->sssplit_age_AGE_MIG_PART_time = WAIT(TIME(POINTS(AGE_MIG_PART)[poDerivedStates->sssplit_age_AGE_MIG_PART] - POINTS(AGE_MIG_PART)[poDerivedStates->sssplit_age_AGE_MIG_PART - 1]));
		}
		else
		{		poDerivedStates->sssplit_age_AGE_MIG_PART_time = TIME_INFINITE;
		}
	}
	schedule_event_time = TIME_INFINITE;
}

double Person::EventTime( int nEventNum, int *pnEventInfo )
{
	double	event_time = double( TIME_INFINITE );

	gdRuntimeTime = (double) time;
	gnRuntimeActorId = actor_id;
	switch( nEventNum ) {
		case 0:
			gszRuntimeEventName = _T("Person.BackMigrationEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeBackMigrationEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeBackMigrationEvent"), event_time);
			break;
		case 1:
			gszRuntimeEventName = _T("Person.BirthdayEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeBirthdayEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeBirthdayEvent"), event_time);
			break;
		case 6:
			gszRuntimeEventName = _T("Person.DeathAtMaxLifespanEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeDeathAtMaxLifespanEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeDeathAtMaxLifespanEvent"), event_time);
			break;
		case 7:
			gszRuntimeEventName = _T("Person.EmigrationEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeEmigrationEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeEmigrationEvent"), event_time);
			break;
		case 10:
			gszRuntimeEventName = _T("Person.FertilityBaseEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeFertilityBaseEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeFertilityBaseEvent"), event_time);
			break;
		case 11:
			gszRuntimeEventName = _T("Person.FertilityDetailedEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeFertilityDetailedEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeFertilityDetailedEvent"), event_time);
			break;
		case 12:
			gszRuntimeEventName = _T("Person.FirstImmigrationFromPoolEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeFirstImmigrationFromPoolEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeFirstImmigrationFromPoolEvent"), event_time);
			break;
		case 13:
			gszRuntimeEventName = _T("Person.FirstImmigrationFromScratchEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeFirstImmigrationFromScratchEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeFirstImmigrationFromScratchEvent"), event_time);
			break;
		case 15:
			gszRuntimeEventName = _T("Person.FirstUnionFormationEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeFirstUnionFormationEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeFirstUnionFormationEvent"), event_time);
			break;
		case 19:
			gszRuntimeEventName = _T("Person.LeavingHomeEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeLeavingHomeEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeLeavingHomeEvent"), event_time);
			break;
		case 20:
			gszRuntimeEventName = _T("Person.MigrationEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeMigrationEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeMigrationEvent"), event_time);
			break;
		case 21:
			gszRuntimeEventName = _T("Person.MortalityBaseEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeMortalityBaseEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeMortalityBaseEvent"), event_time);
			break;
		case 22:
			gszRuntimeEventName = _T("Person.MortalityDetailedEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeMortalityDetailedEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeMortalityDetailedEvent"), event_time);
			break;
		case 24:
			gszRuntimeEventName = _T("Person.SetAliveEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeSetAliveEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeSetAliveEvent"), event_time);
			break;
		case 27:
			gszRuntimeEventName = _T("Person.SetPreschoolYearsEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeSetPreschoolYearsEvent() );
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("timeSetPreschoolYearsEvent"), event_time);
			break;
		case 23:
			gszRuntimeEventName = _T("ScheduleEvent(time)");
			event_time = (double) timeScheduleEvent();
			break;
		case -1:
		default:;
	}
	return (TIME) event_time;
}

void Person::Implement( int nEventNum, int nEventInfo )
{
	gdThreadEvents += 1;
	gdRuntimeTime = (double) time;
	gnRuntimeActorId = actor_id;
	switch( nEventNum )
	{
		case 0:
			CHKSUM(time, 0);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.BackMigrationEvent"), time);
			gszRuntimeEventName = _T("Person.BackMigrationEvent(implement)");
			BackMigrationEvent();
			break;
		case 1:
			CHKSUM(time, 1);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.BirthdayEvent"), time);
			gszRuntimeEventName = _T("Person.BirthdayEvent(implement)");
			BirthdayEvent();
			break;
		case 6:
			CHKSUM(time, 2);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.DeathAtMaxLifespanEvent"), time);
			gszRuntimeEventName = _T("Person.DeathAtMaxLifespanEvent(implement)");
			DeathAtMaxLifespanEvent();
			break;
		case 7:
			CHKSUM(time, 3);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.EmigrationEvent"), time);
			gszRuntimeEventName = _T("Person.EmigrationEvent(implement)");
			EmigrationEvent();
			break;
		case 10:
			CHKSUM(time, 4);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.FertilityBaseEvent"), time);
			gszRuntimeEventName = _T("Person.FertilityBaseEvent(implement)");
			FertilityBaseEvent();
			break;
		case 11:
			CHKSUM(time, 5);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.FertilityDetailedEvent"), time);
			gszRuntimeEventName = _T("Person.FertilityDetailedEvent(implement)");
			FertilityDetailedEvent();
			break;
		case 12:
			CHKSUM(time, 6);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.FirstImmigrationFromPoolEvent"), time);
			gszRuntimeEventName = _T("Person.FirstImmigrationFromPoolEvent(implement)");
			FirstImmigrationFromPoolEvent();
			break;
		case 13:
			CHKSUM(time, 7);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.FirstImmigrationFromScratchEvent"), time);
			gszRuntimeEventName = _T("Person.FirstImmigrationFromScratchEvent(implement)");
			FirstImmigrationFromScratchEvent();
			break;
		case 15:
			CHKSUM(time, 8);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.FirstUnionFormationEvent"), time);
			gszRuntimeEventName = _T("Person.FirstUnionFormationEvent(implement)");
			FirstUnionFormationEvent();
			break;
		case 19:
			CHKSUM(time, 9);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.LeavingHomeEvent"), time);
			gszRuntimeEventName = _T("Person.LeavingHomeEvent(implement)");
			LeavingHomeEvent();
			break;
		case 20:
			CHKSUM(time, 10);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.MigrationEvent"), time);
			gszRuntimeEventName = _T("Person.MigrationEvent(implement)");
			MigrationEvent();
			break;
		case 21:
			CHKSUM(time, 11);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.MortalityBaseEvent"), time);
			gszRuntimeEventName = _T("Person.MortalityBaseEvent(implement)");
			MortalityBaseEvent();
			break;
		case 22:
			CHKSUM(time, 12);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.MortalityDetailedEvent"), time);
			gszRuntimeEventName = _T("Person.MortalityDetailedEvent(implement)");
			MortalityDetailedEvent();
			break;
		case 24:
			CHKSUM(time, 13);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.SetAliveEvent"), time);
			gszRuntimeEventName = _T("Person.SetAliveEvent(implement)");
			SetAliveEvent();
			break;
		case 27:
			CHKSUM(time, 14);
			TRCEVNT(case_seed, _T("Person"), actor_id, _T("Person.SetPreschoolYearsEvent"), time);
			gszRuntimeEventName = _T("Person.SetPreschoolYearsEvent(implement)");
			SetPreschoolYearsEvent();
			break;
		case 23:
			//CHKSUM(time, 15);
			ScheduleEvent();
			break;
		case -1:
		default:;
	}
}

double Person::Set_actor_subsample_weight( double value )
{
	return poStateFunctions->Set_actor_subsample_weight( value);
}

double Person::Set_actor_weight( double value )
{
	return poStateFunctions->Set_actor_weight( value);
}

TIME PersonStateFunctions::Set_age( TIME value, bool bStartFunction )
{

	if ( poParent->age != value || bStartFunction )
	{
		if (bStartFunction)
		{
			poDerivedStates->Set_sssplit_age_AGE_MIG_PART( SPLIT(value, AGE_MIG_PART) );
			if ( poDerivedStates->sssplit_age_AGE_MIG_PART < 12 )
			{
				poDerivedStates->sssplit_age_AGE_MIG_PART_time = poParent->time +  POINTS(AGE_MIG_PART)[poDerivedStates->sssplit_age_AGE_MIG_PART] - value ;

			}
			else
			{
				poDerivedStates->sssplit_age_AGE_MIG_PART_time = TIME_INFINITE;
			}
			if (poParent->schedule_event_time > poDerivedStates->sssplit_age_AGE_MIG_PART_time)
			{
				poParent->schedule_event_time = poDerivedStates->sssplit_age_AGE_MIG_PART_time;
			}
		}


		TIME	OldValue;
		OldValue = poParent->age;
		poParent->age = value;

	}
	return poParent->age = value;
}

ALL_YEAR_RANGE PersonStateFunctions::Set_calendar_year( ALL_YEAR_RANGE value )
{

	CHKLMT(value, calendar_year, r_min_ALL_YEAR_RANGE, r_max_ALL_YEAR_RANGE);
	if ( value < r_min_ALL_YEAR_RANGE ) {
		value = r_min_ALL_YEAR_RANGE;
	}
	else if ( value > r_max_ALL_YEAR_RANGE ) {
		value = r_max_ALL_YEAR_RANGE;
	}
	if ( poParent->calendar_year != value )
	{
		bool bScheduleTimeChange = false;	//indicates that schedule event should be recalculated
		if (	(poParent->schedule_event_time == poDerivedStates->trigger_changes_calendar_year_time))
		{
			//force timeScheduleEvent to calculate min time
			poParent->schedule_event_time = TIME_INFINITE;
		//indicate that schedule event needs to be recalculated
			bScheduleTimeChange = true;
		}
#line 55 "../code/TablesPopulation.mpp"
		Set_table_filter_tabBirthsYearPlace( WITHIN(SIM_YEAR_RANGE,value) );
#line 5580 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 13 "../code/TablesFamily.mpp"
		Set_table_filter_tabHavingSpouse( WITHIN(SIM_YEAR_RANGE,value)&& poParent->is_resident );
#line 5583 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 66 "../code/TablesPopulation.mpp"
		Set_table_filter_tabMigrationOriginDestination( WITHIN(SIM_YEAR_RANGE,value)&& poDerivedStates->trigger_changes_geo&& poParent->in_selected_year );
#line 5586 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 45 "../code/TablesPopulation.mpp"
		Set_table_filter_tabPopulationYearPlace( WITHIN(SIM_YEAR_RANGE,value) );
#line 5589 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		ALL_YEAR_RANGE	OldValue;
		OldValue = poParent->calendar_year;
		poParent->calendar_year = value;

#line 23 "../code/TablesPopulation.mpp"
		Set_all_year( COERCE(ALL_YEAR_RANGE,poParent->calendar_year) );
#line 5598 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 43 "../code/EducationPrimaryPlanning.mpp"
		Set_educ_one_required_rooms( ( poParent->is_educ_one_student&& poParent->is_resident )? 1.0/ gprParam->Educ1Infrastructure[CHKBND(SI_ROOMS, 2, "Educ1Infrastructure", 0)][CHKBND(poParent->geo_nat, 12, "Educ1Infrastructure", 1)][CHKBND(RANGE_POS(SIM_YEAR_RANGE,poParent->calendar_year), 151, "Educ1Infrastructure", 2)]: 0.0 );
#line 5601 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 39 "../code/EducationPrimaryPlanning.mpp"
		Set_educ_one_required_teachers( ( poParent->is_educ_one_student&& poParent->is_resident )? 1.0/ gprParam->Educ1Infrastructure[CHKBND(SI_TEACHERS, 2, "Educ1Infrastructure", 0)][CHKBND(poParent->geo_nat, 12, "Educ1Infrastructure", 1)][CHKBND(RANGE_POS(SIM_YEAR_RANGE,poParent->calendar_year), 151, "Educ1Infrastructure", 2)]: 0.0 );
#line 5604 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 18 "../code/TablesEducation.mpp"
		Set_in_projected_time( ( poParent->calendar_year>= MIN(SIM_YEAR_RANGE) ) );
#line 5607 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 36 "../code/_ContextABC.mpp"
		Set_in_sel_years( poParent->calendar_year== 2000|| poParent->calendar_year== 2010|| poParent->calendar_year== 2020|| poParent->calendar_year== 2030|| poParent->calendar_year== 2040|| poParent->calendar_year== 2050|| poParent->calendar_year== 2060 );
#line 5610 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 34 "../code/TablesPopulation.mpp"
		Set_in_selected_year( ( poParent->selected_year10!= SY40|| poParent->calendar_year== MIN(SIM_YEAR_RANGE)+ 40 ) );
#line 5613 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 47 "../code/FertilityGeneral.mpp"
		Set_is_fertile( ( poParent->sex== FEMALE&& WITHIN(FERTILE_AGE_RANGE,poParent->integer_age)&& ( poParent->creation_type== CT_BIRTH|| ( poParent->creation_type== CT_START&& poParent->calendar_year>= MIN(SIM_YEAR_RANGE) )|| ( poParent->creation_type== CT_POOL&& poParent->ever_resident )|| ( poParent->creation_type== CT_SCRATCH&& poParent->ever_resident ) ) )? TRUE: FALSE );
#line 5616 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 45 "../code/MortalityGeneral.mpp"
		Set_is_mortal( ( poParent->creation_type== CT_BIRTH|| ( poParent->creation_type== CT_START&& poParent->calendar_year>= MIN(SIM_YEAR_RANGE) )|| ( poParent->creation_type== CT_POOL&& poParent->ever_resident )|| ( poParent->creation_type== CT_SCRATCH&& poParent->ever_resident ) )? TRUE: FALSE );
#line 5619 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 44 "../code/_ContextABC.mpp"
		Set_sel_years( ( poParent->calendar_year== 2000 )? TSY_2000: ( poParent->calendar_year== 2010 )? TSY_2010: ( poParent->calendar_year== 2020 )? TSY_2020: ( poParent->calendar_year== 2030 )? TSY_2030: ( poParent->calendar_year== 2040 )? TSY_2040: ( poParent->calendar_year== 2050 )? TSY_2050: TSY_2060 );
#line 5622 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 31 "../code/TablesPopulation.mpp"
		Set_selected_year10( ( poParent->calendar_year== MIN(SIM_YEAR_RANGE) )? SY00: ( poParent->calendar_year== MIN(SIM_YEAR_RANGE)+ 10 )? SY10: ( poParent->calendar_year== MIN(SIM_YEAR_RANGE)+ 20 )? SY20: ( poParent->calendar_year== MIN(SIM_YEAR_RANGE)+ 30 )? SY30: SY40 );
#line 5625 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 24 "../code/TablesPopulation.mpp"
		Set_sim_year( COERCE(SIM_YEAR_RANGE,poParent->calendar_year) );
#line 5628 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		poDerivedStates->Set_trigger_changes_calendar_year( TRUE );
		poDerivedStates->trigger_changes_calendar_year_time = poParent->__time;
		if ( poParent->schedule_event_time > poDerivedStates->trigger_changes_calendar_year_time)
		{
			poParent->schedule_event_time = poDerivedStates->trigger_changes_calendar_year_time;
			//indicate that schedule event should be recalculated
			bScheduleTimeChange = true;
		}
		if ( poParent->poBackMigrationEvent->cRecalc == 0 ) {
			poParent->poBackMigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poBackMigrationEvent );
		};
		if ( poParent->poFertilityBaseEvent->cRecalc == 0 ) {
			poParent->poFertilityBaseEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFertilityBaseEvent );
		};
		if ( poParent->poFertilityDetailedEvent->cRecalc == 0 ) {
			poParent->poFertilityDetailedEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFertilityDetailedEvent );
		};
		if ( poParent->poFirstUnionFormationEvent->cRecalc == 0 ) {
			poParent->poFirstUnionFormationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstUnionFormationEvent );
		};
		if ( poParent->poMigrationEvent->cRecalc == 0 ) {
			poParent->poMigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMigrationEvent );
		};
		if ( poParent->poMortalityBaseEvent->cRecalc == 0 ) {
			poParent->poMortalityBaseEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMortalityBaseEvent );
		};
		if ( poParent->poMortalityDetailedEvent->cRecalc == 0 ) {
			poParent->poMortalityDetailedEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMortalityDetailedEvent );
		};

		if (bScheduleTimeChange && poParent->poScheduleEvent->cRecalc == 0)
		{
			poParent->poScheduleEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent(poParent->poScheduleEvent);
		}
	}
	return poParent->calendar_year = value;
}

CREATION_TYPE PersonStateFunctions::Set_creation_type( CREATION_TYPE value )
{

	CHKLMT(value, creation_type, 0, rc_size_CREATION_TYPE - 1);
	if ( value < 0 ) {
		value = (CREATION_TYPE) 0;
	}
	else if ( value >= rc_size_CREATION_TYPE ) {
		value = (CREATION_TYPE) ( rc_size_CREATION_TYPE - 1 );
	}
	if ( poParent->creation_type != value )
	{
#line 32 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleByTypeImmiScratch( poParent->is_alive&& value== CT_SCRATCH&& poParent->sex== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->is_mortal ) );
#line 5689 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleImmiScratch( poParent->is_alive&& value== CT_SCRATCH&& poParent->sex== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->is_mortal ) );
#line 5692 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 36 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch( poParent->is_alive&& value== CT_SCRATCH&& poParent->sex== FEMALE&& poParent->has_spouse&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->is_mortal ) );
#line 5695 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 10 "../code/ImmigrationFromScratch.mpp"
		Set_actorset_filter_asPotentialImmigrantMothers( poParent->is_alive&& poParent->sex== FEMALE&& value== CT_SCRATCH&& !( poParent->ever_resident )&& poParent->family_role!= FR_CHILD&& poParent->ever_union );
#line 5698 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 18 "../code/ChildVaccination.mpp"
		Set_actorset_filter_asResidentBabies( poParent->is_alive&& poParent->integer_age== 0&& poParent->is_resident&& value== CT_BIRTH );
#line 5701 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 35 "../code/EducationPrimaryTransmission.mpp"
		Set_actorset_filter_asSimBornAge0( poParent->is_alive&& value== CT_BIRTH&& poParent->integer_age== 0 );
#line 5704 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 46 "../code/TablesVaccination.mpp"
		Set_table_filter_TabChildVaccination( poParent->in_projected_time&& value== CT_BIRTH&& poParent->is_resident&& poParent->integer_age== 0 );
#line 5707 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 115 "../code/TablesEducation.mpp"
		Set_table_filter_TabEducFateByGroup( poParent->in_projected_time&& poDerivedStates->trigger_changes_calendar_year&& poParent->integer_age== 0&& value== CT_BIRTH );
#line 5710 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 31 "../code/TablesVaccination.mpp"
		Set_table_filter_TabImmunizationChildren( poParent->in_projected_time&& poParent->integer_age== 0&& value== CT_BIRTH&& poParent->is_resident );
#line 5713 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 19 "../code/TablesVaccination.mpp"
		Set_table_filter_TabPrenatCare( poParent->in_projected_time&& poParent->integer_age== 0&& value== CT_BIRTH&& poParent->is_resident );
#line 5716 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 27 "../code/TablesStunting.mpp"
		Set_table_filter_tabStuntingSexRegMotherYob( poParent->in_projected_time&& poParent->is_resident&& value== CT_BIRTH&& poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE );
#line 5719 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 17 "../code/TablesStunting.mpp"
		Set_table_filter_tabStuntingSexRegYob( poParent->in_projected_time&& poParent->is_resident&& value== CT_BIRTH&& poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE );
#line 5722 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		CREATION_TYPE	OldValue;
		OldValue = poParent->creation_type;
		poParent->creation_type = value;

#line 82 "../code/HumanCapitalIndex.mpp"
		Set_in_hci_sample( ( poParent->creation_type== CT_BIRTH&& poParent->geo_birth< SIZE(GEO_NAT) )? TRUE: FALSE );
#line 5731 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 47 "../code/FertilityGeneral.mpp"
		Set_is_fertile( ( poParent->sex== FEMALE&& WITHIN(FERTILE_AGE_RANGE,poParent->integer_age)&& ( poParent->creation_type== CT_BIRTH|| ( poParent->creation_type== CT_START&& poParent->calendar_year>= MIN(SIM_YEAR_RANGE) )|| ( poParent->creation_type== CT_POOL&& poParent->ever_resident )|| ( poParent->creation_type== CT_SCRATCH&& poParent->ever_resident ) ) )? TRUE: FALSE );
#line 5734 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 45 "../code/MortalityGeneral.mpp"
		Set_is_mortal( ( poParent->creation_type== CT_BIRTH|| ( poParent->creation_type== CT_START&& poParent->calendar_year>= MIN(SIM_YEAR_RANGE) )|| ( poParent->creation_type== CT_POOL&& poParent->ever_resident )|| ( poParent->creation_type== CT_SCRATCH&& poParent->ever_resident ) )? TRUE: FALSE );
#line 5737 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		if ( poParent->poFirstImmigrationFromPoolEvent->cRecalc == 0 ) {
			poParent->poFirstImmigrationFromPoolEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstImmigrationFromPoolEvent );
		};
		if ( poParent->poFirstImmigrationFromScratchEvent->cRecalc == 0 ) {
			poParent->poFirstImmigrationFromScratchEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstImmigrationFromScratchEvent );
		};
		if ( poParent->poFirstUnionFormationEvent->cRecalc == 0 ) {
			poParent->poFirstUnionFormationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstUnionFormationEvent );
		};
		if ( poParent->poSetPreschoolYearsEvent->cRecalc == 0 ) {
			poParent->poSetPreschoolYearsEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poSetPreschoolYearsEvent );
		};
	}
	return poParent->creation_type = value;
}

EDUC_ONE_LEVEL PersonStateFunctions::Set_educ_mother( EDUC_ONE_LEVEL value )
{

	CHKLMT(value, educ_mother, 0, rc_size_EDUC_ONE_LEVEL - 1);
	if ( value < 0 ) {
		value = (EDUC_ONE_LEVEL) 0;
	}
	else if ( value >= rc_size_EDUC_ONE_LEVEL ) {
		value = (EDUC_ONE_LEVEL) ( rc_size_EDUC_ONE_LEVEL - 1 );
	}
	if ( poParent->educ_mother != value )
	{


		EDUC_ONE_LEVEL	OldValue;
		OldValue = poParent->educ_mother;
		poParent->educ_mother = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[21] >= 0 && ( poParent->tableinfo.tab21dim_filter || poParent->table_filter_tabStuntingSexRegMotherYob ) ) {
				if ( poParent->tableinfo.tab21events == -1 || poParent->tableinfo.tab21events < poParent->__events ) {
					if ( poParent->tableinfo.tab21events != -1 ) {
						poParent->tableinfo.tab21open0 = (real) poParent->CloseUpdate( 21, 0, poParent->tableinfo.tab21open0, poParent->tableinfo.tab21close0, poParent->tableinfo.tab21dim_filter, poParent->tableinfo.tab21cell, poParent->tableinfo.tab21total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						Set_table_unit_tabStuntingSexRegMotherYob(poParent->table_unit_tabStuntingSexRegMotherYob+1);
						poParent->tableinfo.tab21close1 = (real) poParent->table_unit_tabStuntingSexRegMotherYob;
						poParent->tableinfo.tab21open1 = (real) poParent->CloseUpdate( 21, 1, poParent->tableinfo.tab21open1, poParent->tableinfo.tab21close1, poParent->tableinfo.tab21dim_filter, poParent->tableinfo.tab21cell, poParent->tableinfo.tab21total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab21events = poParent->__events;
					poParent->tableinfo.tab21close0 = (real) poParent->is_stunted_value_in;
					poParent->tableinfo.tab21close1 = (real) poParent->table_unit_tabStuntingSexRegMotherYob;
					poParent->tableinfo.tab21cell = (index) poParent->TableCellIndex( 21, 4,
						poParent->sex,
						OldValue,
						poParent->region_nat,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab21total_cell = (index) poParent->tableinfo.tab21cell;
					}
					else {
						poParent->tableinfo.tab21total_cell = (index) poParent->TableTotalCellIndex( 21, 4,
						poParent->sex,
						OldValue,
						poParent->region_nat,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab21dim_filter = poParent->table_filter_tabStuntingSexRegMotherYob;
				}
			}
			if ( gplTableIndexes[22] >= 0 && ( poParent->tableinfo.tab22dim_filter || poParent->table_filter_TabPrenatCare ) ) {
				if ( poParent->tableinfo.tab22events == -1 || poParent->tableinfo.tab22events < poParent->__events ) {
					if ( poParent->tableinfo.tab22events != -1 ) {
						poParent->tableinfo.tab22open0 = (real) poParent->CloseUpdate( 22, 0, poParent->tableinfo.tab22open0, poParent->tableinfo.tab22close0, poParent->tableinfo.tab22dim_filter, poParent->tableinfo.tab22cell, poParent->tableinfo.tab22total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab22open1 = (real) poParent->CloseUpdate( 22, 1, poParent->tableinfo.tab22open1, poParent->tableinfo.tab22close1, poParent->tableinfo.tab22dim_filter, poParent->tableinfo.tab22cell, poParent->tableinfo.tab22total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab22events = poParent->__events;
					poParent->tableinfo.tab22close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab22close1 = (real) poDerivedStates->Get_duration_got_prenat_care_GPC_YES();
					poParent->tableinfo.tab22cell = (index) poParent->TableCellIndex( 22, 3,
						poParent->region_nat,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab22total_cell = (index) poParent->tableinfo.tab22cell;
					}
					else {
						poParent->tableinfo.tab22total_cell = (index) poParent->TableTotalCellIndex( 22, 3,
						poParent->region_nat,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab22dim_filter = poParent->table_filter_TabPrenatCare;
				}
			}
			if ( gplTableIndexes[23] >= 0 && ( poParent->tableinfo.tab23dim_filter || poParent->table_filter_TabImmunizationChildren ) ) {
				if ( poParent->tableinfo.tab23events == -1 || poParent->tableinfo.tab23events < poParent->__events ) {
					if ( poParent->tableinfo.tab23events != -1 ) {
						poParent->tableinfo.tab23open0 = (real) poParent->CloseUpdate( 23, 0, poParent->tableinfo.tab23open0, poParent->tableinfo.tab23close0, poParent->tableinfo.tab23dim_filter, poParent->tableinfo.tab23cell, poParent->tableinfo.tab23total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab23open1 = (real) poParent->CloseUpdate( 23, 1, poParent->tableinfo.tab23open1, poParent->tableinfo.tab23close1, poParent->tableinfo.tab23dim_filter, poParent->tableinfo.tab23cell, poParent->tableinfo.tab23total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab23events = poParent->__events;
					poParent->tableinfo.tab23close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab23close1 = (real) poDerivedStates->Get_duration_is_immunized_TRUE();
					poParent->tableinfo.tab23cell = (index) poParent->TableCellIndex( 23, 6,
						poParent->sex,
						poParent->ethnicity_short,
						poParent->got_prenat_care,
						poParent->region_nat,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab23total_cell = (index) poParent->tableinfo.tab23cell;
					}
					else {
						poParent->tableinfo.tab23total_cell = (index) poParent->TableTotalCellIndex( 23, 6,
						poParent->sex,
						poParent->ethnicity_short,
						poParent->got_prenat_care,
						poParent->region_nat,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab23dim_filter = poParent->table_filter_TabImmunizationChildren;
				}
			}
		}
#line 79 "../code/MortalityDetailed.mpp"
		Set_child_mortality_group( ( poParent->mother_age_at_birth>= 17&& poParent->educ_mother== EOL_HIGH )? CMG_00: ( poParent->mother_age_at_birth>= 17&& poParent->educ_mother== EOL_MEDIUM )? CMG_01: ( poParent->mother_age_at_birth>= 17&& poParent->educ_mother== EOL_LOW )? CMG_02: ( poParent->mother_age_at_birth>= 15&& poParent->educ_mother== EOL_HIGH )? CMG_03: ( poParent->mother_age_at_birth>= 15&& poParent->educ_mother== EOL_MEDIUM )? CMG_04: ( poParent->mother_age_at_birth>= 15&& poParent->educ_mother== EOL_LOW )? CMG_05: ( poParent->educ_mother== EOL_HIGH )? CMG_06: ( poParent->educ_mother== EOL_MEDIUM )? CMG_07: CMG_08 );
#line 5865 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 100 "../code/EducationPrimaryTransmission.mpp"
		Set_educ1_group( ( poParent->educ_mother== EOL_LOW )? E1G_00: ( poParent->educ_mother== EOL_MEDIUM )? E1G_01: E1G_02 );
#line 5868 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->educ_mother = value;
}

EDUC_ONE_ENTRY_AGE PersonStateFunctions::Set_educ_one_entry_age( EDUC_ONE_ENTRY_AGE value )
{

	CHKLMT(value, educ_one_entry_age, r_min_EDUC_ONE_ENTRY_AGE, r_max_EDUC_ONE_ENTRY_AGE);
	if ( value < r_min_EDUC_ONE_ENTRY_AGE ) {
		value = r_min_EDUC_ONE_ENTRY_AGE;
	}
	else if ( value > r_max_EDUC_ONE_ENTRY_AGE ) {
		value = r_max_EDUC_ONE_ENTRY_AGE;
	}
	if ( poParent->educ_one_entry_age != value )
	{


		EDUC_ONE_ENTRY_AGE	OldValue;
		OldValue = poParent->educ_one_entry_age;
		poParent->educ_one_entry_age = value;

	}
	return poParent->educ_one_entry_age = value;
}

EDUC_ONE_LEVEL PersonStateFunctions::Set_educ_one_fate( EDUC_ONE_LEVEL value )
{

	CHKLMT(value, educ_one_fate, 0, rc_size_EDUC_ONE_LEVEL - 1);
	if ( value < 0 ) {
		value = (EDUC_ONE_LEVEL) 0;
	}
	else if ( value >= rc_size_EDUC_ONE_LEVEL ) {
		value = (EDUC_ONE_LEVEL) ( rc_size_EDUC_ONE_LEVEL - 1 );
	}
	if ( poParent->educ_one_fate != value )
	{
#line 25 "../code/EducationPrimaryTracking.mpp"
		Set_actorset_filter_asPotentialSchoolOneStudent( poParent->is_alive&& poParent->educ_one_status!= EOS_OUT&& value!= EOL_LOW );
#line 5909 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		if ( poParent->educ_one_fate == EOL_HIGH ) {
			poDerivedStates->Set_duration_educ_one_fate_EOL_HIGH( poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH() );
		}
		else if ( value == EOL_HIGH ) {
			poDerivedStates->duration_educ_one_fate_EOL_HIGH_time = poParent->__time;
		}
		if ( poParent->educ_one_fate == EOL_LOW ) {
			poDerivedStates->Set_duration_educ_one_fate_EOL_LOW( poDerivedStates->Get_duration_educ_one_fate_EOL_LOW() );
		}
		else if ( value == EOL_LOW ) {
			poDerivedStates->duration_educ_one_fate_EOL_LOW_time = poParent->__time;
		}
		if ( poParent->educ_one_fate == EOL_MEDIUM ) {
			poDerivedStates->Set_duration_educ_one_fate_EOL_MEDIUM( poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM() );
		}
		else if ( value == EOL_MEDIUM ) {
			poDerivedStates->duration_educ_one_fate_EOL_MEDIUM_time = poParent->__time;
		}


		EDUC_ONE_LEVEL	OldValue;
		OldValue = poParent->educ_one_fate;
		poParent->educ_one_fate = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[3] >= 0 && ( poParent->tableinfo.tab3dim_filter || poParent->table_filter_TabEducFateByGroup ) ) {
				if ( poParent->tableinfo.tab3events == -1 || poParent->tableinfo.tab3events < poParent->__events ) {
					if ( poParent->tableinfo.tab3events != -1 ) {
						Set_table_unit_TabEducFateByGroup(poParent->table_unit_TabEducFateByGroup+1);
						poParent->tableinfo.tab3close0 = (real) poParent->table_unit_TabEducFateByGroup;
						poParent->tableinfo.tab3open0 = (real) poParent->CloseUpdate( 3, 0, poParent->tableinfo.tab3open0, poParent->tableinfo.tab3close0, poParent->tableinfo.tab3dim_filter, poParent->tableinfo.tab3cell, poParent->tableinfo.tab3total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab3events = poParent->__events;
					poParent->tableinfo.tab3close0 = (real) poParent->table_unit_TabEducFateByGroup;
					poParent->tableinfo.tab3cell = (index) poParent->TableCellIndex( 3, 3,
						poParent->educ1_group,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE,
						OldValue );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab3total_cell = (index) poParent->tableinfo.tab3cell;
					}
					else {
						poParent->tableinfo.tab3total_cell = (index) poParent->TableTotalCellIndex( 3, 3,
						poParent->educ1_group,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE,
						OldValue );
					}
					poParent->tableinfo.tab3dim_filter = poParent->table_filter_TabEducFateByGroup;
				}
			}
			if ( gplTableIndexes[6] >= 0 && ( poParent->tableinfo.tab6dim_filter || poParent->table_filter_TabEducFateDistrYob ) ) {
				if ( poParent->tableinfo.tab6events == -1 || poParent->tableinfo.tab6events < poParent->__events ) {
					if ( poParent->tableinfo.tab6events != -1 ) {
						Set_table_unit_TabEducFateDistrYob(poParent->table_unit_TabEducFateDistrYob+1);
						poParent->tableinfo.tab6close0 = (real) poParent->table_unit_TabEducFateDistrYob;
						poParent->tableinfo.tab6open0 = (real) poParent->CloseUpdate( 6, 0, poParent->tableinfo.tab6open0, poParent->tableinfo.tab6close0, poParent->tableinfo.tab6dim_filter, poParent->tableinfo.tab6cell, poParent->tableinfo.tab6total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab6events = poParent->__events;
					poParent->tableinfo.tab6close0 = (real) poParent->table_unit_TabEducFateDistrYob;
					poParent->tableinfo.tab6cell = (index) poParent->TableCellIndex( 6, 3,
						poDerivedStates->aggregate_region_birth_REGION_NAT,
						poParent->year_of_birth - r_min_ALL_YEAR_RANGE,
						OldValue );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab6total_cell = (index) poParent->tableinfo.tab6cell;
					}
					else {
						poParent->tableinfo.tab6total_cell = (index) poParent->TableTotalCellIndex( 6, 3,
						poDerivedStates->aggregate_region_birth_REGION_NAT,
						poParent->year_of_birth - r_min_ALL_YEAR_RANGE,
						OldValue );
					}
					poParent->tableinfo.tab6dim_filter = poParent->table_filter_TabEducFateDistrYob;
				}
			}
		}
#line 22 "../code/TablesEducation.mpp"
		Set_is_aded0( ( !( poParent->is_child )&& !( poParent->is_old )&& poParent->educ_one_fate== EOL_LOW ) );
#line 5988 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 23 "../code/TablesEducation.mpp"
		Set_is_aded1( ( !( poParent->is_child )&& !( poParent->is_old )&& poParent->educ_one_fate== EOL_MEDIUM ) );
#line 5991 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 24 "../code/TablesEducation.mpp"
		Set_is_aded2( ( !( poParent->is_child )&& !( poParent->is_old )&& poParent->educ_one_fate== EOL_HIGH ) );
#line 5994 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 88 "../code/FamilyFirstUnion.mpp"
		Set_union1_group( ( poParent->educ_one_fate== EOL_LOW )? U1G_00: ( poParent->educ_one_fate== EOL_MEDIUM )? U1G_01: U1G_02 );
#line 5997 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		if ( poParent->lSpouse != NULL ) {
			if ( poParent->lSpouse->bUpdates ) {
				poParent->lSpouse->UpdateTime();
				poParent->lSpouse->poStateFunctions->Set_lSpouse_educ_one_fate( poParent->educ_one_fate );
			}
		}
		if ( poParent->poFertilityDetailedEvent->cRecalc == 0 ) {
			poParent->poFertilityDetailedEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFertilityDetailedEvent );
		};
	}
	return poParent->educ_one_fate = value;
}

EDUC_ONE_GRADE0 PersonStateFunctions::Set_educ_one_grade_attended( EDUC_ONE_GRADE0 value )
{

	CHKLMT(value, educ_one_grade_attended, r_min_EDUC_ONE_GRADE0, r_max_EDUC_ONE_GRADE0);
	if ( value < r_min_EDUC_ONE_GRADE0 ) {
		value = r_min_EDUC_ONE_GRADE0;
	}
	else if ( value > r_max_EDUC_ONE_GRADE0 ) {
		value = r_max_EDUC_ONE_GRADE0;
	}
	if ( poParent->educ_one_grade_attended != value )
	{
		if ( value == 1 ) {
			poDerivedStates->Set_entrances_educ_one_grade_attended_1( poDerivedStates->entrances_educ_one_grade_attended_1 + 1 );
		}


		EDUC_ONE_GRADE0	OldValue;
		OldValue = poParent->educ_one_grade_attended;
		poParent->educ_one_grade_attended = value;

#line 107 "../code/EducationPrimaryTracking.mpp"
		Set_educ_one_grade( COERCE(EDUC_ONE_GRADE,poParent->educ_one_grade_attended) );
#line 6035 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 112 "../code/EducationPrimaryTracking.mpp"
		Set_educ_one_level( ( poParent->educ_one_grade_passed== MAX(EDUC_ONE_GRADE) )? EOL_HIGH: ( poParent->educ_one_grade_attended== 0 )? EOL_LOW: EOL_MEDIUM );
#line 6038 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 104 "../code/EducationPrimaryTracking.mpp"
		Set_is_educ_one_student( ( poParent->educ_one_status== EOS_ATTEND|| ( poParent->educ_one_status== EOS_WAIT&& poParent->educ_one_grade_attended> 0 ) )? TRUE: FALSE );
#line 6041 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		if ( poParent->poSetPreschoolYearsEvent->cRecalc == 0 ) {
			poParent->poSetPreschoolYearsEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poSetPreschoolYearsEvent );
		};
	}
	return poParent->educ_one_grade_attended = value;
}

EDUC_ONE_GRADE PersonStateFunctions::Set_educ_one_grade_fate( EDUC_ONE_GRADE value )
{

	CHKLMT(value, educ_one_grade_fate, r_min_EDUC_ONE_GRADE, r_max_EDUC_ONE_GRADE);
	if ( value < r_min_EDUC_ONE_GRADE ) {
		value = r_min_EDUC_ONE_GRADE;
	}
	else if ( value > r_max_EDUC_ONE_GRADE ) {
		value = r_max_EDUC_ONE_GRADE;
	}
	if ( poParent->educ_one_grade_fate != value )
	{


		EDUC_ONE_GRADE	OldValue;
		OldValue = poParent->educ_one_grade_fate;
		poParent->educ_one_grade_fate = value;

	}
	return poParent->educ_one_grade_fate = value;
}

EDUC_ONE_GRADE0 PersonStateFunctions::Set_educ_one_grade_passed( EDUC_ONE_GRADE0 value )
{

	CHKLMT(value, educ_one_grade_passed, r_min_EDUC_ONE_GRADE0, r_max_EDUC_ONE_GRADE0);
	if ( value < r_min_EDUC_ONE_GRADE0 ) {
		value = r_min_EDUC_ONE_GRADE0;
	}
	else if ( value > r_max_EDUC_ONE_GRADE0 ) {
		value = r_max_EDUC_ONE_GRADE0;
	}
	if ( poParent->educ_one_grade_passed != value )
	{
#line 26 "../code/EducationSecondaryBase.mpp"
		Set_actorset_filter_asPotentialSchoolTwoStudent( value== MAX(EDUC_ONE_GRADE)&& poParent->educ_two_status!= ETS_OUT&& poParent->educ_two_delay<= gprParam->Educ2AllowedDelays );
#line 6086 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		EDUC_ONE_GRADE0	OldValue;
		OldValue = poParent->educ_one_grade_passed;
		poParent->educ_one_grade_passed = value;

#line 112 "../code/EducationPrimaryTracking.mpp"
		Set_educ_one_level( ( poParent->educ_one_grade_passed== MAX(EDUC_ONE_GRADE) )? EOL_HIGH: ( poParent->educ_one_grade_attended== 0 )? EOL_LOW: EOL_MEDIUM );
#line 6095 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 25 "../code/TablesEducation.mpp"
		Set_has_passed_primary( ( poParent->educ_one_grade_passed== MAX(EDUC_ONE_GRADE) ) );
#line 6098 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 85 "../code/HumanCapitalIndex.mpp"
		poParent->years_of_schooling = poParent->years_preschool+ poParent->educ_one_grade_passed+ poParent->educ_two_grade_passed;
#line 6101 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->educ_one_grade_passed = value;
}

EDUC_ONE_STATUS PersonStateFunctions::Set_educ_one_status( EDUC_ONE_STATUS value )
{

	CHKLMT(value, educ_one_status, 0, rc_size_EDUC_ONE_STATUS - 1);
	if ( value < 0 ) {
		value = (EDUC_ONE_STATUS) 0;
	}
	else if ( value >= rc_size_EDUC_ONE_STATUS ) {
		value = (EDUC_ONE_STATUS) ( rc_size_EDUC_ONE_STATUS - 1 );
	}
	if ( poParent->educ_one_status != value )
	{
#line 25 "../code/EducationPrimaryTracking.mpp"
		Set_actorset_filter_asPotentialSchoolOneStudent( poParent->is_alive&& value!= EOS_OUT&& poParent->educ_one_fate!= EOL_LOW );
#line 6120 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		EDUC_ONE_STATUS	OldValue;
		OldValue = poParent->educ_one_status;
		poParent->educ_one_status = value;

#line 104 "../code/EducationPrimaryTracking.mpp"
		Set_is_educ_one_student( ( poParent->educ_one_status== EOS_ATTEND|| ( poParent->educ_one_status== EOS_WAIT&& poParent->educ_one_grade_attended> 0 ) )? TRUE: FALSE );
#line 6129 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/TablesEducation.mpp"
		Set_tab_out_of_primary( ( poParent->educ_one_status!= EOS_WAIT&& poParent->educ_one_status!= EOS_ATTEND&& !( poParent->has_passed_primary ) ) );
#line 6132 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->educ_one_status = value;
}

logical PersonStateFunctions::Set_educ_one_to_process( logical value )
{

	value = value != 0;
	if ( poParent->educ_one_to_process != value )
	{
#line 28 "../code/EducationPrimaryTracking.mpp"
		Set_actorset_filter_asPotentialSchoolOneStudentToProcess( poParent->is_alive&& value );
#line 6145 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		logical	OldValue;
		OldValue = poParent->educ_one_to_process;
		poParent->educ_one_to_process = value;

	}
	return poParent->educ_one_to_process = value;
}

int PersonStateFunctions::Set_educ_two_delay( int value )
{

	if ( poParent->educ_two_delay != value )
	{
#line 26 "../code/EducationSecondaryBase.mpp"
		Set_actorset_filter_asPotentialSchoolTwoStudent( poParent->educ_one_grade_passed== MAX(EDUC_ONE_GRADE)&& poParent->educ_two_status!= ETS_OUT&& value<= gprParam->Educ2AllowedDelays );
#line 6163 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		int	OldValue;
		OldValue = poParent->educ_two_delay;
		poParent->educ_two_delay = value;

	}
	return poParent->educ_two_delay = value;
}

EDUC_TWO_GRADE0 PersonStateFunctions::Set_educ_two_grade_attended( EDUC_TWO_GRADE0 value )
{

	CHKLMT(value, educ_two_grade_attended, r_min_EDUC_TWO_GRADE0, r_max_EDUC_TWO_GRADE0);
	if ( value < r_min_EDUC_TWO_GRADE0 ) {
		value = r_min_EDUC_TWO_GRADE0;
	}
	else if ( value > r_max_EDUC_TWO_GRADE0 ) {
		value = r_max_EDUC_TWO_GRADE0;
	}
	if ( poParent->educ_two_grade_attended != value )
	{


		EDUC_TWO_GRADE0	OldValue;
		OldValue = poParent->educ_two_grade_attended;
		poParent->educ_two_grade_attended = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[12] >= 0 && ( poParent->tableinfo.tab12dim_filter || poParent->table_filter_TabSchool2TrackTab ) ) {
				if ( poParent->tableinfo.tab12events == -1 || poParent->tableinfo.tab12events < poParent->__events ) {
					if ( poParent->tableinfo.tab12events != -1 ) {
						Set_table_unit_TabSchool2TrackTab(poParent->table_unit_TabSchool2TrackTab+1);
						poParent->tableinfo.tab12close0 = (real) poParent->table_unit_TabSchool2TrackTab;
						poParent->tableinfo.tab12open0 = (real) poParent->CloseUpdate( 12, 0, poParent->tableinfo.tab12open0, poParent->tableinfo.tab12close0, poParent->tableinfo.tab12dim_filter, poParent->tableinfo.tab12cell, poParent->tableinfo.tab12total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab12events = poParent->__events;
					poParent->tableinfo.tab12close0 = (real) poParent->table_unit_TabSchool2TrackTab;
					poParent->tableinfo.tab12cell = (index) poParent->TableCellIndex( 12, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						OldValue - r_min_EDUC_TWO_GRADE0 );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab12total_cell = (index) poParent->tableinfo.tab12cell;
					}
					else {
						poParent->tableinfo.tab12total_cell = (index) poParent->TableTotalCellIndex( 12, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						OldValue - r_min_EDUC_TWO_GRADE0 );
					}
					poParent->tableinfo.tab12dim_filter = poParent->table_filter_TabSchool2TrackTab;
				}
			}
		}
	}
	return poParent->educ_two_grade_attended = value;
}

EDUC_TWO_GRADE0 PersonStateFunctions::Set_educ_two_grade_passed( EDUC_TWO_GRADE0 value )
{

	CHKLMT(value, educ_two_grade_passed, r_min_EDUC_TWO_GRADE0, r_max_EDUC_TWO_GRADE0);
	if ( value < r_min_EDUC_TWO_GRADE0 ) {
		value = r_min_EDUC_TWO_GRADE0;
	}
	else if ( value > r_max_EDUC_TWO_GRADE0 ) {
		value = r_max_EDUC_TWO_GRADE0;
	}
	if ( poParent->educ_two_grade_passed != value )
	{


		EDUC_TWO_GRADE0	OldValue;
		OldValue = poParent->educ_two_grade_passed;
		poParent->educ_two_grade_passed = value;

#line 113 "../code/EducationSecondaryBase.mpp"
		Set_educ_two_level( ( poParent->educ_two_status== ETS_NEVER )? ETL_NEVER: ( poParent->educ_two_grade_passed== MAX(EDUC_TWO_GRADE) )? ETL_INCOMPLETE: ETL_GRADUATE );
#line 6241 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 85 "../code/HumanCapitalIndex.mpp"
		poParent->years_of_schooling = poParent->years_preschool+ poParent->educ_one_grade_passed+ poParent->educ_two_grade_passed;
#line 6244 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->educ_two_grade_passed = value;
}

EDUC_TWO_STATUS PersonStateFunctions::Set_educ_two_status( EDUC_TWO_STATUS value )
{

	CHKLMT(value, educ_two_status, 0, rc_size_EDUC_TWO_STATUS - 1);
	if ( value < 0 ) {
		value = (EDUC_TWO_STATUS) 0;
	}
	else if ( value >= rc_size_EDUC_TWO_STATUS ) {
		value = (EDUC_TWO_STATUS) ( rc_size_EDUC_TWO_STATUS - 1 );
	}
	if ( poParent->educ_two_status != value )
	{
#line 26 "../code/EducationSecondaryBase.mpp"
		Set_actorset_filter_asPotentialSchoolTwoStudent( poParent->educ_one_grade_passed== MAX(EDUC_ONE_GRADE)&& value!= ETS_OUT&& poParent->educ_two_delay<= gprParam->Educ2AllowedDelays );
#line 6263 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 221 "../code/TablesEducation.mpp"
		Set_table_filter_TabSchool2TrackTab( poDerivedStates->trigger_changes_calendar_year&& poParent->in_projected_time&& poParent->is_resident&& value== ETS_ATTEND );
#line 6266 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		EDUC_TWO_STATUS	OldValue;
		OldValue = poParent->educ_two_status;
		poParent->educ_two_status = value;

#line 113 "../code/EducationSecondaryBase.mpp"
		Set_educ_two_level( ( poParent->educ_two_status== ETS_NEVER )? ETL_NEVER: ( poParent->educ_two_grade_passed== MAX(EDUC_TWO_GRADE) )? ETL_INCOMPLETE: ETL_GRADUATE );
#line 6275 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->educ_two_status = value;
}

logical PersonStateFunctions::Set_educ_two_to_process( logical value )
{

	value = value != 0;
	if ( poParent->educ_two_to_process != value )
	{
#line 29 "../code/EducationSecondaryBase.mpp"
		Set_actorset_filter_asPotentialSchoolTwoStudentToProcess( poParent->is_alive&& value );
#line 6288 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		logical	OldValue;
		OldValue = poParent->educ_two_to_process;
		poParent->educ_two_to_process = value;

	}
	return poParent->educ_two_to_process = value;
}

ETHNICITY PersonStateFunctions::Set_ethnicity( ETHNICITY value )
{

	CHKLMT(value, ethnicity, 0, rc_size_ETHNICITY - 1);
	if ( value < 0 ) {
		value = (ETHNICITY) 0;
	}
	else if ( value >= rc_size_ETHNICITY ) {
		value = (ETHNICITY) ( rc_size_ETHNICITY - 1 );
	}
	if ( poParent->ethnicity != value )
	{
		poDerivedStates->Set_aggregate_ethnicity_ETHNICITY_SHORT( ETHNICITY_To_ETHNICITY_SHORT( value ) );


		ETHNICITY	OldValue;
		OldValue = poParent->ethnicity;
		poParent->ethnicity = value;

	}
	return poParent->ethnicity = value;
}

logical PersonStateFunctions::Set_ever_resident( logical value )
{

	value = value != 0;
	if ( poParent->ever_resident != value )
	{
#line 10 "../code/ImmigrationFromScratch.mpp"
		Set_actorset_filter_asPotentialImmigrantMothers( poParent->is_alive&& poParent->sex== FEMALE&& poParent->creation_type== CT_SCRATCH&& !( value )&& poParent->family_role!= FR_CHILD&& poParent->ever_union );
#line 6330 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		logical	OldValue;
		OldValue = poParent->ever_resident;
		poParent->ever_resident = value;

#line 47 "../code/FertilityGeneral.mpp"
		Set_is_fertile( ( poParent->sex== FEMALE&& WITHIN(FERTILE_AGE_RANGE,poParent->integer_age)&& ( poParent->creation_type== CT_BIRTH|| ( poParent->creation_type== CT_START&& poParent->calendar_year>= MIN(SIM_YEAR_RANGE) )|| ( poParent->creation_type== CT_POOL&& poParent->ever_resident )|| ( poParent->creation_type== CT_SCRATCH&& poParent->ever_resident ) ) )? TRUE: FALSE );
#line 6339 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 45 "../code/MortalityGeneral.mpp"
		Set_is_mortal( ( poParent->creation_type== CT_BIRTH|| ( poParent->creation_type== CT_START&& poParent->calendar_year>= MIN(SIM_YEAR_RANGE) )|| ( poParent->creation_type== CT_POOL&& poParent->ever_resident )|| ( poParent->creation_type== CT_SCRATCH&& poParent->ever_resident ) )? TRUE: FALSE );
#line 6342 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		if ( poParent->poBackMigrationEvent->cRecalc == 0 ) {
			poParent->poBackMigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poBackMigrationEvent );
		};
		if ( poParent->poFirstImmigrationFromPoolEvent->cRecalc == 0 ) {
			poParent->poFirstImmigrationFromPoolEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstImmigrationFromPoolEvent );
		};
		if ( poParent->poFirstImmigrationFromScratchEvent->cRecalc == 0 ) {
			poParent->poFirstImmigrationFromScratchEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstImmigrationFromScratchEvent );
		};
	}
	return poParent->ever_resident = value;
}

logical PersonStateFunctions::Set_ever_union( logical value )
{

	value = value != 0;
	if ( poParent->ever_union != value )
	{
#line 10 "../code/ImmigrationFromScratch.mpp"
		Set_actorset_filter_asPotentialImmigrantMothers( poParent->is_alive&& poParent->sex== FEMALE&& poParent->creation_type== CT_SCRATCH&& !( poParent->ever_resident )&& poParent->family_role!= FR_CHILD&& value );
#line 6367 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		if ( poParent->ever_union == TRUE ) {
			poDerivedStates->Set_duration_ever_union_TRUE( poDerivedStates->Get_duration_ever_union_TRUE() );
		}
		else if ( value == TRUE ) {
			poDerivedStates->duration_ever_union_TRUE_time = poParent->__time;
		}


		logical	OldValue;
		OldValue = poParent->ever_union;
		poParent->ever_union = value;

		if ( poParent->poFertilityDetailedEvent->cRecalc == 0 ) {
			poParent->poFertilityDetailedEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFertilityDetailedEvent );
		};
		if ( poParent->poFirstUnionFormationEvent->cRecalc == 0 ) {
			poParent->poFirstUnionFormationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstUnionFormationEvent );
		};
	}
	return poParent->ever_union = value;
}

GEO PersonStateFunctions::Set_geo( GEO value )
{

	CHKLMT(value, geo, 0, rc_size_GEO - 1);
	if ( value < 0 ) {
		value = (GEO) 0;
	}
	else if ( value >= rc_size_GEO ) {
		value = (GEO) ( rc_size_GEO - 1 );
	}
	if ( poParent->geo != value )
	{
		bool bScheduleTimeChange = false;	//indicates that schedule event should be recalculated
		if (	(poParent->schedule_event_time == poDerivedStates->trigger_changes_geo_time))
		{
			//force timeScheduleEvent to calculate min time
			poParent->schedule_event_time = TIME_INFINITE;
		//indicate that schedule event needs to be recalculated
			bScheduleTimeChange = true;
		}
		poDerivedStates->Set_aggregate_geo_GEO_NAT( GEO_To_GEO_NAT( value ) );
		poDerivedStates->Set_aggregate_geo_REGION_INT( GEO_To_REGION_INT( value ) );


		GEO	OldValue;
		OldValue = poParent->geo;
		poParent->geo = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[0] >= 0 && ( poParent->tableinfo.tab0dim_filter || poParent->table_filter_tabPreSchool ) ) {
				if ( poParent->tableinfo.tab0events == -1 || poParent->tableinfo.tab0events < poParent->__events ) {
					if ( poParent->tableinfo.tab0events != -1 ) {
						poParent->tableinfo.tab0open0 = (real) poParent->CloseUpdate( 0, 0, poParent->tableinfo.tab0open0, poParent->tableinfo.tab0close0, poParent->tableinfo.tab0dim_filter, poParent->tableinfo.tab0cell, poParent->tableinfo.tab0total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab0open1 = (real) poParent->CloseUpdate( 0, 1, poParent->tableinfo.tab0open1, poParent->tableinfo.tab0close1, poParent->tableinfo.tab0dim_filter, poParent->tableinfo.tab0cell, poParent->tableinfo.tab0total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab0open2 = (real) poParent->CloseUpdate( 0, 2, poParent->tableinfo.tab0open2, poParent->tableinfo.tab0close2, poParent->tableinfo.tab0dim_filter, poParent->tableinfo.tab0cell, poParent->tableinfo.tab0total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab0open3 = (real) poParent->CloseUpdate( 0, 3, poParent->tableinfo.tab0open3, poParent->tableinfo.tab0close3, poParent->tableinfo.tab0dim_filter, poParent->tableinfo.tab0cell, poParent->tableinfo.tab0total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab0events = poParent->__events;
					poParent->tableinfo.tab0close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab0close1 = (real) poDerivedStates->Get_duration_years_preschool_0();
					poParent->tableinfo.tab0close2 = (real) poDerivedStates->Get_duration_years_preschool_1();
					poParent->tableinfo.tab0close3 = (real) poDerivedStates->Get_duration_years_preschool_2();
					poParent->tableinfo.tab0cell = (index) poParent->TableCellIndex( 0, 2,
						OldValue,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab0total_cell = (index) poParent->tableinfo.tab0cell;
					}
					else {
						poParent->tableinfo.tab0total_cell = (index) poParent->TableTotalCellIndex( 0, 2,
						OldValue,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab0dim_filter = poParent->table_filter_tabPreSchool;
				}
			}
			if ( gplTableIndexes[17] >= 0 && ( poParent->tableinfo.tab17dim_filter || poParent->table_filter_tabPopulationYearPlace ) ) {
				if ( poParent->tableinfo.tab17events == -1 || poParent->tableinfo.tab17events < poParent->__events ) {
					if ( poParent->tableinfo.tab17events != -1 ) {
						poParent->tableinfo.tab17open0 = (real) poParent->CloseUpdate( 17, 0, poParent->tableinfo.tab17open0, poParent->tableinfo.tab17close0, poParent->tableinfo.tab17dim_filter, poParent->tableinfo.tab17cell, poParent->tableinfo.tab17total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab17events = poParent->__events;
					poParent->tableinfo.tab17close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab17cell = (index) poParent->TableCellIndex( 17, 2,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab17total_cell = (index) poParent->tableinfo.tab17cell;
					}
					else {
						poParent->tableinfo.tab17total_cell = (index) poParent->TableTotalCellIndex( 17, 2,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab17dim_filter = poParent->table_filter_tabPopulationYearPlace;
				}
			}
			if ( gplTableIndexes[18] >= 0 && ( poParent->tableinfo.tab18dim_filter || poParent->table_filter_tabBirthsYearPlace ) ) {
				if ( poParent->tableinfo.tab18events == -1 || poParent->tableinfo.tab18events < poParent->__events ) {
					if ( poParent->tableinfo.tab18events != -1 ) {
						poParent->tableinfo.tab18open0 = (real) poParent->CloseUpdate( 18, 0, poParent->tableinfo.tab18open0, poParent->tableinfo.tab18close0, poParent->tableinfo.tab18dim_filter, poParent->tableinfo.tab18cell, poParent->tableinfo.tab18total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab18events = poParent->__events;
					poParent->tableinfo.tab18close0 = (real) poDerivedStates->transitions_is_alive_FALSE_TRUE;
					poParent->tableinfo.tab18cell = (index) poParent->TableCellIndex( 18, 2,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab18total_cell = (index) poParent->tableinfo.tab18cell;
					}
					else {
						poParent->tableinfo.tab18total_cell = (index) poParent->TableTotalCellIndex( 18, 2,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab18dim_filter = poParent->table_filter_tabBirthsYearPlace;
				}
			}
			if ( gplTableIndexes[19] >= 0 && ( poParent->tableinfo.tab19dim_filter || poParent->table_filter_tabMigrationOriginDestination ) ) {
				if ( poParent->tableinfo.tab19events == -1 || poParent->tableinfo.tab19events < poParent->__events ) {
					if ( poParent->tableinfo.tab19events != -1 ) {
						poParent->tableinfo.tab19open0 = (real) poParent->CloseUpdate( 19, 0, poParent->tableinfo.tab19open0, poParent->tableinfo.tab19close0, poParent->tableinfo.tab19dim_filter, poParent->tableinfo.tab19cell, poParent->tableinfo.tab19total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab19events = poParent->__events;
					poParent->tableinfo.tab19close0 = (real) poParent->is_alive_value_in;
					poParent->tableinfo.tab19cell = (index) poParent->TableCellIndex( 19, 3,
						poParent->selected_year10,
						poParent->geo_prev,
						OldValue );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab19total_cell = (index) poParent->tableinfo.tab19cell;
					}
					else {
						poParent->tableinfo.tab19total_cell = (index) poParent->TableTotalCellIndex( 19, 3,
						poParent->selected_year10,
						poParent->geo_prev,
						OldValue );
					}
					poParent->tableinfo.tab19dim_filter = poParent->table_filter_tabMigrationOriginDestination;
				}
			}
		}
#line 79 "../code/ActorPerson.mpp"
		Set_is_resident( ( poParent->geo< SIZE(GEO_NAT) ) );
#line 6516 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		poDerivedStates->Set_trigger_changes_geo( TRUE );
		poDerivedStates->trigger_changes_geo_time = poParent->__time;
		if ( poParent->schedule_event_time > poDerivedStates->trigger_changes_geo_time)
		{
			poParent->schedule_event_time = poDerivedStates->trigger_changes_geo_time;
			//indicate that schedule event should be recalculated
			bScheduleTimeChange = true;
		}
		// actorset dimension or order state
		if (poParent->actorset_filter_asAvailableMale)
		{
			asAvailableMale[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][OldValue]->Remove(poParent->poasAvailableMale);
			poParent->poasAvailableMale = asAvailableMale[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->geo]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asAvailableMaleByType)
		{
			asAvailableMaleByType[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->spouse_group][OldValue]->Remove(poParent->poasAvailableMaleByType);
			poParent->poasAvailableMaleByType = asAvailableMaleByType[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->spouse_group][poParent->geo]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asAvailableMaleByTypeImmiScratch)
		{
			asAvailableMaleByTypeImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->spouse_group][OldValue][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->Remove(poParent->poasAvailableMaleByTypeImmiScratch);
			poParent->poasAvailableMaleByTypeImmiScratch = asAvailableMaleByTypeImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->spouse_group][poParent->geo][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asAvailableMaleImmiScratch)
		{
			asAvailableMaleImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][OldValue][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->Remove(poParent->poasAvailableMaleImmiScratch);
			poParent->poasAvailableMaleImmiScratch = asAvailableMaleImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->geo][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asPotentialImmigrantMothers)
		{
			asPotentialImmigrantMothers[RANGE_POS(AGE_RANGE, poParent->integer_age)][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)][OldValue]->Remove(poParent->poasPotentialImmigrantMothers);
			poParent->poasPotentialImmigrantMothers = asPotentialImmigrantMothers[RANGE_POS(AGE_RANGE, poParent->integer_age)][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)][poParent->geo]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asSimBornAge0)
		{
			asSimBornAge0[poParent->sex][OldValue][poParent->educ1_group]->Remove(poParent->poasSimBornAge0);
			poParent->poasSimBornAge0 = asSimBornAge0[poParent->sex][poParent->geo][poParent->educ1_group]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asWantToMove)
		{
			asWantToMove[OldValue][poParent->geo_want_to_move][poParent->sex][RANGE_POS(AGE_MIG, poParent->age_mig)]->Remove(poParent->poasWantToMove);
			poParent->poasWantToMove = asWantToMove[poParent->geo][poParent->geo_want_to_move][poParent->sex][RANGE_POS(AGE_MIG, poParent->age_mig)]->NewNode(poParent);
		}
		if ( poParent->poEmigrationEvent->cRecalc == 0 ) {
			poParent->poEmigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poEmigrationEvent );
		};
		if ( poParent->poFertilityDetailedEvent->cRecalc == 0 ) {
			poParent->poFertilityDetailedEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFertilityDetailedEvent );
		};
		if ( poParent->poMigrationEvent->cRecalc == 0 ) {
			poParent->poMigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMigrationEvent );
		};

		if (bScheduleTimeChange && poParent->poScheduleEvent->cRecalc == 0)
		{
			poParent->poScheduleEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent(poParent->poScheduleEvent);
		}
	}
	return poParent->geo = value;
}

GEO PersonStateFunctions::Set_geo_birth( GEO value )
{

	CHKLMT(value, geo_birth, 0, rc_size_GEO - 1);
	if ( value < 0 ) {
		value = (GEO) 0;
	}
	else if ( value >= rc_size_GEO ) {
		value = (GEO) ( rc_size_GEO - 1 );
	}
	if ( poParent->geo_birth != value )
	{
		poDerivedStates->Set_aggregate_geo_birth_GEO_NAT( GEO_To_GEO_NAT( value ) );
		poDerivedStates->Set_aggregate_geo_birth_REGION_INT( GEO_To_REGION_INT( value ) );


		GEO	OldValue;
		OldValue = poParent->geo_birth;
		poParent->geo_birth = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[1] >= 0 && ( poParent->tableinfo.tab1dim_filter || poParent->table_filter_tabEducationFateGeobirYob ) ) {
				if ( poParent->tableinfo.tab1events == -1 || poParent->tableinfo.tab1events < poParent->__events ) {
					if ( poParent->tableinfo.tab1events != -1 ) {
						poParent->tableinfo.tab1open0 = (real) poParent->CloseUpdate( 1, 0, poParent->tableinfo.tab1open0, poParent->tableinfo.tab1close0, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab1open1 = (real) poParent->CloseUpdate( 1, 1, poParent->tableinfo.tab1open1, poParent->tableinfo.tab1close1, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab1open2 = (real) poParent->CloseUpdate( 1, 2, poParent->tableinfo.tab1open2, poParent->tableinfo.tab1close2, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab1open3 = (real) poParent->CloseUpdate( 1, 3, poParent->tableinfo.tab1open3, poParent->tableinfo.tab1close3, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab1events = poParent->__events;
					poParent->tableinfo.tab1close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab1close1 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH();
					poParent->tableinfo.tab1close2 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_LOW();
					poParent->tableinfo.tab1close3 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM();
					poParent->tableinfo.tab1cell = (index) poParent->TableCellIndex( 1, 3,
						poParent->sex,
						OldValue,
						poParent->year_of_birth - r_min_ALL_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab1total_cell = (index) poParent->tableinfo.tab1cell;
					}
					else {
						poParent->tableinfo.tab1total_cell = (index) poParent->TableTotalCellIndex( 1, 3,
						poParent->sex,
						OldValue,
						poParent->year_of_birth - r_min_ALL_YEAR_RANGE );
					}
					poParent->tableinfo.tab1dim_filter = poParent->table_filter_tabEducationFateGeobirYob;
				}
			}
			if ( gplTableIndexes[8] >= 0 && ( poParent->tableinfo.tab8dim_filter || poParent->table_filter_TabEduc15ByDistrictBirth ) ) {
				if ( poParent->tableinfo.tab8events == -1 || poParent->tableinfo.tab8events < poParent->__events ) {
					if ( poParent->tableinfo.tab8events != -1 ) {
						poParent->tableinfo.tab8open0 = (real) poParent->CloseUpdate( 8, 0, poParent->tableinfo.tab8open0, poParent->tableinfo.tab8close0, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab8open1 = (real) poParent->CloseUpdate( 8, 1, poParent->tableinfo.tab8open1, poParent->tableinfo.tab8close1, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab8open2 = (real) poParent->CloseUpdate( 8, 2, poParent->tableinfo.tab8open2, poParent->tableinfo.tab8close2, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab8open3 = (real) poParent->CloseUpdate( 8, 3, poParent->tableinfo.tab8open3, poParent->tableinfo.tab8close3, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab8events = poParent->__events;
					poParent->tableinfo.tab8close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab8close1 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH();
					poParent->tableinfo.tab8close2 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_LOW();
					poParent->tableinfo.tab8close3 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM();
					poParent->tableinfo.tab8cell = (index) poParent->TableCellIndex( 8, 3,
						OldValue,
						poParent->sex,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab8total_cell = (index) poParent->tableinfo.tab8cell;
					}
					else {
						poParent->tableinfo.tab8total_cell = (index) poParent->TableTotalCellIndex( 8, 3,
						OldValue,
						poParent->sex,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab8dim_filter = poParent->table_filter_TabEduc15ByDistrictBirth;
				}
			}
		}
#line 82 "../code/HumanCapitalIndex.mpp"
		Set_in_hci_sample( ( poParent->creation_type== CT_BIRTH&& poParent->geo_birth< SIZE(GEO_NAT) )? TRUE: FALSE );
#line 6665 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		// actorset dimension or order state
		if (poParent->actorset_filter_asResidentsAge0SexGeo)
		{
			asResidentsAge0SexGeo[poParent->sex][OldValue]->Remove(poParent->poasResidentsAge0SexGeo);
			poParent->poasResidentsAge0SexGeo = asResidentsAge0SexGeo[poParent->sex][poParent->geo_birth]->NewNode(poParent);
		}
	}
	return poParent->geo_birth = value;
}

GEO PersonStateFunctions::Set_geo_prev( GEO value )
{

	CHKLMT(value, geo_prev, 0, rc_size_GEO - 1);
	if ( value < 0 ) {
		value = (GEO) 0;
	}
	else if ( value >= rc_size_GEO ) {
		value = (GEO) ( rc_size_GEO - 1 );
	}
	if ( poParent->geo_prev != value )
	{


		GEO	OldValue;
		OldValue = poParent->geo_prev;
		poParent->geo_prev = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[19] >= 0 && ( poParent->tableinfo.tab19dim_filter || poParent->table_filter_tabMigrationOriginDestination ) ) {
				if ( poParent->tableinfo.tab19events == -1 || poParent->tableinfo.tab19events < poParent->__events ) {
					if ( poParent->tableinfo.tab19events != -1 ) {
						poParent->tableinfo.tab19open0 = (real) poParent->CloseUpdate( 19, 0, poParent->tableinfo.tab19open0, poParent->tableinfo.tab19close0, poParent->tableinfo.tab19dim_filter, poParent->tableinfo.tab19cell, poParent->tableinfo.tab19total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab19events = poParent->__events;
					poParent->tableinfo.tab19close0 = (real) poParent->is_alive_value_in;
					poParent->tableinfo.tab19cell = (index) poParent->TableCellIndex( 19, 3,
						poParent->selected_year10,
						OldValue,
						poParent->geo );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab19total_cell = (index) poParent->tableinfo.tab19cell;
					}
					else {
						poParent->tableinfo.tab19total_cell = (index) poParent->TableTotalCellIndex( 19, 3,
						poParent->selected_year10,
						OldValue,
						poParent->geo );
					}
					poParent->tableinfo.tab19dim_filter = poParent->table_filter_tabMigrationOriginDestination;
				}
			}
		}
	}
	return poParent->geo_prev = value;
}

GEO PersonStateFunctions::Set_geo_want_to_move( GEO value )
{

	CHKLMT(value, geo_want_to_move, 0, rc_size_GEO - 1);
	if ( value < 0 ) {
		value = (GEO) 0;
	}
	else if ( value >= rc_size_GEO ) {
		value = (GEO) ( rc_size_GEO - 1 );
	}
	if ( poParent->geo_want_to_move != value )
	{


		GEO	OldValue;
		OldValue = poParent->geo_want_to_move;
		poParent->geo_want_to_move = value;

		// actorset dimension or order state
		if (poParent->actorset_filter_asWantToMove)
		{
			asWantToMove[poParent->geo][OldValue][poParent->sex][RANGE_POS(AGE_MIG, poParent->age_mig)]->Remove(poParent->poasWantToMove);
			poParent->poasWantToMove = asWantToMove[poParent->geo][poParent->geo_want_to_move][poParent->sex][RANGE_POS(AGE_MIG, poParent->age_mig)]->NewNode(poParent);
		}
	}
	return poParent->geo_want_to_move = value;
}

GOT_PRENAT_CARE PersonStateFunctions::Set_got_prenat_care( GOT_PRENAT_CARE value )
{

	CHKLMT(value, got_prenat_care, 0, rc_size_GOT_PRENAT_CARE - 1);
	if ( value < 0 ) {
		value = (GOT_PRENAT_CARE) 0;
	}
	else if ( value >= rc_size_GOT_PRENAT_CARE ) {
		value = (GOT_PRENAT_CARE) ( rc_size_GOT_PRENAT_CARE - 1 );
	}
	if ( poParent->got_prenat_care != value )
	{
		if ( poParent->got_prenat_care == GPC_YES ) {
			poDerivedStates->Set_duration_got_prenat_care_GPC_YES( poDerivedStates->Get_duration_got_prenat_care_GPC_YES() );
		}
		else if ( value == GPC_YES ) {
			poDerivedStates->duration_got_prenat_care_GPC_YES_time = poParent->__time;
		}


		GOT_PRENAT_CARE	OldValue;
		OldValue = poParent->got_prenat_care;
		poParent->got_prenat_care = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[23] >= 0 && ( poParent->tableinfo.tab23dim_filter || poParent->table_filter_TabImmunizationChildren ) ) {
				if ( poParent->tableinfo.tab23events == -1 || poParent->tableinfo.tab23events < poParent->__events ) {
					if ( poParent->tableinfo.tab23events != -1 ) {
						poParent->tableinfo.tab23open0 = (real) poParent->CloseUpdate( 23, 0, poParent->tableinfo.tab23open0, poParent->tableinfo.tab23close0, poParent->tableinfo.tab23dim_filter, poParent->tableinfo.tab23cell, poParent->tableinfo.tab23total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab23open1 = (real) poParent->CloseUpdate( 23, 1, poParent->tableinfo.tab23open1, poParent->tableinfo.tab23close1, poParent->tableinfo.tab23dim_filter, poParent->tableinfo.tab23cell, poParent->tableinfo.tab23total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab23events = poParent->__events;
					poParent->tableinfo.tab23close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab23close1 = (real) poDerivedStates->Get_duration_is_immunized_TRUE();
					poParent->tableinfo.tab23cell = (index) poParent->TableCellIndex( 23, 6,
						poParent->sex,
						poParent->ethnicity_short,
						OldValue,
						poParent->region_nat,
						poParent->educ_mother,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab23total_cell = (index) poParent->tableinfo.tab23cell;
					}
					else {
						poParent->tableinfo.tab23total_cell = (index) poParent->TableTotalCellIndex( 23, 6,
						poParent->sex,
						poParent->ethnicity_short,
						OldValue,
						poParent->region_nat,
						poParent->educ_mother,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab23dim_filter = poParent->table_filter_TabImmunizationChildren;
				}
			}
		}
	}
	return poParent->got_prenat_care = value;
}

AGE_RANGE PersonStateFunctions::Set_integer_age( AGE_RANGE value )
{

	CHKLMT(value, integer_age, r_min_AGE_RANGE, r_max_AGE_RANGE);
	if ( value < r_min_AGE_RANGE ) {
		value = r_min_AGE_RANGE;
	}
	else if ( value > r_max_AGE_RANGE ) {
		value = r_max_AGE_RANGE;
	}
	if ( poParent->integer_age != value )
	{
#line 14 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMale( poParent->is_alive&& poParent->sex== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,value)&& poParent->is_mortal );
#line 6826 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 18 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleByType( poParent->is_alive&& poParent->sex== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,value)&& poParent->is_mortal );
#line 6829 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 32 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleByTypeImmiScratch( poParent->is_alive&& poParent->creation_type== CT_SCRATCH&& poParent->sex== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,value)&& !( poParent->is_mortal ) );
#line 6832 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleImmiScratch( poParent->is_alive&& poParent->creation_type== CT_SCRATCH&& poParent->sex== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,value)&& !( poParent->is_mortal ) );
#line 6835 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 22 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAge( poParent->is_alive&& poParent->sex== FEMALE&& poParent->has_spouse&& WITHIN(FEMALE_SPOUSE_AGE,value)&& poParent->is_mortal );
#line 6838 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 36 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch( poParent->is_alive&& poParent->creation_type== CT_SCRATCH&& poParent->sex== FEMALE&& poParent->has_spouse&& WITHIN(FEMALE_SPOUSE_AGE,value)&& !( poParent->is_mortal ) );
#line 6841 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 18 "../code/ChildVaccination.mpp"
		Set_actorset_filter_asResidentBabies( poParent->is_alive&& value== 0&& poParent->is_resident&& poParent->creation_type== CT_BIRTH );
#line 6844 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 30 "../code/EducationPrimaryBase.mpp"
		Set_actorset_filter_asResidentsAge0Sex( poParent->is_alive&& value== 0&& poParent->is_resident );
#line 6847 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 27 "../code/EducationPrimaryBase.mpp"
		Set_actorset_filter_asResidentsAge0SexGeo( poParent->is_alive&& poParent->is_resident&& value== 0 );
#line 6850 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 35 "../code/EducationPrimaryTransmission.mpp"
		Set_actorset_filter_asSimBornAge0( poParent->is_alive&& poParent->creation_type== CT_BIRTH&& value== 0 );
#line 6853 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_actorset_filter_asWomenNoChildren( poParent->is_alive&& poParent->sex== FEMALE&& WITHIN(FEMALE_SPOUSE_AGE,value)&& !( poParent->lives_with_dependent_child )&& poParent->is_mortal&& !( poParent->is_blocked_from_marriage ) );
#line 6856 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 23 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_actorset_filter_asWomenWithChildren( poParent->is_alive&& poParent->sex== FEMALE&& WITHIN(FEMALE_SPOUSE_AGE,value)&& poParent->lives_with_dependent_child&& poParent->is_mortal&& !( poParent->is_blocked_from_marriage ) );
#line 6859 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 46 "../code/TablesVaccination.mpp"
		Set_table_filter_TabChildVaccination( poParent->in_projected_time&& poParent->creation_type== CT_BIRTH&& poParent->is_resident&& value== 0 );
#line 6862 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 165 "../code/TablesEducation.mpp"
		Set_table_filter_TabEduc15ByDistrict( value== 15&& poParent->in_projected_time&& poParent->is_resident );
#line 6865 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 178 "../code/TablesEducation.mpp"
		Set_table_filter_TabEduc15ByDistrictBirth( value== 15&& poParent->in_projected_time&& poParent->is_resident );
#line 6868 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 88 "../code/TablesEducation.mpp"
		Set_table_filter_tabEducationFateGeobirYob( value== 1 );
#line 6871 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 115 "../code/TablesEducation.mpp"
		Set_table_filter_TabEducFateByGroup( poParent->in_projected_time&& poDerivedStates->trigger_changes_calendar_year&& value== 0&& poParent->creation_type== CT_BIRTH );
#line 6874 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 154 "../code/TablesEducation.mpp"
		Set_table_filter_TabEducFateDistrYob( poDerivedStates->trigger_changes_calendar_year&& value== 0&& poParent->is_resident );
#line 6877 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 31 "../code/TablesVaccination.mpp"
		Set_table_filter_TabImmunizationChildren( poParent->in_projected_time&& value== 0&& poParent->creation_type== CT_BIRTH&& poParent->is_resident );
#line 6880 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 19 "../code/TablesVaccination.mpp"
		Set_table_filter_TabPrenatCare( poParent->in_projected_time&& value== 0&& poParent->creation_type== CT_BIRTH&& poParent->is_resident );
#line 6883 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 76 "../code/TablesEducation.mpp"
		Set_table_filter_tabPreSchool( WITHIN(SIM_YEAR_RANGE,poParent->year_of_birth)&& value== 10 );
#line 6886 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 211 "../code/TablesEducation.mpp"
		Set_table_filter_TabPrimSchoolOutOfSchool9to11( poParent->in_projected_time&& poParent->is_resident&& value>= 9&& value<= 11 );
#line 6889 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 231 "../code/TablesEducation.mpp"
		Set_table_filter_TabSchool2AttainmentsTab( poParent->in_projected_time&& poParent->is_resident&& value>= 20&& value< 25 );
#line 6892 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		if ( poParent->integer_age == 14 && value == 15 ) {
			poDerivedStates->Set_transitions_integer_age_14_15( poDerivedStates->transitions_integer_age_14_15 + 1 );
		}


		AGE_RANGE	OldValue;
		OldValue = poParent->integer_age;
		poParent->integer_age = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[14] >= 0 && ( poParent->tableinfo.tab14dim_filter || poParent->table_filter_tabHavingSpouse ) ) {
				if ( poParent->tableinfo.tab14events == -1 || poParent->tableinfo.tab14events < poParent->__events ) {
					if ( poParent->tableinfo.tab14events != -1 ) {
						poParent->tableinfo.tab14open0 = (real) poParent->CloseUpdate( 14, 0, poParent->tableinfo.tab14open0, poParent->tableinfo.tab14close0, poParent->tableinfo.tab14dim_filter, poParent->tableinfo.tab14cell, poParent->tableinfo.tab14total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab14open1 = (real) poParent->CloseUpdate( 14, 1, poParent->tableinfo.tab14open1, poParent->tableinfo.tab14close1, poParent->tableinfo.tab14dim_filter, poParent->tableinfo.tab14cell, poParent->tableinfo.tab14total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab14open2 = (real) poParent->CloseUpdate( 14, 2, poParent->tableinfo.tab14open2, poParent->tableinfo.tab14close2, poParent->tableinfo.tab14dim_filter, poParent->tableinfo.tab14cell, poParent->tableinfo.tab14total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab14events = poParent->__events;
					poParent->tableinfo.tab14close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab14close1 = (real) poDerivedStates->Get_duration_ever_union_TRUE();
					poParent->tableinfo.tab14close2 = (real) poDerivedStates->Get_duration_has_spouse_TRUE();
					poParent->tableinfo.tab14cell = (index) poParent->TableCellIndex( 14, 3,
						poParent->sex,
						OldValue - r_min_AGE_RANGE,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab14total_cell = (index) poParent->tableinfo.tab14cell;
					}
					else {
						poParent->tableinfo.tab14total_cell = (index) poParent->TableTotalCellIndex( 14, 3,
						poParent->sex,
						OldValue - r_min_AGE_RANGE,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab14dim_filter = poParent->table_filter_tabHavingSpouse;
				}
			}
		}
#line 96 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_age_last_birth( ( poParent->lives_with_dependent_child )? poParent->integer_age- poParent->age_youngest_child: TIME_INFINITE );
#line 6933 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 67 "../code/MortalityDetailed.mpp"
		Set_child_mortality( ( poParent->integer_age== 0&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_0): ( poParent->integer_age== 1&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_1): ( poParent->integer_age== 2&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_2): ( poParent->integer_age== 3&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_3): ( poParent->integer_age== 4&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_4): ( poParent->integer_age== 0&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_0): ( poParent->integer_age== 1&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_1): ( poParent->integer_age== 2&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_2): ( poParent->integer_age== 3&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_3): ( poParent->integer_age== 4&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_4): 0.0 );
#line 6936 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 81 "../code/FamilyPartnerMatching.mpp"
		Set_female_spouse_age( ( poParent->sex== FEMALE )? COERCE(FEMALE_SPOUSE_AGE,poParent->integer_age): MAX(FEMALE_SPOUSE_AGE) );
#line 6939 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 49 "../code/FertilityGeneral.mpp"
		Set_fertile_age( COERCE(FERTILE_AGE_RANGE,poParent->integer_age) );
#line 6942 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 103 "../code/HumanCapitalIndex.mpp"
		poParent->is_adult = ( poParent->integer_age>= 15 )? TRUE: FALSE;
#line 6945 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 20 "../code/TablesEducation.mpp"
		Set_is_child( ( poParent->integer_age< 15 ) );
#line 6948 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 47 "../code/FertilityGeneral.mpp"
		Set_is_fertile( ( poParent->sex== FEMALE&& WITHIN(FERTILE_AGE_RANGE,poParent->integer_age)&& ( poParent->creation_type== CT_BIRTH|| ( poParent->creation_type== CT_START&& poParent->calendar_year>= MIN(SIM_YEAR_RANGE) )|| ( poParent->creation_type== CT_POOL&& poParent->ever_resident )|| ( poParent->creation_type== CT_SCRATCH&& poParent->ever_resident ) ) )? TRUE: FALSE );
#line 6951 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 21 "../code/TablesEducation.mpp"
		Set_is_old( ( poParent->integer_age>= 60 ) );
#line 6954 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 78 "../code/FamilyPartnerMatching.mpp"
		Set_male_age( ( poParent->sex== MALE )? COERCE(MALE_SPOUSE_AGE,poParent->integer_age): MAX(MALE_SPOUSE_AGE) );
#line 6957 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		{
			int	nNewGrp;

			nNewGrp = Split( (double) poParent->integer_age, poDerivedStates->split_integer_age_TAB_AGEC5, 20, p_points_TAB_AGEC5 );
			if ( nNewGrp != poDerivedStates->split_integer_age_TAB_AGEC5 ) {
				poDerivedStates->Set_split_integer_age_TAB_AGEC5( nNewGrp );
			}
		}
#line 97 "../code/HumanCapitalIndex.mpp"
		Set_survived_early_years( ( poParent->integer_age>= 5 )? TRUE: FALSE );
#line 6968 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		if ( poParent->lSpouse != NULL ) {
			if ( poParent->lSpouse->bUpdates ) {
				poParent->lSpouse->UpdateTime();
				poParent->lSpouse->poStateFunctions->Set_lSpouse_integer_age( poParent->integer_age );
			}
		}
		if ( poParent->lHHMother != NULL ) {
			poParent->lHHMother->UpdateTime();
			poParent->lHHMother->poStateFunctions->Set_min_over_mlHHMotherChildren_integer_age( poParent->lHHMother->mlHHMotherChildren->MinOver_integer_age( poParent, poParent->integer_age ) );
		}
		// actorset dimension or order state
		if (poParent->actorset_filter_asPotentialImmigrantMothers)
		{
			asPotentialImmigrantMothers[RANGE_POS(AGE_RANGE,OldValue)][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)][poParent->geo]->Remove(poParent->poasPotentialImmigrantMothers);
			poParent->poasPotentialImmigrantMothers = asPotentialImmigrantMothers[RANGE_POS(AGE_RANGE, poParent->integer_age)][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)][poParent->geo]->NewNode(poParent);
		}
		if ( poParent->poBirthdayEvent->cRecalc == 0 ) {
			poParent->poBirthdayEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poBirthdayEvent );
		};
		if ( poParent->poFertilityBaseEvent->cRecalc == 0 ) {
			poParent->poFertilityBaseEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFertilityBaseEvent );
		};
		if ( poParent->poFertilityDetailedEvent->cRecalc == 0 ) {
			poParent->poFertilityDetailedEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFertilityDetailedEvent );
		};
		if ( poParent->poFirstUnionFormationEvent->cRecalc == 0 ) {
			poParent->poFirstUnionFormationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstUnionFormationEvent );
		};
		if ( poParent->poMortalityBaseEvent->cRecalc == 0 ) {
			poParent->poMortalityBaseEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMortalityBaseEvent );
		};
		if ( poParent->poMortalityDetailedEvent->cRecalc == 0 ) {
			poParent->poMortalityDetailedEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMortalityDetailedEvent );
		};
	}
	return poParent->integer_age = value;
}

logical PersonStateFunctions::Set_is_alive( logical value )
{

	value = value != 0;
	if ( poParent->is_alive != value )
	{
#line 40 "../code/ActorPerson.mpp"
		Set_actorset_filter_asAllPerson( value );
#line 7021 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 14 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMale( value&& poParent->sex== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& poParent->is_mortal );
#line 7024 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 18 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleByType( value&& poParent->sex== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& poParent->is_mortal );
#line 7027 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 32 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleByTypeImmiScratch( value&& poParent->creation_type== CT_SCRATCH&& poParent->sex== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->is_mortal ) );
#line 7030 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleImmiScratch( value&& poParent->creation_type== CT_SCRATCH&& poParent->sex== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->is_mortal ) );
#line 7033 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 22 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAge( value&& poParent->sex== FEMALE&& poParent->has_spouse&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& poParent->is_mortal );
#line 7036 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 36 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch( value&& poParent->creation_type== CT_SCRATCH&& poParent->sex== FEMALE&& poParent->has_spouse&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->is_mortal ) );
#line 7039 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 10 "../code/ImmigrationFromScratch.mpp"
		Set_actorset_filter_asPotentialImmigrantMothers( value&& poParent->sex== FEMALE&& poParent->creation_type== CT_SCRATCH&& !( poParent->ever_resident )&& poParent->family_role!= FR_CHILD&& poParent->ever_union );
#line 7042 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 25 "../code/EducationPrimaryTracking.mpp"
		Set_actorset_filter_asPotentialSchoolOneStudent( value&& poParent->educ_one_status!= EOS_OUT&& poParent->educ_one_fate!= EOL_LOW );
#line 7045 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/EducationPrimaryTracking.mpp"
		Set_actorset_filter_asPotentialSchoolOneStudentToProcess( value&& poParent->educ_one_to_process );
#line 7048 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 29 "../code/EducationSecondaryBase.mpp"
		Set_actorset_filter_asPotentialSchoolTwoStudentToProcess( value&& poParent->educ_two_to_process );
#line 7051 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 18 "../code/ChildVaccination.mpp"
		Set_actorset_filter_asResidentBabies( value&& poParent->integer_age== 0&& poParent->is_resident&& poParent->creation_type== CT_BIRTH );
#line 7054 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 30 "../code/EducationPrimaryBase.mpp"
		Set_actorset_filter_asResidentsAge0Sex( value&& poParent->integer_age== 0&& poParent->is_resident );
#line 7057 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 27 "../code/EducationPrimaryBase.mpp"
		Set_actorset_filter_asResidentsAge0SexGeo( value&& poParent->is_resident&& poParent->integer_age== 0 );
#line 7060 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 35 "../code/EducationPrimaryTransmission.mpp"
		Set_actorset_filter_asSimBornAge0( value&& poParent->creation_type== CT_BIRTH&& poParent->integer_age== 0 );
#line 7063 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 11 "../code/MigrationGeneral.mpp"
		Set_actorset_filter_asWantToMove( value );
#line 7066 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_actorset_filter_asWomenNoChildren( value&& poParent->sex== FEMALE&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->lives_with_dependent_child )&& poParent->is_mortal&& !( poParent->is_blocked_from_marriage ) );
#line 7069 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 23 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_actorset_filter_asWomenWithChildren( value&& poParent->sex== FEMALE&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& poParent->lives_with_dependent_child&& poParent->is_mortal&& !( poParent->is_blocked_from_marriage ) );
#line 7072 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		if ( poParent->is_alive == FALSE && value == TRUE ) {
			poDerivedStates->Set_transitions_is_alive_FALSE_TRUE( poDerivedStates->transitions_is_alive_FALSE_TRUE + 1 );
		}
		if ( poParent->is_alive == TRUE && value == FALSE ) {
			poDerivedStates->Set_transitions_is_alive_TRUE_FALSE( poDerivedStates->transitions_is_alive_TRUE_FALSE + 1 );
		}
		if ( poParent->is_alive == TRUE && value == FALSE ) {
			poDerivedStates->Set_value_at_transitions_is_alive_TRUE_FALSE_adult_survival( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_adult_survival + poParent->adult_survival );
		}
		if ( poParent->is_alive == TRUE && value == FALSE ) {
			poDerivedStates->Set_value_at_transitions_is_alive_TRUE_FALSE_ind_hci( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_ind_hci + poParent->ind_hci );
		}
		if ( poParent->is_alive == TRUE && value == FALSE ) {
			poDerivedStates->Set_value_at_transitions_is_alive_TRUE_FALSE_is_stunted( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_is_stunted + poParent->is_stunted );
		}
		if ( poParent->is_alive == TRUE && value == FALSE ) {
			poDerivedStates->Set_value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling + poParent->quality_of_schooling );
		}
		if ( poParent->is_alive == TRUE && value == FALSE ) {
			poDerivedStates->Set_value_at_transitions_is_alive_TRUE_FALSE_survived_early_years( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_survived_early_years + poParent->survived_early_years );
		}
		if ( poParent->is_alive == TRUE && value == FALSE ) {
			poDerivedStates->Set_value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling( poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling + poParent->years_of_schooling );
		}


		logical	OldValue;
		OldValue = poParent->is_alive;
		poParent->is_alive = value;

#line 73 "../code/TablesPopulation.mpp"
		Set_is_alive_value_in( poParent->is_alive );
#line 7105 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		if ( poParent->lHHFather != NULL ) {
			poParent->lHHFather->UpdateTime();
			poParent->lHHFather->poStateFunctions->Set_sum_over_mlHHFatherChildren_is_alive( poParent->lHHFather->sum_over_mlHHFatherChildren_is_alive + poParent->is_alive - OldValue);
		}
		if ( poParent->lHHMother != NULL ) {
			poParent->lHHMother->UpdateTime();
			poParent->lHHMother->poStateFunctions->Set_sum_over_mlHHMotherChildren_is_alive( poParent->lHHMother->sum_over_mlHHMotherChildren_is_alive + poParent->is_alive - OldValue);
		}
	}
	return poParent->is_alive = value;
}

logical PersonStateFunctions::Set_is_blocked_from_marriage( logical value )
{

	value = value != 0;
	if ( poParent->is_blocked_from_marriage != value )
	{
#line 28 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_actorset_filter_asWomenNoChildren( poParent->is_alive&& poParent->sex== FEMALE&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->lives_with_dependent_child )&& poParent->is_mortal&& !( value ) );
#line 7126 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 23 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_actorset_filter_asWomenWithChildren( poParent->is_alive&& poParent->sex== FEMALE&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& poParent->lives_with_dependent_child&& poParent->is_mortal&& !( value ) );
#line 7129 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		logical	OldValue;
		OldValue = poParent->is_blocked_from_marriage;
		poParent->is_blocked_from_marriage = value;

	}
	return poParent->is_blocked_from_marriage = value;
}

logical PersonStateFunctions::Set_is_immunized( logical value )
{

	value = value != 0;
	if ( poParent->is_immunized != value )
	{
		if ( poParent->is_immunized == TRUE ) {
			poDerivedStates->Set_duration_is_immunized_TRUE( poDerivedStates->Get_duration_is_immunized_TRUE() );
		}
		else if ( value == TRUE ) {
			poDerivedStates->duration_is_immunized_TRUE_time = poParent->__time;
		}


		logical	OldValue;
		OldValue = poParent->is_immunized;
		poParent->is_immunized = value;

	}
	return poParent->is_immunized = value;
}

logical PersonStateFunctions::Set_is_ready_for_birthtables( logical value )
{

	value = value != 0;
	if ( poParent->is_ready_for_birthtables != value )
	{
		bool bScheduleTimeChange = false;	//indicates that schedule event should be recalculated
		if (	(poParent->schedule_event_time == poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE_time))
		{
			//force timeScheduleEvent to calculate min time
			poParent->schedule_event_time = TIME_INFINITE;
		//indicate that schedule event needs to be recalculated
			bScheduleTimeChange = true;
		}


		logical	OldValue;
		OldValue = poParent->is_ready_for_birthtables;
		poParent->is_ready_for_birthtables = value;

		if ( value == TRUE ) 
		{
			poDerivedStates->Set_trigger_entrances_is_ready_for_birthtables_TRUE( TRUE );
			poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE_time = poParent->__time;
			if ( poParent->schedule_event_time > poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE_time)
			{
				poParent->schedule_event_time = poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE_time;
				//indicate that schedule event should be recalculated
				bScheduleTimeChange = true;
			}
		}

		if (bScheduleTimeChange && poParent->poScheduleEvent->cRecalc == 0)
		{
			poParent->poScheduleEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent(poParent->poScheduleEvent);
		}
	}
	return poParent->is_ready_for_birthtables = value;
}

logical PersonStateFunctions::Set_is_stunted( logical value )
{

	value = value != 0;
	if ( poParent->is_stunted != value )
	{


		logical	OldValue;
		OldValue = poParent->is_stunted;
		poParent->is_stunted = value;

#line 23 "../code/TablesStunting.mpp"
		Set_is_stunted_value_in( poParent->is_stunted );
#line 7217 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->is_stunted = value;
}

double PersonStateFunctions::Set_mother_age_at_birth( double value )
{

	if ( poParent->mother_age_at_birth != value )
	{


		double	OldValue;
		OldValue = poParent->mother_age_at_birth;
		poParent->mother_age_at_birth = value;

#line 79 "../code/MortalityDetailed.mpp"
		Set_child_mortality_group( ( poParent->mother_age_at_birth>= 17&& poParent->educ_mother== EOL_HIGH )? CMG_00: ( poParent->mother_age_at_birth>= 17&& poParent->educ_mother== EOL_MEDIUM )? CMG_01: ( poParent->mother_age_at_birth>= 17&& poParent->educ_mother== EOL_LOW )? CMG_02: ( poParent->mother_age_at_birth>= 15&& poParent->educ_mother== EOL_HIGH )? CMG_03: ( poParent->mother_age_at_birth>= 15&& poParent->educ_mother== EOL_MEDIUM )? CMG_04: ( poParent->mother_age_at_birth>= 15&& poParent->educ_mother== EOL_LOW )? CMG_05: ( poParent->educ_mother== EOL_HIGH )? CMG_06: ( poParent->educ_mother== EOL_MEDIUM )? CMG_07: CMG_08 );
#line 7235 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->mother_age_at_birth = value;
}

PARITY_RANGE PersonStateFunctions::Set_parity( PARITY_RANGE value )
{

	CHKLMT(value, parity, r_min_PARITY_RANGE, r_max_PARITY_RANGE);
	if ( value < r_min_PARITY_RANGE ) {
		value = r_min_PARITY_RANGE;
	}
	else if ( value > r_max_PARITY_RANGE ) {
		value = r_max_PARITY_RANGE;
	}
	if ( poParent->parity != value )
	{


		PARITY_RANGE	OldValue;
		OldValue = poParent->parity;
		poParent->parity = value;

		if ( poParent->poFertilityDetailedEvent->cRecalc == 0 ) {
			poParent->poFertilityDetailedEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFertilityDetailedEvent );
		};
	}
	return poParent->parity = value;
}

logical PersonStateFunctions::Set_preschool_is_decided( logical value )
{

	value = value != 0;
	if ( poParent->preschool_is_decided != value )
	{


		logical	OldValue;
		OldValue = poParent->preschool_is_decided;
		poParent->preschool_is_decided = value;

		if ( poParent->poSetPreschoolYearsEvent->cRecalc == 0 ) {
			poParent->poSetPreschoolYearsEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poSetPreschoolYearsEvent );
		};
	}
	return poParent->preschool_is_decided = value;
}

logical PersonStateFunctions::Set_ready_to_set_alive( logical value )
{

	value = value != 0;
	if ( poParent->ready_to_set_alive != value )
	{


		logical	OldValue;
		OldValue = poParent->ready_to_set_alive;
		poParent->ready_to_set_alive = value;

		if ( poParent->poSetAliveEvent->cRecalc == 0 ) {
			poParent->poSetAliveEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poSetAliveEvent );
		};
	}
	return poParent->ready_to_set_alive = value;
}

SEX PersonStateFunctions::Set_sex( SEX value )
{

	CHKLMT(value, sex, 0, rc_size_SEX - 1);
	if ( value < 0 ) {
		value = (SEX) 0;
	}
	else if ( value >= rc_size_SEX ) {
		value = (SEX) ( rc_size_SEX - 1 );
	}
	if ( poParent->sex != value )
	{
#line 14 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMale( poParent->is_alive&& value== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& poParent->is_mortal );
#line 7320 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 18 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleByType( poParent->is_alive&& value== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& poParent->is_mortal );
#line 7323 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 32 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleByTypeImmiScratch( poParent->is_alive&& poParent->creation_type== CT_SCRATCH&& value== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->is_mortal ) );
#line 7326 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleImmiScratch( poParent->is_alive&& poParent->creation_type== CT_SCRATCH&& value== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->is_mortal ) );
#line 7329 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 22 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAge( poParent->is_alive&& value== FEMALE&& poParent->has_spouse&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& poParent->is_mortal );
#line 7332 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 36 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch( poParent->is_alive&& poParent->creation_type== CT_SCRATCH&& value== FEMALE&& poParent->has_spouse&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->is_mortal ) );
#line 7335 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 10 "../code/ImmigrationFromScratch.mpp"
		Set_actorset_filter_asPotentialImmigrantMothers( poParent->is_alive&& value== FEMALE&& poParent->creation_type== CT_SCRATCH&& !( poParent->ever_resident )&& poParent->family_role!= FR_CHILD&& poParent->ever_union );
#line 7338 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_actorset_filter_asWomenNoChildren( poParent->is_alive&& value== FEMALE&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->lives_with_dependent_child )&& poParent->is_mortal&& !( poParent->is_blocked_from_marriage ) );
#line 7341 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 23 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_actorset_filter_asWomenWithChildren( poParent->is_alive&& value== FEMALE&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& poParent->lives_with_dependent_child&& poParent->is_mortal&& !( poParent->is_blocked_from_marriage ) );
#line 7344 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		SEX	OldValue;
		OldValue = poParent->sex;
		poParent->sex = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[1] >= 0 && ( poParent->tableinfo.tab1dim_filter || poParent->table_filter_tabEducationFateGeobirYob ) ) {
				if ( poParent->tableinfo.tab1events == -1 || poParent->tableinfo.tab1events < poParent->__events ) {
					if ( poParent->tableinfo.tab1events != -1 ) {
						poParent->tableinfo.tab1open0 = (real) poParent->CloseUpdate( 1, 0, poParent->tableinfo.tab1open0, poParent->tableinfo.tab1close0, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab1open1 = (real) poParent->CloseUpdate( 1, 1, poParent->tableinfo.tab1open1, poParent->tableinfo.tab1close1, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab1open2 = (real) poParent->CloseUpdate( 1, 2, poParent->tableinfo.tab1open2, poParent->tableinfo.tab1close2, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab1open3 = (real) poParent->CloseUpdate( 1, 3, poParent->tableinfo.tab1open3, poParent->tableinfo.tab1close3, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab1events = poParent->__events;
					poParent->tableinfo.tab1close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab1close1 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH();
					poParent->tableinfo.tab1close2 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_LOW();
					poParent->tableinfo.tab1close3 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM();
					poParent->tableinfo.tab1cell = (index) poParent->TableCellIndex( 1, 3,
						OldValue,
						poParent->geo_birth,
						poParent->year_of_birth - r_min_ALL_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab1total_cell = (index) poParent->tableinfo.tab1cell;
					}
					else {
						poParent->tableinfo.tab1total_cell = (index) poParent->TableTotalCellIndex( 1, 3,
						OldValue,
						poParent->geo_birth,
						poParent->year_of_birth - r_min_ALL_YEAR_RANGE );
					}
					poParent->tableinfo.tab1dim_filter = poParent->table_filter_tabEducationFateGeobirYob;
				}
			}
			if ( gplTableIndexes[2] >= 0 && ( poParent->tableinfo.tab2dim_filter || poParent->table_filter_PopPyramidByEduc ) ) {
				if ( poParent->tableinfo.tab2events == -1 || poParent->tableinfo.tab2events < poParent->__events ) {
					if ( poParent->tableinfo.tab2events != -1 ) {
						poParent->tableinfo.tab2open0 = (real) poParent->CloseUpdate( 2, 0, poParent->tableinfo.tab2open0, poParent->tableinfo.tab2close0, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab2open1 = (real) poParent->CloseUpdate( 2, 1, poParent->tableinfo.tab2open1, poParent->tableinfo.tab2close1, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab2open2 = (real) poParent->CloseUpdate( 2, 2, poParent->tableinfo.tab2open2, poParent->tableinfo.tab2close2, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab2events = poParent->__events;
					poParent->tableinfo.tab2close0 = (real) poDerivedStates->Get_duration_educ_one_level_EOL_HIGH();
					poParent->tableinfo.tab2close1 = (real) poDerivedStates->Get_duration_educ_one_level_EOL_LOW();
					poParent->tableinfo.tab2close2 = (real) poDerivedStates->Get_duration_educ_one_level_EOL_MEDIUM();
					poParent->tableinfo.tab2cell = (index) poParent->TableCellIndex( 2, 4,
						poParent->sel_years,
						OldValue,
						poParent->geo_nat,
						poDerivedStates->split_integer_age_TAB_AGEC5 );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab2total_cell = (index) poParent->tableinfo.tab2cell;
					}
					else {
						poParent->tableinfo.tab2total_cell = (index) poParent->TableTotalCellIndex( 2, 4,
						poParent->sel_years,
						OldValue,
						poParent->geo_nat,
						poDerivedStates->split_integer_age_TAB_AGEC5 );
					}
					poParent->tableinfo.tab2dim_filter = poParent->table_filter_PopPyramidByEduc;
				}
			}
			if ( gplTableIndexes[5] >= 0 && ( poParent->tableinfo.tab5dim_filter || poParent->table_filter_TabPopProvAgeEducSex ) ) {
				if ( poParent->tableinfo.tab5events == -1 || poParent->tableinfo.tab5events < poParent->__events ) {
					if ( poParent->tableinfo.tab5events != -1 ) {
						poParent->tableinfo.tab5open0 = (real) poParent->CloseUpdate( 5, 0, poParent->tableinfo.tab5open0, poParent->tableinfo.tab5close0, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open1 = (real) poParent->CloseUpdate( 5, 1, poParent->tableinfo.tab5open1, poParent->tableinfo.tab5close1, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open2 = (real) poParent->CloseUpdate( 5, 2, poParent->tableinfo.tab5open2, poParent->tableinfo.tab5close2, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open3 = (real) poParent->CloseUpdate( 5, 3, poParent->tableinfo.tab5open3, poParent->tableinfo.tab5close3, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open4 = (real) poParent->CloseUpdate( 5, 4, poParent->tableinfo.tab5open4, poParent->tableinfo.tab5close4, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab5events = poParent->__events;
					poParent->tableinfo.tab5close0 = (real) poDerivedStates->Get_duration_is_aded0_TRUE();
					poParent->tableinfo.tab5close1 = (real) poDerivedStates->Get_duration_is_aded1_TRUE();
					poParent->tableinfo.tab5close2 = (real) poDerivedStates->Get_duration_is_aded2_TRUE();
					poParent->tableinfo.tab5close3 = (real) poDerivedStates->Get_duration_is_child_TRUE();
					poParent->tableinfo.tab5close4 = (real) poDerivedStates->Get_duration_is_old_TRUE();
					poParent->tableinfo.tab5cell = (index) poParent->TableCellIndex( 5, 3,
						poParent->geo_nat,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab5total_cell = (index) poParent->tableinfo.tab5cell;
					}
					else {
						poParent->tableinfo.tab5total_cell = (index) poParent->TableTotalCellIndex( 5, 3,
						poParent->geo_nat,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab5dim_filter = poParent->table_filter_TabPopProvAgeEducSex;
				}
			}
			if ( gplTableIndexes[7] >= 0 && ( poParent->tableinfo.tab7dim_filter || poParent->table_filter_TabEduc15ByDistrict ) ) {
				if ( poParent->tableinfo.tab7events == -1 || poParent->tableinfo.tab7events < poParent->__events ) {
					if ( poParent->tableinfo.tab7events != -1 ) {
						poParent->tableinfo.tab7open0 = (real) poParent->CloseUpdate( 7, 0, poParent->tableinfo.tab7open0, poParent->tableinfo.tab7close0, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab7open1 = (real) poParent->CloseUpdate( 7, 1, poParent->tableinfo.tab7open1, poParent->tableinfo.tab7close1, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab7open2 = (real) poParent->CloseUpdate( 7, 2, poParent->tableinfo.tab7open2, poParent->tableinfo.tab7close2, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab7open3 = (real) poParent->CloseUpdate( 7, 3, poParent->tableinfo.tab7open3, poParent->tableinfo.tab7close3, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab7events = poParent->__events;
					poParent->tableinfo.tab7close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab7close1 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH();
					poParent->tableinfo.tab7close2 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_LOW();
					poParent->tableinfo.tab7close3 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM();
					poParent->tableinfo.tab7cell = (index) poParent->TableCellIndex( 7, 3,
						poParent->geo_nat,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab7total_cell = (index) poParent->tableinfo.tab7cell;
					}
					else {
						poParent->tableinfo.tab7total_cell = (index) poParent->TableTotalCellIndex( 7, 3,
						poParent->geo_nat,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab7dim_filter = poParent->table_filter_TabEduc15ByDistrict;
				}
			}
			if ( gplTableIndexes[8] >= 0 && ( poParent->tableinfo.tab8dim_filter || poParent->table_filter_TabEduc15ByDistrictBirth ) ) {
				if ( poParent->tableinfo.tab8events == -1 || poParent->tableinfo.tab8events < poParent->__events ) {
					if ( poParent->tableinfo.tab8events != -1 ) {
						poParent->tableinfo.tab8open0 = (real) poParent->CloseUpdate( 8, 0, poParent->tableinfo.tab8open0, poParent->tableinfo.tab8close0, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab8open1 = (real) poParent->CloseUpdate( 8, 1, poParent->tableinfo.tab8open1, poParent->tableinfo.tab8close1, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab8open2 = (real) poParent->CloseUpdate( 8, 2, poParent->tableinfo.tab8open2, poParent->tableinfo.tab8close2, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab8open3 = (real) poParent->CloseUpdate( 8, 3, poParent->tableinfo.tab8open3, poParent->tableinfo.tab8close3, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab8events = poParent->__events;
					poParent->tableinfo.tab8close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab8close1 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH();
					poParent->tableinfo.tab8close2 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_LOW();
					poParent->tableinfo.tab8close3 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM();
					poParent->tableinfo.tab8cell = (index) poParent->TableCellIndex( 8, 3,
						poParent->geo_birth,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab8total_cell = (index) poParent->tableinfo.tab8cell;
					}
					else {
						poParent->tableinfo.tab8total_cell = (index) poParent->TableTotalCellIndex( 8, 3,
						poParent->geo_birth,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab8dim_filter = poParent->table_filter_TabEduc15ByDistrictBirth;
				}
			}
			if ( gplTableIndexes[13] >= 0 && ( poParent->tableinfo.tab13dim_filter || poParent->table_filter_TabSchool2AttainmentsTab ) ) {
				if ( poParent->tableinfo.tab13events == -1 || poParent->tableinfo.tab13events < poParent->__events ) {
					if ( poParent->tableinfo.tab13events != -1 ) {
						poParent->tableinfo.tab13open0 = (real) poParent->CloseUpdate( 13, 0, poParent->tableinfo.tab13open0, poParent->tableinfo.tab13close0, poParent->tableinfo.tab13dim_filter, poParent->tableinfo.tab13cell, poParent->tableinfo.tab13total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab13open1 = (real) poParent->CloseUpdate( 13, 1, poParent->tableinfo.tab13open1, poParent->tableinfo.tab13close1, poParent->tableinfo.tab13dim_filter, poParent->tableinfo.tab13cell, poParent->tableinfo.tab13total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab13open2 = (real) poParent->CloseUpdate( 13, 2, poParent->tableinfo.tab13open2, poParent->tableinfo.tab13close2, poParent->tableinfo.tab13dim_filter, poParent->tableinfo.tab13cell, poParent->tableinfo.tab13total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab13open3 = (real) poParent->CloseUpdate( 13, 3, poParent->tableinfo.tab13open3, poParent->tableinfo.tab13close3, poParent->tableinfo.tab13dim_filter, poParent->tableinfo.tab13cell, poParent->tableinfo.tab13total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab13events = poParent->__events;
					poParent->tableinfo.tab13close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab13close1 = (real) poDerivedStates->Get_duration_educ_two_level_ETL_GRADUATE();
					poParent->tableinfo.tab13close2 = (real) poDerivedStates->Get_duration_educ_two_level_ETL_INCOMPLETE();
					poParent->tableinfo.tab13close3 = (real) poDerivedStates->Get_duration_educ_two_level_ETL_NEVER();
					poParent->tableinfo.tab13cell = (index) poParent->TableCellIndex( 13, 2,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab13total_cell = (index) poParent->tableinfo.tab13cell;
					}
					else {
						poParent->tableinfo.tab13total_cell = (index) poParent->TableTotalCellIndex( 13, 2,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab13dim_filter = poParent->table_filter_TabSchool2AttainmentsTab;
				}
			}
			if ( gplTableIndexes[14] >= 0 && ( poParent->tableinfo.tab14dim_filter || poParent->table_filter_tabHavingSpouse ) ) {
				if ( poParent->tableinfo.tab14events == -1 || poParent->tableinfo.tab14events < poParent->__events ) {
					if ( poParent->tableinfo.tab14events != -1 ) {
						poParent->tableinfo.tab14open0 = (real) poParent->CloseUpdate( 14, 0, poParent->tableinfo.tab14open0, poParent->tableinfo.tab14close0, poParent->tableinfo.tab14dim_filter, poParent->tableinfo.tab14cell, poParent->tableinfo.tab14total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab14open1 = (real) poParent->CloseUpdate( 14, 1, poParent->tableinfo.tab14open1, poParent->tableinfo.tab14close1, poParent->tableinfo.tab14dim_filter, poParent->tableinfo.tab14cell, poParent->tableinfo.tab14total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab14open2 = (real) poParent->CloseUpdate( 14, 2, poParent->tableinfo.tab14open2, poParent->tableinfo.tab14close2, poParent->tableinfo.tab14dim_filter, poParent->tableinfo.tab14cell, poParent->tableinfo.tab14total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab14events = poParent->__events;
					poParent->tableinfo.tab14close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab14close1 = (real) poDerivedStates->Get_duration_ever_union_TRUE();
					poParent->tableinfo.tab14close2 = (real) poDerivedStates->Get_duration_has_spouse_TRUE();
					poParent->tableinfo.tab14cell = (index) poParent->TableCellIndex( 14, 3,
						OldValue,
						poParent->integer_age - r_min_AGE_RANGE,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab14total_cell = (index) poParent->tableinfo.tab14cell;
					}
					else {
						poParent->tableinfo.tab14total_cell = (index) poParent->TableTotalCellIndex( 14, 3,
						OldValue,
						poParent->integer_age - r_min_AGE_RANGE,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab14dim_filter = poParent->table_filter_tabHavingSpouse;
				}
			}
			if ( gplTableIndexes[15] >= 0 && ( poParent->tableinfo.tab15dim_filter || poParent->table_filter_tabHCI ) ) {
				if ( poParent->tableinfo.tab15events == -1 || poParent->tableinfo.tab15events < poParent->__events ) {
					if ( poParent->tableinfo.tab15events != -1 ) {
						poParent->tableinfo.tab15open0 = (real) poParent->CloseUpdate( 15, 0, poParent->tableinfo.tab15open0, poParent->tableinfo.tab15close0, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open1 = (real) poParent->CloseUpdate( 15, 1, poParent->tableinfo.tab15open1, poParent->tableinfo.tab15close1, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open2 = (real) poParent->CloseUpdate( 15, 2, poParent->tableinfo.tab15open2, poParent->tableinfo.tab15close2, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open3 = (real) poParent->CloseUpdate( 15, 3, poParent->tableinfo.tab15open3, poParent->tableinfo.tab15close3, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open4 = (real) poParent->CloseUpdate( 15, 4, poParent->tableinfo.tab15open4, poParent->tableinfo.tab15close4, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open5 = (real) poParent->CloseUpdate( 15, 5, poParent->tableinfo.tab15open5, poParent->tableinfo.tab15close5, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open6 = (real) poParent->CloseUpdate( 15, 6, poParent->tableinfo.tab15open6, poParent->tableinfo.tab15close6, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open7 = (real) poParent->CloseUpdate( 15, 7, poParent->tableinfo.tab15open7, poParent->tableinfo.tab15close7, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open8 = (real) poParent->CloseUpdate( 15, 8, poParent->tableinfo.tab15open8, poParent->tableinfo.tab15close8, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab15events = poParent->__events;
					poParent->tableinfo.tab15close0 = (real) poDerivedStates->transitions_integer_age_14_15;
					poParent->tableinfo.tab15close1 = (real) poDerivedStates->transitions_is_alive_TRUE_FALSE;
					poParent->tableinfo.tab15close2 = (real) poDerivedStates->transitions_survived_early_years_FALSE_TRUE;
					poParent->tableinfo.tab15close3 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_adult_survival;
					poParent->tableinfo.tab15close4 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_ind_hci;
					poParent->tableinfo.tab15close5 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_is_stunted;
					poParent->tableinfo.tab15close6 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling;
					poParent->tableinfo.tab15close7 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_survived_early_years;
					poParent->tableinfo.tab15close8 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling;
					poParent->tableinfo.tab15cell = (index) poParent->TableCellIndex( 15, 3,
						OldValue,
						poDerivedStates->aggregate_region_birth_REGION_NAT,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab15total_cell = (index) poParent->tableinfo.tab15cell;
					}
					else {
						poParent->tableinfo.tab15total_cell = (index) poParent->TableTotalCellIndex( 15, 3,
						OldValue,
						poDerivedStates->aggregate_region_birth_REGION_NAT,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab15dim_filter = poParent->table_filter_tabHCI;
				}
			}
			if ( gplTableIndexes[16] >= 0 && ( poParent->tableinfo.tab16dim_filter || poParent->table_filter_TabHCIDistrict ) ) {
				if ( poParent->tableinfo.tab16events == -1 || poParent->tableinfo.tab16events < poParent->__events ) {
					if ( poParent->tableinfo.tab16events != -1 ) {
						poParent->tableinfo.tab16open0 = (real) poParent->CloseUpdate( 16, 0, poParent->tableinfo.tab16open0, poParent->tableinfo.tab16close0, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open1 = (real) poParent->CloseUpdate( 16, 1, poParent->tableinfo.tab16open1, poParent->tableinfo.tab16close1, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open2 = (real) poParent->CloseUpdate( 16, 2, poParent->tableinfo.tab16open2, poParent->tableinfo.tab16close2, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open3 = (real) poParent->CloseUpdate( 16, 3, poParent->tableinfo.tab16open3, poParent->tableinfo.tab16close3, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open4 = (real) poParent->CloseUpdate( 16, 4, poParent->tableinfo.tab16open4, poParent->tableinfo.tab16close4, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open5 = (real) poParent->CloseUpdate( 16, 5, poParent->tableinfo.tab16open5, poParent->tableinfo.tab16close5, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open6 = (real) poParent->CloseUpdate( 16, 6, poParent->tableinfo.tab16open6, poParent->tableinfo.tab16close6, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open7 = (real) poParent->CloseUpdate( 16, 7, poParent->tableinfo.tab16open7, poParent->tableinfo.tab16close7, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open8 = (real) poParent->CloseUpdate( 16, 8, poParent->tableinfo.tab16open8, poParent->tableinfo.tab16close8, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab16events = poParent->__events;
					poParent->tableinfo.tab16close0 = (real) poDerivedStates->transitions_integer_age_14_15;
					poParent->tableinfo.tab16close1 = (real) poDerivedStates->transitions_is_alive_TRUE_FALSE;
					poParent->tableinfo.tab16close2 = (real) poDerivedStates->transitions_survived_early_years_FALSE_TRUE;
					poParent->tableinfo.tab16close3 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_adult_survival;
					poParent->tableinfo.tab16close4 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_ind_hci;
					poParent->tableinfo.tab16close5 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_is_stunted;
					poParent->tableinfo.tab16close6 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling;
					poParent->tableinfo.tab16close7 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_survived_early_years;
					poParent->tableinfo.tab16close8 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling;
					poParent->tableinfo.tab16cell = (index) poParent->TableCellIndex( 16, 3,
						OldValue,
						poDerivedStates->aggregate_geo_birth_GEO_NAT,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab16total_cell = (index) poParent->tableinfo.tab16cell;
					}
					else {
						poParent->tableinfo.tab16total_cell = (index) poParent->TableTotalCellIndex( 16, 3,
						OldValue,
						poDerivedStates->aggregate_geo_birth_GEO_NAT,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab16dim_filter = poParent->table_filter_TabHCIDistrict;
				}
			}
			if ( gplTableIndexes[21] >= 0 && ( poParent->tableinfo.tab21dim_filter || poParent->table_filter_tabStuntingSexRegMotherYob ) ) {
				if ( poParent->tableinfo.tab21events == -1 || poParent->tableinfo.tab21events < poParent->__events ) {
					if ( poParent->tableinfo.tab21events != -1 ) {
						poParent->tableinfo.tab21open0 = (real) poParent->CloseUpdate( 21, 0, poParent->tableinfo.tab21open0, poParent->tableinfo.tab21close0, poParent->tableinfo.tab21dim_filter, poParent->tableinfo.tab21cell, poParent->tableinfo.tab21total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						Set_table_unit_tabStuntingSexRegMotherYob(poParent->table_unit_tabStuntingSexRegMotherYob+1);
						poParent->tableinfo.tab21close1 = (real) poParent->table_unit_tabStuntingSexRegMotherYob;
						poParent->tableinfo.tab21open1 = (real) poParent->CloseUpdate( 21, 1, poParent->tableinfo.tab21open1, poParent->tableinfo.tab21close1, poParent->tableinfo.tab21dim_filter, poParent->tableinfo.tab21cell, poParent->tableinfo.tab21total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab21events = poParent->__events;
					poParent->tableinfo.tab21close0 = (real) poParent->is_stunted_value_in;
					poParent->tableinfo.tab21close1 = (real) poParent->table_unit_tabStuntingSexRegMotherYob;
					poParent->tableinfo.tab21cell = (index) poParent->TableCellIndex( 21, 4,
						OldValue,
						poParent->educ_mother,
						poParent->region_nat,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab21total_cell = (index) poParent->tableinfo.tab21cell;
					}
					else {
						poParent->tableinfo.tab21total_cell = (index) poParent->TableTotalCellIndex( 21, 4,
						OldValue,
						poParent->educ_mother,
						poParent->region_nat,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab21dim_filter = poParent->table_filter_tabStuntingSexRegMotherYob;
				}
			}
			if ( gplTableIndexes[23] >= 0 && ( poParent->tableinfo.tab23dim_filter || poParent->table_filter_TabImmunizationChildren ) ) {
				if ( poParent->tableinfo.tab23events == -1 || poParent->tableinfo.tab23events < poParent->__events ) {
					if ( poParent->tableinfo.tab23events != -1 ) {
						poParent->tableinfo.tab23open0 = (real) poParent->CloseUpdate( 23, 0, poParent->tableinfo.tab23open0, poParent->tableinfo.tab23close0, poParent->tableinfo.tab23dim_filter, poParent->tableinfo.tab23cell, poParent->tableinfo.tab23total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab23open1 = (real) poParent->CloseUpdate( 23, 1, poParent->tableinfo.tab23open1, poParent->tableinfo.tab23close1, poParent->tableinfo.tab23dim_filter, poParent->tableinfo.tab23cell, poParent->tableinfo.tab23total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab23events = poParent->__events;
					poParent->tableinfo.tab23close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab23close1 = (real) poDerivedStates->Get_duration_is_immunized_TRUE();
					poParent->tableinfo.tab23cell = (index) poParent->TableCellIndex( 23, 6,
						OldValue,
						poParent->ethnicity_short,
						poParent->got_prenat_care,
						poParent->region_nat,
						poParent->educ_mother,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab23total_cell = (index) poParent->tableinfo.tab23cell;
					}
					else {
						poParent->tableinfo.tab23total_cell = (index) poParent->TableTotalCellIndex( 23, 6,
						OldValue,
						poParent->ethnicity_short,
						poParent->got_prenat_care,
						poParent->region_nat,
						poParent->educ_mother,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab23dim_filter = poParent->table_filter_TabImmunizationChildren;
				}
			}
		}
#line 85 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_age_youngest_child( ( poParent->sex== FEMALE&& poParent->count_mlHHMotherChildren> 0 )? double(poParent->min_over_mlHHMotherChildren_integer_age): TIME_INFINITE );
#line 7694 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 67 "../code/MortalityDetailed.mpp"
		Set_child_mortality( ( poParent->integer_age== 0&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_0): ( poParent->integer_age== 1&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_1): ( poParent->integer_age== 2&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_2): ( poParent->integer_age== 3&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_3): ( poParent->integer_age== 4&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_4): ( poParent->integer_age== 0&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_0): ( poParent->integer_age== 1&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_1): ( poParent->integer_age== 2&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_2): ( poParent->integer_age== 3&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_3): ( poParent->integer_age== 4&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_4): 0.0 );
#line 7697 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 75 "../code/FamilyGeneral.mpp"
		poParent->children_in_household = ( poParent->sex== FEMALE )? poParent->sum_over_mlHHMotherChildren_is_alive: poParent->sum_over_mlHHFatherChildren_is_alive;
#line 7700 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 71 "../code/FamilyGeneral.mpp"
		Set_family_role( ( poParent->lSpouse&& poParent->sex== MALE )? FR_SPOUSE: ( poParent->lHHFather|| poParent->lHHMother )? FR_CHILD: FR_HEAD );
#line 7703 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 81 "../code/FamilyPartnerMatching.mpp"
		Set_female_spouse_age( ( poParent->sex== FEMALE )? COERCE(FEMALE_SPOUSE_AGE,poParent->integer_age): MAX(FEMALE_SPOUSE_AGE) );
#line 7706 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 47 "../code/FertilityGeneral.mpp"
		Set_is_fertile( ( poParent->sex== FEMALE&& WITHIN(FERTILE_AGE_RANGE,poParent->integer_age)&& ( poParent->creation_type== CT_BIRTH|| ( poParent->creation_type== CT_START&& poParent->calendar_year>= MIN(SIM_YEAR_RANGE) )|| ( poParent->creation_type== CT_POOL&& poParent->ever_resident )|| ( poParent->creation_type== CT_SCRATCH&& poParent->ever_resident ) ) )? TRUE: FALSE );
#line 7709 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 92 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_lives_with_dependent_child( ( ( poParent->sex== FEMALE&& poParent->count_mlHHMotherChildren> 0 )|| ( poParent->sex== MALE&& poParent->count_mlHHFatherChildren> 0 ) )? TRUE: FALSE );
#line 7712 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 78 "../code/FamilyPartnerMatching.mpp"
		Set_male_age( ( poParent->sex== MALE )? COERCE(MALE_SPOUSE_AGE,poParent->integer_age): MAX(MALE_SPOUSE_AGE) );
#line 7715 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 75 "../code/FamilyPartnerMatching.mpp"
		Set_male_spouse_age( ( poParent->has_spouse&& poParent->sex== FEMALE )? COERCE(MALE_SPOUSE_AGE,poParent->lSpouse_integer_age): MAX(MALE_SPOUSE_AGE) );
#line 7718 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		// actorset dimension or order state
		if (poParent->actorset_filter_asResidentsAge0Sex)
		{
			asResidentsAge0Sex[OldValue]->Remove(poParent->poasResidentsAge0Sex);
			poParent->poasResidentsAge0Sex = asResidentsAge0Sex[poParent->sex]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asResidentsAge0SexGeo)
		{
			asResidentsAge0SexGeo[OldValue][poParent->geo_birth]->Remove(poParent->poasResidentsAge0SexGeo);
			poParent->poasResidentsAge0SexGeo = asResidentsAge0SexGeo[poParent->sex][poParent->geo_birth]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asSimBornAge0)
		{
			asSimBornAge0[OldValue][poParent->geo][poParent->educ1_group]->Remove(poParent->poasSimBornAge0);
			poParent->poasSimBornAge0 = asSimBornAge0[poParent->sex][poParent->geo][poParent->educ1_group]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asWantToMove)
		{
			asWantToMove[poParent->geo][poParent->geo_want_to_move][OldValue][RANGE_POS(AGE_MIG, poParent->age_mig)]->Remove(poParent->poasWantToMove);
			poParent->poasWantToMove = asWantToMove[poParent->geo][poParent->geo_want_to_move][poParent->sex][RANGE_POS(AGE_MIG, poParent->age_mig)]->NewNode(poParent);
		}
		if ( poParent->poEmigrationEvent->cRecalc == 0 ) {
			poParent->poEmigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poEmigrationEvent );
		};
		if ( poParent->poFirstUnionFormationEvent->cRecalc == 0 ) {
			poParent->poFirstUnionFormationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstUnionFormationEvent );
		};
		if ( poParent->poMigrationEvent->cRecalc == 0 ) {
			poParent->poMigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMigrationEvent );
		};
		if ( poParent->poMortalityBaseEvent->cRecalc == 0 ) {
			poParent->poMortalityBaseEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMortalityBaseEvent );
		};
		if ( poParent->poMortalityDetailedEvent->cRecalc == 0 ) {
			poParent->poMortalityDetailedEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMortalityDetailedEvent );
		};
	}
	return poParent->sex = value;
}

counter PersonStateFunctions::Set_table_unit_TabEducFateByGroup( counter value )
{

	if ( poParent->table_unit_TabEducFateByGroup != value )
	{


		counter	OldValue;
		OldValue = poParent->table_unit_TabEducFateByGroup;
		poParent->table_unit_TabEducFateByGroup = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[3] >= 0 && poParent->tableinfo.tab3events == poParent->__events ) {
				poParent->tableinfo.tab3close0 = (real) poParent->table_unit_TabEducFateByGroup;
			}
		}
	}
	return poParent->table_unit_TabEducFateByGroup = value;
}

counter PersonStateFunctions::Set_table_unit_TabEducFateDistrYob( counter value )
{

	if ( poParent->table_unit_TabEducFateDistrYob != value )
	{


		counter	OldValue;
		OldValue = poParent->table_unit_TabEducFateDistrYob;
		poParent->table_unit_TabEducFateDistrYob = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[6] >= 0 && poParent->tableinfo.tab6events == poParent->__events ) {
				poParent->tableinfo.tab6close0 = (real) poParent->table_unit_TabEducFateDistrYob;
			}
		}
	}
	return poParent->table_unit_TabEducFateDistrYob = value;
}

counter PersonStateFunctions::Set_table_unit_TabPrimarySchoolPlanning( counter value )
{

	if ( poParent->table_unit_TabPrimarySchoolPlanning != value )
	{


		counter	OldValue;
		OldValue = poParent->table_unit_TabPrimarySchoolPlanning;
		poParent->table_unit_TabPrimarySchoolPlanning = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[4] >= 0 && poParent->tableinfo.tab4events == poParent->__events ) {
				poParent->tableinfo.tab4close2 = (real) poParent->table_unit_TabPrimarySchoolPlanning;
			}
		}
	}
	return poParent->table_unit_TabPrimarySchoolPlanning = value;
}

counter PersonStateFunctions::Set_table_unit_TabSchool2TrackTab( counter value )
{

	if ( poParent->table_unit_TabSchool2TrackTab != value )
	{


		counter	OldValue;
		OldValue = poParent->table_unit_TabSchool2TrackTab;
		poParent->table_unit_TabSchool2TrackTab = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[12] >= 0 && poParent->tableinfo.tab12events == poParent->__events ) {
				poParent->tableinfo.tab12close0 = (real) poParent->table_unit_TabSchool2TrackTab;
			}
		}
	}
	return poParent->table_unit_TabSchool2TrackTab = value;
}

counter PersonStateFunctions::Set_table_unit_tabStuntingSexRegMotherYob( counter value )
{

	if ( poParent->table_unit_tabStuntingSexRegMotherYob != value )
	{


		counter	OldValue;
		OldValue = poParent->table_unit_tabStuntingSexRegMotherYob;
		poParent->table_unit_tabStuntingSexRegMotherYob = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[21] >= 0 && poParent->tableinfo.tab21events == poParent->__events ) {
				poParent->tableinfo.tab21close1 = (real) poParent->table_unit_tabStuntingSexRegMotherYob;
			}
		}
	}
	return poParent->table_unit_tabStuntingSexRegMotherYob = value;
}

counter PersonStateFunctions::Set_table_unit_tabStuntingSexRegYob( counter value )
{

	if ( poParent->table_unit_tabStuntingSexRegYob != value )
	{


		counter	OldValue;
		OldValue = poParent->table_unit_tabStuntingSexRegYob;
		poParent->table_unit_tabStuntingSexRegYob = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[20] >= 0 && poParent->tableinfo.tab20events == poParent->__events ) {
				poParent->tableinfo.tab20close1 = (real) poParent->table_unit_tabStuntingSexRegYob;
			}
		}
	}
	return poParent->table_unit_tabStuntingSexRegYob = value;
}

logical PersonStateFunctions::Set_this_parity_spell( logical value )
{

	value = value != 0;
	if ( poParent->this_parity_spell != value )
	{
		bool bScheduleTimeChange = false;	//indicates that schedule event should be recalculated
		if (	(poParent->schedule_event_time == poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time))
		{
			//force timeScheduleEvent to calculate min time
			poParent->schedule_event_time = TIME_INFINITE;
		//indicate that schedule event needs to be recalculated
			bScheduleTimeChange = true;
		}
		if ( poParent->this_parity_spell == TRUE ) {
			poDerivedStates->completed_spell_duration_this_parity_spell_TRUE = poDerivedStates->Get_active_spell_duration_this_parity_spell_TRUE();
			poDerivedStates->active_spell_duration_this_parity_spell_TRUE = (TIME) 0;
		}
		else if ( value == TRUE ) {
			poDerivedStates->active_spell_duration_this_parity_spell_TRUE_time = poParent->__time;
		}
		poDerivedStates->Set_sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY( 0 );
		if ( value == TRUE ) {
			poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time = WAIT1(poParent->time, POINTS(DUR_TIME_IN_PARITY)[0] );
			if ( (!bScheduleTimeChange || poParent->schedule_event_time != TIME_INFINITE) && poParent->schedule_event_time > poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time)
			{
				poParent->schedule_event_time = poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time;
				//indicate that schedule event should be recalculated
				bScheduleTimeChange = true;
			}
		}
		else
		{
			poDerivedStates->sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY_time = TIME_INFINITE;
		}


		logical	OldValue;
		OldValue = poParent->this_parity_spell;
		poParent->this_parity_spell = value;


		if (bScheduleTimeChange && poParent->poScheduleEvent->cRecalc == 0)
		{
			poParent->poScheduleEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent(poParent->poScheduleEvent);
		}
	}
	return poParent->this_parity_spell = value;
}

TIME PersonStateFunctions::Set_time_last_move( TIME value )
{

	if ( poParent->time_last_move != value )
	{


		TIME	OldValue;
		OldValue = poParent->time_last_move;
		poParent->time_last_move = value;

		if ( poParent->poMigrationEvent->cRecalc == 0 ) {
			poParent->poMigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMigrationEvent );
		};
	}
	return poParent->time_last_move = value;
}

TIME PersonStateFunctions::Set_time_next_birthday( TIME value )
{

	if ( poParent->time_next_birthday != value )
	{


		TIME	OldValue;
		OldValue = poParent->time_next_birthday;
		poParent->time_next_birthday = value;

		if ( poParent->poBirthdayEvent->cRecalc == 0 ) {
			poParent->poBirthdayEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poBirthdayEvent );
		};
	}
	return poParent->time_next_birthday = value;
}

TIME PersonStateFunctions::Set_time_of_birth( TIME value )
{

	if ( poParent->time_of_birth != value )
	{


		TIME	OldValue;
		OldValue = poParent->time_of_birth;
		poParent->time_of_birth = value;

#line 88 "../code/ActorPerson.mpp"
		Set_year_of_birth( COERCE(ALL_YEAR_RANGE,int(poParent->time_of_birth)) );
#line 7986 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		if ( poParent->poBirthdayEvent->cRecalc == 0 ) {
			poParent->poBirthdayEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poBirthdayEvent );
		};
		if ( poParent->poDeathAtMaxLifespanEvent->cRecalc == 0 ) {
			poParent->poDeathAtMaxLifespanEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poDeathAtMaxLifespanEvent );
		};
		if ( poParent->poLeavingHomeEvent->cRecalc == 0 ) {
			poParent->poLeavingHomeEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poLeavingHomeEvent );
		};
	}
	return poParent->time_of_birth = value;
}

TIME PersonStateFunctions::Set_time_of_first_immigration( TIME value )
{

	if ( poParent->time_of_first_immigration != value )
	{


		TIME	OldValue;
		OldValue = poParent->time_of_first_immigration;
		poParent->time_of_first_immigration = value;

#line 52 "../code/ImmigrationFromScratch.mpp"
		Set_year_of_first_immigration( COERCE(SIM_YEAR_RANGE,int(poParent->time_of_first_immigration)) );
#line 8016 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		if ( poParent->poFirstImmigrationFromPoolEvent->cRecalc == 0 ) {
			poParent->poFirstImmigrationFromPoolEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstImmigrationFromPoolEvent );
		};
		if ( poParent->poFirstImmigrationFromScratchEvent->cRecalc == 0 ) {
			poParent->poFirstImmigrationFromScratchEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstImmigrationFromScratchEvent );
		};
	}
	return poParent->time_of_first_immigration = value;
}

int PersonStateFunctions::Set_years_preschool( int value )
{

	if ( poParent->years_preschool != value )
	{
		if ( poParent->years_preschool == 0 ) {
			poDerivedStates->Set_duration_years_preschool_0( poDerivedStates->Get_duration_years_preschool_0() );
		}
		else if ( value == 0 ) {
			poDerivedStates->duration_years_preschool_0_time = poParent->__time;
		}
		if ( poParent->years_preschool == 1 ) {
			poDerivedStates->Set_duration_years_preschool_1( poDerivedStates->Get_duration_years_preschool_1() );
		}
		else if ( value == 1 ) {
			poDerivedStates->duration_years_preschool_1_time = poParent->__time;
		}
		if ( poParent->years_preschool == 2 ) {
			poDerivedStates->Set_duration_years_preschool_2( poDerivedStates->Get_duration_years_preschool_2() );
		}
		else if ( value == 2 ) {
			poDerivedStates->duration_years_preschool_2_time = poParent->__time;
		}


		int	OldValue;
		OldValue = poParent->years_preschool;
		poParent->years_preschool = value;

#line 85 "../code/HumanCapitalIndex.mpp"
		poParent->years_of_schooling = poParent->years_preschool+ poParent->educ_one_grade_passed+ poParent->educ_two_grade_passed;
#line 8060 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->years_preschool = value;
}

logical PersonStateFunctions::Set_actorset_filter_asAllFertilePersonsForFertilityAlignment( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asAllFertilePersonsForFertilityAlignment != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asAllFertilePersonsForFertilityAlignment;
		poParent->actorset_filter_asAllFertilePersonsForFertilityAlignment = value;

		// actorset filter
		if (poParent->actorset_filter_asAllFertilePersonsForFertilityAlignment)
		{
			poParent->poasAllFertilePersonsForFertilityAlignment = asAllFertilePersonsForFertilityAlignment[poParent->is_resident]->NewNode(poParent);
		}
		else
		{
			asAllFertilePersonsForFertilityAlignment[poParent->is_resident]->Remove(poParent->poasAllFertilePersonsForFertilityAlignment);
			poParent->poasAllFertilePersonsForFertilityAlignment = NULL;
		}
	}
	return poParent->actorset_filter_asAllFertilePersonsForFertilityAlignment = value;
}

logical PersonStateFunctions::Set_actorset_filter_asAllFertilePersonsForFertilityAlignmentByAge( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asAllFertilePersonsForFertilityAlignmentByAge != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asAllFertilePersonsForFertilityAlignmentByAge;
		poParent->actorset_filter_asAllFertilePersonsForFertilityAlignmentByAge = value;

		// actorset filter
		if (poParent->actorset_filter_asAllFertilePersonsForFertilityAlignmentByAge)
		{
			poParent->poasAllFertilePersonsForFertilityAlignmentByAge = asAllFertilePersonsForFertilityAlignmentByAge[RANGE_POS(FERTILE_AGE_RANGE, poParent->fertile_age)][poParent->is_resident]->NewNode(poParent);
		}
		else
		{
			asAllFertilePersonsForFertilityAlignmentByAge[RANGE_POS(FERTILE_AGE_RANGE, poParent->fertile_age)][poParent->is_resident]->Remove(poParent->poasAllFertilePersonsForFertilityAlignmentByAge);
			poParent->poasAllFertilePersonsForFertilityAlignmentByAge = NULL;
		}
	}
	return poParent->actorset_filter_asAllFertilePersonsForFertilityAlignmentByAge = value;
}

logical PersonStateFunctions::Set_actorset_filter_asAllPerson( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asAllPerson != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asAllPerson;
		poParent->actorset_filter_asAllPerson = value;

		// actorset filter
		if (poParent->actorset_filter_asAllPerson)
		{
			poParent->poasAllPerson = asAllPerson->NewNode(poParent);
		}
		else
		{
			asAllPerson->Remove(poParent->poasAllPerson);
			poParent->poasAllPerson = NULL;
		}
	}
	return poParent->actorset_filter_asAllPerson = value;
}

logical PersonStateFunctions::Set_actorset_filter_asAvailableMale( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asAvailableMale != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asAvailableMale;
		poParent->actorset_filter_asAvailableMale = value;

		// actorset filter
		if (poParent->actorset_filter_asAvailableMale)
		{
			poParent->poasAvailableMale = asAvailableMale[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->geo]->NewNode(poParent);
		}
		else
		{
			asAvailableMale[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->geo]->Remove(poParent->poasAvailableMale);
			poParent->poasAvailableMale = NULL;
		}
	}
	return poParent->actorset_filter_asAvailableMale = value;
}

logical PersonStateFunctions::Set_actorset_filter_asAvailableMaleByType( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asAvailableMaleByType != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asAvailableMaleByType;
		poParent->actorset_filter_asAvailableMaleByType = value;

		// actorset filter
		if (poParent->actorset_filter_asAvailableMaleByType)
		{
			poParent->poasAvailableMaleByType = asAvailableMaleByType[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->spouse_group][poParent->geo]->NewNode(poParent);
		}
		else
		{
			asAvailableMaleByType[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->spouse_group][poParent->geo]->Remove(poParent->poasAvailableMaleByType);
			poParent->poasAvailableMaleByType = NULL;
		}
	}
	return poParent->actorset_filter_asAvailableMaleByType = value;
}

logical PersonStateFunctions::Set_actorset_filter_asAvailableMaleByTypeImmiScratch( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asAvailableMaleByTypeImmiScratch != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asAvailableMaleByTypeImmiScratch;
		poParent->actorset_filter_asAvailableMaleByTypeImmiScratch = value;

		// actorset filter
		if (poParent->actorset_filter_asAvailableMaleByTypeImmiScratch)
		{
			poParent->poasAvailableMaleByTypeImmiScratch = asAvailableMaleByTypeImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->spouse_group][poParent->geo][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->NewNode(poParent);
		}
		else
		{
			asAvailableMaleByTypeImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->spouse_group][poParent->geo][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->Remove(poParent->poasAvailableMaleByTypeImmiScratch);
			poParent->poasAvailableMaleByTypeImmiScratch = NULL;
		}
	}
	return poParent->actorset_filter_asAvailableMaleByTypeImmiScratch = value;
}

logical PersonStateFunctions::Set_actorset_filter_asAvailableMaleImmiScratch( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asAvailableMaleImmiScratch != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asAvailableMaleImmiScratch;
		poParent->actorset_filter_asAvailableMaleImmiScratch = value;

		// actorset filter
		if (poParent->actorset_filter_asAvailableMaleImmiScratch)
		{
			poParent->poasAvailableMaleImmiScratch = asAvailableMaleImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->geo][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->NewNode(poParent);
		}
		else
		{
			asAvailableMaleImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->geo][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->Remove(poParent->poasAvailableMaleImmiScratch);
			poParent->poasAvailableMaleImmiScratch = NULL;
		}
	}
	return poParent->actorset_filter_asAvailableMaleImmiScratch = value;
}

logical PersonStateFunctions::Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAge( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAge != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAge;
		poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAge = value;

		// actorset filter
		if (poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAge)
		{
			poParent->poasFemaleInUnionByAgeAndPartnerAge = asFemaleInUnionByAgeAndPartnerAge[RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][RANGE_POS(MALE_SPOUSE_AGE, poParent->male_spouse_age)]->NewNode(poParent);
		}
		else
		{
			asFemaleInUnionByAgeAndPartnerAge[RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][RANGE_POS(MALE_SPOUSE_AGE, poParent->male_spouse_age)]->Remove(poParent->poasFemaleInUnionByAgeAndPartnerAge);
			poParent->poasFemaleInUnionByAgeAndPartnerAge = NULL;
		}
	}
	return poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAge = value;
}

logical PersonStateFunctions::Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch;
		poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch = value;

		// actorset filter
		if (poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch)
		{
			poParent->poasFemaleInUnionByAgeAndPartnerAgeImmiScratch = asFemaleInUnionByAgeAndPartnerAgeImmiScratch[RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][RANGE_POS(MALE_SPOUSE_AGE, poParent->male_spouse_age)][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->NewNode(poParent);
		}
		else
		{
			asFemaleInUnionByAgeAndPartnerAgeImmiScratch[RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][RANGE_POS(MALE_SPOUSE_AGE, poParent->male_spouse_age)][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->Remove(poParent->poasFemaleInUnionByAgeAndPartnerAgeImmiScratch);
			poParent->poasFemaleInUnionByAgeAndPartnerAgeImmiScratch = NULL;
		}
	}
	return poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch = value;
}

logical PersonStateFunctions::Set_actorset_filter_asPotentialImmigrantMothers( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asPotentialImmigrantMothers != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asPotentialImmigrantMothers;
		poParent->actorset_filter_asPotentialImmigrantMothers = value;

		// actorset filter
		if (poParent->actorset_filter_asPotentialImmigrantMothers)
		{
			poParent->poasPotentialImmigrantMothers = asPotentialImmigrantMothers[RANGE_POS(AGE_RANGE, poParent->integer_age)][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)][poParent->geo]->NewNode(poParent);
		}
		else
		{
			asPotentialImmigrantMothers[RANGE_POS(AGE_RANGE, poParent->integer_age)][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)][poParent->geo]->Remove(poParent->poasPotentialImmigrantMothers);
			poParent->poasPotentialImmigrantMothers = NULL;
		}
	}
	return poParent->actorset_filter_asPotentialImmigrantMothers = value;
}

logical PersonStateFunctions::Set_actorset_filter_asPotentialSchoolOneStudent( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asPotentialSchoolOneStudent != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asPotentialSchoolOneStudent;
		poParent->actorset_filter_asPotentialSchoolOneStudent = value;

		// actorset filter
		if (poParent->actorset_filter_asPotentialSchoolOneStudent)
		{
			poParent->poasPotentialSchoolOneStudent = asPotentialSchoolOneStudent->NewNode(poParent);
		}
		else
		{
			asPotentialSchoolOneStudent->Remove(poParent->poasPotentialSchoolOneStudent);
			poParent->poasPotentialSchoolOneStudent = NULL;
		}
	}
	return poParent->actorset_filter_asPotentialSchoolOneStudent = value;
}

logical PersonStateFunctions::Set_actorset_filter_asPotentialSchoolOneStudentToProcess( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asPotentialSchoolOneStudentToProcess != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asPotentialSchoolOneStudentToProcess;
		poParent->actorset_filter_asPotentialSchoolOneStudentToProcess = value;

		// actorset filter
		if (poParent->actorset_filter_asPotentialSchoolOneStudentToProcess)
		{
			poParent->poasPotentialSchoolOneStudentToProcess = asPotentialSchoolOneStudentToProcess->NewNode(poParent);
		}
		else
		{
			asPotentialSchoolOneStudentToProcess->Remove(poParent->poasPotentialSchoolOneStudentToProcess);
			poParent->poasPotentialSchoolOneStudentToProcess = NULL;
		}
	}
	return poParent->actorset_filter_asPotentialSchoolOneStudentToProcess = value;
}

logical PersonStateFunctions::Set_actorset_filter_asPotentialSchoolTwoStudent( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asPotentialSchoolTwoStudent != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asPotentialSchoolTwoStudent;
		poParent->actorset_filter_asPotentialSchoolTwoStudent = value;

		// actorset filter
		if (poParent->actorset_filter_asPotentialSchoolTwoStudent)
		{
			poParent->poasPotentialSchoolTwoStudent = asPotentialSchoolTwoStudent->NewNode(poParent);
		}
		else
		{
			asPotentialSchoolTwoStudent->Remove(poParent->poasPotentialSchoolTwoStudent);
			poParent->poasPotentialSchoolTwoStudent = NULL;
		}
	}
	return poParent->actorset_filter_asPotentialSchoolTwoStudent = value;
}

logical PersonStateFunctions::Set_actorset_filter_asPotentialSchoolTwoStudentToProcess( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asPotentialSchoolTwoStudentToProcess != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asPotentialSchoolTwoStudentToProcess;
		poParent->actorset_filter_asPotentialSchoolTwoStudentToProcess = value;

		// actorset filter
		if (poParent->actorset_filter_asPotentialSchoolTwoStudentToProcess)
		{
			poParent->poasPotentialSchoolTwoStudentToProcess = asPotentialSchoolTwoStudentToProcess->NewNode(poParent);
		}
		else
		{
			asPotentialSchoolTwoStudentToProcess->Remove(poParent->poasPotentialSchoolTwoStudentToProcess);
			poParent->poasPotentialSchoolTwoStudentToProcess = NULL;
		}
	}
	return poParent->actorset_filter_asPotentialSchoolTwoStudentToProcess = value;
}

logical PersonStateFunctions::Set_actorset_filter_asResidentBabies( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asResidentBabies != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asResidentBabies;
		poParent->actorset_filter_asResidentBabies = value;

		// actorset filter
		if (poParent->actorset_filter_asResidentBabies)
		{
			poParent->poasResidentBabies = asResidentBabies->NewNode(poParent);
		}
		else
		{
			asResidentBabies->Remove(poParent->poasResidentBabies);
			poParent->poasResidentBabies = NULL;
		}
	}
	return poParent->actorset_filter_asResidentBabies = value;
}

logical PersonStateFunctions::Set_actorset_filter_asResidentsAge0Sex( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asResidentsAge0Sex != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asResidentsAge0Sex;
		poParent->actorset_filter_asResidentsAge0Sex = value;

		// actorset filter
		if (poParent->actorset_filter_asResidentsAge0Sex)
		{
			poParent->poasResidentsAge0Sex = asResidentsAge0Sex[poParent->sex]->NewNode(poParent);
		}
		else
		{
			asResidentsAge0Sex[poParent->sex]->Remove(poParent->poasResidentsAge0Sex);
			poParent->poasResidentsAge0Sex = NULL;
		}
	}
	return poParent->actorset_filter_asResidentsAge0Sex = value;
}

logical PersonStateFunctions::Set_actorset_filter_asResidentsAge0SexGeo( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asResidentsAge0SexGeo != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asResidentsAge0SexGeo;
		poParent->actorset_filter_asResidentsAge0SexGeo = value;

		// actorset filter
		if (poParent->actorset_filter_asResidentsAge0SexGeo)
		{
			poParent->poasResidentsAge0SexGeo = asResidentsAge0SexGeo[poParent->sex][poParent->geo_birth]->NewNode(poParent);
		}
		else
		{
			asResidentsAge0SexGeo[poParent->sex][poParent->geo_birth]->Remove(poParent->poasResidentsAge0SexGeo);
			poParent->poasResidentsAge0SexGeo = NULL;
		}
	}
	return poParent->actorset_filter_asResidentsAge0SexGeo = value;
}

logical PersonStateFunctions::Set_actorset_filter_asSimBornAge0( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asSimBornAge0 != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asSimBornAge0;
		poParent->actorset_filter_asSimBornAge0 = value;

		// actorset filter
		if (poParent->actorset_filter_asSimBornAge0)
		{
			poParent->poasSimBornAge0 = asSimBornAge0[poParent->sex][poParent->geo][poParent->educ1_group]->NewNode(poParent);
		}
		else
		{
			asSimBornAge0[poParent->sex][poParent->geo][poParent->educ1_group]->Remove(poParent->poasSimBornAge0);
			poParent->poasSimBornAge0 = NULL;
		}
	}
	return poParent->actorset_filter_asSimBornAge0 = value;
}

logical PersonStateFunctions::Set_actorset_filter_asWantToMove( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asWantToMove != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asWantToMove;
		poParent->actorset_filter_asWantToMove = value;

		// actorset filter
		if (poParent->actorset_filter_asWantToMove)
		{
			poParent->poasWantToMove = asWantToMove[poParent->geo][poParent->geo_want_to_move][poParent->sex][RANGE_POS(AGE_MIG, poParent->age_mig)]->NewNode(poParent);
		}
		else
		{
			asWantToMove[poParent->geo][poParent->geo_want_to_move][poParent->sex][RANGE_POS(AGE_MIG, poParent->age_mig)]->Remove(poParent->poasWantToMove);
			poParent->poasWantToMove = NULL;
		}
	}
	return poParent->actorset_filter_asWantToMove = value;
}

logical PersonStateFunctions::Set_actorset_filter_asWomenNoChildren( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asWomenNoChildren != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asWomenNoChildren;
		poParent->actorset_filter_asWomenNoChildren = value;

		// actorset filter
		if (poParent->actorset_filter_asWomenNoChildren)
		{
			poParent->poasWomenNoChildren = asWomenNoChildren[poParent->union1_group][RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][poParent->has_spouse]->NewNode(poParent);
		}
		else
		{
			asWomenNoChildren[poParent->union1_group][RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][poParent->has_spouse]->Remove(poParent->poasWomenNoChildren);
			poParent->poasWomenNoChildren = NULL;
		}
	}
	return poParent->actorset_filter_asWomenNoChildren = value;
}

logical PersonStateFunctions::Set_actorset_filter_asWomenWithChildren( logical value )
{

	value = value != 0;
	if ( poParent->actorset_filter_asWomenWithChildren != value )
	{


		logical	OldValue;
		OldValue = poParent->actorset_filter_asWomenWithChildren;
		poParent->actorset_filter_asWomenWithChildren = value;

		// actorset filter
		if (poParent->actorset_filter_asWomenWithChildren)
		{
			poParent->poasWomenWithChildren = asWomenWithChildren[poParent->union1_group][poParent->child_agegr][poParent->moth_agegr][poParent->has_spouse]->NewNode(poParent);
		}
		else
		{
			asWomenWithChildren[poParent->union1_group][poParent->child_agegr][poParent->moth_agegr][poParent->has_spouse]->Remove(poParent->poasWomenWithChildren);
			poParent->poasWomenWithChildren = NULL;
		}
	}
	return poParent->actorset_filter_asWomenWithChildren = value;
}

double PersonStateFunctions::Set_age_last_birth( double value )
{

	if ( poParent->age_last_birth != value )
	{


		double	OldValue;
		OldValue = poParent->age_last_birth;
		poParent->age_last_birth = value;

		{
			int	nNewGrp;

			nNewGrp = Split( (double) poParent->age_last_birth, poDerivedStates->split_age_last_birth_MOTH_AGEGR_PART, 5, p_points_MOTH_AGEGR_PART );
			if ( nNewGrp != poDerivedStates->split_age_last_birth_MOTH_AGEGR_PART ) {
				poDerivedStates->Set_split_age_last_birth_MOTH_AGEGR_PART( nNewGrp );
			}
		}
	}
	return poParent->age_last_birth = value;
}

AGE_MIG PersonStateFunctions::Set_age_mig( AGE_MIG value )
{

	CHKLMT(value, age_mig, r_min_AGE_MIG, r_max_AGE_MIG);
	if ( value < r_min_AGE_MIG ) {
		value = r_min_AGE_MIG;
	}
	else if ( value > r_max_AGE_MIG ) {
		value = r_max_AGE_MIG;
	}
	if ( poParent->age_mig != value )
	{


		AGE_MIG	OldValue;
		OldValue = poParent->age_mig;
		poParent->age_mig = value;

		// actorset dimension or order state
		if (poParent->actorset_filter_asWantToMove)
		{
			asWantToMove[poParent->geo][poParent->geo_want_to_move][poParent->sex][RANGE_POS(AGE_MIG,OldValue)]->Remove(poParent->poasWantToMove);
			poParent->poasWantToMove = asWantToMove[poParent->geo][poParent->geo_want_to_move][poParent->sex][RANGE_POS(AGE_MIG, poParent->age_mig)]->NewNode(poParent);
		}
		if ( poParent->poEmigrationEvent->cRecalc == 0 ) {
			poParent->poEmigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poEmigrationEvent );
		};
		if ( poParent->poMigrationEvent->cRecalc == 0 ) {
			poParent->poMigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMigrationEvent );
		};
	}
	return poParent->age_mig = value;
}

double PersonStateFunctions::Set_age_youngest_child( double value )
{

	if ( poParent->age_youngest_child != value )
	{


		double	OldValue;
		OldValue = poParent->age_youngest_child;
		poParent->age_youngest_child = value;

#line 96 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_age_last_birth( ( poParent->lives_with_dependent_child )? poParent->integer_age- poParent->age_youngest_child: TIME_INFINITE );
#line 8683 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		{
			int	nNewGrp;

			nNewGrp = Split( (double) poParent->age_youngest_child, poDerivedStates->split_age_youngest_child_CHILD_AGEGR_PART, 6, p_points_CHILD_AGEGR_PART );
			if ( nNewGrp != poDerivedStates->split_age_youngest_child_CHILD_AGEGR_PART ) {
				poDerivedStates->Set_split_age_youngest_child_CHILD_AGEGR_PART( nNewGrp );
			}
		}
	}
	return poParent->age_youngest_child = value;
}

ALL_YEAR_RANGE PersonStateFunctions::Set_all_year( ALL_YEAR_RANGE value )
{

	CHKLMT(value, all_year, r_min_ALL_YEAR_RANGE, r_max_ALL_YEAR_RANGE);
	if ( value < r_min_ALL_YEAR_RANGE ) {
		value = r_min_ALL_YEAR_RANGE;
	}
	else if ( value > r_max_ALL_YEAR_RANGE ) {
		value = r_max_ALL_YEAR_RANGE;
	}
	if ( poParent->all_year != value )
	{


		ALL_YEAR_RANGE	OldValue;
		OldValue = poParent->all_year;
		poParent->all_year = value;

	}
	return poParent->all_year = value;
}

CHILD_AGEGR PersonStateFunctions::Set_child_agegr( CHILD_AGEGR value )
{

	CHKLMT(value, child_agegr, 0, rc_size_CHILD_AGEGR - 1);
	if ( value < 0 ) {
		value = (CHILD_AGEGR) 0;
	}
	else if ( value >= rc_size_CHILD_AGEGR ) {
		value = (CHILD_AGEGR) ( rc_size_CHILD_AGEGR - 1 );
	}
	if ( poParent->child_agegr != value )
	{


		CHILD_AGEGR	OldValue;
		OldValue = poParent->child_agegr;
		poParent->child_agegr = value;

		// actorset dimension or order state
		if (poParent->actorset_filter_asWomenWithChildren)
		{
			asWomenWithChildren[poParent->union1_group][OldValue][poParent->moth_agegr][poParent->has_spouse]->Remove(poParent->poasWomenWithChildren);
			poParent->poasWomenWithChildren = asWomenWithChildren[poParent->union1_group][poParent->child_agegr][poParent->moth_agegr][poParent->has_spouse]->NewNode(poParent);
		}
	}
	return poParent->child_agegr = value;
}

int PersonStateFunctions::Set_child_agegr_part( int value )
{

	if ( poParent->child_agegr_part != value )
	{


		int	OldValue;
		OldValue = poParent->child_agegr_part;
		poParent->child_agegr_part = value;

#line 114 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_child_agegr( ( poParent->child_agegr_part== 0 )? CA00: ( poParent->child_agegr_part== 1 )? CA01: ( poParent->child_agegr_part== 2 )? CA03: ( poParent->child_agegr_part== 3 )? CA06: ( poParent->child_agegr_part== 4 )? CA09: ( poParent->child_agegr_part== 5 )? CA12: CA15 );
#line 8759 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->child_agegr_part = value;
}

double PersonStateFunctions::Set_child_mortality( double value )
{

	if ( poParent->child_mortality != value )
	{


		double	OldValue;
		OldValue = poParent->child_mortality;
		poParent->child_mortality = value;

		if ( poParent->poMortalityDetailedEvent->cRecalc == 0 ) {
			poParent->poMortalityDetailedEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMortalityDetailedEvent );
		};
	}
	return poParent->child_mortality = value;
}

CHILD_MORTALITY_GROUP PersonStateFunctions::Set_child_mortality_group( CHILD_MORTALITY_GROUP value )
{

	CHKLMT(value, child_mortality_group, 0, rc_size_CHILD_MORTALITY_GROUP - 1);
	if ( value < 0 ) {
		value = (CHILD_MORTALITY_GROUP) 0;
	}
	else if ( value >= rc_size_CHILD_MORTALITY_GROUP ) {
		value = (CHILD_MORTALITY_GROUP) ( rc_size_CHILD_MORTALITY_GROUP - 1 );
	}
	if ( poParent->child_mortality_group != value )
	{


		CHILD_MORTALITY_GROUP	OldValue;
		OldValue = poParent->child_mortality_group;
		poParent->child_mortality_group = value;

		if ( poParent->poMortalityDetailedEvent->cRecalc == 0 ) {
			poParent->poMortalityDetailedEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMortalityDetailedEvent );
		};
	}
	return poParent->child_mortality_group = value;
}

EDUC1_GROUP PersonStateFunctions::Set_educ1_group( EDUC1_GROUP value )
{

	CHKLMT(value, educ1_group, 0, rc_size_EDUC1_GROUP - 1);
	if ( value < 0 ) {
		value = (EDUC1_GROUP) 0;
	}
	else if ( value >= rc_size_EDUC1_GROUP ) {
		value = (EDUC1_GROUP) ( rc_size_EDUC1_GROUP - 1 );
	}
	if ( poParent->educ1_group != value )
	{


		EDUC1_GROUP	OldValue;
		OldValue = poParent->educ1_group;
		poParent->educ1_group = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[3] >= 0 && ( poParent->tableinfo.tab3dim_filter || poParent->table_filter_TabEducFateByGroup ) ) {
				if ( poParent->tableinfo.tab3events == -1 || poParent->tableinfo.tab3events < poParent->__events ) {
					if ( poParent->tableinfo.tab3events != -1 ) {
						Set_table_unit_TabEducFateByGroup(poParent->table_unit_TabEducFateByGroup+1);
						poParent->tableinfo.tab3close0 = (real) poParent->table_unit_TabEducFateByGroup;
						poParent->tableinfo.tab3open0 = (real) poParent->CloseUpdate( 3, 0, poParent->tableinfo.tab3open0, poParent->tableinfo.tab3close0, poParent->tableinfo.tab3dim_filter, poParent->tableinfo.tab3cell, poParent->tableinfo.tab3total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab3events = poParent->__events;
					poParent->tableinfo.tab3close0 = (real) poParent->table_unit_TabEducFateByGroup;
					poParent->tableinfo.tab3cell = (index) poParent->TableCellIndex( 3, 3,
						OldValue,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE,
						poParent->educ_one_fate );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab3total_cell = (index) poParent->tableinfo.tab3cell;
					}
					else {
						poParent->tableinfo.tab3total_cell = (index) poParent->TableTotalCellIndex( 3, 3,
						OldValue,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE,
						poParent->educ_one_fate );
					}
					poParent->tableinfo.tab3dim_filter = poParent->table_filter_TabEducFateByGroup;
				}
			}
		}
		// actorset dimension or order state
		if (poParent->actorset_filter_asSimBornAge0)
		{
			asSimBornAge0[poParent->sex][poParent->geo][OldValue]->Remove(poParent->poasSimBornAge0);
			poParent->poasSimBornAge0 = asSimBornAge0[poParent->sex][poParent->geo][poParent->educ1_group]->NewNode(poParent);
		}
	}
	return poParent->educ1_group = value;
}

EDUC_ONE_GRADE PersonStateFunctions::Set_educ_one_grade( EDUC_ONE_GRADE value )
{

	CHKLMT(value, educ_one_grade, r_min_EDUC_ONE_GRADE, r_max_EDUC_ONE_GRADE);
	if ( value < r_min_EDUC_ONE_GRADE ) {
		value = r_min_EDUC_ONE_GRADE;
	}
	else if ( value > r_max_EDUC_ONE_GRADE ) {
		value = r_max_EDUC_ONE_GRADE;
	}
	if ( poParent->educ_one_grade != value )
	{


		EDUC_ONE_GRADE	OldValue;
		OldValue = poParent->educ_one_grade;
		poParent->educ_one_grade = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[4] >= 0 && ( poParent->tableinfo.tab4dim_filter || poParent->table_filter_TabPrimarySchoolPlanning ) ) {
				if ( poParent->tableinfo.tab4events == -1 || poParent->tableinfo.tab4events < poParent->__events ) {
					if ( poParent->tableinfo.tab4events != -1 ) {
						poParent->tableinfo.tab4open0 = (real) poParent->CloseUpdate( 4, 0, poParent->tableinfo.tab4open0, poParent->tableinfo.tab4close0, poParent->tableinfo.tab4dim_filter, poParent->tableinfo.tab4cell, poParent->tableinfo.tab4total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab4open1 = (real) poParent->CloseUpdate( 4, 1, poParent->tableinfo.tab4open1, poParent->tableinfo.tab4close1, poParent->tableinfo.tab4dim_filter, poParent->tableinfo.tab4cell, poParent->tableinfo.tab4total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						Set_table_unit_TabPrimarySchoolPlanning(poParent->table_unit_TabPrimarySchoolPlanning+1);
						poParent->tableinfo.tab4close2 = (real) poParent->table_unit_TabPrimarySchoolPlanning;
						poParent->tableinfo.tab4open2 = (real) poParent->CloseUpdate( 4, 2, poParent->tableinfo.tab4open2, poParent->tableinfo.tab4close2, poParent->tableinfo.tab4dim_filter, poParent->tableinfo.tab4cell, poParent->tableinfo.tab4total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab4events = poParent->__events;
					poParent->tableinfo.tab4close0 = (real) poParent->educ_one_required_rooms_value_out;
					poParent->tableinfo.tab4close1 = (real) poParent->educ_one_required_teachers_value_out;
					poParent->tableinfo.tab4close2 = (real) poParent->table_unit_TabPrimarySchoolPlanning;
					poParent->tableinfo.tab4cell = (index) poParent->TableCellIndex( 4, 3,
						poParent->geo_nat,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						OldValue - r_min_EDUC_ONE_GRADE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab4total_cell = (index) poParent->tableinfo.tab4cell;
					}
					else {
						poParent->tableinfo.tab4total_cell = (index) poParent->TableTotalCellIndex( 4, 3,
						poParent->geo_nat,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						OldValue - r_min_EDUC_ONE_GRADE );
					}
					poParent->tableinfo.tab4dim_filter = poParent->table_filter_TabPrimarySchoolPlanning;
				}
			}
		}
	}
	return poParent->educ_one_grade = value;
}

EDUC_ONE_LEVEL PersonStateFunctions::Set_educ_one_level( EDUC_ONE_LEVEL value )
{

	CHKLMT(value, educ_one_level, 0, rc_size_EDUC_ONE_LEVEL - 1);
	if ( value < 0 ) {
		value = (EDUC_ONE_LEVEL) 0;
	}
	else if ( value >= rc_size_EDUC_ONE_LEVEL ) {
		value = (EDUC_ONE_LEVEL) ( rc_size_EDUC_ONE_LEVEL - 1 );
	}
	if ( poParent->educ_one_level != value )
	{
		if ( poParent->educ_one_level == EOL_HIGH ) {
			poDerivedStates->Set_duration_educ_one_level_EOL_HIGH( poDerivedStates->Get_duration_educ_one_level_EOL_HIGH() );
		}
		else if ( value == EOL_HIGH ) {
			poDerivedStates->duration_educ_one_level_EOL_HIGH_time = poParent->__time;
		}
		if ( poParent->educ_one_level == EOL_LOW ) {
			poDerivedStates->Set_duration_educ_one_level_EOL_LOW( poDerivedStates->Get_duration_educ_one_level_EOL_LOW() );
		}
		else if ( value == EOL_LOW ) {
			poDerivedStates->duration_educ_one_level_EOL_LOW_time = poParent->__time;
		}
		if ( poParent->educ_one_level == EOL_MEDIUM ) {
			poDerivedStates->Set_duration_educ_one_level_EOL_MEDIUM( poDerivedStates->Get_duration_educ_one_level_EOL_MEDIUM() );
		}
		else if ( value == EOL_MEDIUM ) {
			poDerivedStates->duration_educ_one_level_EOL_MEDIUM_time = poParent->__time;
		}


		EDUC_ONE_LEVEL	OldValue;
		OldValue = poParent->educ_one_level;
		poParent->educ_one_level = value;

	}
	return poParent->educ_one_level = value;
}

double PersonStateFunctions::Set_educ_one_required_rooms( double value )
{

	if ( poParent->educ_one_required_rooms != value )
	{


		double	OldValue;
		OldValue = poParent->educ_one_required_rooms;
		poParent->educ_one_required_rooms = value;

#line 135 "../code/TablesEducation.mpp"
		Set_educ_one_required_rooms_value_out( poParent->educ_one_required_rooms );
#line 8970 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->educ_one_required_rooms = value;
}

double PersonStateFunctions::Set_educ_one_required_rooms_value_out( double value )
{

	if ( poParent->educ_one_required_rooms_value_out != value )
	{


		double	OldValue;
		OldValue = poParent->educ_one_required_rooms_value_out;
		poParent->educ_one_required_rooms_value_out = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[4] >= 0 && poParent->tableinfo.tab4events == poParent->__events ) {
				poParent->tableinfo.tab4close0 = (real) poParent->educ_one_required_rooms_value_out;
			}
		}
	}
	return poParent->educ_one_required_rooms_value_out = value;
}

double PersonStateFunctions::Set_educ_one_required_teachers( double value )
{

	if ( poParent->educ_one_required_teachers != value )
	{


		double	OldValue;
		OldValue = poParent->educ_one_required_teachers;
		poParent->educ_one_required_teachers = value;

#line 135 "../code/TablesEducation.mpp"
		Set_educ_one_required_teachers_value_out( poParent->educ_one_required_teachers );
#line 9008 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->educ_one_required_teachers = value;
}

double PersonStateFunctions::Set_educ_one_required_teachers_value_out( double value )
{

	if ( poParent->educ_one_required_teachers_value_out != value )
	{


		double	OldValue;
		OldValue = poParent->educ_one_required_teachers_value_out;
		poParent->educ_one_required_teachers_value_out = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[4] >= 0 && poParent->tableinfo.tab4events == poParent->__events ) {
				poParent->tableinfo.tab4close1 = (real) poParent->educ_one_required_teachers_value_out;
			}
		}
	}
	return poParent->educ_one_required_teachers_value_out = value;
}

EDUC_TWO_LEVEL PersonStateFunctions::Set_educ_two_level( EDUC_TWO_LEVEL value )
{

	CHKLMT(value, educ_two_level, 0, rc_size_EDUC_TWO_LEVEL - 1);
	if ( value < 0 ) {
		value = (EDUC_TWO_LEVEL) 0;
	}
	else if ( value >= rc_size_EDUC_TWO_LEVEL ) {
		value = (EDUC_TWO_LEVEL) ( rc_size_EDUC_TWO_LEVEL - 1 );
	}
	if ( poParent->educ_two_level != value )
	{
		if ( poParent->educ_two_level == ETL_GRADUATE ) {
			poDerivedStates->Set_duration_educ_two_level_ETL_GRADUATE( poDerivedStates->Get_duration_educ_two_level_ETL_GRADUATE() );
		}
		else if ( value == ETL_GRADUATE ) {
			poDerivedStates->duration_educ_two_level_ETL_GRADUATE_time = poParent->__time;
		}
		if ( poParent->educ_two_level == ETL_INCOMPLETE ) {
			poDerivedStates->Set_duration_educ_two_level_ETL_INCOMPLETE( poDerivedStates->Get_duration_educ_two_level_ETL_INCOMPLETE() );
		}
		else if ( value == ETL_INCOMPLETE ) {
			poDerivedStates->duration_educ_two_level_ETL_INCOMPLETE_time = poParent->__time;
		}
		if ( poParent->educ_two_level == ETL_NEVER ) {
			poDerivedStates->Set_duration_educ_two_level_ETL_NEVER( poDerivedStates->Get_duration_educ_two_level_ETL_NEVER() );
		}
		else if ( value == ETL_NEVER ) {
			poDerivedStates->duration_educ_two_level_ETL_NEVER_time = poParent->__time;
		}


		EDUC_TWO_LEVEL	OldValue;
		OldValue = poParent->educ_two_level;
		poParent->educ_two_level = value;

	}
	return poParent->educ_two_level = value;
}

ETHNICITY_SHORT PersonStateFunctions::Set_ethnicity_short( ETHNICITY_SHORT value )
{

	CHKLMT(value, ethnicity_short, 0, rc_size_ETHNICITY_SHORT - 1);
	if ( value < 0 ) {
		value = (ETHNICITY_SHORT) 0;
	}
	else if ( value >= rc_size_ETHNICITY_SHORT ) {
		value = (ETHNICITY_SHORT) ( rc_size_ETHNICITY_SHORT - 1 );
	}
	if ( poParent->ethnicity_short != value )
	{


		ETHNICITY_SHORT	OldValue;
		OldValue = poParent->ethnicity_short;
		poParent->ethnicity_short = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[23] >= 0 && ( poParent->tableinfo.tab23dim_filter || poParent->table_filter_TabImmunizationChildren ) ) {
				if ( poParent->tableinfo.tab23events == -1 || poParent->tableinfo.tab23events < poParent->__events ) {
					if ( poParent->tableinfo.tab23events != -1 ) {
						poParent->tableinfo.tab23open0 = (real) poParent->CloseUpdate( 23, 0, poParent->tableinfo.tab23open0, poParent->tableinfo.tab23close0, poParent->tableinfo.tab23dim_filter, poParent->tableinfo.tab23cell, poParent->tableinfo.tab23total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab23open1 = (real) poParent->CloseUpdate( 23, 1, poParent->tableinfo.tab23open1, poParent->tableinfo.tab23close1, poParent->tableinfo.tab23dim_filter, poParent->tableinfo.tab23cell, poParent->tableinfo.tab23total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab23events = poParent->__events;
					poParent->tableinfo.tab23close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab23close1 = (real) poDerivedStates->Get_duration_is_immunized_TRUE();
					poParent->tableinfo.tab23cell = (index) poParent->TableCellIndex( 23, 6,
						poParent->sex,
						OldValue,
						poParent->got_prenat_care,
						poParent->region_nat,
						poParent->educ_mother,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab23total_cell = (index) poParent->tableinfo.tab23cell;
					}
					else {
						poParent->tableinfo.tab23total_cell = (index) poParent->TableTotalCellIndex( 23, 6,
						poParent->sex,
						OldValue,
						poParent->got_prenat_care,
						poParent->region_nat,
						poParent->educ_mother,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab23dim_filter = poParent->table_filter_TabImmunizationChildren;
				}
			}
		}
	}
	return poParent->ethnicity_short = value;
}

FAM_ROLE PersonStateFunctions::Set_family_role( FAM_ROLE value )
{

	CHKLMT(value, family_role, 0, rc_size_FAM_ROLE - 1);
	if ( value < 0 ) {
		value = (FAM_ROLE) 0;
	}
	else if ( value >= rc_size_FAM_ROLE ) {
		value = (FAM_ROLE) ( rc_size_FAM_ROLE - 1 );
	}
	if ( poParent->family_role != value )
	{
#line 10 "../code/ImmigrationFromScratch.mpp"
		Set_actorset_filter_asPotentialImmigrantMothers( poParent->is_alive&& poParent->sex== FEMALE&& poParent->creation_type== CT_SCRATCH&& !( poParent->ever_resident )&& value!= FR_CHILD&& poParent->ever_union );
#line 9142 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		FAM_ROLE	OldValue;
		OldValue = poParent->family_role;
		poParent->family_role = value;

		if ( poParent->poFirstImmigrationFromScratchEvent->cRecalc == 0 ) {
			poParent->poFirstImmigrationFromScratchEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstImmigrationFromScratchEvent );
		};
		if ( poParent->poLeavingHomeEvent->cRecalc == 0 ) {
			poParent->poLeavingHomeEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poLeavingHomeEvent );
		};
	}
	return poParent->family_role = value;
}

FEMALE_SPOUSE_AGE PersonStateFunctions::Set_female_spouse_age( FEMALE_SPOUSE_AGE value )
{

	CHKLMT(value, female_spouse_age, r_min_FEMALE_SPOUSE_AGE, r_max_FEMALE_SPOUSE_AGE);
	if ( value < r_min_FEMALE_SPOUSE_AGE ) {
		value = r_min_FEMALE_SPOUSE_AGE;
	}
	else if ( value > r_max_FEMALE_SPOUSE_AGE ) {
		value = r_max_FEMALE_SPOUSE_AGE;
	}
	if ( poParent->female_spouse_age != value )
	{


		FEMALE_SPOUSE_AGE	OldValue;
		OldValue = poParent->female_spouse_age;
		poParent->female_spouse_age = value;

		// actorset dimension or order state
		if (poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAge)
		{
			asFemaleInUnionByAgeAndPartnerAge[RANGE_POS(FEMALE_SPOUSE_AGE,OldValue)][RANGE_POS(MALE_SPOUSE_AGE, poParent->male_spouse_age)]->Remove(poParent->poasFemaleInUnionByAgeAndPartnerAge);
			poParent->poasFemaleInUnionByAgeAndPartnerAge = asFemaleInUnionByAgeAndPartnerAge[RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][RANGE_POS(MALE_SPOUSE_AGE, poParent->male_spouse_age)]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch)
		{
			asFemaleInUnionByAgeAndPartnerAgeImmiScratch[RANGE_POS(FEMALE_SPOUSE_AGE,OldValue)][RANGE_POS(MALE_SPOUSE_AGE, poParent->male_spouse_age)][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->Remove(poParent->poasFemaleInUnionByAgeAndPartnerAgeImmiScratch);
			poParent->poasFemaleInUnionByAgeAndPartnerAgeImmiScratch = asFemaleInUnionByAgeAndPartnerAgeImmiScratch[RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][RANGE_POS(MALE_SPOUSE_AGE, poParent->male_spouse_age)][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asWomenNoChildren)
		{
			asWomenNoChildren[poParent->union1_group][RANGE_POS(FEMALE_SPOUSE_AGE,OldValue)][poParent->has_spouse]->Remove(poParent->poasWomenNoChildren);
			poParent->poasWomenNoChildren = asWomenNoChildren[poParent->union1_group][RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][poParent->has_spouse]->NewNode(poParent);
		}
	}
	return poParent->female_spouse_age = value;
}

FERTILE_AGE_RANGE PersonStateFunctions::Set_fertile_age( FERTILE_AGE_RANGE value )
{

	CHKLMT(value, fertile_age, r_min_FERTILE_AGE_RANGE, r_max_FERTILE_AGE_RANGE);
	if ( value < r_min_FERTILE_AGE_RANGE ) {
		value = r_min_FERTILE_AGE_RANGE;
	}
	else if ( value > r_max_FERTILE_AGE_RANGE ) {
		value = r_max_FERTILE_AGE_RANGE;
	}
	if ( poParent->fertile_age != value )
	{


		FERTILE_AGE_RANGE	OldValue;
		OldValue = poParent->fertile_age;
		poParent->fertile_age = value;

		// actorset dimension or order state
		if (poParent->actorset_filter_asAllFertilePersonsForFertilityAlignmentByAge)
		{
			asAllFertilePersonsForFertilityAlignmentByAge[RANGE_POS(FERTILE_AGE_RANGE,OldValue)][poParent->is_resident]->Remove(poParent->poasAllFertilePersonsForFertilityAlignmentByAge);
			poParent->poasAllFertilePersonsForFertilityAlignmentByAge = asAllFertilePersonsForFertilityAlignmentByAge[RANGE_POS(FERTILE_AGE_RANGE, poParent->fertile_age)][poParent->is_resident]->NewNode(poParent);
		}
	}
	return poParent->fertile_age = value;
}

GEO_NAT PersonStateFunctions::Set_geo_nat( GEO_NAT value )
{

	CHKLMT(value, geo_nat, 0, rc_size_GEO_NAT - 1);
	if ( value < 0 ) {
		value = (GEO_NAT) 0;
	}
	else if ( value >= rc_size_GEO_NAT ) {
		value = (GEO_NAT) ( rc_size_GEO_NAT - 1 );
	}
	if ( poParent->geo_nat != value )
	{


		GEO_NAT	OldValue;
		OldValue = poParent->geo_nat;
		poParent->geo_nat = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[2] >= 0 && ( poParent->tableinfo.tab2dim_filter || poParent->table_filter_PopPyramidByEduc ) ) {
				if ( poParent->tableinfo.tab2events == -1 || poParent->tableinfo.tab2events < poParent->__events ) {
					if ( poParent->tableinfo.tab2events != -1 ) {
						poParent->tableinfo.tab2open0 = (real) poParent->CloseUpdate( 2, 0, poParent->tableinfo.tab2open0, poParent->tableinfo.tab2close0, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab2open1 = (real) poParent->CloseUpdate( 2, 1, poParent->tableinfo.tab2open1, poParent->tableinfo.tab2close1, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab2open2 = (real) poParent->CloseUpdate( 2, 2, poParent->tableinfo.tab2open2, poParent->tableinfo.tab2close2, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab2events = poParent->__events;
					poParent->tableinfo.tab2close0 = (real) poDerivedStates->Get_duration_educ_one_level_EOL_HIGH();
					poParent->tableinfo.tab2close1 = (real) poDerivedStates->Get_duration_educ_one_level_EOL_LOW();
					poParent->tableinfo.tab2close2 = (real) poDerivedStates->Get_duration_educ_one_level_EOL_MEDIUM();
					poParent->tableinfo.tab2cell = (index) poParent->TableCellIndex( 2, 4,
						poParent->sel_years,
						poParent->sex,
						OldValue,
						poDerivedStates->split_integer_age_TAB_AGEC5 );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab2total_cell = (index) poParent->tableinfo.tab2cell;
					}
					else {
						poParent->tableinfo.tab2total_cell = (index) poParent->TableTotalCellIndex( 2, 4,
						poParent->sel_years,
						poParent->sex,
						OldValue,
						poDerivedStates->split_integer_age_TAB_AGEC5 );
					}
					poParent->tableinfo.tab2dim_filter = poParent->table_filter_PopPyramidByEduc;
				}
			}
			if ( gplTableIndexes[4] >= 0 && ( poParent->tableinfo.tab4dim_filter || poParent->table_filter_TabPrimarySchoolPlanning ) ) {
				if ( poParent->tableinfo.tab4events == -1 || poParent->tableinfo.tab4events < poParent->__events ) {
					if ( poParent->tableinfo.tab4events != -1 ) {
						poParent->tableinfo.tab4open0 = (real) poParent->CloseUpdate( 4, 0, poParent->tableinfo.tab4open0, poParent->tableinfo.tab4close0, poParent->tableinfo.tab4dim_filter, poParent->tableinfo.tab4cell, poParent->tableinfo.tab4total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab4open1 = (real) poParent->CloseUpdate( 4, 1, poParent->tableinfo.tab4open1, poParent->tableinfo.tab4close1, poParent->tableinfo.tab4dim_filter, poParent->tableinfo.tab4cell, poParent->tableinfo.tab4total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						Set_table_unit_TabPrimarySchoolPlanning(poParent->table_unit_TabPrimarySchoolPlanning+1);
						poParent->tableinfo.tab4close2 = (real) poParent->table_unit_TabPrimarySchoolPlanning;
						poParent->tableinfo.tab4open2 = (real) poParent->CloseUpdate( 4, 2, poParent->tableinfo.tab4open2, poParent->tableinfo.tab4close2, poParent->tableinfo.tab4dim_filter, poParent->tableinfo.tab4cell, poParent->tableinfo.tab4total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab4events = poParent->__events;
					poParent->tableinfo.tab4close0 = (real) poParent->educ_one_required_rooms_value_out;
					poParent->tableinfo.tab4close1 = (real) poParent->educ_one_required_teachers_value_out;
					poParent->tableinfo.tab4close2 = (real) poParent->table_unit_TabPrimarySchoolPlanning;
					poParent->tableinfo.tab4cell = (index) poParent->TableCellIndex( 4, 3,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						poParent->educ_one_grade - r_min_EDUC_ONE_GRADE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab4total_cell = (index) poParent->tableinfo.tab4cell;
					}
					else {
						poParent->tableinfo.tab4total_cell = (index) poParent->TableTotalCellIndex( 4, 3,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						poParent->educ_one_grade - r_min_EDUC_ONE_GRADE );
					}
					poParent->tableinfo.tab4dim_filter = poParent->table_filter_TabPrimarySchoolPlanning;
				}
			}
			if ( gplTableIndexes[5] >= 0 && ( poParent->tableinfo.tab5dim_filter || poParent->table_filter_TabPopProvAgeEducSex ) ) {
				if ( poParent->tableinfo.tab5events == -1 || poParent->tableinfo.tab5events < poParent->__events ) {
					if ( poParent->tableinfo.tab5events != -1 ) {
						poParent->tableinfo.tab5open0 = (real) poParent->CloseUpdate( 5, 0, poParent->tableinfo.tab5open0, poParent->tableinfo.tab5close0, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open1 = (real) poParent->CloseUpdate( 5, 1, poParent->tableinfo.tab5open1, poParent->tableinfo.tab5close1, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open2 = (real) poParent->CloseUpdate( 5, 2, poParent->tableinfo.tab5open2, poParent->tableinfo.tab5close2, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open3 = (real) poParent->CloseUpdate( 5, 3, poParent->tableinfo.tab5open3, poParent->tableinfo.tab5close3, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open4 = (real) poParent->CloseUpdate( 5, 4, poParent->tableinfo.tab5open4, poParent->tableinfo.tab5close4, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab5events = poParent->__events;
					poParent->tableinfo.tab5close0 = (real) poDerivedStates->Get_duration_is_aded0_TRUE();
					poParent->tableinfo.tab5close1 = (real) poDerivedStates->Get_duration_is_aded1_TRUE();
					poParent->tableinfo.tab5close2 = (real) poDerivedStates->Get_duration_is_aded2_TRUE();
					poParent->tableinfo.tab5close3 = (real) poDerivedStates->Get_duration_is_child_TRUE();
					poParent->tableinfo.tab5close4 = (real) poDerivedStates->Get_duration_is_old_TRUE();
					poParent->tableinfo.tab5cell = (index) poParent->TableCellIndex( 5, 3,
						OldValue,
						poParent->sex,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab5total_cell = (index) poParent->tableinfo.tab5cell;
					}
					else {
						poParent->tableinfo.tab5total_cell = (index) poParent->TableTotalCellIndex( 5, 3,
						OldValue,
						poParent->sex,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab5dim_filter = poParent->table_filter_TabPopProvAgeEducSex;
				}
			}
			if ( gplTableIndexes[7] >= 0 && ( poParent->tableinfo.tab7dim_filter || poParent->table_filter_TabEduc15ByDistrict ) ) {
				if ( poParent->tableinfo.tab7events == -1 || poParent->tableinfo.tab7events < poParent->__events ) {
					if ( poParent->tableinfo.tab7events != -1 ) {
						poParent->tableinfo.tab7open0 = (real) poParent->CloseUpdate( 7, 0, poParent->tableinfo.tab7open0, poParent->tableinfo.tab7close0, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab7open1 = (real) poParent->CloseUpdate( 7, 1, poParent->tableinfo.tab7open1, poParent->tableinfo.tab7close1, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab7open2 = (real) poParent->CloseUpdate( 7, 2, poParent->tableinfo.tab7open2, poParent->tableinfo.tab7close2, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab7open3 = (real) poParent->CloseUpdate( 7, 3, poParent->tableinfo.tab7open3, poParent->tableinfo.tab7close3, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab7events = poParent->__events;
					poParent->tableinfo.tab7close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab7close1 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH();
					poParent->tableinfo.tab7close2 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_LOW();
					poParent->tableinfo.tab7close3 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM();
					poParent->tableinfo.tab7cell = (index) poParent->TableCellIndex( 7, 3,
						OldValue,
						poParent->sex,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab7total_cell = (index) poParent->tableinfo.tab7cell;
					}
					else {
						poParent->tableinfo.tab7total_cell = (index) poParent->TableTotalCellIndex( 7, 3,
						OldValue,
						poParent->sex,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab7dim_filter = poParent->table_filter_TabEduc15ByDistrict;
				}
			}
			if ( gplTableIndexes[9] >= 0 && ( poParent->tableinfo.tab9dim_filter || poParent->table_filter_TabPrimSchoolEntries ) ) {
				if ( poParent->tableinfo.tab9events == -1 || poParent->tableinfo.tab9events < poParent->__events ) {
					if ( poParent->tableinfo.tab9events != -1 ) {
						poParent->tableinfo.tab9open0 = (real) poParent->CloseUpdate( 9, 0, poParent->tableinfo.tab9open0, poParent->tableinfo.tab9close0, poParent->tableinfo.tab9dim_filter, poParent->tableinfo.tab9cell, poParent->tableinfo.tab9total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab9events = poParent->__events;
					poParent->tableinfo.tab9close0 = (real) poDerivedStates->entrances_educ_one_grade_attended_1;
					poParent->tableinfo.tab9cell = (index) poParent->TableCellIndex( 9, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						OldValue );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab9total_cell = (index) poParent->tableinfo.tab9cell;
					}
					else {
						poParent->tableinfo.tab9total_cell = (index) poParent->TableTotalCellIndex( 9, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						OldValue );
					}
					poParent->tableinfo.tab9dim_filter = poParent->table_filter_TabPrimSchoolEntries;
				}
			}
			if ( gplTableIndexes[10] >= 0 && ( poParent->tableinfo.tab10dim_filter || poParent->table_filter_TabPrimSchoolGraduations ) ) {
				if ( poParent->tableinfo.tab10events == -1 || poParent->tableinfo.tab10events < poParent->__events ) {
					if ( poParent->tableinfo.tab10events != -1 ) {
						poParent->tableinfo.tab10open0 = (real) poParent->CloseUpdate( 10, 0, poParent->tableinfo.tab10open0, poParent->tableinfo.tab10close0, poParent->tableinfo.tab10dim_filter, poParent->tableinfo.tab10cell, poParent->tableinfo.tab10total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab10events = poParent->__events;
					poParent->tableinfo.tab10close0 = (real) poDerivedStates->entrances_has_passed_primary_TRUE;
					poParent->tableinfo.tab10cell = (index) poParent->TableCellIndex( 10, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						OldValue );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab10total_cell = (index) poParent->tableinfo.tab10cell;
					}
					else {
						poParent->tableinfo.tab10total_cell = (index) poParent->TableTotalCellIndex( 10, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						OldValue );
					}
					poParent->tableinfo.tab10dim_filter = poParent->table_filter_TabPrimSchoolGraduations;
				}
			}
			if ( gplTableIndexes[11] >= 0 && ( poParent->tableinfo.tab11dim_filter || poParent->table_filter_TabPrimSchoolOutOfSchool9to11 ) ) {
				if ( poParent->tableinfo.tab11events == -1 || poParent->tableinfo.tab11events < poParent->__events ) {
					if ( poParent->tableinfo.tab11events != -1 ) {
						poParent->tableinfo.tab11open0 = (real) poParent->CloseUpdate( 11, 0, poParent->tableinfo.tab11open0, poParent->tableinfo.tab11close0, poParent->tableinfo.tab11dim_filter, poParent->tableinfo.tab11cell, poParent->tableinfo.tab11total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab11open1 = (real) poParent->CloseUpdate( 11, 1, poParent->tableinfo.tab11open1, poParent->tableinfo.tab11close1, poParent->tableinfo.tab11dim_filter, poParent->tableinfo.tab11cell, poParent->tableinfo.tab11total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab11events = poParent->__events;
					poParent->tableinfo.tab11close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab11close1 = (real) poDerivedStates->Get_duration_tab_out_of_primary_TRUE();
					poParent->tableinfo.tab11cell = (index) poParent->TableCellIndex( 11, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						OldValue );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab11total_cell = (index) poParent->tableinfo.tab11cell;
					}
					else {
						poParent->tableinfo.tab11total_cell = (index) poParent->TableTotalCellIndex( 11, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						OldValue );
					}
					poParent->tableinfo.tab11dim_filter = poParent->table_filter_TabPrimSchoolOutOfSchool9to11;
				}
			}
		}
#line 43 "../code/EducationPrimaryPlanning.mpp"
		Set_educ_one_required_rooms( ( poParent->is_educ_one_student&& poParent->is_resident )? 1.0/ gprParam->Educ1Infrastructure[CHKBND(SI_ROOMS, 2, "Educ1Infrastructure", 0)][CHKBND(poParent->geo_nat, 12, "Educ1Infrastructure", 1)][CHKBND(RANGE_POS(SIM_YEAR_RANGE,poParent->calendar_year), 151, "Educ1Infrastructure", 2)]: 0.0 );
#line 9432 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 39 "../code/EducationPrimaryPlanning.mpp"
		Set_educ_one_required_teachers( ( poParent->is_educ_one_student&& poParent->is_resident )? 1.0/ gprParam->Educ1Infrastructure[CHKBND(SI_TEACHERS, 2, "Educ1Infrastructure", 0)][CHKBND(poParent->geo_nat, 12, "Educ1Infrastructure", 1)][CHKBND(RANGE_POS(SIM_YEAR_RANGE,poParent->calendar_year), 151, "Educ1Infrastructure", 2)]: 0.0 );
#line 9435 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->geo_nat = value;
}

logical PersonStateFunctions::Set_has_passed_primary( logical value )
{

	value = value != 0;
	if ( poParent->has_passed_primary != value )
	{
		if ( value == TRUE ) {
			poDerivedStates->Set_entrances_has_passed_primary_TRUE( poDerivedStates->entrances_has_passed_primary_TRUE + 1 );
		}


		logical	OldValue;
		OldValue = poParent->has_passed_primary;
		poParent->has_passed_primary = value;

#line 28 "../code/TablesEducation.mpp"
		Set_tab_out_of_primary( ( poParent->educ_one_status!= EOS_WAIT&& poParent->educ_one_status!= EOS_ATTEND&& !( poParent->has_passed_primary ) ) );
#line 9457 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->has_passed_primary = value;
}

logical PersonStateFunctions::Set_has_spouse( logical value )
{

	value = value != 0;
	if ( poParent->has_spouse != value )
	{
#line 14 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMale( poParent->is_alive&& poParent->sex== MALE&& !( value )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& poParent->is_mortal );
#line 9470 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 18 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleByType( poParent->is_alive&& poParent->sex== MALE&& !( value )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& poParent->is_mortal );
#line 9473 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 32 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleByTypeImmiScratch( poParent->is_alive&& poParent->creation_type== CT_SCRATCH&& poParent->sex== MALE&& !( value )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->is_mortal ) );
#line 9476 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleImmiScratch( poParent->is_alive&& poParent->creation_type== CT_SCRATCH&& poParent->sex== MALE&& !( value )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->is_mortal ) );
#line 9479 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 22 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAge( poParent->is_alive&& poParent->sex== FEMALE&& value&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& poParent->is_mortal );
#line 9482 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 36 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch( poParent->is_alive&& poParent->creation_type== CT_SCRATCH&& poParent->sex== FEMALE&& value&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->is_mortal ) );
#line 9485 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		if ( poParent->has_spouse == TRUE ) {
			poDerivedStates->Set_duration_has_spouse_TRUE( poDerivedStates->Get_duration_has_spouse_TRUE() );
		}
		else if ( value == TRUE ) {
			poDerivedStates->duration_has_spouse_TRUE_time = poParent->__time;
		}


		logical	OldValue;
		OldValue = poParent->has_spouse;
		poParent->has_spouse = value;

#line 75 "../code/FamilyPartnerMatching.mpp"
		Set_male_spouse_age( ( poParent->has_spouse&& poParent->sex== FEMALE )? COERCE(MALE_SPOUSE_AGE,poParent->lSpouse_integer_age): MAX(MALE_SPOUSE_AGE) );
#line 9500 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 87 "../code/FamilyPartnerMatching.mpp"
		Set_spouse_group( ( poParent->has_spouse&& poParent->lSpouse_educ_one_fate== EOL_HIGH )? SMG_02: ( poParent->has_spouse&& poParent->lSpouse_educ_one_fate== EOL_MEDIUM )? SMG_01: SMG_00 );
#line 9503 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		// actorset dimension or order state
		if (poParent->actorset_filter_asWomenNoChildren)
		{
			asWomenNoChildren[poParent->union1_group][RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][OldValue]->Remove(poParent->poasWomenNoChildren);
			poParent->poasWomenNoChildren = asWomenNoChildren[poParent->union1_group][RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][poParent->has_spouse]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asWomenWithChildren)
		{
			asWomenWithChildren[poParent->union1_group][poParent->child_agegr][poParent->moth_agegr][OldValue]->Remove(poParent->poasWomenWithChildren);
			poParent->poasWomenWithChildren = asWomenWithChildren[poParent->union1_group][poParent->child_agegr][poParent->moth_agegr][poParent->has_spouse]->NewNode(poParent);
		}
	}
	return poParent->has_spouse = value;
}

logical PersonStateFunctions::Set_in_hci_sample( logical value )
{

	value = value != 0;
	if ( poParent->in_hci_sample != value )
	{
#line 17 "../code/TablesHumanCapital.mpp"
		Set_table_filter_tabHCI( poParent->in_projected_time&& value&& poParent->year_of_birth>= MIN(SIM_YEAR_RANGE) );
#line 9527 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 45 "../code/TablesHumanCapital.mpp"
		Set_table_filter_TabHCIDistrict( poParent->in_projected_time&& value&& poParent->year_of_birth>= MIN(SIM_YEAR_RANGE) );
#line 9530 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		logical	OldValue;
		OldValue = poParent->in_hci_sample;
		poParent->in_hci_sample = value;

	}
	return poParent->in_hci_sample = value;
}

logical PersonStateFunctions::Set_in_projected_time( logical value )
{

	value = value != 0;
	if ( poParent->in_projected_time != value )
	{
#line 46 "../code/TablesVaccination.mpp"
		Set_table_filter_TabChildVaccination( value&& poParent->creation_type== CT_BIRTH&& poParent->is_resident&& poParent->integer_age== 0 );
#line 9549 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 165 "../code/TablesEducation.mpp"
		Set_table_filter_TabEduc15ByDistrict( poParent->integer_age== 15&& value&& poParent->is_resident );
#line 9552 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 178 "../code/TablesEducation.mpp"
		Set_table_filter_TabEduc15ByDistrictBirth( poParent->integer_age== 15&& value&& poParent->is_resident );
#line 9555 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 115 "../code/TablesEducation.mpp"
		Set_table_filter_TabEducFateByGroup( value&& poDerivedStates->trigger_changes_calendar_year&& poParent->integer_age== 0&& poParent->creation_type== CT_BIRTH );
#line 9558 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 17 "../code/TablesHumanCapital.mpp"
		Set_table_filter_tabHCI( value&& poParent->in_hci_sample&& poParent->year_of_birth>= MIN(SIM_YEAR_RANGE) );
#line 9561 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 45 "../code/TablesHumanCapital.mpp"
		Set_table_filter_TabHCIDistrict( value&& poParent->in_hci_sample&& poParent->year_of_birth>= MIN(SIM_YEAR_RANGE) );
#line 9564 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 31 "../code/TablesVaccination.mpp"
		Set_table_filter_TabImmunizationChildren( value&& poParent->integer_age== 0&& poParent->creation_type== CT_BIRTH&& poParent->is_resident );
#line 9567 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 139 "../code/TablesEducation.mpp"
		Set_table_filter_TabPopProvAgeEducSex( value&& poParent->is_resident );
#line 9570 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 19 "../code/TablesVaccination.mpp"
		Set_table_filter_TabPrenatCare( value&& poParent->integer_age== 0&& poParent->creation_type== CT_BIRTH&& poParent->is_resident );
#line 9573 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 126 "../code/TablesEducation.mpp"
		Set_table_filter_TabPrimarySchoolPlanning( value&& poDerivedStates->trigger_changes_calendar_year&& poParent->is_educ_one_student&& poParent->is_resident );
#line 9576 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 191 "../code/TablesEducation.mpp"
		Set_table_filter_TabPrimSchoolEntries( value&& poParent->is_resident );
#line 9579 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 201 "../code/TablesEducation.mpp"
		Set_table_filter_TabPrimSchoolGraduations( value&& poParent->is_resident );
#line 9582 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 211 "../code/TablesEducation.mpp"
		Set_table_filter_TabPrimSchoolOutOfSchool9to11( value&& poParent->is_resident&& poParent->integer_age>= 9&& poParent->integer_age<= 11 );
#line 9585 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 231 "../code/TablesEducation.mpp"
		Set_table_filter_TabSchool2AttainmentsTab( value&& poParent->is_resident&& poParent->integer_age>= 20&& poParent->integer_age< 25 );
#line 9588 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 221 "../code/TablesEducation.mpp"
		Set_table_filter_TabSchool2TrackTab( poDerivedStates->trigger_changes_calendar_year&& value&& poParent->is_resident&& poParent->educ_two_status== ETS_ATTEND );
#line 9591 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 27 "../code/TablesStunting.mpp"
		Set_table_filter_tabStuntingSexRegMotherYob( value&& poParent->is_resident&& poParent->creation_type== CT_BIRTH&& poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE );
#line 9594 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 17 "../code/TablesStunting.mpp"
		Set_table_filter_tabStuntingSexRegYob( value&& poParent->is_resident&& poParent->creation_type== CT_BIRTH&& poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE );
#line 9597 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		logical	OldValue;
		OldValue = poParent->in_projected_time;
		poParent->in_projected_time = value;

	}
	return poParent->in_projected_time = value;
}

logical PersonStateFunctions::Set_in_sel_years( logical value )
{

	value = value != 0;
	if ( poParent->in_sel_years != value )
	{
#line 101 "../code/TablesEducation.mpp"
		Set_table_filter_PopPyramidByEduc( value&& poParent->is_resident );
#line 9616 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		logical	OldValue;
		OldValue = poParent->in_sel_years;
		poParent->in_sel_years = value;

	}
	return poParent->in_sel_years = value;
}

logical PersonStateFunctions::Set_in_selected_year( logical value )
{

	value = value != 0;
	if ( poParent->in_selected_year != value )
	{
#line 66 "../code/TablesPopulation.mpp"
		Set_table_filter_tabMigrationOriginDestination( WITHIN(SIM_YEAR_RANGE,poParent->calendar_year)&& poDerivedStates->trigger_changes_geo&& value );
#line 9635 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		logical	OldValue;
		OldValue = poParent->in_selected_year;
		poParent->in_selected_year = value;

	}
	return poParent->in_selected_year = value;
}

logical PersonStateFunctions::Set_is_aded0( logical value )
{

	value = value != 0;
	if ( poParent->is_aded0 != value )
	{
		if ( poParent->is_aded0 == TRUE ) {
			poDerivedStates->Set_duration_is_aded0_TRUE( poDerivedStates->Get_duration_is_aded0_TRUE() );
		}
		else if ( value == TRUE ) {
			poDerivedStates->duration_is_aded0_TRUE_time = poParent->__time;
		}


		logical	OldValue;
		OldValue = poParent->is_aded0;
		poParent->is_aded0 = value;

	}
	return poParent->is_aded0 = value;
}

logical PersonStateFunctions::Set_is_aded1( logical value )
{

	value = value != 0;
	if ( poParent->is_aded1 != value )
	{
		if ( poParent->is_aded1 == TRUE ) {
			poDerivedStates->Set_duration_is_aded1_TRUE( poDerivedStates->Get_duration_is_aded1_TRUE() );
		}
		else if ( value == TRUE ) {
			poDerivedStates->duration_is_aded1_TRUE_time = poParent->__time;
		}


		logical	OldValue;
		OldValue = poParent->is_aded1;
		poParent->is_aded1 = value;

	}
	return poParent->is_aded1 = value;
}

logical PersonStateFunctions::Set_is_aded2( logical value )
{

	value = value != 0;
	if ( poParent->is_aded2 != value )
	{
		if ( poParent->is_aded2 == TRUE ) {
			poDerivedStates->Set_duration_is_aded2_TRUE( poDerivedStates->Get_duration_is_aded2_TRUE() );
		}
		else if ( value == TRUE ) {
			poDerivedStates->duration_is_aded2_TRUE_time = poParent->__time;
		}


		logical	OldValue;
		OldValue = poParent->is_aded2;
		poParent->is_aded2 = value;

	}
	return poParent->is_aded2 = value;
}

logical PersonStateFunctions::Set_is_alive_value_in( logical value )
{

	value = value != 0;
	if ( poParent->is_alive_value_in != value )
	{


		logical	OldValue;
		OldValue = poParent->is_alive_value_in;
		poParent->is_alive_value_in = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[19] >= 0 && poParent->tableinfo.tab19events == poParent->__events ) {
				poParent->tableinfo.tab19close0 = (real) poParent->is_alive_value_in;
			}
		}
	}
	return poParent->is_alive_value_in = value;
}

logical PersonStateFunctions::Set_is_child( logical value )
{

	value = value != 0;
	if ( poParent->is_child != value )
	{
		if ( poParent->is_child == TRUE ) {
			poDerivedStates->Set_duration_is_child_TRUE( poDerivedStates->Get_duration_is_child_TRUE() );
		}
		else if ( value == TRUE ) {
			poDerivedStates->duration_is_child_TRUE_time = poParent->__time;
		}


		logical	OldValue;
		OldValue = poParent->is_child;
		poParent->is_child = value;

#line 22 "../code/TablesEducation.mpp"
		Set_is_aded0( ( !( poParent->is_child )&& !( poParent->is_old )&& poParent->educ_one_fate== EOL_LOW ) );
#line 9753 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 23 "../code/TablesEducation.mpp"
		Set_is_aded1( ( !( poParent->is_child )&& !( poParent->is_old )&& poParent->educ_one_fate== EOL_MEDIUM ) );
#line 9756 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 24 "../code/TablesEducation.mpp"
		Set_is_aded2( ( !( poParent->is_child )&& !( poParent->is_old )&& poParent->educ_one_fate== EOL_HIGH ) );
#line 9759 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->is_child = value;
}

logical PersonStateFunctions::Set_is_educ_one_student( logical value )
{

	value = value != 0;
	if ( poParent->is_educ_one_student != value )
	{
#line 126 "../code/TablesEducation.mpp"
		Set_table_filter_TabPrimarySchoolPlanning( poParent->in_projected_time&& poDerivedStates->trigger_changes_calendar_year&& value&& poParent->is_resident );
#line 9772 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		logical	OldValue;
		OldValue = poParent->is_educ_one_student;
		poParent->is_educ_one_student = value;

#line 43 "../code/EducationPrimaryPlanning.mpp"
		Set_educ_one_required_rooms( ( poParent->is_educ_one_student&& poParent->is_resident )? 1.0/ gprParam->Educ1Infrastructure[CHKBND(SI_ROOMS, 2, "Educ1Infrastructure", 0)][CHKBND(poParent->geo_nat, 12, "Educ1Infrastructure", 1)][CHKBND(RANGE_POS(SIM_YEAR_RANGE,poParent->calendar_year), 151, "Educ1Infrastructure", 2)]: 0.0 );
#line 9781 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 39 "../code/EducationPrimaryPlanning.mpp"
		Set_educ_one_required_teachers( ( poParent->is_educ_one_student&& poParent->is_resident )? 1.0/ gprParam->Educ1Infrastructure[CHKBND(SI_TEACHERS, 2, "Educ1Infrastructure", 0)][CHKBND(poParent->geo_nat, 12, "Educ1Infrastructure", 1)][CHKBND(RANGE_POS(SIM_YEAR_RANGE,poParent->calendar_year), 151, "Educ1Infrastructure", 2)]: 0.0 );
#line 9784 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->is_educ_one_student = value;
}

logical PersonStateFunctions::Set_is_fertile( logical value )
{

	value = value != 0;
	if ( poParent->is_fertile != value )
	{
#line 31 "../code/FertilityDetailed.mpp"
		Set_actorset_filter_asAllFertilePersonsForFertilityAlignment( gprParam->FertilityModel== FEM_ALIGNED_TOTAL&& value );
#line 9797 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 27 "../code/FertilityDetailed.mpp"
		Set_actorset_filter_asAllFertilePersonsForFertilityAlignmentByAge( gprParam->FertilityModel== FEM_ALIGNED_AGE&& value );
#line 9800 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		logical	OldValue;
		OldValue = poParent->is_fertile;
		poParent->is_fertile = value;

		if ( poParent->poFertilityBaseEvent->cRecalc == 0 ) {
			poParent->poFertilityBaseEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFertilityBaseEvent );
		};
		if ( poParent->poFertilityDetailedEvent->cRecalc == 0 ) {
			poParent->poFertilityDetailedEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFertilityDetailedEvent );
		};
	}
	return poParent->is_fertile = value;
}

logical PersonStateFunctions::Set_is_mortal( logical value )
{

	value = value != 0;
	if ( poParent->is_mortal != value )
	{
#line 14 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMale( poParent->is_alive&& poParent->sex== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& value );
#line 9827 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 18 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleByType( poParent->is_alive&& poParent->sex== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& value );
#line 9830 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 32 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleByTypeImmiScratch( poParent->is_alive&& poParent->creation_type== CT_SCRATCH&& poParent->sex== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& !( value ) );
#line 9833 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asAvailableMaleImmiScratch( poParent->is_alive&& poParent->creation_type== CT_SCRATCH&& poParent->sex== MALE&& !( poParent->has_spouse )&& WITHIN(MALE_SPOUSE_AGE,poParent->integer_age)&& !( value ) );
#line 9836 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 22 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAge( poParent->is_alive&& poParent->sex== FEMALE&& poParent->has_spouse&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& value );
#line 9839 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 36 "../code/FamilyPartnerMatching.mpp"
		Set_actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch( poParent->is_alive&& poParent->creation_type== CT_SCRATCH&& poParent->sex== FEMALE&& poParent->has_spouse&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& !( value ) );
#line 9842 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 28 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_actorset_filter_asWomenNoChildren( poParent->is_alive&& poParent->sex== FEMALE&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& !( poParent->lives_with_dependent_child )&& value&& !( poParent->is_blocked_from_marriage ) );
#line 9845 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 23 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_actorset_filter_asWomenWithChildren( poParent->is_alive&& poParent->sex== FEMALE&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& poParent->lives_with_dependent_child&& value&& !( poParent->is_blocked_from_marriage ) );
#line 9848 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		logical	OldValue;
		OldValue = poParent->is_mortal;
		poParent->is_mortal = value;

		if ( poParent->lStartValues != NULL ) {
			if ( poParent->lStartValues->bUpdates ) {
				poParent->lStartValues->UpdateTime();
				poParent->lStartValues->poStateFunctions->Set_lPersonStartpop_is_mortal( poParent->is_mortal );
			}
		}
		if ( poParent->poEmigrationEvent->cRecalc == 0 ) {
			poParent->poEmigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poEmigrationEvent );
		};
		if ( poParent->poFirstUnionFormationEvent->cRecalc == 0 ) {
			poParent->poFirstUnionFormationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstUnionFormationEvent );
		};
		if ( poParent->poMigrationEvent->cRecalc == 0 ) {
			poParent->poMigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMigrationEvent );
		};
		if ( poParent->poMortalityBaseEvent->cRecalc == 0 ) {
			poParent->poMortalityBaseEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMortalityBaseEvent );
		};
		if ( poParent->poMortalityDetailedEvent->cRecalc == 0 ) {
			poParent->poMortalityDetailedEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMortalityDetailedEvent );
		};
	}
	return poParent->is_mortal = value;
}

logical PersonStateFunctions::Set_is_old( logical value )
{

	value = value != 0;
	if ( poParent->is_old != value )
	{
		if ( poParent->is_old == TRUE ) {
			poDerivedStates->Set_duration_is_old_TRUE( poDerivedStates->Get_duration_is_old_TRUE() );
		}
		else if ( value == TRUE ) {
			poDerivedStates->duration_is_old_TRUE_time = poParent->__time;
		}


		logical	OldValue;
		OldValue = poParent->is_old;
		poParent->is_old = value;

#line 22 "../code/TablesEducation.mpp"
		Set_is_aded0( ( !( poParent->is_child )&& !( poParent->is_old )&& poParent->educ_one_fate== EOL_LOW ) );
#line 9905 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 23 "../code/TablesEducation.mpp"
		Set_is_aded1( ( !( poParent->is_child )&& !( poParent->is_old )&& poParent->educ_one_fate== EOL_MEDIUM ) );
#line 9908 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 24 "../code/TablesEducation.mpp"
		Set_is_aded2( ( !( poParent->is_child )&& !( poParent->is_old )&& poParent->educ_one_fate== EOL_HIGH ) );
#line 9911 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->is_old = value;
}

logical PersonStateFunctions::Set_is_resident( logical value )
{

	value = value != 0;
	if ( poParent->is_resident != value )
	{
#line 18 "../code/ChildVaccination.mpp"
		Set_actorset_filter_asResidentBabies( poParent->is_alive&& poParent->integer_age== 0&& value&& poParent->creation_type== CT_BIRTH );
#line 9924 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 30 "../code/EducationPrimaryBase.mpp"
		Set_actorset_filter_asResidentsAge0Sex( poParent->is_alive&& poParent->integer_age== 0&& value );
#line 9927 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 27 "../code/EducationPrimaryBase.mpp"
		Set_actorset_filter_asResidentsAge0SexGeo( poParent->is_alive&& value&& poParent->integer_age== 0 );
#line 9930 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 101 "../code/TablesEducation.mpp"
		Set_table_filter_PopPyramidByEduc( poParent->in_sel_years&& value );
#line 9933 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 46 "../code/TablesVaccination.mpp"
		Set_table_filter_TabChildVaccination( poParent->in_projected_time&& poParent->creation_type== CT_BIRTH&& value&& poParent->integer_age== 0 );
#line 9936 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 165 "../code/TablesEducation.mpp"
		Set_table_filter_TabEduc15ByDistrict( poParent->integer_age== 15&& poParent->in_projected_time&& value );
#line 9939 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 178 "../code/TablesEducation.mpp"
		Set_table_filter_TabEduc15ByDistrictBirth( poParent->integer_age== 15&& poParent->in_projected_time&& value );
#line 9942 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 154 "../code/TablesEducation.mpp"
		Set_table_filter_TabEducFateDistrYob( poDerivedStates->trigger_changes_calendar_year&& poParent->integer_age== 0&& value );
#line 9945 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 13 "../code/TablesFamily.mpp"
		Set_table_filter_tabHavingSpouse( WITHIN(SIM_YEAR_RANGE,poParent->calendar_year)&& value );
#line 9948 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 31 "../code/TablesVaccination.mpp"
		Set_table_filter_TabImmunizationChildren( poParent->in_projected_time&& poParent->integer_age== 0&& poParent->creation_type== CT_BIRTH&& value );
#line 9951 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 139 "../code/TablesEducation.mpp"
		Set_table_filter_TabPopProvAgeEducSex( poParent->in_projected_time&& value );
#line 9954 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 19 "../code/TablesVaccination.mpp"
		Set_table_filter_TabPrenatCare( poParent->in_projected_time&& poParent->integer_age== 0&& poParent->creation_type== CT_BIRTH&& value );
#line 9957 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 126 "../code/TablesEducation.mpp"
		Set_table_filter_TabPrimarySchoolPlanning( poParent->in_projected_time&& poDerivedStates->trigger_changes_calendar_year&& poParent->is_educ_one_student&& value );
#line 9960 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 191 "../code/TablesEducation.mpp"
		Set_table_filter_TabPrimSchoolEntries( poParent->in_projected_time&& value );
#line 9963 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 201 "../code/TablesEducation.mpp"
		Set_table_filter_TabPrimSchoolGraduations( poParent->in_projected_time&& value );
#line 9966 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 211 "../code/TablesEducation.mpp"
		Set_table_filter_TabPrimSchoolOutOfSchool9to11( poParent->in_projected_time&& value&& poParent->integer_age>= 9&& poParent->integer_age<= 11 );
#line 9969 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 231 "../code/TablesEducation.mpp"
		Set_table_filter_TabSchool2AttainmentsTab( poParent->in_projected_time&& value&& poParent->integer_age>= 20&& poParent->integer_age< 25 );
#line 9972 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 221 "../code/TablesEducation.mpp"
		Set_table_filter_TabSchool2TrackTab( poDerivedStates->trigger_changes_calendar_year&& poParent->in_projected_time&& value&& poParent->educ_two_status== ETS_ATTEND );
#line 9975 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 27 "../code/TablesStunting.mpp"
		Set_table_filter_tabStuntingSexRegMotherYob( poParent->in_projected_time&& value&& poParent->creation_type== CT_BIRTH&& poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE );
#line 9978 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 17 "../code/TablesStunting.mpp"
		Set_table_filter_tabStuntingSexRegYob( poParent->in_projected_time&& value&& poParent->creation_type== CT_BIRTH&& poDerivedStates->trigger_entrances_is_ready_for_birthtables_TRUE );
#line 9981 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		logical	OldValue;
		OldValue = poParent->is_resident;
		poParent->is_resident = value;

#line 43 "../code/EducationPrimaryPlanning.mpp"
		Set_educ_one_required_rooms( ( poParent->is_educ_one_student&& poParent->is_resident )? 1.0/ gprParam->Educ1Infrastructure[CHKBND(SI_ROOMS, 2, "Educ1Infrastructure", 0)][CHKBND(poParent->geo_nat, 12, "Educ1Infrastructure", 1)][CHKBND(RANGE_POS(SIM_YEAR_RANGE,poParent->calendar_year), 151, "Educ1Infrastructure", 2)]: 0.0 );
#line 9990 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 39 "../code/EducationPrimaryPlanning.mpp"
		Set_educ_one_required_teachers( ( poParent->is_educ_one_student&& poParent->is_resident )? 1.0/ gprParam->Educ1Infrastructure[CHKBND(SI_TEACHERS, 2, "Educ1Infrastructure", 0)][CHKBND(poParent->geo_nat, 12, "Educ1Infrastructure", 1)][CHKBND(RANGE_POS(SIM_YEAR_RANGE,poParent->calendar_year), 151, "Educ1Infrastructure", 2)]: 0.0 );
#line 9993 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
		// actorset dimension or order state
		if (poParent->actorset_filter_asAllFertilePersonsForFertilityAlignment)
		{
			asAllFertilePersonsForFertilityAlignment[OldValue]->Remove(poParent->poasAllFertilePersonsForFertilityAlignment);
			poParent->poasAllFertilePersonsForFertilityAlignment = asAllFertilePersonsForFertilityAlignment[poParent->is_resident]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asAllFertilePersonsForFertilityAlignmentByAge)
		{
			asAllFertilePersonsForFertilityAlignmentByAge[RANGE_POS(FERTILE_AGE_RANGE, poParent->fertile_age)][OldValue]->Remove(poParent->poasAllFertilePersonsForFertilityAlignmentByAge);
			poParent->poasAllFertilePersonsForFertilityAlignmentByAge = asAllFertilePersonsForFertilityAlignmentByAge[RANGE_POS(FERTILE_AGE_RANGE, poParent->fertile_age)][poParent->is_resident]->NewNode(poParent);
		}
		if ( poParent->poBackMigrationEvent->cRecalc == 0 ) {
			poParent->poBackMigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poBackMigrationEvent );
		};
		if ( poParent->poEmigrationEvent->cRecalc == 0 ) {
			poParent->poEmigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poEmigrationEvent );
		};
		if ( poParent->poMigrationEvent->cRecalc == 0 ) {
			poParent->poMigrationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poMigrationEvent );
		};
		if ( poParent->poSetPreschoolYearsEvent->cRecalc == 0 ) {
			poParent->poSetPreschoolYearsEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poSetPreschoolYearsEvent );
		};
	}
	return poParent->is_resident = value;
}

logical PersonStateFunctions::Set_is_stunted_value_in( logical value )
{

	value = value != 0;
	if ( poParent->is_stunted_value_in != value )
	{


		logical	OldValue;
		OldValue = poParent->is_stunted_value_in;
		poParent->is_stunted_value_in = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[20] >= 0 && poParent->tableinfo.tab20events == poParent->__events ) {
				poParent->tableinfo.tab20close0 = (real) poParent->is_stunted_value_in;
			}
			if ( gplTableIndexes[21] >= 0 && poParent->tableinfo.tab21events == poParent->__events ) {
				poParent->tableinfo.tab21close0 = (real) poParent->is_stunted_value_in;
			}
		}
	}
	return poParent->is_stunted_value_in = value;
}

logical PersonStateFunctions::Set_lives_with_dependent_child( logical value )
{

	value = value != 0;
	if ( poParent->lives_with_dependent_child != value )
	{
#line 28 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_actorset_filter_asWomenNoChildren( poParent->is_alive&& poParent->sex== FEMALE&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& !( value )&& poParent->is_mortal&& !( poParent->is_blocked_from_marriage ) );
#line 10057 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 23 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_actorset_filter_asWomenWithChildren( poParent->is_alive&& poParent->sex== FEMALE&& WITHIN(FEMALE_SPOUSE_AGE,poParent->integer_age)&& value&& poParent->is_mortal&& !( poParent->is_blocked_from_marriage ) );
#line 10060 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		logical	OldValue;
		OldValue = poParent->lives_with_dependent_child;
		poParent->lives_with_dependent_child = value;

#line 96 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_age_last_birth( ( poParent->lives_with_dependent_child )? poParent->integer_age- poParent->age_youngest_child: TIME_INFINITE );
#line 10069 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lives_with_dependent_child = value;
}

MALE_SPOUSE_AGE PersonStateFunctions::Set_male_age( MALE_SPOUSE_AGE value )
{

	CHKLMT(value, male_age, r_min_MALE_SPOUSE_AGE, r_max_MALE_SPOUSE_AGE);
	if ( value < r_min_MALE_SPOUSE_AGE ) {
		value = r_min_MALE_SPOUSE_AGE;
	}
	else if ( value > r_max_MALE_SPOUSE_AGE ) {
		value = r_max_MALE_SPOUSE_AGE;
	}
	if ( poParent->male_age != value )
	{


		MALE_SPOUSE_AGE	OldValue;
		OldValue = poParent->male_age;
		poParent->male_age = value;

		// actorset dimension or order state
		if (poParent->actorset_filter_asAvailableMale)
		{
			asAvailableMale[RANGE_POS(MALE_SPOUSE_AGE,OldValue)][poParent->geo]->Remove(poParent->poasAvailableMale);
			poParent->poasAvailableMale = asAvailableMale[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->geo]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asAvailableMaleByType)
		{
			asAvailableMaleByType[RANGE_POS(MALE_SPOUSE_AGE,OldValue)][poParent->spouse_group][poParent->geo]->Remove(poParent->poasAvailableMaleByType);
			poParent->poasAvailableMaleByType = asAvailableMaleByType[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->spouse_group][poParent->geo]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asAvailableMaleByTypeImmiScratch)
		{
			asAvailableMaleByTypeImmiScratch[RANGE_POS(MALE_SPOUSE_AGE,OldValue)][poParent->spouse_group][poParent->geo][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->Remove(poParent->poasAvailableMaleByTypeImmiScratch);
			poParent->poasAvailableMaleByTypeImmiScratch = asAvailableMaleByTypeImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->spouse_group][poParent->geo][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asAvailableMaleImmiScratch)
		{
			asAvailableMaleImmiScratch[RANGE_POS(MALE_SPOUSE_AGE,OldValue)][poParent->geo][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->Remove(poParent->poasAvailableMaleImmiScratch);
			poParent->poasAvailableMaleImmiScratch = asAvailableMaleImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->geo][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->NewNode(poParent);
		}
	}
	return poParent->male_age = value;
}

MALE_SPOUSE_AGE PersonStateFunctions::Set_male_spouse_age( MALE_SPOUSE_AGE value )
{

	CHKLMT(value, male_spouse_age, r_min_MALE_SPOUSE_AGE, r_max_MALE_SPOUSE_AGE);
	if ( value < r_min_MALE_SPOUSE_AGE ) {
		value = r_min_MALE_SPOUSE_AGE;
	}
	else if ( value > r_max_MALE_SPOUSE_AGE ) {
		value = r_max_MALE_SPOUSE_AGE;
	}
	if ( poParent->male_spouse_age != value )
	{


		MALE_SPOUSE_AGE	OldValue;
		OldValue = poParent->male_spouse_age;
		poParent->male_spouse_age = value;

		// actorset dimension or order state
		if (poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAge)
		{
			asFemaleInUnionByAgeAndPartnerAge[RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][RANGE_POS(MALE_SPOUSE_AGE,OldValue)]->Remove(poParent->poasFemaleInUnionByAgeAndPartnerAge);
			poParent->poasFemaleInUnionByAgeAndPartnerAge = asFemaleInUnionByAgeAndPartnerAge[RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][RANGE_POS(MALE_SPOUSE_AGE, poParent->male_spouse_age)]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch)
		{
			asFemaleInUnionByAgeAndPartnerAgeImmiScratch[RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][RANGE_POS(MALE_SPOUSE_AGE,OldValue)][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->Remove(poParent->poasFemaleInUnionByAgeAndPartnerAgeImmiScratch);
			poParent->poasFemaleInUnionByAgeAndPartnerAgeImmiScratch = asFemaleInUnionByAgeAndPartnerAgeImmiScratch[RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][RANGE_POS(MALE_SPOUSE_AGE, poParent->male_spouse_age)][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->NewNode(poParent);
		}
	}
	return poParent->male_spouse_age = value;
}

MOTH_AGEGR PersonStateFunctions::Set_moth_agegr( MOTH_AGEGR value )
{

	CHKLMT(value, moth_agegr, 0, rc_size_MOTH_AGEGR - 1);
	if ( value < 0 ) {
		value = (MOTH_AGEGR) 0;
	}
	else if ( value >= rc_size_MOTH_AGEGR ) {
		value = (MOTH_AGEGR) ( rc_size_MOTH_AGEGR - 1 );
	}
	if ( poParent->moth_agegr != value )
	{


		MOTH_AGEGR	OldValue;
		OldValue = poParent->moth_agegr;
		poParent->moth_agegr = value;

		// actorset dimension or order state
		if (poParent->actorset_filter_asWomenWithChildren)
		{
			asWomenWithChildren[poParent->union1_group][poParent->child_agegr][OldValue][poParent->has_spouse]->Remove(poParent->poasWomenWithChildren);
			poParent->poasWomenWithChildren = asWomenWithChildren[poParent->union1_group][poParent->child_agegr][poParent->moth_agegr][poParent->has_spouse]->NewNode(poParent);
		}
	}
	return poParent->moth_agegr = value;
}

int PersonStateFunctions::Set_moth_agegr_part( int value )
{

	if ( poParent->moth_agegr_part != value )
	{


		int	OldValue;
		OldValue = poParent->moth_agegr_part;
		poParent->moth_agegr_part = value;

#line 106 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_moth_agegr( ( poParent->moth_agegr_part== 0 )? CMA20: ( poParent->moth_agegr_part== 1 )? CMA25: ( poParent->moth_agegr_part== 2 )? CMA30: ( poParent->moth_agegr_part== 3 )? CMA35: ( poParent->moth_agegr_part== 4 )? CMA40: CMA40P );
#line 10191 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->moth_agegr_part = value;
}

REGION_INT PersonStateFunctions::Set_region_birth( REGION_INT value )
{

	CHKLMT(value, region_birth, 0, rc_size_REGION_INT - 1);
	if ( value < 0 ) {
		value = (REGION_INT) 0;
	}
	else if ( value >= rc_size_REGION_INT ) {
		value = (REGION_INT) ( rc_size_REGION_INT - 1 );
	}
	if ( poParent->region_birth != value )
	{
		poDerivedStates->Set_aggregate_region_birth_REGION_NAT( REGION_INT_To_REGION_NAT( value ) );


		REGION_INT	OldValue;
		OldValue = poParent->region_birth;
		poParent->region_birth = value;

	}
	return poParent->region_birth = value;
}

REGION_INT PersonStateFunctions::Set_region_int( REGION_INT value )
{

	CHKLMT(value, region_int, 0, rc_size_REGION_INT - 1);
	if ( value < 0 ) {
		value = (REGION_INT) 0;
	}
	else if ( value >= rc_size_REGION_INT ) {
		value = (REGION_INT) ( rc_size_REGION_INT - 1 );
	}
	if ( poParent->region_int != value )
	{
		poDerivedStates->Set_aggregate_region_int_REGION_NAT( REGION_INT_To_REGION_NAT( value ) );


		REGION_INT	OldValue;
		OldValue = poParent->region_int;
		poParent->region_int = value;

	}
	return poParent->region_int = value;
}

REGION_NAT PersonStateFunctions::Set_region_nat( REGION_NAT value )
{

	CHKLMT(value, region_nat, 0, rc_size_REGION_NAT - 1);
	if ( value < 0 ) {
		value = (REGION_NAT) 0;
	}
	else if ( value >= rc_size_REGION_NAT ) {
		value = (REGION_NAT) ( rc_size_REGION_NAT - 1 );
	}
	if ( poParent->region_nat != value )
	{


		REGION_NAT	OldValue;
		OldValue = poParent->region_nat;
		poParent->region_nat = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[20] >= 0 && ( poParent->tableinfo.tab20dim_filter || poParent->table_filter_tabStuntingSexRegYob ) ) {
				if ( poParent->tableinfo.tab20events == -1 || poParent->tableinfo.tab20events < poParent->__events ) {
					if ( poParent->tableinfo.tab20events != -1 ) {
						poParent->tableinfo.tab20open0 = (real) poParent->CloseUpdate( 20, 0, poParent->tableinfo.tab20open0, poParent->tableinfo.tab20close0, poParent->tableinfo.tab20dim_filter, poParent->tableinfo.tab20cell, poParent->tableinfo.tab20total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						Set_table_unit_tabStuntingSexRegYob(poParent->table_unit_tabStuntingSexRegYob+1);
						poParent->tableinfo.tab20close1 = (real) poParent->table_unit_tabStuntingSexRegYob;
						poParent->tableinfo.tab20open1 = (real) poParent->CloseUpdate( 20, 1, poParent->tableinfo.tab20open1, poParent->tableinfo.tab20close1, poParent->tableinfo.tab20dim_filter, poParent->tableinfo.tab20cell, poParent->tableinfo.tab20total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab20events = poParent->__events;
					poParent->tableinfo.tab20close0 = (real) poParent->is_stunted_value_in;
					poParent->tableinfo.tab20close1 = (real) poParent->table_unit_tabStuntingSexRegYob;
					poParent->tableinfo.tab20cell = (index) poParent->TableCellIndex( 20, 2,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab20total_cell = (index) poParent->tableinfo.tab20cell;
					}
					else {
						poParent->tableinfo.tab20total_cell = (index) poParent->TableTotalCellIndex( 20, 2,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab20dim_filter = poParent->table_filter_tabStuntingSexRegYob;
				}
			}
			if ( gplTableIndexes[21] >= 0 && ( poParent->tableinfo.tab21dim_filter || poParent->table_filter_tabStuntingSexRegMotherYob ) ) {
				if ( poParent->tableinfo.tab21events == -1 || poParent->tableinfo.tab21events < poParent->__events ) {
					if ( poParent->tableinfo.tab21events != -1 ) {
						poParent->tableinfo.tab21open0 = (real) poParent->CloseUpdate( 21, 0, poParent->tableinfo.tab21open0, poParent->tableinfo.tab21close0, poParent->tableinfo.tab21dim_filter, poParent->tableinfo.tab21cell, poParent->tableinfo.tab21total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						Set_table_unit_tabStuntingSexRegMotherYob(poParent->table_unit_tabStuntingSexRegMotherYob+1);
						poParent->tableinfo.tab21close1 = (real) poParent->table_unit_tabStuntingSexRegMotherYob;
						poParent->tableinfo.tab21open1 = (real) poParent->CloseUpdate( 21, 1, poParent->tableinfo.tab21open1, poParent->tableinfo.tab21close1, poParent->tableinfo.tab21dim_filter, poParent->tableinfo.tab21cell, poParent->tableinfo.tab21total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab21events = poParent->__events;
					poParent->tableinfo.tab21close0 = (real) poParent->is_stunted_value_in;
					poParent->tableinfo.tab21close1 = (real) poParent->table_unit_tabStuntingSexRegMotherYob;
					poParent->tableinfo.tab21cell = (index) poParent->TableCellIndex( 21, 4,
						poParent->sex,
						poParent->educ_mother,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab21total_cell = (index) poParent->tableinfo.tab21cell;
					}
					else {
						poParent->tableinfo.tab21total_cell = (index) poParent->TableTotalCellIndex( 21, 4,
						poParent->sex,
						poParent->educ_mother,
						OldValue,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab21dim_filter = poParent->table_filter_tabStuntingSexRegMotherYob;
				}
			}
			if ( gplTableIndexes[22] >= 0 && ( poParent->tableinfo.tab22dim_filter || poParent->table_filter_TabPrenatCare ) ) {
				if ( poParent->tableinfo.tab22events == -1 || poParent->tableinfo.tab22events < poParent->__events ) {
					if ( poParent->tableinfo.tab22events != -1 ) {
						poParent->tableinfo.tab22open0 = (real) poParent->CloseUpdate( 22, 0, poParent->tableinfo.tab22open0, poParent->tableinfo.tab22close0, poParent->tableinfo.tab22dim_filter, poParent->tableinfo.tab22cell, poParent->tableinfo.tab22total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab22open1 = (real) poParent->CloseUpdate( 22, 1, poParent->tableinfo.tab22open1, poParent->tableinfo.tab22close1, poParent->tableinfo.tab22dim_filter, poParent->tableinfo.tab22cell, poParent->tableinfo.tab22total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab22events = poParent->__events;
					poParent->tableinfo.tab22close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab22close1 = (real) poDerivedStates->Get_duration_got_prenat_care_GPC_YES();
					poParent->tableinfo.tab22cell = (index) poParent->TableCellIndex( 22, 3,
						OldValue,
						poParent->educ_mother,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab22total_cell = (index) poParent->tableinfo.tab22cell;
					}
					else {
						poParent->tableinfo.tab22total_cell = (index) poParent->TableTotalCellIndex( 22, 3,
						OldValue,
						poParent->educ_mother,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab22dim_filter = poParent->table_filter_TabPrenatCare;
				}
			}
			if ( gplTableIndexes[23] >= 0 && ( poParent->tableinfo.tab23dim_filter || poParent->table_filter_TabImmunizationChildren ) ) {
				if ( poParent->tableinfo.tab23events == -1 || poParent->tableinfo.tab23events < poParent->__events ) {
					if ( poParent->tableinfo.tab23events != -1 ) {
						poParent->tableinfo.tab23open0 = (real) poParent->CloseUpdate( 23, 0, poParent->tableinfo.tab23open0, poParent->tableinfo.tab23close0, poParent->tableinfo.tab23dim_filter, poParent->tableinfo.tab23cell, poParent->tableinfo.tab23total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab23open1 = (real) poParent->CloseUpdate( 23, 1, poParent->tableinfo.tab23open1, poParent->tableinfo.tab23close1, poParent->tableinfo.tab23dim_filter, poParent->tableinfo.tab23cell, poParent->tableinfo.tab23total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab23events = poParent->__events;
					poParent->tableinfo.tab23close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab23close1 = (real) poDerivedStates->Get_duration_is_immunized_TRUE();
					poParent->tableinfo.tab23cell = (index) poParent->TableCellIndex( 23, 6,
						poParent->sex,
						poParent->ethnicity_short,
						poParent->got_prenat_care,
						OldValue,
						poParent->educ_mother,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab23total_cell = (index) poParent->tableinfo.tab23cell;
					}
					else {
						poParent->tableinfo.tab23total_cell = (index) poParent->TableTotalCellIndex( 23, 6,
						poParent->sex,
						poParent->ethnicity_short,
						poParent->got_prenat_care,
						OldValue,
						poParent->educ_mother,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab23dim_filter = poParent->table_filter_TabImmunizationChildren;
				}
			}
			if ( gplTableIndexes[24] >= 0 && ( poParent->tableinfo.tab24dim_filter || poParent->table_filter_TabChildVaccination ) ) {
				if ( poParent->tableinfo.tab24events == -1 || poParent->tableinfo.tab24events < poParent->__events ) {
					if ( poParent->tableinfo.tab24events != -1 ) {
						poParent->tableinfo.tab24open0 = (real) poParent->CloseUpdate( 24, 0, poParent->tableinfo.tab24open0, poParent->tableinfo.tab24close0, poParent->tableinfo.tab24dim_filter, poParent->tableinfo.tab24cell, poParent->tableinfo.tab24total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab24open1 = (real) poParent->CloseUpdate( 24, 1, poParent->tableinfo.tab24open1, poParent->tableinfo.tab24close1, poParent->tableinfo.tab24dim_filter, poParent->tableinfo.tab24cell, poParent->tableinfo.tab24total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab24events = poParent->__events;
					poParent->tableinfo.tab24close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab24close1 = (real) poDerivedStates->Get_duration_is_immunized_TRUE();
					poParent->tableinfo.tab24cell = (index) poParent->TableCellIndex( 24, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						OldValue );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab24total_cell = (index) poParent->tableinfo.tab24cell;
					}
					else {
						poParent->tableinfo.tab24total_cell = (index) poParent->TableTotalCellIndex( 24, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						OldValue );
					}
					poParent->tableinfo.tab24dim_filter = poParent->table_filter_TabChildVaccination;
				}
			}
		}
	}
	return poParent->region_nat = value;
}

TAB_SEL_YEARS PersonStateFunctions::Set_sel_years( TAB_SEL_YEARS value )
{

	CHKLMT(value, sel_years, 0, rc_size_TAB_SEL_YEARS - 1);
	if ( value < 0 ) {
		value = (TAB_SEL_YEARS) 0;
	}
	else if ( value >= rc_size_TAB_SEL_YEARS ) {
		value = (TAB_SEL_YEARS) ( rc_size_TAB_SEL_YEARS - 1 );
	}
	if ( poParent->sel_years != value )
	{


		TAB_SEL_YEARS	OldValue;
		OldValue = poParent->sel_years;
		poParent->sel_years = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[2] >= 0 && ( poParent->tableinfo.tab2dim_filter || poParent->table_filter_PopPyramidByEduc ) ) {
				if ( poParent->tableinfo.tab2events == -1 || poParent->tableinfo.tab2events < poParent->__events ) {
					if ( poParent->tableinfo.tab2events != -1 ) {
						poParent->tableinfo.tab2open0 = (real) poParent->CloseUpdate( 2, 0, poParent->tableinfo.tab2open0, poParent->tableinfo.tab2close0, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab2open1 = (real) poParent->CloseUpdate( 2, 1, poParent->tableinfo.tab2open1, poParent->tableinfo.tab2close1, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab2open2 = (real) poParent->CloseUpdate( 2, 2, poParent->tableinfo.tab2open2, poParent->tableinfo.tab2close2, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab2events = poParent->__events;
					poParent->tableinfo.tab2close0 = (real) poDerivedStates->Get_duration_educ_one_level_EOL_HIGH();
					poParent->tableinfo.tab2close1 = (real) poDerivedStates->Get_duration_educ_one_level_EOL_LOW();
					poParent->tableinfo.tab2close2 = (real) poDerivedStates->Get_duration_educ_one_level_EOL_MEDIUM();
					poParent->tableinfo.tab2cell = (index) poParent->TableCellIndex( 2, 4,
						OldValue,
						poParent->sex,
						poParent->geo_nat,
						poDerivedStates->split_integer_age_TAB_AGEC5 );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab2total_cell = (index) poParent->tableinfo.tab2cell;
					}
					else {
						poParent->tableinfo.tab2total_cell = (index) poParent->TableTotalCellIndex( 2, 4,
						OldValue,
						poParent->sex,
						poParent->geo_nat,
						poDerivedStates->split_integer_age_TAB_AGEC5 );
					}
					poParent->tableinfo.tab2dim_filter = poParent->table_filter_PopPyramidByEduc;
				}
			}
		}
	}
	return poParent->sel_years = value;
}

SELECTED_YEAR10 PersonStateFunctions::Set_selected_year10( SELECTED_YEAR10 value )
{

	CHKLMT(value, selected_year10, 0, rc_size_SELECTED_YEAR10 - 1);
	if ( value < 0 ) {
		value = (SELECTED_YEAR10) 0;
	}
	else if ( value >= rc_size_SELECTED_YEAR10 ) {
		value = (SELECTED_YEAR10) ( rc_size_SELECTED_YEAR10 - 1 );
	}
	if ( poParent->selected_year10 != value )
	{


		SELECTED_YEAR10	OldValue;
		OldValue = poParent->selected_year10;
		poParent->selected_year10 = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[19] >= 0 && ( poParent->tableinfo.tab19dim_filter || poParent->table_filter_tabMigrationOriginDestination ) ) {
				if ( poParent->tableinfo.tab19events == -1 || poParent->tableinfo.tab19events < poParent->__events ) {
					if ( poParent->tableinfo.tab19events != -1 ) {
						poParent->tableinfo.tab19open0 = (real) poParent->CloseUpdate( 19, 0, poParent->tableinfo.tab19open0, poParent->tableinfo.tab19close0, poParent->tableinfo.tab19dim_filter, poParent->tableinfo.tab19cell, poParent->tableinfo.tab19total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab19events = poParent->__events;
					poParent->tableinfo.tab19close0 = (real) poParent->is_alive_value_in;
					poParent->tableinfo.tab19cell = (index) poParent->TableCellIndex( 19, 3,
						OldValue,
						poParent->geo_prev,
						poParent->geo );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab19total_cell = (index) poParent->tableinfo.tab19cell;
					}
					else {
						poParent->tableinfo.tab19total_cell = (index) poParent->TableTotalCellIndex( 19, 3,
						OldValue,
						poParent->geo_prev,
						poParent->geo );
					}
					poParent->tableinfo.tab19dim_filter = poParent->table_filter_tabMigrationOriginDestination;
				}
			}
		}
#line 34 "../code/TablesPopulation.mpp"
		Set_in_selected_year( ( poParent->selected_year10!= SY40|| poParent->calendar_year== MIN(SIM_YEAR_RANGE)+ 40 ) );
#line 10497 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->selected_year10 = value;
}

SIM_YEAR_RANGE PersonStateFunctions::Set_sim_year( SIM_YEAR_RANGE value )
{

	CHKLMT(value, sim_year, r_min_SIM_YEAR_RANGE, r_max_SIM_YEAR_RANGE);
	if ( value < r_min_SIM_YEAR_RANGE ) {
		value = r_min_SIM_YEAR_RANGE;
	}
	else if ( value > r_max_SIM_YEAR_RANGE ) {
		value = r_max_SIM_YEAR_RANGE;
	}
	if ( poParent->sim_year != value )
	{


		SIM_YEAR_RANGE	OldValue;
		OldValue = poParent->sim_year;
		poParent->sim_year = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[4] >= 0 && ( poParent->tableinfo.tab4dim_filter || poParent->table_filter_TabPrimarySchoolPlanning ) ) {
				if ( poParent->tableinfo.tab4events == -1 || poParent->tableinfo.tab4events < poParent->__events ) {
					if ( poParent->tableinfo.tab4events != -1 ) {
						poParent->tableinfo.tab4open0 = (real) poParent->CloseUpdate( 4, 0, poParent->tableinfo.tab4open0, poParent->tableinfo.tab4close0, poParent->tableinfo.tab4dim_filter, poParent->tableinfo.tab4cell, poParent->tableinfo.tab4total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab4open1 = (real) poParent->CloseUpdate( 4, 1, poParent->tableinfo.tab4open1, poParent->tableinfo.tab4close1, poParent->tableinfo.tab4dim_filter, poParent->tableinfo.tab4cell, poParent->tableinfo.tab4total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						Set_table_unit_TabPrimarySchoolPlanning(poParent->table_unit_TabPrimarySchoolPlanning+1);
						poParent->tableinfo.tab4close2 = (real) poParent->table_unit_TabPrimarySchoolPlanning;
						poParent->tableinfo.tab4open2 = (real) poParent->CloseUpdate( 4, 2, poParent->tableinfo.tab4open2, poParent->tableinfo.tab4close2, poParent->tableinfo.tab4dim_filter, poParent->tableinfo.tab4cell, poParent->tableinfo.tab4total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab4events = poParent->__events;
					poParent->tableinfo.tab4close0 = (real) poParent->educ_one_required_rooms_value_out;
					poParent->tableinfo.tab4close1 = (real) poParent->educ_one_required_teachers_value_out;
					poParent->tableinfo.tab4close2 = (real) poParent->table_unit_TabPrimarySchoolPlanning;
					poParent->tableinfo.tab4cell = (index) poParent->TableCellIndex( 4, 3,
						poParent->geo_nat,
						OldValue - r_min_SIM_YEAR_RANGE,
						poParent->educ_one_grade - r_min_EDUC_ONE_GRADE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab4total_cell = (index) poParent->tableinfo.tab4cell;
					}
					else {
						poParent->tableinfo.tab4total_cell = (index) poParent->TableTotalCellIndex( 4, 3,
						poParent->geo_nat,
						OldValue - r_min_SIM_YEAR_RANGE,
						poParent->educ_one_grade - r_min_EDUC_ONE_GRADE );
					}
					poParent->tableinfo.tab4dim_filter = poParent->table_filter_TabPrimarySchoolPlanning;
				}
			}
			if ( gplTableIndexes[5] >= 0 && ( poParent->tableinfo.tab5dim_filter || poParent->table_filter_TabPopProvAgeEducSex ) ) {
				if ( poParent->tableinfo.tab5events == -1 || poParent->tableinfo.tab5events < poParent->__events ) {
					if ( poParent->tableinfo.tab5events != -1 ) {
						poParent->tableinfo.tab5open0 = (real) poParent->CloseUpdate( 5, 0, poParent->tableinfo.tab5open0, poParent->tableinfo.tab5close0, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open1 = (real) poParent->CloseUpdate( 5, 1, poParent->tableinfo.tab5open1, poParent->tableinfo.tab5close1, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open2 = (real) poParent->CloseUpdate( 5, 2, poParent->tableinfo.tab5open2, poParent->tableinfo.tab5close2, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open3 = (real) poParent->CloseUpdate( 5, 3, poParent->tableinfo.tab5open3, poParent->tableinfo.tab5close3, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open4 = (real) poParent->CloseUpdate( 5, 4, poParent->tableinfo.tab5open4, poParent->tableinfo.tab5close4, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab5events = poParent->__events;
					poParent->tableinfo.tab5close0 = (real) poDerivedStates->Get_duration_is_aded0_TRUE();
					poParent->tableinfo.tab5close1 = (real) poDerivedStates->Get_duration_is_aded1_TRUE();
					poParent->tableinfo.tab5close2 = (real) poDerivedStates->Get_duration_is_aded2_TRUE();
					poParent->tableinfo.tab5close3 = (real) poDerivedStates->Get_duration_is_child_TRUE();
					poParent->tableinfo.tab5close4 = (real) poDerivedStates->Get_duration_is_old_TRUE();
					poParent->tableinfo.tab5cell = (index) poParent->TableCellIndex( 5, 3,
						poParent->geo_nat,
						poParent->sex,
						OldValue - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab5total_cell = (index) poParent->tableinfo.tab5cell;
					}
					else {
						poParent->tableinfo.tab5total_cell = (index) poParent->TableTotalCellIndex( 5, 3,
						poParent->geo_nat,
						poParent->sex,
						OldValue - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab5dim_filter = poParent->table_filter_TabPopProvAgeEducSex;
				}
			}
			if ( gplTableIndexes[7] >= 0 && ( poParent->tableinfo.tab7dim_filter || poParent->table_filter_TabEduc15ByDistrict ) ) {
				if ( poParent->tableinfo.tab7events == -1 || poParent->tableinfo.tab7events < poParent->__events ) {
					if ( poParent->tableinfo.tab7events != -1 ) {
						poParent->tableinfo.tab7open0 = (real) poParent->CloseUpdate( 7, 0, poParent->tableinfo.tab7open0, poParent->tableinfo.tab7close0, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab7open1 = (real) poParent->CloseUpdate( 7, 1, poParent->tableinfo.tab7open1, poParent->tableinfo.tab7close1, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab7open2 = (real) poParent->CloseUpdate( 7, 2, poParent->tableinfo.tab7open2, poParent->tableinfo.tab7close2, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab7open3 = (real) poParent->CloseUpdate( 7, 3, poParent->tableinfo.tab7open3, poParent->tableinfo.tab7close3, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab7events = poParent->__events;
					poParent->tableinfo.tab7close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab7close1 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH();
					poParent->tableinfo.tab7close2 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_LOW();
					poParent->tableinfo.tab7close3 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM();
					poParent->tableinfo.tab7cell = (index) poParent->TableCellIndex( 7, 3,
						poParent->geo_nat,
						poParent->sex,
						OldValue - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab7total_cell = (index) poParent->tableinfo.tab7cell;
					}
					else {
						poParent->tableinfo.tab7total_cell = (index) poParent->TableTotalCellIndex( 7, 3,
						poParent->geo_nat,
						poParent->sex,
						OldValue - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab7dim_filter = poParent->table_filter_TabEduc15ByDistrict;
				}
			}
			if ( gplTableIndexes[8] >= 0 && ( poParent->tableinfo.tab8dim_filter || poParent->table_filter_TabEduc15ByDistrictBirth ) ) {
				if ( poParent->tableinfo.tab8events == -1 || poParent->tableinfo.tab8events < poParent->__events ) {
					if ( poParent->tableinfo.tab8events != -1 ) {
						poParent->tableinfo.tab8open0 = (real) poParent->CloseUpdate( 8, 0, poParent->tableinfo.tab8open0, poParent->tableinfo.tab8close0, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab8open1 = (real) poParent->CloseUpdate( 8, 1, poParent->tableinfo.tab8open1, poParent->tableinfo.tab8close1, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab8open2 = (real) poParent->CloseUpdate( 8, 2, poParent->tableinfo.tab8open2, poParent->tableinfo.tab8close2, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab8open3 = (real) poParent->CloseUpdate( 8, 3, poParent->tableinfo.tab8open3, poParent->tableinfo.tab8close3, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab8events = poParent->__events;
					poParent->tableinfo.tab8close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab8close1 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH();
					poParent->tableinfo.tab8close2 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_LOW();
					poParent->tableinfo.tab8close3 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM();
					poParent->tableinfo.tab8cell = (index) poParent->TableCellIndex( 8, 3,
						poParent->geo_birth,
						poParent->sex,
						OldValue - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab8total_cell = (index) poParent->tableinfo.tab8cell;
					}
					else {
						poParent->tableinfo.tab8total_cell = (index) poParent->TableTotalCellIndex( 8, 3,
						poParent->geo_birth,
						poParent->sex,
						OldValue - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab8dim_filter = poParent->table_filter_TabEduc15ByDistrictBirth;
				}
			}
			if ( gplTableIndexes[9] >= 0 && ( poParent->tableinfo.tab9dim_filter || poParent->table_filter_TabPrimSchoolEntries ) ) {
				if ( poParent->tableinfo.tab9events == -1 || poParent->tableinfo.tab9events < poParent->__events ) {
					if ( poParent->tableinfo.tab9events != -1 ) {
						poParent->tableinfo.tab9open0 = (real) poParent->CloseUpdate( 9, 0, poParent->tableinfo.tab9open0, poParent->tableinfo.tab9close0, poParent->tableinfo.tab9dim_filter, poParent->tableinfo.tab9cell, poParent->tableinfo.tab9total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab9events = poParent->__events;
					poParent->tableinfo.tab9close0 = (real) poDerivedStates->entrances_educ_one_grade_attended_1;
					poParent->tableinfo.tab9cell = (index) poParent->TableCellIndex( 9, 2,
						OldValue - r_min_SIM_YEAR_RANGE,
						poParent->geo_nat );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab9total_cell = (index) poParent->tableinfo.tab9cell;
					}
					else {
						poParent->tableinfo.tab9total_cell = (index) poParent->TableTotalCellIndex( 9, 2,
						OldValue - r_min_SIM_YEAR_RANGE,
						poParent->geo_nat );
					}
					poParent->tableinfo.tab9dim_filter = poParent->table_filter_TabPrimSchoolEntries;
				}
			}
			if ( gplTableIndexes[10] >= 0 && ( poParent->tableinfo.tab10dim_filter || poParent->table_filter_TabPrimSchoolGraduations ) ) {
				if ( poParent->tableinfo.tab10events == -1 || poParent->tableinfo.tab10events < poParent->__events ) {
					if ( poParent->tableinfo.tab10events != -1 ) {
						poParent->tableinfo.tab10open0 = (real) poParent->CloseUpdate( 10, 0, poParent->tableinfo.tab10open0, poParent->tableinfo.tab10close0, poParent->tableinfo.tab10dim_filter, poParent->tableinfo.tab10cell, poParent->tableinfo.tab10total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab10events = poParent->__events;
					poParent->tableinfo.tab10close0 = (real) poDerivedStates->entrances_has_passed_primary_TRUE;
					poParent->tableinfo.tab10cell = (index) poParent->TableCellIndex( 10, 2,
						OldValue - r_min_SIM_YEAR_RANGE,
						poParent->geo_nat );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab10total_cell = (index) poParent->tableinfo.tab10cell;
					}
					else {
						poParent->tableinfo.tab10total_cell = (index) poParent->TableTotalCellIndex( 10, 2,
						OldValue - r_min_SIM_YEAR_RANGE,
						poParent->geo_nat );
					}
					poParent->tableinfo.tab10dim_filter = poParent->table_filter_TabPrimSchoolGraduations;
				}
			}
			if ( gplTableIndexes[11] >= 0 && ( poParent->tableinfo.tab11dim_filter || poParent->table_filter_TabPrimSchoolOutOfSchool9to11 ) ) {
				if ( poParent->tableinfo.tab11events == -1 || poParent->tableinfo.tab11events < poParent->__events ) {
					if ( poParent->tableinfo.tab11events != -1 ) {
						poParent->tableinfo.tab11open0 = (real) poParent->CloseUpdate( 11, 0, poParent->tableinfo.tab11open0, poParent->tableinfo.tab11close0, poParent->tableinfo.tab11dim_filter, poParent->tableinfo.tab11cell, poParent->tableinfo.tab11total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab11open1 = (real) poParent->CloseUpdate( 11, 1, poParent->tableinfo.tab11open1, poParent->tableinfo.tab11close1, poParent->tableinfo.tab11dim_filter, poParent->tableinfo.tab11cell, poParent->tableinfo.tab11total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab11events = poParent->__events;
					poParent->tableinfo.tab11close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab11close1 = (real) poDerivedStates->Get_duration_tab_out_of_primary_TRUE();
					poParent->tableinfo.tab11cell = (index) poParent->TableCellIndex( 11, 2,
						OldValue - r_min_SIM_YEAR_RANGE,
						poParent->geo_nat );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab11total_cell = (index) poParent->tableinfo.tab11cell;
					}
					else {
						poParent->tableinfo.tab11total_cell = (index) poParent->TableTotalCellIndex( 11, 2,
						OldValue - r_min_SIM_YEAR_RANGE,
						poParent->geo_nat );
					}
					poParent->tableinfo.tab11dim_filter = poParent->table_filter_TabPrimSchoolOutOfSchool9to11;
				}
			}
			if ( gplTableIndexes[12] >= 0 && ( poParent->tableinfo.tab12dim_filter || poParent->table_filter_TabSchool2TrackTab ) ) {
				if ( poParent->tableinfo.tab12events == -1 || poParent->tableinfo.tab12events < poParent->__events ) {
					if ( poParent->tableinfo.tab12events != -1 ) {
						Set_table_unit_TabSchool2TrackTab(poParent->table_unit_TabSchool2TrackTab+1);
						poParent->tableinfo.tab12close0 = (real) poParent->table_unit_TabSchool2TrackTab;
						poParent->tableinfo.tab12open0 = (real) poParent->CloseUpdate( 12, 0, poParent->tableinfo.tab12open0, poParent->tableinfo.tab12close0, poParent->tableinfo.tab12dim_filter, poParent->tableinfo.tab12cell, poParent->tableinfo.tab12total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab12events = poParent->__events;
					poParent->tableinfo.tab12close0 = (real) poParent->table_unit_TabSchool2TrackTab;
					poParent->tableinfo.tab12cell = (index) poParent->TableCellIndex( 12, 2,
						OldValue - r_min_SIM_YEAR_RANGE,
						poParent->educ_two_grade_attended - r_min_EDUC_TWO_GRADE0 );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab12total_cell = (index) poParent->tableinfo.tab12cell;
					}
					else {
						poParent->tableinfo.tab12total_cell = (index) poParent->TableTotalCellIndex( 12, 2,
						OldValue - r_min_SIM_YEAR_RANGE,
						poParent->educ_two_grade_attended - r_min_EDUC_TWO_GRADE0 );
					}
					poParent->tableinfo.tab12dim_filter = poParent->table_filter_TabSchool2TrackTab;
				}
			}
			if ( gplTableIndexes[13] >= 0 && ( poParent->tableinfo.tab13dim_filter || poParent->table_filter_TabSchool2AttainmentsTab ) ) {
				if ( poParent->tableinfo.tab13events == -1 || poParent->tableinfo.tab13events < poParent->__events ) {
					if ( poParent->tableinfo.tab13events != -1 ) {
						poParent->tableinfo.tab13open0 = (real) poParent->CloseUpdate( 13, 0, poParent->tableinfo.tab13open0, poParent->tableinfo.tab13close0, poParent->tableinfo.tab13dim_filter, poParent->tableinfo.tab13cell, poParent->tableinfo.tab13total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab13open1 = (real) poParent->CloseUpdate( 13, 1, poParent->tableinfo.tab13open1, poParent->tableinfo.tab13close1, poParent->tableinfo.tab13dim_filter, poParent->tableinfo.tab13cell, poParent->tableinfo.tab13total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab13open2 = (real) poParent->CloseUpdate( 13, 2, poParent->tableinfo.tab13open2, poParent->tableinfo.tab13close2, poParent->tableinfo.tab13dim_filter, poParent->tableinfo.tab13cell, poParent->tableinfo.tab13total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab13open3 = (real) poParent->CloseUpdate( 13, 3, poParent->tableinfo.tab13open3, poParent->tableinfo.tab13close3, poParent->tableinfo.tab13dim_filter, poParent->tableinfo.tab13cell, poParent->tableinfo.tab13total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab13events = poParent->__events;
					poParent->tableinfo.tab13close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab13close1 = (real) poDerivedStates->Get_duration_educ_two_level_ETL_GRADUATE();
					poParent->tableinfo.tab13close2 = (real) poDerivedStates->Get_duration_educ_two_level_ETL_INCOMPLETE();
					poParent->tableinfo.tab13close3 = (real) poDerivedStates->Get_duration_educ_two_level_ETL_NEVER();
					poParent->tableinfo.tab13cell = (index) poParent->TableCellIndex( 13, 2,
						poParent->sex,
						OldValue - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab13total_cell = (index) poParent->tableinfo.tab13cell;
					}
					else {
						poParent->tableinfo.tab13total_cell = (index) poParent->TableTotalCellIndex( 13, 2,
						poParent->sex,
						OldValue - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab13dim_filter = poParent->table_filter_TabSchool2AttainmentsTab;
				}
			}
			if ( gplTableIndexes[14] >= 0 && ( poParent->tableinfo.tab14dim_filter || poParent->table_filter_tabHavingSpouse ) ) {
				if ( poParent->tableinfo.tab14events == -1 || poParent->tableinfo.tab14events < poParent->__events ) {
					if ( poParent->tableinfo.tab14events != -1 ) {
						poParent->tableinfo.tab14open0 = (real) poParent->CloseUpdate( 14, 0, poParent->tableinfo.tab14open0, poParent->tableinfo.tab14close0, poParent->tableinfo.tab14dim_filter, poParent->tableinfo.tab14cell, poParent->tableinfo.tab14total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab14open1 = (real) poParent->CloseUpdate( 14, 1, poParent->tableinfo.tab14open1, poParent->tableinfo.tab14close1, poParent->tableinfo.tab14dim_filter, poParent->tableinfo.tab14cell, poParent->tableinfo.tab14total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab14open2 = (real) poParent->CloseUpdate( 14, 2, poParent->tableinfo.tab14open2, poParent->tableinfo.tab14close2, poParent->tableinfo.tab14dim_filter, poParent->tableinfo.tab14cell, poParent->tableinfo.tab14total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab14events = poParent->__events;
					poParent->tableinfo.tab14close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab14close1 = (real) poDerivedStates->Get_duration_ever_union_TRUE();
					poParent->tableinfo.tab14close2 = (real) poDerivedStates->Get_duration_has_spouse_TRUE();
					poParent->tableinfo.tab14cell = (index) poParent->TableCellIndex( 14, 3,
						poParent->sex,
						poParent->integer_age - r_min_AGE_RANGE,
						OldValue - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab14total_cell = (index) poParent->tableinfo.tab14cell;
					}
					else {
						poParent->tableinfo.tab14total_cell = (index) poParent->TableTotalCellIndex( 14, 3,
						poParent->sex,
						poParent->integer_age - r_min_AGE_RANGE,
						OldValue - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab14dim_filter = poParent->table_filter_tabHavingSpouse;
				}
			}
			if ( gplTableIndexes[17] >= 0 && ( poParent->tableinfo.tab17dim_filter || poParent->table_filter_tabPopulationYearPlace ) ) {
				if ( poParent->tableinfo.tab17events == -1 || poParent->tableinfo.tab17events < poParent->__events ) {
					if ( poParent->tableinfo.tab17events != -1 ) {
						poParent->tableinfo.tab17open0 = (real) poParent->CloseUpdate( 17, 0, poParent->tableinfo.tab17open0, poParent->tableinfo.tab17close0, poParent->tableinfo.tab17dim_filter, poParent->tableinfo.tab17cell, poParent->tableinfo.tab17total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab17events = poParent->__events;
					poParent->tableinfo.tab17close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab17cell = (index) poParent->TableCellIndex( 17, 2,
						poParent->geo,
						OldValue - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab17total_cell = (index) poParent->tableinfo.tab17cell;
					}
					else {
						poParent->tableinfo.tab17total_cell = (index) poParent->TableTotalCellIndex( 17, 2,
						poParent->geo,
						OldValue - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab17dim_filter = poParent->table_filter_tabPopulationYearPlace;
				}
			}
			if ( gplTableIndexes[18] >= 0 && ( poParent->tableinfo.tab18dim_filter || poParent->table_filter_tabBirthsYearPlace ) ) {
				if ( poParent->tableinfo.tab18events == -1 || poParent->tableinfo.tab18events < poParent->__events ) {
					if ( poParent->tableinfo.tab18events != -1 ) {
						poParent->tableinfo.tab18open0 = (real) poParent->CloseUpdate( 18, 0, poParent->tableinfo.tab18open0, poParent->tableinfo.tab18close0, poParent->tableinfo.tab18dim_filter, poParent->tableinfo.tab18cell, poParent->tableinfo.tab18total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab18events = poParent->__events;
					poParent->tableinfo.tab18close0 = (real) poDerivedStates->transitions_is_alive_FALSE_TRUE;
					poParent->tableinfo.tab18cell = (index) poParent->TableCellIndex( 18, 2,
						poParent->geo,
						OldValue - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab18total_cell = (index) poParent->tableinfo.tab18cell;
					}
					else {
						poParent->tableinfo.tab18total_cell = (index) poParent->TableTotalCellIndex( 18, 2,
						poParent->geo,
						OldValue - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab18dim_filter = poParent->table_filter_tabBirthsYearPlace;
				}
			}
			if ( gplTableIndexes[20] >= 0 && ( poParent->tableinfo.tab20dim_filter || poParent->table_filter_tabStuntingSexRegYob ) ) {
				if ( poParent->tableinfo.tab20events == -1 || poParent->tableinfo.tab20events < poParent->__events ) {
					if ( poParent->tableinfo.tab20events != -1 ) {
						poParent->tableinfo.tab20open0 = (real) poParent->CloseUpdate( 20, 0, poParent->tableinfo.tab20open0, poParent->tableinfo.tab20close0, poParent->tableinfo.tab20dim_filter, poParent->tableinfo.tab20cell, poParent->tableinfo.tab20total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						Set_table_unit_tabStuntingSexRegYob(poParent->table_unit_tabStuntingSexRegYob+1);
						poParent->tableinfo.tab20close1 = (real) poParent->table_unit_tabStuntingSexRegYob;
						poParent->tableinfo.tab20open1 = (real) poParent->CloseUpdate( 20, 1, poParent->tableinfo.tab20open1, poParent->tableinfo.tab20close1, poParent->tableinfo.tab20dim_filter, poParent->tableinfo.tab20cell, poParent->tableinfo.tab20total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab20events = poParent->__events;
					poParent->tableinfo.tab20close0 = (real) poParent->is_stunted_value_in;
					poParent->tableinfo.tab20close1 = (real) poParent->table_unit_tabStuntingSexRegYob;
					poParent->tableinfo.tab20cell = (index) poParent->TableCellIndex( 20, 2,
						poParent->region_nat,
						OldValue - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab20total_cell = (index) poParent->tableinfo.tab20cell;
					}
					else {
						poParent->tableinfo.tab20total_cell = (index) poParent->TableTotalCellIndex( 20, 2,
						poParent->region_nat,
						OldValue - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab20dim_filter = poParent->table_filter_tabStuntingSexRegYob;
				}
			}
			if ( gplTableIndexes[21] >= 0 && ( poParent->tableinfo.tab21dim_filter || poParent->table_filter_tabStuntingSexRegMotherYob ) ) {
				if ( poParent->tableinfo.tab21events == -1 || poParent->tableinfo.tab21events < poParent->__events ) {
					if ( poParent->tableinfo.tab21events != -1 ) {
						poParent->tableinfo.tab21open0 = (real) poParent->CloseUpdate( 21, 0, poParent->tableinfo.tab21open0, poParent->tableinfo.tab21close0, poParent->tableinfo.tab21dim_filter, poParent->tableinfo.tab21cell, poParent->tableinfo.tab21total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						Set_table_unit_tabStuntingSexRegMotherYob(poParent->table_unit_tabStuntingSexRegMotherYob+1);
						poParent->tableinfo.tab21close1 = (real) poParent->table_unit_tabStuntingSexRegMotherYob;
						poParent->tableinfo.tab21open1 = (real) poParent->CloseUpdate( 21, 1, poParent->tableinfo.tab21open1, poParent->tableinfo.tab21close1, poParent->tableinfo.tab21dim_filter, poParent->tableinfo.tab21cell, poParent->tableinfo.tab21total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab21events = poParent->__events;
					poParent->tableinfo.tab21close0 = (real) poParent->is_stunted_value_in;
					poParent->tableinfo.tab21close1 = (real) poParent->table_unit_tabStuntingSexRegMotherYob;
					poParent->tableinfo.tab21cell = (index) poParent->TableCellIndex( 21, 4,
						poParent->sex,
						poParent->educ_mother,
						poParent->region_nat,
						OldValue - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab21total_cell = (index) poParent->tableinfo.tab21cell;
					}
					else {
						poParent->tableinfo.tab21total_cell = (index) poParent->TableTotalCellIndex( 21, 4,
						poParent->sex,
						poParent->educ_mother,
						poParent->region_nat,
						OldValue - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab21dim_filter = poParent->table_filter_tabStuntingSexRegMotherYob;
				}
			}
			if ( gplTableIndexes[22] >= 0 && ( poParent->tableinfo.tab22dim_filter || poParent->table_filter_TabPrenatCare ) ) {
				if ( poParent->tableinfo.tab22events == -1 || poParent->tableinfo.tab22events < poParent->__events ) {
					if ( poParent->tableinfo.tab22events != -1 ) {
						poParent->tableinfo.tab22open0 = (real) poParent->CloseUpdate( 22, 0, poParent->tableinfo.tab22open0, poParent->tableinfo.tab22close0, poParent->tableinfo.tab22dim_filter, poParent->tableinfo.tab22cell, poParent->tableinfo.tab22total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab22open1 = (real) poParent->CloseUpdate( 22, 1, poParent->tableinfo.tab22open1, poParent->tableinfo.tab22close1, poParent->tableinfo.tab22dim_filter, poParent->tableinfo.tab22cell, poParent->tableinfo.tab22total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab22events = poParent->__events;
					poParent->tableinfo.tab22close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab22close1 = (real) poDerivedStates->Get_duration_got_prenat_care_GPC_YES();
					poParent->tableinfo.tab22cell = (index) poParent->TableCellIndex( 22, 3,
						poParent->region_nat,
						poParent->educ_mother,
						OldValue - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab22total_cell = (index) poParent->tableinfo.tab22cell;
					}
					else {
						poParent->tableinfo.tab22total_cell = (index) poParent->TableTotalCellIndex( 22, 3,
						poParent->region_nat,
						poParent->educ_mother,
						OldValue - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab22dim_filter = poParent->table_filter_TabPrenatCare;
				}
			}
			if ( gplTableIndexes[23] >= 0 && ( poParent->tableinfo.tab23dim_filter || poParent->table_filter_TabImmunizationChildren ) ) {
				if ( poParent->tableinfo.tab23events == -1 || poParent->tableinfo.tab23events < poParent->__events ) {
					if ( poParent->tableinfo.tab23events != -1 ) {
						poParent->tableinfo.tab23open0 = (real) poParent->CloseUpdate( 23, 0, poParent->tableinfo.tab23open0, poParent->tableinfo.tab23close0, poParent->tableinfo.tab23dim_filter, poParent->tableinfo.tab23cell, poParent->tableinfo.tab23total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab23open1 = (real) poParent->CloseUpdate( 23, 1, poParent->tableinfo.tab23open1, poParent->tableinfo.tab23close1, poParent->tableinfo.tab23dim_filter, poParent->tableinfo.tab23cell, poParent->tableinfo.tab23total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab23events = poParent->__events;
					poParent->tableinfo.tab23close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab23close1 = (real) poDerivedStates->Get_duration_is_immunized_TRUE();
					poParent->tableinfo.tab23cell = (index) poParent->TableCellIndex( 23, 6,
						poParent->sex,
						poParent->ethnicity_short,
						poParent->got_prenat_care,
						poParent->region_nat,
						poParent->educ_mother,
						OldValue - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab23total_cell = (index) poParent->tableinfo.tab23cell;
					}
					else {
						poParent->tableinfo.tab23total_cell = (index) poParent->TableTotalCellIndex( 23, 6,
						poParent->sex,
						poParent->ethnicity_short,
						poParent->got_prenat_care,
						poParent->region_nat,
						poParent->educ_mother,
						OldValue - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab23dim_filter = poParent->table_filter_TabImmunizationChildren;
				}
			}
			if ( gplTableIndexes[24] >= 0 && ( poParent->tableinfo.tab24dim_filter || poParent->table_filter_TabChildVaccination ) ) {
				if ( poParent->tableinfo.tab24events == -1 || poParent->tableinfo.tab24events < poParent->__events ) {
					if ( poParent->tableinfo.tab24events != -1 ) {
						poParent->tableinfo.tab24open0 = (real) poParent->CloseUpdate( 24, 0, poParent->tableinfo.tab24open0, poParent->tableinfo.tab24close0, poParent->tableinfo.tab24dim_filter, poParent->tableinfo.tab24cell, poParent->tableinfo.tab24total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab24open1 = (real) poParent->CloseUpdate( 24, 1, poParent->tableinfo.tab24open1, poParent->tableinfo.tab24close1, poParent->tableinfo.tab24dim_filter, poParent->tableinfo.tab24cell, poParent->tableinfo.tab24total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab24events = poParent->__events;
					poParent->tableinfo.tab24close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab24close1 = (real) poDerivedStates->Get_duration_is_immunized_TRUE();
					poParent->tableinfo.tab24cell = (index) poParent->TableCellIndex( 24, 2,
						OldValue - r_min_SIM_YEAR_RANGE,
						poParent->region_nat );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab24total_cell = (index) poParent->tableinfo.tab24cell;
					}
					else {
						poParent->tableinfo.tab24total_cell = (index) poParent->TableTotalCellIndex( 24, 2,
						OldValue - r_min_SIM_YEAR_RANGE,
						poParent->region_nat );
					}
					poParent->tableinfo.tab24dim_filter = poParent->table_filter_TabChildVaccination;
				}
			}
		}
	}
	return poParent->sim_year = value;
}

SPOUSE_GROUP PersonStateFunctions::Set_spouse_group( SPOUSE_GROUP value )
{

	CHKLMT(value, spouse_group, 0, rc_size_SPOUSE_GROUP - 1);
	if ( value < 0 ) {
		value = (SPOUSE_GROUP) 0;
	}
	else if ( value >= rc_size_SPOUSE_GROUP ) {
		value = (SPOUSE_GROUP) ( rc_size_SPOUSE_GROUP - 1 );
	}
	if ( poParent->spouse_group != value )
	{


		SPOUSE_GROUP	OldValue;
		OldValue = poParent->spouse_group;
		poParent->spouse_group = value;

		// actorset dimension or order state
		if (poParent->actorset_filter_asAvailableMaleByType)
		{
			asAvailableMaleByType[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][OldValue][poParent->geo]->Remove(poParent->poasAvailableMaleByType);
			poParent->poasAvailableMaleByType = asAvailableMaleByType[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->spouse_group][poParent->geo]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asAvailableMaleByTypeImmiScratch)
		{
			asAvailableMaleByTypeImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][OldValue][poParent->geo][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->Remove(poParent->poasAvailableMaleByTypeImmiScratch);
			poParent->poasAvailableMaleByTypeImmiScratch = asAvailableMaleByTypeImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->spouse_group][poParent->geo][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->NewNode(poParent);
		}
	}
	return poParent->spouse_group = value;
}

logical PersonStateFunctions::Set_survived_early_years( logical value )
{

	value = value != 0;
	if ( poParent->survived_early_years != value )
	{
		if ( poParent->survived_early_years == FALSE && value == TRUE ) {
			poDerivedStates->Set_transitions_survived_early_years_FALSE_TRUE( poDerivedStates->transitions_survived_early_years_FALSE_TRUE + 1 );
		}


		logical	OldValue;
		OldValue = poParent->survived_early_years;
		poParent->survived_early_years = value;

	}
	return poParent->survived_early_years = value;
}

logical PersonStateFunctions::Set_tab_out_of_primary( logical value )
{

	value = value != 0;
	if ( poParent->tab_out_of_primary != value )
	{
		if ( poParent->tab_out_of_primary == TRUE ) {
			poDerivedStates->Set_duration_tab_out_of_primary_TRUE( poDerivedStates->Get_duration_tab_out_of_primary_TRUE() );
		}
		else if ( value == TRUE ) {
			poDerivedStates->duration_tab_out_of_primary_TRUE_time = poParent->__time;
		}


		logical	OldValue;
		OldValue = poParent->tab_out_of_primary;
		poParent->tab_out_of_primary = value;

	}
	return poParent->tab_out_of_primary = value;
}

SIM_YEAR_RANGE PersonStateFunctions::Set_tab_sim_yob( SIM_YEAR_RANGE value )
{

	CHKLMT(value, tab_sim_yob, r_min_SIM_YEAR_RANGE, r_max_SIM_YEAR_RANGE);
	if ( value < r_min_SIM_YEAR_RANGE ) {
		value = r_min_SIM_YEAR_RANGE;
	}
	else if ( value > r_max_SIM_YEAR_RANGE ) {
		value = r_max_SIM_YEAR_RANGE;
	}
	if ( poParent->tab_sim_yob != value )
	{


		SIM_YEAR_RANGE	OldValue;
		OldValue = poParent->tab_sim_yob;
		poParent->tab_sim_yob = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[0] >= 0 && ( poParent->tableinfo.tab0dim_filter || poParent->table_filter_tabPreSchool ) ) {
				if ( poParent->tableinfo.tab0events == -1 || poParent->tableinfo.tab0events < poParent->__events ) {
					if ( poParent->tableinfo.tab0events != -1 ) {
						poParent->tableinfo.tab0open0 = (real) poParent->CloseUpdate( 0, 0, poParent->tableinfo.tab0open0, poParent->tableinfo.tab0close0, poParent->tableinfo.tab0dim_filter, poParent->tableinfo.tab0cell, poParent->tableinfo.tab0total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab0open1 = (real) poParent->CloseUpdate( 0, 1, poParent->tableinfo.tab0open1, poParent->tableinfo.tab0close1, poParent->tableinfo.tab0dim_filter, poParent->tableinfo.tab0cell, poParent->tableinfo.tab0total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab0open2 = (real) poParent->CloseUpdate( 0, 2, poParent->tableinfo.tab0open2, poParent->tableinfo.tab0close2, poParent->tableinfo.tab0dim_filter, poParent->tableinfo.tab0cell, poParent->tableinfo.tab0total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab0open3 = (real) poParent->CloseUpdate( 0, 3, poParent->tableinfo.tab0open3, poParent->tableinfo.tab0close3, poParent->tableinfo.tab0dim_filter, poParent->tableinfo.tab0cell, poParent->tableinfo.tab0total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab0events = poParent->__events;
					poParent->tableinfo.tab0close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab0close1 = (real) poDerivedStates->Get_duration_years_preschool_0();
					poParent->tableinfo.tab0close2 = (real) poDerivedStates->Get_duration_years_preschool_1();
					poParent->tableinfo.tab0close3 = (real) poDerivedStates->Get_duration_years_preschool_2();
					poParent->tableinfo.tab0cell = (index) poParent->TableCellIndex( 0, 2,
						poParent->geo,
						OldValue - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab0total_cell = (index) poParent->tableinfo.tab0cell;
					}
					else {
						poParent->tableinfo.tab0total_cell = (index) poParent->TableTotalCellIndex( 0, 2,
						poParent->geo,
						OldValue - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab0dim_filter = poParent->table_filter_tabPreSchool;
				}
			}
			if ( gplTableIndexes[3] >= 0 && ( poParent->tableinfo.tab3dim_filter || poParent->table_filter_TabEducFateByGroup ) ) {
				if ( poParent->tableinfo.tab3events == -1 || poParent->tableinfo.tab3events < poParent->__events ) {
					if ( poParent->tableinfo.tab3events != -1 ) {
						Set_table_unit_TabEducFateByGroup(poParent->table_unit_TabEducFateByGroup+1);
						poParent->tableinfo.tab3close0 = (real) poParent->table_unit_TabEducFateByGroup;
						poParent->tableinfo.tab3open0 = (real) poParent->CloseUpdate( 3, 0, poParent->tableinfo.tab3open0, poParent->tableinfo.tab3close0, poParent->tableinfo.tab3dim_filter, poParent->tableinfo.tab3cell, poParent->tableinfo.tab3total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab3events = poParent->__events;
					poParent->tableinfo.tab3close0 = (real) poParent->table_unit_TabEducFateByGroup;
					poParent->tableinfo.tab3cell = (index) poParent->TableCellIndex( 3, 3,
						poParent->educ1_group,
						OldValue - r_min_SIM_YEAR_RANGE,
						poParent->educ_one_fate );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab3total_cell = (index) poParent->tableinfo.tab3cell;
					}
					else {
						poParent->tableinfo.tab3total_cell = (index) poParent->TableTotalCellIndex( 3, 3,
						poParent->educ1_group,
						OldValue - r_min_SIM_YEAR_RANGE,
						poParent->educ_one_fate );
					}
					poParent->tableinfo.tab3dim_filter = poParent->table_filter_TabEducFateByGroup;
				}
			}
			if ( gplTableIndexes[15] >= 0 && ( poParent->tableinfo.tab15dim_filter || poParent->table_filter_tabHCI ) ) {
				if ( poParent->tableinfo.tab15events == -1 || poParent->tableinfo.tab15events < poParent->__events ) {
					if ( poParent->tableinfo.tab15events != -1 ) {
						poParent->tableinfo.tab15open0 = (real) poParent->CloseUpdate( 15, 0, poParent->tableinfo.tab15open0, poParent->tableinfo.tab15close0, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open1 = (real) poParent->CloseUpdate( 15, 1, poParent->tableinfo.tab15open1, poParent->tableinfo.tab15close1, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open2 = (real) poParent->CloseUpdate( 15, 2, poParent->tableinfo.tab15open2, poParent->tableinfo.tab15close2, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open3 = (real) poParent->CloseUpdate( 15, 3, poParent->tableinfo.tab15open3, poParent->tableinfo.tab15close3, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open4 = (real) poParent->CloseUpdate( 15, 4, poParent->tableinfo.tab15open4, poParent->tableinfo.tab15close4, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open5 = (real) poParent->CloseUpdate( 15, 5, poParent->tableinfo.tab15open5, poParent->tableinfo.tab15close5, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open6 = (real) poParent->CloseUpdate( 15, 6, poParent->tableinfo.tab15open6, poParent->tableinfo.tab15close6, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open7 = (real) poParent->CloseUpdate( 15, 7, poParent->tableinfo.tab15open7, poParent->tableinfo.tab15close7, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open8 = (real) poParent->CloseUpdate( 15, 8, poParent->tableinfo.tab15open8, poParent->tableinfo.tab15close8, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab15events = poParent->__events;
					poParent->tableinfo.tab15close0 = (real) poDerivedStates->transitions_integer_age_14_15;
					poParent->tableinfo.tab15close1 = (real) poDerivedStates->transitions_is_alive_TRUE_FALSE;
					poParent->tableinfo.tab15close2 = (real) poDerivedStates->transitions_survived_early_years_FALSE_TRUE;
					poParent->tableinfo.tab15close3 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_adult_survival;
					poParent->tableinfo.tab15close4 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_ind_hci;
					poParent->tableinfo.tab15close5 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_is_stunted;
					poParent->tableinfo.tab15close6 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling;
					poParent->tableinfo.tab15close7 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_survived_early_years;
					poParent->tableinfo.tab15close8 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling;
					poParent->tableinfo.tab15cell = (index) poParent->TableCellIndex( 15, 3,
						poParent->sex,
						poDerivedStates->aggregate_region_birth_REGION_NAT,
						OldValue - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab15total_cell = (index) poParent->tableinfo.tab15cell;
					}
					else {
						poParent->tableinfo.tab15total_cell = (index) poParent->TableTotalCellIndex( 15, 3,
						poParent->sex,
						poDerivedStates->aggregate_region_birth_REGION_NAT,
						OldValue - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab15dim_filter = poParent->table_filter_tabHCI;
				}
			}
			if ( gplTableIndexes[16] >= 0 && ( poParent->tableinfo.tab16dim_filter || poParent->table_filter_TabHCIDistrict ) ) {
				if ( poParent->tableinfo.tab16events == -1 || poParent->tableinfo.tab16events < poParent->__events ) {
					if ( poParent->tableinfo.tab16events != -1 ) {
						poParent->tableinfo.tab16open0 = (real) poParent->CloseUpdate( 16, 0, poParent->tableinfo.tab16open0, poParent->tableinfo.tab16close0, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open1 = (real) poParent->CloseUpdate( 16, 1, poParent->tableinfo.tab16open1, poParent->tableinfo.tab16close1, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open2 = (real) poParent->CloseUpdate( 16, 2, poParent->tableinfo.tab16open2, poParent->tableinfo.tab16close2, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open3 = (real) poParent->CloseUpdate( 16, 3, poParent->tableinfo.tab16open3, poParent->tableinfo.tab16close3, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open4 = (real) poParent->CloseUpdate( 16, 4, poParent->tableinfo.tab16open4, poParent->tableinfo.tab16close4, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open5 = (real) poParent->CloseUpdate( 16, 5, poParent->tableinfo.tab16open5, poParent->tableinfo.tab16close5, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open6 = (real) poParent->CloseUpdate( 16, 6, poParent->tableinfo.tab16open6, poParent->tableinfo.tab16close6, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open7 = (real) poParent->CloseUpdate( 16, 7, poParent->tableinfo.tab16open7, poParent->tableinfo.tab16close7, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open8 = (real) poParent->CloseUpdate( 16, 8, poParent->tableinfo.tab16open8, poParent->tableinfo.tab16close8, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab16events = poParent->__events;
					poParent->tableinfo.tab16close0 = (real) poDerivedStates->transitions_integer_age_14_15;
					poParent->tableinfo.tab16close1 = (real) poDerivedStates->transitions_is_alive_TRUE_FALSE;
					poParent->tableinfo.tab16close2 = (real) poDerivedStates->transitions_survived_early_years_FALSE_TRUE;
					poParent->tableinfo.tab16close3 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_adult_survival;
					poParent->tableinfo.tab16close4 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_ind_hci;
					poParent->tableinfo.tab16close5 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_is_stunted;
					poParent->tableinfo.tab16close6 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling;
					poParent->tableinfo.tab16close7 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_survived_early_years;
					poParent->tableinfo.tab16close8 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling;
					poParent->tableinfo.tab16cell = (index) poParent->TableCellIndex( 16, 3,
						poParent->sex,
						poDerivedStates->aggregate_geo_birth_GEO_NAT,
						OldValue - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab16total_cell = (index) poParent->tableinfo.tab16cell;
					}
					else {
						poParent->tableinfo.tab16total_cell = (index) poParent->TableTotalCellIndex( 16, 3,
						poParent->sex,
						poDerivedStates->aggregate_geo_birth_GEO_NAT,
						OldValue - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab16dim_filter = poParent->table_filter_TabHCIDistrict;
				}
			}
		}
	}
	return poParent->tab_sim_yob = value;
}

logical PersonStateFunctions::Set_table_filter_PopPyramidByEduc( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_PopPyramidByEduc != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_PopPyramidByEduc;
		poParent->table_filter_PopPyramidByEduc = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[2] >= 0 ) {
				if ( poParent->tableinfo.tab2events == -1 || poParent->tableinfo.tab2events < poParent->__events ) {
					if ( poParent->tableinfo.tab2events != -1 ) {
						poParent->tableinfo.tab2open0 = (real) poParent->CloseUpdate( 2, 0, poParent->tableinfo.tab2open0, poParent->tableinfo.tab2close0, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab2open1 = (real) poParent->CloseUpdate( 2, 1, poParent->tableinfo.tab2open1, poParent->tableinfo.tab2close1, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab2open2 = (real) poParent->CloseUpdate( 2, 2, poParent->tableinfo.tab2open2, poParent->tableinfo.tab2close2, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab2events = poParent->__events;
					poParent->tableinfo.tab2close0 = (real) poDerivedStates->Get_duration_educ_one_level_EOL_HIGH();
					poParent->tableinfo.tab2close1 = (real) poDerivedStates->Get_duration_educ_one_level_EOL_LOW();
					poParent->tableinfo.tab2close2 = (real) poDerivedStates->Get_duration_educ_one_level_EOL_MEDIUM();
					poParent->tableinfo.tab2cell = (index) poParent->TableCellIndex( 2, 4,
						poParent->sel_years,
						poParent->sex,
						poParent->geo_nat,
						poDerivedStates->split_integer_age_TAB_AGEC5 );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab2total_cell = (index) poParent->tableinfo.tab2cell;
					}
					else {
						poParent->tableinfo.tab2total_cell = (index) poParent->TableTotalCellIndex( 2, 4,
						poParent->sel_years,
						poParent->sex,
						poParent->geo_nat,
						poDerivedStates->split_integer_age_TAB_AGEC5 );
					}
					poParent->tableinfo.tab2dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_PopPyramidByEduc = value;
}

logical PersonStateFunctions::Set_table_filter_tabBirthsYearPlace( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_tabBirthsYearPlace != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_tabBirthsYearPlace;
		poParent->table_filter_tabBirthsYearPlace = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[18] >= 0 ) {
				if ( poParent->tableinfo.tab18events == -1 || poParent->tableinfo.tab18events < poParent->__events ) {
					if ( poParent->tableinfo.tab18events != -1 ) {
						poParent->tableinfo.tab18open0 = (real) poParent->CloseUpdate( 18, 0, poParent->tableinfo.tab18open0, poParent->tableinfo.tab18close0, poParent->tableinfo.tab18dim_filter, poParent->tableinfo.tab18cell, poParent->tableinfo.tab18total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab18events = poParent->__events;
					poParent->tableinfo.tab18close0 = (real) poDerivedStates->transitions_is_alive_FALSE_TRUE;
					poParent->tableinfo.tab18cell = (index) poParent->TableCellIndex( 18, 2,
						poParent->geo,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab18total_cell = (index) poParent->tableinfo.tab18cell;
					}
					else {
						poParent->tableinfo.tab18total_cell = (index) poParent->TableTotalCellIndex( 18, 2,
						poParent->geo,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab18dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_tabBirthsYearPlace = value;
}

logical PersonStateFunctions::Set_table_filter_TabChildVaccination( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabChildVaccination != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabChildVaccination;
		poParent->table_filter_TabChildVaccination = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[24] >= 0 ) {
				if ( poParent->tableinfo.tab24events == -1 || poParent->tableinfo.tab24events < poParent->__events ) {
					if ( poParent->tableinfo.tab24events != -1 ) {
						poParent->tableinfo.tab24open0 = (real) poParent->CloseUpdate( 24, 0, poParent->tableinfo.tab24open0, poParent->tableinfo.tab24close0, poParent->tableinfo.tab24dim_filter, poParent->tableinfo.tab24cell, poParent->tableinfo.tab24total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab24open1 = (real) poParent->CloseUpdate( 24, 1, poParent->tableinfo.tab24open1, poParent->tableinfo.tab24close1, poParent->tableinfo.tab24dim_filter, poParent->tableinfo.tab24cell, poParent->tableinfo.tab24total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab24events = poParent->__events;
					poParent->tableinfo.tab24close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab24close1 = (real) poDerivedStates->Get_duration_is_immunized_TRUE();
					poParent->tableinfo.tab24cell = (index) poParent->TableCellIndex( 24, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						poParent->region_nat );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab24total_cell = (index) poParent->tableinfo.tab24cell;
					}
					else {
						poParent->tableinfo.tab24total_cell = (index) poParent->TableTotalCellIndex( 24, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						poParent->region_nat );
					}
					poParent->tableinfo.tab24dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabChildVaccination = value;
}

logical PersonStateFunctions::Set_table_filter_TabEduc15ByDistrict( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabEduc15ByDistrict != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabEduc15ByDistrict;
		poParent->table_filter_TabEduc15ByDistrict = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[7] >= 0 ) {
				if ( poParent->tableinfo.tab7events == -1 || poParent->tableinfo.tab7events < poParent->__events ) {
					if ( poParent->tableinfo.tab7events != -1 ) {
						poParent->tableinfo.tab7open0 = (real) poParent->CloseUpdate( 7, 0, poParent->tableinfo.tab7open0, poParent->tableinfo.tab7close0, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab7open1 = (real) poParent->CloseUpdate( 7, 1, poParent->tableinfo.tab7open1, poParent->tableinfo.tab7close1, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab7open2 = (real) poParent->CloseUpdate( 7, 2, poParent->tableinfo.tab7open2, poParent->tableinfo.tab7close2, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab7open3 = (real) poParent->CloseUpdate( 7, 3, poParent->tableinfo.tab7open3, poParent->tableinfo.tab7close3, poParent->tableinfo.tab7dim_filter, poParent->tableinfo.tab7cell, poParent->tableinfo.tab7total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab7events = poParent->__events;
					poParent->tableinfo.tab7close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab7close1 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH();
					poParent->tableinfo.tab7close2 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_LOW();
					poParent->tableinfo.tab7close3 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM();
					poParent->tableinfo.tab7cell = (index) poParent->TableCellIndex( 7, 3,
						poParent->geo_nat,
						poParent->sex,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab7total_cell = (index) poParent->tableinfo.tab7cell;
					}
					else {
						poParent->tableinfo.tab7total_cell = (index) poParent->TableTotalCellIndex( 7, 3,
						poParent->geo_nat,
						poParent->sex,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab7dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabEduc15ByDistrict = value;
}

logical PersonStateFunctions::Set_table_filter_TabEduc15ByDistrictBirth( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabEduc15ByDistrictBirth != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabEduc15ByDistrictBirth;
		poParent->table_filter_TabEduc15ByDistrictBirth = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[8] >= 0 ) {
				if ( poParent->tableinfo.tab8events == -1 || poParent->tableinfo.tab8events < poParent->__events ) {
					if ( poParent->tableinfo.tab8events != -1 ) {
						poParent->tableinfo.tab8open0 = (real) poParent->CloseUpdate( 8, 0, poParent->tableinfo.tab8open0, poParent->tableinfo.tab8close0, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab8open1 = (real) poParent->CloseUpdate( 8, 1, poParent->tableinfo.tab8open1, poParent->tableinfo.tab8close1, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab8open2 = (real) poParent->CloseUpdate( 8, 2, poParent->tableinfo.tab8open2, poParent->tableinfo.tab8close2, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab8open3 = (real) poParent->CloseUpdate( 8, 3, poParent->tableinfo.tab8open3, poParent->tableinfo.tab8close3, poParent->tableinfo.tab8dim_filter, poParent->tableinfo.tab8cell, poParent->tableinfo.tab8total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab8events = poParent->__events;
					poParent->tableinfo.tab8close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab8close1 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH();
					poParent->tableinfo.tab8close2 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_LOW();
					poParent->tableinfo.tab8close3 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM();
					poParent->tableinfo.tab8cell = (index) poParent->TableCellIndex( 8, 3,
						poParent->geo_birth,
						poParent->sex,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab8total_cell = (index) poParent->tableinfo.tab8cell;
					}
					else {
						poParent->tableinfo.tab8total_cell = (index) poParent->TableTotalCellIndex( 8, 3,
						poParent->geo_birth,
						poParent->sex,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab8dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabEduc15ByDistrictBirth = value;
}

logical PersonStateFunctions::Set_table_filter_tabEducationFateGeobirYob( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_tabEducationFateGeobirYob != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_tabEducationFateGeobirYob;
		poParent->table_filter_tabEducationFateGeobirYob = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[1] >= 0 ) {
				if ( poParent->tableinfo.tab1events == -1 || poParent->tableinfo.tab1events < poParent->__events ) {
					if ( poParent->tableinfo.tab1events != -1 ) {
						poParent->tableinfo.tab1open0 = (real) poParent->CloseUpdate( 1, 0, poParent->tableinfo.tab1open0, poParent->tableinfo.tab1close0, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab1open1 = (real) poParent->CloseUpdate( 1, 1, poParent->tableinfo.tab1open1, poParent->tableinfo.tab1close1, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab1open2 = (real) poParent->CloseUpdate( 1, 2, poParent->tableinfo.tab1open2, poParent->tableinfo.tab1close2, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab1open3 = (real) poParent->CloseUpdate( 1, 3, poParent->tableinfo.tab1open3, poParent->tableinfo.tab1close3, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab1events = poParent->__events;
					poParent->tableinfo.tab1close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab1close1 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH();
					poParent->tableinfo.tab1close2 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_LOW();
					poParent->tableinfo.tab1close3 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM();
					poParent->tableinfo.tab1cell = (index) poParent->TableCellIndex( 1, 3,
						poParent->sex,
						poParent->geo_birth,
						poParent->year_of_birth - r_min_ALL_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab1total_cell = (index) poParent->tableinfo.tab1cell;
					}
					else {
						poParent->tableinfo.tab1total_cell = (index) poParent->TableTotalCellIndex( 1, 3,
						poParent->sex,
						poParent->geo_birth,
						poParent->year_of_birth - r_min_ALL_YEAR_RANGE );
					}
					poParent->tableinfo.tab1dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_tabEducationFateGeobirYob = value;
}

logical PersonStateFunctions::Set_table_filter_TabEducFateByGroup( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabEducFateByGroup != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabEducFateByGroup;
		poParent->table_filter_TabEducFateByGroup = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[3] >= 0 ) {
				if ( poParent->tableinfo.tab3events == -1 || poParent->tableinfo.tab3events < poParent->__events ) {
					if ( poParent->tableinfo.tab3events != -1 ) {
						Set_table_unit_TabEducFateByGroup(poParent->table_unit_TabEducFateByGroup+1);
						poParent->tableinfo.tab3close0 = (real) poParent->table_unit_TabEducFateByGroup;
						poParent->tableinfo.tab3open0 = (real) poParent->CloseUpdate( 3, 0, poParent->tableinfo.tab3open0, poParent->tableinfo.tab3close0, poParent->tableinfo.tab3dim_filter, poParent->tableinfo.tab3cell, poParent->tableinfo.tab3total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab3events = poParent->__events;
					poParent->tableinfo.tab3close0 = (real) poParent->table_unit_TabEducFateByGroup;
					poParent->tableinfo.tab3cell = (index) poParent->TableCellIndex( 3, 3,
						poParent->educ1_group,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE,
						poParent->educ_one_fate );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab3total_cell = (index) poParent->tableinfo.tab3cell;
					}
					else {
						poParent->tableinfo.tab3total_cell = (index) poParent->TableTotalCellIndex( 3, 3,
						poParent->educ1_group,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE,
						poParent->educ_one_fate );
					}
					poParent->tableinfo.tab3dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabEducFateByGroup = value;
}

logical PersonStateFunctions::Set_table_filter_TabEducFateDistrYob( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabEducFateDistrYob != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabEducFateDistrYob;
		poParent->table_filter_TabEducFateDistrYob = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[6] >= 0 ) {
				if ( poParent->tableinfo.tab6events == -1 || poParent->tableinfo.tab6events < poParent->__events ) {
					if ( poParent->tableinfo.tab6events != -1 ) {
						Set_table_unit_TabEducFateDistrYob(poParent->table_unit_TabEducFateDistrYob+1);
						poParent->tableinfo.tab6close0 = (real) poParent->table_unit_TabEducFateDistrYob;
						poParent->tableinfo.tab6open0 = (real) poParent->CloseUpdate( 6, 0, poParent->tableinfo.tab6open0, poParent->tableinfo.tab6close0, poParent->tableinfo.tab6dim_filter, poParent->tableinfo.tab6cell, poParent->tableinfo.tab6total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab6events = poParent->__events;
					poParent->tableinfo.tab6close0 = (real) poParent->table_unit_TabEducFateDistrYob;
					poParent->tableinfo.tab6cell = (index) poParent->TableCellIndex( 6, 3,
						poDerivedStates->aggregate_region_birth_REGION_NAT,
						poParent->year_of_birth - r_min_ALL_YEAR_RANGE,
						poParent->educ_one_fate );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab6total_cell = (index) poParent->tableinfo.tab6cell;
					}
					else {
						poParent->tableinfo.tab6total_cell = (index) poParent->TableTotalCellIndex( 6, 3,
						poDerivedStates->aggregate_region_birth_REGION_NAT,
						poParent->year_of_birth - r_min_ALL_YEAR_RANGE,
						poParent->educ_one_fate );
					}
					poParent->tableinfo.tab6dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabEducFateDistrYob = value;
}

logical PersonStateFunctions::Set_table_filter_tabHavingSpouse( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_tabHavingSpouse != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_tabHavingSpouse;
		poParent->table_filter_tabHavingSpouse = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[14] >= 0 ) {
				if ( poParent->tableinfo.tab14events == -1 || poParent->tableinfo.tab14events < poParent->__events ) {
					if ( poParent->tableinfo.tab14events != -1 ) {
						poParent->tableinfo.tab14open0 = (real) poParent->CloseUpdate( 14, 0, poParent->tableinfo.tab14open0, poParent->tableinfo.tab14close0, poParent->tableinfo.tab14dim_filter, poParent->tableinfo.tab14cell, poParent->tableinfo.tab14total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab14open1 = (real) poParent->CloseUpdate( 14, 1, poParent->tableinfo.tab14open1, poParent->tableinfo.tab14close1, poParent->tableinfo.tab14dim_filter, poParent->tableinfo.tab14cell, poParent->tableinfo.tab14total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab14open2 = (real) poParent->CloseUpdate( 14, 2, poParent->tableinfo.tab14open2, poParent->tableinfo.tab14close2, poParent->tableinfo.tab14dim_filter, poParent->tableinfo.tab14cell, poParent->tableinfo.tab14total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab14events = poParent->__events;
					poParent->tableinfo.tab14close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab14close1 = (real) poDerivedStates->Get_duration_ever_union_TRUE();
					poParent->tableinfo.tab14close2 = (real) poDerivedStates->Get_duration_has_spouse_TRUE();
					poParent->tableinfo.tab14cell = (index) poParent->TableCellIndex( 14, 3,
						poParent->sex,
						poParent->integer_age - r_min_AGE_RANGE,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab14total_cell = (index) poParent->tableinfo.tab14cell;
					}
					else {
						poParent->tableinfo.tab14total_cell = (index) poParent->TableTotalCellIndex( 14, 3,
						poParent->sex,
						poParent->integer_age - r_min_AGE_RANGE,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab14dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_tabHavingSpouse = value;
}

logical PersonStateFunctions::Set_table_filter_tabHCI( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_tabHCI != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_tabHCI;
		poParent->table_filter_tabHCI = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[15] >= 0 ) {
				if ( poParent->tableinfo.tab15events == -1 || poParent->tableinfo.tab15events < poParent->__events ) {
					if ( poParent->tableinfo.tab15events != -1 ) {
						poParent->tableinfo.tab15open0 = (real) poParent->CloseUpdate( 15, 0, poParent->tableinfo.tab15open0, poParent->tableinfo.tab15close0, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open1 = (real) poParent->CloseUpdate( 15, 1, poParent->tableinfo.tab15open1, poParent->tableinfo.tab15close1, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open2 = (real) poParent->CloseUpdate( 15, 2, poParent->tableinfo.tab15open2, poParent->tableinfo.tab15close2, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open3 = (real) poParent->CloseUpdate( 15, 3, poParent->tableinfo.tab15open3, poParent->tableinfo.tab15close3, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open4 = (real) poParent->CloseUpdate( 15, 4, poParent->tableinfo.tab15open4, poParent->tableinfo.tab15close4, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open5 = (real) poParent->CloseUpdate( 15, 5, poParent->tableinfo.tab15open5, poParent->tableinfo.tab15close5, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open6 = (real) poParent->CloseUpdate( 15, 6, poParent->tableinfo.tab15open6, poParent->tableinfo.tab15close6, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open7 = (real) poParent->CloseUpdate( 15, 7, poParent->tableinfo.tab15open7, poParent->tableinfo.tab15close7, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open8 = (real) poParent->CloseUpdate( 15, 8, poParent->tableinfo.tab15open8, poParent->tableinfo.tab15close8, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab15events = poParent->__events;
					poParent->tableinfo.tab15close0 = (real) poDerivedStates->transitions_integer_age_14_15;
					poParent->tableinfo.tab15close1 = (real) poDerivedStates->transitions_is_alive_TRUE_FALSE;
					poParent->tableinfo.tab15close2 = (real) poDerivedStates->transitions_survived_early_years_FALSE_TRUE;
					poParent->tableinfo.tab15close3 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_adult_survival;
					poParent->tableinfo.tab15close4 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_ind_hci;
					poParent->tableinfo.tab15close5 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_is_stunted;
					poParent->tableinfo.tab15close6 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling;
					poParent->tableinfo.tab15close7 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_survived_early_years;
					poParent->tableinfo.tab15close8 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling;
					poParent->tableinfo.tab15cell = (index) poParent->TableCellIndex( 15, 3,
						poParent->sex,
						poDerivedStates->aggregate_region_birth_REGION_NAT,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab15total_cell = (index) poParent->tableinfo.tab15cell;
					}
					else {
						poParent->tableinfo.tab15total_cell = (index) poParent->TableTotalCellIndex( 15, 3,
						poParent->sex,
						poDerivedStates->aggregate_region_birth_REGION_NAT,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab15dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_tabHCI = value;
}

logical PersonStateFunctions::Set_table_filter_TabHCIDistrict( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabHCIDistrict != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabHCIDistrict;
		poParent->table_filter_TabHCIDistrict = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[16] >= 0 ) {
				if ( poParent->tableinfo.tab16events == -1 || poParent->tableinfo.tab16events < poParent->__events ) {
					if ( poParent->tableinfo.tab16events != -1 ) {
						poParent->tableinfo.tab16open0 = (real) poParent->CloseUpdate( 16, 0, poParent->tableinfo.tab16open0, poParent->tableinfo.tab16close0, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open1 = (real) poParent->CloseUpdate( 16, 1, poParent->tableinfo.tab16open1, poParent->tableinfo.tab16close1, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open2 = (real) poParent->CloseUpdate( 16, 2, poParent->tableinfo.tab16open2, poParent->tableinfo.tab16close2, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open3 = (real) poParent->CloseUpdate( 16, 3, poParent->tableinfo.tab16open3, poParent->tableinfo.tab16close3, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open4 = (real) poParent->CloseUpdate( 16, 4, poParent->tableinfo.tab16open4, poParent->tableinfo.tab16close4, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open5 = (real) poParent->CloseUpdate( 16, 5, poParent->tableinfo.tab16open5, poParent->tableinfo.tab16close5, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open6 = (real) poParent->CloseUpdate( 16, 6, poParent->tableinfo.tab16open6, poParent->tableinfo.tab16close6, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open7 = (real) poParent->CloseUpdate( 16, 7, poParent->tableinfo.tab16open7, poParent->tableinfo.tab16close7, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open8 = (real) poParent->CloseUpdate( 16, 8, poParent->tableinfo.tab16open8, poParent->tableinfo.tab16close8, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab16events = poParent->__events;
					poParent->tableinfo.tab16close0 = (real) poDerivedStates->transitions_integer_age_14_15;
					poParent->tableinfo.tab16close1 = (real) poDerivedStates->transitions_is_alive_TRUE_FALSE;
					poParent->tableinfo.tab16close2 = (real) poDerivedStates->transitions_survived_early_years_FALSE_TRUE;
					poParent->tableinfo.tab16close3 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_adult_survival;
					poParent->tableinfo.tab16close4 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_ind_hci;
					poParent->tableinfo.tab16close5 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_is_stunted;
					poParent->tableinfo.tab16close6 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling;
					poParent->tableinfo.tab16close7 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_survived_early_years;
					poParent->tableinfo.tab16close8 = (real) poDerivedStates->value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling;
					poParent->tableinfo.tab16cell = (index) poParent->TableCellIndex( 16, 3,
						poParent->sex,
						poDerivedStates->aggregate_geo_birth_GEO_NAT,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab16total_cell = (index) poParent->tableinfo.tab16cell;
					}
					else {
						poParent->tableinfo.tab16total_cell = (index) poParent->TableTotalCellIndex( 16, 3,
						poParent->sex,
						poDerivedStates->aggregate_geo_birth_GEO_NAT,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab16dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabHCIDistrict = value;
}

logical PersonStateFunctions::Set_table_filter_TabImmunizationChildren( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabImmunizationChildren != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabImmunizationChildren;
		poParent->table_filter_TabImmunizationChildren = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[23] >= 0 ) {
				if ( poParent->tableinfo.tab23events == -1 || poParent->tableinfo.tab23events < poParent->__events ) {
					if ( poParent->tableinfo.tab23events != -1 ) {
						poParent->tableinfo.tab23open0 = (real) poParent->CloseUpdate( 23, 0, poParent->tableinfo.tab23open0, poParent->tableinfo.tab23close0, poParent->tableinfo.tab23dim_filter, poParent->tableinfo.tab23cell, poParent->tableinfo.tab23total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab23open1 = (real) poParent->CloseUpdate( 23, 1, poParent->tableinfo.tab23open1, poParent->tableinfo.tab23close1, poParent->tableinfo.tab23dim_filter, poParent->tableinfo.tab23cell, poParent->tableinfo.tab23total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab23events = poParent->__events;
					poParent->tableinfo.tab23close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab23close1 = (real) poDerivedStates->Get_duration_is_immunized_TRUE();
					poParent->tableinfo.tab23cell = (index) poParent->TableCellIndex( 23, 6,
						poParent->sex,
						poParent->ethnicity_short,
						poParent->got_prenat_care,
						poParent->region_nat,
						poParent->educ_mother,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab23total_cell = (index) poParent->tableinfo.tab23cell;
					}
					else {
						poParent->tableinfo.tab23total_cell = (index) poParent->TableTotalCellIndex( 23, 6,
						poParent->sex,
						poParent->ethnicity_short,
						poParent->got_prenat_care,
						poParent->region_nat,
						poParent->educ_mother,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab23dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabImmunizationChildren = value;
}

logical PersonStateFunctions::Set_table_filter_tabMigrationOriginDestination( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_tabMigrationOriginDestination != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_tabMigrationOriginDestination;
		poParent->table_filter_tabMigrationOriginDestination = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[19] >= 0 ) {
				if ( poParent->tableinfo.tab19events == -1 || poParent->tableinfo.tab19events < poParent->__events ) {
					if ( poParent->tableinfo.tab19events != -1 ) {
						poParent->tableinfo.tab19open0 = (real) poParent->CloseUpdate( 19, 0, poParent->tableinfo.tab19open0, poParent->tableinfo.tab19close0, poParent->tableinfo.tab19dim_filter, poParent->tableinfo.tab19cell, poParent->tableinfo.tab19total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab19events = poParent->__events;
					poParent->tableinfo.tab19close0 = (real) poParent->is_alive_value_in;
					poParent->tableinfo.tab19cell = (index) poParent->TableCellIndex( 19, 3,
						poParent->selected_year10,
						poParent->geo_prev,
						poParent->geo );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab19total_cell = (index) poParent->tableinfo.tab19cell;
					}
					else {
						poParent->tableinfo.tab19total_cell = (index) poParent->TableTotalCellIndex( 19, 3,
						poParent->selected_year10,
						poParent->geo_prev,
						poParent->geo );
					}
					poParent->tableinfo.tab19dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_tabMigrationOriginDestination = value;
}

logical PersonStateFunctions::Set_table_filter_TabPopProvAgeEducSex( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabPopProvAgeEducSex != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabPopProvAgeEducSex;
		poParent->table_filter_TabPopProvAgeEducSex = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[5] >= 0 ) {
				if ( poParent->tableinfo.tab5events == -1 || poParent->tableinfo.tab5events < poParent->__events ) {
					if ( poParent->tableinfo.tab5events != -1 ) {
						poParent->tableinfo.tab5open0 = (real) poParent->CloseUpdate( 5, 0, poParent->tableinfo.tab5open0, poParent->tableinfo.tab5close0, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open1 = (real) poParent->CloseUpdate( 5, 1, poParent->tableinfo.tab5open1, poParent->tableinfo.tab5close1, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open2 = (real) poParent->CloseUpdate( 5, 2, poParent->tableinfo.tab5open2, poParent->tableinfo.tab5close2, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open3 = (real) poParent->CloseUpdate( 5, 3, poParent->tableinfo.tab5open3, poParent->tableinfo.tab5close3, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab5open4 = (real) poParent->CloseUpdate( 5, 4, poParent->tableinfo.tab5open4, poParent->tableinfo.tab5close4, poParent->tableinfo.tab5dim_filter, poParent->tableinfo.tab5cell, poParent->tableinfo.tab5total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab5events = poParent->__events;
					poParent->tableinfo.tab5close0 = (real) poDerivedStates->Get_duration_is_aded0_TRUE();
					poParent->tableinfo.tab5close1 = (real) poDerivedStates->Get_duration_is_aded1_TRUE();
					poParent->tableinfo.tab5close2 = (real) poDerivedStates->Get_duration_is_aded2_TRUE();
					poParent->tableinfo.tab5close3 = (real) poDerivedStates->Get_duration_is_child_TRUE();
					poParent->tableinfo.tab5close4 = (real) poDerivedStates->Get_duration_is_old_TRUE();
					poParent->tableinfo.tab5cell = (index) poParent->TableCellIndex( 5, 3,
						poParent->geo_nat,
						poParent->sex,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab5total_cell = (index) poParent->tableinfo.tab5cell;
					}
					else {
						poParent->tableinfo.tab5total_cell = (index) poParent->TableTotalCellIndex( 5, 3,
						poParent->geo_nat,
						poParent->sex,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab5dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabPopProvAgeEducSex = value;
}

logical PersonStateFunctions::Set_table_filter_tabPopulationYearPlace( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_tabPopulationYearPlace != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_tabPopulationYearPlace;
		poParent->table_filter_tabPopulationYearPlace = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[17] >= 0 ) {
				if ( poParent->tableinfo.tab17events == -1 || poParent->tableinfo.tab17events < poParent->__events ) {
					if ( poParent->tableinfo.tab17events != -1 ) {
						poParent->tableinfo.tab17open0 = (real) poParent->CloseUpdate( 17, 0, poParent->tableinfo.tab17open0, poParent->tableinfo.tab17close0, poParent->tableinfo.tab17dim_filter, poParent->tableinfo.tab17cell, poParent->tableinfo.tab17total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab17events = poParent->__events;
					poParent->tableinfo.tab17close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab17cell = (index) poParent->TableCellIndex( 17, 2,
						poParent->geo,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab17total_cell = (index) poParent->tableinfo.tab17cell;
					}
					else {
						poParent->tableinfo.tab17total_cell = (index) poParent->TableTotalCellIndex( 17, 2,
						poParent->geo,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab17dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_tabPopulationYearPlace = value;
}

logical PersonStateFunctions::Set_table_filter_TabPrenatCare( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabPrenatCare != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabPrenatCare;
		poParent->table_filter_TabPrenatCare = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[22] >= 0 ) {
				if ( poParent->tableinfo.tab22events == -1 || poParent->tableinfo.tab22events < poParent->__events ) {
					if ( poParent->tableinfo.tab22events != -1 ) {
						poParent->tableinfo.tab22open0 = (real) poParent->CloseUpdate( 22, 0, poParent->tableinfo.tab22open0, poParent->tableinfo.tab22close0, poParent->tableinfo.tab22dim_filter, poParent->tableinfo.tab22cell, poParent->tableinfo.tab22total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab22open1 = (real) poParent->CloseUpdate( 22, 1, poParent->tableinfo.tab22open1, poParent->tableinfo.tab22close1, poParent->tableinfo.tab22dim_filter, poParent->tableinfo.tab22cell, poParent->tableinfo.tab22total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab22events = poParent->__events;
					poParent->tableinfo.tab22close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab22close1 = (real) poDerivedStates->Get_duration_got_prenat_care_GPC_YES();
					poParent->tableinfo.tab22cell = (index) poParent->TableCellIndex( 22, 3,
						poParent->region_nat,
						poParent->educ_mother,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab22total_cell = (index) poParent->tableinfo.tab22cell;
					}
					else {
						poParent->tableinfo.tab22total_cell = (index) poParent->TableTotalCellIndex( 22, 3,
						poParent->region_nat,
						poParent->educ_mother,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab22dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabPrenatCare = value;
}

logical PersonStateFunctions::Set_table_filter_tabPreSchool( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_tabPreSchool != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_tabPreSchool;
		poParent->table_filter_tabPreSchool = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[0] >= 0 ) {
				if ( poParent->tableinfo.tab0events == -1 || poParent->tableinfo.tab0events < poParent->__events ) {
					if ( poParent->tableinfo.tab0events != -1 ) {
						poParent->tableinfo.tab0open0 = (real) poParent->CloseUpdate( 0, 0, poParent->tableinfo.tab0open0, poParent->tableinfo.tab0close0, poParent->tableinfo.tab0dim_filter, poParent->tableinfo.tab0cell, poParent->tableinfo.tab0total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab0open1 = (real) poParent->CloseUpdate( 0, 1, poParent->tableinfo.tab0open1, poParent->tableinfo.tab0close1, poParent->tableinfo.tab0dim_filter, poParent->tableinfo.tab0cell, poParent->tableinfo.tab0total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab0open2 = (real) poParent->CloseUpdate( 0, 2, poParent->tableinfo.tab0open2, poParent->tableinfo.tab0close2, poParent->tableinfo.tab0dim_filter, poParent->tableinfo.tab0cell, poParent->tableinfo.tab0total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab0open3 = (real) poParent->CloseUpdate( 0, 3, poParent->tableinfo.tab0open3, poParent->tableinfo.tab0close3, poParent->tableinfo.tab0dim_filter, poParent->tableinfo.tab0cell, poParent->tableinfo.tab0total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab0events = poParent->__events;
					poParent->tableinfo.tab0close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab0close1 = (real) poDerivedStates->Get_duration_years_preschool_0();
					poParent->tableinfo.tab0close2 = (real) poDerivedStates->Get_duration_years_preschool_1();
					poParent->tableinfo.tab0close3 = (real) poDerivedStates->Get_duration_years_preschool_2();
					poParent->tableinfo.tab0cell = (index) poParent->TableCellIndex( 0, 2,
						poParent->geo,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab0total_cell = (index) poParent->tableinfo.tab0cell;
					}
					else {
						poParent->tableinfo.tab0total_cell = (index) poParent->TableTotalCellIndex( 0, 2,
						poParent->geo,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab0dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_tabPreSchool = value;
}

logical PersonStateFunctions::Set_table_filter_TabPrimarySchoolPlanning( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabPrimarySchoolPlanning != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabPrimarySchoolPlanning;
		poParent->table_filter_TabPrimarySchoolPlanning = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[4] >= 0 ) {
				if ( poParent->tableinfo.tab4events == -1 || poParent->tableinfo.tab4events < poParent->__events ) {
					if ( poParent->tableinfo.tab4events != -1 ) {
						poParent->tableinfo.tab4open0 = (real) poParent->CloseUpdate( 4, 0, poParent->tableinfo.tab4open0, poParent->tableinfo.tab4close0, poParent->tableinfo.tab4dim_filter, poParent->tableinfo.tab4cell, poParent->tableinfo.tab4total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab4open1 = (real) poParent->CloseUpdate( 4, 1, poParent->tableinfo.tab4open1, poParent->tableinfo.tab4close1, poParent->tableinfo.tab4dim_filter, poParent->tableinfo.tab4cell, poParent->tableinfo.tab4total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						Set_table_unit_TabPrimarySchoolPlanning(poParent->table_unit_TabPrimarySchoolPlanning+1);
						poParent->tableinfo.tab4close2 = (real) poParent->table_unit_TabPrimarySchoolPlanning;
						poParent->tableinfo.tab4open2 = (real) poParent->CloseUpdate( 4, 2, poParent->tableinfo.tab4open2, poParent->tableinfo.tab4close2, poParent->tableinfo.tab4dim_filter, poParent->tableinfo.tab4cell, poParent->tableinfo.tab4total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab4events = poParent->__events;
					poParent->tableinfo.tab4close0 = (real) poParent->educ_one_required_rooms_value_out;
					poParent->tableinfo.tab4close1 = (real) poParent->educ_one_required_teachers_value_out;
					poParent->tableinfo.tab4close2 = (real) poParent->table_unit_TabPrimarySchoolPlanning;
					poParent->tableinfo.tab4cell = (index) poParent->TableCellIndex( 4, 3,
						poParent->geo_nat,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						poParent->educ_one_grade - r_min_EDUC_ONE_GRADE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab4total_cell = (index) poParent->tableinfo.tab4cell;
					}
					else {
						poParent->tableinfo.tab4total_cell = (index) poParent->TableTotalCellIndex( 4, 3,
						poParent->geo_nat,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						poParent->educ_one_grade - r_min_EDUC_ONE_GRADE );
					}
					poParent->tableinfo.tab4dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabPrimarySchoolPlanning = value;
}

logical PersonStateFunctions::Set_table_filter_TabPrimSchoolEntries( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabPrimSchoolEntries != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabPrimSchoolEntries;
		poParent->table_filter_TabPrimSchoolEntries = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[9] >= 0 ) {
				if ( poParent->tableinfo.tab9events == -1 || poParent->tableinfo.tab9events < poParent->__events ) {
					if ( poParent->tableinfo.tab9events != -1 ) {
						poParent->tableinfo.tab9open0 = (real) poParent->CloseUpdate( 9, 0, poParent->tableinfo.tab9open0, poParent->tableinfo.tab9close0, poParent->tableinfo.tab9dim_filter, poParent->tableinfo.tab9cell, poParent->tableinfo.tab9total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab9events = poParent->__events;
					poParent->tableinfo.tab9close0 = (real) poDerivedStates->entrances_educ_one_grade_attended_1;
					poParent->tableinfo.tab9cell = (index) poParent->TableCellIndex( 9, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						poParent->geo_nat );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab9total_cell = (index) poParent->tableinfo.tab9cell;
					}
					else {
						poParent->tableinfo.tab9total_cell = (index) poParent->TableTotalCellIndex( 9, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						poParent->geo_nat );
					}
					poParent->tableinfo.tab9dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabPrimSchoolEntries = value;
}

logical PersonStateFunctions::Set_table_filter_TabPrimSchoolGraduations( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabPrimSchoolGraduations != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabPrimSchoolGraduations;
		poParent->table_filter_TabPrimSchoolGraduations = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[10] >= 0 ) {
				if ( poParent->tableinfo.tab10events == -1 || poParent->tableinfo.tab10events < poParent->__events ) {
					if ( poParent->tableinfo.tab10events != -1 ) {
						poParent->tableinfo.tab10open0 = (real) poParent->CloseUpdate( 10, 0, poParent->tableinfo.tab10open0, poParent->tableinfo.tab10close0, poParent->tableinfo.tab10dim_filter, poParent->tableinfo.tab10cell, poParent->tableinfo.tab10total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab10events = poParent->__events;
					poParent->tableinfo.tab10close0 = (real) poDerivedStates->entrances_has_passed_primary_TRUE;
					poParent->tableinfo.tab10cell = (index) poParent->TableCellIndex( 10, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						poParent->geo_nat );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab10total_cell = (index) poParent->tableinfo.tab10cell;
					}
					else {
						poParent->tableinfo.tab10total_cell = (index) poParent->TableTotalCellIndex( 10, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						poParent->geo_nat );
					}
					poParent->tableinfo.tab10dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabPrimSchoolGraduations = value;
}

logical PersonStateFunctions::Set_table_filter_TabPrimSchoolOutOfSchool9to11( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabPrimSchoolOutOfSchool9to11 != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabPrimSchoolOutOfSchool9to11;
		poParent->table_filter_TabPrimSchoolOutOfSchool9to11 = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[11] >= 0 ) {
				if ( poParent->tableinfo.tab11events == -1 || poParent->tableinfo.tab11events < poParent->__events ) {
					if ( poParent->tableinfo.tab11events != -1 ) {
						poParent->tableinfo.tab11open0 = (real) poParent->CloseUpdate( 11, 0, poParent->tableinfo.tab11open0, poParent->tableinfo.tab11close0, poParent->tableinfo.tab11dim_filter, poParent->tableinfo.tab11cell, poParent->tableinfo.tab11total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab11open1 = (real) poParent->CloseUpdate( 11, 1, poParent->tableinfo.tab11open1, poParent->tableinfo.tab11close1, poParent->tableinfo.tab11dim_filter, poParent->tableinfo.tab11cell, poParent->tableinfo.tab11total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab11events = poParent->__events;
					poParent->tableinfo.tab11close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab11close1 = (real) poDerivedStates->Get_duration_tab_out_of_primary_TRUE();
					poParent->tableinfo.tab11cell = (index) poParent->TableCellIndex( 11, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						poParent->geo_nat );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab11total_cell = (index) poParent->tableinfo.tab11cell;
					}
					else {
						poParent->tableinfo.tab11total_cell = (index) poParent->TableTotalCellIndex( 11, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						poParent->geo_nat );
					}
					poParent->tableinfo.tab11dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabPrimSchoolOutOfSchool9to11 = value;
}

logical PersonStateFunctions::Set_table_filter_TabSchool2AttainmentsTab( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabSchool2AttainmentsTab != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabSchool2AttainmentsTab;
		poParent->table_filter_TabSchool2AttainmentsTab = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[13] >= 0 ) {
				if ( poParent->tableinfo.tab13events == -1 || poParent->tableinfo.tab13events < poParent->__events ) {
					if ( poParent->tableinfo.tab13events != -1 ) {
						poParent->tableinfo.tab13open0 = (real) poParent->CloseUpdate( 13, 0, poParent->tableinfo.tab13open0, poParent->tableinfo.tab13close0, poParent->tableinfo.tab13dim_filter, poParent->tableinfo.tab13cell, poParent->tableinfo.tab13total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab13open1 = (real) poParent->CloseUpdate( 13, 1, poParent->tableinfo.tab13open1, poParent->tableinfo.tab13close1, poParent->tableinfo.tab13dim_filter, poParent->tableinfo.tab13cell, poParent->tableinfo.tab13total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab13open2 = (real) poParent->CloseUpdate( 13, 2, poParent->tableinfo.tab13open2, poParent->tableinfo.tab13close2, poParent->tableinfo.tab13dim_filter, poParent->tableinfo.tab13cell, poParent->tableinfo.tab13total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab13open3 = (real) poParent->CloseUpdate( 13, 3, poParent->tableinfo.tab13open3, poParent->tableinfo.tab13close3, poParent->tableinfo.tab13dim_filter, poParent->tableinfo.tab13cell, poParent->tableinfo.tab13total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab13events = poParent->__events;
					poParent->tableinfo.tab13close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab13close1 = (real) poDerivedStates->Get_duration_educ_two_level_ETL_GRADUATE();
					poParent->tableinfo.tab13close2 = (real) poDerivedStates->Get_duration_educ_two_level_ETL_INCOMPLETE();
					poParent->tableinfo.tab13close3 = (real) poDerivedStates->Get_duration_educ_two_level_ETL_NEVER();
					poParent->tableinfo.tab13cell = (index) poParent->TableCellIndex( 13, 2,
						poParent->sex,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab13total_cell = (index) poParent->tableinfo.tab13cell;
					}
					else {
						poParent->tableinfo.tab13total_cell = (index) poParent->TableTotalCellIndex( 13, 2,
						poParent->sex,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab13dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabSchool2AttainmentsTab = value;
}

logical PersonStateFunctions::Set_table_filter_TabSchool2TrackTab( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_TabSchool2TrackTab != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_TabSchool2TrackTab;
		poParent->table_filter_TabSchool2TrackTab = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[12] >= 0 ) {
				if ( poParent->tableinfo.tab12events == -1 || poParent->tableinfo.tab12events < poParent->__events ) {
					if ( poParent->tableinfo.tab12events != -1 ) {
						Set_table_unit_TabSchool2TrackTab(poParent->table_unit_TabSchool2TrackTab+1);
						poParent->tableinfo.tab12close0 = (real) poParent->table_unit_TabSchool2TrackTab;
						poParent->tableinfo.tab12open0 = (real) poParent->CloseUpdate( 12, 0, poParent->tableinfo.tab12open0, poParent->tableinfo.tab12close0, poParent->tableinfo.tab12dim_filter, poParent->tableinfo.tab12cell, poParent->tableinfo.tab12total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab12events = poParent->__events;
					poParent->tableinfo.tab12close0 = (real) poParent->table_unit_TabSchool2TrackTab;
					poParent->tableinfo.tab12cell = (index) poParent->TableCellIndex( 12, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						poParent->educ_two_grade_attended - r_min_EDUC_TWO_GRADE0 );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab12total_cell = (index) poParent->tableinfo.tab12cell;
					}
					else {
						poParent->tableinfo.tab12total_cell = (index) poParent->TableTotalCellIndex( 12, 2,
						poParent->sim_year - r_min_SIM_YEAR_RANGE,
						poParent->educ_two_grade_attended - r_min_EDUC_TWO_GRADE0 );
					}
					poParent->tableinfo.tab12dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_TabSchool2TrackTab = value;
}

logical PersonStateFunctions::Set_table_filter_tabStuntingSexRegMotherYob( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_tabStuntingSexRegMotherYob != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_tabStuntingSexRegMotherYob;
		poParent->table_filter_tabStuntingSexRegMotherYob = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[21] >= 0 ) {
				if ( poParent->tableinfo.tab21events == -1 || poParent->tableinfo.tab21events < poParent->__events ) {
					if ( poParent->tableinfo.tab21events != -1 ) {
						poParent->tableinfo.tab21open0 = (real) poParent->CloseUpdate( 21, 0, poParent->tableinfo.tab21open0, poParent->tableinfo.tab21close0, poParent->tableinfo.tab21dim_filter, poParent->tableinfo.tab21cell, poParent->tableinfo.tab21total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						Set_table_unit_tabStuntingSexRegMotherYob(poParent->table_unit_tabStuntingSexRegMotherYob+1);
						poParent->tableinfo.tab21close1 = (real) poParent->table_unit_tabStuntingSexRegMotherYob;
						poParent->tableinfo.tab21open1 = (real) poParent->CloseUpdate( 21, 1, poParent->tableinfo.tab21open1, poParent->tableinfo.tab21close1, poParent->tableinfo.tab21dim_filter, poParent->tableinfo.tab21cell, poParent->tableinfo.tab21total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab21events = poParent->__events;
					poParent->tableinfo.tab21close0 = (real) poParent->is_stunted_value_in;
					poParent->tableinfo.tab21close1 = (real) poParent->table_unit_tabStuntingSexRegMotherYob;
					poParent->tableinfo.tab21cell = (index) poParent->TableCellIndex( 21, 4,
						poParent->sex,
						poParent->educ_mother,
						poParent->region_nat,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab21total_cell = (index) poParent->tableinfo.tab21cell;
					}
					else {
						poParent->tableinfo.tab21total_cell = (index) poParent->TableTotalCellIndex( 21, 4,
						poParent->sex,
						poParent->educ_mother,
						poParent->region_nat,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab21dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_tabStuntingSexRegMotherYob = value;
}

logical PersonStateFunctions::Set_table_filter_tabStuntingSexRegYob( logical value )
{

	value = value != 0;
	if ( poParent->table_filter_tabStuntingSexRegYob != value )
	{


		logical	OldValue;
		OldValue = poParent->table_filter_tabStuntingSexRegYob;
		poParent->table_filter_tabStuntingSexRegYob = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[20] >= 0 ) {
				if ( poParent->tableinfo.tab20events == -1 || poParent->tableinfo.tab20events < poParent->__events ) {
					if ( poParent->tableinfo.tab20events != -1 ) {
						poParent->tableinfo.tab20open0 = (real) poParent->CloseUpdate( 20, 0, poParent->tableinfo.tab20open0, poParent->tableinfo.tab20close0, poParent->tableinfo.tab20dim_filter, poParent->tableinfo.tab20cell, poParent->tableinfo.tab20total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						Set_table_unit_tabStuntingSexRegYob(poParent->table_unit_tabStuntingSexRegYob+1);
						poParent->tableinfo.tab20close1 = (real) poParent->table_unit_tabStuntingSexRegYob;
						poParent->tableinfo.tab20open1 = (real) poParent->CloseUpdate( 20, 1, poParent->tableinfo.tab20open1, poParent->tableinfo.tab20close1, poParent->tableinfo.tab20dim_filter, poParent->tableinfo.tab20cell, poParent->tableinfo.tab20total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab20events = poParent->__events;
					poParent->tableinfo.tab20close0 = (real) poParent->is_stunted_value_in;
					poParent->tableinfo.tab20close1 = (real) poParent->table_unit_tabStuntingSexRegYob;
					poParent->tableinfo.tab20cell = (index) poParent->TableCellIndex( 20, 2,
						poParent->region_nat,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab20total_cell = (index) poParent->tableinfo.tab20cell;
					}
					else {
						poParent->tableinfo.tab20total_cell = (index) poParent->TableTotalCellIndex( 20, 2,
						poParent->region_nat,
						poParent->sim_year - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab20dim_filter = OldValue;
				}
			}
		}
	}
	return poParent->table_filter_tabStuntingSexRegYob = value;
}

int PersonStateFunctions::Set_time_in_parity( int value )
{

	if ( poParent->time_in_parity != value )
	{


		int	OldValue;
		OldValue = poParent->time_in_parity;
		poParent->time_in_parity = value;

		if ( poParent->poFertilityDetailedEvent->cRecalc == 0 ) {
			poParent->poFertilityDetailedEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFertilityDetailedEvent );
		};
	}
	return poParent->time_in_parity = value;
}

UNION1_GROUP PersonStateFunctions::Set_union1_group( UNION1_GROUP value )
{

	CHKLMT(value, union1_group, 0, rc_size_UNION1_GROUP - 1);
	if ( value < 0 ) {
		value = (UNION1_GROUP) 0;
	}
	else if ( value >= rc_size_UNION1_GROUP ) {
		value = (UNION1_GROUP) ( rc_size_UNION1_GROUP - 1 );
	}
	if ( poParent->union1_group != value )
	{


		UNION1_GROUP	OldValue;
		OldValue = poParent->union1_group;
		poParent->union1_group = value;

		// actorset dimension or order state
		if (poParent->actorset_filter_asWomenNoChildren)
		{
			asWomenNoChildren[OldValue][RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][poParent->has_spouse]->Remove(poParent->poasWomenNoChildren);
			poParent->poasWomenNoChildren = asWomenNoChildren[poParent->union1_group][RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][poParent->has_spouse]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asWomenWithChildren)
		{
			asWomenWithChildren[OldValue][poParent->child_agegr][poParent->moth_agegr][poParent->has_spouse]->Remove(poParent->poasWomenWithChildren);
			poParent->poasWomenWithChildren = asWomenWithChildren[poParent->union1_group][poParent->child_agegr][poParent->moth_agegr][poParent->has_spouse]->NewNode(poParent);
		}
		if ( poParent->poFirstUnionFormationEvent->cRecalc == 0 ) {
			poParent->poFirstUnionFormationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstUnionFormationEvent );
		};
	}
	return poParent->union1_group = value;
}

ALL_YEAR_RANGE PersonStateFunctions::Set_year_of_birth( ALL_YEAR_RANGE value )
{

	CHKLMT(value, year_of_birth, r_min_ALL_YEAR_RANGE, r_max_ALL_YEAR_RANGE);
	if ( value < r_min_ALL_YEAR_RANGE ) {
		value = r_min_ALL_YEAR_RANGE;
	}
	else if ( value > r_max_ALL_YEAR_RANGE ) {
		value = r_max_ALL_YEAR_RANGE;
	}
	if ( poParent->year_of_birth != value )
	{
#line 17 "../code/TablesHumanCapital.mpp"
		Set_table_filter_tabHCI( poParent->in_projected_time&& poParent->in_hci_sample&& value>= MIN(SIM_YEAR_RANGE) );
#line 12380 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 45 "../code/TablesHumanCapital.mpp"
		Set_table_filter_TabHCIDistrict( poParent->in_projected_time&& poParent->in_hci_sample&& value>= MIN(SIM_YEAR_RANGE) );
#line 12383 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 76 "../code/TablesEducation.mpp"
		Set_table_filter_tabPreSchool( WITHIN(SIM_YEAR_RANGE,value)&& poParent->integer_age== 10 );
#line 12386 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"


		ALL_YEAR_RANGE	OldValue;
		OldValue = poParent->year_of_birth;
		poParent->year_of_birth = value;

		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[1] >= 0 && ( poParent->tableinfo.tab1dim_filter || poParent->table_filter_tabEducationFateGeobirYob ) ) {
				if ( poParent->tableinfo.tab1events == -1 || poParent->tableinfo.tab1events < poParent->__events ) {
					if ( poParent->tableinfo.tab1events != -1 ) {
						poParent->tableinfo.tab1open0 = (real) poParent->CloseUpdate( 1, 0, poParent->tableinfo.tab1open0, poParent->tableinfo.tab1close0, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab1open1 = (real) poParent->CloseUpdate( 1, 1, poParent->tableinfo.tab1open1, poParent->tableinfo.tab1close1, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab1open2 = (real) poParent->CloseUpdate( 1, 2, poParent->tableinfo.tab1open2, poParent->tableinfo.tab1close2, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab1open3 = (real) poParent->CloseUpdate( 1, 3, poParent->tableinfo.tab1open3, poParent->tableinfo.tab1close3, poParent->tableinfo.tab1dim_filter, poParent->tableinfo.tab1cell, poParent->tableinfo.tab1total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab1events = poParent->__events;
					poParent->tableinfo.tab1close0 = (real) poDerivedStates->Get_duration_();
					poParent->tableinfo.tab1close1 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_HIGH();
					poParent->tableinfo.tab1close2 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_LOW();
					poParent->tableinfo.tab1close3 = (real) poDerivedStates->Get_duration_educ_one_fate_EOL_MEDIUM();
					poParent->tableinfo.tab1cell = (index) poParent->TableCellIndex( 1, 3,
						poParent->sex,
						poParent->geo_birth,
						OldValue - r_min_ALL_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab1total_cell = (index) poParent->tableinfo.tab1cell;
					}
					else {
						poParent->tableinfo.tab1total_cell = (index) poParent->TableTotalCellIndex( 1, 3,
						poParent->sex,
						poParent->geo_birth,
						OldValue - r_min_ALL_YEAR_RANGE );
					}
					poParent->tableinfo.tab1dim_filter = poParent->table_filter_tabEducationFateGeobirYob;
				}
			}
			if ( gplTableIndexes[6] >= 0 && ( poParent->tableinfo.tab6dim_filter || poParent->table_filter_TabEducFateDistrYob ) ) {
				if ( poParent->tableinfo.tab6events == -1 || poParent->tableinfo.tab6events < poParent->__events ) {
					if ( poParent->tableinfo.tab6events != -1 ) {
						Set_table_unit_TabEducFateDistrYob(poParent->table_unit_TabEducFateDistrYob+1);
						poParent->tableinfo.tab6close0 = (real) poParent->table_unit_TabEducFateDistrYob;
						poParent->tableinfo.tab6open0 = (real) poParent->CloseUpdate( 6, 0, poParent->tableinfo.tab6open0, poParent->tableinfo.tab6close0, poParent->tableinfo.tab6dim_filter, poParent->tableinfo.tab6cell, poParent->tableinfo.tab6total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab6events = poParent->__events;
					poParent->tableinfo.tab6close0 = (real) poParent->table_unit_TabEducFateDistrYob;
					poParent->tableinfo.tab6cell = (index) poParent->TableCellIndex( 6, 3,
						poDerivedStates->aggregate_region_birth_REGION_NAT,
						OldValue - r_min_ALL_YEAR_RANGE,
						poParent->educ_one_fate );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab6total_cell = (index) poParent->tableinfo.tab6cell;
					}
					else {
						poParent->tableinfo.tab6total_cell = (index) poParent->TableTotalCellIndex( 6, 3,
						poDerivedStates->aggregate_region_birth_REGION_NAT,
						OldValue - r_min_ALL_YEAR_RANGE,
						poParent->educ_one_fate );
					}
					poParent->tableinfo.tab6dim_filter = poParent->table_filter_TabEducFateDistrYob;
				}
			}
		}
#line 19 "../code/TablesEducation.mpp"
		Set_tab_sim_yob( COERCE(SIM_YEAR_RANGE,poParent->year_of_birth) );
#line 12451 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->year_of_birth = value;
}

SIM_YEAR_RANGE PersonStateFunctions::Set_year_of_first_immigration( SIM_YEAR_RANGE value )
{

	CHKLMT(value, year_of_first_immigration, r_min_SIM_YEAR_RANGE, r_max_SIM_YEAR_RANGE);
	if ( value < r_min_SIM_YEAR_RANGE ) {
		value = r_min_SIM_YEAR_RANGE;
	}
	else if ( value > r_max_SIM_YEAR_RANGE ) {
		value = r_max_SIM_YEAR_RANGE;
	}
	if ( poParent->year_of_first_immigration != value )
	{


		SIM_YEAR_RANGE	OldValue;
		OldValue = poParent->year_of_first_immigration;
		poParent->year_of_first_immigration = value;

		// actorset dimension or order state
		if (poParent->actorset_filter_asAvailableMaleByTypeImmiScratch)
		{
			asAvailableMaleByTypeImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->spouse_group][poParent->geo][RANGE_POS(SIM_YEAR_RANGE,OldValue)]->Remove(poParent->poasAvailableMaleByTypeImmiScratch);
			poParent->poasAvailableMaleByTypeImmiScratch = asAvailableMaleByTypeImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->spouse_group][poParent->geo][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asAvailableMaleImmiScratch)
		{
			asAvailableMaleImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->geo][RANGE_POS(SIM_YEAR_RANGE,OldValue)]->Remove(poParent->poasAvailableMaleImmiScratch);
			poParent->poasAvailableMaleImmiScratch = asAvailableMaleImmiScratch[RANGE_POS(MALE_SPOUSE_AGE, poParent->male_age)][poParent->geo][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asFemaleInUnionByAgeAndPartnerAgeImmiScratch)
		{
			asFemaleInUnionByAgeAndPartnerAgeImmiScratch[RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][RANGE_POS(MALE_SPOUSE_AGE, poParent->male_spouse_age)][RANGE_POS(SIM_YEAR_RANGE,OldValue)]->Remove(poParent->poasFemaleInUnionByAgeAndPartnerAgeImmiScratch);
			poParent->poasFemaleInUnionByAgeAndPartnerAgeImmiScratch = asFemaleInUnionByAgeAndPartnerAgeImmiScratch[RANGE_POS(FEMALE_SPOUSE_AGE, poParent->female_spouse_age)][RANGE_POS(MALE_SPOUSE_AGE, poParent->male_spouse_age)][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)]->NewNode(poParent);
		}
		if (poParent->actorset_filter_asPotentialImmigrantMothers)
		{
			asPotentialImmigrantMothers[RANGE_POS(AGE_RANGE, poParent->integer_age)][RANGE_POS(SIM_YEAR_RANGE,OldValue)][poParent->geo]->Remove(poParent->poasPotentialImmigrantMothers);
			poParent->poasPotentialImmigrantMothers = asPotentialImmigrantMothers[RANGE_POS(AGE_RANGE, poParent->integer_age)][RANGE_POS(SIM_YEAR_RANGE, poParent->year_of_first_immigration)][poParent->geo]->NewNode(poParent);
		}
	}
	return poParent->year_of_first_immigration = value;
}

counter PersonStateFunctions::Set_count_mlHHFatherChildren( counter value )
{

	if ( poParent->count_mlHHFatherChildren != value )
	{


		counter	OldValue;
		OldValue = poParent->count_mlHHFatherChildren;
		poParent->count_mlHHFatherChildren = value;

#line 92 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_lives_with_dependent_child( ( ( poParent->sex== FEMALE&& poParent->count_mlHHMotherChildren> 0 )|| ( poParent->sex== MALE&& poParent->count_mlHHFatherChildren> 0 ) )? TRUE: FALSE );
#line 12512 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->count_mlHHFatherChildren = value;
}

counter PersonStateFunctions::Set_count_mlHHMotherChildren( counter value )
{

	if ( poParent->count_mlHHMotherChildren != value )
	{


		counter	OldValue;
		OldValue = poParent->count_mlHHMotherChildren;
		poParent->count_mlHHMotherChildren = value;

#line 85 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_age_youngest_child( ( poParent->sex== FEMALE&& poParent->count_mlHHMotherChildren> 0 )? double(poParent->min_over_mlHHMotherChildren_integer_age): TIME_INFINITE );
#line 12530 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 92 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_lives_with_dependent_child( ( ( poParent->sex== FEMALE&& poParent->count_mlHHMotherChildren> 0 )|| ( poParent->sex== MALE&& poParent->count_mlHHFatherChildren> 0 ) )? TRUE: FALSE );
#line 12533 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->count_mlHHMotherChildren = value;
}

double PersonStateFunctions::Set_lCalibrator_mort_female_0( double value )
{

	if ( poParent->lCalibrator_mort_female_0 != value )
	{


		double	OldValue;
		OldValue = poParent->lCalibrator_mort_female_0;
		poParent->lCalibrator_mort_female_0 = value;

#line 67 "../code/MortalityDetailed.mpp"
		Set_child_mortality( ( poParent->integer_age== 0&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_0): ( poParent->integer_age== 1&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_1): ( poParent->integer_age== 2&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_2): ( poParent->integer_age== 3&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_3): ( poParent->integer_age== 4&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_4): ( poParent->integer_age== 0&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_0): ( poParent->integer_age== 1&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_1): ( poParent->integer_age== 2&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_2): ( poParent->integer_age== 3&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_3): ( poParent->integer_age== 4&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_4): 0.0 );
#line 12551 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lCalibrator_mort_female_0 = value;
}

double PersonStateFunctions::Set_lCalibrator_mort_female_1( double value )
{

	if ( poParent->lCalibrator_mort_female_1 != value )
	{


		double	OldValue;
		OldValue = poParent->lCalibrator_mort_female_1;
		poParent->lCalibrator_mort_female_1 = value;

#line 67 "../code/MortalityDetailed.mpp"
		Set_child_mortality( ( poParent->integer_age== 0&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_0): ( poParent->integer_age== 1&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_1): ( poParent->integer_age== 2&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_2): ( poParent->integer_age== 3&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_3): ( poParent->integer_age== 4&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_4): ( poParent->integer_age== 0&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_0): ( poParent->integer_age== 1&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_1): ( poParent->integer_age== 2&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_2): ( poParent->integer_age== 3&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_3): ( poParent->integer_age== 4&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_4): 0.0 );
#line 12569 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lCalibrator_mort_female_1 = value;
}

double PersonStateFunctions::Set_lCalibrator_mort_female_2( double value )
{

	if ( poParent->lCalibrator_mort_female_2 != value )
	{


		double	OldValue;
		OldValue = poParent->lCalibrator_mort_female_2;
		poParent->lCalibrator_mort_female_2 = value;

#line 67 "../code/MortalityDetailed.mpp"
		Set_child_mortality( ( poParent->integer_age== 0&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_0): ( poParent->integer_age== 1&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_1): ( poParent->integer_age== 2&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_2): ( poParent->integer_age== 3&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_3): ( poParent->integer_age== 4&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_4): ( poParent->integer_age== 0&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_0): ( poParent->integer_age== 1&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_1): ( poParent->integer_age== 2&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_2): ( poParent->integer_age== 3&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_3): ( poParent->integer_age== 4&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_4): 0.0 );
#line 12587 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lCalibrator_mort_female_2 = value;
}

double PersonStateFunctions::Set_lCalibrator_mort_female_3( double value )
{

	if ( poParent->lCalibrator_mort_female_3 != value )
	{


		double	OldValue;
		OldValue = poParent->lCalibrator_mort_female_3;
		poParent->lCalibrator_mort_female_3 = value;

#line 67 "../code/MortalityDetailed.mpp"
		Set_child_mortality( ( poParent->integer_age== 0&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_0): ( poParent->integer_age== 1&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_1): ( poParent->integer_age== 2&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_2): ( poParent->integer_age== 3&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_3): ( poParent->integer_age== 4&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_4): ( poParent->integer_age== 0&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_0): ( poParent->integer_age== 1&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_1): ( poParent->integer_age== 2&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_2): ( poParent->integer_age== 3&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_3): ( poParent->integer_age== 4&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_4): 0.0 );
#line 12605 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lCalibrator_mort_female_3 = value;
}

double PersonStateFunctions::Set_lCalibrator_mort_female_4( double value )
{

	if ( poParent->lCalibrator_mort_female_4 != value )
	{


		double	OldValue;
		OldValue = poParent->lCalibrator_mort_female_4;
		poParent->lCalibrator_mort_female_4 = value;

#line 67 "../code/MortalityDetailed.mpp"
		Set_child_mortality( ( poParent->integer_age== 0&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_0): ( poParent->integer_age== 1&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_1): ( poParent->integer_age== 2&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_2): ( poParent->integer_age== 3&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_3): ( poParent->integer_age== 4&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_4): ( poParent->integer_age== 0&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_0): ( poParent->integer_age== 1&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_1): ( poParent->integer_age== 2&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_2): ( poParent->integer_age== 3&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_3): ( poParent->integer_age== 4&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_4): 0.0 );
#line 12623 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lCalibrator_mort_female_4 = value;
}

double PersonStateFunctions::Set_lCalibrator_mort_male_0( double value )
{

	if ( poParent->lCalibrator_mort_male_0 != value )
	{


		double	OldValue;
		OldValue = poParent->lCalibrator_mort_male_0;
		poParent->lCalibrator_mort_male_0 = value;

#line 67 "../code/MortalityDetailed.mpp"
		Set_child_mortality( ( poParent->integer_age== 0&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_0): ( poParent->integer_age== 1&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_1): ( poParent->integer_age== 2&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_2): ( poParent->integer_age== 3&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_3): ( poParent->integer_age== 4&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_4): ( poParent->integer_age== 0&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_0): ( poParent->integer_age== 1&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_1): ( poParent->integer_age== 2&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_2): ( poParent->integer_age== 3&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_3): ( poParent->integer_age== 4&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_4): 0.0 );
#line 12641 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lCalibrator_mort_male_0 = value;
}

double PersonStateFunctions::Set_lCalibrator_mort_male_1( double value )
{

	if ( poParent->lCalibrator_mort_male_1 != value )
	{


		double	OldValue;
		OldValue = poParent->lCalibrator_mort_male_1;
		poParent->lCalibrator_mort_male_1 = value;

#line 67 "../code/MortalityDetailed.mpp"
		Set_child_mortality( ( poParent->integer_age== 0&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_0): ( poParent->integer_age== 1&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_1): ( poParent->integer_age== 2&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_2): ( poParent->integer_age== 3&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_3): ( poParent->integer_age== 4&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_4): ( poParent->integer_age== 0&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_0): ( poParent->integer_age== 1&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_1): ( poParent->integer_age== 2&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_2): ( poParent->integer_age== 3&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_3): ( poParent->integer_age== 4&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_4): 0.0 );
#line 12659 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lCalibrator_mort_male_1 = value;
}

double PersonStateFunctions::Set_lCalibrator_mort_male_2( double value )
{

	if ( poParent->lCalibrator_mort_male_2 != value )
	{


		double	OldValue;
		OldValue = poParent->lCalibrator_mort_male_2;
		poParent->lCalibrator_mort_male_2 = value;

#line 67 "../code/MortalityDetailed.mpp"
		Set_child_mortality( ( poParent->integer_age== 0&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_0): ( poParent->integer_age== 1&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_1): ( poParent->integer_age== 2&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_2): ( poParent->integer_age== 3&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_3): ( poParent->integer_age== 4&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_4): ( poParent->integer_age== 0&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_0): ( poParent->integer_age== 1&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_1): ( poParent->integer_age== 2&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_2): ( poParent->integer_age== 3&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_3): ( poParent->integer_age== 4&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_4): 0.0 );
#line 12677 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lCalibrator_mort_male_2 = value;
}

double PersonStateFunctions::Set_lCalibrator_mort_male_3( double value )
{

	if ( poParent->lCalibrator_mort_male_3 != value )
	{


		double	OldValue;
		OldValue = poParent->lCalibrator_mort_male_3;
		poParent->lCalibrator_mort_male_3 = value;

#line 67 "../code/MortalityDetailed.mpp"
		Set_child_mortality( ( poParent->integer_age== 0&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_0): ( poParent->integer_age== 1&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_1): ( poParent->integer_age== 2&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_2): ( poParent->integer_age== 3&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_3): ( poParent->integer_age== 4&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_4): ( poParent->integer_age== 0&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_0): ( poParent->integer_age== 1&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_1): ( poParent->integer_age== 2&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_2): ( poParent->integer_age== 3&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_3): ( poParent->integer_age== 4&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_4): 0.0 );
#line 12695 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lCalibrator_mort_male_3 = value;
}

double PersonStateFunctions::Set_lCalibrator_mort_male_4( double value )
{

	if ( poParent->lCalibrator_mort_male_4 != value )
	{


		double	OldValue;
		OldValue = poParent->lCalibrator_mort_male_4;
		poParent->lCalibrator_mort_male_4 = value;

#line 67 "../code/MortalityDetailed.mpp"
		Set_child_mortality( ( poParent->integer_age== 0&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_0): ( poParent->integer_age== 1&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_1): ( poParent->integer_age== 2&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_2): ( poParent->integer_age== 3&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_3): ( poParent->integer_age== 4&& poParent->sex== MALE )? double(poParent->lCalibrator_mort_male_4): ( poParent->integer_age== 0&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_0): ( poParent->integer_age== 1&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_1): ( poParent->integer_age== 2&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_2): ( poParent->integer_age== 3&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_3): ( poParent->integer_age== 4&& poParent->sex== FEMALE )? double(poParent->lCalibrator_mort_female_4): 0.0 );
#line 12713 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lCalibrator_mort_male_4 = value;
}

EDUC_ONE_LEVEL PersonStateFunctions::Set_lSpouse_educ_one_fate( EDUC_ONE_LEVEL value )
{

	CHKLMT(value, lSpouse_educ_one_fate, 0, rc_size_EDUC_ONE_LEVEL - 1);
	if ( value < 0 ) {
		value = (EDUC_ONE_LEVEL) 0;
	}
	else if ( value >= rc_size_EDUC_ONE_LEVEL ) {
		value = (EDUC_ONE_LEVEL) ( rc_size_EDUC_ONE_LEVEL - 1 );
	}
	if ( poParent->lSpouse_educ_one_fate != value )
	{


		EDUC_ONE_LEVEL	OldValue;
		OldValue = poParent->lSpouse_educ_one_fate;
		poParent->lSpouse_educ_one_fate = value;

#line 87 "../code/FamilyPartnerMatching.mpp"
		Set_spouse_group( ( poParent->has_spouse&& poParent->lSpouse_educ_one_fate== EOL_HIGH )? SMG_02: ( poParent->has_spouse&& poParent->lSpouse_educ_one_fate== EOL_MEDIUM )? SMG_01: SMG_00 );
#line 12738 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lSpouse_educ_one_fate = value;
}

AGE_RANGE PersonStateFunctions::Set_lSpouse_integer_age( AGE_RANGE value )
{

	CHKLMT(value, lSpouse_integer_age, r_min_AGE_RANGE, r_max_AGE_RANGE);
	if ( value < r_min_AGE_RANGE ) {
		value = r_min_AGE_RANGE;
	}
	else if ( value > r_max_AGE_RANGE ) {
		value = r_max_AGE_RANGE;
	}
	if ( poParent->lSpouse_integer_age != value )
	{


		AGE_RANGE	OldValue;
		OldValue = poParent->lSpouse_integer_age;
		poParent->lSpouse_integer_age = value;

#line 75 "../code/FamilyPartnerMatching.mpp"
		Set_male_spouse_age( ( poParent->has_spouse&& poParent->sex== FEMALE )? COERCE(MALE_SPOUSE_AGE,poParent->lSpouse_integer_age): MAX(MALE_SPOUSE_AGE) );
#line 12763 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lSpouse_integer_age = value;
}

AGE_RANGE PersonStateFunctions::Set_min_over_mlHHMotherChildren_integer_age( AGE_RANGE value )
{

	CHKLMT(value, min_over_mlHHMotherChildren_integer_age, r_min_AGE_RANGE, r_max_AGE_RANGE);
	if ( value < r_min_AGE_RANGE ) {
		value = r_min_AGE_RANGE;
	}
	else if ( value > r_max_AGE_RANGE ) {
		value = r_max_AGE_RANGE;
	}
	if ( poParent->min_over_mlHHMotherChildren_integer_age != value )
	{


		AGE_RANGE	OldValue;
		OldValue = poParent->min_over_mlHHMotherChildren_integer_age;
		poParent->min_over_mlHHMotherChildren_integer_age = value;

#line 85 "../code/FamilyFemalePartnershipStatus.mpp"
		Set_age_youngest_child( ( poParent->sex== FEMALE&& poParent->count_mlHHMotherChildren> 0 )? double(poParent->min_over_mlHHMotherChildren_integer_age): TIME_INFINITE );
#line 12788 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->min_over_mlHHMotherChildren_integer_age = value;
}

integer PersonStateFunctions::Set_sum_over_mlHHFatherChildren_is_alive( integer value )
{

	if ( poParent->sum_over_mlHHFatherChildren_is_alive != value )
	{


		integer	OldValue;
		OldValue = poParent->sum_over_mlHHFatherChildren_is_alive;
		poParent->sum_over_mlHHFatherChildren_is_alive = value;

#line 75 "../code/FamilyGeneral.mpp"
		poParent->children_in_household = ( poParent->sex== FEMALE )? poParent->sum_over_mlHHMotherChildren_is_alive: poParent->sum_over_mlHHFatherChildren_is_alive;
#line 12806 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->sum_over_mlHHFatherChildren_is_alive = value;
}

integer PersonStateFunctions::Set_sum_over_mlHHMotherChildren_is_alive( integer value )
{

	if ( poParent->sum_over_mlHHMotherChildren_is_alive != value )
	{


		integer	OldValue;
		OldValue = poParent->sum_over_mlHHMotherChildren_is_alive;
		poParent->sum_over_mlHHMotherChildren_is_alive = value;

#line 75 "../code/FamilyGeneral.mpp"
		poParent->children_in_household = ( poParent->sex== FEMALE )? poParent->sum_over_mlHHMotherChildren_is_alive: poParent->sum_over_mlHHFatherChildren_is_alive;
#line 12824 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->sum_over_mlHHMotherChildren_is_alive = value;
}

ETHNICITY_SHORT PersonDerivedStates::Set_aggregate_ethnicity_ETHNICITY_SHORT( ETHNICITY_SHORT value )
{

	if ( value < 0 ) {
		value = (ETHNICITY_SHORT) 0;
	}
	else if ( value >= rc_size_ETHNICITY_SHORT ) {
		value = (ETHNICITY_SHORT) ( rc_size_ETHNICITY_SHORT - 1 );
	}
	if ( aggregate_ethnicity_ETHNICITY_SHORT != value ) {

		ETHNICITY_SHORT OldValue;
		OldValue = aggregate_ethnicity_ETHNICITY_SHORT;
		aggregate_ethnicity_ETHNICITY_SHORT = value;
		poStateFunctions->Set_ethnicity_short( aggregate_ethnicity_ETHNICITY_SHORT );
	}
	return aggregate_ethnicity_ETHNICITY_SHORT = value;
}

GEO_NAT PersonDerivedStates::Set_aggregate_geo_birth_GEO_NAT( GEO_NAT value )
{

	if ( value < 0 ) {
		value = (GEO_NAT) 0;
	}
	else if ( value >= rc_size_GEO_NAT ) {
		value = (GEO_NAT) ( rc_size_GEO_NAT - 1 );
	}
	if ( aggregate_geo_birth_GEO_NAT != value ) {

GEO_NAT		OldValue;
		OldValue = aggregate_geo_birth_GEO_NAT;
		aggregate_geo_birth_GEO_NAT = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[16] >= 0 && ( poParent->tableinfo.tab16dim_filter || poParent->table_filter_TabHCIDistrict ) ) {
				if ( poParent->tableinfo.tab16events == -1 || poParent->tableinfo.tab16events < poParent->__events ) {
					if ( poParent->tableinfo.tab16events != -1 ) {
						poParent->tableinfo.tab16open0 = (real) poParent->CloseUpdate( 16, 0, poParent->tableinfo.tab16open0, poParent->tableinfo.tab16close0, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open1 = (real) poParent->CloseUpdate( 16, 1, poParent->tableinfo.tab16open1, poParent->tableinfo.tab16close1, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open2 = (real) poParent->CloseUpdate( 16, 2, poParent->tableinfo.tab16open2, poParent->tableinfo.tab16close2, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open3 = (real) poParent->CloseUpdate( 16, 3, poParent->tableinfo.tab16open3, poParent->tableinfo.tab16close3, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open4 = (real) poParent->CloseUpdate( 16, 4, poParent->tableinfo.tab16open4, poParent->tableinfo.tab16close4, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open5 = (real) poParent->CloseUpdate( 16, 5, poParent->tableinfo.tab16open5, poParent->tableinfo.tab16close5, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open6 = (real) poParent->CloseUpdate( 16, 6, poParent->tableinfo.tab16open6, poParent->tableinfo.tab16close6, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open7 = (real) poParent->CloseUpdate( 16, 7, poParent->tableinfo.tab16open7, poParent->tableinfo.tab16close7, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab16open8 = (real) poParent->CloseUpdate( 16, 8, poParent->tableinfo.tab16open8, poParent->tableinfo.tab16close8, poParent->tableinfo.tab16dim_filter, poParent->tableinfo.tab16cell, poParent->tableinfo.tab16total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab16events = poParent->__events;
					poParent->tableinfo.tab16close0 = (real) transitions_integer_age_14_15;
					poParent->tableinfo.tab16close1 = (real) transitions_is_alive_TRUE_FALSE;
					poParent->tableinfo.tab16close2 = (real) transitions_survived_early_years_FALSE_TRUE;
					poParent->tableinfo.tab16close3 = (real) value_at_transitions_is_alive_TRUE_FALSE_adult_survival;
					poParent->tableinfo.tab16close4 = (real) value_at_transitions_is_alive_TRUE_FALSE_ind_hci;
					poParent->tableinfo.tab16close5 = (real) value_at_transitions_is_alive_TRUE_FALSE_is_stunted;
					poParent->tableinfo.tab16close6 = (real) value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling;
					poParent->tableinfo.tab16close7 = (real) value_at_transitions_is_alive_TRUE_FALSE_survived_early_years;
					poParent->tableinfo.tab16close8 = (real) value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling;
					poParent->tableinfo.tab16cell = (index) poParent->TableCellIndex( 16, 3,
						poParent->sex,
						OldValue,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab16total_cell = (index) poParent->tableinfo.tab16cell;
					}
					else {
						poParent->tableinfo.tab16total_cell = (index) poParent->TableTotalCellIndex( 16, 3,
						poParent->sex,
						OldValue,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab16dim_filter = poParent->table_filter_TabHCIDistrict;
				}
			}
		}
	}
	return aggregate_geo_birth_GEO_NAT = value;
}

REGION_INT PersonDerivedStates::Set_aggregate_geo_birth_REGION_INT( REGION_INT value )
{

	if ( value < 0 ) {
		value = (REGION_INT) 0;
	}
	else if ( value >= rc_size_REGION_INT ) {
		value = (REGION_INT) ( rc_size_REGION_INT - 1 );
	}
	if ( aggregate_geo_birth_REGION_INT != value ) {

		REGION_INT OldValue;
		OldValue = aggregate_geo_birth_REGION_INT;
		aggregate_geo_birth_REGION_INT = value;
		poStateFunctions->Set_region_birth( aggregate_geo_birth_REGION_INT );
	}
	return aggregate_geo_birth_REGION_INT = value;
}

GEO_NAT PersonDerivedStates::Set_aggregate_geo_GEO_NAT( GEO_NAT value )
{

	if ( value < 0 ) {
		value = (GEO_NAT) 0;
	}
	else if ( value >= rc_size_GEO_NAT ) {
		value = (GEO_NAT) ( rc_size_GEO_NAT - 1 );
	}
	if ( aggregate_geo_GEO_NAT != value ) {

		GEO_NAT OldValue;
		OldValue = aggregate_geo_GEO_NAT;
		aggregate_geo_GEO_NAT = value;
		poStateFunctions->Set_geo_nat( aggregate_geo_GEO_NAT );
	}
	return aggregate_geo_GEO_NAT = value;
}

REGION_INT PersonDerivedStates::Set_aggregate_geo_REGION_INT( REGION_INT value )
{

	if ( value < 0 ) {
		value = (REGION_INT) 0;
	}
	else if ( value >= rc_size_REGION_INT ) {
		value = (REGION_INT) ( rc_size_REGION_INT - 1 );
	}
	if ( aggregate_geo_REGION_INT != value ) {

		REGION_INT OldValue;
		OldValue = aggregate_geo_REGION_INT;
		aggregate_geo_REGION_INT = value;
		poStateFunctions->Set_region_int( aggregate_geo_REGION_INT );
	}
	return aggregate_geo_REGION_INT = value;
}

REGION_NAT PersonDerivedStates::Set_aggregate_region_birth_REGION_NAT( REGION_NAT value )
{

	if ( value < 0 ) {
		value = (REGION_NAT) 0;
	}
	else if ( value >= rc_size_REGION_NAT ) {
		value = (REGION_NAT) ( rc_size_REGION_NAT - 1 );
	}
	if ( aggregate_region_birth_REGION_NAT != value ) {

REGION_NAT		OldValue;
		OldValue = aggregate_region_birth_REGION_NAT;
		aggregate_region_birth_REGION_NAT = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[6] >= 0 && ( poParent->tableinfo.tab6dim_filter || poParent->table_filter_TabEducFateDistrYob ) ) {
				if ( poParent->tableinfo.tab6events == -1 || poParent->tableinfo.tab6events < poParent->__events ) {
					if ( poParent->tableinfo.tab6events != -1 ) {
						poStateFunctions->Set_table_unit_TabEducFateDistrYob(poParent->table_unit_TabEducFateDistrYob+1);
						poParent->tableinfo.tab6close0 = (real) poParent->table_unit_TabEducFateDistrYob;
						poParent->tableinfo.tab6open0 = (real) poParent->CloseUpdate( 6, 0, poParent->tableinfo.tab6open0, poParent->tableinfo.tab6close0, poParent->tableinfo.tab6dim_filter, poParent->tableinfo.tab6cell, poParent->tableinfo.tab6total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab6events = poParent->__events;
					poParent->tableinfo.tab6close0 = (real) poParent->table_unit_TabEducFateDistrYob;
					poParent->tableinfo.tab6cell = (index) poParent->TableCellIndex( 6, 3,
						OldValue,
						poParent->year_of_birth - r_min_ALL_YEAR_RANGE,
						poParent->educ_one_fate );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab6total_cell = (index) poParent->tableinfo.tab6cell;
					}
					else {
						poParent->tableinfo.tab6total_cell = (index) poParent->TableTotalCellIndex( 6, 3,
						OldValue,
						poParent->year_of_birth - r_min_ALL_YEAR_RANGE,
						poParent->educ_one_fate );
					}
					poParent->tableinfo.tab6dim_filter = poParent->table_filter_TabEducFateDistrYob;
				}
			}
			if ( gplTableIndexes[15] >= 0 && ( poParent->tableinfo.tab15dim_filter || poParent->table_filter_tabHCI ) ) {
				if ( poParent->tableinfo.tab15events == -1 || poParent->tableinfo.tab15events < poParent->__events ) {
					if ( poParent->tableinfo.tab15events != -1 ) {
						poParent->tableinfo.tab15open0 = (real) poParent->CloseUpdate( 15, 0, poParent->tableinfo.tab15open0, poParent->tableinfo.tab15close0, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open1 = (real) poParent->CloseUpdate( 15, 1, poParent->tableinfo.tab15open1, poParent->tableinfo.tab15close1, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open2 = (real) poParent->CloseUpdate( 15, 2, poParent->tableinfo.tab15open2, poParent->tableinfo.tab15close2, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open3 = (real) poParent->CloseUpdate( 15, 3, poParent->tableinfo.tab15open3, poParent->tableinfo.tab15close3, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open4 = (real) poParent->CloseUpdate( 15, 4, poParent->tableinfo.tab15open4, poParent->tableinfo.tab15close4, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open5 = (real) poParent->CloseUpdate( 15, 5, poParent->tableinfo.tab15open5, poParent->tableinfo.tab15close5, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open6 = (real) poParent->CloseUpdate( 15, 6, poParent->tableinfo.tab15open6, poParent->tableinfo.tab15close6, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open7 = (real) poParent->CloseUpdate( 15, 7, poParent->tableinfo.tab15open7, poParent->tableinfo.tab15close7, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab15open8 = (real) poParent->CloseUpdate( 15, 8, poParent->tableinfo.tab15open8, poParent->tableinfo.tab15close8, poParent->tableinfo.tab15dim_filter, poParent->tableinfo.tab15cell, poParent->tableinfo.tab15total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab15events = poParent->__events;
					poParent->tableinfo.tab15close0 = (real) transitions_integer_age_14_15;
					poParent->tableinfo.tab15close1 = (real) transitions_is_alive_TRUE_FALSE;
					poParent->tableinfo.tab15close2 = (real) transitions_survived_early_years_FALSE_TRUE;
					poParent->tableinfo.tab15close3 = (real) value_at_transitions_is_alive_TRUE_FALSE_adult_survival;
					poParent->tableinfo.tab15close4 = (real) value_at_transitions_is_alive_TRUE_FALSE_ind_hci;
					poParent->tableinfo.tab15close5 = (real) value_at_transitions_is_alive_TRUE_FALSE_is_stunted;
					poParent->tableinfo.tab15close6 = (real) value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling;
					poParent->tableinfo.tab15close7 = (real) value_at_transitions_is_alive_TRUE_FALSE_survived_early_years;
					poParent->tableinfo.tab15close8 = (real) value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling;
					poParent->tableinfo.tab15cell = (index) poParent->TableCellIndex( 15, 3,
						poParent->sex,
						OldValue,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab15total_cell = (index) poParent->tableinfo.tab15cell;
					}
					else {
						poParent->tableinfo.tab15total_cell = (index) poParent->TableTotalCellIndex( 15, 3,
						poParent->sex,
						OldValue,
						poParent->tab_sim_yob - r_min_SIM_YEAR_RANGE );
					}
					poParent->tableinfo.tab15dim_filter = poParent->table_filter_tabHCI;
				}
			}
		}
	}
	return aggregate_region_birth_REGION_NAT = value;
}

REGION_NAT PersonDerivedStates::Set_aggregate_region_int_REGION_NAT( REGION_NAT value )
{

	if ( value < 0 ) {
		value = (REGION_NAT) 0;
	}
	else if ( value >= rc_size_REGION_NAT ) {
		value = (REGION_NAT) ( rc_size_REGION_NAT - 1 );
	}
	if ( aggregate_region_int_REGION_NAT != value ) {

		REGION_NAT OldValue;
		OldValue = aggregate_region_int_REGION_NAT;
		aggregate_region_int_REGION_NAT = value;
		poStateFunctions->Set_region_nat( aggregate_region_int_REGION_NAT );
	}
	return aggregate_region_int_REGION_NAT = value;
}

TIME PersonDerivedStates::Set_duration_( TIME value )
{

	if ( duration_ != value ) {

TIME		OldValue;
		OldValue = duration_;
		duration_ = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[0] >= 0 ) {
				if ( poParent->tableinfo.tab0events == poParent->__events ) {
					poParent->tableinfo.tab0close0 = (real) duration_;
				}
			}
			if ( gplTableIndexes[1] >= 0 ) {
				if ( poParent->tableinfo.tab1events == poParent->__events ) {
					poParent->tableinfo.tab1close0 = (real) duration_;
				}
			}
			if ( gplTableIndexes[7] >= 0 ) {
				if ( poParent->tableinfo.tab7events == poParent->__events ) {
					poParent->tableinfo.tab7close0 = (real) duration_;
				}
			}
			if ( gplTableIndexes[8] >= 0 ) {
				if ( poParent->tableinfo.tab8events == poParent->__events ) {
					poParent->tableinfo.tab8close0 = (real) duration_;
				}
			}
			if ( gplTableIndexes[11] >= 0 ) {
				if ( poParent->tableinfo.tab11events == poParent->__events ) {
					poParent->tableinfo.tab11close0 = (real) duration_;
				}
			}
			if ( gplTableIndexes[13] >= 0 ) {
				if ( poParent->tableinfo.tab13events == poParent->__events ) {
					poParent->tableinfo.tab13close0 = (real) duration_;
				}
			}
			if ( gplTableIndexes[14] >= 0 ) {
				if ( poParent->tableinfo.tab14events == poParent->__events ) {
					poParent->tableinfo.tab14close0 = (real) duration_;
				}
			}
			if ( gplTableIndexes[17] >= 0 ) {
				if ( poParent->tableinfo.tab17events == poParent->__events ) {
					poParent->tableinfo.tab17close0 = (real) duration_;
				}
			}
			if ( gplTableIndexes[22] >= 0 ) {
				if ( poParent->tableinfo.tab22events == poParent->__events ) {
					poParent->tableinfo.tab22close0 = (real) duration_;
				}
			}
			if ( gplTableIndexes[23] >= 0 ) {
				if ( poParent->tableinfo.tab23events == poParent->__events ) {
					poParent->tableinfo.tab23close0 = (real) duration_;
				}
			}
			if ( gplTableIndexes[24] >= 0 ) {
				if ( poParent->tableinfo.tab24events == poParent->__events ) {
					poParent->tableinfo.tab24close0 = (real) duration_;
				}
			}
		}
	}
	return duration_ = value;
}

TIME PersonDerivedStates::Set_duration_educ_one_fate_EOL_HIGH( TIME value )
{

	if ( duration_educ_one_fate_EOL_HIGH != value ) {

TIME		OldValue;
		OldValue = duration_educ_one_fate_EOL_HIGH;
		duration_educ_one_fate_EOL_HIGH = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[1] >= 0 ) {
				if ( poParent->tableinfo.tab1events == poParent->__events ) {
					poParent->tableinfo.tab1close1 = (real) duration_educ_one_fate_EOL_HIGH;
				}
			}
			if ( gplTableIndexes[7] >= 0 ) {
				if ( poParent->tableinfo.tab7events == poParent->__events ) {
					poParent->tableinfo.tab7close1 = (real) duration_educ_one_fate_EOL_HIGH;
				}
			}
			if ( gplTableIndexes[8] >= 0 ) {
				if ( poParent->tableinfo.tab8events == poParent->__events ) {
					poParent->tableinfo.tab8close1 = (real) duration_educ_one_fate_EOL_HIGH;
				}
			}
		}
	}
	return duration_educ_one_fate_EOL_HIGH = value;
}

TIME PersonDerivedStates::Set_duration_educ_one_fate_EOL_LOW( TIME value )
{

	if ( duration_educ_one_fate_EOL_LOW != value ) {

TIME		OldValue;
		OldValue = duration_educ_one_fate_EOL_LOW;
		duration_educ_one_fate_EOL_LOW = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[1] >= 0 ) {
				if ( poParent->tableinfo.tab1events == poParent->__events ) {
					poParent->tableinfo.tab1close2 = (real) duration_educ_one_fate_EOL_LOW;
				}
			}
			if ( gplTableIndexes[7] >= 0 ) {
				if ( poParent->tableinfo.tab7events == poParent->__events ) {
					poParent->tableinfo.tab7close2 = (real) duration_educ_one_fate_EOL_LOW;
				}
			}
			if ( gplTableIndexes[8] >= 0 ) {
				if ( poParent->tableinfo.tab8events == poParent->__events ) {
					poParent->tableinfo.tab8close2 = (real) duration_educ_one_fate_EOL_LOW;
				}
			}
		}
	}
	return duration_educ_one_fate_EOL_LOW = value;
}

TIME PersonDerivedStates::Set_duration_educ_one_fate_EOL_MEDIUM( TIME value )
{

	if ( duration_educ_one_fate_EOL_MEDIUM != value ) {

TIME		OldValue;
		OldValue = duration_educ_one_fate_EOL_MEDIUM;
		duration_educ_one_fate_EOL_MEDIUM = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[1] >= 0 ) {
				if ( poParent->tableinfo.tab1events == poParent->__events ) {
					poParent->tableinfo.tab1close3 = (real) duration_educ_one_fate_EOL_MEDIUM;
				}
			}
			if ( gplTableIndexes[7] >= 0 ) {
				if ( poParent->tableinfo.tab7events == poParent->__events ) {
					poParent->tableinfo.tab7close3 = (real) duration_educ_one_fate_EOL_MEDIUM;
				}
			}
			if ( gplTableIndexes[8] >= 0 ) {
				if ( poParent->tableinfo.tab8events == poParent->__events ) {
					poParent->tableinfo.tab8close3 = (real) duration_educ_one_fate_EOL_MEDIUM;
				}
			}
		}
	}
	return duration_educ_one_fate_EOL_MEDIUM = value;
}

TIME PersonDerivedStates::Set_duration_educ_one_level_EOL_HIGH( TIME value )
{

	if ( duration_educ_one_level_EOL_HIGH != value ) {

TIME		OldValue;
		OldValue = duration_educ_one_level_EOL_HIGH;
		duration_educ_one_level_EOL_HIGH = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[2] >= 0 ) {
				if ( poParent->tableinfo.tab2events == poParent->__events ) {
					poParent->tableinfo.tab2close0 = (real) duration_educ_one_level_EOL_HIGH;
				}
			}
		}
	}
	return duration_educ_one_level_EOL_HIGH = value;
}

TIME PersonDerivedStates::Set_duration_educ_one_level_EOL_LOW( TIME value )
{

	if ( duration_educ_one_level_EOL_LOW != value ) {

TIME		OldValue;
		OldValue = duration_educ_one_level_EOL_LOW;
		duration_educ_one_level_EOL_LOW = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[2] >= 0 ) {
				if ( poParent->tableinfo.tab2events == poParent->__events ) {
					poParent->tableinfo.tab2close1 = (real) duration_educ_one_level_EOL_LOW;
				}
			}
		}
	}
	return duration_educ_one_level_EOL_LOW = value;
}

TIME PersonDerivedStates::Set_duration_educ_one_level_EOL_MEDIUM( TIME value )
{

	if ( duration_educ_one_level_EOL_MEDIUM != value ) {

TIME		OldValue;
		OldValue = duration_educ_one_level_EOL_MEDIUM;
		duration_educ_one_level_EOL_MEDIUM = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[2] >= 0 ) {
				if ( poParent->tableinfo.tab2events == poParent->__events ) {
					poParent->tableinfo.tab2close2 = (real) duration_educ_one_level_EOL_MEDIUM;
				}
			}
		}
	}
	return duration_educ_one_level_EOL_MEDIUM = value;
}

TIME PersonDerivedStates::Set_duration_educ_two_level_ETL_GRADUATE( TIME value )
{

	if ( duration_educ_two_level_ETL_GRADUATE != value ) {

TIME		OldValue;
		OldValue = duration_educ_two_level_ETL_GRADUATE;
		duration_educ_two_level_ETL_GRADUATE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[13] >= 0 ) {
				if ( poParent->tableinfo.tab13events == poParent->__events ) {
					poParent->tableinfo.tab13close1 = (real) duration_educ_two_level_ETL_GRADUATE;
				}
			}
		}
	}
	return duration_educ_two_level_ETL_GRADUATE = value;
}

TIME PersonDerivedStates::Set_duration_educ_two_level_ETL_INCOMPLETE( TIME value )
{

	if ( duration_educ_two_level_ETL_INCOMPLETE != value ) {

TIME		OldValue;
		OldValue = duration_educ_two_level_ETL_INCOMPLETE;
		duration_educ_two_level_ETL_INCOMPLETE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[13] >= 0 ) {
				if ( poParent->tableinfo.tab13events == poParent->__events ) {
					poParent->tableinfo.tab13close2 = (real) duration_educ_two_level_ETL_INCOMPLETE;
				}
			}
		}
	}
	return duration_educ_two_level_ETL_INCOMPLETE = value;
}

TIME PersonDerivedStates::Set_duration_educ_two_level_ETL_NEVER( TIME value )
{

	if ( duration_educ_two_level_ETL_NEVER != value ) {

TIME		OldValue;
		OldValue = duration_educ_two_level_ETL_NEVER;
		duration_educ_two_level_ETL_NEVER = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[13] >= 0 ) {
				if ( poParent->tableinfo.tab13events == poParent->__events ) {
					poParent->tableinfo.tab13close3 = (real) duration_educ_two_level_ETL_NEVER;
				}
			}
		}
	}
	return duration_educ_two_level_ETL_NEVER = value;
}

TIME PersonDerivedStates::Set_duration_ever_union_TRUE( TIME value )
{

	if ( duration_ever_union_TRUE != value ) {

TIME		OldValue;
		OldValue = duration_ever_union_TRUE;
		duration_ever_union_TRUE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[14] >= 0 ) {
				if ( poParent->tableinfo.tab14events == poParent->__events ) {
					poParent->tableinfo.tab14close1 = (real) duration_ever_union_TRUE;
				}
			}
		}
	}
	return duration_ever_union_TRUE = value;
}

TIME PersonDerivedStates::Set_duration_got_prenat_care_GPC_YES( TIME value )
{

	if ( duration_got_prenat_care_GPC_YES != value ) {

TIME		OldValue;
		OldValue = duration_got_prenat_care_GPC_YES;
		duration_got_prenat_care_GPC_YES = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[22] >= 0 ) {
				if ( poParent->tableinfo.tab22events == poParent->__events ) {
					poParent->tableinfo.tab22close1 = (real) duration_got_prenat_care_GPC_YES;
				}
			}
		}
	}
	return duration_got_prenat_care_GPC_YES = value;
}

TIME PersonDerivedStates::Set_duration_has_spouse_TRUE( TIME value )
{

	if ( duration_has_spouse_TRUE != value ) {

TIME		OldValue;
		OldValue = duration_has_spouse_TRUE;
		duration_has_spouse_TRUE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[14] >= 0 ) {
				if ( poParent->tableinfo.tab14events == poParent->__events ) {
					poParent->tableinfo.tab14close2 = (real) duration_has_spouse_TRUE;
				}
			}
		}
	}
	return duration_has_spouse_TRUE = value;
}

TIME PersonDerivedStates::Set_duration_is_aded0_TRUE( TIME value )
{

	if ( duration_is_aded0_TRUE != value ) {

TIME		OldValue;
		OldValue = duration_is_aded0_TRUE;
		duration_is_aded0_TRUE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[5] >= 0 ) {
				if ( poParent->tableinfo.tab5events == poParent->__events ) {
					poParent->tableinfo.tab5close0 = (real) duration_is_aded0_TRUE;
				}
			}
		}
	}
	return duration_is_aded0_TRUE = value;
}

TIME PersonDerivedStates::Set_duration_is_aded1_TRUE( TIME value )
{

	if ( duration_is_aded1_TRUE != value ) {

TIME		OldValue;
		OldValue = duration_is_aded1_TRUE;
		duration_is_aded1_TRUE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[5] >= 0 ) {
				if ( poParent->tableinfo.tab5events == poParent->__events ) {
					poParent->tableinfo.tab5close1 = (real) duration_is_aded1_TRUE;
				}
			}
		}
	}
	return duration_is_aded1_TRUE = value;
}

TIME PersonDerivedStates::Set_duration_is_aded2_TRUE( TIME value )
{

	if ( duration_is_aded2_TRUE != value ) {

TIME		OldValue;
		OldValue = duration_is_aded2_TRUE;
		duration_is_aded2_TRUE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[5] >= 0 ) {
				if ( poParent->tableinfo.tab5events == poParent->__events ) {
					poParent->tableinfo.tab5close2 = (real) duration_is_aded2_TRUE;
				}
			}
		}
	}
	return duration_is_aded2_TRUE = value;
}

TIME PersonDerivedStates::Set_duration_is_child_TRUE( TIME value )
{

	if ( duration_is_child_TRUE != value ) {

TIME		OldValue;
		OldValue = duration_is_child_TRUE;
		duration_is_child_TRUE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[5] >= 0 ) {
				if ( poParent->tableinfo.tab5events == poParent->__events ) {
					poParent->tableinfo.tab5close3 = (real) duration_is_child_TRUE;
				}
			}
		}
	}
	return duration_is_child_TRUE = value;
}

TIME PersonDerivedStates::Set_duration_is_immunized_TRUE( TIME value )
{

	if ( duration_is_immunized_TRUE != value ) {

TIME		OldValue;
		OldValue = duration_is_immunized_TRUE;
		duration_is_immunized_TRUE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[23] >= 0 ) {
				if ( poParent->tableinfo.tab23events == poParent->__events ) {
					poParent->tableinfo.tab23close1 = (real) duration_is_immunized_TRUE;
				}
			}
			if ( gplTableIndexes[24] >= 0 ) {
				if ( poParent->tableinfo.tab24events == poParent->__events ) {
					poParent->tableinfo.tab24close1 = (real) duration_is_immunized_TRUE;
				}
			}
		}
	}
	return duration_is_immunized_TRUE = value;
}

TIME PersonDerivedStates::Set_duration_is_old_TRUE( TIME value )
{

	if ( duration_is_old_TRUE != value ) {

TIME		OldValue;
		OldValue = duration_is_old_TRUE;
		duration_is_old_TRUE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[5] >= 0 ) {
				if ( poParent->tableinfo.tab5events == poParent->__events ) {
					poParent->tableinfo.tab5close4 = (real) duration_is_old_TRUE;
				}
			}
		}
	}
	return duration_is_old_TRUE = value;
}

TIME PersonDerivedStates::Set_duration_tab_out_of_primary_TRUE( TIME value )
{

	if ( duration_tab_out_of_primary_TRUE != value ) {

TIME		OldValue;
		OldValue = duration_tab_out_of_primary_TRUE;
		duration_tab_out_of_primary_TRUE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[11] >= 0 ) {
				if ( poParent->tableinfo.tab11events == poParent->__events ) {
					poParent->tableinfo.tab11close1 = (real) duration_tab_out_of_primary_TRUE;
				}
			}
		}
	}
	return duration_tab_out_of_primary_TRUE = value;
}

TIME PersonDerivedStates::Set_duration_years_preschool_0( TIME value )
{

	if ( duration_years_preschool_0 != value ) {

TIME		OldValue;
		OldValue = duration_years_preschool_0;
		duration_years_preschool_0 = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[0] >= 0 ) {
				if ( poParent->tableinfo.tab0events == poParent->__events ) {
					poParent->tableinfo.tab0close1 = (real) duration_years_preschool_0;
				}
			}
		}
	}
	return duration_years_preschool_0 = value;
}

TIME PersonDerivedStates::Set_duration_years_preschool_1( TIME value )
{

	if ( duration_years_preschool_1 != value ) {

TIME		OldValue;
		OldValue = duration_years_preschool_1;
		duration_years_preschool_1 = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[0] >= 0 ) {
				if ( poParent->tableinfo.tab0events == poParent->__events ) {
					poParent->tableinfo.tab0close2 = (real) duration_years_preschool_1;
				}
			}
		}
	}
	return duration_years_preschool_1 = value;
}

TIME PersonDerivedStates::Set_duration_years_preschool_2( TIME value )
{

	if ( duration_years_preschool_2 != value ) {

TIME		OldValue;
		OldValue = duration_years_preschool_2;
		duration_years_preschool_2 = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[0] >= 0 ) {
				if ( poParent->tableinfo.tab0events == poParent->__events ) {
					poParent->tableinfo.tab0close3 = (real) duration_years_preschool_2;
				}
			}
		}
	}
	return duration_years_preschool_2 = value;
}

counter PersonDerivedStates::Set_entrances_educ_one_grade_attended_1( counter value )
{

	if ( entrances_educ_one_grade_attended_1 != value ) {

counter		OldValue;
		OldValue = entrances_educ_one_grade_attended_1;
		entrances_educ_one_grade_attended_1 = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[9] >= 0 ) {
				if ( poParent->tableinfo.tab9events == poParent->__events ) {
					poParent->tableinfo.tab9close0 = (real) entrances_educ_one_grade_attended_1;
				}
			}
		}
	}
	return entrances_educ_one_grade_attended_1 = value;
}

counter PersonDerivedStates::Set_entrances_has_passed_primary_TRUE( counter value )
{

	if ( entrances_has_passed_primary_TRUE != value ) {

counter		OldValue;
		OldValue = entrances_has_passed_primary_TRUE;
		entrances_has_passed_primary_TRUE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[10] >= 0 ) {
				if ( poParent->tableinfo.tab10events == poParent->__events ) {
					poParent->tableinfo.tab10close0 = (real) entrances_has_passed_primary_TRUE;
				}
			}
		}
	}
	return entrances_has_passed_primary_TRUE = value;
}

MOTH_AGEGR_PART PersonDerivedStates::Set_split_age_last_birth_MOTH_AGEGR_PART( MOTH_AGEGR_PART value )
{

	if ( split_age_last_birth_MOTH_AGEGR_PART != value ) {

		MOTH_AGEGR_PART OldValue;
		OldValue = split_age_last_birth_MOTH_AGEGR_PART;
		split_age_last_birth_MOTH_AGEGR_PART = value;
		poStateFunctions->Set_moth_agegr_part( split_age_last_birth_MOTH_AGEGR_PART );
	}
	return split_age_last_birth_MOTH_AGEGR_PART = value;
}

CHILD_AGEGR_PART PersonDerivedStates::Set_split_age_youngest_child_CHILD_AGEGR_PART( CHILD_AGEGR_PART value )
{

	if ( split_age_youngest_child_CHILD_AGEGR_PART != value ) {

		CHILD_AGEGR_PART OldValue;
		OldValue = split_age_youngest_child_CHILD_AGEGR_PART;
		split_age_youngest_child_CHILD_AGEGR_PART = value;
		poStateFunctions->Set_child_agegr_part( split_age_youngest_child_CHILD_AGEGR_PART );
	}
	return split_age_youngest_child_CHILD_AGEGR_PART = value;
}

TAB_AGEC5 PersonDerivedStates::Set_split_integer_age_TAB_AGEC5( TAB_AGEC5 value )
{

	if ( split_integer_age_TAB_AGEC5 != value ) {

TAB_AGEC5		OldValue;
		OldValue = split_integer_age_TAB_AGEC5;
		split_integer_age_TAB_AGEC5 = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[2] >= 0 && ( poParent->tableinfo.tab2dim_filter || poParent->table_filter_PopPyramidByEduc ) ) {
				if ( poParent->tableinfo.tab2events == -1 || poParent->tableinfo.tab2events < poParent->__events ) {
					if ( poParent->tableinfo.tab2events != -1 ) {
						poParent->tableinfo.tab2open0 = (real) poParent->CloseUpdate( 2, 0, poParent->tableinfo.tab2open0, poParent->tableinfo.tab2close0, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab2open1 = (real) poParent->CloseUpdate( 2, 1, poParent->tableinfo.tab2open1, poParent->tableinfo.tab2close1, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
						poParent->tableinfo.tab2open2 = (real) poParent->CloseUpdate( 2, 2, poParent->tableinfo.tab2open2, poParent->tableinfo.tab2close2, poParent->tableinfo.tab2dim_filter, poParent->tableinfo.tab2cell, poParent->tableinfo.tab2total_cell, poParent->actor_weight, poParent->actor_subsample_weight );
					}
					poParent->tableinfo.tab2events = poParent->__events;
					poParent->tableinfo.tab2close0 = (real) Get_duration_educ_one_level_EOL_HIGH();
					poParent->tableinfo.tab2close1 = (real) Get_duration_educ_one_level_EOL_LOW();
					poParent->tableinfo.tab2close2 = (real) Get_duration_educ_one_level_EOL_MEDIUM();
					poParent->tableinfo.tab2cell = (index) poParent->TableCellIndex( 2, 4,
						poParent->sel_years,
						poParent->sex,
						poParent->geo_nat,
						OldValue );
					if ( gpModelApp->m_nSubSamples == 1 ) {
						poParent->tableinfo.tab2total_cell = (index) poParent->tableinfo.tab2cell;
					}
					else {
						poParent->tableinfo.tab2total_cell = (index) poParent->TableTotalCellIndex( 2, 4,
						poParent->sel_years,
						poParent->sex,
						poParent->geo_nat,
						OldValue );
					}
					poParent->tableinfo.tab2dim_filter = poParent->table_filter_PopPyramidByEduc;
				}
			}
		}
	}
	return split_integer_age_TAB_AGEC5 = value;
}

DUR_TIME_IN_PARITY PersonDerivedStates::Set_sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY( DUR_TIME_IN_PARITY value )
{

	if ( sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY != value ) {

		DUR_TIME_IN_PARITY OldValue;
		OldValue = sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY;
		sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY = value;
		poStateFunctions->Set_time_in_parity( sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY );
	}
	return sssplit_active_spell_duration_this_parity_spell_TRUE_DUR_TIME_IN_PARITY = value;
}

AGE_MIG_PART PersonDerivedStates::Set_sssplit_age_AGE_MIG_PART( AGE_MIG_PART value )
{

	if ( sssplit_age_AGE_MIG_PART != value ) {

		AGE_MIG_PART OldValue;
		OldValue = sssplit_age_AGE_MIG_PART;
		sssplit_age_AGE_MIG_PART = value;
		poStateFunctions->Set_age_mig( COERCE(AGE_MIG,sssplit_age_AGE_MIG_PART) );
	}
	return sssplit_age_AGE_MIG_PART = value;
}

counter PersonDerivedStates::Set_transitions_integer_age_14_15( counter value )
{

	if ( transitions_integer_age_14_15 != value ) {

counter		OldValue;
		OldValue = transitions_integer_age_14_15;
		transitions_integer_age_14_15 = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[15] >= 0 ) {
				if ( poParent->tableinfo.tab15events == poParent->__events ) {
					poParent->tableinfo.tab15close0 = (real) transitions_integer_age_14_15;
				}
			}
			if ( gplTableIndexes[16] >= 0 ) {
				if ( poParent->tableinfo.tab16events == poParent->__events ) {
					poParent->tableinfo.tab16close0 = (real) transitions_integer_age_14_15;
				}
			}
		}
	}
	return transitions_integer_age_14_15 = value;
}

counter PersonDerivedStates::Set_transitions_is_alive_FALSE_TRUE( counter value )
{

	if ( transitions_is_alive_FALSE_TRUE != value ) {

counter		OldValue;
		OldValue = transitions_is_alive_FALSE_TRUE;
		transitions_is_alive_FALSE_TRUE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[18] >= 0 ) {
				if ( poParent->tableinfo.tab18events == poParent->__events ) {
					poParent->tableinfo.tab18close0 = (real) transitions_is_alive_FALSE_TRUE;
				}
			}
		}
	}
	return transitions_is_alive_FALSE_TRUE = value;
}

counter PersonDerivedStates::Set_transitions_is_alive_TRUE_FALSE( counter value )
{

	if ( transitions_is_alive_TRUE_FALSE != value ) {

counter		OldValue;
		OldValue = transitions_is_alive_TRUE_FALSE;
		transitions_is_alive_TRUE_FALSE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[15] >= 0 ) {
				if ( poParent->tableinfo.tab15events == poParent->__events ) {
					poParent->tableinfo.tab15close1 = (real) transitions_is_alive_TRUE_FALSE;
				}
			}
			if ( gplTableIndexes[16] >= 0 ) {
				if ( poParent->tableinfo.tab16events == poParent->__events ) {
					poParent->tableinfo.tab16close1 = (real) transitions_is_alive_TRUE_FALSE;
				}
			}
		}
	}
	return transitions_is_alive_TRUE_FALSE = value;
}

counter PersonDerivedStates::Set_transitions_survived_early_years_FALSE_TRUE( counter value )
{

	if ( transitions_survived_early_years_FALSE_TRUE != value ) {

counter		OldValue;
		OldValue = transitions_survived_early_years_FALSE_TRUE;
		transitions_survived_early_years_FALSE_TRUE = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[15] >= 0 ) {
				if ( poParent->tableinfo.tab15events == poParent->__events ) {
					poParent->tableinfo.tab15close2 = (real) transitions_survived_early_years_FALSE_TRUE;
				}
			}
			if ( gplTableIndexes[16] >= 0 ) {
				if ( poParent->tableinfo.tab16events == poParent->__events ) {
					poParent->tableinfo.tab16close2 = (real) transitions_survived_early_years_FALSE_TRUE;
				}
			}
		}
	}
	return transitions_survived_early_years_FALSE_TRUE = value;
}

logical PersonDerivedStates::Set_trigger_changes_calendar_year( logical value )
{

	value = value != 0;
	if ( trigger_changes_calendar_year != value ) {
		poStateFunctions->Set_table_filter_TabEducFateByGroup( poParent->in_projected_time&& value&& poParent->integer_age== 0&& poParent->creation_type== CT_BIRTH );
		poStateFunctions->Set_table_filter_TabEducFateDistrYob( value&& poParent->integer_age== 0&& poParent->is_resident );
		poStateFunctions->Set_table_filter_TabPrimarySchoolPlanning( poParent->in_projected_time&& value&& poParent->is_educ_one_student&& poParent->is_resident );
		poStateFunctions->Set_table_filter_TabSchool2TrackTab( value&& poParent->in_projected_time&& poParent->is_resident&& poParent->educ_two_status== ETS_ATTEND );

		logical OldValue;
		OldValue = trigger_changes_calendar_year;
		trigger_changes_calendar_year = value;
	}
	return trigger_changes_calendar_year = value;
}

logical PersonDerivedStates::Set_trigger_changes_geo( logical value )
{

	value = value != 0;
	if ( trigger_changes_geo != value ) {
		poStateFunctions->Set_table_filter_tabMigrationOriginDestination( WITHIN(SIM_YEAR_RANGE,poParent->calendar_year)&& value&& poParent->in_selected_year );

		logical OldValue;
		OldValue = trigger_changes_geo;
		trigger_changes_geo = value;
	}
	return trigger_changes_geo = value;
}

logical PersonDerivedStates::Set_trigger_entrances_is_ready_for_birthtables_TRUE( logical value )
{

	value = value != 0;
	if ( trigger_entrances_is_ready_for_birthtables_TRUE != value ) {
		poStateFunctions->Set_table_filter_tabStuntingSexRegMotherYob( poParent->in_projected_time&& poParent->is_resident&& poParent->creation_type== CT_BIRTH&& value );
		poStateFunctions->Set_table_filter_tabStuntingSexRegYob( poParent->in_projected_time&& poParent->is_resident&& poParent->creation_type== CT_BIRTH&& value );

		logical OldValue;
		OldValue = trigger_entrances_is_ready_for_birthtables_TRUE;
		trigger_entrances_is_ready_for_birthtables_TRUE = value;
	}
	return trigger_entrances_is_ready_for_birthtables_TRUE = value;
}

real PersonDerivedStates::Set_value_at_transitions_is_alive_TRUE_FALSE_adult_survival( real value )
{

	if ( value_at_transitions_is_alive_TRUE_FALSE_adult_survival != value ) {

real		OldValue;
		OldValue = value_at_transitions_is_alive_TRUE_FALSE_adult_survival;
		value_at_transitions_is_alive_TRUE_FALSE_adult_survival = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[15] >= 0 ) {
				if ( poParent->tableinfo.tab15events == poParent->__events ) {
					poParent->tableinfo.tab15close3 = (real) value_at_transitions_is_alive_TRUE_FALSE_adult_survival;
				}
			}
			if ( gplTableIndexes[16] >= 0 ) {
				if ( poParent->tableinfo.tab16events == poParent->__events ) {
					poParent->tableinfo.tab16close3 = (real) value_at_transitions_is_alive_TRUE_FALSE_adult_survival;
				}
			}
		}
	}
	return value_at_transitions_is_alive_TRUE_FALSE_adult_survival = value;
}

real PersonDerivedStates::Set_value_at_transitions_is_alive_TRUE_FALSE_ind_hci( real value )
{

	if ( value_at_transitions_is_alive_TRUE_FALSE_ind_hci != value ) {

real		OldValue;
		OldValue = value_at_transitions_is_alive_TRUE_FALSE_ind_hci;
		value_at_transitions_is_alive_TRUE_FALSE_ind_hci = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[15] >= 0 ) {
				if ( poParent->tableinfo.tab15events == poParent->__events ) {
					poParent->tableinfo.tab15close4 = (real) value_at_transitions_is_alive_TRUE_FALSE_ind_hci;
				}
			}
			if ( gplTableIndexes[16] >= 0 ) {
				if ( poParent->tableinfo.tab16events == poParent->__events ) {
					poParent->tableinfo.tab16close4 = (real) value_at_transitions_is_alive_TRUE_FALSE_ind_hci;
				}
			}
		}
	}
	return value_at_transitions_is_alive_TRUE_FALSE_ind_hci = value;
}

integer PersonDerivedStates::Set_value_at_transitions_is_alive_TRUE_FALSE_is_stunted( integer value )
{

	if ( value_at_transitions_is_alive_TRUE_FALSE_is_stunted != value ) {

integer		OldValue;
		OldValue = value_at_transitions_is_alive_TRUE_FALSE_is_stunted;
		value_at_transitions_is_alive_TRUE_FALSE_is_stunted = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[15] >= 0 ) {
				if ( poParent->tableinfo.tab15events == poParent->__events ) {
					poParent->tableinfo.tab15close5 = (real) value_at_transitions_is_alive_TRUE_FALSE_is_stunted;
				}
			}
			if ( gplTableIndexes[16] >= 0 ) {
				if ( poParent->tableinfo.tab16events == poParent->__events ) {
					poParent->tableinfo.tab16close5 = (real) value_at_transitions_is_alive_TRUE_FALSE_is_stunted;
				}
			}
		}
	}
	return value_at_transitions_is_alive_TRUE_FALSE_is_stunted = value;
}

real PersonDerivedStates::Set_value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling( real value )
{

	if ( value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling != value ) {

real		OldValue;
		OldValue = value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling;
		value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[15] >= 0 ) {
				if ( poParent->tableinfo.tab15events == poParent->__events ) {
					poParent->tableinfo.tab15close6 = (real) value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling;
				}
			}
			if ( gplTableIndexes[16] >= 0 ) {
				if ( poParent->tableinfo.tab16events == poParent->__events ) {
					poParent->tableinfo.tab16close6 = (real) value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling;
				}
			}
		}
	}
	return value_at_transitions_is_alive_TRUE_FALSE_quality_of_schooling = value;
}

integer PersonDerivedStates::Set_value_at_transitions_is_alive_TRUE_FALSE_survived_early_years( integer value )
{

	if ( value_at_transitions_is_alive_TRUE_FALSE_survived_early_years != value ) {

integer		OldValue;
		OldValue = value_at_transitions_is_alive_TRUE_FALSE_survived_early_years;
		value_at_transitions_is_alive_TRUE_FALSE_survived_early_years = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[15] >= 0 ) {
				if ( poParent->tableinfo.tab15events == poParent->__events ) {
					poParent->tableinfo.tab15close7 = (real) value_at_transitions_is_alive_TRUE_FALSE_survived_early_years;
				}
			}
			if ( gplTableIndexes[16] >= 0 ) {
				if ( poParent->tableinfo.tab16events == poParent->__events ) {
					poParent->tableinfo.tab16close7 = (real) value_at_transitions_is_alive_TRUE_FALSE_survived_early_years;
				}
			}
		}
	}
	return value_at_transitions_is_alive_TRUE_FALSE_survived_early_years = value;
}

integer PersonDerivedStates::Set_value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling( integer value )
{

	if ( value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling != value ) {

integer		OldValue;
		OldValue = value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling;
		value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling = value;
		if ( poParent->bUpdates ) {
			if ( gplTableIndexes[15] >= 0 ) {
				if ( poParent->tableinfo.tab15events == poParent->__events ) {
					poParent->tableinfo.tab15close8 = (real) value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling;
				}
			}
			if ( gplTableIndexes[16] >= 0 ) {
				if ( poParent->tableinfo.tab16events == poParent->__events ) {
					poParent->tableinfo.tab16close8 = (real) value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling;
				}
			}
		}
	}
	return value_at_transitions_is_alive_TRUE_FALSE_years_of_schooling = value;
}

void *StartpopValues::operator new( size_t count )
{
	StartpopValues *poActor;
	if ( gpoEmptyStartpopValues != NULL ) {
		poActor = gpoEmptyStartpopValues;
		gpoEmptyStartpopValues = (StartpopValues *) gpoEmptyStartpopValues->poNext;
	}
	else {
		poActor = ::new StartpopValues( TRUE ) ;
	}
	poActor->__finished = false;
	return poActor;
}

void StartpopValues::Set_lPersonStartpop( Person *value )
{
	Person *poOldLink;

	if ( lPersonStartpop != value ) {
		if (lPersonStartpop != NULL)
		{
			lPersonStartpop->UpdateTime();
		}
		if (value != NULL)
		{
			value->UpdateTime();
		}
		poOldLink = lPersonStartpop;
		lPersonStartpop = value;
		UpdateLinkedStates();
		if ( poOldLink != NULL && poOldLink->lStartValues == this) {
			poOldLink->Set_lStartValues( NULL );
		}
		if ( value != NULL && value->lStartValues != this) {
			value->Set_lStartValues( this );
		}
	}
}

void StartpopValues::InitActor()
{
	poDerivedStates = new StartpopValuesDerivedStates(); 
	poDerivedStates->poParent = this;
	poStateFunctions = new StartpopValuesStateFunctions();
	poStateFunctions->poParent = this;
	poStateFunctions->poDerivedStates = poDerivedStates;
	poDerivedStates->poStateFunctions = poStateFunctions;
	bUpdates = FALSE;
	nActorNumber = 4;
	InitializeStates();
	lPersonStartpop = NULL;
}

void StartpopValues::DeleteActor()
{
	delete poDerivedStates;
	delete poStateFunctions;
}

#pragma optimize( "", off )

void StartpopValues::InitializeStates()
{
	// initialization of states
	actor_id = (long) 0;
	actor_subsample_weight = (double) 0;
	actor_weight = (double) 0;
	age = (TIME) 0;
	case_id = (long) 0;
	case_seed = (double) 0;
#line 22 "../code/ActorStartpopValues.mpp"
	is_activated = FALSE != 0;
#line 14075 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 97 "../code/FamilyFirstUnion.mpp"
	is_set_first_union_clock = FALSE != 0;
#line 14078 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 60 "../code/FertilityGeneral.mpp"
	is_set_last_birth_clock = FALSE != 0;
#line 14081 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 66 "../code/MigrationGeneral.mpp"
	is_set_residential_move_clock = FALSE != 0;
#line 14084 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 70 "../code/MigrationGeneral.mpp"
	scheduled_time_first_move = (TIME) TIME_INFINITE;
#line 14087 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 71 "../code/MigrationGeneral.mpp"
	scheduled_time_last_move = (TIME) TIME_INFINITE;
#line 14090 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	time = (TIME) gdThreadCurrentTime;
#line 100 "../code/FamilyFirstUnion.mpp"
	time_first_union_imputation = (TIME) TIME_INFINITE;
#line 14094 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
#line 63 "../code/FertilityGeneral.mpp"
	time_last_birth_imputation = (TIME) TIME_INFINITE;
#line 14097 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	is_removeable = FALSE;
	lPersonStartpop_is_mortal = FALSE;
	__finished = FALSE;
	__time = (TIME) gdThreadCurrentTime;
	__events = 0;
}

#pragma optimize( "", on )

void StartpopValues::UpdateLinkedStates()
{
	if ( lPersonStartpop != NULL ) {
		poStateFunctions->Set_lPersonStartpop_is_mortal( lPersonStartpop->is_mortal );
	}
	else {
		poStateFunctions->Set_lPersonStartpop_is_mortal( logical( 0 ) );
	}
}

void StartpopValues::StartSpecial1()
{
	nSubSample = GetSubSample();
	InitializeStates();
	case_seed = GetCaseSeed();
	case_id = GetCaseID();
	lActorId = actor_id = GetObjectID();
	GetCaseWeight( &actor_weight, &actor_subsample_weight );
	lPersonStartpop = NULL;

	// initialization of events
	poFirstMoveEvent = gpoEventQueue->NewEvent( this, 14, 0 );
	poFirstUnionImputationEvent = gpoEventQueue->NewEvent( this, 16, 0 );
	poLastBirthImputationEvent = gpoEventQueue->NewEvent( this, 17, 0 );
	poLastMoveEvent = gpoEventQueue->NewEvent( this, 18, 0 );
	poSetFirstUnionClock = gpoEventQueue->NewEvent( this, 25, 0 );
	poSetLastBirthClock = gpoEventQueue->NewEvent( this, 26, 0 );
	poSetResidentialMoveClock = gpoEventQueue->NewEvent( this, 28, 0 );
	poSuicideEvent = gpoEventQueue->NewEvent( this, 31, 0 );
	//Initalize derived states
#line 28 "../code/ActorStartpopValues.mpp"
	poStateFunctions->Set_is_removeable( logical( ( lPersonStartpop_is_mortal ) ) );
#line 14139 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	CountActor(4);
}

#pragma optimize( "", off )

void StartpopValues::StartSpecial2()
{
	// initialization of expression and history states
	poStateFunctions->Set_time(TIME( CoarsenMantissa( time ) ), true);
	poStateFunctions->Set_age(TIME( CoarsenMantissa( age ) ), true);
	gdThreadCurrentTime = time;
	__time = time;
	UpdateLinkedStates();
	poPrev = NULL;
	poNext = gpoEventQueue->poHeadActor;
	if ( gpoEventQueue->poHeadActor != NULL ) {
		gpoEventQueue->poHeadActor->poPrev = this;
	}
	gpoEventQueue->poHeadActor = this;
	bUpdates = TRUE;
}

void StartpopValues::FinishSpecial()
{
	CountFinishedActor(4);
	bUpdates = FALSE;
	gpoEventQueue->DestroyEvent( poFirstMoveEvent );
	gpoEventQueue->DestroyEvent( poFirstUnionImputationEvent );
	gpoEventQueue->DestroyEvent( poLastBirthImputationEvent );
	gpoEventQueue->DestroyEvent( poLastMoveEvent );
	gpoEventQueue->DestroyEvent( poSetFirstUnionClock );
	gpoEventQueue->DestroyEvent( poSetLastBirthClock );
	gpoEventQueue->DestroyEvent( poSetResidentialMoveClock );
	gpoEventQueue->DestroyEvent( poSuicideEvent );
	Set_lPersonStartpop( NULL );

	// Remove actor from master list of active actors
	if ( poPrev != NULL )
	{
		poPrev->poNext = poNext;
	}
	else
	{
		if (gpoEventQueue->poHeadActor != NULL && gpoEventQueue->poHeadActor == this)
		{
			gpoEventQueue->poHeadActor = poNext;
		}
	}
	if ( poNext != NULL )
	{
		poNext->poPrev = poPrev;
	}
	poPrev = NULL;

	// Add actor to head of empty list of actors of this type
	poNext = gpoEmptyStartpopValues;
	gpoEmptyStartpopValues = this;
}

#pragma optimize( "", on )

void StartpopValues::WaitUntil( double target_time )
{
	//don't increase gdEventsForTabulation if special events because they are all treated as one event for tabulation purposes
	if (!gpoEventQueue->m_bSpecialEvents)
	{
		gdEventsForTabulation++;
	}
	__events = gdEventsForTabulation;
	if ( time < (TIME) target_time )
	{
		__time = (TIME) target_time;
		WaitSpecial( (TIME) target_time - time );
		poStateFunctions->Set_time( (TIME) target_time );
	}
}

StartpopValues *StartpopValues::UpdateTime(  )
{
	__events = gdEventsForTabulation;
	if ( time < (TIME) gdThreadCurrentTime ) {
		__time = (TIME) gdThreadCurrentTime;
		WaitSpecial( (TIME) gdThreadCurrentTime - time );
		poStateFunctions->Set_time( (TIME) gdThreadCurrentTime );
	}
	return this;
}

void StartpopValues::WaitUntilThisActor( double target_time )
{
	WaitUntil(target_time);
}

void StartpopValues::WaitSpecial( TIME wait_time )
{

	// increment the age
		poStateFunctions->Set_age( age + wait_time );
}

double StartpopValues::EventTime( int nEventNum, int *pnEventInfo )
{
	double	event_time = double( TIME_INFINITE );

	gdRuntimeTime = (double) time;
	gnRuntimeActorId = actor_id;
	switch( nEventNum ) {
		case 14:
			gszRuntimeEventName = _T("StartpopValues.FirstMoveEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeFirstMoveEvent() );
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("timeFirstMoveEvent"), event_time);
			break;
		case 16:
			gszRuntimeEventName = _T("StartpopValues.FirstUnionImputationEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeFirstUnionImputationEvent() );
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("timeFirstUnionImputationEvent"), event_time);
			break;
		case 17:
			gszRuntimeEventName = _T("StartpopValues.LastBirthImputationEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeLastBirthImputationEvent() );
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("timeLastBirthImputationEvent"), event_time);
			break;
		case 18:
			gszRuntimeEventName = _T("StartpopValues.LastMoveEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeLastMoveEvent() );
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("timeLastMoveEvent"), event_time);
			break;
		case 25:
			gszRuntimeEventName = _T("StartpopValues.SetFirstUnionClock(time)");
			event_time = CoarsenMantissa( (TIME) timeSetFirstUnionClock() );
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("timeSetFirstUnionClock"), event_time);
			break;
		case 26:
			gszRuntimeEventName = _T("StartpopValues.SetLastBirthClock(time)");
			event_time = CoarsenMantissa( (TIME) timeSetLastBirthClock() );
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("timeSetLastBirthClock"), event_time);
			break;
		case 28:
			gszRuntimeEventName = _T("StartpopValues.SetResidentialMoveClock(time)");
			event_time = CoarsenMantissa( (TIME) timeSetResidentialMoveClock() );
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("timeSetResidentialMoveClock"), event_time);
			break;
		case 31:
			gszRuntimeEventName = _T("StartpopValues.SuicideEvent(time)");
			event_time = CoarsenMantissa( (TIME) timeSuicideEvent() );
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("timeSuicideEvent"), event_time);
			break;
		case -1:
		default:;
	}
	return (TIME) event_time;
}

void StartpopValues::Implement( int nEventNum, int nEventInfo )
{
	gdThreadEvents += 1;
	gdRuntimeTime = (double) time;
	gnRuntimeActorId = actor_id;
	switch( nEventNum )
	{
		case 14:
			CHKSUM(time, 0);
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("StartpopValues.FirstMoveEvent"), time);
			gszRuntimeEventName = _T("StartpopValues.FirstMoveEvent(implement)");
			FirstMoveEvent();
			break;
		case 16:
			CHKSUM(time, 1);
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("StartpopValues.FirstUnionImputationEvent"), time);
			gszRuntimeEventName = _T("StartpopValues.FirstUnionImputationEvent(implement)");
			FirstUnionImputationEvent();
			break;
		case 17:
			CHKSUM(time, 2);
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("StartpopValues.LastBirthImputationEvent"), time);
			gszRuntimeEventName = _T("StartpopValues.LastBirthImputationEvent(implement)");
			LastBirthImputationEvent();
			break;
		case 18:
			CHKSUM(time, 3);
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("StartpopValues.LastMoveEvent"), time);
			gszRuntimeEventName = _T("StartpopValues.LastMoveEvent(implement)");
			LastMoveEvent();
			break;
		case 25:
			CHKSUM(time, 4);
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("StartpopValues.SetFirstUnionClock"), time);
			gszRuntimeEventName = _T("StartpopValues.SetFirstUnionClock(implement)");
			SetFirstUnionClock();
			break;
		case 26:
			CHKSUM(time, 5);
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("StartpopValues.SetLastBirthClock"), time);
			gszRuntimeEventName = _T("StartpopValues.SetLastBirthClock(implement)");
			SetLastBirthClock();
			break;
		case 28:
			CHKSUM(time, 6);
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("StartpopValues.SetResidentialMoveClock"), time);
			gszRuntimeEventName = _T("StartpopValues.SetResidentialMoveClock(implement)");
			SetResidentialMoveClock();
			break;
		case 31:
			CHKSUM(time, 7);
			TRCEVNT(case_seed, _T("StartpopValues"), actor_id, _T("StartpopValues.SuicideEvent"), time);
			gszRuntimeEventName = _T("StartpopValues.SuicideEvent(implement)");
			SuicideEvent();
			break;
		case -1:
		default:;
	}
}

double StartpopValues::Set_actor_subsample_weight( double value )
{
	return poStateFunctions->Set_actor_subsample_weight( value);
}

double StartpopValues::Set_actor_weight( double value )
{
	return poStateFunctions->Set_actor_weight( value);
}

logical StartpopValuesStateFunctions::Set_is_activated( logical value )
{

	value = value != 0;
	if ( poParent->is_activated != value )
	{


		logical	OldValue;
		OldValue = poParent->is_activated;
		poParent->is_activated = value;

		if ( poParent->poSetFirstUnionClock->cRecalc == 0 ) {
			poParent->poSetFirstUnionClock->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poSetFirstUnionClock );
		};
		if ( poParent->poSetLastBirthClock->cRecalc == 0 ) {
			poParent->poSetLastBirthClock->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poSetLastBirthClock );
		};
		if ( poParent->poSetResidentialMoveClock->cRecalc == 0 ) {
			poParent->poSetResidentialMoveClock->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poSetResidentialMoveClock );
		};
	}
	return poParent->is_activated = value;
}

logical StartpopValuesStateFunctions::Set_is_set_first_union_clock( logical value )
{

	value = value != 0;
	if ( poParent->is_set_first_union_clock != value )
	{


		logical	OldValue;
		OldValue = poParent->is_set_first_union_clock;
		poParent->is_set_first_union_clock = value;

		if ( poParent->poSetFirstUnionClock->cRecalc == 0 ) {
			poParent->poSetFirstUnionClock->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poSetFirstUnionClock );
		};
	}
	return poParent->is_set_first_union_clock = value;
}

logical StartpopValuesStateFunctions::Set_is_set_last_birth_clock( logical value )
{

	value = value != 0;
	if ( poParent->is_set_last_birth_clock != value )
	{


		logical	OldValue;
		OldValue = poParent->is_set_last_birth_clock;
		poParent->is_set_last_birth_clock = value;

		if ( poParent->poSetLastBirthClock->cRecalc == 0 ) {
			poParent->poSetLastBirthClock->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poSetLastBirthClock );
		};
	}
	return poParent->is_set_last_birth_clock = value;
}

logical StartpopValuesStateFunctions::Set_is_set_residential_move_clock( logical value )
{

	value = value != 0;
	if ( poParent->is_set_residential_move_clock != value )
	{


		logical	OldValue;
		OldValue = poParent->is_set_residential_move_clock;
		poParent->is_set_residential_move_clock = value;

		if ( poParent->poSetResidentialMoveClock->cRecalc == 0 ) {
			poParent->poSetResidentialMoveClock->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poSetResidentialMoveClock );
		};
	}
	return poParent->is_set_residential_move_clock = value;
}

TIME StartpopValuesStateFunctions::Set_scheduled_time_first_move( TIME value )
{

	if ( poParent->scheduled_time_first_move != value )
	{


		TIME	OldValue;
		OldValue = poParent->scheduled_time_first_move;
		poParent->scheduled_time_first_move = value;

		if ( poParent->poFirstMoveEvent->cRecalc == 0 ) {
			poParent->poFirstMoveEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstMoveEvent );
		};
	}
	return poParent->scheduled_time_first_move = value;
}

TIME StartpopValuesStateFunctions::Set_scheduled_time_last_move( TIME value )
{

	if ( poParent->scheduled_time_last_move != value )
	{


		TIME	OldValue;
		OldValue = poParent->scheduled_time_last_move;
		poParent->scheduled_time_last_move = value;

		if ( poParent->poLastMoveEvent->cRecalc == 0 ) {
			poParent->poLastMoveEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poLastMoveEvent );
		};
	}
	return poParent->scheduled_time_last_move = value;
}

TIME StartpopValuesStateFunctions::Set_time_first_union_imputation( TIME value )
{

	if ( poParent->time_first_union_imputation != value )
	{


		TIME	OldValue;
		OldValue = poParent->time_first_union_imputation;
		poParent->time_first_union_imputation = value;

		if ( poParent->poFirstUnionImputationEvent->cRecalc == 0 ) {
			poParent->poFirstUnionImputationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poFirstUnionImputationEvent );
		};
	}
	return poParent->time_first_union_imputation = value;
}

TIME StartpopValuesStateFunctions::Set_time_last_birth_imputation( TIME value )
{

	if ( poParent->time_last_birth_imputation != value )
	{


		TIME	OldValue;
		OldValue = poParent->time_last_birth_imputation;
		poParent->time_last_birth_imputation = value;

		if ( poParent->poLastBirthImputationEvent->cRecalc == 0 ) {
			poParent->poLastBirthImputationEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poLastBirthImputationEvent );
		};
	}
	return poParent->time_last_birth_imputation = value;
}

logical StartpopValuesStateFunctions::Set_is_removeable( logical value )
{

	value = value != 0;
	if ( poParent->is_removeable != value )
	{


		logical	OldValue;
		OldValue = poParent->is_removeable;
		poParent->is_removeable = value;

		if ( poParent->poSuicideEvent->cRecalc == 0 ) {
			poParent->poSuicideEvent->cRecalc = 1;
			gpoEventQueue->RecalculateEvent( poParent->poSuicideEvent );
		};
	}
	return poParent->is_removeable = value;
}

logical StartpopValuesStateFunctions::Set_lPersonStartpop_is_mortal( logical value )
{

	value = value != 0;
	if ( poParent->lPersonStartpop_is_mortal != value )
	{


		logical	OldValue;
		OldValue = poParent->lPersonStartpop_is_mortal;
		poParent->lPersonStartpop_is_mortal = value;

#line 28 "../code/ActorStartpopValues.mpp"
		Set_is_removeable( ( poParent->lPersonStartpop_is_mortal ) );
#line 14561 "C:/Users/Martin/Dropbox/openM++/models/Dynamis-Pop-3-06-ABC/modgen/src/ACTORS.CPP"
	}
	return poParent->lPersonStartpop_is_mortal = value;
}
Thread Calibrator	*gpoEmptyCalibrator = NULL;

void DeleteAllCalibratorActors()
{
}
Thread Clock	*gpoEmptyClock = NULL;

void DeleteAllClockActors()
{
}
Thread Observation	*gpoEmptyObservation = NULL;

void DeleteAllObservationActors()
{
}
Thread Person	*gpoEmptyPerson = NULL;

void DeleteAllPersonActors()
{
}
Thread StartpopValues	*gpoEmptyStartpopValues = NULL;

void DeleteAllStartpopValuesActors()
{
}

void Calibrator_mlCalibratorToPerson::Add( Person *prPerson )
{
	BOOL	bFound = FALSE;
	int	nIndex;
	int	nNewPosition = -1;
	int	nSize;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		if ( arData[nIndex] == prPerson ) {
			bFound = TRUE;
			break;
		}
		if ( arData[nIndex] == NULL && nNewPosition == -1 ) {
			nNewPosition = nIndex;
		}
	}
	if ( !bFound ) {
		prPerson->UpdateTime();
		if ( nNewPosition != -1 ) {
			arData[nNewPosition] = prPerson;
		}
		else {
			arData.Add( prPerson );
		}
		UpdateReferences();
		prPerson->Set_lCalibrator( prParentCalibrator );
	}
}

void Calibrator_mlCalibratorToPerson::Remove( Person *prPerson )
{
	int	nIndex;
	int	nSize;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		if ( arData[nIndex] == prPerson ) {
			prPerson->UpdateTime();
			arData[nIndex] = NULL;
			UpdateReferences();
			if ( prPerson->lCalibrator == prParentCalibrator ) {
				prPerson->Set_lCalibrator( NULL );
			}
			break;
		}
	}
}

void Calibrator_mlCalibratorToPerson::RemoveAll()
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prPerson->UpdateTime();
			arData[nIndex] = NULL;
			if ( prPerson->lCalibrator == prParentCalibrator ) {
				prPerson->Set_lCalibrator( NULL );
			}
		}
	}
	UpdateReferences();
}

void Calibrator_mlCalibratorToPerson::Finish()
{
	int	nSize = 0;
	nSize = arData.GetSize();
	RemoveAll();
	SetMaxActorLinked(0, 1, nSize);
}

void Calibrator_mlCalibratorToPerson::FinishAll()
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prPerson->UpdateTime();
			prPerson->Finish();
			arData[nIndex] = NULL;
		}
	}
	UpdateReferences();
}

Person *Calibrator_mlCalibratorToPerson::GetNext( int nInitPos, int *pnPos )
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;
	Person	*prFoundPerson = NULL;

	*pnPos = -1;
	nSize = arData.GetSize();
	for ( nIndex = nInitPos; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prFoundPerson = prPerson;
			*pnPos = nIndex;
			break;
		}
	}
	return prFoundPerson;
}

void Calibrator_mlCalibratorToPerson::UpdateReferences()
{
}

void Person_mlBioFatherChildren::Add( Person *prPerson )
{
	BOOL	bFound = FALSE;
	int	nIndex;
	int	nNewPosition = -1;
	int	nSize;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		if ( arData[nIndex] == prPerson ) {
			bFound = TRUE;
			break;
		}
		if ( arData[nIndex] == NULL && nNewPosition == -1 ) {
			nNewPosition = nIndex;
		}
	}
	if ( !bFound ) {
		prPerson->UpdateTime();
		if ( nNewPosition != -1 ) {
			arData[nNewPosition] = prPerson;
		}
		else {
			arData.Add( prPerson );
		}
		UpdateReferences();
		prPerson->Set_lBioFather( prParentPerson );
	}
}

void Person_mlBioFatherChildren::Remove( Person *prPerson )
{
	int	nIndex;
	int	nSize;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		if ( arData[nIndex] == prPerson ) {
			prPerson->UpdateTime();
			arData[nIndex] = NULL;
			UpdateReferences();
			if ( prPerson->lBioFather == prParentPerson ) {
				prPerson->Set_lBioFather( NULL );
			}
			break;
		}
	}
}

void Person_mlBioFatherChildren::RemoveAll()
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prPerson->UpdateTime();
			arData[nIndex] = NULL;
			if ( prPerson->lBioFather == prParentPerson ) {
				prPerson->Set_lBioFather( NULL );
			}
		}
	}
	UpdateReferences();
}

void Person_mlBioFatherChildren::Finish()
{
	int	nSize = 0;
	nSize = arData.GetSize();
	RemoveAll();
	SetMaxActorLinked(3, 7, nSize);
}

void Person_mlBioFatherChildren::FinishAll()
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prPerson->UpdateTime();
			prPerson->Finish();
			arData[nIndex] = NULL;
		}
	}
	UpdateReferences();
}

Person *Person_mlBioFatherChildren::GetNext( int nInitPos, int *pnPos )
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;
	Person	*prFoundPerson = NULL;

	*pnPos = -1;
	nSize = arData.GetSize();
	for ( nIndex = nInitPos; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prFoundPerson = prPerson;
			*pnPos = nIndex;
			break;
		}
	}
	return prFoundPerson;
}

void Person_mlBioFatherChildren::UpdateReferences()
{
}

void Person_mlBioMotherChildren::Add( Person *prPerson )
{
	BOOL	bFound = FALSE;
	int	nIndex;
	int	nNewPosition = -1;
	int	nSize;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		if ( arData[nIndex] == prPerson ) {
			bFound = TRUE;
			break;
		}
		if ( arData[nIndex] == NULL && nNewPosition == -1 ) {
			nNewPosition = nIndex;
		}
	}
	if ( !bFound ) {
		prPerson->UpdateTime();
		if ( nNewPosition != -1 ) {
			arData[nNewPosition] = prPerson;
		}
		else {
			arData.Add( prPerson );
		}
		UpdateReferences();
		prPerson->Set_lBioMother( prParentPerson );
	}
}

void Person_mlBioMotherChildren::Remove( Person *prPerson )
{
	int	nIndex;
	int	nSize;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		if ( arData[nIndex] == prPerson ) {
			prPerson->UpdateTime();
			arData[nIndex] = NULL;
			UpdateReferences();
			if ( prPerson->lBioMother == prParentPerson ) {
				prPerson->Set_lBioMother( NULL );
			}
			break;
		}
	}
}

void Person_mlBioMotherChildren::RemoveAll()
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prPerson->UpdateTime();
			arData[nIndex] = NULL;
			if ( prPerson->lBioMother == prParentPerson ) {
				prPerson->Set_lBioMother( NULL );
			}
		}
	}
	UpdateReferences();
}

void Person_mlBioMotherChildren::Finish()
{
	int	nSize = 0;
	nSize = arData.GetSize();
	RemoveAll();
	SetMaxActorLinked(3, 8, nSize);
}

void Person_mlBioMotherChildren::FinishAll()
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prPerson->UpdateTime();
			prPerson->Finish();
			arData[nIndex] = NULL;
		}
	}
	UpdateReferences();
}

Person *Person_mlBioMotherChildren::GetNext( int nInitPos, int *pnPos )
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;
	Person	*prFoundPerson = NULL;

	*pnPos = -1;
	nSize = arData.GetSize();
	for ( nIndex = nInitPos; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prFoundPerson = prPerson;
			*pnPos = nIndex;
			break;
		}
	}
	return prFoundPerson;
}

void Person_mlBioMotherChildren::UpdateReferences()
{
}

void Person_mlHHFatherChildren::Add( Person *prPerson )
{
	BOOL	bFound = FALSE;
	int	nIndex;
	int	nNewPosition = -1;
	int	nSize;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		if ( arData[nIndex] == prPerson ) {
			bFound = TRUE;
			break;
		}
		if ( arData[nIndex] == NULL && nNewPosition == -1 ) {
			nNewPosition = nIndex;
		}
	}
	if ( !bFound ) {
		prPerson->UpdateTime();
		if ( nNewPosition != -1 ) {
			arData[nNewPosition] = prPerson;
		}
		else {
			arData.Add( prPerson );
		}
		UpdateReferences();
		prPerson->Set_lHHFather( prParentPerson );
	}
}

void Person_mlHHFatherChildren::Remove( Person *prPerson )
{
	int	nIndex;
	int	nSize;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		if ( arData[nIndex] == prPerson ) {
			prPerson->UpdateTime();
			arData[nIndex] = NULL;
			UpdateReferences();
			if ( prPerson->lHHFather == prParentPerson ) {
				prPerson->Set_lHHFather( NULL );
			}
			break;
		}
	}
}

void Person_mlHHFatherChildren::RemoveAll()
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prPerson->UpdateTime();
			arData[nIndex] = NULL;
			if ( prPerson->lHHFather == prParentPerson ) {
				prPerson->Set_lHHFather( NULL );
			}
		}
	}
	UpdateReferences();
}

void Person_mlHHFatherChildren::Finish()
{
	int	nSize = 0;
	nSize = arData.GetSize();
	RemoveAll();
	SetMaxActorLinked(3, 9, nSize);
}

void Person_mlHHFatherChildren::FinishAll()
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prPerson->UpdateTime();
			prPerson->Finish();
			arData[nIndex] = NULL;
		}
	}
	UpdateReferences();
}

Person *Person_mlHHFatherChildren::GetNext( int nInitPos, int *pnPos )
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;
	Person	*prFoundPerson = NULL;

	*pnPos = -1;
	nSize = arData.GetSize();
	for ( nIndex = nInitPos; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prFoundPerson = prPerson;
			*pnPos = nIndex;
			break;
		}
	}
	return prFoundPerson;
}

counter Person_mlHHFatherChildren::Count()
{
	int	nIndex;
	int	nResult = 0;
	int	nSize;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		if(arData[nIndex] != NULL)
		{
			nResult++;
		}
	}
	return nResult;
}

integer Person_mlHHFatherChildren::SumOver_is_alive()
{
	Person	*prPerson;
	int	nIndex;
	int	nSize;
	integer	result;

	result = (integer) 0;
	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			result += prPerson->is_alive;
		}
	}
	return result;
}

void Person_mlHHFatherChildren::UpdateReferences()
{
	prParentPerson->poStateFunctions->Set_count_mlHHFatherChildren( Count() );
	prParentPerson->poStateFunctions->Set_sum_over_mlHHFatherChildren_is_alive( SumOver_is_alive() );
}

void Person_mlHHMotherChildren::Add( Person *prPerson )
{
	BOOL	bFound = FALSE;
	int	nIndex;
	int	nNewPosition = -1;
	int	nSize;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		if ( arData[nIndex] == prPerson ) {
			bFound = TRUE;
			break;
		}
		if ( arData[nIndex] == NULL && nNewPosition == -1 ) {
			nNewPosition = nIndex;
		}
	}
	if ( !bFound ) {
		prPerson->UpdateTime();
		if ( nNewPosition != -1 ) {
			arData[nNewPosition] = prPerson;
		}
		else {
			arData.Add( prPerson );
		}
		UpdateReferences();
		prPerson->Set_lHHMother( prParentPerson );
	}
}

void Person_mlHHMotherChildren::Remove( Person *prPerson )
{
	int	nIndex;
	int	nSize;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		if ( arData[nIndex] == prPerson ) {
			prPerson->UpdateTime();
			arData[nIndex] = NULL;
			UpdateReferences();
			if ( prPerson->lHHMother == prParentPerson ) {
				prPerson->Set_lHHMother( NULL );
			}
			break;
		}
	}
}

void Person_mlHHMotherChildren::RemoveAll()
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prPerson->UpdateTime();
			arData[nIndex] = NULL;
			if ( prPerson->lHHMother == prParentPerson ) {
				prPerson->Set_lHHMother( NULL );
			}
		}
	}
	UpdateReferences();
}

void Person_mlHHMotherChildren::Finish()
{
	int	nSize = 0;
	nSize = arData.GetSize();
	RemoveAll();
	SetMaxActorLinked(3, 10, nSize);
}

void Person_mlHHMotherChildren::FinishAll()
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prPerson->UpdateTime();
			prPerson->Finish();
			arData[nIndex] = NULL;
		}
	}
	UpdateReferences();
}

Person *Person_mlHHMotherChildren::GetNext( int nInitPos, int *pnPos )
{
	int	nIndex;
	int	nSize;
	Person	*prPerson;
	Person	*prFoundPerson = NULL;

	*pnPos = -1;
	nSize = arData.GetSize();
	for ( nIndex = nInitPos; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			prFoundPerson = prPerson;
			*pnPos = nIndex;
			break;
		}
	}
	return prFoundPerson;
}

counter Person_mlHHMotherChildren::Count()
{
	int	nIndex;
	int	nResult = 0;
	int	nSize;

	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		if(arData[nIndex] != NULL)
		{
			nResult++;
		}
	}
	return nResult;
}

AGE_RANGE Person_mlHHMotherChildren::MinOver_integer_age( Person *prSpecialPerson, AGE_RANGE value )
{
	BOOL	bFound = FALSE;
	Person	*prPerson;
	int	nIndex;
	int	nSize;
	AGE_RANGE	current_value;
	AGE_RANGE	result;

	result = (AGE_RANGE) r_min_AGE_RANGE;
	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			current_value = ( prPerson != prSpecialPerson ) ? prPerson->integer_age : value;
			result = bFound ? min( result, current_value ) : current_value;
			bFound = TRUE;
		}
	}
	return result;
}

integer Person_mlHHMotherChildren::SumOver_is_alive()
{
	Person	*prPerson;
	int	nIndex;
	int	nSize;
	integer	result;

	result = (integer) 0;
	nSize = arData.GetSize();
	for ( nIndex = 0; nIndex < nSize; nIndex++ ) {
		prPerson = arData[nIndex];
		if ( prPerson != NULL ) {
			result += prPerson->is_alive;
		}
	}
	return result;
}

void Person_mlHHMotherChildren::UpdateReferences()
{
	prParentPerson->poStateFunctions->Set_count_mlHHMotherChildren( Count() );
	prParentPerson->poStateFunctions->Set_min_over_mlHHMotherChildren_integer_age( MinOver_integer_age( NULL, 0 ) );
	prParentPerson->poStateFunctions->Set_sum_over_mlHHMotherChildren_is_alive( SumOver_is_alive() );
}

int	asAllFertilePersonsForFertilityAlignmentActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asAllFertilePersonsForFertilityAlignmentActorSetTree	*asAllFertilePersonsForFertilityAlignment[2];

int	asAllFertilePersonsForFertilityAlignmentByAgeActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asAllFertilePersonsForFertilityAlignmentByAgeActorSetTree	*asAllFertilePersonsForFertilityAlignmentByAge[rc_size_FERTILE_AGE_RANGE][2];

int	asAllPersonActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asAllPersonActorSetTree	*asAllPerson;

int	asAvailableMaleActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asAvailableMaleActorSetTree	*asAvailableMale[rc_size_MALE_SPOUSE_AGE][rc_size_GEO];

int	asAvailableMaleByTypeActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asAvailableMaleByTypeActorSetTree	*asAvailableMaleByType[rc_size_MALE_SPOUSE_AGE][rc_size_SPOUSE_GROUP][rc_size_GEO];

int	asAvailableMaleByTypeImmiScratchActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asAvailableMaleByTypeImmiScratchActorSetTree	*asAvailableMaleByTypeImmiScratch[rc_size_MALE_SPOUSE_AGE][rc_size_SPOUSE_GROUP][rc_size_GEO][rc_size_SIM_YEAR_RANGE];

int	asAvailableMaleImmiScratchActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asAvailableMaleImmiScratchActorSetTree	*asAvailableMaleImmiScratch[rc_size_MALE_SPOUSE_AGE][rc_size_GEO][rc_size_SIM_YEAR_RANGE];

int	asCalibratorActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Calibrator	*poCalibrator1 = NULL;
	Calibrator	*poCalibrator2 = NULL;

	poCalibrator1 = (Calibrator *) poActor1;
	poCalibrator2 = (Calibrator *) poActor2;
	return (poCalibrator1->actor_id == poCalibrator2->actor_id ? 0 :
			poCalibrator1->actor_id < poCalibrator2->actor_id ? -1 : 1);
}
Thread asCalibratorActorSetTree	*asCalibrator;

int	asClockActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Clock	*poClock1 = NULL;
	Clock	*poClock2 = NULL;

	poClock1 = (Clock *) poActor1;
	poClock2 = (Clock *) poActor2;
	return (poClock1->actor_id == poClock2->actor_id ? 0 :
			poClock1->actor_id < poClock2->actor_id ? -1 : 1);
}
Thread asClockActorSetTree	*asClock;

int	asFemaleInUnionByAgeAndPartnerAgeActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asFemaleInUnionByAgeAndPartnerAgeActorSetTree	*asFemaleInUnionByAgeAndPartnerAge[rc_size_FEMALE_SPOUSE_AGE][rc_size_MALE_SPOUSE_AGE];

int	asFemaleInUnionByAgeAndPartnerAgeImmiScratchActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asFemaleInUnionByAgeAndPartnerAgeImmiScratchActorSetTree	*asFemaleInUnionByAgeAndPartnerAgeImmiScratch[rc_size_FEMALE_SPOUSE_AGE][rc_size_MALE_SPOUSE_AGE][rc_size_SIM_YEAR_RANGE];

int	asObservationAllActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Observation	*poObservation1 = NULL;
	Observation	*poObservation2 = NULL;

	poObservation1 = (Observation *) poActor1;
	poObservation2 = (Observation *) poActor2;
	return (poObservation1->actor_id == poObservation2->actor_id ? 0 :
			poObservation1->actor_id < poObservation2->actor_id ? -1 : 1);
}
Thread asObservationAllActorSetTree	*asObservationAll[rc_size_POP_POOL];

int	asObservationByFamActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Observation	*poObservation1 = NULL;
	Observation	*poObservation2 = NULL;

	poObservation1 = (Observation *) poActor1;
	poObservation2 = (Observation *) poActor2;
	return (poObservation1->actor_id == poObservation2->actor_id ? 0 :
			poObservation1->actor_id < poObservation2->actor_id ? -1 : 1);
}
Thread asObservationByFamActorSetTree	*asObservationByFam[rc_size_FAM_ID];

int	asObservationByFamOldestActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Observation	*poObservation1 = NULL;
	Observation	*poObservation2 = NULL;

	poObservation1 = (Observation *) poActor1;
	poObservation2 = (Observation *) poActor2;
	return (poObservation1->actor_id == poObservation2->actor_id ? 0 :
			poObservation1->actor_id < poObservation2->actor_id ? -1 : 1);
}
Thread asObservationByFamOldestActorSetTree	*asObservationByFamOldest[rc_size_FAM_ID];

int	asObservationHeadsActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Observation	*poObservation1 = NULL;
	Observation	*poObservation2 = NULL;

	poObservation1 = (Observation *) poActor1;
	poObservation2 = (Observation *) poActor2;
	return (poObservation1->actor_id == poObservation2->actor_id ? 0 :
			poObservation1->actor_id < poObservation2->actor_id ? -1 : 1);
}
Thread asObservationHeadsActorSetTree	*asObservationHeads[rc_size_POP_POOL];

int	asObservationNonHeadsActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Observation	*poObservation1 = NULL;
	Observation	*poObservation2 = NULL;

	poObservation1 = (Observation *) poActor1;
	poObservation2 = (Observation *) poActor2;
	return (poObservation1->actor_id == poObservation2->actor_id ? 0 :
			poObservation1->actor_id < poObservation2->actor_id ? -1 : 1);
}
Thread asObservationNonHeadsActorSetTree	*asObservationNonHeads[rc_size_FAM_ID][rc_size_POP_POOL];

int	asObservationsActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Observation	*poObservation1 = NULL;
	Observation	*poObservation2 = NULL;

	poObservation1 = (Observation *) poActor1;
	poObservation2 = (Observation *) poActor2;
	return (poObservation1->actor_id == poObservation2->actor_id ? 0 :
			poObservation1->actor_id < poObservation2->actor_id ? -1 : 1);
}
Thread asObservationsActorSetTree	*asObservations;

int	asPotentialImmigrantMothersActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asPotentialImmigrantMothersActorSetTree	*asPotentialImmigrantMothers[rc_size_AGE_RANGE][rc_size_SIM_YEAR_RANGE][rc_size_GEO];

int	asPotentialSchoolOneStudentActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asPotentialSchoolOneStudentActorSetTree	*asPotentialSchoolOneStudent;

int	asPotentialSchoolOneStudentToProcessActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asPotentialSchoolOneStudentToProcessActorSetTree	*asPotentialSchoolOneStudentToProcess;

int	asPotentialSchoolTwoStudentActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asPotentialSchoolTwoStudentActorSetTree	*asPotentialSchoolTwoStudent;

int	asPotentialSchoolTwoStudentToProcessActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asPotentialSchoolTwoStudentToProcessActorSetTree	*asPotentialSchoolTwoStudentToProcess;

int	asResidentBabiesActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asResidentBabiesActorSetTree	*asResidentBabies;

int	asResidentsAge0SexActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asResidentsAge0SexActorSetTree	*asResidentsAge0Sex[rc_size_SEX];

int	asResidentsAge0SexGeoActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asResidentsAge0SexGeoActorSetTree	*asResidentsAge0SexGeo[rc_size_SEX][rc_size_GEO];

int	asSimBornAge0ActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asSimBornAge0ActorSetTree	*asSimBornAge0[rc_size_SEX][rc_size_GEO][rc_size_EDUC1_GROUP];

int	asSimulatedObservationHeadsActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Observation	*poObservation1 = NULL;
	Observation	*poObservation2 = NULL;

	poObservation1 = (Observation *) poActor1;
	poObservation2 = (Observation *) poActor2;
	return (poObservation1->actor_id == poObservation2->actor_id ? 0 :
			poObservation1->actor_id < poObservation2->actor_id ? -1 : 1);
}
Thread asSimulatedObservationHeadsActorSetTree	*asSimulatedObservationHeads[rc_size_POP_POOL];

int	asWantToMoveActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asWantToMoveActorSetTree	*asWantToMove[rc_size_GEO][rc_size_GEO][rc_size_SEX][rc_size_AGE_MIG];

int	asWomenNoChildrenActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asWomenNoChildrenActorSetTree	*asWomenNoChildren[rc_size_UNION1_GROUP][rc_size_FEMALE_SPOUSE_AGE][2];

int	asWomenWithChildrenActorSetTree::ActorCompare(ActorClass *poActor1, ActorClass *poActor2)
	{
	Person	*poPerson1 = NULL;
	Person	*poPerson2 = NULL;

	poPerson1 = (Person *) poActor1;
	poPerson2 = (Person *) poActor2;
	return (poPerson1->actor_id == poPerson2->actor_id ? 0 :
			poPerson1->actor_id < poPerson2->actor_id ? -1 : 1);
}
Thread asWomenWithChildrenActorSetTree	*asWomenWithChildren[rc_size_UNION1_GROUP][rc_size_CHILD_AGEGR][rc_size_MOTH_AGEGR][2];

void PreSimulation() {
	PreSimulation0();
	PreSimulation1();
	PreSimulation2();
	PreSimulation3();
	PreSimulation4();
}

void PostSimulation() {
}

void UserTables() {
}

BOOL ValidateParameters( SCENARIO_EVENT eEvent ) {
	BOOL bProceed = TRUE;

	 return bProceed;
}

void InitActorSets()
{
	int	nIndex0;
	int	nIndex1;
	int	nIndex2;
	int	nIndex3;

	for (nIndex0 = 0; nIndex0 < 2; nIndex0++)
	{
		asAllFertilePersonsForFertilityAlignment[nIndex0] = new asAllFertilePersonsForFertilityAlignmentActorSetTree();
	}
	for (nIndex0 = 0; nIndex0 < 38; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 2; nIndex1++)
		{
			asAllFertilePersonsForFertilityAlignmentByAge[nIndex0][nIndex1] = new asAllFertilePersonsForFertilityAlignmentByAgeActorSetTree();
		}
	}
	asAllPerson = new asAllPersonActorSetTree();
	for (nIndex0 = 0; nIndex0 < 86; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 13; nIndex1++)
		{
			asAvailableMale[nIndex0][nIndex1] = new asAvailableMaleActorSetTree();
		}
	}
	for (nIndex0 = 0; nIndex0 < 86; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 3; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 13; nIndex2++)
			{
				asAvailableMaleByType[nIndex0][nIndex1][nIndex2] = new asAvailableMaleByTypeActorSetTree();
			}
		}
	}
	for (nIndex0 = 0; nIndex0 < 86; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 3; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 13; nIndex2++)
			{
				for (nIndex3 = 0; nIndex3 < 151; nIndex3++)
				{
					asAvailableMaleByTypeImmiScratch[nIndex0][nIndex1][nIndex2][nIndex3] = new asAvailableMaleByTypeImmiScratchActorSetTree();
				}
			}
		}
	}
	for (nIndex0 = 0; nIndex0 < 86; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 13; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 151; nIndex2++)
			{
				asAvailableMaleImmiScratch[nIndex0][nIndex1][nIndex2] = new asAvailableMaleImmiScratchActorSetTree();
			}
		}
	}
	asCalibrator = new asCalibratorActorSetTree();
	asClock = new asClockActorSetTree();
	for (nIndex0 = 0; nIndex0 < 51; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 86; nIndex1++)
		{
			asFemaleInUnionByAgeAndPartnerAge[nIndex0][nIndex1] = new asFemaleInUnionByAgeAndPartnerAgeActorSetTree();
		}
	}
	for (nIndex0 = 0; nIndex0 < 51; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 86; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 151; nIndex2++)
			{
				asFemaleInUnionByAgeAndPartnerAgeImmiScratch[nIndex0][nIndex1][nIndex2] = new asFemaleInUnionByAgeAndPartnerAgeImmiScratchActorSetTree();
			}
		}
	}
	for (nIndex0 = 0; nIndex0 < 4; nIndex0++)
	{
		asObservationAll[nIndex0] = new asObservationAllActorSetTree();
	}
	for (nIndex0 = 0; nIndex0 < 2000001; nIndex0++)
	{
		asObservationByFam[nIndex0] = new asObservationByFamActorSetTree();
	}
	for (nIndex0 = 0; nIndex0 < 2000001; nIndex0++)
	{
		asObservationByFamOldest[nIndex0] = new asObservationByFamOldestActorSetTree();
	}
	for (nIndex0 = 0; nIndex0 < 4; nIndex0++)
	{
		asObservationHeads[nIndex0] = new asObservationHeadsActorSetTree();
	}
	for (nIndex0 = 0; nIndex0 < 2000001; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 4; nIndex1++)
		{
			asObservationNonHeads[nIndex0][nIndex1] = new asObservationNonHeadsActorSetTree();
		}
	}
	asObservations = new asObservationsActorSetTree();
	for (nIndex0 = 0; nIndex0 < 101; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 151; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 13; nIndex2++)
			{
				asPotentialImmigrantMothers[nIndex0][nIndex1][nIndex2] = new asPotentialImmigrantMothersActorSetTree();
			}
		}
	}
	asPotentialSchoolOneStudent = new asPotentialSchoolOneStudentActorSetTree();
	asPotentialSchoolOneStudentToProcess = new asPotentialSchoolOneStudentToProcessActorSetTree();
	asPotentialSchoolTwoStudent = new asPotentialSchoolTwoStudentActorSetTree();
	asPotentialSchoolTwoStudentToProcess = new asPotentialSchoolTwoStudentToProcessActorSetTree();
	asResidentBabies = new asResidentBabiesActorSetTree();
	for (nIndex0 = 0; nIndex0 < 2; nIndex0++)
	{
		asResidentsAge0Sex[nIndex0] = new asResidentsAge0SexActorSetTree();
	}
	for (nIndex0 = 0; nIndex0 < 2; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 13; nIndex1++)
		{
			asResidentsAge0SexGeo[nIndex0][nIndex1] = new asResidentsAge0SexGeoActorSetTree();
		}
	}
	for (nIndex0 = 0; nIndex0 < 2; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 13; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 3; nIndex2++)
			{
				asSimBornAge0[nIndex0][nIndex1][nIndex2] = new asSimBornAge0ActorSetTree();
			}
		}
	}
	for (nIndex0 = 0; nIndex0 < 4; nIndex0++)
	{
		asSimulatedObservationHeads[nIndex0] = new asSimulatedObservationHeadsActorSetTree();
	}
	for (nIndex0 = 0; nIndex0 < 13; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 13; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 2; nIndex2++)
			{
				for (nIndex3 = 0; nIndex3 < 13; nIndex3++)
				{
					asWantToMove[nIndex0][nIndex1][nIndex2][nIndex3] = new asWantToMoveActorSetTree();
				}
			}
		}
	}
	for (nIndex0 = 0; nIndex0 < 3; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 51; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 2; nIndex2++)
			{
				asWomenNoChildren[nIndex0][nIndex1][nIndex2] = new asWomenNoChildrenActorSetTree();
			}
		}
	}
	for (nIndex0 = 0; nIndex0 < 3; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 7; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 6; nIndex2++)
			{
				for (nIndex3 = 0; nIndex3 < 2; nIndex3++)
				{
					asWomenWithChildren[nIndex0][nIndex1][nIndex2][nIndex3] = new asWomenWithChildrenActorSetTree();
				}
			}
		}
	}
}

void	DeleteActorSets()
{
	int	nIndex0;
	int	nIndex1;
	int	nIndex2;
	int	nIndex3;

	for (nIndex0 = 0; nIndex0 < 2; nIndex0++)
	{
		delete asAllFertilePersonsForFertilityAlignment[nIndex0];
		asAllFertilePersonsForFertilityAlignment[nIndex0] = NULL;
	}
	for (nIndex0 = 0; nIndex0 < 38; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 2; nIndex1++)
		{
			delete asAllFertilePersonsForFertilityAlignmentByAge[nIndex0][nIndex1];
			asAllFertilePersonsForFertilityAlignmentByAge[nIndex0][nIndex1] = NULL;
		}
	}
	delete asAllPerson;
	asAllPerson = NULL;
	for (nIndex0 = 0; nIndex0 < 86; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 13; nIndex1++)
		{
			delete asAvailableMale[nIndex0][nIndex1];
			asAvailableMale[nIndex0][nIndex1] = NULL;
		}
	}
	for (nIndex0 = 0; nIndex0 < 86; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 3; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 13; nIndex2++)
			{
				delete asAvailableMaleByType[nIndex0][nIndex1][nIndex2];
				asAvailableMaleByType[nIndex0][nIndex1][nIndex2] = NULL;
			}
		}
	}
	for (nIndex0 = 0; nIndex0 < 86; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 3; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 13; nIndex2++)
			{
				for (nIndex3 = 0; nIndex3 < 151; nIndex3++)
				{
					delete asAvailableMaleByTypeImmiScratch[nIndex0][nIndex1][nIndex2][nIndex3];
					asAvailableMaleByTypeImmiScratch[nIndex0][nIndex1][nIndex2][nIndex3] = NULL;
				}
			}
		}
	}
	for (nIndex0 = 0; nIndex0 < 86; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 13; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 151; nIndex2++)
			{
				delete asAvailableMaleImmiScratch[nIndex0][nIndex1][nIndex2];
				asAvailableMaleImmiScratch[nIndex0][nIndex1][nIndex2] = NULL;
			}
		}
	}
	delete asCalibrator;
	asCalibrator = NULL;
	delete asClock;
	asClock = NULL;
	for (nIndex0 = 0; nIndex0 < 51; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 86; nIndex1++)
		{
			delete asFemaleInUnionByAgeAndPartnerAge[nIndex0][nIndex1];
			asFemaleInUnionByAgeAndPartnerAge[nIndex0][nIndex1] = NULL;
		}
	}
	for (nIndex0 = 0; nIndex0 < 51; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 86; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 151; nIndex2++)
			{
				delete asFemaleInUnionByAgeAndPartnerAgeImmiScratch[nIndex0][nIndex1][nIndex2];
				asFemaleInUnionByAgeAndPartnerAgeImmiScratch[nIndex0][nIndex1][nIndex2] = NULL;
			}
		}
	}
	for (nIndex0 = 0; nIndex0 < 4; nIndex0++)
	{
		delete asObservationAll[nIndex0];
		asObservationAll[nIndex0] = NULL;
	}
	for (nIndex0 = 0; nIndex0 < 2000001; nIndex0++)
	{
		delete asObservationByFam[nIndex0];
		asObservationByFam[nIndex0] = NULL;
	}
	for (nIndex0 = 0; nIndex0 < 2000001; nIndex0++)
	{
		delete asObservationByFamOldest[nIndex0];
		asObservationByFamOldest[nIndex0] = NULL;
	}
	for (nIndex0 = 0; nIndex0 < 4; nIndex0++)
	{
		delete asObservationHeads[nIndex0];
		asObservationHeads[nIndex0] = NULL;
	}
	for (nIndex0 = 0; nIndex0 < 2000001; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 4; nIndex1++)
		{
			delete asObservationNonHeads[nIndex0][nIndex1];
			asObservationNonHeads[nIndex0][nIndex1] = NULL;
		}
	}
	delete asObservations;
	asObservations = NULL;
	for (nIndex0 = 0; nIndex0 < 101; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 151; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 13; nIndex2++)
			{
				delete asPotentialImmigrantMothers[nIndex0][nIndex1][nIndex2];
				asPotentialImmigrantMothers[nIndex0][nIndex1][nIndex2] = NULL;
			}
		}
	}
	delete asPotentialSchoolOneStudent;
	asPotentialSchoolOneStudent = NULL;
	delete asPotentialSchoolOneStudentToProcess;
	asPotentialSchoolOneStudentToProcess = NULL;
	delete asPotentialSchoolTwoStudent;
	asPotentialSchoolTwoStudent = NULL;
	delete asPotentialSchoolTwoStudentToProcess;
	asPotentialSchoolTwoStudentToProcess = NULL;
	delete asResidentBabies;
	asResidentBabies = NULL;
	for (nIndex0 = 0; nIndex0 < 2; nIndex0++)
	{
		delete asResidentsAge0Sex[nIndex0];
		asResidentsAge0Sex[nIndex0] = NULL;
	}
	for (nIndex0 = 0; nIndex0 < 2; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 13; nIndex1++)
		{
			delete asResidentsAge0SexGeo[nIndex0][nIndex1];
			asResidentsAge0SexGeo[nIndex0][nIndex1] = NULL;
		}
	}
	for (nIndex0 = 0; nIndex0 < 2; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 13; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 3; nIndex2++)
			{
				delete asSimBornAge0[nIndex0][nIndex1][nIndex2];
				asSimBornAge0[nIndex0][nIndex1][nIndex2] = NULL;
			}
		}
	}
	for (nIndex0 = 0; nIndex0 < 4; nIndex0++)
	{
		delete asSimulatedObservationHeads[nIndex0];
		asSimulatedObservationHeads[nIndex0] = NULL;
	}
	for (nIndex0 = 0; nIndex0 < 13; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 13; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 2; nIndex2++)
			{
				for (nIndex3 = 0; nIndex3 < 13; nIndex3++)
				{
					delete asWantToMove[nIndex0][nIndex1][nIndex2][nIndex3];
					asWantToMove[nIndex0][nIndex1][nIndex2][nIndex3] = NULL;
				}
			}
		}
	}
	for (nIndex0 = 0; nIndex0 < 3; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 51; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 2; nIndex2++)
			{
				delete asWomenNoChildren[nIndex0][nIndex1][nIndex2];
				asWomenNoChildren[nIndex0][nIndex1][nIndex2] = NULL;
			}
		}
	}
	for (nIndex0 = 0; nIndex0 < 3; nIndex0++)
	{
		for (nIndex1 = 0; nIndex1 < 7; nIndex1++)
		{
			for (nIndex2 = 0; nIndex2 < 6; nIndex2++)
			{
				for (nIndex3 = 0; nIndex3 < 2; nIndex3++)
				{
					delete asWomenWithChildren[nIndex0][nIndex1][nIndex2][nIndex3];
					asWomenWithChildren[nIndex0][nIndex1][nIndex2][nIndex3] = NULL;
				}
			}
		}
	}
}

} // namespace 
