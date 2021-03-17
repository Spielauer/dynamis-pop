--
-- drop compatibility views for model: Dynamis-Pop-3-05-ABC
-- model digest:   408768e630b6131383a26171a0e0e299
-- script created: 2021-03-16 22:42:10.122
--
-- Dear user:
--   this part of database is optional and NOT used by openM++
--   if you want it for any reason please enjoy else just ignore it
-- Or other words:
--   if you don't know what is this then you don't need it
--

--
-- drop input parameters compatibility views
--
DROP VIEW AgeImmiSearchMother;
DROP VIEW AgeImmigrantsScratch;
DROP VIEW AgeLeavingHome;
DROP VIEW AgeOfImmigrantMother;
DROP VIEW AgeSpecificFertility;
DROP VIEW BackMigrationHazard;
DROP VIEW BirthTrends;
DROP VIEW ChildMortalityBaseRisk;
DROP VIEW ChildMortalityRelativeRisks;
DROP VIEW ChildMortalityTrend;
DROP VIEW ChildVaccinationOdds;
DROP VIEW Educ1FirstCohortRefinedModel;
DROP VIEW Educ1GradOdds;
DROP VIEW Educ1Infrastructure;
DROP VIEW Educ1Model;
DROP VIEW Educ1StartOdds;
DROP VIEW Educ2AllowedDelays;
DROP VIEW Educ2DelayedProgressionIntake;
DROP VIEW Educ2DelayedRepetitionIntake;
DROP VIEW Educ2DirectProgressionIntake;
DROP VIEW Educ2DirectRepetitionIntake;
DROP VIEW Educ2PeriodSuccess;
DROP VIEW EducOneDropoutGrade;
DROP VIEW EducOneEntryAge;
DROP VIEW EducTrans1;
DROP VIEW EducTrans2;
DROP VIEW EmigrationDestination;
DROP VIEW EmigrationRatesDistrict;
DROP VIEW EndSchoolOneYear;
DROP VIEW EthnicTransmission;
DROP VIEW EthnicityImmigrantsScratch;
DROP VIEW FertilityModel;
DROP VIEW FirstBirthRates;
DROP VIEW HCICoefficients;
DROP VIEW HigherOrderBirthsPara;
DROP VIEW ImmiPoolDestination;
DROP VIEW ImmiPoolSize;
DROP VIEW ImmiScratchDestination;
DROP VIEW InUnionProbNoChildren;
DROP VIEW InUnionProbWithChildren;
DROP VIEW LifeExpectancy;
DROP VIEW MicroDataInputFile;
DROP VIEW MigrationDestination;
DROP VIEW MigrationProbability;
DROP VIEW MigrationTryKeepingFamiliesTogether;
DROP VIEW ModelBackmigration;
DROP VIEW ModelEmigration;
DROP VIEW ModelImmigrationFromPools;
DROP VIEW ModelImmigrationFromScratch;
DROP VIEW ModelMigration;
DROP VIEW MortalityModel;
DROP VIEW MortalityTable;
DROP VIEW NumberImmigrantsFromScratch;
DROP VIEW PartnerAgeDistribution;
DROP VIEW PartnerCharacteristicDistribution;
DROP VIEW PreNatalCareOdds;
DROP VIEW PreSchoolAttendance;
DROP VIEW ProbStayWithMother;
DROP VIEW ProportionStunting;
DROP VIEW SchoolOneInterruptionRate;
DROP VIEW SchoolOneRepetitionRate;
DROP VIEW SchoolQuality;
DROP VIEW SexRatio;
DROP VIEW SimulationEnd;
DROP VIEW SimulationSeed;
DROP VIEW StartPopSampleSize;
DROP VIEW StartSchoolOneYear;
DROP VIEW TotalFertilityRate;
DROP VIEW Union1Choice;
DROP VIEW Union1ParametersCMN;
DROP VIEW Union1ParametersHazards;

--
-- drop output tables compatibility views
--
DROP VIEW PopPyramidByEduc;
DROP VIEW TabChildVaccination;
DROP VIEW TabEduc15ByDistrict;
DROP VIEW TabEduc15ByDistrictBirth;
DROP VIEW TabEducFateByGroup;
DROP VIEW TabEducFateDistrYob;
DROP VIEW TabHCIDistrict;
DROP VIEW TabImmunizationChildren;
DROP VIEW TabPopProvAgeEducSex;
DROP VIEW TabPrenatCare;
DROP VIEW TabPrimSchoolEntries;
DROP VIEW TabPrimSchoolGraduations;
DROP VIEW TabPrimSchoolOutOfSchool9to11;
DROP VIEW TabPrimarySchoolPlanning;
DROP VIEW TabSchool2AttainmentsTab;
DROP VIEW TabSchool2TrackTab;
DROP VIEW tabBirthsYearPlace;
DROP VIEW tabEducationFateGeobirYob;
DROP VIEW tabHCI;
DROP VIEW tabHavingSpouse;
DROP VIEW tabMigrationOriginDestination;
DROP VIEW tabPopulationYearPlace;
DROP VIEW tabPreSchool;
DROP VIEW tabStuntingSexRegMotherYob;
DROP VIEW tabStuntingSexRegYob;

