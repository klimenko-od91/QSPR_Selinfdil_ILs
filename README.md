# QSPR_Selinfdil_ILs
This repository contains an archive with scripts and data files that are needed to reproduce the QSPR model development, validation and virtual screening for Selectivity at infinite dilution of ionic liquids. The repository is created by demand from Journal of cheminformatics.

The scripts that belong to the same cross-validation (CV) scheme and Property type (log10 [S] OR big IDAC flag) must be run sequentially.
The CV split scripts are differentiated with 80,50,20 suffices.
The property scripts either have BIDAC prefix for big IDAC flag modeling or no prefix at all in case of log10 [S] modeling.

-------------------------------------------The example of running sequential scripts for 80% CV split in log10 [S] modeling------------------------------------------------------

First, Selectivity_QSPR_server_80.R script must be executed in order to build models and get training set and CV results.
Next, Selectivity_prediction_server_80.R script is run to predict optimization set. Then, use test_set_stat_80.R script to get the optimization set statistics for all models and get training set and CV results for the best model according to the decision function. Test set prediction is done using Selectivity_prediction_external_80.R script.
After that, use external_test_stats_80.R to get the statistical assessment for the external test set. 

IDAC flag models must be made after the log10 [S] modeling.

The applicability domain scripts are the same for all schemes. There is one for the optimization set (AD_calculation.R) and one for the test set (AD_calculation_external.R).
The applicability domain of the particular set (test, external test or combinatorial library) must be calculated before the  prediction is made.

--------------------------------------------Virtual screening with the best models-----------------------------------------------------------------------------------------------
Next, comblib_prediction_server_andod.R is used to predict the log10 [S] of Combinatorial library for aniline/n-dodecane system.
comblib_prediction_server_BIDAC_andod.R is used to predict the big IDAC flag of Combinatorial library for aniline/n-dodecane system.
The applicability domain assessment is done using AD_calculation_andod.R
