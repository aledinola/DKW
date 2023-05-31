Instructions to reproduce the output of "Rescue Policies for Small Businesses during the COVID-19 Recession".
Authors: Alessandro Di Nola, Leo Kaas, Haomin Wang.

*********************************************************************************************************
*** Notes on SOFTWARE, OPERATING SYSTEM, and COMPUTATION TIME ***
*********************************************************************************************************
- All programs were run on Windows 10 on a computer with Intel-core-i7 CPU 3.00 GHz and 64 GB RAM.
- Matlab codes are run in Matlab R2022b and make use of the following Matlab toolboxes: 
  Parallel Computing Toolbox
  Optimization Toolbox
  Global Optimization Toolbox
- A few subroutines are coded in Fortran and compiled using the MEX utility in Matlab. 
- The Fortran codes make use of the OpenMP library. 
- The runtime for the code in the folders containing Matlab codes (baseline, conditional and alternative) was of about 37 seconds for the steady-state version of the model, and about 3 hours for the transition after a pandemic shock. 

************************* DATA ********************
Folder "Data" contains public data and codes to compute moments and other statistics. 

Data sources:
1) SUSB: see subfolder "SUSB". The data are downloaded from https://www.census.gov/programs-surveys/susb/data/tables.html.
2) BDS: see subfolder "BDS". The data are downloaded from https://www.census.gov/data/datasets/time-series/econ/bds/bds-datasets.html.
3) FRED: see subfolder "FRED". The data are downloaded from https://fred.stlouisfed.org.
4) Business expense data: see folder "other_data". The data are downloaded from https://washingtonpost.com/wp-srv/special/business/costofrunningabusiness.html.
5) Small firm output data: see folder "other_data". The data are from Bloom, Fletcher and Yeh (2021) "The Impact of COVID-19" on US Firms."
6) Kauffman Firm Survey. We use KFS confidential microdata file. To access the dataset, an application is required to be made to the Kauffman Foundation. See link: https://www.kauffman.org/entrepreneurship/research/kauffman-firm-survey/.

Code to create moments and other statistics
1) Run prog/main.do in STATA to generate moments and statistics. Outputs are stored in subfolder "moments."
	- data_moments.txt: a text file containing data targets for the steady state calibration.
	- FRED_Graphs and FRED_Graphs_ee are Figures 1a and 1b in the paper.
2) The program prog/main.do uses raw data saved in subfolders SUSB, BDS, FRED and other_data as well as statistics computed based on the KFS data. 
Since the KFS data is confidential, we provide only the STATA code that is used to generate the KFS moments. See file "prog/main_KFS.do."


************************* MAIN RESULTS (Baseline grant, laissez-faire, targeted grant) ********************
Codes to generate main results are in folder "baseline"

CODE FOR MODEL SOLUTION 

To produce steady-state and transition results and save to corresponding mat files.
1) To reproduce the results for the steady-state economy:
1.1) Open 'main.m', set do_calib=0, do_save=1. Open 'set_parameters.m' and set par.do_howard = 1. 
1.2) Run 'main.m' This will run the steady-state economy and save results in a mat file called 'ss.mat'
     in the subfolder 'mat'.
     Overall, the replication of the steady-state takes approximately 30 seconds.

2) To reproduce the results for the transition after a pandemic shock:
2.1) Open 'main.m', set do_calib=1, do_save=1. Open 'set_parameters.m' and set par.do_howard = 0.
     There are several options for different type of policy:
     grant_flag=1 and grant_target=0 ==> Baseline PPP grant policy (untargeted) ==> Output saved in 'grant_baseline.mat'
     grant_flag=0 ==> Counterfactual NO grant policy (laissez-faire) ==> Output saved in 'nogrant.mat'
     grant_flag=5 and grant_target=2 ==> Counterfactual targeted grant ==> Output saved in 'grant_targslim.mat'
     grant_flag=6 and grant_target=2 ==> Counterfactual targeted grant (large grant) ==> Output saved in 'grant_targslim_large.mat'
     grant_flag=7 and grant_target=2 ==> Counterfactual targeted grant (small grant) ==> Output saved in 'grant_targslim_small.mat'
     Overall, the replication of the transition takes approximately three hours.
     All mat files are saved in the 'mat' subfolder.

CODE FOR FIGURES AND TABLES 

3)   Run 'main_tables.m'. This will generate the following tables in tex format, saved in subfolder 'tables':
     - exo_parameters.tex:     correpsonding to Table 1 "External parameters" in the paper
     - parameters.tex:         correpsonding to Table 2 "Internal parameters" in the paper
     - moments.tex:            corresponding to Table 3 "Model fit" in the paper
     - comp_parameters.tex:    Table with Computational parameters, not present in the paper
     - steady_state.tex:       Table with some steady-state values (not in the paper)
     - tran_shocks.tex:        Calibrated pandemic shock parameters (Table 4)
     - transition_moments.tex: Pandemic impact (Table 5)
     - cost_grants.tex:        Cost per job saved in small firms (Table 6)
     - cost_grants_robust.tex  Cost per job saved in small firms (Table 9)

4) Run 'main_plots.m'. This will generate the following figures saved in subfolder 'figures':
     - FirmSize:      Figure 2a Firm Shares
     - empShare:      Figure 2b Employment shares  
     - exit_type_ss:  Figure 3a Exit decisions
     - entry_type_ss: Figure 3b Entry decisions
     - inv_type_ss:   Figure 4a Investment decision
     - borrow_pol_ss: Figure 4b Borrowing decision
     - model_mom_trans: Figure 5a Aggregate variables, adjustment to the pandemic shock in the baseline calibration
     - model_mom_trans_bis: Figure 5b Entry and exit rates, adjustment to the pandemic shock in the baseline calibration
     - irf_ave_bk:    Firm indebtness measured as the average debt-capital ratio (FIgure 10)
     - emp_exit_rate_lbins: Impact of the grant on the exit rate and empl. shares by firm size bins (Figure 8)
     - zombie.tex:    Zombie firms (Tables 7 and 10), saved in 'tables'
     - capadj_rate_ss.tex: Four types of capital adjustment rates in the steady-state (Table 8)
     - cum_capadj_decomp_sr: Short-run cumulative excess investment by types of capital adjustment (Figure 9a) 
     - cum_capadj_decomp_mr: Medium-run cumulative excess investment by types of capital adjustment (Figure 9b) 
     - cum_capadj_decomp_lr: Long-run cumulative excess investment by types of capital adjustment (Figure 9c) 
     - cum_capadj_decomp_all: Q1-Q40 cumulative excess investment by types of capital adjustment (Figure 9d) 
     - welfare_analysis.tex: see Section 5.7 (Welfare analysis), last paragraph 

5) Run 'make_plots_compare.m'. This will generate the impulse responses in Figures 6-7-17.

6) Run 'cum_impact_compare.m'. This will generate bar figures for short and long-run effects in Figures 11-12.
 

************************* CONDITIONAL GRANT *************************************************************
The programs that generates results of the conditional grant are in the folder 'conditional'. 

1) To reproduce the results for the steady-state economy:
1.1) Open 'main.m', set do_calib=0, do_save=1. Open 'set_parameters.m' and set par.do_howard = 1. 
1.2) Run 'main.m' This will run the steady-state economy and save results in a mat file called 'ss.mat'
     in the subfolder 'mat'.
     Overall, the replication of the steady-state takes approximately 30 seconds.
 
2) To reproduce the results for the transition after a pandemic shock:
2.1) Open 'main.m', set do_calib=1, do_save=1. Open 'set_parameters.m' and set par.do_howard = 0.
     There are several options for different type of policy:
     grant_flag=1 and grant_target=0 ==> Baseline PPP grant policy (untargeted) ==> Output saved in 'grant_baseline.mat'
     grant_flag=0 ==> Counterfactual NO grant policy (laissez-faire) ==> Output saved in 'nogrant.mat'
     All mat files are saved in the 'mat' subfolder.

3) Run 'make_plots_compare.m'. This will generate the impulse responses in Figures 13-14 (in Appendix C).

4) Run 'cum_impact_compare.m'. This will generate bar figures for short and long-run effects in Figures 15-16 (in Appendix C).

************************* ALTERNATIVE CALIBRATION *******************************************************
The programs that generates results for the alternative calibration of the productivity process (Appendix E) 
are in the folder 'alternative'. 

1) To reproduce the results for the steady-state economy:
1.1) Open 'main.m', set do_calib=0, do_save=1. Open 'set_parameters.m' and set par.do_howard = 1. 
1.2) Run 'main.m' This will run the steady-state economy and save results in a mat file called 'ss.mat'
     in the subfolder 'mat'.
     Overall, the replication of the steady-state takes approximately 30 seconds.

2) To reproduce the results for the transition after a pandemic shock:
2.1) Open 'main.m', set do_calib=1, do_save=1. Open 'set_parameters.m' and set par.do_howard = 0.
     There are several options for different type of policy:
     grant_flag=1 and grant_target=0 ==> Baseline PPP grant policy (untargeted) ==> Output saved in 'grant_baseline.mat'
     grant_flag=0 ==> Counterfactual NO grant policy (laissez-faire) ==> Output saved in 'nogrant.mat'
     All mat files are saved in the 'mat' subfolder.

3) Run 'main_tables.m'. This gives you Table 11 (Model Fit, alternative calibration) in Appendix E.
4) Run 'make_plots_compare.m'. This will generate the impulse responses in Figures 18-19 (in Appendix E).

************************* MEX FILES INSTRUCTIONS ********************************************************

The codes calls two MEX functions:
- sub_vfi_onestep_mex
- myinterp1q

The source files (written in Fortran) are stored in subfolder 'tools\mex'. Before running the codes, the user should 
compile each MEX file into an executable (machine dependent). We write here the instructions for Windows and MAC.

Step 1: Change directory to 'tools\mex' in the matlab command window.
Step 2: Compile
- WINDOWS
mex -v OPTIMFLAGS="/O3 /Qopenmp /Qprec-div- /QxHost /DNDEBUG" -R2018a sub_vfi_onestep_mex.f90 -output sub_vfi_onestep_mex
mex -v OPTIMFLAGS="/O3 /Qprec-div- /QxHost /DNDEBUG" -R2018a myinterp1q.f90 -output myinterp1q
- MAC
mex -v FOPTIMFLAGS="-O3 -qopenmp "  sub_vfi_onestep_mex.f90 -output sub_vfi_onestep_mex
mex -v FOPTIMFLAGS="-O3"  myinterp1q.f90 -output myinterp1q

************************* CODE REFERENCES ********************************************************
We make use of a few other Matlab functions written by others:

FUNCTION      AUTHOR(S)
paretojo      In Hwan Jo and Tatsuro Senga (publicly available on ideas/repec Replication code RED 2019)
bddparetocdf  In Hwan Jo and Tatsuro Senga (publicly available on ideas/repec Replication code RED 2019)
fminsearchcon John D'Errico (publicly available on Mathworks website)
markovapprox  Eva Carceles-Poveda

In Hwan Jo & Tatsuro Senga (2019). 
"Code and data files for "Aggregate Consequences of Credit Subsidy Policies: Firm Dynamics and Misallocation"," 
Computer Codes 17-402, Review of Economic Dynamics.

John D'Errico (2023). fminsearchbnd, fminsearchcon (https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon), MATLAB Central File Exchange. Retrieved May 29, 2023.



