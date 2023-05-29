

clear all 
set more off, perm
macro drop _all
set mem 5000m
set linesize 120
set maxvar 10000
*****************************
* Directories
*****************************

* KFS Data Directory
*global raw "K:\Datasets\Common\kfs8"
global raw "K:\Datasets\Common\KFS Manual and Data"

* Working Directories
global prog "H:\DKW\prog"
global work "H:\DKW\work"
global log "H:\DKW\log"
global output "H:\DKW\output"

*****************************
* 0. Create dataset for analysis
*****************************
* Use the logically imputed long form dataset
use mprid c2_owners* c3a_owner_operators* c1z2_legal_status* c5_num_employees*  ///
	tot_assets_? tot_liab_? tot_debt_? tot_equity_? naics_code* f29_assetval_cash* ///
	f16a_rev* f17a_total* f18a_wage_exp* g10c_net_worth_01* wgt* cswgt_final* ///
	classf* using "${raw}\KFS8_LI.dta",clear
capture drop tot_equity_all
* rename some variables
forvalues y = 0/7 {
	ren c2_owners_`y' n_owners_`y'
	ren c3a_owner_operators_`y' n_actowners_`y'
	ren c1z2_legal_status_`y' lfo_`y'
	ren c5_num_employees_`y' size_`y'
	ren f16a_rev_amt_`y' tot_rev_`y'
	ren f17a_total_exp_amt_`y' tot_exp_`y'
	ren f18a_wage_exp_amt_`y' wage_exp_`y'
	ren g10c_net_worth_01_`y' networth_`y' /* net worth of the first owner, only available from 2008 */
	ren f29_assetval_cash_`y' cash_assets_`y'
	if `y'>0 {
		ren wgt_`y'_long wgt_long_`y'
	}	
}

* reshape 
reshape long lfo_ naics_code_ n_owners_ n_actowners_ size_ tot_rev_ tot_exp_ wage_exp_ networth_ classf_ tot_liab_ tot_debt_ ///
	tot_assets_ cash_assets_  tot_equity_  wgt_long_ cswgt_final_, i(mprid) j(y)

gen year = 2004+y
drop y

foreach var in lfo naics_code n_owners n_actowners size tot_rev tot_exp wage_exp networth classf tot_assets tot_liab tot_debt ///
	cash_assets tot_equity cswgt_final wgt_long {
		ren `var'_ `var'
		
	}

	
* Create new variables

* Gen ind_fin: businesses in finance
gen ind_fin =  (naics_code>=520000 & naics_code<=529999)
	
* Real asset and net debt
gen RealAsset = tot_assets - cash_assets 
gen NetDebt = tot_debt - cash_assets
gen VA = tot_rev - (tot_exp-wage_exp)

* Gen size_cat to indicate small firms
gen size_cat = .
replace size_cat = 0 if size==0
replace size_cat = 1 if size>0 & size<500
replace size_cat = 2 if size>=500 & size<.
replace size_cat = -1 if classf ~= 6
* Gen lead var of size_cat
sort mprid year
by mprid: gen size_cat_lead = size_cat[_n+1] if year+1==year[_n+1]

label define size_cat -1 "inactive" 0 "size = 0" 1 "size 1-499" 2 "size 500+"
label value size_cat  size_cat_lead size_cat

compress
save ${work}\KFS8_clean.dta,replace



*****************************
* 1. Compute Moments
*	1.1 Debt to Asset Ratio of entrants and incumbents
*	1.2 capital to revenue ratio
*	1.3 Autocorrelation of employment and revenue
*	1.4 Labor Share
*	1.5 firm dynamics: transition rates
*	1.6 Capital adjustment, investment
*	1.7 Debt and savings distribution
*****************************
cap log close
log using "${log}\DKW_1_crdata.txt",replace nomsg text

use ${work}\KFS8_clean.dta,clear

quietly {
*	1.1 Debt to Asset Ratio of entrants and incumbents
*  a. agg debt to agg asset ratio of all firms
capture drop asset_weight
gen asset_weight = RealAsset*cswgt_final if size_cat==1
capture drop debt_weight
gen debt_weight = max(NetDebt,0)*cswgt_final if size_cat==1
capture drop aggDebt aggRealAsset
egen aggDebt = total(debt_weight)
egen aggRealAsset = total(asset_weight)
capture drop agg_debt_asset_ratio
gen agg_debt_asset_ratio = aggDebt/aggRealAsset
sum agg_debt_asset_ratio
* number of observations
count if debt_weight<. & asset_weight <.
capture drop Nobs_agg_debt_asset_ratio
gen Nobs_agg_debt_asset_ratio = r(N)

* by year
capture drop aggDebt_year
bysort year:   egen aggDebt_year = total(debt_weight)
capture drop aggRealAsset_year
bysort year:   egen aggRealAsset_year = total(asset_weight)
capture drop agg_debt_asset_ratio_year
gen agg_debt_asset_ratio_year = aggDebt_year / aggRealAsset_year
* number of observations
capture drop Nobs_agg_debt_asset_ratio_year
gen Nobs_agg_debt_asset_ratio_year = .
forvalues y = 2004/2011 {
    count if debt_weight<. & asset_weight<. & year == `y'
	replace Nobs_agg_debt_asset_ratio_year = r(N) if year == `y'
}

noisily display "Ave debt asset ratio = " agg_debt_asset_ratio  "	No. observations = " Nobs_agg_debt_asset_ratio
noisily display "Ave debt asset ratio by year and No. of observations: " 
noisily tabstat  agg_debt_asset_ratio_year Nobs_agg_debt_asset_ratio_year,s(mean) by(year)  nototal



* b. agg debt-to-asset ratio of those with debt and fraction of firms with debt
* Fraction of firms with net debt
capture drop has_NetDebt
gen has_NetDebt = (NetDebt>0 & NetDebt<.)
sum has_NetDebt [aw = cswgt_final ] if size_cat==1
noisily display "ave has_NetDebt = " r(mean)  "	No. observations = " r(N)

noisily tabstat has_NetDebt [aw = cswgt_final ] if size_cat==1,by(year) s(mean N) 


* agg Debt to payroll ratio conditional on positive debt
capture drop payroll_weight
gen payroll_weight = wage_exp*cswgt_final if size_cat==1 & has_NetDebt ==1
capture drop debt_weight
gen debt_weight = NetDebt*cswgt_final if size_cat==1 & has_NetDebt ==1
capture drop aggDebt 
capture drop aggPayroll
egen aggDebt = total(debt_weight)
egen aggPayroll = total(payroll_weight)
capture drop agg_debt_payroll_ratio_cond 
gen agg_debt_payroll_ratio_cond = aggDebt/aggPayroll
* number of observations 
count if payroll_weight<. & debt_weight<.
capture drop Nobs_agg_debt_payroll_ratio_cond
gen Nobs_agg_debt_payroll_ratio_cond = r(N)


* by year
capture drop aggDebt_year
bysort year: egen aggDebt_year = total(debt_weight)
capture drop aggPayroll_year
bysort year: egen aggPayroll_year = total(payroll_weight)
capture drop agg_debt_payroll_ratio_cond_year
gen agg_debt_payroll_ratio_cond_year = aggDebt_year / aggPayroll_year
* number of observations
capture drop Nobs_agg_debt_payroll_cond_year
gen Nobs_agg_debt_payroll_cond_year = .
forvalues y = 2004/2011 {
    count if payroll_weight<. & debt_weight<. & year== `y'
	replace Nobs_agg_debt_payroll_cond_year = r(N) if year == `y'
}

noisily display "ave debt payroll ratio (conditional) = " agg_debt_payroll_ratio_cond "	No. observations = " Nobs_agg_debt_payroll_ratio_cond
noisily display "ave debt payroll ratio (conditional) by year and No. observations:"
noisily tabstat  agg_debt_payroll_ratio_cond_year Nobs_agg_debt_payroll_cond_year,s(mean) by(year) 


* agg Financial asset to payroll ratio conditional on netDebt<0
capture drop payroll_weight
gen payroll_weight = wage_exp*cswgt_final if size_cat==1 & has_NetDebt ==0
capture drop debt_weight
gen debt_weight = NetDebt*cswgt_final if size_cat==1 & has_NetDebt ==0
capture drop aggDebt 
capture drop aggPayroll
egen aggDebt = total(debt_weight)
egen aggPayroll = total(payroll_weight)
capture drop agg_debt_payroll_ratio_cond 
gen agg_debt_payroll_ratio_cond = aggDebt/aggPayroll
* Number of observations
count if payroll_weight<. & debt_weight<.
capture drop Nobs_agg_debt_payroll
gen Nobs_agg_debt_payroll = r(N)

* by year
capture drop aggDebt_year
bysort year: egen aggDebt_year = total(debt_weight)
capture drop aggPayroll_year
bysort year: egen aggPayroll_year = total(payroll_weight)
capture drop agg_debt_payroll_ratio_cond_year
gen agg_debt_payroll_ratio_cond_year = aggDebt_year / aggPayroll_year
* number of observations
capture drop Nobs_agg_debt_payroll_cond_year
gen Nobs_agg_debt_payroll_cond_year = .
forvalues y = 2004/2011 {
    count if payroll_weight<. & debt_weight<. & year== `y'
	replace Nobs_agg_debt_payroll_cond_year = r(N) if year == `y'
}

noisily display "ave financial asset payroll ratio (conditional) = " agg_debt_payroll_ratio_cond "	No. observations = " Nobs_agg_debt_payroll
noisily display "ave financial asset to payroll ratio (conditional) by year and No. observations:"
noisily tabstat  agg_debt_payroll_ratio_cond_year Nobs_agg_debt_payroll_cond_year,s(mean) by(year) 


} /*quietly*/

*	1.2 capital to revenue ratio

quietly{
* a. agg k to agg rev ratio
capture drop k_weight
gen k_weight = RealAsset*cswgt_final if size_cat==1
capture drop rev_weight
gen rev_weight = tot_rev * cswgt_final if size_cat==1
capture drop aggK
egen aggK = total(k_weight)
capture drop aggRev
egen aggRev = total(rev_weight)
capture drop agg_k_rev_ratio 
gen agg_k_rev_ratio= aggK/aggRev
* Number of observations
count if k_weight<. & rev_weight<.
capture drop Nobs_k_rev 
gen Nobs_k_rev = r(N)

* by year
capture drop aggK_year
bysort year: egen aggK_year = total(k_weight)
capture drop aggRev_year
bysort year: egen aggRev_year = total(rev_weight)
capture drop agg_k_rev_ratio_year
gen agg_k_rev_ratio_year = aggK_year / aggRev_year
* number of observations
capture drop Nobs_k_rev_year
gen Nobs_k_rev_year = .
forvalues y = 2004/2011 {
    count if k_weight<. & rev_weight<. & year == `y'
	replace Nobs_k_rev_year = r(N) if year == `y'
}

noisily display "ave K rev ratio = " agg_k_rev_ratio "	No. observations = " Nobs_k_rev
noisily display "ave K rev ratio by year and No. observations"
noisily tabstat  agg_k_rev_ratio_year Nobs_k_rev_year,s(mean) by(year)

* b. agg k to agg VA ratio
capture drop k_weight
gen k_weight = RealAsset*cswgt_final if size_cat==1
capture drop va_weight
gen va_weight = VA * cswgt_final if size_cat==1
capture drop aggK
egen aggK = total(k_weight)
capture drop aggVA
egen aggVA = total(va_weight)
capture drop agg_k_va_ratio 
gen agg_k_va_ratio= aggK/aggVA
* number of observations
count if k_weight<. & va_weight<.
capture drop Nobs_k_va
gen Nobs_k_va = r(N)

* by year
capture drop aggK_year
bysort year: egen aggK_year = total(k_weight)
capture drop aggVA_year
bysort year: egen aggVA_year = total(va_weight)
capture drop agg_k_va_ratio_year
gen agg_k_va_ratio_year = aggK_year / aggVA_year
* Number of observations
capture drop Nobs_k_va_year 
gen Nobs_k_va_year = .
forvalues y = 2004/2011 {
    count if k_weight<. & va_weight<. & year == `y'
	replace Nobs_k_va_year = r(N) if year == `y'
}

noisily display "ave K VA ratio = " agg_k_va_ratio "	No. observations = " Nobs_k_va
noisily display "ave k VA ratio by year and No. observations"
noisily tabstat  agg_k_va_ratio_year Nobs_k_va_year,s(mean) by(year)


* b. agg k to agg payroll ratio
capture drop k_weight
gen k_weight = RealAsset*cswgt_final if size_cat==1
capture drop wage_weight
gen wage_weight = wage_exp * cswgt_final if size_cat==1
capture drop aggK
egen aggK = total(k_weight)
capture drop aggW
egen aggW = total(wage_weight)
capture drop agg_k_w_ratio 
gen agg_k_w_ratio= aggK/aggW
* number of observations
count if k_weight<. & wage_weight<.
capture drop Nobs_k_w
gen Nobs_k_w = r(N)

* by year
capture drop aggK_year
bysort year: egen aggK_year = total(k_weight)
capture drop aggW_year
bysort year: egen aggW_year = total(wage_weight)
capture drop agg_k_w_ratio_year
gen agg_k_w_ratio_year = aggK_year / aggW_year
* number of observations
capture drop Nobs_k_w_year
gen Nobs_k_w_year = .
forvalues y = 2004/2011 {
    count if k_weight<. & wage_weight<. & year == `y'
	replace Nobs_k_w_year = r(N) if year == `y'
}

noisily display "ave K payroll ratio = " agg_k_w_ratio "	No. observations = " Nobs_k_w
noisily display "ave K payroll ratio by year and No. observations"
noisily tabstat  agg_k_w_ratio_year Nobs_k_w_year,s(mean) by(year)


} /*quietly*/

*	1.3 Autocorrelation of employment and revenue

capture drop size0 tot_rev0
sort mprid year
by mprid: gen size0 = size[_n-1] if year==year[_n-1]+1
by mprid: gen tot_rev0  = tot_rev[_n-1] if year==year[_n-1]+1
* log revenue
capture drop logRev
gen logRev = log(tot_rev) if tot_rev>0
capture drop logRev0
by mprid: gen logRev0 = logRev[_n-1] if year==year[_n-1]+1
* log va
capture drop logVA
gen logVA = log(VA) if VA>0
capture drop logVA0
by mprid: gen logVA0 = logVA[_n-1] if year==year[_n-1]+1


quietly{
* Correlation between size and lagged size
corr size size0 [aw = cswgt_final ] if size0>0 & size>0 
noisily display "autocorr emp all years = " r(rho)	"	No. observations = " r(N)
forvalues y = 2005/2011 {
    corr size size0 [aw = cswgt_final ] if size0>0 & size>0 & year == `y'	
	noisily display "autocorr emp `y'  = " r(rho)	"	No. observations = " r(N)

}

* Correlation between log revenue and lagged log revenue
corr logRev logRev0 [aw = cswgt_final ] if size0>0 & size>0
noisily display "autocorr log rev all years = " r(rho)	"	No. observations = " r(N)

forvalues y = 2005/2011 {
	corr logRev logRev0 [aw = cswgt_final ] if size0>0 & size>0 & year == `y'	
	noisily display "autocorr log rev `y' = " r(rho)	"	No. observations = " r(N)
}

* Correlation between log revenue and lagged log revenue
corr logRev logRev0 [aw = cswgt_final ] if size0>0 & size>0
noisily display "autocorr log VA all years = " r(rho)	"	No. observations = " r(N)

forvalues y = 2005/2011 {
	corr logVA logVA0 [aw = cswgt_final ] if size0>0 & size>0 & year == `y'	
	noisily display "autocorr log VA `y' = " r(rho)	"	No. observations = " r(N)
}
}/*quietly*/


log close
