
* "Rescue Policies for Small Businesses During the Covid-19 Recession"
* By Alessandro Di Nola, Leo Kaas and Haomin Wang

clear all
set more off, perm
macro drop _all
set mem 5000m
set matsize 1000

set linesize 120   



**************************************************************
************* 		Directories 		**********************
**************************************************************

* Main directory
* Change [PATH] to the directory of the replication package
global path 	"[PATH]"

* Data directories
global SUSB 	${path}/DATA/SUSB
global BDS 	${path}/DATA/BDS
global FRED ${path}/DATA/FRED
global other ${path}/DATA/other_data
* Program directory
global prog 	${path}/DATA/prog
* Intermediate dataset directory
global work 	${path}/DATA/work
* Output directory
global moments ${path}/DATA/moments


**************************************************************
************* 	A.	Read and clean SUSB data 	**************
* Excel tables are downloaded here: 
* https://www.census.gov/programs-surveys/susb/data/tables.html
**************************************************************
* A.1 Read raw SUSB tables by year and combine them into one file
* A.2 Change variable names and add labels, create indicators for covid-impacted sectors
**************************************************************


			
* A.1 Read raw SUSB tables by year and combine them into one file
*2010
import excel using ${SUSB}/us_naicssector_lfo_2010.xls, ///
			cellrange(A8:K1380)  clear
save ${work}/SUSB_empdist_2010.dta,replace
* 2011
import excel using ${SUSB}/us_naicssector_lfo_2011.xls, ///
			cellrange(A8:K1377)  clear
save ${work}/SUSB_empdist_2011.dta,replace
* 2012
import excel using ${SUSB}/us_naicssector_lfo_2012.xls, ///
			cellrange(A9:M1378)  clear
save ${work}/SUSB_empdist_2012.dta,replace
* 2013
import excel using ${SUSB}/us_naicssector_lfo_2013.xlsx, ///
			cellrange(A8:K1374)  clear
save ${work}/SUSB_empdist_2013.dta,replace
* 2014
import excel using ${SUSB}/us_naicssector_lfo_2014.xlsx, ///
			cellrange(A9:K1383)  clear
save ${work}/SUSB_empdist_2014.dta,replace
* 2015
import excel using ${SUSB}/us_naicssector_lfo_2015.xlsx, ///
			cellrange(A9:K1379)  clear
save ${work}/SUSB_empdist_2015.dta,replace
* 2016
import excel using ${SUSB}/us_naicssector_lfo_2016.xlsx, ///
			cellrange(A10:K1380)  clear
save ${work}/SUSB_empdist_2016.dta,replace
* 2017
import excel using ${SUSB}/us_naicssector_lfo_2017.xlsx, ///
			cellrange(A10:K1355)  clear
save ${work}/SUSB_empdist_2017.dta,replace

* Combine data from all years
forvalues year = 2010/2017 {
		use ${work}/SUSB_empdist_`year'.dta,clear
		* Gen year
		display "year `year'"
		gen year = `year'

		* Save file
		if (`year'>2010) {
			append using "${work}/SUSB_empdist.dta"
		}
		compress
		save ${work}/SUSB_empdist.dta,replace
		rm ${work}/SUSB_empdist_`year'.dta
}

* A.2 Change variable names and add labels, create indicators for covid-impacted sectors
use ${work}/SUSB_empdist.dta,clear
* Gen lfo: legal form
gen lfo = substr(A,1,1)
destring lfo, replace
drop A
* NAICS code
gen naics = substr(B,1,2)
destring naics,replace force
replace naics = 0 if naics ==.
drop B C
* Employment size categories
gen emp_c = substr(D,1,1)
destring emp_c,replace
drop D


* Number of firms and establishment
ren E N_firm
ren F N_estab
* Number of total employment
ren G tot_emp
* Annual Payroll
ren J tot_payroll
* Annual Receipts
ren L tot_receipts
* Clean up
drop H I K M
order year lfo naics emp_c N_firm N_estab tot_emp tot_payroll tot_receipts
* Labels
label define legal_form 1 "Total" 2 "Corporation" 3 "S-Corp" 4 "Partnership" ///
			5 "Sole Proprietorship" 6 "Non-Profit" 7 "Government" 8 "Other"
label define naics 0 "All" 11 "Agriculture, forestry, fishing and hunting" ///
	21 "Mining, quarrying, and oil and gas extraction" 22 "Utilities" ///
	23 "Construction" 31 "Manufacturing" 42 "Wholesale trade" 44 "Retail trade" ///
	48 "Transportation and warehousing" 51 "Information" 52 "Finance and insurance" ///
	53 "Real estate and rental and leasing" 54 "Professional, scientific, and technical services" ///
	55 "Management of companies and enterprises" 56 "Administrative and support and waste management and remediation services" ///
	61 "Educational services" 62 "Health care and social assistance" 71 "Arts, entertainment, and recreation" /// 
	72 "Accommodation and food services" 81 "Other services (except public administration)" ///
	99 "Industries not classified"

label define emp_c 1 "Total" 2 "0-4" 3 "5-9" 4 "10-19" 5 "<20" 6 "20-99" 7 "100-499" 8 "<500" 9 "500+"

label var lfo "Legal Form of Org."
label value lfo legal_form
label var naics "NAICS Code"
label value naics naics
label var emp_c "Employment Size Category"
label value emp_c emp_c
label var N_firm "Number of Firms"
label var N_estab "Number of Establishment"
label var tot_emp "Number of total employment"
label var tot_payroll "Annual Payroll ($1,000)"
label var tot_receipts "Annual Receipts ($1,000)"

* Define sectors impacted by covid-shock
capture drop impacted
gen impacted = inlist(naics, 44,71,72,81)
replace impacted = . if naics ==0
tab impacted [aw = tot_emp] if inlist(lfo, 2,3,4,5) & emp_c == 1 /* & year ==2017*/

* Drop emp categories: total,  <20 , <500
drop if inlist(emp_c,1,5,8)

compress
save ${work}/SUSB_empdist.dta,replace



**************************************************************
*************  B. Extract Moments from SUSB data 	**********
* B.1 Average employment for small and non-small businesses
* B.2 Relative size of small businesses and impacted businesses in terms of employment and receipts (only 2012)
* B.3 Employment and firm share by firm size category
* B.4 Labor share (payroll expense divided by revenue)
**************************************************************

use ${work}/SUSB_empdist.dta,clear

* B.1 Average employment for small and non-small businesses

* Drop rows for lfo = all, government, non-profit, and other.
drop if inlist(lfo,1, 6, 7, 8)
* Drop rows such that naics category is all
drop if naics ==0
* Gen small: indicator for firms with 0-499 employees
capture drop small
gen small = (inlist(emp_c,2,3,4,6,7))
label define small 1 "0-499 Employees" 0 "500+ Employees"
label value small small

quietly {
cap log close                                           
log using ${moments}/data_moments.txt, text replace   nomsg 
 log off

* Display average firm size in terms of employment by small and non-small firms
*	two definitions of small firms: small_nc and small_all

	preserve
		bysort year small: egen emp = total(tot_emp)
		bysort year small: egen firms = total(N_firm)
		gen ave_emp = emp/firms
		keep year small ave_emp
		duplicates drop 
		sort year small
		export delimited using "${moments}/SUSB_avefirmsize" , delim(",") nolabel replace 
		log on
		sum ave_emp if small == 1
		noisily display "avefirmsize	" r(mean) 
		log off
	restore


* B.2 Relative size of small businesses and impacted businesses in terms of employment and receipts (only 2012)

* rel size by employment
	preserve
		capture drop emp
		bysort year small impacted: egen emp = total(tot_emp)
		capture drop emp_year
		bysort year: egen emp_year = total(tot_emp)
		capture drop emp_share
		gen emp_share = emp/emp_year
		keep year small impacted emp emp_share
		duplicates drop 
		sort year small impacted	
		export delimited using "${moments}/SUSB_empshare" , delim(",") nolabel replace 
		log on
		sum emp_share if small==1 & impacted==1
		noisily display "empshare_small_imp	" r(mean) 
		sum emp_share if small==1 & impacted==0
		noisily display "empshare_small_unimp	" r(mean) 
		log off
	restore
	


* rel size by receipts
	preserve
		keep if year == 2012 /* receipts only availabe in 2012*/
		capture drop receipts
		bysort small impacted: egen receipts = total(tot_receipts)
		capture drop receipts_year
		egen receipts_year = total(tot_receipts)
		capture drop receipts_share
		gen receipts_share = receipts/receipts_year
		keep year small impacted receipts receipts_share
		duplicates drop 
		sort year small impacted	
		export delimited using "${moments}/SUSB_receiptshare" , delim(",") nolabel replace 
		log on
		sum receipts_share if small==1 & impacted==1
		noisily display "revshare_small_imp	" r(mean) 
		sum receipts_share if small==1 & impacted==0
		noisily display "revshare_small_unimp	" r(mean) 
		log off
	restore

* share of impacted firms
	preserve
		capture drop firms
		bysort year small impacted: egen firms = total(N_firm)
		capture drop firms_year
		bysort year: egen firms_year = total(N_firm)
		capture drop firms_share
		gen firms_share = firms/firms_year
		keep year small impacted firms firms_share
		duplicates drop 
		sort year small impacted	
		export delimited using "${moments}/SUSB_firmshare" , delim(",") nolabel replace 
		* Share of small firms that are impacted
		log on
		sum firms_share if small==1 & impacted==1
		local imp = r(mean) 
		sum firms_share if small==1 & impacted==0
		local unimp = r(mean)
		local share_imp = `imp'/(`imp'+`unimp')
		noisily display "share_impacted		`share_imp'" 
		log off
	restore

* B.3 Employment and firm share by firm size category


	preserve
		
		keep if  small ==1 
		bysort year emp_c: egen tot_emp_c = total(tot_emp) 
		bysort year emp_c: egen tot_firms = total(N_firm) 
		bysort year emp_c: egen tot_rev = total(tot_receipts)
		bysort year: egen emp = total(tot_emp) 
		bysort year: egen firms = total(N_firm)
		bysort year: egen rev = total(tot_receipts)
		gen emp_share = tot_emp_c/emp
		gen firm_share = tot_firms/firms
		gen rev_share = tot_rev/rev
		keep year emp_c tot_emp_c emp_share tot_firms firm_share rev_share
		duplicates drop 
		sort year emp_c	
		export delimited using "${moments}/SUSB_empdist" , delim(",") nolabel replace 
		log on
		foreach var in emp firm rev {
			* 0-4 employees
			sum `var'_share if emp_c ==2
			noisily display "`var'share_0_4		" r(mean) 
			* 5-9 employees
			sum `var'_share if emp_c ==3
			noisily display "`var'share_5_9		" r(mean)
			* 10-19 employees
			sum `var'_share if emp_c ==4
			noisily display "`var'share_10_19	" r(mean) 
			* 20-99 employees
			sum `var'_share if emp_c ==6
			noisily display "`var'share_20_99	" r(mean) 
			* 100-499 employees
			sum `var'_share if emp_c ==7
			noisily display "`var'share_100_499	" r(mean) 
		}
		log off
	restore
	
* B.4 payroll-to-sales ratio
	preserve
		keep if  small ==1 & year == 2012
		egen payroll_all = total(tot_payroll)
		egen receipts_all = total(tot_receipts)
		gen labor_share = payroll_all/receipts_all
		
		log on
		sum labor_share
		noisily display "payroll_sales_ratio		" r(mean) 
		
		log off
	restore
	
* End of quietly
log close
}

**************************************************************
************* 	C.	Read and clean BDS data 	**************
* Excel tables are downloaded here: 
* https://www.census.gov/data/datasets/time-series/econ/bds/bds-datasets.html
**************************************************************
* C.1 Read and clean BDS table 
* C.2 Compute moments from BDS
* C.2.a firm exit rate of small firms for all firms and by firm size
* C.2.b firm entry rate: entry rate = firms of age zero / total num of firms
* C.3 Average firm size of new entrants
* C.4 job creation and job destruction rate for continuing establishments and by firm size
* C.5 Weighted average of KFS moments (weighted by firm age distribution)
**************************************************************

* C.1 Read and clean BDS table 

*import delim using "${BDS}/bds2018_vcnaics3_fsize.csv",delim(",") varn(noname)  rowr(2) clear

import delim using "${BDS}/bds2018_fage_fsize.csv",delim(",") varn(noname)  rowr(2) clear


* Clean up dataset
* year
ren v1 year
* firm age 
gen fage = substr(v2,1,1)
drop v2
replace fage = "1" if fage=="a"
replace fage = "2" if fage=="b"
replace fage = "3" if fage=="c"
replace fage = "4" if fage=="d"
replace fage = "5" if fage=="e"
replace fage = "6" if fage=="f"
replace fage = "7" if fage=="g"
replace fage = "8" if fage=="h"
replace fage = "9" if fage=="i"
replace fage = "10" if fage=="j"
replace fage = "11" if fage=="k"
replace fage = "12" if fage=="l"
destring fage,replace
* firm size category
gen emp_c_str = substr(v3,1,1)
*keep if inlist(emp_c_str,"a","b","c", "d","e")
gen emp_c=.
replace emp_c = 1 if emp_c_str=="a"
replace emp_c = 2 if emp_c_str=="b"
replace emp_c = 3 if emp_c_str=="c"
replace emp_c = 4 if emp_c_str=="d"
replace emp_c = 5 if emp_c_str=="e"
replace emp_c = 6 if emp_c_str=="f"
replace emp_c = 7 if emp_c_str=="g"
replace emp_c = 8 if emp_c_str=="h"
replace emp_c = 9 if emp_c_str=="i"
replace emp_c = 10 if emp_c_str=="j"


compress
drop emp_c_str v3
* Number of firms and establishments
destring v4,replace force
ren v4 N_firms
destring v5,replace force
ren v5 N_estab
* total number of employment 
destring v6,replace force
ren v6 tot_emp
* denom: average employment between t and t-1
destring v7,replace force
ren v7 denom
* drop establishment entry and exit variables
drop v8 v9 v10 v11
* Total job creation and job destruction
destring v12,replace force
ren v12 job_creation
destring v17,replace force
ren v17 job_destruction
* Job creation and destruction of continuers 
destring v14,replace force
ren v14 job_creation_continuers
destring v19,replace force
ren v19 job_destruction_continuers
* drop variables related to job creation and job destruction 
*	because these are based on establishment concepts
drop v13 v15 v16 v18 v20-v24
* Count of firms that have exited in their entirety during the period
destring v25,replace force
ren v25 firmdeath_firms
* drop firmdeath_estabs
drop v26
* Count of employment associated with firm deaths.
destring v27, replace force
ren v27 firmdeath_emp

* Labels

label define emp_c 1 "1-4" 2 "5-9" 3 "10-19" 4 "20-99" 5 "100-499" 6 "500 to 999" ///
	7 "1000 to 2499" 8 "2500 to 4999" 9 "5000 to 9999" 10 "10000+"
label var emp_c "Employment Size Category"
label value emp_c emp_c
label var N_firm "Number of Firms"
label var N_estab "Number of Establishment"
label var denom "average of employment for t and t-1"
label var fage "Firm Age Category"
label define fage 1 "0" 2 "1" 3 "2" 4 "3" 5 "4" 6 "5"  7 "6-10" 8 "11-15" ///
	9 "16-20" 10 "21-25" 11 "26+" 12 "Left Censored"
label value fage fage

order year fage emp_c

compress
save ${work}/BDS_fage_fsize.dta,replace

* C.2.a firm exit rate of small firms for all firms and by firm size
use ${work}/BDS_fage_fsize.dta,clear
keep if emp_c>=1 & emp_c<=5
quietly{

cap log close                                           
log using ${moments}/data_moments.txt, text  append  nomsg 
log off


preserve
	keep if year>=2010
	* by year only  
	bysort year: egen tot_firms = total(N_firms)
	bysort year: egen tot_firmdeaths = total(firmdeath_firms)
	gen firm_exitrate = tot_firmdeaths/tot_firms

	* by year, firm age category
	bysort year fage: egen tot_firms_byage = total(N_firms)
	bysort year fage: egen tot_firmdeaths_byage = total(firmdeath_firms)
	gen firm_exitrate_byage = tot_firmdeaths_byage/tot_firms_byage
	
	keep year fage firm_exitrate firm_exitrate_byage
	duplicates drop
	export delimited using "${moments}/BDS_firmexit" , delim(",") nolabel replace 
	log on
	sum firm_exitrate
	noisily display "exitrate		" r(mean) 
	sum firm_exitrate_byage if fage == 2 
	noisily display "exitrate_age0		" r(mean)
	log off	
restore


preserve
	keep if year>=2010

	* by year, firm size category
	* Four categories: 0-10,10-20,20-100,100+
	recode emp_c (1 2 =1) (3=2) (4=3) (5=4),gen(emp_c_4)
	label define emp_c_4 1 "1-9" 2 "10-19" 3 "20-99" 4 "100-499"
	label value emp_c_4 emp_c_4
	bysort year emp_c_4: egen tot_firms_bysize = total(N_firms)
	bysort year emp_c_4: egen tot_firmdeaths_bysize = total(firmdeath_firms)
	gen firm_exitrate_bysize = tot_firmdeaths_bysize/tot_firms_bysize
	
	keep year emp_c_4 firm_exitrate_bysize
	duplicates drop
	export delimited using "${moments}/BDS_firmexit_size" , delim(",") nolabel replace 
	log on
	sum firm_exitrate if emp_c_4 ==1
	noisily display "exitrate_0_9		" r(mean) 
	sum firm_exitrate if emp_c_4 ==2
	noisily display "exitrate_10_19		" r(mean) 
	sum firm_exitrate if emp_c_4 ==3
	noisily display "exitrate_20_99		" r(mean) 
	sum firm_exitrate if emp_c_4 ==4
	noisily display "exitrate_100_499		" r(mean) 
	log off	
restore

* C.2.b firm entry rate: entry rate = firms of age zero / total num of firms

preserve
	keep if year>=2010
	gen N_newfirms = 0
	replace N_newfirms = N_firms if fage ==1
	* by year only  
	bysort year: egen tot_firms = total(N_firms)
	bysort year: egen tot_newfirms = total(N_newfirms)  
	gen firm_entryrate = tot_newfirms/tot_firms
	
	keep year firm_entryrate 
	duplicates drop
	
	export delimited using "${moments}/BDS_firmentry" , delim(",") nolabel replace 
	log on
	sum firm_entryrate
	noisily display "entryrate		" r(mean) 
	log off	
restore

* C.3 Average firm size by firm age
preserve
	keep if year>=2010
	* all small firms  
	bysort year fage: egen emp = total(tot_emp)
	bysort year fage: egen firms = total(N_firms)
	gen ave_emp = emp/firms

	keep year fage ave_emp
	duplicates drop
	export delimited using "${moments}/BDS_firmsize_age" , delim(",") nolabel replace 
	log on
	forvalues a = 1/6 {
		sum ave_emp if fage==`a'
		local aa = `a' - 1
		noisily display "avefirmsize_age`aa'		" r(mean)
	}
	log off	
	
restore


* C.4 job creation and job destruction rate for continuing establishments and by firm size
preserve
	keep if year>=2010
	drop if fage ==1 /*drop new firms*/
	* by year only  
	bysort year: egen tot_denom = total(denom)
	bysort year: egen tot_job_creation_cont = total(job_creation_continuers)
	bysort year: egen tot_job_destruction_cont = total(job_destruction_continuers)
	gen jcrate_cont = tot_job_creation_cont/tot_denom
	gen jdrate_cont = tot_job_destruction_cont/tot_denom

	* by year and firm size category
	bysort year emp_c: egen tot_denom_c = total(denom)
	bysort year emp_c: egen tot_job_creation_cont_c = total(job_creation_continuers)
	bysort year emp_c: egen tot_job_destruction_cont_c = total(job_destruction_continuers)
	gen jcrate_cont_c = tot_job_creation_cont_c/tot_denom_c
	gen jdrate_cont_c = tot_job_destruction_cont_c/tot_denom_c

	keep year emp_c jcrate_cont jdrate_cont jcrate_cont_c jdrate_cont_c
	order year emp_c jcrate_cont jdrate_cont jcrate_cont_c jdrate_cont_c
	duplicates drop
	export delimited using "${moments}/BDS_jobcrcontuing" , delim(",") nolabel replace 
	
	log on
	sum jcrate_cont
	noisily display "jcr		" r(mean) 
	sum jdrate_cont
	noisily display "jdr		" r(mean) 
	log off	
	
	
	
restore

* C.5 Weighted average of KFS moments (weighted by firm age distribution)
* Data in C.5 come from the Kauffman Firm Survey. See code main_KFS.do.
preserve
	keep if year>=2010 & fage<12
	* Import debt-to-asset ratio by firm age from KFS results mannually
	capture drop daratio
	gen daratio = .
	replace daratio = 0.1915053 if fage == 1 /* age 0 */
	replace daratio = 0.1840776 if fage == 2 /* age 1 */
	replace daratio = 0.0863333 if fage == 3 /* age 2 */
	replace daratio = 0.1069989 if fage == 4 /* age 3 */
	replace daratio = 0.344448 if fage == 5 /* age 4 */
	replace daratio = 0.055098 if fage == 6 /* age 5 */
	replace daratio = (0.055098+0.0441066+0.0083223)/3 if fage >6 & fage<. /* age 6+ */
	
	* Import frac of firms with net debt by firm age from KFS results mannually
	capture drop hasDebt
	gen hasDebt = .
	replace hasDebt = 0.4202558 if fage == 1 /* age 0 */
	replace hasDebt = 0.3733265 if fage == 2 /* age 1 */
	replace hasDebt = 0.3727061 if fage == 3 /* age 2 */
	replace hasDebt = 0.3598956 if fage == 4 /* age 3 */
	replace hasDebt = 0.3854402 if fage == 5 /* age 4 */
	replace hasDebt = 0.3278793 if fage == 6 /* age 5 */
	replace hasDebt = (0.3278793+0.294475+0.2802074)/3 if fage >6 & fage<. /* age 6+ */
	
	* Import debt to payroll ratio conditional on having debt by firm age from KFS results mannually
	capture drop debt_payroll_cond
	gen debt_payroll_cond = .
	replace debt_payroll_cond = 2.939828 if fage == 1 /* age 0 */
	replace debt_payroll_cond = 4.169449 if fage == 2 /* age 1 */
	replace debt_payroll_cond = 0.5912598 if fage == 3 /* age 2 */
	replace debt_payroll_cond = 1.387434 if fage == 4 /* age 3 */
	replace debt_payroll_cond = 4.331478 if fage == 5 /* age 4 */
	replace debt_payroll_cond = 0.9035707 if fage == 6 /* age 5 */
	replace debt_payroll_cond = (0.9035707+0.796663+0.5896994)/3 if fage >6 & fage<. /* age 6+ */
	
	* Import cash to payroll ratio conditional on having no net debt by firm age from KFS results mannually
	capture drop cash_payroll_cond
	gen cash_payroll_cond = .
	replace cash_payroll_cond = 1.061513 if fage == 1 /* age 0 */
	replace cash_payroll_cond = 1.025575 if fage == 2 /* age 1 */
	replace cash_payroll_cond = 0.2601712 if fage == 3 /* age 2 */
	replace cash_payroll_cond = 1.695405 if fage == 4 /* age 3 */
	replace cash_payroll_cond = 0.4012887 if fage == 5 /* age 4 */
	replace cash_payroll_cond = 4.53265 if fage == 6 /* age 5 */
	replace cash_payroll_cond = (4.53265+0.8841256+0.1238163)/3 if fage >6 & fage<. /* age 6+ */


	* Import capital to payroll  ratio by firm age from KFS results mannually
	capture drop k_w_ratio
	gen k_w_ratio = .
	replace k_w_ratio = 7.838887 if fage == 1
	replace k_w_ratio = 6.340726 if fage == 2
	replace k_w_ratio = 2.920279 if fage == 3
	replace k_w_ratio = 4.324912 if fage == 4
	replace k_w_ratio = 5.099328 if fage == 5
	replace k_w_ratio = 3.739924 if fage == 6
	replace k_w_ratio = (3.739924+4.236708+4.542665)/3 if fage >6 & fage<.
	
	* Import employment autocorrelation by firm age from KFS results mannually
	capture drop empcorr
	gen empcorr = .
	replace empcorr = 0.6372023 if fage == 1
	replace empcorr = 0.81197774 if fage == 2
	replace empcorr = 0.76093911 if fage == 3
	replace empcorr = 0.79982172 if fage == 4
	replace empcorr = 0.92514354 if fage == 5
	replace empcorr = 0.90937235 if fage == 6
	replace empcorr = (0.90937235+0.9338262)/2 if fage >6 & fage<.
	
	
	* Calculate average KFS moments after weighting
	
	
	log on
	* debt-asset-ratio
	sum daratio [aw = N_firms],d	
	noisily display "debt_asset_all		" r(mean) 
	sum daratio [aw = N_firms] if fage==1,d	
	noisily display "debt_asset_entrants		" r(mean) 
	* has net debt
	sum hasDebt [aw = N_firms],d	
	noisily display "hasNetDebt		" r(mean) 
	* has net debt entrants
	sum hasDebt [aw = N_firms] if fage == 1,d	
	noisily display "hasNetDebt_age0		" r(mean) 
	* debt_payroll_cond
	sum debt_payroll_cond [aw = N_firms],d	
	noisily display "debt_payroll_cond		" r(mean) 
	* cash_payroll_cond
	sum cash_payroll_cond [aw = N_firms],d	
	noisily display "cash_payroll_cond		" r(mean) 
	* capital-payroll ratio
	sum k_w_ratio [aw = N_firms],d	
	noisily display "k_payroll_ratio		" r(mean) 
	* autocorrelation of employment
	sum empcorr [aw = N_firms],d	
	noisily display "autocorr_emp		" r(mean) 
	log off	
restore	


log close
}



**************************************************************
* 	D.1.	Read and compute fixed expenses 	
* Table prepared by washington post.
*	D.2. 	Moments from Khan and Thomas (2013 JPE) 
**************************************************************
* 	D.1.	Read and compute fixed expenses 	
import delim using "${other}/business_expenses.txt",delim(",") varnames(1)  clear
foreach var in overhead rent cost_of_sales payroll {
	replace `var' = subinstr(`var', "%", "",.) 
	destring `var',replace
	replace `var' = `var'/100
}
compress
save ${work}/business_expenses.dta,replace

use ${work}/business_expenses.dta,clear
* gen fixed cost = overhead - payroll
gen fixed = overhead-payroll
gen fixed_payroll = fixed/payroll
quietly{

	cap log close                                           
	log using ${moments}/data_moments.txt, text  append  nomsg 
	* compute average of fixed cost weighted by yearly revenue
	sum fixed [aw = yearly_revenue]
	noisily display "fixedcost_to_rev		" r(mean) 
	sum fixed_payroll [aw = yearly_revenue]
	noisily display "fixedcost_to_payroll		" r(mean) 
	log close
}

*	D.2. 	Moments from Khan and Thomas (2013 JPE) 
quietly{

	cap log close                                           
	log using ${moments}/data_moments.txt, text  append  nomsg 
	noisily display "stdev_invrate		" 0.337 
	noisily display "scor_invrate		" 0.058 
	noisily display "freq_lumpinv		" 0.186 
	noisily display "mean_invrate		" 0.122 

	log close
}

**************************************************************
************* 	E.	Read and compute impact on small firm output using data from Bloom et al (2021) 	
**************************************************************

import delim using "${other}/bloom_etal_sales.txt",delim(",") varnames(1) clear
save ${work}/bloom_etal_sales.dta,replace
* Get weights from SUSB

use ${work}/SUSB_empdist.dta,clear
drop if inlist(lfo,1,6,7,8)
keep if naics == 0
keep if inlist(emp_c,2,3,4,6,7)
recode emp_c (2=1) (3=2) (4=3) (6 7 = 4),gen(size)
* Weight based on emp
bysort year size: egen weight_emp = total(N_firm)
bysort year: egen N_firm_all = total(N_firm)
replace weight_emp  = weight_emp /N_firm_all
* weight based on sales
bysort year size: egen weight_sales = total(tot_receipts)
bysort year: egen tot_receipts_all = total(tot_receipts)
replace weight_sales = weight_sales/tot_receipts_all
keep year size weight*
duplicates drop
keep if year == 2012
drop year
* merge weight to bloom_etal_sales.dta
merge 1:m size using ${work}/bloom_etal_sales.dta
drop _merge
save ${work}/bloom_etal_sales.dta,replace


* Quarterly drop in sales, 100 = 2020Q1
use  ${work}/bloom_etal_sales.dta,clear
gen rel_sales = 1-drop_in_sales
sort size quarter
by size: gen rel_sales_q1 = rel_sales/rel_sales[1]
gen drop_rel_q1 = rel_sales_q1-1
* average relative drop in sales, weighted
cap log close                                           
log using ${moments}/output_small.txt, text replace   nomsg 
tabstat  drop_rel_q1 [aw = weight_sales] ,by(quarter)
log close



**************************************************************
*********** 	F.	Plot Response to the Pandemic  ***********
**************************************************************


**** Load data
import delim using "${FRED}/fredgraph.csv",delim(",") varnames(noname) ///
	colrange(1:8) rowrange(2:7)  clear

* date
capture drop date
gen date = date(v1,"MDY",2050)
format date %td
gen quarter = qofd(date)
format quarter %tq

* year 
gen year = year(date)

* change variable names
ren v2 consumption
ren v3 investment
ren v4 employment
ren v5 GDP
ren v6 output
ren v7 exit_rate
ren v8 entry_rate

* gen yearly averages
foreach var in GDP output consumption investment employment exit_rate entry_rate {
bysort year: egen `var'_y = mean(`var')
}

* Plot 1: output, consumption, investment, employment
set scheme s2color
twoway line output consumption investment employment quarter if quarter>=241 & quarter<=244, ///
	lp("l" "_" "-" "_-") lw("thick" "thick" "thick" "thick") xlabel(241(1)244) ///
	plotregion(fcolor(white)) graphregion(fcolor(white)) ///
	ytitle("Index (Pre-pandemic = 100)") xtitle("") name(fred,replace) ///
	legend(order(1 "Output" 2 "Consumption" 3 "Private Investment" 4 "Employment-Population Ratio") ) 
graph export ${moments}/FRED_Graphs.eps,name(fred) replace
graph export ${moments}/FRED_Graphs.png,name(fred) replace

* Plot 2: Entry and Exit
set scheme s2color
twoway line entry_rate exit_rate quarter if quarter>=241 & quarter<=244, ///
	lp("l" "_") lw("thick" "thick")  xlabel(241(1)244) ///
	plotregion(fcolor(white)) graphregion(fcolor(white)) ///
	ytitle("Index (Pre-pandemic = 100)") xtitle("") name(fred,replace) ///
	legend(order(1 "Entry Rate" 2 "Exit Rate") ) 
graph export ${moments}/FRED_Graphs_ee.eps,name(fred) replace
graph export ${moments}/FRED_Graphs_ee.png,name(fred) replace


