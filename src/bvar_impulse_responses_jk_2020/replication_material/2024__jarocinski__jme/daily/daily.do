* Prepare data for local projection estimation

* FRED from text
import delimited using daily_Daily.txt, varnames(1) clear
rename date temp
gen date = date(temp,"YMD"), before(temp)
format date %td
rename dgs3mo gs3mo
keep date gs3mo
save fred_daily1.dta, replace

import delimited using daily_Daily_7_Day.txt, varnames(1) clear
rename date temp
gen date = date(temp,"YMD"), before(temp)
format date %td
rename dff ff
keep date ff
save fred_daily2.dta, replace

import delimited using "daily_Daily_Close.txt", varnames(1) clear
rename date temp
gen date = date(temp,"YMD"), before(temp)
format date %td
rename bamlc0a2caa bofaml_us_aa_oas
rename bamlc0a2caaey bofaml_us_aa_yld
rename bamlc0a4cbbb bofaml_us_bbb_oas
rename bamlc0a4cbbbey bofaml_us_bbb_yld
rename bamlh0a0hym2 bofaml_us_hyld_oas
rename bamlhe00ehyioas bofaml_ea_hyld_oas
keep date bofaml_us*
save fred_daily3.dta, replace

import delimited using daily_Weekly_Ending_Friday.txt, varnames(1) clear
rename date temp
gen date = date(temp,"YMD"), before(temp)
format date %td
tsset date
tsfill
ipolate nfci date, gen(nfci_lin)
keep date nfci_lin
save fred_daily4.dta, replace

* BIS EUR per USD, from 1974, the same as the ECB official rate from 1999
* https://data.bis.org/topics/XRU/BIS,WS_XRU,1.0/D.XM.EUR.A
import delimited using bis_dp_search_export_20230915-113252.csv, rowrange(5) varnames(nonames) colrange(7:8) clear
gen date = date(v1,"YMD"), before (v1)
format date %td
gen eurusd = 100*ln(v2)
keep date eurusd
save bis_eurusd.dta, replace


* SP500 (from Haver)
import excel using haver_stock.xlsx, sheet("data") first clear
gen date = date(Jan,"DMY",2050), before(Jan)
format date %td
drop Jan
drop if year(date)<1970
drop if missing(date)
gen sp500 = 100*ln(SP500DAILY)
keep date sp500
save haver_stock.dta, replace

* Zero coupon yields by Gurkaynak, Sack, Wright (2006) from
* https://www.federalreserve.gov/data/nominal-yield-curve.htm
import delimited using feds200628.csv, rowrange(11) varnames(10) clear
rename date datestr
gen date = date(datestr, "YMD"), before(datestr)
format date %td
drop datestr
keep date sveny01 sveny02 sveny05 sveny10
destring, replace force
save feds200628.dta, replace

* TIPS curve from Gurkaynak, Sack and Wright (2008)
* https://www.federalreserve.gov/pubs/feds/2008/200805/200805abs.html
import delimited using "feds200805.csv", rowrange(20) varnames(19) clear
rename date datestr
gen date = date(datestr, "YMD"), before(datestr)
format date %td
drop datestr
keep date bkeven02-bkeven05 bkeven10 bkeven1f4 bkeven1f9 bkeven5f05
destring, replace force
save feds200805.dta, replace

* Broad, AFE, EME Nominal Dollar Indices extended back to 1973
use broad/broad.dta, clear
replace broad = 100*ln(broad)
replace afe = 100*ln(afe)
replace eme = 100*ln(eme)
save broad.dta, replace

* Shadow Fed funds rate by Leo Krippner from https://www.ljkmfa.com/
import excel using "SSR_Estimates_20210726.xlsx", sheet("D. Daily SSR series") cellrange(A7) case(l) firstrow clear
gen date = date(a, "DMY"), before(a)
format date %td
drop if missing(date)
rename usssr krippner
keep date krippner
save SSR_Estimates_20210726.dta, replace


* merge
use fred_daily1, clear
erase fred_daily1.dta
merge 1:1 date using fred_daily2.dta, nogenerate
erase fred_daily2.dta
merge 1:1 date using fred_daily3.dta, nogenerate
erase fred_daily3.dta
merge 1:1 date using fred_daily4.dta, nogenerate
erase fred_daily4.dta
merge 1:1 date using bis_eurusd.dta, nogenerate
erase bis_eurusd.dta
merge 1:1 date using haver_stock.dta, nogenerate
erase haver_stock.dta
merge 1:1 date using feds200628.dta, nogenerate
erase feds200628.dta
merge 1:1 date using feds200805.dta, nogenerate
erase feds200805.dta
merge 1:1 date using broad.dta, nogenerate
erase broad.dta
merge 1:1 date using SSR_Estimates_20210726.dta, nogenerate
erase SSR_Estimates_20210726.dta

sort date
keep if year(date)>=1990

drop if dow(date)==0 | dow(date)==6
drop if month(date)==1 & day(date)==1
drop if month(date)==12 & day(date)==25

gen krippnerlong = krippner
replace krippnerlong = ff if date<td(3jan1995)

export delimited using daily_raw.csv, replace

quietly ds date, not
local all_vars `r(varlist)'

* fill missing values
foreach var in `all_vars' {
replace `var' = `var'[_n-1] if `var'==. & `var'[_n+2]!=.
}

* transformations
foreach var in `all_vars' {
	*gen tl1`var' = `var'[_n] - `var'[_n-1]
	gen f1l1`var' = `var'[_n+1] - `var'[_n-1]
	gen f2l1`var' = `var'[_n+2] - `var'[_n-1]
	gen f3l1`var' = `var'[_n+3] - `var'[_n-1]
	gen f4l1`var' = `var'[_n+4] - `var'[_n-1]
	gen f5l1`var' = `var'[_n+5] - `var'[_n-1]
	gen f10l1`var' = `var'[_n+10] - `var'[_n-1]
	gen f15l1`var' = `var'[_n+15] - `var'[_n-1]
	gen f20l1`var' = `var'[_n+20] - `var'[_n-1]
	gen f25l1`var' = `var'[_n+25] - `var'[_n-1]
}

keep if year(date)<=2019

export delimited using daily.csv, replace


