* broad afe eme (back to 2006): https://www.federalreserve.gov/datadownload/Build.aspx?rel=H10
import delimited using "FRB_H10.csv", rowrange(7) varnames(6) clear
gen date = date(timeperiod, "YMD"), before(timeperiod)
format date %td
drop timeperiod
destring, replace force
rename jrxwtfb_nb broad
rename jrxwtfn_nb afe
rename jrxwtfo_nb eme
save "FRB_H10.dta", replace

* data back to 1973:
* https://www.federalreserve.gov/econres/notes/feds-notes/revisions-to-the-federal-reserve-dollar-indexes-20190115.html
* https://www.federalreserve.gov/econres/notes/ifdp-notes/IFDP_Note_Data_Appendix.xlsx
import excel using IFDP_Note_Data_Appendix.xlsx, sheet("Nominal Dollar Indexes Daily") firstrow case(lower) clear
rename period date
rename nominalbroaddaily broad_h
rename nominalafedaily afe_h
rename nominalemedaily eme_h
destring, replace force

* merge
merge 1:1 date using "FRB_H10.dta", nogenerate
erase "FRB_H10.dta"
sort date

* are they similar? yes
//graph twoway (line broad* date) if date>td(02jan2006)

replace broad = broad_h if date<td(02jan2006)
drop broad_h
replace afe = afe_h if date<td(02jan2006)
drop afe_h
replace eme = eme_h if date<td(02jan2006)
drop eme_h

save broad_afe_eme.dta, replace
