clear all

* Set up the results file
tempname handle
postfile `handle' ate att atu ols prob_treat prob_control n tymons_ate tymons_att tymons_atu w1 w0 delta using results.dta, replace

* Loop through the iterations
forvalues i = 4/4996 {
    clear
    set seed 5150
    set obs 5000
    gen treat = 0 
    replace treat = 1 in `i'/5000

    * Imbalanced covariates
    gen 	age = rnormal(25, 2.5) if treat == 1
    replace age = rnormal(30, 3) if treat == 0
    gen 	gpa = rnormal(2.3, 0.75) if treat == 0
    replace gpa = rnormal(1.76, 0.5) if treat == 1

    * Re-center the covariates
    su age, meanonly
    replace age = age - r(mean)

    su gpa, meanonly
    replace gpa = gpa - r(mean)

    * Quadratics and interaction
    quietly gen age_sq = age^2
    quietly gen gpa_sq = gpa^2
    quietly gen interaction = gpa * age

    * Modeling potential outcomes
    quietly gen y0 = 15000 + 10.25*age + -10.5*age_sq + 1000*gpa + -10.5*gpa_sq + 500*interaction + rnormal(0, 5)
    quietly gen y1 = y0 + 2500 + 100 * age + 1100 * gpa
    quietly gen treatment_effect = y1 - y0

    * Calculate ATE, ATT, and ATU
    su treatment_effect, meanonly
    local ate = r(mean)
    su treatment_effect if treat == 1, meanonly
    local att = r(mean)
    su treatment_effect if treat == 0, meanonly
    local atu = r(mean)

    * Generate earnings variable
    quietly gen earnings = treat * y1 + (1 - treat) * y0

	* Get the weights, OLS coefficient and ALPE 
	quietly hettreatreg age gpa age_sq gpa_sq interaction, o(earnings) t(treat) vce(robust)
	local ols `e(ols1)'
	local prob_treat `e(p1)'
	local prob_control `e(p0)'
	local n `e(N)'
	local tymons_ate `e(ate)'
	local tymons_att `e(att)'
	local tymons_atu `e(atu)'
	local w1 `e(w1)'
	local w0 `e(w0)'
	local delta `e(delta)'
	

    * Post the results to the results file
    post `handle' (`ate') (`att') (`atu') (`ols') (`prob_treat') (`prob_control') (`n') (`tymons_ate') (`tymons_att') (`tymons_atu') (`w1') (`w0') (`delta') 

	}

* Close the postfile
postclose `handle'

* Use the results
use results.dta, clear

gsort prob_treat
gen id=_n

* Simple graph
twoway (line w1 prob_treat), ytitle("OLS weight on ATT") xtitle("Share of units in treatment group") title("Weighted Average Interpretation of Sloczynsi OLS Theorem") subtitle("Treatment Group Size vs. ATT Weights") note("5000 OLS regressions varying treatment share from 0.001 to 0.9994") legend(position(6))

#delimit cr


