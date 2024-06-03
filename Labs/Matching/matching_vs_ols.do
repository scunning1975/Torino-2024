********************************************
* name: matching_vs_ols.do
* author: scott cunningham (baylor)
* description: illustrating bias in matching and OLS and their respective corrections
* last updated: april 4 2024
********************************************

	clear
	capture log close
    drop _all 
	set seed 100
	set obs 5000
	gen 	treat = 0 
	replace treat = 1 in 2501/5000
	
	* Poor pre-treatment fit
	gen 	age = rnormal(40,2.5) 		if treat==1
	replace age = rnormal(27,3) 		if treat==0
	gen 	gpa = rnormal(2.2,0.75) 	if treat==0
	replace gpa = rnormal(1.6,0.5) 		if treat==1

	* Visualize the imbalance
	twoway (histogram age if treat==1,  color(green)) ///
       (histogram age if treat==0,  ///
	   fcolor(none) lcolor(black)), legend(order(1 "Treated" 2 "Not treated" ))

	twoway (histogram gpa if treat==1,  color(blue)) ///
       (histogram gpa if treat==0,  ///
	   fcolor(none) lcolor(black)), legend(order(1 "Treated" 2 "Not treated" ))

	* Recenter the covariates (needed to generate potential outcomes)
	su age
	replace age = age - `r(mean)'

	su gpa
	replace gpa = gpa - `r(mean)'

	* Polynomials and interactions
	gen age_sq 		= age^2
	gen gpa_sq 		= gpa^2
	gen interaction	= gpa*age

	* Generate the potential outcomes
	gen y0 = 15000 + 10.25*age + -10.5*age_sq + 1000*gpa + -10.5*gpa_sq + 500*interaction + rnormal(0,5)
	label variable y0 "Earnings if you do not get a masters"
	
	gen y1 = y0 + 2500 + 100 * age + 4000 * gpa
	label variable y1 "Earnings if you get a masters"

	* Individual treatment effects
	gen delta = y1 - y0
	label variable delta "Causal effect of getting a masters on earnings"

	* Aggregate causal parameters
	su delta // ATE = 2500
	su delta if treat==1 // ATT = 1936.37
	local att = r(mean)
	scalar att = `att'
	gen att = `att'

	* Use the switching equation to assign y1 or y0 to a unit based on their treatment status
	gen earnings = treat*y1 + (1-treat)*y0

	* Nearest neighbor matching without bias adjustment
	teffects nnmatch (earnings age age_sq gpa gpa_sq interaction) (treat), atet nn(1) metric(euclidean) gen(match)
	
	* Nearest neighbor matching with bias adjustment
	teffects nnmatch (earnings age age_sq gpa gpa_sq interaction) (treat), atet nn(1) metric(euclidean) biasadj(age age_sq gpa gpa_sq interaction)
	
	* OLS regressions - this regression assumes "constant treatment effects" and the coefficient on treat is trying to estimate the ATE. 
	regress earnings age age_sq gpa gpa_sq interaction treat, robust
	
	* Regression adjustment
	teffects ra (earnings age age_sq gpa gpa_sq interaction) (treat), atet
	
	
	* Regression adjustment "the long way". At its core, regression adjustment "interacts" the treatment variable with all of the covariates. 	
	
	* First step is run the fully interacted regression.
	#delimit ;
	
	regress earnings 	i.treat##c.age 
						i.treat##c.age_sq
						i.treat##c.gpa 
						i.treat##c.gpa_sq					
						i.treat##c.age##c.gpa;
	#delimit cr					
	
	local ate2=_b[1.treat]
	scalar ate2 = `ate2'
	gen ate2=`ate2'
	
	* Second step is save each of those treatment coefficients. 
	
	local treat_coef 		= _b[1.treat] // 1
	local age_treat_coef 	= _b[1.treat#c.age] // 2
	local agesq_treat_coef 	= _b[1.treat#c.age_sq] // 3
	local gpa_treat_coef 	= _b[1.treat#c.gpa] // 4
	local gpasq_treat_coef 	= _b[1.treat#c.gpa_sq] // 5
	local age_gpa_coef 		= _b[1.treat#c.age#c.gpa] // 6
	
	scalar 	treat_coef = `treat_coef'
	gen 	treat_coef_var = `treat_coef' // 1

	scalar 	age_treat_coef = `age_treat_coef'
	gen 	age_treat_coef_var = `age_treat_coef' // 2
	
	scalar 	agesq_treat_coef = `agesq_treat_coef'
	gen 	agesq_treat_coeff_var = `agesq_treat_coef' // 3

	scalar 	gpa_treat_coef = `gpa_treat_coef'
	gen 	gpa_treat_coef_var = `gpa_treat_coef' // 4

	scalar 	gpasq_treat_coef = `gpasq_treat_coef'
	gen 	gpasq_treat_coef_var = `gpasq_treat_coef' // 5

	scalar 	age_gpa_coef = `age_gpa_coef'
	gen 	age_gpa_coef_var = `age_gpa_coef' // 6
		
	* Step three will now calculate the average value of each of the covariates (for example age or gpa or the interaction of age and gpa) for the treatment group only
	
	su 	age if treat==1
	local mean_age = `r(mean)'
	gen mean_age = `mean_age'
	
	su age_sq if treat==1
	local mean_agesq = `r(mean)'
	gen mean_agesq = `mean_agesq'
	
	su gpa if treat==1
	local mean_gpa = `r(mean)'
	gen mean_gpa = `mean_gpa'
	
	su gpa_sq if treat==1
	local mean_gpasq = `r(mean)'
	gen mean_gpasq = `mean_gpasq'
	
	su interaction if treat==1
	local mean_agegpa = `r(mean)'
	gen mean_agegpa = `mean_agegpa'

	* Step four is calculate or estimate the ATT using all that information (coefficient x sample means for treatment group added together)
	
gen treat4 = 	treat_coef_var + /// 1
                age_treat_coef_var * mean_age + /// 2
                agesq_treat_coeff_var * mean_agesq + /// 3
                gpa_treat_coef_var * mean_gpa + /// 4
                gpasq_treat_coef_var * mean_gpasq + /// 5
                age_gpa_coef_var * mean_agegpa  /// 6
				
sum delta if treat==1
sum treat4 // regression adjustment "the long way"

	* Regression adjustment "the medium way" (Kitagawa-Oaxaca-Blinder method)
	
	regress earnings age age_sq gpa gpa_sq c.age##c.gpa if treat == 0
	predict y0_hat, xb
 
	regress earnings age age_sq gpa gpa_sq c.age##c.gpa if treat == 1
	predict y1_hat, xb
 
	gen te_hat = y1_hat - y0_hat
	sum te_hat if treat == 1
	
	* Regression adjustment using teffects ("the short way")
	teffects ra (earnings age age_sq gpa gpa_sq interaction) (treat), atet
				
capture log close
exit

