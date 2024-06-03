* Assignment is not based on the propensity score; rather it is two separate samples with different imbalances in covariates which is oddly enough not the same
    clear 
    drop _all 
	set obs 5000
	gen 	treat = 0 
	replace treat = 1 in 2500/5000

	ssc install hettreatreg, replace
	
	* Poor pre-treatment fit
	gen 	age = rnormal(25,2.5) 		if treat==1
	replace age = rnormal(30,3) 		if treat==0
	gen 	gpa = rnormal(2.3,0.75) 	if treat==0
	replace gpa = rnormal(1.76,0.5) 	if treat==1

	su age
	replace age = age - `r(mean)'

	su gpa
	replace gpa = gpa - `r(mean)'

	* All combinations 
	gen age_sq 		= age^2
	gen gpa_sq 		= gpa^2
	gen interaction	= gpa*age
	gen agegpa		= age*gpa	 

	gen y0 = 15000 + 10.25*age + -10.5*age_sq + 1000*gpa + -10.5*gpa_sq + 2000*interaction + rnormal(0,5)
	gen y1 = y0 + 2500 + 100 * age + 1000*gpa
	gen delta = y1 - y0

	su delta 			 // ATE = 2500
	su delta if treat==1 // ATT = 1991
	su delta if treat==0 // ATT = 3009

	* Switching equation
	gen earnings = treat*y1 + (1-treat)*y0

	* OLS with additive controls under exogeneity
	reg earnings treat age gpa age_sq gpa_sq interaction, robust

	su delta
	su delta if treat==1
	su delta if treat==0

	* Decomposition of OLS using Tymon's OLS theorem
	hettreatreg age gpa age_sq gpa_sq interaction, o(earnings) t(treat)

	* Or use RA in teffects to estimate the ATT

	teffects ra (earnings age gpa age_sq gpa_sq interaction) (treat), atet
	su delta if treat==1


	* Or use nearest neighbor in teffects to estimate the ATT

	teffects nnmatch (earnings age gpa age_sq gpa_sq interaction) (treat), biasadj(age gpa age_sq gpa_sq interaction) atet
	su delta if treat==1


	* Or use propensity score matching in teffects to estimate the ATT

	teffects psmatch (earnings) (treat age gpa age_sq gpa_sq interaction), atet
	su delta if treat==1
	
