********************************************************************************
* name: matching.do
* description: illustrating nearest neighbor matching with and without bias correction/adjustment
* author: scott cunningham (baylor)
* last updated: june 3, 2024
********************************************************************************

    clear 
	capture log close
    drop _all 
	set seed 5000
	set obs 5000
	gen 	treat = 0 
	replace treat = 1 in 2501/5000
	
	* Poor pre-treatment fit
	gen 	age = rnormal(25,2.5) 		if treat==1
	replace age = rnormal(30,3) 		if treat==0
	gen 	gpa = rnormal(2.3,0.75) 	if treat==0
	replace gpa = rnormal(1.76,0.5) 	if treat==1
	
	twoway (histogram age if treat==1,  color(green)) ///
       (histogram age if treat==0,  ///
	   fcolor(none) lcolor(black)), legend(order(1 "Treated" 2 "Not treated" ))
	   
	twoway (histogram gpa if treat==1,  color(blue)) ///
       (histogram gpa if treat==0,  ///
	   fcolor(none) lcolor(black)), legend(order(1 "Treated" 2 "Not treated" ))

	su age
	replace age = age - `r(mean)'

	su gpa
	replace gpa = gpa - `r(mean)'

	* All combinations 
	gen age_sq 		= age^2
	gen gpa_sq 		= gpa^2
	gen interaction	= gpa*age
	gen agegpa		= age*gpa	 

	gen y0 = 15000 + 10.25*age + -10.5*age_sq + 1000*gpa + -10.5*gpa_sq + 500*interaction + rnormal(0,5)
	gen y1 = y0 + 2500 + 100 * age + 1000 * gpa
	gen delta = y1 - y0

	su delta // ATE = 2500
	su delta if treat==1 // ATT = 1977
	local att = r(mean)
	scalar att = `att'
	gen att = `att'

	gen earnings = treat*y1 + (1-treat)*y0
	
********************************************************************************
* Estimation
********************************************************************************
	
	
	** Baseline regression model that assumes constant treatment effects
	regress earnings age gpa age_sq gpa_sq agegpa treat, robust

	** Nearest neighbor (mahanalobis distance minimization) matching (ATE, ATT)
	* Matching model 1 for ATE
	teffects nnmatch (earnings age gpa age_sq gpa_sq agegpa) (treat), ate nn(1) metric(maha) 

	* Matching model 2 for ATT
	teffects nnmatch (earnings age gpa age_sq gpa_sq agegpa) (treat), atet nn(1) metric(maha) 
	
	* Matching Model 3 for ATE with bias adjustment
	teffects nnmatch (earnings age gpa age_sq gpa_sq agegpa) (treat), ate nn(1) metric(maha) biasadj(age age_sq gpa gpa_sq agegpa)

	* Matching model 4 for ATT with bias adjustment
	teffects nnmatch (earnings age gpa age_sq gpa_sq agegpa) (treat), atet nn(1) metric(maha) biasadj(age age_sq gpa gpa_sq agegpa)
	
** Suggestion: remove covariates in bias adjustment formulas for the outcome regression model to help yourself better understand the importance of model misspecification in the bias adjustment 	

	
  ** IPW method for estimating ATT
  logit treat age age_sq gpa gpa_sq interaction
  
  * predict propensity score
  predict pscore

  * inverse propensity score weights (ATT)
  gen inv_ps_weight = treat + (1-treat) * pscore/(1-pscore)

  * Estimation with IPW weights
  reg earnings i.treat [aw=inv_ps_weight], robust
  su att
  
  * Trim
  drop if pscore<0.1 | pscore>0.9
  reg earnings i.treat [aw=inv_ps_weight], robust
  su att
  su delta if treat==1 // parameter changed
  
  