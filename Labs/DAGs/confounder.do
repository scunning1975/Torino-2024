********************************************************************************
* name: confounder.do
* author: scott cunningham (baylor)
* description: simple depiction of confounder
* last updated: june 4, 2024
********************************************************************************

clear
capture log close
set seed 10501

set obs 10000

* Confounder triangle: 
* D -> Y
* X -> D
* X -> Y

* Generate Covariates (Baseline values)
  gen x = rnormal(35, 10)
  sum x
  replace x = x-`r(mean)'

* Treatment probability increases with age and decrease with gpa
  gen prob = 0.3 + 0.3 * (x > 0) 
  gen treat = runiform() < prob

* Data generating process for error  
  gen          e = rnormal(0, 5)
  gen     y0 = 100 + 100 * x + e 

* Covariate-based treatment effect with constant treatment effects
  gen         y1 = y0 + 1000 
  
* Aggregate causal parameters
  gen delta = y1-y0
  sum delta, meanonly
  gen ate = `r(mean)'
  su delta if treat==1, meanonly
  gen att = `r(mean)'
  su delta if treat==0, meanonly
  gen atu = `r(mean)'
  
  su ate att atu // constant treatment effects

* Switching equation   
  gen earnings = treat*y1 + (1-treat)*y0
  
* Estimation with and without confounder controls  
  reg earnings treat, robust
  reg earnings treat x, robust

capture log close
exit


  
  
  