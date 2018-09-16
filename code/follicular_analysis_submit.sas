dm 'log;clear;output;clear;';

libname data "C:\akira\data";

/*create necessary indicator variables
for event of interest and competing risk*/
/*create agegroup dichotomized at 65*/
/*event of interest: failure from disease, either no response to
treatment or relapse*/
/*competing risk: death without failure*/
/*time to event: dftime*/

data follicular;
	set data.follicular;
	a = (age > 65)+ 0;
	/*we name evcens as the indicator for event of interest*/
	if resp = 'NR' or relsite ^= ' ' then evcens = 1;
	else evcens = 0;
	/*we name crcens as the indicator for competing risk*/
	if resp = 'CR' and relsite = ' ' and stat = 1 then crcens = 1;
	else crcens = 0;
	/*we also define a general indicator for all types of events*/
	/*cens = 0(cesoring), 1(event of interest), 2(competing risk)*/
	cens = evcens + 2*crcens;
run;

/*estimates and tests without covariates*/

/*we ask whether the type-specific hazard functions are the same for all
event types*/

/*first a quick look at the frequencies*/
title "a quick look at the event freuencies";
proc freq data=follicular;
	table cens;
run;

/*a formal test*/
/*we have 541 - 193 = 348 all types of events*/
/*so the expected value given equal type-specific hazard function is
348/2 = 174*/
/*a Pearson chi square statistic is then (272 - 174)^2/174  with degree of freedom 1*/
/*p value <0.0001 and hence reject the null that type-specific hazard function is the same
for different events*/

/*examine if the type specific hazard functions are proportional*/
data failure;
	set follicular;
	event = (cens = 1);
	type = 1;
run;

data death;
	set follicular;
	event = (cens = 2);
	type = 2;
run;

data combine;
	set failure death;
run;

title "examine if the type specific hazard functions are proportional";
proc lifetest data=combine plots = lls;
	time dftime*event(0);
	strata type;
run;	

/*comments to be deleted*/
/*curve for death is lower than the curve for failure, i.e.
the event time to death is longer than the failure*/


/*Examine the smoothed hazard plots using the kernel smoothing option*/
title "examine the smoothed hazard plots using the kernel smoothing option";
proc lifetest data=combine plots = H(BW = 10);
	time dftime*event(0);
	strata type;
run;

/*comments to be deleted*/
/*the sharp increase in the hazard for death can be disregarded, since standard error becomes large
at larger points in time*/
/*as expected, the hazard for death before relapse is lower than the failure to treatment*/

/*a parametric test proposed by Cox and Oakes(1984) for testing the proportional hazard hypothesis*/
title "a parametric test proposed by Cox and Oakes for testing the proportional hazard hypothesis";
proc logistic data=follicular;
	where cens NE 0;
	/*only two event types, do not really need glogit option here*/
	model cens = dftime / link = glogit;
run;

/*comments to be deleted*/
/*effect of event time is highly significant, indicating a rejection of the proportionality hypothesis*/
/*it tells us that the hazard for treatment failur increases more slowly with time than
the hazard for death before relapse. Specifically the ratio decrease by about 100(1 - exp(-0.1817)) = 16.7%
each year. This matches the plot above*/


/*Covariate effects via cox models*/

/*look at whether the effects of covariates are the same or different across event types by fitting a Cox
model to each type.*/

data follicular2;
	set follicular;
	if ch = 'Y' then chemo = 1;
	else chemo = 0;
/*in fact all rt = 'Y' so we could not test its effect*/
/*	if rt = 'Y' then radio = 1;*/
/*	else radio = 0;*/
run;

/*using age as a continuous covariate*/
/*proc phreg data=follicular2;*/
/*	class clinstg(ref = '1') chemo(ref = '0');*/
/*	model dftime*cens(0) = age hgb clinstg chemo;*/
/*run;*/

/*taking a as an indicator of age > 65*/
/*here we treat all event types the same*/
title "covariates effect via cox model: treating all event types as the same";
proc phreg data=follicular2;
	class a(ref = '0') clinstg(ref = '1') chemo(ref = '0');
	model dftime*cens(0) = a hgb clinstg chemo;
run;
/*comments to be deleted:*/
/*age group and clinical stage have  highly significant effect on treatment failure.
chemotherapy is marginally significant. Patients who are older than 65 at entry have about twice the risk
of having a treatment failure. Patients who are at clinical stage 2 at entry have about 50% more risk 
of having a treatment failure(no resopnse, relapse or death before relapse). Patients having chemotherapy 
have about 25% reduced risk of having experiencing a treatment failure. The global test for regression coefficients
suggest that there are defintely some effect from these covariates(p value around 0.001).*/

title "consider only event of interest";
proc phreg data=follicular2;
	class a(ref = '0') clinstg(ref = '1') chemo(ref = '0');
	model dftime*cens(0, 2) = a hgb clinstg chemo;
run;
/*or equivalently*/
proc phreg data=follicular2;
	class a(ref = '0') clinstg(ref = '1') chemo(ref = '0');
	model dftime*evcens(0) = a hgb clinstg chemo;
run;

/*comments to be deleted:*/
/*age group and clinical stage have highly sigificant effect on event of interest, i.e.
no response or relapse. chemotherapy is marginally significant. Patients who are older than 65 have about
50% more risk of having no response or relapse, patients who are at clinical stage 2 have about 60% more risk of 
having no response or relapse. Patients having chemotherapy have about 30% reduced risk of having no response or 
relapse. The global test for regression coefficients still tell covariates effect are significant (p value <0.001)*/

title "consider the competing risk as event";
proc phreg data=follicular2;
	class a(ref = '0') clinstg(ref = '1') chemo(ref = '0');
	model dftime*cens(0, 1) = a hgb clinstg chemo;
run;
/*or equivalently*/
proc phreg data=follicular2;
	class a(ref = '0') clinstg(ref = '1') chemo(ref = '0');
	model dftime*crcens(0) = a hgb clinstg chemo;
run;
/*comments to be deleted*/
/*only age group has significant effect on death before relapse. It tells that patients above age of 65
have about 6.5 times the risk of death before relapse than patients who are under 65.*/
/*This is not surprising to see. Since we are looking at death before relapse, covariates reflecting
the condition of  cell-lymphoma should not be directly contributing to the death. Age is a significant here because
older peole in general sustain high risk of death due to other causes.*/

/*comments to be deleted*/
/*coefficients can differ greatly across different event types. But perhaps these differences are merely the result of 
random variation.[cite allison here]. What we need is a test of null hypothesis beta_j = beta for all j, where beta_j
is the vector of coefficients for event type j.A test statistic is readly constructed from output given by proc phreg.
For each model, proc phreg reports -2*log-likelihood(for the model with covariates). For the three models we just estimated
the values are:
All types combined: 3888.993
event of interest(no response or relapse): 3140.689
competing risk(death before relapse): 719.480
The sume of event types is 3140.689 + 719.480 = 3860.169
The difference between the sum and the all types combined is 
3888.993 - 3860.169 = 28.824, which is the likelihood ratio chi square statistic for the null hypothesis
degrees of freedom is (4 + 4)(from 2 seperate models) - (4)(from the combined model) = 4*/

title "test for beta_j = beta cross all j";
data test;
	p = 1 - CDF('CHISQUARE', 28.824, 4);
run;

proc print data=test;
run;

/*comments to be deleted:*/
/*p value <0.0001 and we reject the null hypothesis that the coefficients are all equal across event types.*/
/*the results are not surprising because death before relapse is unlikely to have the same determinants as 
no response or relapse.*/
/*three quick comments here:
1: Above test is valid based on the breslow method of handling ties. If other methods are used, for example, efron
we need to use the combine dataset created earlier, and stratify on the event type. detail see allison
2. If somehow the above test ended up fail to reject, the the next test we would be interested in is to see whether 
beta = 0 or not. Then we can go back to look at the global test for beta = 0 when we consider all event types together.
3. we can construct test statistics for hypothesis about coeffieicnets for specific covariates.

For age group:
(0.41908- 1.87145)^2/([0.13567]^2 + [0.26520]^2) = 23.77102
For hgb::
(0.00230 - 0.00162)^2/((0.00406)^2 + (0.00812)^2)= 0.0056
For clinical stage:
(0.48191 - 0.21985)^2/((0.13095)^2 + (0.28336)^2)= 0.7048
For chemotherapy:
(-0.32457+ 0.09558)^2/((0.16639)^2 + (0.35714)^2)= 0.3378
*/
/*quick comment here for the test above, some may be tempted to ask that we may need a covariance term on the denominator
of test statistic, since we do not have independent groups here (we are using all 541 patients in both models).
However it is not an issue here, since the likelihood function fctors into distinct likelihood for each event
type, and the parameter estimates for each event type are asymptotically independent of the parameter estimates for all 
other event types.(cite Allison)*/
title "test for individual regression coefficient between models";
data test2;
	p1 = 1 - CDF('CHISQUARE', 23.7102, 1);
	p2 = 1 - CDF('CHISQUARE', 0.0056, 1);
	p3 = 1 - CDF('CHISQUARE', 0.7048, 1);
	p4 = 1 - CDF('CHISQUARE', 0.3378, 1);
run;

proc print data=test2;
run;
/*comments to delete*/
/*we have significant statistical justification for conclude that the coefficients for age group are different
between the two models*/
/*for other three covariates, we do not have enough evidence to support the conclusion that the coefficients are
different.(p values: <0.0001, 0.94, 0.40, 0.56)*/
/*we may want to further test whether the three covariates are equal to 0, given that we found no reason to reject that
they are different in the two models*/
/*For hgb:
0.3196 + 0.0399 = 0.3595
For clinical stage:
13.5429 + 0.6020 = 14.1449
For chemotherapy:
3.8049 + 0.0716 = 3.8765*/
title "test for individual regression coefficients to see if they are 0";
data test3;
	p1 = 1 - CDF('CHISQUARE', 0.3595, 2);
	p2 = 1 - CDF('CHISQUARE', 14.1449, 2);
	p3 = 1 - CDF('CHISQUARE', 3.8765, 2);
run;
proc print data=test3;
run;
/*the p values are (0.84, <0.001, 0.14)*/
/*so we fail to reject that the coefficients for hgb and chemotherapy are 0, and reject the null hypothesis that
the coefficient for clinical stage is 0.*/



/*Accelerated Failure Time Models*/

/*exclude covariates that are not significant in any of those previous cox model*/

/*exponential model*/
title "exponential model for event of interest";
proc lifereg data=follicular2;
	model dftime*evcens(0) = a clinstg/dist = exponential;
run;
title "exponential model for competing risk";
proc lifereg data=follicular2;
	model dftime*crcens(0) = a clinstg/dist = exponential;
run;

/*weibull model*/
title "weibull model for event of interest";
proc lifereg data=follicular2;
	model dftime*evcens(0) = a clinstg/dist = weibull;
run;
title "weibull model for competing risk";
proc lifereg data=follicular2;
	model dftime*crcens(0) = a clinstg/dist = weibull;
run;

/*generalized gamma model*/
title "generalized gamma model for event of interest";
proc lifereg data=follicular2;
	model dftime*evcens(0) = a clinstg/dist = gamma;
run;
title "generalized gamma model for competing risk";
proc lifereg data=follicular2;
	model dftime*crcens(0) = a clinstg/dist = gamma;
run;

/*log normal model*/
title "log normal model for event of interest";
proc lifereg data=follicular2;
	model dftime*evcens(0) = a clinstg/dist = lnormal;
run;
title "log normal model for competing risk";
proc lifereg data=follicular2;
	model dftime*crcens(0) = a clinstg/dist = lnormal;
run;

/*log logistic model*/
title "log logistic model for event of interest";
proc lifereg data=follicular2;
	model dftime*evcens(0) = a clinstg/dist = llogistic;
run;
title "log logistic model for competing risk";
proc lifereg data=follicular2;
	model dftime*crcens(0) = a clinstg/dist = llogistic;
run;

/*comment to be deleted:*/
/*a table for -2log-likelihood:
				event of interest		competing risk
exponential        1825.207          	  436.588						
weibull 		   1673.698               415.013
gamma			   1673.695               413.606
lognormal		   1700.212               425.910
log logistic       1673.949               419.122
*/
/*a table for AIC
				event of interest        competing risk
exponential		   1831.207			      442.588
weibull			   1681.698				  423.013
gamma			   1683.695				  423.606
lognormal		   1708.212				  433.910
log logistic	   1681.949				  427.122
*/	
/*a table for BIC
				 event of iterest		 competing risk
exponential		   1844.087			      455.468
weibull			   1698.872				  440.187
gamma			   1705.162				  445.073
lognormal		   1725.386				  451.083
log logistic	   1699.122				  444.296
*/

/*comments to be deleted
For events of interest(relapse or no response to treatment)
exponential model should be rejected. The test on scale = 1 has a highly significant p value(<0.0000001)
, also the confidence interval for the scale parameter from the output of weibull model is (1.6147, 1.9965) which does
not include 1.
The goodness of fit test between weibull and generalized gamma has a chi-square test statistic with degree of freedom 1:
(1673.698 - 1673.695) = 0.003, which gives a p value of 0.96 
and it highly insignificant. Among all 5 models the weibull model has the smallest AIC
and BIC value as well. So we will be in favor of using Weibull model here.*/
/*comments to be deleted
For competing risk(death before relapse)
exponential model should also be reject. The test on sclae = 1 has a highly significant p value(<0.0000001),
also the confidence interval for the scale parameter form the output of weibull model is (0.5451, 0.7621), which does not
include 1. 
The goodness of fit test between weibull and generalized gamma has a chi-square test statistic with degree of freedom 1:
(415.013 - 413.606) = 1.407, which gives a p value of 0.24 and is not significant. So between weibull and generalized gamma
we are in favor of Weibull. Also weibull model has the smallest AIC and BIC among all five models, so clearly we will be
in favor of weibul model here.*/

/*For the event of interest scale parameter estimate is 1.796, using the transformation 1/(1.796)- 1 = -0.44, we got 
the coefficient of log t in the equivalent proportional hazards model. This result indicate that the hazard of treatment 
failure(relapse or no response) decreases with time . On the othe hand, for death before relapse, we got
1/0.6445 - 1 = 0.55 and it indicates that the hazard of death without relapse increases with time.*/

/*the next step is to construct tests of hypothesis about equality of coefficients across different type of events
we could have used the same approach as in cox model to test on the individual coefficients. But we take a different 
approach here instead. For theoretical model, see 225 Allison*/

title "test of difference on regression coefficients between models";
data follicular3;
	set follicular2;
	ldftime = log(dftime);
run;
proc logistic data=follicular3;
	where cens = 1 OR cens = 2;
	class a(ref = '0')/param = glm;
	/*specify the death as reference category, so the event of interest is in the event category*/
	model cens(ref = '2') = ldftime a clinstg;
run;
/*comments to be deleted*/
/*we can interpret each of the coefficients in the output as an estimate of the difference between corresponding
coefficients in the two Weibull models. There are highly sigificant differences in the coefficients for age group, 
log of time from diagnosis to first failure(no response, relapse or death). The test for ldftime is equivalent to
a test of whether the scale parameters are the same in the accelerated failure time version of the model
The null hypothesis that all the corresponding coefficients are equal in the two Weibull models is equivalent to the hypothesis that all the 
coefficients in the implied logit model are 0. That hypothesis is rejected if we look at the output.*/

/*Other Approach*/

/*the issue of non-independence between different types of events*/

/*1. Conditional Process Approach*/
/*P(T = t, J = j) = P(T = t)P(J = j|T = t)*/

/*Step 1. For the second term (the conditional probability), we could use the logistic model from above*/
/*Given that event happens (either no response or relapse, or death before response), the patients above age 65
only has about 25% odds of experiencing event of interest, */

/*Step 2. We still need a model for timing of event*/
/*here we also removed chemo and hgb which are not significant in some of those previous models*/
proc phreg data=follicular2;
	class a(ref = '0') clinstg(ref = '1');
	model dftime*cens(0) = a clinstg;
run;
/*comments to be deleted:*/
/*from the output, we can interpret that patients who are above age 65 have about twice as much risk of 
experiencing treatment failure(by means of either no response or relapse, or death before relapse).
also, patients who are at clinical stage 2 have about 50% more risk of experiencing treatment failure.*/

/*2. Cumulative Incidence Functions*/
%CUMINCID(data = follicular2, time = dftime, status = cens, event = 1, compete = 2, censored = 0)

