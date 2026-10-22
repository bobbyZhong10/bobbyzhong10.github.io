---
title: "Maximum Likelihood"
subtitle: "Likelihood, Information, and Misspecification"
seriesline: "Foundations of Information Systems Economics · Chapter 6"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 6 · Maximum Likelihood"
---

## Introduction

Sorrel is a fictional meal-kit subscription platform: users choose a monthly fee tier (USD 45, 60, 75, or 90) and receive ingredients and recipes every month until they cancel. The growth team wants to answer a pricing question: for every USD 1 added to the monthly fee, how much faster do users cancel? The data at hand are subscription records exported 36 months after the platform launched. On the day of the export, 43.3% of these users were still subscribed, and all we know about them is that they have subscribed "at least this long." Section 1 shows that if these unfinished subscriptions are treated as finished and tenure is then regressed on the fee, the fee effect comes out nearly a third too small.

This chapter is about maximum likelihood. Chapter 5 already placed it in the framework of M-estimation: take $q(w,\theta) = -\ln f(y \mid x;\theta)$, and consistency and asymptotic normality are special cases of Theorems 2 and 3 of Chapter 5, which this chapter does not repeat. Instead, it covers four things that are unique to the likelihood. First, the likelihood knows how the data were observed: observations that are only partly seen, such as censored or grouped ones, can enter the likelihood as they are, provided we can write down the probability of observing them; there is no need to discard them or to pretend they were seen in full. Second, under correct specification the information matrix equality collapses the sandwich of Chapter 5 into the inverse of the information matrix, MLE has the smallest variance in a large class of estimators, and the three tests, Wald, LR, and LM, are asymptotically equivalent. Third, under misspecification the likelihood estimates the model closest to the true distribution in the sense of KL divergence, and then only the sandwich can be trusted. Fourth, the likelihood can be written down only in part (partial MLE) or estimated in two steps (two-step MLE).

The chapter proceeds as follows. Section 2 sets out the hazard model of tenure and the estimands of the four data blocks. Section 3, the most detailed in the chapter, takes up in turn the consistency of conditional MLE, the likelihood under censoring, the information matrix equality and asymptotic normality, efficiency, the trinity of tests, QMLE and the Vuong test, partial MLE and grouped durations, and two-step MLE. Section 4 covers algorithms, variance estimation, and software, Section 5 examines three anchoring papers, Section 6 runs the whole toolkit on Sorrel's data, and Section 7 discusses failure modes. The theory of conditional MLE, partial MLE, and the trinity established in this chapter returns in Chapters 8 and 12.

## 1 Censoring Shrinks the Fee Effect

::: {.case}
Sorrel's block C1 has 5,000 subscribers who signed up over 36 months, with the number of signups growing 4% a month. The four fee tiers hold 30%, 35%, 20%, and 15% of subscribers, and about three in ten users joined through a friend's referral code. The data were exported at the end of month 36: a user who signed up in month $\mathrm{signup}_i$ and was still subscribed at the export is observed for only $37 - \mathrm{signup}_i$ months. Of the 5,000, 2,834 are observed to cancel, and the remaining 43.3% are still subscribed at the export.

Because this is a simulation, we know the answer. Raising the monthly fee by USD 1 multiplies the hazard of cancellation by 1.0202, a rise of 2.02%; moving up one tier (USD 15) raises it by 35.0%; and joining through a referral makes the hazard lower by 33.0%. The hazard also declines slowly with tenure. For non-referred users in the USD 60 tier, median tenure is 8.02 months.
:::

The analyst's first instinct is a regression: take logs of each user's observed duration and run OLS on the fee and a referral indicator, counting users still subscribed at the export as having subscribed for $37 - \mathrm{signup}_i$ months. The fee coefficient is -0.01512 (robust standard error 0.00119), read as "for every USD 1 added to the monthly fee, tenure falls about 1.5%." Section 2.1 shows that in Sorrel's model the true effect of the fee on log duration is $-0.02/0.9 = -0.02222$. The naive regression recovers only 0.680 of the truth, and its 95% interval [-0.01744, -0.01279] lies far from it. Running the same regression on only the 2,834 users who have already canceled does even worse. This is not sampling error: run the same regression on the latent durations, which only a simulation can see, and the coefficient is -0.02319, with an interval that covers the truth.

![Estimates and 95% intervals of the fee coefficient $\beta_{\mathrm{fee}}$ (log hazard per dollar); the vertical dashed line is the truth, 0.02. Triangles are the two OLS regressions that ignore censoring (treating censored observations as exits, or using only users observed to cancel), converted to the hazard scale through the Weibull log-linear form of Section 2.1, $\hat\kappa = \pi/\sqrt{6s^2}$ and $\hat\beta_{\mathrm{fee}} = -\hat\kappa\,\hat b_{\mathrm{fee}}$. Circles are the two estimates that handle censoring through the likelihood (the continuous-time Weibull MLE on block C1, and the binary regression on the first 24 months of block C3, grouped by month, which is the grouped-duration MLE of Section 3.8). The square is the infeasible OLS on the latent durations $t^*$, converted by the same formulas.](assets/fig/fig_06_censoring.svg)

The figure puts the same estimates on the hazard scale: the naive regression gives 0.01533 against a truth of 0.02. Change the estimator and the problem disappears. In the Weibull MLE on the same 5,000 users, censored users are not treated as exits; in the likelihood they contribute only the probability of having subscribed "at least this long." The hazard coefficient of the fee is 0.02056 (standard error 0.00120), and the interval [0.01820, 0.02292] covers the truth, 0.02. Converted to log duration it is -0.02323, with an interval that likewise covers -0.02222. Moving up one fee tier multiplies the hazard by 1.361, against a true 1.350. This chapter answers three questions: why the naive regression is wrong, and which assumption lets the likelihood recover the information that was cut off; when the likelihood's built-in standard errors can be trusted; and what MLE estimates when the Weibull is misspecified.

## 2 The Economic Model and the Estimand

### 2.1 Tenure, the Hazard, and Censoring

Let $t_i^*$ be the latent tenure of subscriber $i$ in months, that is, how long they would stay subscribed if the data were never exported and we could watch indefinitely, and let $c_i$ be the censoring time, the number of months from signup to the data export. What we observe is $t_i = \min(t_i^*, c_i)$ and $\delta_i = \mathbf{1}\{t_i^* \le c_i\}$. Here $\delta_i = 1$ means that the cancellation was seen and $\delta_i = 0$ that the user was still subscribed at the export, a situation called right censoring.

The central object of a duration model is the hazard: the rate at which someone who has stayed subscribed until month $t$ cancels in the next short stretch of time, $h(t \mid x) = f(t \mid x)/S(t \mid x)$, where $f$ is the conditional density of $t^*$ and $S(t \mid x) = \Pr(t^* > t \mid x)$ is the survival function. Here $h$ is the hazard; it has nothing to do with the variance function $\mathrm{Var}(y \mid x) = \sigma^2 h(x)$ of Section 3.4 of Chapter 5 or the $h(x)$ in it. The hazard and the survival function determine each other: $S(t \mid x) = \exp\{-\int_0^t h(\tilde t \mid x)\,d\tilde t\}$ (Wooldridge 2010 equation (22.7)). Why not model the mean duration directly? Because the question the platform cares about is conditional to begin with: "for a user who has already stayed half a year, how likely is it that they leave next month?" Moreover, as Section 3.2 shows, the hazard makes censoring easy to handle.

The model for C1 is a Weibull proportional hazard model, with $x_i = (1, \mathrm{fee}_i, \mathrm{ref}_i)'$ and $\beta = (\beta_0, \beta_{\mathrm{fee}}, \beta_{\mathrm{ref}})'$:

$$h(t \mid x) = \exp(x'\beta)\,\kappa\,t^{\kappa - 1}, \qquad S(t \mid x) = \exp\{-\exp(x'\beta)\,t^{\kappa}\}.$$

The factor $\exp(x'\beta)$ scales the hazards of different people proportionally, and $\beta_{\mathrm{fee}}$ is the semi-elasticity of the hazard with respect to the fee, with a true value of 0.02, a rise of 2.02% per dollar. The shape parameter is $\kappa$: $\kappa = 1$ gives the exponential distribution, whose hazard does not change with duration; with $\kappa < 1$ the hazard falls with duration, which is called negative duration dependence (Wooldridge 2010 Example 22.3, p. 987). In C1, $\kappa = 0.9$: the longer users stay, the less likely they are to leave, but the decline is slow. Median duration is $(\ln 2/\exp(x'\beta))^{1/\kappa}$, which for non-referred users in the USD 60 tier is 8.02 months, and the probability of staying a full 12 months is 0.369.

The same model can also be written as a regression for duration. Let $\varepsilon_i = \ln\{\exp(x_i'\beta)\,t_i^{*\kappa}\}$; it is the log of an exponential random variable with mean one and is independent of $x_i$ (Wooldridge 2010 p. 998), so

$$\ln t_i^* = -\frac{x_i'\beta}{\kappa} + \frac{\varepsilon_i}{\kappa},$$

where $\varepsilon_i$ has mean equal to minus Euler's constant and variance $\pi^2/6$. This is called the accelerated failure time (AFT) form: the effect of the fee on log duration is $-\beta_{\mathrm{fee}}/\kappa = -0.02222$. Without censoring, the OLS slope estimates $-\beta/\kappa$ and the residual variance estimates $\pi^2/(6\kappa^2)$, so one can back out $\hat\kappa = \pi/\sqrt{6s^2}$ and $\hat\beta_{\mathrm{fee}} = -\hat\kappa\,\hat b_{\mathrm{fee}}$; this is how the two OLS rows in the figure of Section 1 were converted.

The chapter's estimand on C1 is the parameter of the conditional density, $\theta = (\beta', \ln\kappa)'$, and the quantity of interest is $\beta_{\mathrm{fee}}$ or its AFT counterpart; the shape parameter enters in logs so that the parameter space has no boundary. Whether $\beta_{\mathrm{fee}}$ can be read as "the causal effect of a price increase" depends on how users choose their tiers. In the DGP the tier is independent of the other determinants of duration, so it can. In real data, people who choose an expensive tier may be different to begin with; that is an identification problem the likelihood cannot solve, and this chapter deals only with estimation and inference given the covariates. The censoring time is set by the signup month: a user who signs up in month $\mathrm{signup}_i$ has $c_i = 37 - \mathrm{signup}_i$. In the DGP the signup month and the latent duration are independent, and the population share of censored observations is 0.420; this independence is the key assumption of Section 3.2.

### 2.2 The Other Three Data Blocks

Block C2 holds 8,000 Sorrel subscribers in another city. Fees, referrals, and the export are the same as in C1, but durations come from a mixture of two kinds of people. One kind is deal-seekers, users who come for the first-box discount; they make up 37.8% of non-referred users and 18.2% of referred users. Their durations follow a log-logistic distribution with median 1.5 months and shape parameter 0.5: many leave as soon as the first box arrives, yet a few stay for years. The other kind is regular users, with median 12 months and shape parameter 5, whose hazard first rises and then falls, peaking in month 15.8 for the USD 60 tier. The fee stretches or shrinks the time scales of both groups together: each additional USD 1 makes the time scale smaller by 1.98%. The analyst still fits the Weibull used for C1. No Weibull is true. For non-referred users in the USD 60 tier, the true hazard falls to its minimum of 0.0187 in month 4, rises to 0.175 in month 13.5, and then declines again; a monotone Weibull hazard cannot draw this curve. The estimand is now the pseudo-true value $\theta_*$, the minimizer of the population expected negative log-likelihood. The pseudo-true shape parameter is $\kappa_* = 0.7586$, obtained by fitting the Weibull to 2,000,000 simulated observations; numerical integration of the population objective gives 0.7588, a difference of about $1.8 \times 10^{-4}$.

Block C3 splits C1's users into months, with one row per user for each month subscribed, a layout called person-period data. Only the first 24 months are used, giving 37,819 rows and 2,784 cancellations; users who cancel after month 24 are recorded as censored at the end of month 24. Because C1's hazard is proportional, one free baseline term per month, $\eta_m$, is just enough, with $\eta_m = \beta_0 + \ln\{m^\kappa - (m-1)^\kappa\}$, and the coefficients of the fee and the referral are still C1's $\beta$.

Block C4 holds 5,000 subscribers who signed up after the app launched. At signup the platform records household size and a standardized score on a short cooking quiz, written $z_i = (1, \mathrm{hh}_i, \mathrm{quiz}_i)'$. What really drives whether users stay is a usage-intensity index $z_i'\gamma$, with $\gamma = (-1.0, 0.25, 0.35)'$, the log of the expected number of recipe check-ins in the app during the first week. The platform sees only one week's check-ins $v_i$, which given $z_i$ follow a Poisson distribution with mean $\exp(z_i'\gamma)$; the average is only 0.816, and 47.1% of users have none at all. The hazard is $\exp(\beta_0 + \beta_{\mathrm{fee}}\mathrm{fee}_i + \beta_{\mathrm{use}}\,z_i'\gamma)\,\kappa\,t^{\kappa - 1}$ with $\beta_{\mathrm{use}} = -0.80$: an index higher by 1 makes the hazard lower by 55.1%. The index must be estimated first, which makes this a two-step MLE.

The estimands of the four data blocks share one form: each is the minimizer of the population expected negative log-likelihood $\mathbb{E}[-\ln f(y \mid x;\theta)]$. When the model is correctly specified, this point is the true parameter; when it is misspecified, it is the member of the model family closest to the true distribution, with closeness measured by KL divergence. Section 3 begins by turning this sentence into theorems.

## 3 Identification and Large-Sample Theory

This section follows the format of Chapter 5: it first uses Sorrel's numbers or figures to show what is to be proved, then states numbered assumptions and theorems, followed by proof sketches of two to six steps; each step names the assumptions it uses, and regularity conditions are only cited to their source. After each theorem, an "In words" paragraph answers four questions: under which assumptions it estimates which quantity, what it estimates when the specification is wrong, who sets the weights, and what Sorrel's numbers say. Assumptions and theorems are numbered within this chapter; the parentheses give the theorem numbers, equation numbers, and page numbers of Wooldridge (2010). Sections 3.1 to 3.4 assume a correctly specified likelihood, Section 3.5 covers testing, Sections 3.6 and 3.7 cover misspecification, and Sections 3.8 and 3.9 cover likelihoods written down only in part and two-step estimation.

### 3.1 Conditional MLE: The Likelihood Is a Loss Function Too

Start with what a single subscriber in C1 contributes to the likelihood. A user whose cancellation is seen ($\delta_i = 1$) contributes the density of "canceling exactly at $t_i$," $f(t_i \mid x_i;\theta)$; a user still subscribed at the export contributes the probability of "having subscribed for at least $c_i$ months," $S(t_i \mid x_i;\theta)$. Rewriting with $f = hS$, the log-likelihood of a single observation is

$$\ell_i(\theta) = \delta_i \ln h(t_i \mid x_i;\theta) - \int_0^{t_i} h(\tilde t \mid x_i;\theta)\,d\tilde t,$$

which under the Weibull is $\ell_i(\theta) = \delta_i[x_i'\beta + \ln\kappa + (\kappa - 1)\ln t_i] - \exp(x_i'\beta)\,t_i^{\kappa}$. The expression reads naturally: an observed cancellation adds the log of the hazard at the moment of cancellation, and everyone subtracts the hazard accumulated over their observation window.

In general, write $w = (y, x)$, let $f(y \mid x;\theta)$ be a model for the conditional density of $y$ given $x$, let $f_0(y \mid x)$ be the true conditional density, and let $\ell_i(\theta) = \ln f(y_i \mid x_i;\theta)$ and $\ell_n(\theta) = \sum_{i=1}^n \ell_i(\theta)$. Following the convention of Chapter 5, take $q(w,\theta) = -\ln f(y \mid x;\theta)$, so that $Q_n(\theta) = -\ell_n(\theta)/n$ and MLE is the minimizer of $Q_n$. The estimator defined this way is called CMLE (conditional MLE; Wooldridge 2010 §13.3). The word "conditional" means that only $y$ given $x$ is modeled, and no parametric form is imposed on the distribution of $x$. Sorrel need not model how users choose their fee tiers. Leaving the regressors unmodeled in this way is the most common form of CMLE in applications, and the binary response models of Chapter 12 belong to the same class.

Why is $\theta_0$ the minimizer of the population objective? Section 3.2 of Chapter 5 proved the unconditional version with Jensen's inequality; the conditional version simply repeats the argument for each $x$:

$$\mathbb{E}[q(w,\theta) \mid x] - \mathbb{E}[q(w,\theta_0) \mid x] = \mathrm{KL}\big(f_0(\cdot \mid x) \,\Vert\, f(\cdot \mid x;\theta)\big) \ge 0,$$

provided $f(\cdot \mid x;\theta_0) = f_0(\cdot \mid x)$. The right side is the conditional KL divergence (Wooldridge 2010 equation (13.9)), which is zero only when the two conditional densities are the same.

::: {.assumption}
**ML.1 (correct specification and identification; Wooldridge 2010 Theorem 13.1 (a) and (b), equations (13.8), (13.11), and (13.14), pp. 473-475)** For all $x$ and $\theta \in \Theta$, $f(\cdot \mid x;\theta)$ is a density; there exists $\theta_0 \in \Theta$ such that $f(\cdot \mid x;\theta_0) = f_0(\cdot \mid x)$ for all $x$; and $\theta_0$ is the unique minimizer of $\mathbb{E}[q(w,\theta)]$ on $\Theta$.

**ML.2 (Wooldridge 2010 Theorem 13.1 (c) to (f), p. 475)** $\{(x_i, y_i)\}$ is a random sample, $\Theta$ is compact, $\ln f(y \mid x;\theta)$ is Borel measurable for each $\theta$ and, for each $w$, continuous on $\Theta$, and $\lvert \ln f(y \mid x;\theta) \rvert \le b(w)$ for all $\theta$, with $\mathbb{E}[b(w)] < \infty$. That is, M.2 of Chapter 5 holds for $q = -\ln f$.
:::

ML.1 is the only assumption that matters, and it has two halves. The first half is correct specification: the true density is in the model family, which in C1 means "duration follows a Weibull with these three covariates." The second half is identification: different values of $\theta$ give different densities, which in C1 holds as long as the fee and the referral indicator are not perfectly collinear and some users are observed to cancel. ML.2 is a technical condition; Wooldridge points out that in practice almost the only things to check are correct specification, identification, and continuity of the log-likelihood (p. 476).

::: {.theorem}
**Theorem 1 (consistency of the CMLE; Wooldridge 2010 Theorem 13.1, p. 475)** Under ML.1 and ML.2, $Q_n(\theta) = -\ell_n(\theta)/n$ has a minimizer $\hat\theta$, and $\hat\theta \xrightarrow{p} \theta_0$.
:::

The proof sketch has three steps. Step one (ML.1): the conditional KL inequality above shows that $\theta_0$ minimizes the conditional expected loss at every $x$ (equations (13.9) to (13.12)). Step two (ML.1): taking expectations over $x$, $\theta_0$ is a minimizer of $\mathbb{E}[q(w,\theta)]$ (equation (13.14)), and uniqueness is M.1 of Chapter 5. Step three (ML.2): Theorem 2 of Chapter 5 (Wooldridge 2010 Theorem 12.2) gives the result.

In words: under ML.1 and ML.2, CMLE estimates the true parameter of the conditional density. C1's fee coefficient, 0.02056, aims at 0.02, and the shape parameter, 0.8848 (standard error 0.0136), aims at 0.9. Under misspecification, the KL inequality of step one still holds, but $f_0$ is not in the model family, so the minimizer becomes the pseudo-true value of Section 3.6. Who sets the weights: $\mathbb{E}[q]$ differs from the average of the conditional KL divergence, weighted by the distribution of $x$, only by a constant that does not depend on $\theta$. Under correct specification the KL divergence at every $x$ is zero at $\theta_0$, so the distribution of $x$ affects only the strength of identification, not where $\theta_0$ lies; this is why CMLE can ignore how $x$ is generated. Under misspecification, the distribution of $x$ helps determine the pseudo-true value.

### 3.2 Censoring: Which Assumption Lets the Likelihood Recover What Was Cut Off

Why is the naive regression of Section 1 wrong? Users on low fees would stay long anyway, so they are the most likely to run into the export date and be cut short; users on high fees leave early and are rarely cut off. Cutting durations short compresses the gap between tiers and pulls the slope toward zero, leaving only 0.680 of the truth. It also compresses the dispersion of durations, so the residual variance shrinks and the shape parameter converted as in Section 2.1 is 1.014, larger than the true 0.9. On the hazard scale the two biases partly cancel, giving 0.01533, or 0.767 of the truth, but the offset is pure coincidence. Using only the users observed to cancel is worse: it systematically keeps only short durations, especially among low-fee users, and the coefficient on the hazard scale is only 0.0107, 0.535 of the truth.

::: {.intuition}
What a user still subscribed at the export tells you is not "they subscribed for 14 months" but "they subscribed for at least 14 months." Treating them as canceling after 14 months invents an ending for them; dropping them pretends they do not exist, and it is the most loyal users who get dropped. The likelihood keeps honest books: this user's contribution to the likelihood is the probability of "surviving past 14 months," $S(14 \mid x)$. Keeping these books straight requires only one thing: the fact that the user was cut off at month 14 must have nothing to do with how long they would otherwise have stayed.
:::

::: {.assumption}
**CENS.1 (independent censoring; Wooldridge 2010 equation (22.22), p. 994)** Given $x_i$, the conditional distribution of the latent duration $t_i^*$ does not depend on the censoring time $c_i$ or on the time of signup.
:::

CENS.1 allows the censoring time to differ across people and to be correlated with the covariates (for example, if referred users tend to sign up later); it rules out only one thing: "given the covariates, the censoring time still carries information about the latent duration." Sorrel's censoring comes from a single export date, so the censoring time is determined by the signup month alone, and CENS.1 holds provided people who signed up in different months do not differ systematically beyond the covariates. If early users are inherently more loyal, the signup month must be put into $x$ (Wooldridge 2010 p. 994). This assumption cannot be tested in the data: the part of the duration that was cut off is the part we cannot see.

::: {.theorem}
**Theorem 2 (the likelihood under independent right censoring; Wooldridge 2010 §22.3.2, equations (22.21) to (22.24), pp. 993-995)** Let $f(t \mid x;\theta)$ be a model for the conditional density of $t^*$ given $x$, correctly specified at $\theta_0$, and let CENS.1 hold. Write $\ell_i(\theta) = \delta_i \ln f(t_i \mid x_i;\theta) + (1 - \delta_i)\ln S(t_i \mid x_i;\theta)$, and suppose that $\theta_0$ is the unique maximizer of the expected log-likelihood of the censored data, $\mathbb{E}[\ell_i(\theta)]$, on $\Theta$ and that ML.2 holds for this $\ell_i$. Then the conditional density of $(t_i, \delta_i)$ given $(x_i, c_i)$ is

$$f(t_i \mid x_i;\theta_0)^{\delta_i}\,S(t_i \mid x_i;\theta_0)^{1 - \delta_i},$$

and the CMLE with log-likelihood $\ell_i(\theta)$ satisfies Theorem 1 and is asymptotically normal under ML.3 of Section 3.3.
:::

The proof sketch has three steps. Step one (CENS.1): if $\delta_i = 1$, what is observed is $t_i = t_i^* \le c_i$, and CENS.1 makes the density of $t^*$ given $(x, c)$ still $f(t \mid x;\theta_0)$. Step two (CENS.1): if $\delta_i = 0$, all that is observed is the event $t_i^* > c_i$, whose probability is $S(c_i \mid x_i;\theta_0)$, and in this case $t_i = c_i$. Step three (identification and ML.2): together the two cases give equation (22.23), the true conditional density of $(t_i, \delta_i)$ given $(x_i, c_i)$; $c_i$ and the time of signup enter only through $t_i$ and need not be modeled separately. Substituting $w_i = (t_i, \delta_i, x_i, c_i)$ into Theorem 1 gives the result (Wooldridge 2010 p. 995). Identification must hold for the likelihood of the censored data: if some parameter affected only the part of the distribution beyond every censoring time, the data could not pin it down; the Weibull has no such problem as long as someone cancels within the observation window.

In words: given CENS.1 and a correctly specified Weibull, the MLE on censored data still estimates the parameters of the conditional distribution of the latent duration $t^*$, the same $\beta$ that would be estimated without censoring; on C1, 0.02056 aims at 0.02. What lets the likelihood recover the information that was cut off is CENS.1: it allows the probability of "surviving past $c_i$" to be written with the same $S$. Misspecification comes in two kinds. If the duration distribution is wrong, MLE estimates the pseudo-true value under this censoring design (Section 3.6). If censoring is not independent, for example if accounts whose payment failed in the month before the export are purged from the data early, equation (22.23) is no longer the density of the data, the quantity the estimator aims at has no clear meaning, and the data themselves cannot diagnose the problem. Who sets the weights: a censored user contributes only the single term $-\int_0^{c_i} h(\tilde t \mid x_i;\theta)\,d\tilde t$, and the shorter the observation window, the less information they carry. Users who sign up late have short censoring times, so the platform's growth path determines how much each user weighs in the likelihood. This affects precision, not bias: 43.3% of C1's observations are censored, yet MLE is still consistent, and in line with this, the interval from this single C1 sample covers the truth.

The same logic also appears outside durations. If an outcome has many zeros, monthly add-on purchases for example, the Tobit model of Tobin (1958) writes $y_i = \max(0, y_i^*)$, $y_i^* = x_i'\beta + \varepsilon_i$, $\varepsilon_i \mid x_i \sim \mathcal{N}(0, \sigma^2)$: an observation at zero contributes $\ln\Phi(-x_i'\beta/\sigma)$, and a positive observation contributes $\ln\phi\big((y_i - x_i'\beta)/\sigma\big) - \ln\sigma$ (Wooldridge 2010 equations (17.19) and (17.20), p. 676). The structure is the same as in equation (22.24): fully seen observations contribute a density, and partly seen ones contribute the probability of falling into the region where they lie. The difference is that Sorrel's censoring points vary across people and are set by an external event, the export date, so CENS.1 is an assumption about sampling; in corner-solution applications such as add-on spending, zero is a genuine choice by the consumer, and whether the model can be trusted rests entirely on two distributional assumptions, normality and homoskedasticity. The robustness of Section 3.7 of Chapter 5, which relies only on a correctly specified mean, is not available here.

### 3.3 The Score, the Information Matrix Equality, and Asymptotic Normality

Start with two standard errors for the same coefficient on C1. Computed from the curvature of the likelihood (the Hessian), the standard error of the fee coefficient is 0.001202; computed from the outer product of the scores (OPG), it is 0.001223, a ratio of 0.983. Chapter 5 noted that these two quantities are, in general, the two slices of bread and the filling of the sandwich, with no reason to be equal; here they are almost the same. The reason is that when the likelihood is correctly specified, they are one and the same quantity.

The score follows the definition of Chapter 5, $s_i(\theta) = \nabla_\theta q(w_i,\theta)' = -\nabla_\theta \ell_i(\theta)'$, which differs from the score of the log-likelihood in Wooldridge's book only by a minus sign; the sign affects neither the outer product $s_is_i'$ nor $\mathbb{E}[s] = 0$. In C1 the score for $\beta$ is

$$s_{\beta,i}(\theta) = -x_i\Big(\delta_i - \int_0^{t_i} h(\tilde t \mid x_i;\theta)\,d\tilde t\Big),$$

where the term in parentheses is "the number of cancellations observed minus the cumulative hazard." At the truth its expectation given $(x_i, c_i)$ is zero: the expected number of times a person cancels within the observation window equals the expected cumulative hazard at the end of their observation.

This zero mean is no quirk of the Weibull. For any correctly specified likelihood, differentiating the identity "the density integrates to one" with respect to the parameter gives the zero mean of the score (Wooldridge 2010 equations (13.20) to (13.22)):

$$0 = \nabla_\theta \int f(y \mid x;\theta)\,dy = \int \nabla_\theta \ln f(y \mid x;\theta)\,f(y \mid x;\theta)\,dy = -\,\mathbb{E}_\theta\big[s(w,\theta)' \mid x\big],$$

where $\mathbb{E}_\theta$ denotes an expectation taken under the model density $f(\cdot \mid x;\theta)$. Differentiating once more, with the product rule, gives $\mathbb{E}_\theta[\nabla^2_\theta q \mid x] = \mathbb{E}_\theta[s\,s' \mid x]$. These two equalities hold at every $\theta$, but only at $\theta_0$ is the model density the true density of the data, and only there is an expectation taken under the model an expectation over the data.

::: {.intuition}
Curvature and the variance of the score look completely unrelated, so why should they be equal? They are really two windows onto the same thing, both measuring how firmly the data pin down $\theta$. One window is the curvature: the more sharply the negative log-likelihood bends at the valley floor, the faster the average loss rises as you move slightly away from the truth, and the stronger the data's opinion about $\theta$. The other window is the fluctuation of the score: each observation "votes" on which direction $\theta$ should move, and the more decisive the votes, the more the samples disagree with one another. When the likelihood is correctly specified, each observation's vote is weighted by its true distribution, and the two windows see the same number, the Fisher information. When the model is wrong, the votes are weighted by the wrong distribution, and each window sees a number of its own.
:::

::: {.assumption}
**ML.3 (conditions of Wooldridge 2010 Theorem 13.2, p. 479)** In addition to ML.1 and ML.2: (a) $\theta_0$ is in the interior of $\Theta$; (b) $\ln f(y \mid x;\cdot)$ is twice continuously differentiable on the interior of $\Theta$; (c) integration and differentiation can be interchanged in equations (13.21) and (13.25); (d) every element of $\nabla^2_\theta \ln f$ is dominated by a function with finite expectation; (e) $A_0 = \mathbb{E}[\nabla^2_\theta q(w,\theta_0)]$ is positive definite.
:::

Condition (c) is the new one. It rules out models in which the support of $y$ changes with the parameter, such as an auction in which the lower bound of the bids depends on the parameter; there MLE may not even converge at rate $\sqrt n$ (Wooldridge 2010 p. 479). The common models for durations, counts, and binary responses all satisfy it.

::: {.theorem}
**Theorem 3 (asymptotic normality of the CMLE and the information matrix equality; Wooldridge 2010 Theorem 13.2, equations (13.20) and (13.26) to (13.29), pp. 477-479)** Under ML.1 through ML.3: (i) $\mathbb{E}[s(w,\theta_0) \mid x] = 0$; (ii) the conditional information matrix equality holds, $\mathbb{E}[\nabla^2_\theta q(w,\theta_0) \mid x] = \mathbb{E}[s(w,\theta_0)\,s(w,\theta_0)' \mid x]$, and taking expectations gives $A_0 = B_0 \equiv I(\theta_0)$; (iii)

$$\sqrt{n}(\hat\theta_{MLE} - \theta_0) \xrightarrow{d} \mathcal{N}\big(0,\ I(\theta_0)^{-1}\big).$$
:::

The proof sketch has three steps. Step one (ML.1, ML.3(b) and (c)): the two differentiations of $\int f(y \mid x;\theta)\,dy = 1$ above (equations (13.21), (13.22), and (13.25)) give, at every $\theta$, $\mathbb{E}_\theta[s(w,\theta) \mid x] = 0$ and $\mathbb{E}_\theta[\nabla^2_\theta q \mid x] = \mathbb{E}_\theta[s\,s' \mid x]$; ML.1 says that the model density at $\theta_0$ is the true density, which yields the conditional forms of (i) and (ii) (equation (13.26)), and taking expectations over $x$ gives equation (13.27). Step two (ML.2 and ML.3): ML.3 together with (i) is M.3 of Chapter 5, and Theorem 3 of Chapter 5 gives $\mathcal{N}(0, A_0^{-1}B_0A_0^{-1})$. Step three ((ii)): substituting $A_0 = B_0 = I(\theta_0)$, the sandwich collapses to $I(\theta_0)^{-1}$ (equations (13.28) and (13.29)).

In words: under ML.1 through ML.3, MLE estimates the true parameter, and its sampling fluctuation is governed by a single matrix, $I(\theta_0)$; the curvature and the variance of the score need not be estimated separately. The collapse rests entirely on one fact, that the model density is the true density: the expectations of step one are taken under the model and equal expectations over the data only under correct specification. Under misspecification, (i) and (ii) can both fail, and we are back to the sandwich of Chapter 5 (Section 3.6 of this chapter). Who sets the weights: $I(\theta_0)$ is the average of the information in the individual observations, weighted by the distribution of $x$ and of the censoring time. In a proportional hazard model, the expected information in the $\beta$ block is $\Pr(\delta_i = 1 \mid x_i, c_i)\,x_ix_i'$, proportional to the probability that the person is observed to cancel. So the cost of censoring is fewer observed cancellations, and with them a smaller $I(\theta_0)$ and larger standard errors, not bias. Person by person, though, someone who subscribed for a long time before being censored may contribute more curvature than someone who canceled very early. On C1 the ratio of the Hessian and OPG standard errors is 0.983 and the sandwich gives 0.001183, so the three differ only by sampling error. The logit of Chapter 12 will show the same agreement on another correctly specified likelihood.

There are at least three ways to estimate $I(\theta_0)$ (Wooldridge 2010 equations (13.30) to (13.32)): the average of the sample Hessians, the average of the outer products of the scores, and the expected Hessian, averaged after first taking the conditional expectation over $y$ given $x$. The third has a closed form in many common models and in some cases performs clearly better than OPG in finite samples. The second needs only first derivatives and is the BHHH previewed in Section 4.1 of Chapter 5; it is justified by the information matrix equality, and its price is that it often performs poorly in finite samples (Wooldridge 2010 p. 480).

### 3.4 Efficiency: The Class of Estimators in Which MLE Is Best

Stay with C1. Run OLS on the complete latent durations, which only the simulation has, convert as in Section 2.1, and the standard error of the fee coefficient is 0.001201; the MLE on the censored data gives 0.001202, almost the same. This does not mean censoring is free: the two estimators do not use the same data, and the theorem below does not compare them. One reason the two are close is that OLS uses only two moments of log duration, the mean and the variance, and not the shape of the distribution of $\varepsilon_i$, so it loses information even on complete data. This raises a question: under correct specification, is there any estimator more precise than MLE?

The classic answer is the Cramér-Rao lower bound: any regular unbiased estimator has a finite-sample variance no smaller than $[n\,I(\theta_0)]^{-1}$ (Hansen 2014 Appendix B.11, Theorem B.11.3), and the asymptotic variance of MLE is exactly $I(\theta_0)^{-1}$. This statement requires unbiasedness, and MLE is generally biased. Wooldridge (2010) §14.4 gives a purely asymptotic version that states clearly "in which class it is best."

::: {.theorem}
**Theorem 4 (efficiency of the CMLE; Wooldridge 2010 Lemma 14.1 and Theorem 14.3, equations (14.48) to (14.53), pp. 538-541)** Let ML.1 through ML.3 hold. Consider all estimators that solve a sample moment equation $n^{-1}\sum_{i=1}^n g(w_i, \hat\theta) = 0$, where the moment function satisfies $\mathbb{E}_\theta[g(w,\theta) \mid x] = 0$ for all $x$ and all $\theta$, together with regularity conditions of the kind in M.3 of Chapter 5. The asymptotic variance of any estimator in this class satisfies

$$\mathrm{Avar}\big(\sqrt{n}(\hat\theta - \theta_0)\big) \geq I(\theta_0)^{-1}.$$

The CMLE belongs to this class (take $g = s$) and attains this bound.
:::

The proof sketch has three steps. Step one (ML.1 and ML.3(c)): differentiating $\mathbb{E}_\theta[g(w,\theta) \mid x] = 0$ with respect to $\theta$, as in step one of Theorem 3, gives the Jacobian of the moment function, $D_0 \equiv \mathbb{E}[\nabla_\theta g(w,\theta_0)] = \mathbb{E}[g\,s']$ (the derivation of equation (13.39)). Step two (Theorem 3 of Chapter 5): the asymptotic variance of this estimator is $D_0^{-1}\mathbb{E}[gg']D_0^{-1\prime}$, and that of MLE is $I(\theta_0)^{-1}$. Step three (step one): the difference between the two equals $D_0^{-1}\,\mathbb{E}\big[(g - D_0\,I(\theta_0)^{-1}s)(g - D_0\,I(\theta_0)^{-1}s)'\big]\,D_0^{-1\prime}$, a variance matrix and hence positive semidefinite (Lemma 14.1 and equation (14.51), with the proportionality constant equal to 1).

In words: under ML.1 through ML.3, MLE is asymptotically the most efficient member of the large class of estimators whose moment functions have conditional mean zero under the model at every parameter value. This class includes the converted OLS on complete durations in C1 (under the Weibull model its two moments have conditional mean zero at every $\theta$). The theorem compares only estimators on the same data; what it says is that if complete durations were available, the MLE on the complete data would be at least as precise as that OLS. Under misspecification the theorem no longer holds, the efficiency of MLE is moot, and what one gets instead is the pseudo-true value of Section 3.6. Who sets the weights: step three shows that any moment function can be decomposed into "a linear combination of the score" plus "noise uncorrelated with the score," and the latter only adds variance; MLE uses the score directly, which amounts to weighting each observation by its information under the true distribution. The class excludes estimators that exploit the distribution of $x$: if the distribution of $x$ also depends on $\theta$, the joint likelihood of $(y, x)$ is more efficient, provided the distribution of $x$ is also correctly specified (Wooldridge 2010 pp. 541-542).

### 3.5 The Trinity of Tests

On C1, test whether Sorrel's durations exhibit duration dependence, that is, $H_0: \kappa = 1$ (the exponential distribution, with a hazard that does not change with duration). The three classic statistics are $\mathrm{W} = 63.73$, $\mathrm{LR} = 67.62$, and $\mathrm{LM} = 63.93$, all far above the 5% critical value of 3.84, and they differ by at most 5.75% of LR. Note that the null hypothesis here is false (the truth is 0.9). The theorem below describes behavior when the null holds or the alternative is very close to it, so the closeness of the three numbers in this sample is merely a fact, not a demonstration of the theorem.

Write the null hypothesis as $H_0: R(\theta_0) = 0$, where $R$ is a $\nu \times 1$ vector of restrictions and $\nabla_\theta R(\theta)$ is its $\nu \times k$ Jacobian; in C1, $R(\theta) = \ln\kappa$ and $\nu = 1$. Let $\hat\theta_r$ be the restricted estimate, and $\ell_u$ and $\ell_r$ the maximized unrestricted and restricted log-likelihoods. The three statistics use, respectively, only the unrestricted estimate, the likelihood values at both estimates, and only the restricted estimate:

$$\begin{aligned}
\mathrm{W} &= R(\hat\theta)'\Big[\nabla_\theta R(\hat\theta)\,\Big(\textstyle\sum_{i=1}^n \nabla^2_\theta q(w_i,\hat\theta)\Big)^{-1}\nabla_\theta R(\hat\theta)'\Big]^{-1}R(\hat\theta), \\
\mathrm{LR} &= 2(\ell_u - \ell_r), \\
\mathrm{LM} &= \Big(\textstyle\sum_{i=1}^n s_i(\hat\theta_r)\Big)'\Big(\textstyle\sum_{i=1}^n \nabla^2_\theta q(w_i,\hat\theta_r)\Big)^{-1}\Big(\textstyle\sum_{i=1}^n s_i(\hat\theta_r)\Big).
\end{aligned}$$

Wald asks how far the unrestricted estimate is from the restriction, LR (likelihood ratio) asks how much the likelihood drops once the restriction is imposed, and LM (Lagrange multiplier, also called the score test) asks how far the score at the restricted estimate is from zero.

::: {.assumption}
**TR.1 (Wooldridge 2010 §12.6.1 and §12.6.2, pp. 420-423)** $R$ is continuously differentiable and $\nabla_\theta R(\theta_0)$ has rank $\nu$; under $H_0$, $\theta_0$ is in the interior of $\Theta$.
:::

The rank condition rules out cases in which some parameter is not identified under the null, such as the mean model of Section 3.2 of Chapter 5, in which $\theta_3 = 0$ leaves no way to determine $\theta_4$; the interior condition rules out nulls on the boundary, such as a test that a variance is zero, where neither Wald nor LR is $\chi^2_\nu$ any longer (Wooldridge 2010 p. 420, p. 481), though LM can still be used in some cases (p. 424).

::: {.theorem}
**Theorem 5 (the trinity of tests; Wooldridge 2010 §12.6, equations (12.63), (12.68), (12.69), (12.76), (12.77), and (12.79), pp. 420-431; §13.6, equations (13.35) and (13.36), pp. 481-482)** (i) Under ML.1 through ML.3, TR.1, and $H_0$, $\mathrm{W}$, $\mathrm{LR}$, and $\mathrm{LM}$ all converge in distribution to $\chi^2_\nu$, and $\mathrm{LR} - \mathrm{LM} \xrightarrow{p} 0$; under local alternatives at a distance of order $1/\sqrt n$ from the null, the three have the same noncentral $\chi^2_\nu$ limit and the same asymptotic local power. (ii) Suppose the likelihood is misspecified, but M.1 through M.3 of Chapter 5 hold at the pseudo-true value $\theta_*$ (QML.1 of Section 3.6), TR.1 holds at $\theta_*$, and the null is written $H_0: R(\theta_*) = 0$. Then in general $A_* \neq B_*$; with the variance in the middle of W replaced by the sandwich and LM replaced by the robust form of equation (12.68), both still converge in distribution to $\chi^2_\nu$, but LR no longer does (Wooldridge 2010 p. 504).
:::

The proof sketch has four steps. Step one (TR.1, Theorem 3 of Chapter 5): the Delta method (a first-order expansion of $R(\hat\theta)$ around $\theta_0$) shows that $\sqrt n R(\hat\theta)$ is asymptotically normal, with variance obtained by wrapping $\nabla R$ around both sides of the asymptotic variance of $\hat\theta$; a quadratic form built from a consistent estimate of this variance gives $\chi^2_\nu$ (equation (12.63)). Under correct specification this asymptotic variance is $I(\theta_0)^{-1}$, which gives the W above. Step two (TR.1): expand the sample score at the restricted estimate around $\theta_0$ (equations (12.66) and (12.67)); the restriction makes $\nabla R\,\sqrt n(\hat\theta_r - \theta_0)$ converge to zero in probability, which gives the robust LM (equation (12.68)), and under $A_0 = B_0$ this simplifies to the form above (equation (12.69)). Step three (ML.1, via Theorem 3 (ii)): take a second-order expansion of $\sum_i q(w_i, \hat\theta_r)$ around $\hat\theta$; the first-order term vanishes by the unrestricted first-order condition (equation (12.76)), and substituting the linear representation of $\sqrt n(\hat\theta_r - \hat\theta)$ shows that LR equals a quadratic form in the score with $A_0^{-1}$ in the middle, plus $o_p(1)$ (equation (12.77)); for this quadratic form to be $\chi^2_\nu$, we need $B_0 = A_0$. Step four (Wooldridge 2010 §12.6.4): repeating the first three steps under local alternatives shows that all three limits are noncentral $\chi^2_\nu$, with the same noncentrality parameter when $A_0 = B_0$ (p. 431).

In words: under ML.1 through ML.3 and TR.1, all three tests concern the same restriction on the true parameter and cannot be ranked asymptotically. Under misspecification, the robust W and the robust LM test a restriction on the pseudo-true value $\theta_*$ and remain valid; LR relies on the information matrix equality, and its limiting distribution then depends on unknown parameters (Wooldridge 2010 p. 504), so the $\chi^2$ table can no longer be used. Who sets the weights: W weights by the curvature at the unrestricted estimate, LM by the curvature at the restricted estimate, and LR lets the likelihood function do its own weighting, which is why LR is invariant to reparameterization while W changes with the way the restriction is written (p. 421). On C1, 63.73, 67.62, and 63.93 are all far above 3.84, and the three reach the same conclusion: Sorrel's durations are not exponential. How to choose among the three in practice, and how Wald's reparameterization trap shows up in data, are left for the demonstrations in Sections 4.1 and 6.6 of Chapter 8.

### 3.6 QMLE: What the Likelihood Estimates When It Is Misspecified

Turn to block C2. The true durations are a mixture of two kinds of people; the analyst still fits a Weibull and gets a shape parameter of 0.7275. The three standard errors differ widely: the model standard error (Hessian) is 0.00934, OPG gives 0.00629, and the sandwich 0.01453; the ratio of the model standard error to OPG is 1.484, and the sandwich is 1.556 times the model standard error. Under correct specification these three numbers should be close, as in C1; that they diverge is itself a sign of misspecification.

![Durations of non-referred subscribers in the USD 60 tier in block C2. Left: the cancellation hazard (log scale); the dark gray solid line is the true mixture, the orange dotted and green dash-dotted lines are the hazards of the two groups, deal-seekers and regular users, each on its own, the blue solid line is the Weibull at the pseudo-true value, and the blue dashed line is the sample QMLE. Right: the share still subscribed, for the true mixture, the Weibull at the pseudo-true value and at the sample QMLE, and the sample's Kaplan-Meier estimate (light gray step line).](assets/fig/fig_06_qmle.svg)

The figure makes it plain: the true hazard falls, rises, and falls again, while the Weibull can draw only a monotone line. It picks a slowly declining curve that overstates the true hazard at the trough in month 4 and understates it at the peak in month 13.5, so the two survival curves cross twice. This line is the pseudo-true value: among all Weibulls, the one closest to the true distribution in KL.

::: {.assumption}
**QML.1 (Wooldridge 2010 §13.11.1, pp. 502-503; White 1982)** $f_0$ is not required to be in the model family; $\theta_*$ is the unique minimizer of $\mathbb{E}[q(w,\theta)]$ on $\Theta$; ML.2 holds; and M.3 of Chapter 5 holds at $\theta_*$ for $q = -\ln f$ (including $\mathbb{E}[s(w,\theta_*)] = 0$ and finite second moments of $s(w,\theta_*)$).
:::

QML.1 replaces the "correct specification" of ML.1 with "there is a unique best approximation," and replaces ML.3 with M.3 of Chapter 5: under misspecification $B_*$ no longer equals $A_*$, and a finite variance of the score must be assumed separately. The interchange of integration and differentiation under the model density in ML.3(c) is no longer needed, since it served only expectations taken under the model.

::: {.theorem}
**Theorem 6 (QMLE; Wooldridge 2010 §13.11.1, equation (13.71), pp. 502-504; White 1982)** Under QML.1: (i) $\hat\theta \xrightarrow{p} \theta_*$, where $\theta_*$ is the minimizer of the expected conditional KL divergence $\mathbb{E}\big[\mathrm{KL}\big(f_0(\cdot \mid x) \,\Vert\, f(\cdot \mid x;\theta)\big)\big]$; (ii) $\sqrt n(\hat\theta - \theta_*) \xrightarrow{d} \mathcal{N}(0, A_*^{-1}B_*A_*^{-1})$, with $A_* = \mathbb{E}[\nabla^2_\theta q(w,\theta_*)]$ and $B_* = \mathbb{E}[s(w,\theta_*)\,s(w,\theta_*)']$, and in general $A_* \neq B_*$; (iii) without further assumptions, the only valid variance estimator is the sandwich with the sample Hessian as its bread (equation (13.71)), LM tests must use the robust form, and LR is not asymptotically $\chi^2$ (p. 504).
:::

The proof sketch has three steps. Step one (QML.1): $\mathbb{E}[q(w,\theta)] = \mathbb{E}[-\ln f_0(y \mid x)] + \mathbb{E}[\mathrm{KL}(f_0(\cdot \mid x) \Vert f(\cdot \mid x;\theta))]$, and the first term does not depend on $\theta$, so the minimizer of the population objective is the minimizer of the expected KL divergence; Theorem 2 of Chapter 5 gives (i). Step two (QML.1): $\theta_*$ is an interior minimizer, so the population first-order condition gives $\mathbb{E}[s(w,\theta_*)] = 0$, and Theorem 3 of Chapter 5, applied at $\theta_*$, gives (ii). Step three (the proof of Theorem 3): step one of Theorem 3 takes expectations under the model density, which equal expectations over the data only when the model density equals the true density. At $\theta_*$ the two differ, the conditional mean $\mathbb{E}[s(w,\theta_*) \mid x]$ is generally not zero, and the information matrix equality has no reason to hold; the sandwich does not collapse, and (iii) follows.

In words: under QML.1, QMLE estimates the best KL approximation to the true distribution within the model family, which is not a structural parameter of any kind. C2's pseudo-true shape, $\kappa_* = 0.7586$, says that "if this city must be summarized by a single Weibull, the best one has a hazard that falls with duration." This falling hazard does not describe regular users, whose hazard rises and then falls; what it mainly summarizes is the steep decline of the population hazard in the first few months. That decline has two sources. The deal-seekers' own hazard falls with duration to begin with (a log-logistic with shape parameter 0.5 has a monotonically decreasing hazard), and this is the main part. Selection then amplifies it: most deal-seekers leave very early, so the users who remain in the first few months are increasingly regular users. The pseudo-true fee coefficient, 0.01832, is likewise only the slope of this best-approximating line; in the true structure the fee makes everyone's time scale smaller by 1.98%, and the two are not the same quantity. Who sets the weights: the KL divergence is weighted by the distribution of $x$ and by the censoring design, and censored users contribute only survival probabilities, so if the same city's data were exported half a year earlier, the pseudo-true value would change; the model family the researcher picks determines the shape of the "approximation." For inference only the sandwich can be trusted, and the direction in which the model standard error errs is not fixed: on C2 the model standard error of the shape parameter is too small (in this sample the sandwich is 1.556 times as large), while that of the fee coefficient is, on the contrary, conservative. The Monte Carlo of Section 6.3 computes both.

The Poisson PML of Section 3.7 of Chapter 5 is a lucky special case of QMLE: the Poisson belongs to the linear exponential family, and as long as the conditional mean is correctly specified, the pseudo-true value equals the true value of the mean parameters (Gourieroux, Monfort and Trognon 1984b). The Weibull duration model has no such property: get the distribution wrong, and the fee coefficient and the shape parameter move together. Wooldridge (2010) pp. 996-997 records a related empirical regularity: Weibull shape parameters estimated without covariates tend to be smaller than those estimated with covariates. The theoretical explanation is unobserved heterogeneity: ignoring multiplicative heterogeneity makes the estimate of the shape parameter asymptotically too low (p. 1006). In C2 this mechanism is compounded by the deal-seekers' own negative duration dependence, which pushes the pseudo-true shape parameter lower still; the anchoring paper of Section 5.2 returns to the heterogeneity side from a marketing angle.

The C2 sample has one more result that must be reported plainly: the sandwich interval for the shape parameter, [0.6990, 0.7560], does not cover the pseudo-true value 0.7586, and the estimate is lower by 2.14 sandwich standard errors. A nominal 95% interval should miss about once in twenty anyway, and one sample says nothing about the quality of a variance formula; coverage has to be read from the Monte Carlo. Across 2,000 replications, the shape-parameter estimates average 0.759076, only 0.00044 away from the pseudo-true value, and the coverage rates of the model-standard-error, OPG, and sandwich intervals are 0.774, 0.5805, and 0.9455.

### 3.7 Model Selection: The Vuong Test

C2 also invites another question: of two models that are both wrong, the Weibull and the log-logistic, which is closer to the true distribution? Both have 4 parameters; the log-likelihood is -16034.97 for the Weibull and -16678.82 for a single log-logistic, so per subscriber the Weibull is ahead by 0.0805 on average. Is this difference sampling noise, or is it real? The two models are not nested, so the $\chi^2$ theory of LR does not apply.

::: {.assumption}
**VU.1 (Wooldridge 2010 equations (13.73) and (13.77), pp. 506-507)** Each of the two models satisfies QML.1, and their log-likelihoods at the pseudo-true values are not identical: their difference is nonzero with positive probability, that is, the two models are strictly non-nested.
:::

::: {.theorem}
**Theorem 7 (the Vuong model selection test; Wooldridge 2010 §13.11.2, equations (13.72) to (13.76), pp. 505-506; Vuong 1989)** Let $\Delta\ell_i$ be the difference, for observation $i$, between the log-likelihoods of the two models, each evaluated at its own QMLE, and let $\hat\sigma_{\Delta\ell}$ be its sample standard deviation. Under VU.1 and the null hypothesis that "the two models have equal expected log-likelihoods," $\mathrm{V}_n = n^{-1/2}\sum_{i=1}^n \Delta\ell_i / \hat\sigma_{\Delta\ell} \xrightarrow{d} \mathcal{N}(0, 1)$.
:::

The proof sketch has three steps. Step one (QML.1): expand each model's summed log-likelihood around its pseudo-true value. The first-order term is "the sum of the sample scores" times $\hat\theta - \theta_*$, which are $O_p(\sqrt n)$ and $O_p(n^{-1/2})$ respectively, so after division by $\sqrt n$ the term goes to zero, and estimating the parameters does not affect the limiting distribution of $n^{-1/2}\sum_i \Delta\ell_i$ (equation (13.72)). Step two (VU.1): this quantity is therefore asymptotically equivalent to a sum of independent and identically distributed terms at the pseudo-true values, with mean zero under the null and a variance $\mathrm{Var}(\Delta\ell_i)$ that is positive because the models are strictly non-nested, so the central limit theorem applies (equation (13.74)). Step three: standardizing by the consistent estimate of equation (13.75) gives equation (13.76).

In words: the Vuong test asks which model is closer to the true distribution in the sense of expected KL divergence, not which model is right; rejecting the null says only that the winning model fits better (Wooldridge 2010 p. 507). Misspecification is its very premise: the null itself requires both models to be wrong, since if one were right it would win by a strict inequality. Who sets the weights: as in Theorem 6, the comparison is weighted by the distribution of $x$ and by the censoring design, so with a different export date the conclusion could in principle change. On C2, $\mathrm{V}_n = 38.09$, and the common logarithm of the two-sided $p$ value from the normal distribution is -316.8: the Weibull wins by an overwhelming margin, yet we know it is still wrong. Nested models (for example, the exponential against the Weibull) should not be compared with the Vuong test, because under the null the two models give the same likelihood and $\Delta\ell_i$ degenerates to zero; use the tests of Section 3.5 instead: the LR test under correct specification, and the robust Wald or LM test when the specification is in doubt (Theorem 5 (ii)).

### 3.8 Partial MLE and Grouped Durations

Block C3 splits C1 into 37,819 monthly rows. A binary regression on these rows with month dummies and the complementary log-log link function (cloglog from here on) gives a fee coefficient of 0.02031; Newton iterations applied directly to the likelihood of the grouped durations give estimates that differ from these by $3.2 \times 10^{-13}$, so the two are the same estimator. This section explains why, and uses the example to introduce partial MLE.

Begin with partial MLE in general. The data have $M$ periods (or $M$ intervals), and $\ell_{im}(\theta) = \ln f_m(y_{im} \mid x_{im};\theta)$ is the log-likelihood of period $m$. The partial log-likelihood, here also called the pooled log-likelihood, is defined as $\ell_i(\theta) = \sum_{m=1}^M \ell_{im}(\theta)$ (Wooldridge 2010 equation (13.44)), and partial MLE maximizes its sum over the sample (equation (13.45)). The difference from CMLE is that the conditional density of each period is correctly specified, but their product need not be the joint density of $(y_{i1}, \ldots, y_{iM})$. The name resembles the partial likelihood of Cox (1975), but the meaning is different.

::: {.assumption}
**PL.1 (Wooldridge 2010 §13.8.1, p. 486)** For each $m$, $f_m(\cdot \mid x_m;\theta_0)$ is the true conditional density of $y_{im}$ given $x_{im}$; $\theta_0$ is the unique minimizer of $\mathbb{E}[-\ell_i(\theta)]$; and ML.2 and ML.3 hold for $q = -\ell_i$.

**PL.2 (dynamic completeness; Wooldridge 2010 equation (13.54), p. 492)** For each $m$, $f_m(\cdot \mid x_m;\theta_0)$ is also the conditional density of $y_{im}$ given $x_{im}$ and the entire past $(y_{i,m-1}, x_{i,m-1}, \ldots, y_{i1}, x_{i1})$.
:::

::: {.theorem}
**Theorem 8 (partial MLE; Wooldridge 2010 §13.8.2 and §13.8.3, equations (13.50) to (13.55), pp. 490-493)** (i) Under PL.1, partial MLE is consistent and $\sqrt n$-asymptotically normal with variance $A_0^{-1}B_0A_0^{-1}$, where $A_0 = \sum_m \mathbb{E}[\nabla^2_\theta q_{im}]$, $B_0 = \mathbb{E}[(\sum_m s_{im})(\sum_m s_{im})']$, $q_{im} = -\ell_{im}$, and $s_{im}$ is its score; $B_0$ is estimated by the outer product of the scores clustered by individual (equation (13.51)). (ii) If the scores of the different periods are serially uncorrelated (equation (13.52)), in particular when PL.2 holds, $B_0 = A_0$, and the standard errors and LR tests computed by treating the pooled data as an ordinary likelihood are all valid.
:::

The proof sketch has three steps. Step one (PL.1): by the KL inequality applied period by period, $\theta_0$ maximizes each $\mathbb{E}[\ell_{im}]$ and hence their sum; Theorems 2 and 3 of Chapter 5 give consistency and the sandwich. Step two (PL.1): the conditional information matrix equality holds for each period separately, so $B_0 = A_0 + \sum_{m \neq m'}\mathbb{E}[s_{im}s_{im'}']$, and whether the sandwich collapses depends only on the cross-period terms. Step three (PL.2): dynamic completeness gives $\mathbb{E}[s_{im} \mid x_{im}, y_{i,m-1}, x_{i,m-1}, \ldots] = 0$ (equation (13.55)), and the score of a period $m' < m$ is a function of this conditioning set, so iterated expectations make the cross terms zero, which is (ii).

In words: under PL.1, partial MLE estimates the common parameter of the per-period conditional densities, even if the dependence structure across periods is not written down at all; correlation across periods enters only the variance, which is why clustering is the default. Under misspecification, for instance if the density of some period is wrong, it becomes a QMLE in the sense of Section 3.6, and again only the sandwich can be used (Wooldridge 2010 p. 504). Who sets the weights: the pooled likelihood gives every person-period row equal weight, so people subscribed for many months weigh more, whereas the cluster variance treats each person as a single independent unit. C3 is a case in which PL.2 holds (Theorem 9 below explains why), and the standard error clustered by subscriber is 1.003 times the usual one: no difference. The panel binary models of Chapter 12 will meet cases in which PL.2 fails, and there clustering is indispensable.

Grouped durations are a textbook application of partial MLE. In the notation of Wooldridge (2010) §22.4, the intervals are $[a_{m-1}, a_m)$ with $a_0 = 0$ and $m = 1, \ldots, M$; in C3, $a_m = m$ months and $M = 24$. The indicator $y_{im} = 1$ means that subscriber $i$ cancels in interval $m$, and $m_i$ is the interval in which their spell ends, where a spell ends with either a cancellation or censoring.

::: {.assumption}
**GD.1 (Wooldridge 2010 equations (22.40) and (22.46), pp. 1012-1013)** The hazard is proportional, $h(t \mid x) = \exp(x'\beta)\,h^{\mathrm{base}}(t)$, with no intercept in $x$; given $x$, $t^*$ is independent of censoring; censoring occurs only at the cut points, and cancellations after $a_M$ are recorded as censored at $a_M$.
:::

::: {.theorem}
**Theorem 9 (the MLE for grouped durations is pooled cloglog; Wooldridge 2010 §22.4.1, equations (22.41) to (22.47), pp. 1012-1013; for the cloglog form, see Prentice and Gloeckler 1978)** Under GD.1, the conditional probability of canceling in interval $m$ is

$$h_m(x) = \Pr(a_{m-1} \le t^* < a_m \mid t^* \ge a_{m-1}, x) = 1 - \exp\{-\exp(\eta_m + x'\beta)\}, \qquad \eta_m = \ln\int_{a_{m-1}}^{a_m} h^{\mathrm{base}}(\tilde t)\,d\tilde t;$$

and the grouped log-likelihood of subscriber $i$ (equation (22.45)) equals exactly the sum of the Bernoulli log-likelihoods over their person-period rows,

$$\sum_{m=1}^{m_i - 1 + \delta_i}\big[y_{im}\ln h_m(x_i) + (1 - y_{im})\ln\{1 - h_m(x_i)\}\big], \qquad y_{im} = \delta_i\,\mathbf{1}\{m = m_i\}.$$

Here $\delta_i$ already records cancellations after $a_M$ as censored, following GD.1. So the pooled cloglog binary regression with interval dummies is the MLE for grouped durations; PL.2 holds, and the usual standard errors are valid.
:::

The proof sketch has three steps. Step one (GD.1): under independent censoring, the conditional probability of surviving interval $m$ is $\exp\{-\int_{a_{m-1}}^{a_m} h(\tilde t \mid x)\,d\tilde t\}$ (equations (22.41) and (22.42)), and the proportional form splits the integral into $\exp(x'\beta)$ times the integral of the baseline hazard, which is the $\exp\{-\exp(\eta_m + x'\beta)\}$ of equation (22.47). Step two (GD.1): equation (22.45) is the log of the probability of "surviving the first $m_i - 1$ intervals," plus, if a cancellation is observed, the log of the conditional probability of canceling in interval $m_i$. The first $m_i - 1$ rows have $y_{im} = 0$ and each contributes $\ln\{1 - h_m\}$, while the row of the cancellation contributes $\ln h_{m_i}$, so the two agree term by term. Step three (PL.2): the density of each row conditions on "surviving to the start of the interval," and the past of a spell has only one possibility, "having survived every earlier interval," so this is also the conditional density given the entire past, and Theorem 8 (ii) applies.

In words: under GD.1, the $\beta$ estimated by pooled cloglog is the $\beta$ of the continuous-time proportional hazard model. C3's fee coefficient, 0.02031, aims at 0.02; compared with C1's continuous-time MLE it is lower by 0.00025. The term $\eta_m$ estimates the log of the integral of the baseline hazard over each interval and needs no shape assumption like the Weibull's, and that is the advantage of the grouped formulation. Under misspecification, if the hazard is not proportional (for example, if the effect of the fee decays with tenure), cloglog estimates a pseudo-true value that is an average of sorts; replacing cloglog with logit gives a different discrete-time model, whose coefficients are no longer semi-elasticities of the continuous-time hazard. Who sets the weights: the precision of $\eta_m$ is set by the number of cancellations in that month. Months 23 and 24 have only 6 cancellations each, the standard error of $\hat\eta_m$ is about 0.415, and the estimates for these two months are mostly noise. In month 21, $\hat\eta_m$ lies more than three standard errors above the truth. A deviation this large is unusual even among 24 estimates, but the baseline hazard of the DGP declines smoothly through month 21 without any jump, so it can only be sampling error and should not be read as a peak in cancellations in month 21.

### 3.9 Two-Step MLE

The first step in block C4 is a Poisson CMLE: regressing the check-ins $v_i$ on $z_i$ gives $\hat\gamma = (-1.0086, 0.2534, 0.3438)$, against true values $(-1.0, 0.25, 0.35)$. The second step puts the fitted index $z_i'\hat\gamma$ into the Weibull hazard as a covariate and gives $\hat\beta_{\mathrm{use}} = -0.8181$, against a truth of -0.80. Treating $z_i'\hat\gamma$ as data, the software reports a standard error of 0.04120; once the first step's estimation error is accounted for, the standard error is 0.05006, 1.215 times as large.

This is a two-step estimator whose second step is an MLE, and its theory is simply Theorem 6 of Section 3.6 of Chapter 5, with $q$ taken to be the negative log-likelihood of the second step (Wooldridge 2010 §13.10.1, p. 499). The first step's influence function comes from the Poisson CMLE, $\psi_i = (\mathbb{E}[\exp(z_i'\gamma_0)\,z_iz_i'])^{-1}z_i\{v_i - \exp(z_i'\gamma_0)\}$ (equation (13.65)). The matrix $\Gamma_0 = \mathbb{E}[\nabla_\gamma s(w,\theta_0,\gamma_0)]$ contains the factor $\beta_{\mathrm{use}}$: if the index did not affect the hazard, the first step's error would have no way to propagate. In C4, $\beta_{\mathrm{use}} = -0.80 \neq 0$; given $z_i$, $v_i$ is independent of the latent duration, so the scores of the two steps are uncorrelated, which is the case of Theorem 6 (ii) of Chapter 5, and the corrected variance cannot be smaller than the naive one. Yet Wooldridge, in §13.10.2, gives a result that seems to say the opposite.

Write the first-step log-likelihood as $\ell_i^{(1)}(\gamma)$, and let $\mathbb{E}_\gamma$ denote an expectation taken under the first-step model's distribution of $v$ given $z$.

::: {.theorem}
**Theorem 10 (efficiency when the first step is a CMLE; Wooldridge 2010 §13.10.2, equations (13.64) to (13.69), pp. 500-501)** Suppose observation $i$ is split into $(v_i, w_i, z_i)$, where $w_i$ denotes only the second-step data other than $v_i$ and $z_i$ (in C4, tenure, the censoring indicator, and the fee); the second-step objective $q(v_i, w_i, z_i, \theta, \gamma)$ contains the first-step data $v_i$ and satisfies TS.1 through TS.3 of Chapter 5; the first step is a correctly specified CMLE of $v_i$ given $z_i$; and (a) conditional independence holds: given $z_i$, $v_i$ is independent of $w_i$ (equation (13.66)); (b) equation (13.67) holds: $\mathbb{E}[\nabla_\gamma s_i(\theta_0,\gamma_0) \mid w_i, z_i] = -\mathbb{E}[s_i\,\nabla_\gamma \ell_i^{(1)}(\gamma_0) \mid w_i, z_i]$. Then $s_i + \Gamma_0\psi_i$ is the residual from the population linear regression of $s_i$ on the first-step score $\nabla_\gamma \ell_i^{(1)}(\gamma_0)'$, and the asymptotic variance is $A_0^{-1}\mathrm{Var}(s_i + \Gamma_0\psi_i)A_0^{-1}$ (equations (13.68) and (13.69)); $B_0 - \mathrm{Var}(s_i + \Gamma_0\psi_i)$ is positive semidefinite, so the two-step estimator that estimates $\gamma$ by CMLE is at least as efficient as plugging in the true $\gamma_0$, and the sandwich standard error $\hat A^{-1}\hat B\hat A^{-1}/n$ computed while ignoring the first step is conservative.
:::

The proof sketch has three steps. Step one (the first step is a CMLE): by the information matrix equality, $\psi_i = \{\mathbb{E}[\nabla_\gamma \ell_i^{(1)\prime}\nabla_\gamma \ell_i^{(1)}]\}^{-1}\nabla_\gamma \ell_i^{(1)}(\gamma_0)'$ (equation (13.65)). Step two ((a) and (b)): taking expectations over $(w_i, z_i)$ gives $\Gamma_0 = -\mathbb{E}[s_i\,\nabla_\gamma \ell_i^{(1)}]$, and substituting into the $s_i + \Gamma_0\psi_i$ of Theorem 6 of Chapter 5 yields the residual of a linear projection. Step three: the variance of a residual is no larger than the variance of the original variable, that is, $\mathrm{Var}(s_i + \Gamma_0\psi_i) \le B_0$; the two are equal if and only if $s_i$ is uncorrelated with the first-step score, in which case $\Gamma_0 = 0$.

In words: under the conditions of Theorem 10, the two-step estimator still estimates $\theta_0$, and far from adding noise, the first step removes the part of the noise in $s_i$ that the first-step score can explain. The key is condition (b): it is equivalent to the derivative of $\mathbb{E}_\gamma[s_i(\theta_0,\gamma) \mid w_i, z_i]$ with respect to $\gamma$ being zero at $\gamma_0$, that is, to the effect of $\gamma$ on the second-step score being exactly offset by the change in the distribution of $v$. Wooldridge derives equation (13.67) from equations (13.66) and (13.39); the derivation of equation (13.39) requires that conditional expectations taken under the model do not change with the parameter (equation (13.40)), and this chapter states that requirement separately as condition (b). Inverse probability weighting is an example: if Sorrel observed durations only for people who completed a survey, with $v_i$ indicating whether the survey was completed, and the second step weighted by the inverse of the estimated response probability, with response depending only on $z_i$, then estimating the response probability would actually be more efficient than using the true probability. C4 does not satisfy (b): the second-step likelihood does not contain $v_i$ at all, and $\gamma$ enters only through the index $z_i'\gamma$, so there is nothing to offset. Under misspecification, the first step converges to $\gamma_*$, the theorem fails, and we are back to the general formula of Theorem 6 of Chapter 5. Who sets the weights: each first-step error $\psi_i$ enters the second step only after being multiplied by $\Gamma_0$; $\Gamma_0$ is proportional to $\beta_{\mathrm{use}}$ and weighted by the cumulative hazard, so people with a large cumulative hazard amplify the first-step error the most. On C4 the corrected standard error is 0.05006 and the naive one 0.04120; across 1,000 Monte Carlo replications the coverage rate of the naive interval is 0.903, that of the corrected interval 0.951, and that of the infeasible estimator using the true index 0.941.

## 4 Estimation: Algorithms, Variances, and Software

### 4.1 Numerical Optimization and Parameterization

MLE generally has no closed-form solution, and the workhorse is still the Newton-Raphson of Section 4.1 of Chapter 5: at the current point $\theta^{(j)}$, use the score and the Hessian to approximate $Q_n$ by a quadratic, and jump to its lowest point. The Weibull's score and Hessian both have analytic expressions, and the Hessian of the $\beta$ block, $\sum_i \exp(x_i'\beta)\,t_i^{\kappa}x_ix_i'$, is positive definite everywhere. On C1, starting from the rate of the exponential distribution, the iteration converges in 5 steps; `survival::survreg` reaches the same point, and the two differ by $8.9 \times 10^{-15}$.

Two practical details deserve attention. The first is parameterization: the shape parameter enters as $\ln\kappa$, so the parameter space has no boundary and the iterations cannot wander into negative values. MLE is invariant to reparameterization, $\hat\kappa = \exp(\widehat{\ln\kappa})$, and the standard error comes from the Delta method, $\mathrm{se}(\hat\kappa) = \hat\kappa\,\mathrm{se}(\widehat{\ln\kappa})$. The second is scale: with the fee measured in dollars, the coefficients are on the order of hundredths while the intercept is in the single digits, so the Hessian has a large condition number, and the default step halving of `maxLik` stalls where the gradient is small but not zero. Switching to the Marquardt modification and judging convergence by the gradient rather than by the function value makes the iterations stable. Measuring the fee in units of USD 10 is another common remedy. After convergence, check that the gradient really is zero and the Hessian positive definite, and re-estimate from a different starting point.

`survreg` reports the AFT form of Section 2.1: it estimates the coefficients of $\ln t$ on $x$ and a scale, and the conversion back to the proportional hazard form is $\kappa = 1/\mathrm{scale}$, $\beta = -\mathrm{coef}/\mathrm{scale}$. On C1 the AFT fee coefficient is -0.02323 and the scale 1.1301, which convert to the hazard coefficient 0.02056 and the shape 0.8848. The coefficients of the two forms have opposite signs, so before reading the results, confirm which form the software uses; Stata's `streg` reports the proportional hazard form by default.

### 4.2 Variance Estimators and Software

| Variance estimator | Bread | Filling | When it is valid |
|---|---|---|---|
| Hessian (equation (13.32), first form) | Average of the sample Hessians | Same as the bread (collapsed) | ML.1 through ML.3 |
| Expected Hessian (equation (13.31)) | Expected Hessian given $x$ | Same as the bread (collapsed) | As above, and the conditional expectation has a closed form |
| OPG / BHHH (equation (13.32), second form) | Average of the outer products of the scores | Same as the bread (collapsed) | As above; often poor in finite samples |
| Sandwich (equation (13.71)) | Average of the sample Hessians | Average of the outer products of the scores | QML.1, including misspecification |
| Cluster sandwich (equation (13.51)) | Average of the sample Hessians | Outer product of the scores summed by individual | Partial MLE with scores correlated across periods |
| Two-step correction | Second-step Hessian | Outer product after adding $\hat\Gamma\hat\psi_i$ | When the first step matters ($\Gamma_0 \neq 0$) |

In software, `maxLik::maxLik(logLik, grad, hess, start, method = "NR")` is a general-purpose MLE engine, and `method = "BHHH"` requires the log-likelihood to return a vector with one entry per observation; `sandwich::sandwich()` gives equation (13.71) for `maxLik` objects; `numDeriv::hessian` and `jacobian` can be used to check hand-written derivatives, and when both packages are loaded, write `maxLik::hessian`, or it will be masked by the function of the same name in `numDeriv`. `lmtest::lrtest(m0, m1)` computes LR, which is meaningful only under correct specification. For grouped durations use `glm(y ~ 0 + factor(m) + x, family = binomial("cloglog"))`, with `sandwich::vcovCL` for the cluster variance. For Tobit, use `censReg::censReg` or `AER::tobit`.

## 5 Anchoring Papers

Each of the three anchoring papers corresponds to one thread of the chapter: what the likelihood estimates under misspecification, how retention in subscription and contractual settings is modeled with a likelihood, and the standard treatment of duration data and hazards in economics. We have not yet found an anchoring application of survival analysis or churn modeling in IS journals, so the second paper comes from marketing and the third is a survey in the Journal of Economic Literature. Each paper is laid out in five parts: paper, method, data, results, and limitations.

### 5.1 White (1982)

::: {.case}
Paper: "Maximum Likelihood Estimation of Misspecified Models," Econometrica. It studies the consequences and the diagnosis of model misspecification when maximum likelihood is used for estimation and inference.

Method: without assuming that the model family contains the true distribution, the paper studies the limit of the QMLE. In White's formulation, the QMLE converges to a well-defined limit, which may or may not be consistent for the parameters the researcher cares about. White characterizes this limit as the best approximation in the KL sense (as paraphrased in Wooldridge 2010 p. 503); this is Theorem 6 of this chapter. Around this limit the estimator is asymptotically normal with a sandwich variance; the conventional Wald, LM, and LR tests break down under misspecification, and more general statistics are needed. The paper also uses the properties of the QMLE and of the information matrix to construct specification tests, the best known of which is the information matrix test: it compares two estimates of the information matrix, one from the Hessian and one from the outer product of the scores, and a large difference signals a problem with the specification.

Data: this is a theoretical paper, with no empirical data.

Results: the paper establishes the consistency, asymptotic normality, and sandwich variance of the QMLE, together with a family of inference methods and specification tests that remain valid under misspecification. On C2 the ratio of the model standard error to OPG is 1.484, which is the kind of divergence the information matrix test formalizes.

Limitations: the pseudo-true value is a best approximation in a statistical sense and need not correspond to any economically meaningful quantity; "inference is possible" does not mean "inference is worthwhile." When computed by an outer-product regression, tests of this kind often perform poorly in finite samples, and Wooldridge (2010) p. 484 recommends forms based on the expected Hessian in specific models instead.
:::

This paper turned "all models are wrong" from a slogan into a workable theory: misspecification does not make the estimator meaningless; it only changes what the estimator estimates, and inference must then switch to the sandwich.

### 5.2 Fader and Hardie (2007)

::: {.case}
Paper: "How to Project Customer Retention," Journal of Interactive Marketing. The question is how a subscription or contractual business can extrapolate a few years of retention data into the future in order to estimate customers' expected lifetimes and value. The numbers and page numbers below follow the journal version.

Method: the shifted-beta-geometric (sBG) model. Each customer cancels in each period with a fixed probability, and this probability follows a beta distribution across customers. The probability of canceling in a given period and the survival function then both have closed forms, and the aggregate retention rate rises with customer tenure even though each person's cancellation probability is constant. This is an illusion created by heterogeneity: those with high cancellation probabilities leave first, and those who remain are more loyal on average. The parameters are estimated by maximum likelihood: the customers who cancel in each of the first 7 years contribute the corresponding cancellation probability, and those still active after 7 years contribute the survival function. This is the grouped-duration likelihood of Section 3.8 of this chapter, with "those still active" as the right-censored observations. The authors carry out all the computations in Excel.

Data: 12 years of retention curves for two segments (regular and high-end) of a subscription business, from the book by Berry and Linoff (2004); they consist only of the percentage still active in each year. In year 7, 24.1% of regular customers and 49.1% of high-end customers are still active; in year 12, the shares are 17.3% and 39.4%.

Results: fitted to the first 7 years only, the sBG projection of year 12 retention underestimates it by only 4% (high-end) and 2% (regular); the paper's text calls these overestimates, but the model curves in its Figure 4 and the values computed from its parameter estimates both lie below the actual values. For comparison, linear and exponential regressions on the same data underestimate the high-end customers' year 12 retention by 81% and 30%, and a quadratic regression overestimates it by 92%; for regular customers, the exponential regression underestimates it by 54%. The estimated mean cancellation probability is 0.37 for regular customers and 0.15 for high-end customers.

Limitations: the model has no covariates, so it can describe only one group and cannot answer "what happens if the monthly fee is USD 1 higher?" That individual cancellation probabilities do not change over time is an assumption; when the data are only an aggregate retention curve, it is hard to separate from "individuals become ever more loyal" without additional structural assumptions. For the high-end customers' first beta parameter, the paper's text (p. 82) gives an estimate of 0.688, while Appendix B (p. 90) gives 0.668; the maximum likelihood estimate consistent with the data in Table 1 is 0.668, so the 0.688 in the text is a typo.
:::

This paper and C2 address the same problem: a rising aggregate retention rate (or a falling aggregate hazard) can come entirely from differences between people, which is how the sBG is specified; in C2, heterogeneity is compounded by the deal-seekers' own falling hazard, and the aggregate curve alone cannot tell the two apart. Fitting a model that ignores heterogeneity yields a pseudo-true value that looks reasonable but extrapolates badly.

### 5.3 Kiefer (1988)

::: {.case}
Paper: "Economic Duration Data and Hazard Functions," Journal of Economic Literature. It is a survey of duration analysis written for economists.

Method: this box is written from Wooldridge's (2010) two citations of the survey. The survey's subject is duration data and hazard functions in economics. Wooldridge (2010) cites two of its specific treatments: handling grouped durations as a sequence of binary outcomes (p. 1010), which is the idea behind Theorem 9 of this chapter, and the name length-biased sampling for the selection problem that arises when only people still in the state at some point in time are sampled (p. 1000).

Data: this is a survey, focused on methods.

Results: as a survey, its contribution is a methodological framework for handling duration data, not an estimate.

Limitations: the survey is fairly old, and methods for unobserved heterogeneity, multiple spells, and similar problems have advanced considerably since, so it should be read alongside more recent material such as Wooldridge (2010) §22; as a survey, it does not defend the identifying assumptions of any particular empirical study.
:::

Kiefer's survey is worth reading because it organizes its discussion of economic duration data around the hazard function, and the two core moves of this chapter are to describe behavior with the hazard and to handle "partly seen" observations with the likelihood.

Taken together: White shows that under misspecification the likelihood estimates a best approximation and inference must switch to the sandwich; Fader and Hardie show that in subscription businesses, a model that ignores heterogeneity yields just such a misleading best approximation; Kiefer's survey is organized around duration data and hazard functions and discusses both grouped durations and length-biased sampling. All three point to the chapter's prescription: write the observation mechanism into the likelihood, leave inference to the sandwich, and put the assumptions about heterogeneity on the table.

## 6 A Full Walkthrough on Sorrel Data

The code below runs under R 4.6.1. It shows only the statements that generate the data, estimate the models, compute variances, and summarize the Monte Carlo, with several statements merged onto one line. Every number in the text comes from the output of an actual run, and the result lines in the code that start with `#` copy what the run printed. `weibull_parts`, `fit_weibull`, `gen_C2`, `gen_C4`, `grouped_data`, `grouped_mle`, `loglogistic_ll_obs`, and `two_step_vcov` are helper functions whose roles are described in the code comments.

### 6.1 The DGP

The key settings of C1 are these: signup months are drawn with 4% monthly growth, the censoring time is the number of months until the export date, and latent durations are obtained by inverting the Weibull survival function.

```r
P1 <- list(n = 5000, months = 36, growth = 1.04, fees = c(45, 60, 75, 90), p_fee = c(0.30, 0.35, 0.20, 0.15),
           p_ref = 0.30, beta = c(-3.44, 0.02, -0.40), kappa = 0.9)
gen_C1 <- function(par, n = par$n) {
  signup <- sample.int(par$months, n, replace = TRUE, prob = par$growth^seq_len(par$months))
  fee <- sample(par$fees, n, replace = TRUE, prob = par$p_fee)
  referred <- rbinom(n, 1, par$p_ref)
  c_months <- par$months - signup + 1                     # censoring time: whole months until the data pull
  eta <- par$beta[1] + par$beta[2] * fee + par$beta[3] * referred
  t_star <- (rexp(n) / exp(eta))^(1 / par$kappa)          # S(t | x) = exp(-exp(x'beta) t^kappa)
  data.frame(fee_usd = fee, referred = referred, c_months = c_months,
             t_months = pmin(t_star, c_months), delta = as.numeric(t_star <= c_months), t_star_true = t_star)
}
set.seed(6101); C1 <- gen_C1(P1)
cat(sprintf("# C1: censored share %.4f; median tenure, 60 USD without referral: %.2f months (truth)\n",
            1 - mean(C1$delta), (log(2) / exp(P1$beta[1] + 60 * P1$beta[2]))^(1 / P1$kappa)))
# C1: censored share 0.4332; median tenure, 60 USD without referral: 8.02 months (truth)
```

C2 and C4 share the same signup and export design. In C2 the probability of being a deal-seeker is `plogis(g0 + g_ref * referred)`, and the durations of both groups are log-logistic, with the median multiplied by the fee's scaling factor; in C4 the check-ins are Poisson, and the hazard contains the true index.

```r
P2 <- modifyList(P1[c("months", "growth", "fees", "p_fee", "p_ref")],
                 list(n = 8000, g0 = -0.5, g_ref = -1.0, med_deal = 1.5, p_deal = 0.5, med_reg = 12, p_reg = 5,
                      theta_fee = 0.02, fee_ref = 60, R_mc = 2000))
P4 <- modifyList(P1[c("months", "growth", "fees", "p_fee")],
                 list(n = 5000, hh_sizes = 1:5, p_hh = c(0.15, 0.35, 0.20, 0.20, 0.10), gamma = c(-1.0, 0.25, 0.35),
                      beta = c(-3.70, 0.02, -0.80), kappa = 0.9, R_mc = 1000))
set.seed(6102); C2 <- gen_C2(P2)       # deal-seekers: log-logistic, median 1.5, shape 0.5; regulars: 12 and 5
set.seed(6104); C4 <- gen_C4(P4)       # v ~ Poisson(exp(z'gamma)); hazard uses the true index z'gamma
```

### 6.2 Block C1: Censoring, MLE, and the Trinity

First run the naive regression, then convert it to the hazard scale as in Section 2.1.

```r
X <- cbind(1, C1$fee_usd, C1$referred); t <- C1$t_months; d <- C1$delta
naive <- lm(log(t) ~ fee_usd + referred, data = C1)      # censored spells treated as exits
b <- unname(coef(naive)[2]); s2 <- mean(residuals(naive)^2)
kap_ols <- pi / sqrt(6 * s2)                             # Var(log t*) = pi^2 / (6 kappa^2) under Weibull
latent <- lm(log(t_star_true) ~ fee_usd + referred, data = C1)   # infeasible: the uncensored t*
cat(sprintf("# naive OLS: slope %.5f (SE %.5f), kappa %.3f, beta_fee %.5f; on t*: slope %.5f; truth -0.02/0.9 = %.5f\n",
            b, sqrt(sandwich::vcovHC(naive, type = "HC0")[2, 2]), kap_ols, -kap_ols * b,
            coef(latent)[[2]], -0.02 / 0.9))
# naive OLS: slope -0.01512 (SE 0.00119), kappa 1.014, beta_fee 0.01533; on t*: slope -0.02319; truth -0.02/0.9 = -0.02222
```

The MLE uses the Newton-Raphson of `maxLik` with analytic derivatives; `info_hat` computes $\hat A$ and $\hat B$ independently from numerical derivatives, and the three standard errors come from these.

```r
nll_weibull_cens <- function(par, t, delta, X, by_obs = FALSE) {   # sum of q_i = -l_i; par = (beta', log kappa)'
  k <- ncol(X); rho <- par[k + 1]; eta <- drop(X %*% par[1:k])
  q <- -(delta * (eta + rho + (exp(rho) - 1) * log(t)) - exp(eta + exp(rho) * log(t)))
  if (by_obs) q else sum(q)
}
wp <- function(p) weibull_parts(p, t, d, X)              # analytic q_i, scores s_i' (rows) and summed Hessian
fm <- maxLik::maxLik(logLik = function(p) -wp(p)$q, grad = function(p) -wp(p)$S, hess = function(p) -wp(p)$H,
                     start = c(log(sum(d) / sum(t)), 0, 0, 0), method = "NR",
                     control = list(tol = -1, reltol = -1, gradtol = 1e-7, iterlim = 200, qac = "marquardt"))
p <- unname(fm$estimate)
info_hat <- function(par, nll, ...) {                    # A_hat: mean Hessian; B_hat: mean outer product of scores
  n <- length(nll(par, ..., by_obs = TRUE))
  A <- numDeriv::hessian(function(q) nll(q, ...), par, method.args = list(eps = 1e-2, d = 0.1, r = 4, v = 2)) / n
  S <- numDeriv::jacobian(function(q) nll(q, ..., by_obs = TRUE), par)
  list(A = (A + t(A)) / 2, B = crossprod(S) / n)
}
ih <- info_hat(p, nll_weibull_cens, t = t, delta = d, X = X); n <- nrow(X); Ai <- solve(ih$A)
se <- sqrt(cbind(hessian = diag(Ai), opg = diag(solve(ih$B)), sandwich = diag(Ai %*% ih$B %*% Ai)) / n)
cat(sprintf("# MLE: beta_fee %.5f, kappa %.4f (%d NR steps); SE of beta_fee: Hessian %.6f, OPG %.6f, sandwich %.6f\n",
            p[2], exp(p[4]), fm$iterations, se[2, 1], se[2, 2], se[2, 3]))
# MLE: beta_fee 0.02056, kappa 0.8848 (5 NR steps); SE of beta_fee: Hessian 0.001202, OPG 0.001223, sandwich 0.001183
sv <- survival::survreg(survival::Surv(t_months, delta) ~ fee_usd + referred, data = C1, dist = "weibull",
                        control = survival::survreg.control(rel.tolerance = 1e-12, maxiter = 100))
ph <- c(-coef(sv) / sv$scale, -log(sv$scale))            # AFT -> PH: kappa = 1 / scale, beta = -coef / scale
cat(sprintf("# survreg: fee %.5f (AFT), scale %.4f; converted to PH, max |difference| from the MLE %.1e\n",
            coef(sv)[[2]], sv$scale, max(abs(ph - p))))
# survreg: fee -0.02323 (AFT), scale 1.1301; converted to PH, max |difference| from the MLE 8.9e-15
```

The trinity: the restricted model fixes $\ln\kappa = 0$, that is, the exponential distribution.

```r
trinity <- function(par_u, par_r, nll, restr, ...) {     # H0: R(theta) = 0; Hessian form of the information
  f <- function(q) nll(q, ...); ma <- list(eps = 1e-2, d = 0.1, r = 4, v = 2)
  r_u <- restr(par_u); Rj <- numDeriv::jacobian(restr, par_u)
  W <- drop(t(r_u) %*% solve(Rj %*% solve(numDeriv::hessian(f, par_u, method.args = ma)) %*% t(Rj), r_u))
  LR <- 2 * (f(par_r) - f(par_u))                        # 2 (l_u - l_r)
  g_r <- numDeriv::grad(f, par_r)                        # sum of the scores at the restricted estimate
  LM <- drop(t(g_r) %*% solve(numDeriv::hessian(f, par_r, method.args = ma), g_r))
  c(W = W, LR = LR, LM = LM)
}
fr <- fit_weibull(t, d, X, log_kappa = 0)                # restricted MLE: kappa = 1, the exponential model
tri <- trinity(p, fr$par, nll_weibull_cens, function(q) q[4], t = t, delta = d, X = X)
cat(sprintf("# H0 kappa = 1: W %.3f, LR %.3f, LM %.3f; largest gap %.2f%% of LR; 5%% critical value %.3f\n",
            tri[["W"]], tri[["LR"]], tri[["LM"]], 100 * max(abs(tri[c("W", "LM")] - tri[["LR"]])) / tri[["LR"]],
            qchisq(0.95, 1)))
# H0 kappa = 1: W 63.733, LR 67.621, LM 63.932; largest gap 5.75% of LR; 5% critical value 3.841
```

The log-likelihood at the MLE is -9888.23, and the estimated median tenure of non-referred users in the USD 60 tier is 8.30 months, against a truth of 8.02.

### 6.3 Block C2: QMLE, Vuong, and Coverage

```r
X2 <- cbind(1, C2$fee_usd, C2$referred)
f2 <- fit_weibull(C2$t_months, C2$delta, X2)             # the same Weibull likelihood, now misspecified
ih2 <- info_hat(f2$par, nll_weibull_cens, t = C2$t_months, delta = C2$delta, X = X2)
k2 <- exp(f2$par[4]); A2i <- solve(ih2$A); n2 <- nrow(X2)
se_k <- k2 * sqrt(c(model = A2i[4, 4], opg = solve(ih2$B)[4, 4], sandwich = (A2i %*% ih2$B %*% A2i)[4, 4]) / n2)
cat(sprintf("# QMLE kappa %.4f; SE model %.5f, OPG %.5f, sandwich %.5f (model / OPG %.3f); sandwich CI [%.4f, %.4f]\n",
            k2, se_k[1], se_k[2], se_k[3], se_k[1] / se_k[2], k2 - qnorm(0.975) * se_k[3], k2 + qnorm(0.975) * se_k[3]))
# QMLE kappa 0.7275; SE model 0.00934, OPG 0.00629, sandwich 0.01453 (model / OPG 1.484); sandwich CI [0.6990, 0.7560]
sl <- survival::survreg(survival::Surv(t_months, delta) ~ fee_usd + referred, data = C2, dist = "loglogistic",
                        control = survival::survreg.control(rel.tolerance = 1e-12, maxiter = 100))
dl <- -nll_weibull_cens(f2$par, C2$t_months, C2$delta, X2, by_obs = TRUE) -
      loglogistic_ll_obs(sl, C2$t_months, C2$delta, X2)  # per-subscriber log-likelihood difference
V_n <- sqrt(n2) * mean(dl) / sqrt(mean((dl - mean(dl))^2))
cat(sprintf("# Vuong: mean difference %.4f, V_n %.2f, log10 p %.1f\n",
            mean(dl), V_n, (log(2) + pnorm(-abs(V_n), log.p = TRUE)) / log(10)))
# Vuong: mean difference 0.0805, V_n 38.09, log10 p -316.8
```

The $p$ value lies far below the normal range of double-precision floating-point numbers, and computing it directly would lose precision, so it is computed in logs. In each Monte Carlo replication we redraw a city of 8,000 subscribers and refit the Weibull, then record how often each of the three intervals covers the pseudo-true value.

```r
k_star <- 0.7586356875395922; b_star <- 0.01832457352899277   # pseudo-true kappa_*, beta_fee,*: 2,000,000 draws
one_rep <- function() {
  dd <- gen_C2(P2); Xr <- cbind(1, dd$fee_usd, dd$referred)
  f <- fit_weibull(dd$t_months, dd$delta, Xr, start = f2$par, engine = "newton")
  k <- exp(f$par[4])
  c(kappa = k, model = k * sqrt(f$V_hess[4, 4]), opg = k * sqrt(f$V_opg[4, 4]), sandwich = k * sqrt(f$V_sand[4, 4]),
    fee = f$par[2], fee_model = sqrt(f$V_hess[2, 2]), fee_opg = sqrt(f$V_opg[2, 2]), fee_sandwich = sqrt(f$V_sand[2, 2]))
}
set.seed(6122); sim <- t(replicate(P2$R_mc, one_rep()))
cover <- function(e, s, v) mean(abs(sim[, e] - v) <= qnorm(0.975) * sim[, s])
cat(sprintf("# 2000 reps: mean kappa_hat %.6f, sd %.5f; coverage model %.4f, OPG %.4f, sandwich %.4f\n",
            mean(sim[, "kappa"]), sd(sim[, "kappa"]), cover("kappa", "model", k_star), cover("kappa", "opg", k_star),
            cover("kappa", "sandwich", k_star)))
# 2000 reps: mean kappa_hat 0.759076, sd 0.01583; coverage model 0.7740, OPG 0.5805, sandwich 0.9455
cat(sprintf("# beta_fee: sd %.6f, mean model SE %.6f; coverage model %.4f, OPG %.4f, sandwich %.4f\n",
            sd(sim[, "fee"]), mean(sim[, "fee_model"]), cover("fee", "fee_model", b_star), cover("fee", "fee_opg", b_star),
            cover("fee", "fee_sandwich", b_star)))
# beta_fee: sd 0.000808, mean model SE 0.000936; coverage model 0.9755, OPG 0.9925, sandwich 0.9570
```

The true standard deviation of the shape-parameter estimates, 0.01583, is 1.62 times the average model standard error, 0.00977, while the average sandwich standard error, 0.01556, is almost equal to it. The fee coefficient is the opposite case: in the same replications its true standard deviation is 0.000808, the model standard error averages 0.000936, which is too large, and its coverage is 0.9755; OPG covers at 0.9925 and the sandwich at 0.957. The direction in which the model standard error errs depends on the parameter, and only the sandwich is right in both places. The Monte Carlo standard error of a coverage rate is about 0.0049.

![Monte Carlo coverage of nominal 95% intervals. The top panel shows three standard errors of block C2's Weibull QMLE for the pseudo-true value $\kappa_*$ (2,000 samples of 8,000 subscribers each); the bottom panel shows four standard errors of block C4's two-step MLE for $\beta_{\mathrm{use}}$ (1,000 samples of 5,000 subscribers each), where the oracle uses the true usage-intensity index and is infeasible. The vertical line is 0.95, and the gray band is the 95% range of Monte Carlo error in the coverage rate for each panel's number of replications.](assets/fig/fig_06_coverage.svg)

### 6.4 Block C3: Grouped Durations

`person_period` expands each subscriber into one row per month subscribed, and `glm` runs the pooled cloglog; `grouped_mle` runs Newton iterations directly on equation (22.45) and never builds these rows.

```r
person_period <- function(t, delta, X, cuts) {           # one row per month entered; y = 1 in the exit month
  M <- length(cuts) - 1
  exit_in <- delta == 1 & t <= cuts[M + 1]               # exits after month M are censored at M
  mi <- ifelse(exit_in, findInterval(t, cuts, left.open = TRUE), pmin(findInterval(t, cuts) - 1, M))
  id <- rep(seq_along(t), mi); m <- sequence(mi)
  data.frame(id = id, m = m, y = as.numeric(exit_in[id] & m == mi[id]), X[id, , drop = FALSE])
}
C3 <- person_period(C1$t_months, C1$delta, cbind(fee_usd = C1$fee_usd, referred = C1$referred), 0:24)
gl <- glm(y ~ 0 + factor(m) + fee_usd + referred, family = binomial(link = "cloglog"), data = C3,
          control = glm.control(epsilon = 1e-14, maxit = 500))
gd <- grouped_data(C1$t_months, C1$delta, 0:24)          # the interval form, Wooldridge (22.45) with (22.47)
gm <- grouped_mle(gd$mi, gd$yi, cbind(C1$fee_usd, C1$referred), M = 24)
se_cl <- sqrt(diag(sandwich::vcovCL(gl, cluster = ~id, type = "HC0", cadjust = FALSE)))[["fee_usd"]]
cat(sprintf("# %d rows, %d exits; max |cloglog - grouped MLE| %.1e; beta_fee %.5f (SE %.5f, clustered / usual %.3f)\n",
            nrow(C3), sum(C3$y), max(abs(unname(coef(gl)) - gm$par)), coef(gl)[["fee_usd"]],
            sqrt(vcov(gl)["fee_usd", "fee_usd"]), se_cl / sqrt(vcov(gl)["fee_usd", "fee_usd"])))
# 37819 rows, 2784 exits; max |cloglog - grouped MLE| 3.2e-13; beta_fee 0.02031 (SE 0.00121, clustered / usual 1.003)
at_risk <- tabulate(C3$m, 24); exits <- tapply(C3$y, C3$m, sum); se_eta <- sqrt(diag(gm$V))
cat(sprintf("# month 1: %d at risk, %d exits; month 24: %d, %d; eta_1 %.3f; eta_21 %.3f (SE %.3f); beta_ref %.4f (SE %.4f)\n",
            at_risk[1], exits[[1]], at_risk[24], exits[[24]], gm$par[1], gm$par[21], se_eta[21], gm$par[26], se_eta[26]))
# month 1: 5000 at risk, 497 exits; month 24: 211, 6; eta_1 -3.470; eta_21 -3.216 (SE 0.181); beta_ref -0.3978 (SE 0.0435)
```

Month 1 has 5,000 subscribers at risk and 497 cancellations; by month 24 only 211 remain at risk, with 6 cancellations. The estimate is $\hat\eta_1 = -3.470$ against a truth of -3.44; for month 21, $\hat\eta_{21} = -3.216$ (standard error 0.181) against a truth of -3.847, which is the sampling deviation described in Section 3.8. The referral coefficient is -0.3978 (standard error 0.0435), against a truth of -0.40.

### 6.5 Block C4: Two-Step MLE

`twostep` runs the Poisson CMLE in the first step and stores the influence function, then fits the Weibull in the second step; `F2` is $\hat\Gamma$ computed analytically, row $i$ of `R1` is $\hat\psi_i'$, and `two_step_vcov`, identical to the function of the same name in Section 6.4 of Chapter 5, replaces the filling with the outer product of $\hat s_i + \hat\Gamma\hat\psi_i$.

```r
twostep <- function(dd, start = NULL) {
  Z <- cbind(1, dd$hh_size, dd$quiz_std); n <- nrow(Z)
  g1 <- glm.fit(Z, dd$usage_wk1, family = poisson(), control = list(epsilon = 1e-12, maxit = 100))   # step 1: CMLE
  gam <- unname(g1$coefficients); uh <- drop(Z %*% gam); mu <- exp(uh)
  R1 <- (Z * (dd$usage_wk1 - mu)) %*% solve(crossprod(Z * mu, Z) / n)    # row i: psi_i' (eq. 13.65)
  X <- cbind(1, dd$fee_usd, uh)                                          # generated regressor z'gamma_hat
  f <- fit_weibull(dd$t_months, dd$delta, X, start = start, engine = "newton")   # step 2: MLE
  p <- f$par; kap <- exp(p[4]); lt <- log(dd$t_months)
  Lam <- exp(drop(X %*% p[1:3]) + kap * lt)                              # cumulative hazard at t_i
  F2 <- rbind(p[3] * crossprod(X * Lam, Z) / n, p[3] * colSums(Z * (kap * lt * Lam)) / n)
  F2[3, ] <- F2[3, ] - colSums(Z * (dd$delta - Lam)) / n                 # Gamma_hat = n^-1 sum grad_gamma s_i
  wpp <- weibull_parts(p, dd$t_months, dd$delta, X); A2 <- wpp$H / n
  c(est = p[3], naive = sqrt(solve(A2)[3, 3] / n), naive_sand = sqrt(f$V_sand[3, 3]),   # naive: z'gamma_hat as data
    corrected = sqrt(two_step_vcov(wpp$S, A2, F2, R1)[3, 3]), gamma = gam)
}
r4 <- twostep(C4)
cat(sprintf("# gamma_hat %.4f %.4f %.4f; beta_u %.4f; SE naive %.5f, corrected %.5f (ratio %.3f)\n",
            r4[["gamma1"]], r4[["gamma2"]], r4[["gamma3"]], r4[["est"]], r4[["naive"]], r4[["corrected"]],
            r4[["corrected"]] / r4[["naive"]]))
# gamma_hat -1.0086 0.2534 0.3438; beta_u -0.8181; SE naive 0.04120, corrected 0.05006 (ratio 1.215)
one_rep4 <- function() {
  dd <- gen_C4(P4)
  r <- twostep(dd, start = c(P4$beta, log(P4$kappa)))
  o <- fit_weibull(dd$t_months, dd$delta, cbind(1, dd$fee_usd, dd$u_true), start = c(P4$beta, log(P4$kappa)),
                   engine = "newton")                                    # oracle: the true index
  c(r[c("est", "naive", "naive_sand", "corrected")], oracle = o$par[3], oracle_se = sqrt(o$V_hess[3, 3]))
}
set.seed(6124); sim4 <- t(replicate(P4$R_mc, one_rep4()))
cv <- function(e, s) mean(abs(sim4[, e] - P4$beta[3]) <= qnorm(0.975) * sim4[, s])
cat(sprintf("# 1000 reps: sd %.4f; mean SE naive %.4f, corrected %.4f; coverage naive %.3f, corrected %.3f, oracle %.3f\n",
            sd(sim4[, "est"]), mean(sim4[, "naive"]), mean(sim4[, "corrected"]),
            cv("est", "naive"), cv("est", "corrected"), cv("oracle", "oracle_se")))
# 1000 reps: sd 0.0481; mean SE naive 0.0405, corrected 0.0485; coverage naive 0.903, corrected 0.951, oracle 0.941
o4 <- fit_weibull(C4$t_months, C4$delta, cbind(1, C4$fee_usd, C4$u_true), engine = "newton")   # oracle, this sample
cat(sprintf("# naive sandwich: SE %.5f, coverage %.3f; oracle in the sample: %.4f (SE %.4f)\n", r4[["naive_sand"]],
            cv("est", "naive_sand"), o4$par[3], sqrt(o4$V_hess[3, 3])))
# naive sandwich: SE 0.04144, coverage 0.905; oracle in the sample: -0.8136 (SE 0.0410)
```

The naive sandwich (still treating the index as data) has a coverage rate of 0.905, almost the same as the naive 0.903: the problem lies not in the specification of the second-step likelihood but in ignoring the first step, and no choice of second-step variance can make up for it. The infeasible estimator with the true index has a standard error of 0.0410 and an interval of [-0.8940, -0.7332].

### 6.6 The Grand Reconciliation

| Block | Estimator | Estimate | 95% interval | Target and true value | Covered |
|---|---|---|---|---|---|
| C1 | OLS, censored as exits (log duration) | -0.01512 | [-0.01744, -0.01279] | -0.02222 | No |
| C1 | OLS, latent durations (log duration, infeasible) | -0.02319 | [-0.02581, -0.02057] | -0.02222 | Yes |
| C1 | Weibull MLE (log duration, `survreg`) | -0.02323 | [-0.02591, -0.02055] | -0.02222 | Yes |
| C1 | OLS, censored as exits (hazard scale) | 0.01533 | [0.01298, 0.01769] | 0.02 | No |
| C1 | OLS, canceled users only (hazard scale) | 0.01070 | [0.00775, 0.01365] | 0.02 | No |
| C1 | Weibull MLE, $\beta_{\mathrm{fee}}$ | 0.02056 | [0.01820, 0.02292] | 0.02 | Yes |
| C1 | Weibull MLE, $\kappa$ | 0.8848 | [0.8583, 0.9114] | 0.9 | Yes |
| C2 | QMLE $\kappa$, model standard error | 0.7275 | [0.7092, 0.7458] | Pseudo-true value 0.7586 | No |
| C2 | QMLE $\kappa$, OPG | 0.7275 | [0.7152, 0.7398] | Pseudo-true value 0.7586 | No |
| C2 | QMLE $\kappa$, sandwich | 0.7275 | [0.6990, 0.7560] | Pseudo-true value 0.7586 | No |
| C2 | QMLE $\beta_{\mathrm{fee}}$, sandwich | 0.01798 | [0.01640, 0.01956] | Pseudo-true value 0.01832 | Yes |
| C3 | Grouped-duration MLE (pooled cloglog), $\beta_{\mathrm{fee}}$ | 0.02031 | [0.01793, 0.02269] | 0.02 | Yes |
| C3 | Same, $\beta_{\mathrm{ref}}$ | -0.3978 | [-0.4830, -0.3126] | -0.40 | Yes |
| C4 | Two-step MLE, naive standard error | -0.8181 | [-0.8989, -0.7373] | -0.80 | Yes |
| C4 | Two-step MLE, corrected standard error | -0.8181 | [-0.9162, -0.7200] | -0.80 | Yes |

The table has two kinds of misses. The three C1 OLS misses are bias and have nothing to do with luck: censoring pushes them away from the truth, and another sample would do the same. The three C2 misses are, first of all, bad luck in this one sample: the shape-parameter estimate is lower by 2.14 sandwich standard errors, so even the sandwich interval misses; which of the three intervals to trust can only be read from the Monte Carlo of Section 6.3. C4's naive interval covers the truth in this sample, but the Monte Carlo says it covers in only 90.3% of replications; a success in a single sample likewise does not show that a variance formula is right.

The reconciliation of this section comes down to one sentence: write the observation mechanism into the likelihood, and the estimator aims at what it should aim at; for inference, either rely on the collapse that correct specification buys, or use the sandwich and the two-step correction honestly.

## 7 Failure Modes and Robustness

Dependent censoring is the most dangerous failure in duration analysis, because it cannot be seen in the data. CENS.1 requires that people who were cut short and people who were not would, with the same covariates, have stayed equally long. Administrative censoring such as an export date usually satisfies it; censoring of the kind "a user moves to a partner platform, the accounts are merged, and the record breaks off" does not, since those cut off may be the very users who were about to leave. There are three practical responses. First, put variables related to censoring (signup month, acquisition channel) among the covariates, so that CENS.1 need hold only given them. Second, separate administrative censoring from exits for other reasons, which often should be modeled as a different outcome (competing risks). Third, run a sensitivity analysis on the censoring mechanism, to see how the conclusions change under different assumptions, such as "censored users leave earlier or later." A closely related trap is a sample of only those still subscribed on some given day: such a sample systematically oversamples long durations, which Kiefer calls length-biased sampling, and the likelihood must be modified for truncation (Wooldridge 2010 §22.3.3).

A misspecified duration distribution is the second kind of risk. The Weibull can draw only a monotone hazard, and C2 shows that it then delivers a plausible-looking pseudo-true value, with the shape parameter and the coefficients shifting together. The response is to make the baseline hazard more flexible: grouped-duration cloglog has one free parameter per interval and keeps only the proportional hazard assumption (Theorem 9); in continuous time there is the proportional hazard model of Cox (1972), which uses the partial likelihood of Cox (1975) to estimate $\beta$ without estimating the baseline hazard. Checking whether the quantities of interest are stable across Weibull, log-logistic, and piecewise-constant specifications is itself a robustness check; non-nested models can be compared with the Vuong test, but it says only which is closer, not which is right.

Unobserved heterogeneity makes the shape of the hazard hard to interpret. A population hazard that falls with duration may mean that everyone becomes less and less inclined to leave, or merely that those inclined to leave left first. C2 has both: the deal-seekers' own hazard falls with duration, and the people remaining in the first few months are increasingly regular users, so the population hazard alone cannot separate the two parts. The model of Fader and Hardie (2007), by contrast, attributes the rise in aggregate retention entirely to the latter. Only under structural assumptions of the proportional hazard kind can the two be separated (Wooldridge 2010 p. 1006, p. 1010); Lancaster (1979) and Heckman and Singer (1984) discuss duration models with heterogeneity and the influence of their distributional assumptions. The discipline in practice: when reporting that "the hazard falls with duration," say clearly whether this describes the population or the behavior of individuals; the latter requires additional assumptions.

Reporting the wrong standard error is the most common problem. The model's built-in standard error is right only when the likelihood is correctly specified, and the direction in which it errs is not fixed: in C2 the model standard error of the shape parameter is too small, while that of the fee coefficient is, on the contrary, too large. The default choice should be the sandwich, which is asymptotically the same as the model standard error under correct specification and remains valid under misspecification. OPG is likewise right only under correct specification, and on C2 its coverage is only 0.5805; even under correct specification it often performs poorly in finite samples and should not be used on its own. The LR test relies on the information matrix equality; when the specification is in doubt, use the robust Wald or LM test instead. With person-period data, each row can be treated as an independent observation only when dynamic completeness holds; otherwise cluster by individual, all the more so when the same person has several spells.

In two-step estimation, the most common mistake is still to treat the first step as known. On C4 the coverage rate of the naive interval is 0.903, and switching to the second-step sandwich gives only 0.905, because the problem lies in the first step. The remedy is the correction of Theorem 6 of Chapter 5, or a bootstrap of the whole two-step procedure that re-estimates the first step in every resample (Section 4.4 of Chapter 8). The conclusion of Theorem 10 that "the first step actually improves efficiency" holds only when the second-step objective contains the first-step data and both conditional independence and condition (b) hold; it cannot be used to justify ignoring the first step.

Seen together, these failure modes show that the robustness of likelihood methods rests on two things. The first is the observation mechanism: the likelihood must be written to match the way censoring, grouping, and truncation arise, and writing it correctly rests on an independence assumption that cannot be tested. The second is the distributional specification: get it right and you obtain efficiency and collapsed standard errors; get it wrong and you obtain a pseudo-true value, in which case you must use the sandwich and ask whether this best approximation is the number you want to report.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Wooldridge (2010, *Econometric Analysis of Cross Section and Panel Data*, 2nd ed.), §13.1 to §13.11 and §22.1 to §22.4. The unified treatment of CMLE, partial MLE, two-step MLE, and QMLE, and of duration analysis; Theorems 1 and 3 of this chapter are its Theorems 13.1 and 13.2, Theorem 2 is §22.3.2, and Theorem 4 is the application of Theorem 14.3 in §14.4.2.
- White (1982, Econometrica). The original source of the QMLE, the pseudo-true value, and the sandwich, and also the origin of the information matrix test.
- Kiefer (1988, Journal of Economic Literature). The classic survey of duration analysis in economics.
- Fader and Hardie (2007, Journal of Interactive Marketing). A grouped-duration likelihood for retention in subscription businesses, and how heterogeneity masquerades as a rise in aggregate retention.

Further reading:

- Vuong (1989, Econometrica). The original paper on model selection tests for non-nested models.
- Engle (1984, *Handbook of Econometrics* Vol. 2). A survey of the asymptotic equivalence of the Wald, LR, and LM tests.
- Lancaster (1979, Econometrica) and Heckman and Singer (1984, Econometrica). Weibull duration models with covariates, and how distributional assumptions about unobserved heterogeneity affect the conclusions.
- Cox (1972, Journal of the Royal Statistical Society Series B) and Cox (1975, Biometrika). The proportional hazard model and the partial likelihood.
- Prentice and Gloeckler (1978, Biometrics). Regression for grouped survival data, and the source of the cloglog form.
- Tobin (1958, Econometrica). The original paper on the Tobit model.
- Schweidel, Fader and Bradlow (2008, Journal of Marketing) and Ascarza and Hardie (2013, Marketing Science). Two extensions of retention models in contractual settings: retention and duration dependence across cohorts, and a joint model of usage and renewal.
- Hansen (2014) Appendix B.11. A concise statement of the asymptotic theory of MLE and the Cramér-Rao lower bound.
:::

::: {.apa-refs}
- Ascarza, E., & Hardie, B. G. S. (2013). A joint model of usage and churn in contractual settings. *Marketing Science, 32*(4), 570-590. https://doi.org/10.1287/mksc.2013.0786
- Berry, M. J. A., & Linoff, G. S. (2004). *Data mining techniques: For marketing, sales, and customer relationship management* (2nd ed.). Wiley.
- Cox, D. R. (1972). Regression models and life-tables. *Journal of the Royal Statistical Society: Series B (Methodological), 34*(2), 187-202. https://doi.org/10.1111/j.2517-6161.1972.tb00899.x
- Cox, D. R. (1975). Partial likelihood. *Biometrika, 62*(2), 269-276. https://doi.org/10.1093/biomet/62.2.269
- Engle, R. F. (1984). Wald, likelihood ratio, and Lagrange multiplier tests in econometrics. In Z. Griliches & M. D. Intriligator (Eds.), *Handbook of econometrics* (Vol. 2, pp. 775-826). North-Holland. https://doi.org/10.1016/S1573-4412(84)02005-5
- Fader, P. S., & Hardie, B. G. S. (2007). How to project customer retention. *Journal of Interactive Marketing, 21*(1), 76-90. https://doi.org/10.1002/dir.20074
- Gourieroux, C., Monfort, A., & Trognon, A. (1984b). Pseudo maximum likelihood methods: Theory. *Econometrica, 52*(3), 681-700. https://doi.org/10.2307/1913471
- Hansen, B. E. (2014). *Econometrics* [Unpublished manuscript]. Department of Economics, University of Wisconsin.
- Heckman, J., & Singer, B. (1984). A method for minimizing the impact of distributional assumptions in econometric models for duration data. *Econometrica, 52*(2), 271-320. https://doi.org/10.2307/1911491
- Kiefer, N. M. (1988). Economic duration data and hazard functions. *Journal of Economic Literature, 26*(2), 646-679.
- Lancaster, T. (1979). Econometric methods for the duration of unemployment. *Econometrica, 47*(4), 939-956. https://doi.org/10.2307/1914140
- Prentice, R. L., & Gloeckler, L. A. (1978). Regression analysis of grouped survival data with application to breast cancer data. *Biometrics, 34*(1), 57-67. https://doi.org/10.2307/2529588
- Schweidel, D. A., Fader, P. S., & Bradlow, E. T. (2008). Understanding service retention within and across cohorts using limited information. *Journal of Marketing, 72*(1), 82-94. https://doi.org/10.1509/jmkg.72.1.082
- Tobin, J. (1958). Estimation of relationships for limited dependent variables. *Econometrica, 26*(1), 24-36. https://doi.org/10.2307/1907382
- Vuong, Q. H. (1989). Likelihood ratio tests for model selection and non-nested hypotheses. *Econometrica, 57*(2), 307-333. https://doi.org/10.2307/1912557
- White, H. (1982). Maximum likelihood estimation of misspecified models. *Econometrica, 50*(1), 1-25. https://doi.org/10.2307/1912526
- Wooldridge, J. M. (2010). *Econometric analysis of cross section and panel data* (2nd ed.). MIT Press.
:::
