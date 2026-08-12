---
title: "Sensitivity, Negative Controls & Partial Identification"
subtitle: "What Remains When Identification Assumptions Are Not Exact"
seriesline: "Foundations of Information Systems Economics · Chapter 15"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 15 · Sensitivity, Negative Controls & Partial Identification"
---

## Introduction

Nova's results table survives the familiar endurance test: swap the controls, drop the outliers, change the clustering level, add fixed effects, and the coefficient on the AI sales assistant stays around $0.696$, with a standard error of only $0.023$. Then the researcher moves an outcome determined before treatment to the left-hand side and still estimates an "effect" of $0.413$. The AI had not yet arrived, and the regression has already improved its past.

You cannot fix this by tacking on one more robustness column. If unobserved managerial ability raises both the probability of adoption and sales, the question you actually need to ask is: how strong would that confounding have to be to push the conclusion across zero, or across a threshold that matters for a decision? Nova's true effect is only $0.200$; the significance is rock solid, and that does nothing to stop the point estimate from being badly wrong.

This chapter therefore turns identification assumptions from a switch into a dial. Unconfoundedness, parallel trends, the exclusion restriction, and continuity at the cutoff are no longer merely declared "plausible"; we let them depart by an interpretable amount and watch how the estimate, the identified interval, and the policy choice move. Negative controls expose associations that should not be there, the robustness value measures how strong an omitted variable would have to be to overturn a conclusion, and partial identification honestly reports a set when the data are not enough to pin down a point. Sensitivity analysis is not an insurance policy stapled onto the original conclusion; it answers which specific failure mechanism, at what magnitude, would be enough to change the action we would take.

## 1 A Significant but Untrustworthy $0.696$

Nova has 8000 merchants. Whether a merchant adopts the AI assistant is not random: merchants with stronger managerial ability are more willing to adopt and also better at selling. The researcher observes baseline performance $X$ but cannot see management capability $U$.

Regressing log sales on adoption and $X$ gives $0.696$, with an SE of $0.023$. The conventional output would report a tiny p-value. If an oracle put $U$ into the regression as well, the coefficient would fall to $0.177$, close to the true effect of $0.200$. Significance does not protect the researcher from omitted-variable bias, because the larger the sample, the more precisely even the wrong conditional mean can be estimated.

Now look at an outcome fixed before treatment: pre-policy complaints. AI adoption cannot reach back and change past complaints, yet the same naive regression returns $0.413$. This is not another treatment effect; it is a diagnostic signal that the design has failed.

::: {.warning}
Sensitivity analysis does not prove that "the result must be right"; it makes public how strong an assumption the conclusion requires. Nor is a negative control an automatic corrector; it is first of all a way to probe residual bias using a relationship that should be zero in theory.
:::

## 2 From Point Identification to a Sensitivity Region

### 2.1 Three Different Kinds of Uncertainty

Sampling uncertainty comes from observing only a finite sample and is expressed through the SE and the confidence interval. Specification uncertainty comes from the many defensible analytic choices. Identification uncertainty comes from the fact that the core untestable assumptions may hold only approximately. This chapter focuses on the third kind, which does not vanish as $N$ grows.

### 2.2 Parameterizing Omitted-Variable Bias

In a linear model, the bias that an unobserved $U$ imparts to the treatment coefficient is jointly determined by two relationships: how much $U$ predicts $D$ given $X$, and how much $U$ predicts $Y$ given $(D,X)$. Cinelli and Hazlett use partial $R^2$ to put both on a comparable scale:

$$R^2_{D\sim U\mid X},\qquad R^2_{Y\sim U\mid D,X}.$$

If either one is close to zero, the bias is limited. Simply saying "there might be unobserved confounding" is too broad; a sensitivity contour asks where the two partial $R^2$ values would have to fall before the adjusted effect reaches zero.

### 2.3 Robustness Value

The robustness value is the smallest equal-strength partial $R^2$ needed to drive the estimate down to some benchmark. Nova's robustness value for the point estimate against zero is $0.287$. Intuitively, an unobserved confounder would have to explain about $28.7\%$ of the residual variance in both adoption and sales, after controlling for $X$, before it could push $0.696$ all the way to zero.

This number does not automatically qualify as "large" or "small." You should compare it against observed benchmark covariates and argue whether unobserved managerial ability could plausibly be stronger than the strongest observed predictor.

### 2.4 Bounds and the Identified Set

If the researcher is willing to assume only that the absolute bias is no larger than $B$, then the treatment effect is no longer a point but belongs to

$$\Theta(B)=[\widehat\tau-B,\widehat\tau+B].$$

Setting $B=0$ recovers point identification; as $B$ grows, the set widens. For Nova's set to include zero for the first time requires $B=0.696$. A wide set is not a failure of the method; it is the amount of information that the available data and assumptions can support.

### 2.5 Breakdown Point and Breakdown Frontier

The breakdown point is the smallest departure from an assumption at which a given conclusion just stops holding. When several assumptions loosen at once, the breakdown frontier traces out the boundary combinations. For example, if we allow both selection on unobservables and treatment effect heterogeneity, which combinations still guarantee that the effect is positive? This representation is closer to the substantive uncertainty than listing dozens of specification coefficients.

To summarize this section: sampling, specification, and identification uncertainty are distinct; partial $R^2$ turns omitted confounding into an interpretable strength; and bounds and the breakdown frontier report what we can still learn once assumptions are loosened.

## 3 Identification: Diagnostics, Calibration, and Sets

### 3.1 Negative-Control Outcome

A negative-control outcome $Y_i^{NC}$ is one that theory says should not be affected by treatment, yet that shares part of the confounding structure with the target outcome. If

$$Y_i^{NC}(1)=Y_i^{NC}(0),$$

but $D$ still significantly predicts $Y^{NC}$ after adjustment, then either exchangeability or the measurement process is off in at least one place. Nova's pre-policy complaints are fixed before adoption, so the coefficient of $0.413$ cannot be interpreted as a causal effect.

::: {.assumption data-label="Negative-control validity"}
The negative-control outcome is not affected by treatment or its versions, but it shares with the target outcome the source of unobserved bias the researcher wants to probe; its measurement error and sample selection must also not mechanically manufacture a treatment association.
:::

The phrase "shared source of bias" is crucial. An old variable with nothing to do with managerial ability, even with a coefficient of zero, cannot rule out management confounding. A negative control must be chosen from domain theory, not picked at random from the database as some pre-treatment outcome.

### 3.2 Negative-Control Exposure

A negative-control exposure is one that theory says does not affect the target outcome, but that shares the confounding with the real exposure. A future treatment, a wrong time window, or a same-family policy known to have no effect can sometimes serve. If it predicts the current outcome, that points to a residual time trend, selection, or measurement leakage.

### 3.3 Calibration by Observed Covariates

Comparing the sensitivity parameter against the strength of an observed covariate is more meaningful than debating $R^2=0.1$ out of thin air. You can report that "an unobserved confounder would need to be twice as strong as baseline sales to overturn the conclusion." But the observed benchmark is not an upper bound, especially when what is truly omitted is management quality, a demand shock, or an algorithmic ranking.

### 3.4 Rosenbaum Bounds

In matched observational studies, Rosenbaum uses $\Gamma$ to denote how much the treatment odds of two observably identical units can differ at most because of unobserved factors. $\Gamma=1$ corresponds to no hidden bias; as $\Gamma$ grows, the randomization-based p-value bounds widen. This directly answers "how much hidden selection would make the significance disappear," but it mainly targets sharp-null testing and is not equivalent to effect-size bounds.

### 3.5 Lee Bounds and Selective Retention

If treatment affects sample selection monotonically, for example if it can only raise and never lower survey response, Lee bounds construct effect bounds by trimming the outcome distribution of the treatment group. Monotonicity is not a free assumption; when a platform treatment can make some users stay while others leave, simple trimming no longer carries a guarantee.

### 3.6 Honest DiD Sensitivity

The pre-trends of Chapter 13 cannot prove parallel trends in the post-period. Rambachan and Roth impose relative-magnitude or smoothness restrictions on the violation of the untreated trend to construct robust confidence sets. You should report how the conclusion changes when the violation is one times, two times, and so on the size of the observed pre-trend, along with the breakdown value at which zero first enters.

### 3.7 Point Identification and Set Identification

The reduced-form partial identification of this chapter differs from the structural partial identification of the game-theoretic part in Chapter 22. Here we start from loosening unconfoundedness, parallel trends, or selection assumptions; the sets in Chapter 22 come from multiple equilibria and an unspecified equilibrium selection. The common principle is: if the data and credible assumptions support only a set, you should not force out a spuriously precise point.

### 3.8 The Sensitivity Parameter Must Correspond to a Specific Failure Mechanism

A generic $B$ can show how the identification region expands, but a substantive analysis must also connect $B$ back to an institutional mechanism. For AI adoption, hidden management quality affects the treatment odds and the outcome residual; for DiD, the departure comes from curvature in the untreated trend; for IV, the departure comes from the instrument's direct effect on the outcome; for attrition, the departure comes from treatment-induced selection. The four cannot share a single phrase like "allow $10\%$ bias."

A better way to report is to write down the causal path that violates the assumption first, then choose the sensitivity parameter that corresponds to that path, and calibrate it in observed benchmarks or institutionally interpretable units. If the parameter has no real-world meaning, the reader cannot judge whether the breakdown value is plausible.

### 3.9 From Statistical Breakdown to Decision Robustness

Whether zero enters the interval is not the only decision threshold. A platform might launch only when the effect exceeds the implementation cost, and a regulator might care whether harm exceeds a tolerated ceiling. Let $c$ be the policy threshold; you should report how far assumptions have to loosen before the lower bound falls below $c$, rather than only asking when the p-value exceeds $0.05$.

This distinction turns sensitivity analysis from "protecting significance" into decision analysis. A small effect that is statistically always positive but drops below the economically meaningful MDE under a tiny departure still leaves the policy conclusion fragile; conversely, an effect whose confidence set includes zero but whose plausible sensitivity region mostly supports enormous gains should not be dismissed as simply "null."

## 4 Estimation: Turning the Sensitivity Parameter into Results

### 4.1 OVB Contour

For each pair $(R^2_{D\sim U\mid X},R^2_{Y\sim U\mid D,X})$, compute the adjusted estimate and draw a contour rather than reporting a single robustness value. The plot can show how the conclusion changes as the confounder predicts treatment or outcome more strongly.

### 4.2 Benchmark Table

Remove observed covariates one at a time, compute each one's partial $R^2$ and induced bias, and express the latent confounder's strength as "if as strong as baseline performance" or "twice as strong." A benchmark must have substantive relevance in advance; you cannot just pick the weak variable that makes the result look robust.

### 4.3 Identified-Region Curve

Increase the allowed absolute bias $B$ continuously from zero and report the lower and upper bounds. Rather than a single threshold, the whole curve lets the reader use their own domain judgment to choose a reasonable $B$.

### 4.4 Negative-Control Estimate

The negative-control coefficient should be paired with the same specification, sample, and inference. If it is clearly nonzero, the main result cannot just note in an appendix that the placebo failed; it should explain which kind of bias the coefficient corresponds to and use its strength for sensitivity calibration.

### 4.5 Inference over Bounds

Estimated bounds also carry sampling uncertainty. You cannot treat the sample lower bound as a known constant. In practice you need confidence regions, bootstrap, or moment-inequality inference suited to partially identified parameters, and you must distinguish the identification region from the confidence region.

### 4.6 Sensitivity Analysis Must Match the Research Design

There is no universal "robustness check" for every identification threat. Different designs rest on different core counterfactuals, so the sensitivity parameter should differ too.

| Research design | Core threat | Appropriate sensitivity target | Question it cannot answer |
|---|---|---|---|
| Selection on observables | Unobserved confounding | partial $R^2$, Rosenbaum $\Gamma$ | treatment versions |
| IV | exclusion, monotonicity | direct-effect bounds, defier share | outcome mismeasurement |
| DiD | differential trend | relative-magnitude or smoothness bound | spillover contamination |
| RDD | manipulation, functional form | donut, bandwidth, bounded sorting | extrapolation beyond the cutoff |
| Attrition | selective observation | Lee bounds, worst-case bounds | construct validity |
| Marketplace experiment | interference | exposure mapping, saturation curve | long-run equilibrium stability |

Lining up many unrelated specifications in a robustness table is not the same as doing sensitivity analysis on the key assumption. A coefficient that stays significant under twenty sets of controls still cannot answer how strong an unobserved confounder would suffice to overturn the conclusion; likewise, changing the bandwidth cannot diagnose the parallel trends of a DiD.

### 4.7 Benchmarking Cannot Mechanically Rely on the "Strongest Covariate"

Observed-covariate benchmarking compares the strength of unobserved confounding against observed variables, and this is the key step for turning an abstract sensitivity parameter onto a substantive scale. But "an unobserved variable cannot be stronger than the strongest observed variable" is not a theorem. Observed covariates are often finely measured, while the management quality, demand shock, or algorithmic targeting score that is truly omitted may strongly predict both treatment and outcome.

A better report gives several benchmarks: the single strongest covariate, a group of theoretically relevant covariates, and the change in the treatment coefficient after removing a particular group of variables. If the conclusion holds only when unobserved confounding is weaker than a very weak demographic variable, then the robustness value, even if nonzero, lacks persuasive force. Conversely, if you would need a confounder that predicts both treatment and outcome more strongly than the complete pre-treatment performance history, the conclusion is more credible.

### 4.8 From Effect Robustness to Policy-Loss Robustness

The decision significance of a research conclusion depends on a threshold, not only on whether the effect crosses zero. Let the net benefit of implementing the policy be

$$\Pi(\tau)=M\tau-C,$$

where $M$ is the affected scale and $C$ is the deployment and opportunity cost. If sensitivity analysis gives $\tau\in[\tau_L(\eta),\tau_U(\eta)]$, then robust decisions come in three kinds: when $M\tau_L-C>0$, implementing is still dominant under the given departure; when $M\tau_U-C<0$, not implementing is still dominant; and the remaining region shows that identification uncertainty is enough to change the decision.

This decision-focused sensitivity changes the emphasis of the report. A small effect that is statistically "still positive" may not be enough to cover the deployment cost; an effect whose interval straddles zero may still be worth piloting when the deployment cost is very low and the downside is bounded. A high-standard empirical report should present both the statistical breakdown frontier and the policy frontier, and state transparently where the cost, scale, and loss function come from.

## 5 Anchor Papers

### 5.1 Cinelli and Hazlett (2020)

::: {.case}
Paper: "Making Sense of Sensitivity: Extending Omitted Variable Bias," Journal of the Royal Statistical Society: Series B.

Method: rewrites omitted-variable bias in linear regression using partial $R^2$, and proposes the robustness value and covariate benchmarking.

Data: the method is illustrated through several empirical examples showing how to convert abstract unobserved confounding into a strength comparable with observed covariates.

Result: the researcher can compute the minimum confounding strength needed to overturn the conclusion without ever observing $U$.

Limitation: the results are tied to a specific linear projection and sensitivity parameterization; they do not prove that a given strength is impossible in a specific institution.
:::

### 5.2 Lipsitch, Tchetgen Tchetgen, and Cohen (2010)

::: {.case}
Paper: "Negative Controls: A Tool for Detecting Confounding and Bias in Observational Studies," Epidemiology.

Method: systematically distinguishes negative-control exposure from negative-control outcome, using relationships that should be zero in theory to probe confounding, selection, and measurement bias.

Data: the paper uses epidemiological designs to illustrate the choice and interpretation of negative controls.

Result: when a valid negative control shows an association, it can veto some of the unbiasedness assumptions that support the main analysis.

Limitation: a zero association cannot prove unbiasedness, especially when the negative control does not share the confounding structure of the main result or has too little power.
:::

### 5.3 Rambachan and Roth (2023)

::: {.case}
Paper: "A More Credible Approach to Parallel Trends," Review of Economic Studies.

Method: allows post-treatment violations to be bounded or smooth relative to the pre-trends, and constructs honest confidence sets and breakdown values.

Data: the paper revisits several classic DiD applications and shows how the conclusions are sensitive to different deviations.

Result: an insignificant pre-trend does not amount to robustness; some conclusions fail under a very small deviation, while others hold over a wide range.

Limitation: relative-magnitude and smoothness bounds still require domain judgment; the method makes that judgment transparent but does not choose it for the researcher.
:::

## 6 A Full Walkthrough on the Nova Data

### 6.1 DGP and the Naive Result

```r
adopt <- rbinom(N,1,plogis(-0.2 + 0.7*x + 0.9*u_management))
log_sales <- 2 + 0.20*adopt + 0.45*x + 0.65*u_management + rnorm(N,sd=0.75)
pre_policy_complaints <- 0.55*u_management + 0.20*x + rnorm(N,sd=0.80)
```

| Estimator or diagnostic | Estimate | SE |
|---|---:|---:|
| Naive adjusted adoption | 0.696 | 0.023 |
| Oracle including management | 0.177 | 0.019 |
| Negative-control outcome | 0.413 | 0.022 |

### 6.2 Robustness Contour

Nova's robustness value is $0.287$. In the figure, the horizontal axis is the partial $R^2$ of $U$ for adoption, the vertical axis is the partial $R^2$ of $U$ for sales, and the color gives the adjusted effect.

![Omitted-variable sensitivity contour. As the unobserved factor predicts both adoption and sales more strongly, the adjusted effect moves from the naive 0.696 toward zero; the plot makes both dimensions of the identification threat public rather than merely swapping controls.](assets/fig/fig_15_robustness_contour.svg)

### 6.3 Identified Region

```r
bias_cap <- seq(0,0.45,by=0.01)
region <- data.frame(lower=b-bias_cap, upper=b+bias_cap)
```

![The identified region as the allowed absolute omitted-variable bias grows. The weaker the assumption, the wider the set; over the range shown (bias allowed up to 0.45) the set never contains zero, and for it to include zero for the first time you must allow the bias to reach the naive coefficient of 0.696.](assets/fig/fig_15_identification_region.svg)

### 6.4 The Full Reconciliation

The large sample estimates the naive association of $0.696$ very precisely, yet it lies far from the true effect of $0.200$. The oracle coefficient of $0.177$ shows that the bias is indeed produced by latent management. The negative control's $0.413$ shows that the same selection process has already contaminated the pre-treatment outcome. Sensitivity analysis will not conjure the oracle back out of thin air, but it tells the reader clearly: what confounding strength it would take to push the result to zero, and which effect values can still be ruled out under weaker bias.

## 7 Failure Modes and Robustness

First, mistaking coefficient stability for causal robustness. Stability after adding more observed controls only shows that those controls did not change the linear projection; it places no restriction on unobserved factors.

Second, treating the robustness value as a quality score. A $0.287$ has no universal "large" or "small"; it must be compared against the plausible confounders in the field.

Third, picking irrelevant negative controls. A placebo that passes easily but does not share the bias structure offers very limited comfort.

Fourth, ignoring a negative-control failure. Continuing to treat the main estimate as unbiased after a placebo turns significant, mentioning it only in a footnote, is inconsistent evidentiary logic.

Fifth, confusing the identification region with the confidence interval. The former expresses assumption uncertainty, the latter sampling uncertainty, and both need to be reported.

Sixth, doing only single-axis sensitivity. In reality several assumptions loosen at once, and the breakdown frontier is more honest than one-at-a-time checks.

Finally, sensitivity analysis cannot substitute for a better design. If you can obtain a random rollout, a valid IV, or a more credible cutoff, improve the assignment mechanism first; sensitivity quantifies the remaining fragility, it does not provide decoration for a weak design.

## 8 Further Reading

::: {.readings}
Required reading:

- Cinelli and Hazlett (2020). Linear OVB, partial $R^2$, and the robustness value.
- Lipsitch, Tchetgen Tchetgen, and Cohen (2010). The basic logic of the negative-control design.
- Rambachan and Roth (2023). Parallel-trends sensitivity and honest confidence sets.

Further reading:

- Rosenbaum (2002). Hidden-bias bounds for matched observational studies.
- Masten and Poirier (2020). Breakdown frontiers and the joint loosening of multiple assumptions.
- Lee (2009). Trimming bounds under monotone sample selection.
:::

::: {.apa-refs}
- Cinelli, C., & Hazlett, C. (2020). Making sense of sensitivity: Extending omitted variable bias. *Journal of the Royal Statistical Society: Series B, 82*(1), 39-67. https://doi.org/10.1111/rssb.12348
- Lee, D. S. (2009). Training, wages, and sample selection: Estimating sharp bounds on treatment effects. *The Review of Economic Studies, 76*(3), 1071-1102. https://doi.org/10.1111/j.1467-937X.2009.00536.x
- Lipsitch, M., Tchetgen Tchetgen, E., & Cohen, T. (2010). Negative controls: A tool for detecting confounding and bias in observational studies. *Epidemiology, 21*(3), 383-388. https://doi.org/10.1097/EDE.0b013e3181d61eeb
- Masten, M. A., & Poirier, A. (2020). Inference on breakdown frontiers. *Quantitative Economics, 11*(1), 41-111. https://doi.org/10.3982/QE1288
- Rambachan, A., & Roth, J. (2023). A more credible approach to parallel trends. *The Review of Economic Studies, 90*(5), 2555-2591. https://doi.org/10.1093/restud/rdad018
- Rosenbaum, P. R. (2002). *Observational studies* (2nd ed.). Springer.
:::
