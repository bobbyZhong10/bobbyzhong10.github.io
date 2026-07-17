---
title: "From an Economic Question to a Research Design"
subtitle: "Units, Treatments, Estimands, and Assignment Mechanisms"
seriesline: "Foundations of Information Systems Economics · Chapter 10"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 10 · From an Economic Question to a Research Design"
---

## Introduction

Aurora randomly opens a recommendation tool early to a batch of merchants. The simplest experimental comparison shows sales rising. The researcher, uneasy, then adds post-adoption engagement to the regression to "control more thoroughly," and the coefficient flips to negative. Randomization did not break, and the software did not miscompute. What broke was the question: engagement is a variable realized after treatment, and once you condition on it, the regression is no longer comparing the two worlds you started with.

Accidents like this rarely surface in the last row of a results table. They are usually planted long before the first line of code. What does "using AI" even mean: gaining access, clicking once, sustained adoption, or being genuinely influenced by the algorithm? Is the outcome sales in the current week, long-run merchant retention, or total platform volume? Is the unit of study the merchant, the consumer, the order, or the market? Once these choices are left vague, a smaller standard error may mean only that the researcher is answering, ever more precisely, a question no one actually asked.

For this reason, this chapter does not rush to pick a regression. It first breaks the research question into unit, treatment, timing, exposure, outcome, estimand, and assignment mechanism. Potential outcomes are there to say which two counterfactuals you mean to compare; a DAG is there to record which paths may be closed and which should not. Neither identifies a causal effect for you automatically. The real test is always the same: does the estimator in front of you still answer the decision problem you started with? Why Aurora's positive effect gets rewritten into a negative one by "more thorough controls" will run through this design language from start to finish.

## 1 An Experimental Effect That Turns from Positive to Negative

The Aurora platform opens a new recommendation tool early to a random half of its merchants. Before launch, the pre-analysis plan states the question clearly: if the tool were opened to all merchants, how would sales change in the next period? In the data, the early-access group and the control group differ in mean log sales by $0.149$, with a standard error of $0.017$. Because assignment is random, this difference targets the intention-to-treat effect directly, and it is close to the known total effect $0.140$ in the DGP.

The analyst then adds baseline quality and region fixed effects, and the estimate becomes $0.157$, with the standard error shrinking to $0.015$. This is still reasonable: both variables are determined before treatment, and adjusting for them absorbs chance imbalance and outcome noise without changing what the treatment means.

Next, the analyst controls for post-launch engagement, on the grounds that "active merchants sell more anyway." The treatment coefficient immediately turns to $-0.099$, and it is highly significant. Does the new tool help sales or hurt them? The randomized experiment has not changed. What changed is the variable the researcher conditioned on and the estimand that goes with it.

::: {.warning}
Adding a control is not an operation that monotonically increases credibility. Whether a variable is realized before or after treatment, whether it is a common cause or a mediator, and whether it is jointly determined by unobservables, matters more than its correlation with the outcome. The regression software will not warn you that you have swapped the question.
:::

## 2 From the Question to the Estimand

### 2.1 First Write Down Eight Design Objects

Before writing a regression equation, fill in the following eight items.

| Design object | Definition in Aurora | What happens if the definition is vague |
|---|---|---|
| Unit | Merchant | User, merchant, and market effects get mixed together |
| Treatment | Early access to the new tool | Confounded with actual use or recommendation exposure |
| Version | Same interface, same ranking model | Different model versions treated as one treatment |
| Timing | Assignment in period 1, measurement in period 2 | Anticipation and lagged effects cannot be separated |
| Exposure | Whether the merchant received access | Clicking into the tool mistaken for assignment |
| Outcome | Next-period log sales | Short-run GMV conflated with long-run retention |
| Population | Eligible merchants in the pilot | No basis for extrapolating to non-participating merchants |
| Assignment | Random early opening within region | Cannot say where the counterfactual comes from |

This table is not administrative paperwork. Change any one cell, and the potential outcomes may change. For example, $Y_i(d)$ can denote merchant $i$'s sales under whether it received access, or its sales under whether it actually used the tool. The former identifies ITT through random assignment; the latter, if noncompliance is present, requires the IV/LATE logic of Chapter 12. The two cannot be collapsed into one estimand just because both use a single binary variable.

### 2.2 Total Effect and Direct Effect

Let $D_i$ denote early access, $M_i(d)$ engagement under treatment status $d$, and $Y_i(d,m)$ the outcome when access and engagement are both set. The platform's original question corresponds to the total effect:

$$\tau^{total} = \mathbb{E}\big[Y_i(1,M_i(1)) - Y_i(0,M_i(0))\big].$$

It lets the tool work through engagement. If engagement is held at a single level, the question becomes the controlled direct effect:

$$\tau^{direct}(m) = \mathbb{E}\big[Y_i(1,m) - Y_i(0,m)\big].$$

Both can be meaningful scientific questions, but they are not the same question. Aurora's structural direct effect is $0.080$, the tool raises engagement by $0.50$ on average, and the marginal effect of engagement on sales is $0.12$, so the total effect is $0.08 + 0.50\times0.12 = 0.140$. Controlling for the mediator is not "estimating the total effect more precisely"; it is an attempt to cut one of the channels through which the effect operates.

::: {.intuition}
First ask through which channels the policy will work, then decide whether to control for the channel variables. If the platform cares about the total consequences of a full rollout, the engagement change induced by treatment is part of the effect, not noise to be scrubbed out.
:::

### 2.3 Estimand, Estimator, and Estimate

Three words must be kept apart. The estimand is the target, for instance $\tau^{total}$; the estimator is the rule mapping a sample to a number, for instance the treated-control mean difference; the estimate is the $0.149$ this sample hands over. Swapping a regression specification may only swap the estimator, or it may quietly swap the estimand. The first discipline of research design is to ask, each time you change the specification, whether it is still aiming at the same object.

The summary of this section is as follows: an empirical question must first fix the unit, treatment version, timing, exposure, outcome, population, and assignment; total and direct effects are different targets; and estimand, estimator, and estimate are not interchangeable.

## 3 Identification: Where Does the Counterfactual Actually Come From

### 3.1 The Assignment Mechanism Is the Engine of the Design

Causal identification does not begin with an outcome regression, it begins with how the treatment was assigned. If $D_i$ is randomly assigned within strata $B_i$ (a stratified-block indicator), then

$$\{Y_i(1),Y_i(0)\}\perp D_i\mid B_i.$$

Aurora opens access early at random within region, so a region-stratified comparison is backed by experimental design. Adding baseline quality can improve precision, but identification does not depend on getting a sales regression right. If early access were instead picked by account managers choosing "the most promising merchants," the assignment mechanism would be different, and observed confounding, IV, RDD, DiD, or some other design would be needed to supply the missing counterfactual.

::: {.assumption data-label="Design-based exchangeability"}
Within the actual randomization level and strata, assignment is independent of all potential outcomes. The analysis must preserve the level, probabilities, and cluster structure used in the randomization.
:::

This assumption constrains treatment assignment, not that the outcome must be linear. It is backed by the implementation protocol and can be checked against balance and an assignment audit to see whether execution drifted, but it cannot be proven from the goodness of fit of an outcome regression.

### 3.2 What Potential Outcomes and DAGs Each Do

Potential outcomes are good at fixing the intervention and the estimand; DAGs are good at displaying the conditional-independence structure among variables. The two are not competing languages. For Aurora, you can draw $D\rightarrow M\rightarrow Y$, $D\rightarrow Y$, $U\rightarrow M$, and $U\rightarrow Y$. This graph reminds us at once that $M$ comes after treatment and is a joint effect of $D$ and $U$.

![Aurora's assignment graph. Early access affects sales directly and also through post-treatment engagement; latent service capability affects both engagement and sales. Controlling for engagement cuts the mediating path and opens a conditional association between early access and latent capability.](assets/fig/fig_10_assignment_dag.svg)

A DAG cannot substitute for economic content. The arrows come from claims about platform institutions, user behavior, and temporal order; they are not facts learned automatically from a correlation matrix. Different plausible DAGs may imply different adjustment sets, and the honest practice is to show the key alternative graphs and state which arrows the conclusion is sensitive to.

### 3.3 Confounder, Mediator, and Collider

A baseline confounder jointly affects $D$ and $Y$ before treatment, and controlling for it appropriately can recover conditional exchangeability. A mediator sits on the path from $D$ to $Y$, and controlling for it takes away part of the total effect. A collider is a variable that two arrows both point into, and conditioning on it manufactures an association between its causes.

In Aurora, engagement is affected by both $D$ and latent service capability $U$, so on the path $D\rightarrow M\leftarrow U\rightarrow Y$ it is a collider. Randomization originally guarantees that $D$ and $U$ are unrelated; conditional on $M$, treated merchants with the same engagement tend to require a lower $U$, so $D$ and $U$ become negatively correlated in the conditioned sample. The result is that the treatment coefficient not only loses the mediated effect but also absorbs collider bias, and it flips from positive to $-0.099$.

::: {.theorem data-label="Bad-control principle"}
If a variable is affected by the treatment, adding it to an outcome regression generally no longer identifies the original total effect; and if it is also a collider between the treatment and an unobserved cause of the outcome, the adjusted treatment coefficient usually does not identify the controlled direct effect either.
:::

### 3.4 Pre-treatment Does Not Mean It Should Be Controlled

"Control only for pre-treatment variables" is a necessary safety screen, but not a sufficient rule. A pre-treatment collider still should not be controlled; a strong predictor close to the treatment that is unrelated to the outcome may lower the precision of some estimators; and in an IV design, the consequences of controlling for the instrument are different again. The reason for choosing controls should come from the assignment and the causal structure, not from stepwise regression, significance, or "throw in everything you can find."

### 3.5 Mechanism Evidence Does Not Mean Mechanism Identification

If the treatment raises engagement and engagement is positively correlated with sales, this is consistent with an engagement channel, but it does not identify the mediated effect. Formal mediation must handle mediator-outcome confounding, post-treatment confounding, and a well-defined mediator intervention. A mechanism prediction being supported by the data says only that the model has not yet been rejected by that prediction; it is not the unique identification of the mechanism against alternative explanations.

### 3.6 Timing, Anticipation, and Treatment Versions

Platform policies often change behavior after the announcement but before the formal launch. If merchants prepare in advance, $Y_i(0)$ is already affected by future assignment in the so-called pre-period, and no anticipation fails. Algorithms also keep iterating, and an "AI assistant" of the same name may correspond to different model versions in different weeks. The researcher should record versions, deployment times, and exposure logs, and make clear whether the estimand is the short-run effect of one version or the policy effect of a family of versions.

### 3.7 Testability and the Design Audit

Randomization probabilities, stratification rules, sample exclusions, treatment crossover, and logging completeness can be audited; the independence of unobserved potential outcomes cannot be fully verified from outcome data. Balance checks are an execution diagnostic, not an item-by-item significance exam. More valuable is reconstructing the assignment code, reconciling exposure timestamps, checking whether attrition is affected by treatment, and locking in the primary estimand while still blind to the outcome.

### 3.8 Attrition, Truncation, and "Who Is Still in the Sample"

Platform outcomes often become unobservable only after treatment. A merchant that exits has no sales, a user that uninstalls has no engagement, and a refused impression has no click. If we define $R_i(d)$ as whether the outcome is recorded under treatment status $d$, what we actually observe is $R_i(D_i)Y_i(D_i)$. Even if $D$ is random, a complete-case comparison compares only between treatment-induced selected populations.

First distinguish three objects. A missing outcome is a unit still in the target population with $Y$ unrecorded; truncation by death is treatment changing whether the unit exists or has a well-defined outcome at all; a sample restriction is the researcher deleting the sample after the fact according to observed behavior. IPW or imputation requires the corresponding missing-at-random condition and cannot automatically handle the case where $R$ depends on the unobserved $Y$. The Lee bounds of Chapter 15 provide weaker conclusions under monotone selection.

### 3.9 Internal Validity, External Validity, and Policy Transport

An experiment's effect on the pilot population can be internally valid without representing a national rollout. Let $S_i=1$ denote entering the study sample; the target-population effect requires transporting the response of $S=1$ to $S=0$. This requires all effect modifiers to be observed and to have overlap, or bounds on the unobserved heterogeneity.

Platform studies especially need to record the eligibility rule, invitation, consent, actual exposure, and analysis inclusion. The study sample is often active users, trackable devices, or firms willing to try AI, precisely the population least suited to unconditional extrapolation. A high-standard conclusion should be written as "valid for which units, which version, which horizon, and which rollout saturation," rather than merely stating in a limitations paragraph that external validity is limited.

The key points of this section can be summarized as follows: identification comes from the assignment mechanism; potential outcomes fix the question and DAGs display the conditional structure; post-treatment controls can simultaneously change the estimand and open a collider path; timing, versions, attrition, and implementation drift all need to enter the design audit; and transport from a pilot sample to a policy population is another set of assumptions that cannot be inferred automatically from internal validity.

## 4 Estimation: Choose the Estimator Only After the Design Is Fixed

### 4.1 Difference in Means

In an executed simple randomized experiment, the most transparent estimator is

$$\hat\tau = \bar Y_1-\bar Y_0.$$

It maps directly to the assignment. Aurora's estimate is $0.149$. If randomization is stratified within region, the analysis should weight by the known assignment probability or add strata fixed effects, and compute standard errors at the actual randomization level.

### 4.2 Baseline Adjustment

Covariates measured in advance that strongly predict the outcome can improve precision. After adding baseline quality and region fixed effects, Aurora obtains $0.157$, with the SE falling from $0.017$ to $0.015$. Randomization guarantees the unadjusted estimator is unbiased; the role of adjustment is mainly precision and finite-sample balance, and it should not be described as "removing endogeneity."

### 4.3 ITT, TOT, and the Exposure Effect

If having access is not the same as using it, then assignment $Z$, actual adoption $D$, and algorithmic exposure $E$ must be kept apart. The ITT of $Z$ on $Y$ is identified by randomization; using $Z$ as an instrument for $D$ identifies the complier LATE under exclusion and monotonicity; and directly comparing $D=1$ with $D=0$ often reintroduces selection. Three binary columns in a data table do not mean three causal effects are all identified by the same experiment.

### 4.4 A Specification Curve Cannot Substitute for Design

Showing several plausible specifications helps clarify whether the limited modeling choices matter, but running hundreds of control combinations does not create exogenous variation. If all specifications rely on the same untrustworthy assignment assumption, their joint stability does not amount to identification. Chapter 15 will parameterize deviations from the core assumption directly, rather than treating "swap a few controls" as a complete sensitivity analysis.

### 4.5 Deriving Distinguishable Empirical Predictions from Economic Mechanisms

Research design cannot stop at "does the treatment work." For IS and OM questions, a more valuable approach is to first write out competing mechanisms, then look for the margin that separates them. Suppose an AI copilot raises customer-service performance; there are at least three explanations: it lowers the handling cost of a single task; it supplies knowledge to low-experience workers; and it changes task allocation, so that easier tasks flow more toward treated workers. All three can raise the average handling volume, yet they give different predictions about quality, heterogeneity, and task composition.

| Mechanism | First-order prediction | Distinguishing additional prediction | Main threat |
|---|---|---|---|
| Cost reduction | Volume completed per unit time rises | Roughly proportional improvement across experience groups | Contemporaneous demand shifts |
| Knowledge supplement | Larger improvement for low-experience workers | Escalation falls, learning curve flattens | Regression to the mean |
| Task reallocation | Treated throughput rises | Untreated task difficulty rises | Interference and selective logging |

The discipline here is that heterogeneity constitutes mechanism evidence only when it is defined before treatment, derived from the mechanism, and paired with the corresponding measurement. Searching for significant results across dozens of subgroups after the fact can only generate hypotheses. Mechanism tests should also avoid controlling for post-treatment variables; if the task mix is itself affected by treatment, putting it into the outcome regression rewrites the total effect into a controlled direct effect that is hard to interpret.

### 4.6 Temporal Structure, Dynamic Paths, and Cumulative Effects

Technology adoption often contains installation, training, active use, and workflow redesign at once. Compressing these dates into a single treatment date conflates the different stages. The research design should distinguish at least three estimands: the short-run disruption right after deployment, the flow effect once use stabilizes, and the cumulative effect that includes the up-front investment cost. If $\tau_k$ denotes the effect in period $k$ after adoption, the cumulative policy value is not whichever $\tau_k$ is most significant, but

$$V_H=\sum_{k=0}^{H}\delta^k\tau_k-C_0,$$

where $C_0$ includes migration, training, and organizational-adjustment costs. The empirical meaning of the productivity J-curve is exactly this: early $\tau_k$ can be negative and later ones positive, and whether $V_H$ is positive depends on the evaluation horizon. The researcher should state in advance whether event time is defined by access, first use, or stable use, and report the adoption rate alongside; otherwise the estimated rollout effect will blend technological capability with adoption behavior.

### 4.7 From Identification Assumptions to an Auditable Design Table

A high-quality paper should let the reader see what fact supports each key assumption, rather than only listing terminology in the methods section. Before estimation, you can build a design table like the following.

| Design object | Question that must be answered | Observable evidence | Part still not testable |
|---|---|---|---|
| Assignment | Who decides when treatment is received | Rollout rule, randomization code, eligibility threshold | Unrecorded manual exceptions |
| Exposure | Does access equal use | Login, feature call, duration | Quality of use and cognitive investment |
| Outcome | Does the metric correspond to the construct | Audit, reconciliation, alternative metric | Unrecorded behavior and long-run welfare |
| Counterfactual | Why the control can represent the untreated state | Balance, pre-trends, placebo | Future untreated potential outcome |
| Interference | Whether units affect one another | Network, boundary flows, saturation response | Unobserved cross-market equilibrium responses |
| Transport | Whether the sample generalizes to the policy population | Overlap, site composition, reweighting | Structural change in a new environment |

The point of this table is not to claim that all assumptions can be tested, but to separate design facts, diagnostic evidence, and untestable assumptions. Rigor comes from exposing the untestable part and stating within what range of deviation the conclusion still holds.

## 5 Anchor Papers

### 5.1 Rubin (1974)

::: {.case}
Paper: "Estimating Causal Effects of Treatments in Randomized and Nonrandomized Studies," Journal of Educational Psychology.

Method: express causal effects through treatment-specific potential outcomes, write causal inference as a missing-data problem, and separate the assignment mechanism from the outcome model.

Data: the paper is a methodological argument and does not rely on any one platform dataset.

Result: the core contribution is not an empirical number but the establishment that "the two potential outcomes of the same unit cannot be observed at once," so a causal design must say how the missing counterfactual is supported by assignment.

Limitation: the framework itself does not automatically supply a valid design. Treatment versions, interference, and dynamic exposure all need extra definition in a concrete application.
:::

### 5.2 Athey and Imbens (2017)

::: {.case}
Paper: "The State of Applied Econometrics: Causality and Policy Evaluation," Journal of Economic Perspectives.

Method: survey experimental, selection-on-observables, IV, RDD, DiD, and synthetic-control designs, stressing that the credibility of modern empirical work comes from an explicit identification strategy.

Data: the review spans applications in education, labor, public economics, and policy evaluation.

Result: the formulas for different estimators are not the heart of the taxonomy; what matters is how assignment or institutional variation constructs the counterfactual.

Limitation: the review is broad in coverage and cannot substitute for a dedicated treatment of interference, measurement, and external validity for each design.
:::

### 5.3 Imbens (2020)

::: {.case}
Paper: "Potential Outcome and Directed Acyclic Graph Approaches to Causality: Relevance for Empirical Practice in Economics," Journal of Economic Literature.

Method: compare potential-outcome and DAG/structural-equation approaches, and discuss the different strengths of the two languages in economic empirical work.

Data: a methodological review that uses several economic designs as examples and reports no single treatment effect.

Result: potential outcomes are especially direct for making the intervention, estimand, and design explicit; DAGs are valuable for expressing complex conditional-independence assumptions. The two can complement each other.

Limitation: any representation relies on substantive assumptions. Drawing the graph completely does not mean the arrows have been verified by the data.
:::

## 6 A Complete Walkthrough on Aurora Data

### 6.1 The DGP

Aurora has 12000 merchants, and the platform randomly opens the tool early to about half of them within 120 regions. The structural direct effect is $0.080$; the tool raises engagement by $0.50$ on average, and the marginal effect of engagement is $0.12$, so the total effect is $0.140$.

```r
set.seed(1010)
N <- 12000
baseline_quality <- rnorm(N)
region <- sample(1:120, N, replace = TRUE)
u_service <- rnorm(N)
early_access <- ave(runif(N), region,
  FUN = function(x) as.integer(x <= median(x)))
engagement <- 0.50*early_access + 0.30*baseline_quality +
  0.60*u_service + rnorm(N, sd=0.70)
log_sales <- 2 + 0.40*baseline_quality + 0.08*early_access +
  0.12*engagement + 0.50*u_service + rnorm(N, sd=0.60)
dat <- data.frame(seller_id = seq_len(N), region, baseline_quality,
  early_access, engagement, log_sales)
```

### 6.2 Three Specifications

```r
fit_diff <- lm(log_sales ~ early_access, data=dat)
fit_pre  <- lm(log_sales ~ early_access + baseline_quality + factor(region), data=dat)
fit_bad  <- lm(log_sales ~ early_access + baseline_quality + engagement +
              factor(region), data=dat)
```

The actual output is as follows.

| Specification | Estimate | SE | Target and question |
|---|---:|---:|---|
| Difference in means | 0.149 | 0.017 | Total effect of random assignment |
| Baseline adjusted | 0.157 | 0.015 | The same total effect, higher precision |
| Post-treatment control | -0.099 | 0.014 | Cuts the mediator and introduces collider bias |

![Estimates of the early-access coefficient under three specifications. The first two center on the true total effect; controlling for post-treatment engagement flips the coefficient negative. The structural direct effect shown by the purple line is not identified by that bad-control regression either.](assets/fig/fig_10_estimand_contrast.svg)

### 6.3 The Full Reconciliation

The conflict in the opening is not that the experiment gives two answers, but that the third regression changed the conditional comparison. Both the difference in means and the baseline adjustment use the original randomization, estimating $0.149$ and $0.157$, close to the true total effect $0.140$. Controlling for engagement gives $-0.099$, which is neither the total effect nor the structural direct effect $0.080$. One more significant control makes the answer worse, not more credible.

## 7 Failure Modes and Robustness

The first kind of failure is that the treatment is unclear. Gaining the tool, opening the tool, accepting a recommendation, and following a recommendation are different interventions. If a paper says "AI adoption" in its motivation while the variable is only access exposure, the conclusion must narrow to an assignment effect.

The second kind of failure is temporal misalignment. Using firm size, review score, or engagement updated after treatment as controls will control a mechanism, manufacture selection, or both at once. Building a table of variable timestamps usually prevents error better than running another set of fixed effects.

The third kind of failure is a mismatch between the unit and the randomization level. If the policy is randomized by city while the regression treats millions of users as independent observations, the standard errors will be severely understated. The unit of analysis can be finer than the randomization unit, but inference must respect the assignment dependence.

The fourth kind of failure is treating balance, placebo, or pre-trend as proof of the assumptions. They can reveal obvious implementation errors, but they cannot verify all unobserved counterfactuals. A research report should separate design evidence, diagnostic evidence, and untestable assumptions.

The fifth kind of failure is mechanism overclaim. Treatment changes the mediator and the mediator predicts the outcome, which provides only mechanism-consistency evidence. Identifying mediation requires extra assignment, sequential ignorability, or other structural constraints.

Last is external validity. Aurora identifies the short-run effect for eligible pilot merchants; it does not automatically answer the equilibrium consequences of all merchants, a mature version, or a whole-market rollout. Population, version, and horizon are all part of the estimand, not a range that the last sentence of the results section may freely widen.

Threading these problems together, the core discipline of research design is a single sentence: first state where the counterfactual comes from, then decide which estimator to summarize it with. When the design is unclear, more controls, fancier machine learning, and smaller p-values will not fill in the missing counterfactual.

## 8 Further Reading

::: {.readings}
Required reading, in the suggested order:

- Rubin (1974). Uses potential outcomes to write the causal question as a missing-counterfactual problem.
- Athey and Imbens (2017). Understands the identification strategy from the whole landscape of applied econometrics.
- Imbens (2020). Compares the applicable boundaries of potential outcomes and DAGs.
- Hernán and Robins (2020). Builds a complete research-design language from the target trial, time-varying treatment, and causal diagrams.

Further reading:

- Angrist and Pischke (2009). Understands empirical design from institutional variation and the experimental ideal.
- Imbens and Rubin (2015). Treats assignment mechanisms, randomization inference, and observational designs systematically.
:::

::: {.apa-refs}
- Angrist, J. D., & Pischke, J.-S. (2009). *Mostly harmless econometrics: An empiricist's companion*. Princeton University Press.
- Athey, S., & Imbens, G. W. (2017). The state of applied econometrics: Causality and policy evaluation. *Journal of Economic Perspectives, 31*(2), 3-32. https://doi.org/10.1257/jep.31.2.3
- Hernán, M. A., & Robins, J. M. (2020). *Causal inference: What if*. Chapman & Hall/CRC.
- Imbens, G. W. (2020). Potential outcome and directed acyclic graph approaches to causality: Relevance for empirical practice in economics. *Journal of Economic Literature, 58*(4), 1129-1179. https://doi.org/10.1257/jel.20191597
- Imbens, G. W., & Rubin, D. B. (2015). *Causal inference for statistics, social, and biomedical sciences: An introduction*. Cambridge University Press. https://doi.org/10.1017/CBO9781139025751
- Rubin, D. B. (1974). Estimating causal effects of treatments in randomized and nonrandomized studies. *Journal of Educational Psychology, 66*(5), 688-701. https://doi.org/10.1037/h0037350
:::
