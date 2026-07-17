---
title: "Interference, Marketplace Experiments & Equilibrium Effects"
subtitle: "Causal Design When One Unit's Treatment Changes Another's Outcome"
seriesline: "Foundations of Information Systems Economics · Chapter 18"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 18 · Interference, Marketplace Experiments & Equilibrium Effects"
---

## Introduction

Meridian randomly waives the delivery fee on half of its orders. The experiment runs flawlessly and the estimate is extremely precise: the treatment effect is $0.101$. On that basis the platform rolls the policy out to everyone, and the average outcome rises by only about $0.040$. This is not an experiment that failed. It is the market answering a question the researcher never asked.

When half of orders are subsidized, treated and control orders share the same pool of couriers. The subsidy delivers a direct benefit of $0.100$ to an individual order, but the higher treated share creates congestion, so every order bears a negative externality. The individual randomized experiment compares two kinds of orders inside the same mixed market. The full rollout changes the state of the entire market. Both are causal effects, but they are not the same causal effect.

Once one person's treatment enters another person's potential outcome, writing down "a randomized experiment" is no longer enough. You must also state the randomization share, the boundary of interaction, the exposure mapping, and the target rollout scale. This chapter uses Meridian to separate direct, indirect, total, and overall effects, and it compares individual, market-cluster, two-stage saturation, and switchback designs. The cluster and switchback designs deliver $0.044$ and $0.043$, closer to the full rollout of $0.040$, but their standard errors are wider too: the effective unit of randomization has moved from the order to the market-period. There is no free lunch here, only an experiment that is more honest, and more expensive, about the right estimand.

## 1 An "Unbiased" $0.101$ That Answers the Wrong Question

Within each market and period, Meridian randomly gives half of orders free delivery. We regress the outcome on individual treatment while controlling for market and period fixed effects, and the coefficient is $0.101$ with an SE of $0.002$. The randomization is executed correctly and the standard error is small.

The platform then pushes the policy to all orders, and the average outcome rises by only $0.040$. Did the experiment fail? No. The individual contrast estimates the direct difference between treated and control orders at the same level of congestion. The full rollout also includes the congestion effect of $-0.060$ that every order bears in common. The experiment is unbiased for its own estimand, but it was used to answer a different policy estimand. Changing the granularity of randomization also changes the effective sample size: 192000 rows of order data are not 192000 independent cluster assignments.

::: {.warning}
Under interference, "a randomized experiment" is not a complete description of the estimand. You must also state the randomization share, the exposure mapping, the cluster boundary, and the target policy saturation.
:::

## 2 Exposure-dependent potential outcomes

### 2.1 From $Y_i(d)$ to $Y_i(d_i,\mathbf d_{-i})$

Without interference, unit $i$'s outcome depends only on its own assignment: $Y_i(d_i)$. With interference,

$$Y_i(\mathbf d)=Y_i(d_i,\mathbf d_{-i}),$$

the dimension of the assignment vector grows with the sample, and it cannot be identified element by element. An exposure mapping compresses the neighbors' assignments into a finite state $G_i=g_i(\mathbf d_{-i})$, written as $Y_i(d_i,G_i)$. In a marketplace, $G_i$ might be the treated share in the same market, the number of nearby treated sellers, or the treatment history over recent periods. Throughout this chapter $g$ always denotes the exposure state and is never used as a group subscript.

### 2.2 Four kinds of effects

When exposure changes from $g_0$ to $g_1$, the direct effect holds neighbors' exposure fixed:

$$\tau_D(g)=\mathbb E[Y_i(1,g)-Y_i(0,g)].$$

The indirect effect holds a unit's own treatment fixed and changes the neighbors' exposure:

$$\tau_I(d)=\mathbb E[Y_i(d,g_1)-Y_i(d,g_0)].$$

The total effect changes both a unit's own treatment and its environment, while the overall effect compares the population mean under two assignment policies. Meridian's direct effect is $0.100$, the congestion effect from zero saturation to full saturation is $-0.060$, so the global rollout effect is $0.040$.

### 2.3 Partial interference

The most common simplification assumes that interference occurs only within a cluster and does not cross clusters. Markets, schools, villages, or geographic cells are often treated as clusters. This assumption must be supported by transaction boundaries and mobility. Couriers taking orders across zones, consumers searching across cities, or sellers operating on multiple platforms all break artificial boundaries.

::: {.assumption data-label="Partial interference"}
Unit $i$'s potential outcome depends only on the assignment vector within its own cluster; there is no treatment spillover between different clusters.
:::

### 2.4 Congestion and business stealing

Marketplace interference usually runs in two directions. Congestion makes treatment compete for a limited supply, so the higher the treated share, the lower the service quality for everyone. Business stealing shifts transactions from control sellers to treated sellers, so the treated-control contrast also contains redistribution. The platform's total GMV may be almost unchanged while the seller-level treatment effect is large.

### 2.5 Short-run and equilibrium effects

Within the experimental window, supply, prices, and participation are usually roughly fixed. A long-run rollout changes driver logins, seller entry, prices, inventory, and user expectations. Even if a cluster design identifies the short-run total effect, it does not automatically identify the long-run equilibrium effect. The latter requires longer experiments, different saturations, institutional extrapolation, or a structural model for support.

To summarize this section: interference makes potential outcomes depend on the assignment environment; direct, indirect, total, and overall effects are not interchangeable; and partial interference and the exposure mapping are dimension-reducing assumptions that must be defended.

## 3 Identification: how design matches exposure

### 3.1 Individual randomization

When about half of orders are treated in every market, individual randomization gives the direct contrast at a fixed saturation. The shared congestion is the same for both groups and cancels in the difference, which is why Meridian estimates $0.100$. If the goal is exactly "who benefits under a half rollout," it is entirely appropriate; if the goal is an all-versus-none policy, it is not enough information.

### 3.2 Cluster randomization

Assigning the whole market-period to treatment or control makes the assignment uniform within a cluster, which lets us compare the total effect of treated markets against control markets. The cost is that the effective unit of randomization drops from the number of orders to the number of market-periods, and the standard error must reflect dependence at the cluster-assignment level.

::: {.assumption data-label="Cluster isolation and exchangeability"}
Cluster assignments are randomized by a known mechanism, there is no spillover between clusters, and the cluster definition is fixed before treatment.
:::

Clusters drawn too small leave residual cross-cluster spillover, while clusters drawn too large cause the effective sample to collapse. Geographic clustering is a bias-variance tradeoff, not a case of bigger is better.

### 3.3 Two-stage randomization

The first stage randomizes clusters to different saturation levels, for example $0$, $0.5$, $1$; the second stage randomizes units within a cluster at the assigned saturation. This lets us estimate the direct effect at a given saturation and also compare the spillover of different saturations onto untreated units. At half saturation, Meridian's spillover onto control is $-0.030$, and the all-versus-none total effect is $0.040$.

### 3.4 Graph-cluster randomization

If interference propagates along a buyer-seller, friendship, or transaction graph, we can group densely connected nodes into the same cluster to reduce cross-group edges. Exposure probabilities are often unequal and require Horvitz-Thompson or Hájek weighting. Whether the graph is defined before treatment or is itself changed by treatment must also be stated.

### 3.5 Switchback experiments

When an entire market must share a single policy, we can switch between treatment and control over time. Each market serves as its own control, which suits ride-hailing pricing, dispatch, and marketplace ranking. Identification requires the assignment schedule to be randomized and requires handling time trends, serial correlation, and carryover.

![Switchback schedule for six markets. Each market switches between treatment and control in two-period blocks, with a randomized starting state; the same market supplies its counterfactual in different periods.](assets/fig/fig_18_switchback.svg)

### 3.6 Carryover and washout

If a subsidy still affects orders, courier positions, or user memory after it ends, the current outcome depends on treatment history. A naive comparison of the current assignment is biased. You should specify a carryover horizon in advance, add washout periods, use longer blocks, or extend the exposure mapping to lagged assignments. Longer blocks reduce carryover bias but reduce the number of switches and the effective sample.

### 3.7 Testability

Randomization code, exposure probabilities, cluster contamination, and treatment delivery can be audited. The absence of any cross-cluster spillover, the adequacy of the chosen exposure mapping, and the isomorphism between long-run equilibrium and the short-run experiment cannot be verified from the current outcome alone. You can provide diagnostic evidence with distance gradients, boundary units, different cluster sizes, and the saturation response.

### 3.8 Design effect and effective sample size

The cost of cluster randomization is not just "the SE is larger." Under the approximation that clusters are of roughly equal size and within-cluster correlation can be summarized by a common ICC, with average cluster size $m$ and intracluster correlation $\rho$, the design effect is

$$DE=1+(m-1)\rho.$$

Relative to individual randomization, the effective sample size is about $N/DE$. Even when $\rho$ is small, the loss can be substantial when $m$ is large. A marketplace design should choose the cluster size to balance bias reduction against variance inflation, and should use pre-period data to estimate the ICC, serial correlation, and cross-boundary traffic for the power calculation.

### 3.9 Supply response and the equilibrium horizon

A short-run switchback holds most of the supply fixed. If the subsidy persists, couriers may increase their online hours, sellers may enter, consumers may form new habits, and the congestion effect will change with the horizon. We can write the policy effect as $\tau(h,s)$, where $h$ is the rollout horizon and $s$ is the saturation. A single two-week experiment at $s=0.5$ identifies only the neighborhood of that support and does not provide the entire response surface.

An extended design can randomize different saturations and lengthen the treatment duration in some clusters, separately identifying short-run congestion and medium-run supply elasticity. If long-run entry cannot be covered by experiment, then "the long-run equilibrium" should be handed explicitly to a model or to external evidence, rather than letting a short-run A/B coefficient overreach.

## 4 Estimation: different designs imply different weights

### 4.1 Individual difference

At a fixed saturation $s$, the treated-control difference estimates $\tau_D(s)$. Ordinary user-level robust SEs ignore shared shocks; at a minimum you should treat correlation at the level of the assignment or interference cluster.

### 4.2 Cluster-level estimator

The most transparent approach first aggregates the outcome to the unit of randomization and then compares treated and control clusters. Regression adjustment can improve precision, but the number of clusters, not the number of rows, drives the leading asymptotic approximation. When there are few clusters, use randomization inference or small-cluster corrections.

### 4.3 Two-stage Horvitz-Thompson

For exposure condition $(d,g)$, weight by its known exposure probability $\pi_i(d,g)$:

$$\widehat\mu(d,g)=\frac{1}{N}\sum_i\frac{\mathbf 1\{D_i=d,G_i=g\}Y_i}{\pi_i(d,g)}.$$

Differences among the various $\widehat\mu(d,g)$ construct the direct and spillover effects. Overlap requires that the target exposure has positive probability for the relevant units.

### 4.4 Switchback regression

We can regress the market-period aggregate outcome on current treatment, market FE, time FE, and pre-specified lag exposures, and infer at the market or randomization-block level. If the schedule is randomized by block, we cannot treat each order as an independent randomization.

### 4.5 Policy value curve

Multiple saturation levels let us estimate the curve of the population outcome as a function of the rollout share. It is better suited than a single A/B coefficient to answering "how far to roll out," but it is credible only within the saturation support the experiment covers; extrapolating from $0.5$ to $1$ still requires a shape assumption.

### 4.6 Saturation design and the response surface

The value of two-stage randomization is not just "correcting the standard error" but actively varying each cluster's treatment saturation. First randomize markets to $s_c\in\{0,0.25,0.5,0.75,1\}$, then assign individuals within a market according to $s_c$, which lets us estimate the response surface of the outcome as a function of own treatment and market saturation:

$$m(d,s)=E[Y_i(d,s_c=s)].$$

Here $m(1,s)-m(0,s)$ describes the direct effect at a given saturation; $m(0,s)-m(0,0)$ describes the spillover borne by untreated units; and $m(1,1)-m(0,0)$ is closer to the full-rollout contrast. If treatment competes for a shared resource, $m(0,s)$ will fall with $s$; if there is knowledge diffusion or a network gain, it may rise. The saturation curve is therefore itself an empirical prediction about the economic mechanism, not merely a technical design choice.

Estimation must use the two-stage assignment probabilities. A treated unit in a high-saturation market cannot be compared with equal weight to a control unit in a low-saturation market, because they represent different exposure cells. For rare $(d,s)$ cells, report positivity and the effective sample size; with no market close to $s=1$, you cannot extrapolate the full-deployment effect from a functional form and call it an experimental estimate.

### 4.7 Bipartite interference: treat one side, measure the outcome on the other

Platforms often have two kinds of nodes: drivers and riders, sellers and buyers, advertisers and consumers. Treatment may be applied on the supply side while the outcome is observed on the demand side. In this case an ordinary social-network exposure mapping is not enough, and a bipartite graph is needed. Let $A_{ij}$ indicate whether supply node $j$ might serve demand node $i$; the demand-side exposure can be written as

$$G_i=\frac{\sum_j A_{ij}D_jw_{ij}}{\sum_j A_{ij}w_{ij}},$$

where $w_{ij}$ can be defined by historical match probabilities or capacity contributions. After AI dispatch is opened to some drivers, a rider's waiting time depends on the treated share among the drivers reachable to that rider, not on the rider's own treatment.

This design requires the assignment graph to be constructed before treatment. If neighbors are defined by actual matching after treatment, the matching itself has already been affected by treatment, producing post-treatment selection. You can build the exposure from a historical window, geographic reachability, or a pre-defined candidate set, and do sensitivity analysis over different graph definitions. In theory you should also distinguish extensive-margin spillover (whether service is obtained) from intensive-margin spillover (waiting time or price).

### 4.8 From a short-run experiment to a market-equilibrium counterfactual

Even if a cluster experiment unbiasedly identifies the short-run market-level effect, it does not automatically equal the long-run equilibrium effect. Supply entry, exit, cross-zone migration, consumer learning, and competitor responses usually take longer. Let the short-run state be $S_0$, with the post-policy transition $S_{t+1}=T(S_t,D,\varepsilon_{t+1})$; the short-run experiment estimates the one-step response given $S_0$, while the long-run policy cares about the average welfare under the new steady-state distribution.

Empirically we can proceed at three levels. First, extend the experiment and plot the duration response, checking whether the effect drifts with exposure age. Second, measure state variables such as active supply, cross-zone mobility, and repeat use, treating them as outcomes of equilibrium adjustment rather than as controls. Third, compare clusters of different sizes and different degrees of boundary openness; if a larger cluster gives an effect closer to the platform-wide counterfactual, that indicates cross-boundary substitution is an important mechanism. If the long-run state is never observed in the experiment, the final extrapolation must depend explicitly on an economic model or a structural assumption and cannot be obtained from cluster-robust SEs alone.

## 5 Anchor papers

### 5.1 Hudgens and Halloran (2008)

::: {.case}
Paper: "Toward Causal Inference With Interference," Journal of the American Statistical Association.

Method: under partial interference and two-stage randomized designs, formally defines direct, indirect, total, and overall effects and gives randomization-based estimators.

Data: illustrated with infectious-disease vaccine and group-intervention scenarios showing how one person's treatment protects others.

Result: different effects correspond to different assignment-policy contrasts and cannot be summarized by a single individual effect.

Limitation: partial interference relies on the correct clusters; cross-group transmission breaks identification.
:::

### 5.2 Aronow and Samii (2017)

::: {.case}
Paper: "Estimating Average Causal Effects Under General Interference," Annals of Applied Statistics.

Method: uses an exposure mapping to compress the high-dimensional assignment vector into a finite set of exposure conditions and builds design-unbiased estimators weighted by exposure probabilities.

Data: the paper demonstrates the method with network and spatial interference examples.

Result: as long as the exposure mapping and assignment probabilities are known, several exposure-specific average effects can be identified.

Limitation: a wrong exposure mapping merges distinct potential outcomes; rare exposures produce extreme weights and low precision.
:::

### 5.3 Holtz et al. (2025)

::: {.case}
Paper: "Reducing Interference Bias in Online Marketplace Experiments Using Cluster Randomization: Evidence from a Pricing Meta-experiment on Airbnb," Management Science.

Method: runs individual-randomized and cluster-randomized pricing experiments simultaneously on Airbnb, directly comparing the two designs.

Data: a platform field experiment in the Airbnb listing marketplace.

Result: at least about $20\%$ of the total-effect estimate from the individual design can be attributed to interference bias, which cluster randomization removes.

Limitation: the cluster design lowers bias but increases variance, and the magnitude of the result depends on that platform's substitution network and cluster construction.
:::

## 6 A full walkthrough on the Meridian data

### 6.1 DGP

Each design has 120 markets, 20 periods, and 80 orders per market-period, for 192000 observations in total. The outcome response is

$$Y_{imt}=\alpha_m+\lambda_t+0.10D_{imt}-0.06S_{mt}+\varepsilon_{imt},$$

where $S_{mt}$ is the market-period treated share. The individual direct effect is $0.100$, and full rollout relative to zero rollout is $0.040$.

### 6.2 Three designs

| Design | Estimate | SE | Corresponding estimand |
|---|---:|---:|---|
| Individual randomization | 0.101 | 0.002 | Half-saturation direct effect |
| Market-period cluster | 0.044 | 0.005 | All-versus-none total effect |
| Switchback | 0.043 | 0.005 | All-versus-none period effect |

![Treatment contrasts for the three randomization designs. The individual design is correct for the direct effect but clearly exceeds the global rollout effect; the cluster and switchback designs change the entire market environment and estimate close to the full-rollout effect.](assets/fig/fig_18_market_interference.svg)

### 6.3 Two-stage estimands

At half saturation, the direct effect is $0.100$; the spillover onto control units relative to zero-saturation markets is $-0.030$; and the all-versus-none total effect is $0.040$. All three numbers are correct within the same DGP, and the difference comes from the policy contrast.

### 6.4 Full reconciliation

The individual experiment's $0.101$ is not statistical bias but an estimand mismatch. The market-period cluster's $0.044$ and the switchback's $0.043$ target the change in the whole-market state and are close to the true value of $0.040$; both have an SE of $0.005$, wider than the individual design's $0.002$, which makes the variance cost of changing the unit of randomization truly visible. A study that only writes "treatment increased outcomes by 10.1%," without stating the saturation and exposure environment, will mis-sell a local direct effect as a rollout effect.

## 7 Failure modes and robustness

First, cluster leakage. When units trade or move between clusters, cluster randomization still has contamination. You should report the share of cross-boundary transactions and do a boundary sensitivity analysis.

Second, the wrong level of inference. Hundreds of thousands of orders are not hundreds of thousands of independent assignments. The SE must respect the market-period or block randomization.

Third, ignored carryover. Switchback blocks that are too short let the previous state contaminate the current outcome. Varying the block length and the washout rule is a necessary sensitivity check.

Fourth, an exposure mapping that is too coarse. Using only the treated-neighbor count may ignore edge weights, capacity, and direction. You should compare theoretically reasonable mappings rather than fitting the best graph to the outcome.

Fifth, business stealing mistaken for total creation. Seller-level gains may be nothing more than transfers from controls. The platform should report both the transaction redistribution and the market-aggregate outcome.

Sixth, extrapolating a short-run experiment to the long-run equilibrium. Supply entry, price response, and user learning usually fall outside a short experiment. You should narrow the conclusion, extend the horizon, or explicitly introduce an equilibrium model.

Finally, cluster randomization is not the default answer. If interference is weak, the individual design is more precise; if the policy can only be implemented market-wide, a switchback is more natural; if you need to separate direct from spillover effects, only a two-stage saturation design provides the required variation. The design choice is determined by the estimand and the interference mechanism.

## 8 Further reading

::: {.readings}
Required reading:

- Hudgens and Halloran (2008). The foundational definitions of direct, indirect, total, and overall effects.
- Aronow and Samii (2017). Exposure mappings and general-interference estimation.
- Holtz et al. (2025). A direct meta-experiment of individual versus cluster designs in a marketplace.

Further reading:

- Baird et al. (2018). The design and analysis of two-stage randomized experiments.
- Bojinov, Simchi-Levi, and Zhao (2023). Switchback design, carryover, and inference.
- Johari et al. (2022). Bias and variance in two-sided marketplace experimentation.
:::

::: {.apa-refs}
- Aronow, P. M., & Samii, C. (2017). Estimating average causal effects under general interference, with application to a social network experiment. *The Annals of Applied Statistics, 11*(4), 1912-1947. https://doi.org/10.1214/16-AOAS1005
- Baird, S., Bohren, J. A., McIntosh, C., & Özler, B. (2018). Optimal design of experiments in the presence of interference. *The Review of Economics and Statistics, 100*(5), 844-860. https://doi.org/10.1162/rest_a_00716
- Bojinov, I., Simchi-Levi, D., & Zhao, J. (2023). Design and analysis of switchback experiments. *Management Science, 69*(7), 3759-3777.
- Holtz, D., Lobel, F., Lobel, R., Liskovich, I., & Aral, S. (2025). Reducing interference bias in online marketplace experiments using cluster randomization: Evidence from a pricing meta-experiment on Airbnb. *Management Science, 71*(1), 390-406. https://doi.org/10.1287/mnsc.2020.01157
- Hudgens, M. G., & Halloran, M. E. (2008). Toward causal inference with interference. *Journal of the American Statistical Association, 103*(482), 832-842. https://doi.org/10.1198/016214508000000292
- Johari, R., Li, H., Liskovich, I., & Weintraub, G. Y. (2022). Experimental design in two-sided platforms: An analysis of bias. *Management Science, 68*(10), 7069-7089.
:::
