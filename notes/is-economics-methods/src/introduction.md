---
title: "Introduction & How to Use"
subtitle: "A Roadmap to the Series"
seriesline: "Foundations of Information Systems Economics · Introduction"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Introduction & How to Use"
---

## Introduction

When you study digital platforms and artificial intelligence, the hardest problem is usually not that you have no method, but that you have too many. A sales jump after an algorithm goes live can be read as a causal effect, or it can be nothing more than selection into adoption. A change in transaction volume after a platform tweaks its rules might come from how consumers respond, or it might come from sellers, competitors, and the recommendation system all adjusting at the same time. Without economic theory, it is hard to say what you should even be comparing. Without an identification design, the prettiest correlation need not carry any causal meaning. And without a grasp of market structure and institutional context, a local effect is a poor guide to questions about competition, welfare, and policy. These troubles usually show up together, yet they are typically scattered across different courses and textbooks.

This series is written for researchers working in the economics of information systems and in empirical operations management, with a focus on the business use of artificial intelligence, digital platforms, information environments, technology adoption, and the market behavior of consumers and firms. The whole book starts from research questions. We first use graduate microeconomic theory to write a real-world dilemma as a problem of choice, equilibrium, and incentives. We then use microeconometrics and causal inference to judge what the data can actually identify. We go on to use empirical industrial organization and structural econometrics to analyze substitution, competitive conduct, and counterfactual policy. And we finally carry these tools back into the specific institutions of information, technology, and platform economics. The goal is not to mold you into a specialist who knows only one method, but to build a research capacity that has both depth and lateral connection: expertise in the economics of information, AI, platforms, and technology, together with the ability to work with economists, marketing scholars, and operations management scholars in a shared language of theory, econometrics, and computation.

## How to Read This Series

The series has twenty-seven numbered modules, grouped into five parts, plus this unnumbered introduction. The map below lays out the whole thing, one sentence per chapter. You do not need to grind through it front to back. The reading paths further down give a few routes tailored to different goals.

**Part I, Microeconomic Foundations**: supplies the primitives (utility, welfare, equilibrium, incentives) that later identification and estimation draw on.

- Chapter 1, Consumer, Producer, and Welfare Theory: UMP/EMP duality, Slutsky, Roy/Shephard, CV/EV and consumer surplus, cost-production duality, and choice under uncertainty.
- Chapter 2, Game Theory: normal and extensive form, Nash equilibrium, Bayesian games, subgame perfection, repeated games, and the folk theorem.
- Chapter 3, Information Economics, Contracts, and Mechanism Design: adverse selection, signaling, screening, moral hazard, the revelation principle, VCG, and auction theory.
- Chapter 4, Pricing, Monopoly, and Imperfect Competition: price discrimination, the Coase conjecture, Cournot/Bertrand/Hotelling, collusion, entry deterrence, and welfare.

**Part II, The Econometric Core**: the machinery of estimation and inference that everything downstream cites.

- Chapter 5, Estimation: MLE, M-estimation, and GMM: the unified extremum-estimator framework, consistency and asymptotic normality, the sandwich variance, GMM efficient weighting, and weak identification.
- Chapter 6, Inference: the Delta Method, the Bootstrap, Clustering, and Testing: the Delta method, bootstrap refinements, HAC, two-step inference, randomization inference, and W/LR/LM.
- Chapter 7, Panel Data and Dynamic Panels: FE/RE, strict and sequential exogeneity, Nickell bias, and Arellano-Bond and Blundell-Bond.
- Chapter 8, Nonparametric and Semiparametric Methods: kernel estimation, local polynomials, series/sieve, Robinson, and the bias-variance tradeoff.

**Part III, The Causal Identification Toolbox**: one module per estimand family, the workhorse of reduced-form causal research.

- Chapter 9, Potential Outcomes, Experiments, and Selection on Observables: the Rubin framework, RCTs and randomization inference, matching, IPW/AIPW, and Roy selection.
- Chapter 10, From an Economic Question to a Research Design: unit, treatment, timing, exposure, estimand, assignment mechanism, and bad controls.
- Chapter 11, Measurement, Digital Traces, and Algorithmically Generated Data: construct validity, logging, selective labels, the validation sample, and prediction-powered inference.
- Chapter 12, IV: LATE, MTE, and the Control Function: the LATE theorem, compliers, modern practice with weak instruments, MTE as the main object, and the control function.
- Chapter 13, DiD and Event Studies: from 2x2 to staggered, the Goodman-Bacon decomposition, and Callaway-Sant'Anna / Sun-Abraham / imputation.
- Chapter 14, RDD and Synthetic Control: continuity-based identification, rdrobust bias correction, fuzzy RD, synthetic control, and augmented SC/SDID.
- Chapter 15, Sensitivity, Negative Controls, and Partial Identification: OVB sensitivity, the robustness value, negative controls, bounds, and the breakdown frontier.
- Chapter 16, Causal ML: DML, Causal Forest, and Policy Learning: Neyman orthogonality, cross-fitting, GATES/CLAN, and policy learning.
- Chapter 17, Digital Experiment Design: power/MDE, CUPED, sequential testing, multiple testing, and experiment governance.
- Chapter 18, Interference, Marketplace Experiments, and Equilibrium Effects: exposure mappings, two-stage randomization, cluster designs, switchback, and rollout effects.

**Part IV, Structural Methods**: consumes the theory from Part I and the identification from Part III, then structures a model so it can do counterfactuals.

- Chapter 19, Discrete Choice: from RUM to Mixed Logit: Gumbel to logit, identification, IIA, nested/GEV, mixed logit, endogeneity, and the control function.
- Chapter 20, Demand Estimation: Berry Inversion and BLP: Berry inversion, the contraction, NFXP and MPEC, instruments, micro-BLP, and identification.
- Chapter 21, The Supply Side: Markups, Conduct, Mergers, and Welfare: the multiproduct Bertrand FOC, markup recovery, GUPPI/UPP, merger simulation, and conduct testing.
- Chapter 22, Entry, Market Structure, and Partial Identification: the Bresnahan-Reiss threshold, Berry 1992, Seim, Ciliberto-Tamer, and moment inequalities.
- Chapter 23, Dynamic Structural Models: Bellman, Rust NFXP, Hotz-Miller CCP, dynamic demand, and an introduction to dynamic games.

**Part V, Topics in the Digital Economy**: brings all three legs down onto the phenomena of the field.

- Chapter 24, Platforms, Network Effects, and Switching Costs: two-sided pricing, multihoming, the competitive bottleneck, compatibility, P2P/commons governance, and state dependence.
- Chapter 25, Consumer Search, Advertising, and Information Goods: the Diamond paradox, search-cost estimation, advertising, experience goods, and Bayesian learning.
- Chapter 26, AI, Algorithms, and the Economics of Data: algorithmic pricing and collusion, AI adoption, text-as-data and inference with generated regressors, and data network effects.
- Chapter 27, Quick Reference and Tools: formula cards, method decision trees, cross-software command tables, and a checklist of common pitfalls.

## Prerequisites

The main text assumes you have a foundation in probability and mathematical statistics, linear algebra, and one semester of asymptotic theory (the law of large numbers, the central limit theorem, Slutsky's theorem). You should be comfortable with expectation, variance, and conditional expectation, and able to read OLS in matrix form. On the computing side, we assume basic R: you can run a regression, read the output, and install a package. The theory chapters (Part I) ask for some comfort with constrained optimization and a bit of real analysis, but we try to give enough intuition at every derivation that a thinner background can still follow along. You do not need to know IO or structural estimation in advance. That is exactly what this series sets out to teach.

## Reading Paths

People with different goals can take different routes. There is no need to grind through in numbered order.

If you mainly do reduced-form causal empirical work, your main line runs from Part II to Part III. First build a solid base in estimation and inference with Chapters 5 and 6, then read the causal language of Chapter 9, the research design of Chapter 10, and the measurement of Chapter 11. After that, pick from Chapters 12 through 18 according to your research question. Add Chapter 7 when you need panels, and Chapter 8 when you need a nonparametric foundation. Finishing this line should leave you able to design, diagnose, and referee an empirical platform study on your own.

If you want to do structural estimation, add Part I and Part IV on top of the above. In the theory chapters, Chapter 1 supplies the consumer and welfare, Chapter 2 supplies games, Chapter 3 supplies information and mechanisms, and Chapter 4 supplies imperfect competition. Then Chapters 19 through 23 build structural capacity in order, and Chapters 19, 20, and 21 have strong dependencies, so read them in sequence.

If you only want to build field intuition quickly, you can start with Part V (Chapters 24 through 26) and go back to the relevant method chapters and the quick reference in Chapter 27 whenever you hit an unfamiliar method.

Whichever line you take, Chapter 27 is worth an early glance. It compresses each method's formulas, use cases, software commands, and common pitfalls into cards, and it is the reference you keep at hand when writing papers and refereeing.

## How This Series Is Written

Every method chapter unfolds in the same eight-beat structure, so that you can move smoothly from one method to another. We open with a real causal question and a number that does not look right, then set up the economic model and the estimand, then pull identification out on its own and work it through (this is the most detailed section of each chapter's analysis), then turn to estimation, anchor to one or two real papers, run the code from scratch on a dedicated simulated case, lay out the failure modes and robustness, and close with annotated further reading and APA references.

A few conventions run through the whole book and are worth knowing up front. Method names, technical terminology, mathematics, and code follow standard usage throughout, because your reading, writing, and refereeing all happen in that vocabulary and the terms need to become muscle memory. Every chapter comes with a concrete, situated running case, and every number in the text comes from R code that was actually run, so it is reproducible rather than made up. In style we look to Mostly Harmless Econometrics: we come in from the problem rather than from the definition, we are sparing with formulas and give each one a sentence of intuition, and we lay out the credibility of an assumption and the places where it cannot be tested rather than hiding them. Where things get hard, we go from shallow to deep, giving intuition, geometry, or a small number first, and only then the formalization.

## Notation

The whole book shares one set of notation, and each chapter links back here instead of redefining it. We write an individual as $i$, a product as $j$, a market or time period as $t$, a firm as $f$, and a cohort (the period of first treatment) as $g$. Potential outcomes are $Y_i(d)$, the treatment indicator is $D_i$, the instrument is $Z_i$, and covariates are $X_i$. Among the target estimands, $\tau$ refers generally to a treatment-effect parameter (ATE, ATT, and LATE each carry a subscript or qualifier), and a structural parameter is $\theta$. The propensity score is $e(x) = P(D=1 \mid X=x)$ (in Chapter 1, $e$ separately denotes the expenditure function $e(p,\bar u)$, which is the standard notation of microeconomic theory, so tell them apart by context). In discrete choice and demand, the utility of individual $i$ for product $j$ in market $t$ is written $u_{ijt} = x_{jt}'\beta_i - \alpha_i p_{jt} + \xi_{jt} + \varepsilon_{ijt}$, where $p$ is price and $\xi$ is an unobserved product characteristic. Cost is $c$ and marginal cost is $mc$, and markups are recovered through the ownership matrix $\mathcal{H}$ (with $\mathcal{H}_{jk}=1$ when products $j$ and $k$ belong to the same firm). Where an individual chapter introduces local notation, it declares that notation explicitly at that point.
