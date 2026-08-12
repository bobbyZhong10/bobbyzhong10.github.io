---
title: "Causal ML: DML, Causal Forests & Policy Learning"
subtitle: "From Regularization Bias to Targeting Policies"
seriesline: "Foundations of Information Systems Economics · Chapter 16"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 16 · Causal ML: DML, Causal Forests & Policy Learning"
---

## Introduction

Halcyon Marketplace has hundreds of possible confounders. Whether a merchant adopts the tool depends on a complicated nonlinear combination of them, and so does the outcome. Linear regression looks hopelessly naive, so the researcher lets machine learning predict treatment and outcome separately. The prediction error falls, but the causal estimate does not automatically improve along with it: the small bias that regularization tolerates in the service of prediction can become a first-order error once it enters a treatment effect.

The point of double/debiased machine learning is not to slap an "artificial intelligence" label on causal inference. It is to redesign the estimating equation so that even when the nuisance functions are estimated by flexible algorithms, their small errors affect the target parameter only at second order. Orthogonalization is responsible for lowering the sensitivity, and cross-fitting is responsible for making sure the same sample is not used both to train and to grade. Drop either one, and no matter how pretty the out-of-sample $R^2$ looks, it cannot certify that the treatment effect is credible.

Halcyon has a second problem: the average effect hides who actually benefits. The AI tool may help small merchants while getting in the way of large merchants who already have mature processes. The causal forest pushes the question from a single ATE toward a CATE that varies with features, and policy learning then asks who should get the tool when the budget is limited. But the richer the heterogeneity, the easier it becomes to tell a story out of noise. The spine of this chapter is therefore one cautious question, asked over and over: is machine learning really helping to estimate nuisances, to discover reproducible differences in effects, or is it only packaging overfitting in more modern clothes?

## 1 An Absurdly Large Number

::: {.case}
Halcyon is a fictional two-sided e-commerce platform. It has switched on an AI business-advisory tool for a subset of its merchants, which we will call Growth Playbook: it reads a merchant's historical data and automatically produces recommendations on assortment, pricing, and ad spend. We tell the story with simulated data, and the advantage is that the data-generating process is fully known, so every estimator can be reconciled against the truth, a teaching condition that real data can never give you. The outcome variable $Y$ is the merchant's log GMV growth over the next 90 days (measured in log points). The treatment $D$ is whether the merchant has switched on Growth Playbook.

The key point is that switching on is not random. Halcyon ran no experiment. Instead it pushed the tool first to merchants who were already growing and relatively large, and merchants themselves are more willing to try new things the more active they are. So adoption is entangled with a whole pile of merchant characteristics: size, tenure on the platform, recent demand momentum, and another dozen or so operating dimensions, 20 covariates $X$ in all. Management's question is blunt: how much extra did switching on Growth Playbook make merchants earn?
:::

Start with the most natural procedure: take the average log GMV growth of merchants who switched the tool on, and subtract the average of those who did not. This difference is 0.517 (SE 0.015). Laid out on a slide for management, it means that merchants who used the tool grew 52 log points more, a number big enough to have the product team banging drums.

The platform's data scientist did not believe it too quickly. She knew the lesson of Chapter 9: the merchants who switched the tool on were the stronger ones to begin with, so this difference is laced with selection bias. So she put all 20 covariates into a regression, all linearly, and controlled them out. The coefficient on $D$ fell from 0.517 to 0.450 (SE 0.014). It came down a bit, but it is still absurdly large. Where is the problem? The covariates are indeed correlated with adoption, but the way they affect adoption and the way they affect the outcome are both nonlinear: size and demand momentum interact, and tenure operates with a kink in it. Linear regression can absorb only the linear projection of these confounders, leaving the curved part in the residual, and the confounding in the residual keeps pushing the growth of adopting merchants upward.

Then bring in the machine learning. She fit a random forest to learn the outcome's dependence on all the covariates, tossed the treatment variable in along with them, and read off a "treatment effect" by differencing the predictions, getting 0.170. She also tried something closer to the textbook: first use a forest to fit away the part of the outcome that the covariates can explain, $\mathbb{E}[Y \mid X]$, and then regress the residual on $D$, getting 0.064 this time. Two procedures that used the same powerful tool, one giving 0.170 and one giving 0.064, differ from each other by more than a factor of two, and, as we will see in a moment, both are wrong.

The thing is, we know the truth. Because it is a simulation, the true sample average treatment effect (SATE) is 0.137. The naive comparison of 0.517 overstates it by nearly threefold, the linear control at 0.450 barely rescues anything, and the two machine-learning procedures, one overstating and one understating, both land outside the confidence interval. Four numbers, four directions, and not one of them credible. Why does bringing in the strongest nonlinear fitter produce a pile of mutually contradictory wrong numbers? The answer is the spine of the first half of this chapter: the regularization that machine learning performs in order to control variance leaves a first-order bias in the causal parameter, unless we specifically reconstruct the estimator to remove it. Section 3 takes this bias apart algebraically, and Section 6 shows that the very same forest, used with orthogonalization plus cross-fitting, recovers 0.137 cleanly.

## 2 The Economic Model and the Estimand

Before touching anything, think through two things: which causal quantity we actually want, and what model to use to write down the relationship among treatment, outcome, and high-dimensional confounding.

Start with the target. The effect of switching the tool on looks like a single number, but it is really a question with at least three levels. The first level is the average: for a merchant drawn at random from all of them, how much growth does switching on bring on average, which is the average treatment effect (ATE). The second level is heterogeneity: does the effect differ for small versus large merchants, which is the CATE, a curve that varies with $X$. The third level is the decision: since the effect differs across merchants and the tool is not free, whom should we give it to, which is a policy. The three levels build on one another, and this chapter unfolds them in order, but they share one set of notation and one identification premise, which we set up here.

Follow the potential outcomes of Chapter 9. Each merchant has two potential outcomes, $Y_i(1)$ is the growth it would realize if it switched the tool on, and $Y_i(0)$ is the growth if it did not, and the individual effect $\tau_i = Y_i(1) - Y_i(0)$ varies from merchant to merchant. The observed outcome is $Y_i = Y_i(0) + \tau_i D_i$. The estimands at the three levels are, respectively,

$$\text{CATE}: \tau(x) = \mathbb{E}[Y_i(1) - Y_i(0) \mid X_i = x], \qquad \text{ATE}: \tau_0 = \mathbb{E}[\tau(X_i)].$$

The policy level is formally defined later, in Section 3.6.

Identification still rests on the two assumptions from Chapter 9, carried over verbatim here without reproving them, only stressing how to read them in high dimensions. The first is unconfoundedness (also called conditional independence, or selection on observables): given the covariates, treatment is independent of the potential outcomes, $\big(Y_i(0), Y_i(1)\big) \perp D_i \mid X_i$. It says that $X$ holds everything that affects both adoption and the outcome. Halcyon's adoption is driven by observable characteristics like size, tenure, and momentum, and as long as they are all in $X$, this assumption holds. The second is overlap (also called positivity): every kind of merchant has a nonzero probability of landing in either treatment state, $0 < e(x) < 1$, where

$$e(x) = \Pr(D_i = 1 \mid X_i = x)$$

is the propensity score. Overlap guarantees that we do not hit a dead zone where some class of merchants all adopt and no control can be found. Halcyon's propensity was clipped to $[0.05, 0.95]$ at generation and its realized values land in $[0.10, 0.95]$, so overlap is ample, a point that Section 7 returns to when it fails.

The real novelty of this chapter is not identification but estimation. Unconfoundedness and overlap are word for word the same as in Chapter 9, and identification was already completed there. The difference is that Chapter 9 tacitly assumed we could write $e(x)$ and the outcome regression $\mathbb{E}[Y \mid X, D]$ correctly with a low-dimensional parametric model, whereas this chapter faces up to the reality that these two functions are high-dimensional, nonlinear, and of unknown form, and hands them to machine learning to estimate. The trouble begins right here: these two functions are nuisance parameters, we do not care about them for their own sake but only treat them as stepping stones toward the causal parameter, yet once a stepping stone is estimated by regularized machine learning, it is slippery to step on.

Write this structure into a model that is convenient to analyze, namely the partially linear model (PLR), also called the Robinson model:

$$Y_i = \theta_0\, D_i + g_0(X_i) + u_i, \qquad \mathbb{E}[u_i \mid X_i, D_i] = 0,$$
$$D_i = m_0(X_i) + v_i, \qquad \mathbb{E}[v_i \mid X_i] = 0.$$

The first equation says the outcome consists of a constant treatment-effect coefficient $\theta_0$ plus an unknown confounding function $g_0(X) = \mathbb{E}[Y(0) \mid X]$; the second says the adoption propensity is determined by the propensity function $m_0(X) = e(X)$, with $v_i$ the exogenous variation in adoption not explained by the covariates. $g_0$ and $m_0$ are the nuisance functions that characterize this model. The estimation in Section 3 does not use $g_0$ directly but instead uses the reduced-form conditional mean of the outcome, $\ell_0(X) = \mathbb{E}[Y \mid X] = \theta_0\, m_0(X) + g_0(X)$, residualizing against it; $\ell_0$ and $m_0$ are the two functions to be handed to machine learning ($g_0$ is only an intermediate quantity that makes the model clear). Here we must flag a detail that will be used repeatedly later: the PLR writes the treatment effect as a constant $\theta_0$, yet Halcyon's effect is plainly the heterogeneous $\tau(X)$. This is not a contradiction. Under heterogeneous effects the PLR's $\theta_0$ is no longer equal to the ATE $\tau_0$, but converges to a particular weighted average, with weights proportional to $\text{Var}(D \mid X) = e(X)(1 - e(X))$:

$$\theta_0 = \frac{\mathbb{E}\big[e(X)(1 - e(X))\, \tau(X)\big]}{\mathbb{E}\big[e(X)(1 - e(X))\big]}.$$

This is an overlap-weighted effect, giving the largest weight to those merchants whose adoption is hardest to predict ($e \approx 1/2$), because it is on them that treatment variation is richest and most informative. In Halcyon this weighted effect is 0.144, close to but not the same as the SATE of 0.137. To get exactly the ATE, you have to switch to an estimator that flattens the weights, namely the AIPW-based interactive regression model (IRM) of Section 4.2. First pin down which quantity you want and what weights it corresponds to, then choose the estimator: this order cannot be reversed here any more than in the earlier chapters.

To summarize: identification follows the unconfoundedness and overlap of Chapter 9, and this chapter's whole challenge is in estimation, because the two nuisance functions are high-dimensional and nonlinear and must be estimated by machine learning; the PLR gives a skeleton convenient for analysis, but bear in mind that under heterogeneous effects it estimates an overlap-weighted effect rather than the ATE, the two being 0.144 and 0.137, respectively, in Halcyon.

## 3 Identification and Debiasing: Regularization Bias, Orthogonalization, Cross-Fitting

This is the most finely analyzed section of the chapter. It answers the cliffhanger of Section 1: why plugging machine-learning predictions straight into a causal estimate goes wrong, and how to fix it. The logic runs in three steps: first the diagnosis (where regularization bias comes from), then the first remedy (how Neyman orthogonality makes the estimate immune to error in the nuisances), then the second remedy (how cross-fitting dismantles the self-contamination of overfitting). After these three steps, the problem shifts from the average effect to heterogeneous effects (the causal forest) and to the decision (policy learning). Throughout, we need only the two identification assumptions of Section 2, and no new counterfactual assumption. This section is about bias in estimation, not about whether identification holds.

### 3.1 Regularization Bias: Why Naive Plug-In Goes Wrong

First state the cause of the disease clearly. Consider the second, failed procedure from Section 1. It is really a machine-learning version of an old and beautiful idea, the Frisch-Waugh-Lovell (FWL) theorem. FWL says that to estimate the partial effect of $D$ on $Y$, you can first regress $Y$ and $D$ each on the covariates, take residuals, and then regress the two residuals on each other, and the slope you get is exactly the coefficient on $D$ in the multiple regression. In the notation of Section 2, that means first estimating $\ell_0(X) = \mathbb{E}[Y \mid X]$ and $m_0$, forming residuals $\tilde Y = Y - \hat\ell(X)$ and $\tilde D = D - \hat m(X)$, and then computing

$$\hat\theta = \frac{\sum_i \tilde D_i\, \tilde Y_i}{\sum_i \tilde D_i^2}.$$

In a linear world this is an identity, and it comes out right however you compute it. The trouble is that the cleanliness of FWL holds only when $\hat\ell$ and $\hat m$ are linear projections, and we are now using a forest to estimate them. A forest is biased: to keep variance down it smooths and shrinks the true function, which is regularization. Trace this bias through to $\hat\theta$ and you see a term that refuses to disappear.

Writing it out makes it clearest. The naive procedure residualizes only $Y$ with $\hat\ell$ and then regresses that residual on the treatment, centering the treatment only around its sample mean $\bar D$ (write $\tilde D_i = D_i - \bar D$), without residualizing it against the covariates, which is exactly the 0.064 procedure of Section 1 and precisely the plug-in that the original DML literature warns against. Then

$$\hat\theta - \theta_0 = \underbrace{\frac{\sum_i \tilde D_i\, u_i}{\sum_i \tilde D_i^2}}_{O_p(n^{-1/2})} + \underbrace{\frac{\sum_i \tilde D_i\, \big(\ell_0(X_i) - \hat\ell(X_i)\big)}{\sum_i \tilde D_i^2}}_{\text{regularization bias}} - \underbrace{\theta_0\, \frac{\sum_i \tilde D_i\, m_0(X_i)}{\sum_i \tilde D_i^2}}_{\text{attenuation}}.$$

The first term on the right is labeled $O_p(n^{-1/2})$, ordinary sampling noise that vanishes at rate $n^{-1/2}$, no problem. The real trouble is in the last two terms, and their root cause is the same: the treatment was not residualized against the covariates, so $\tilde D_i$ still carries $m_0(X_i)$, a function of the covariates. The third term is the deadliest. It contains no $\hat\ell$ at all, so even if the outcome side is learned perfectly ($\hat\ell = \ell_0$) it does not budge, and its probability limit is $-\theta_0\, \mathrm{Var}(m_0) / \mathrm{Var}(D)$, an $O(1)$ attenuation bias that does not vanish as the sample grows, pushing the true effect toward zero by the fraction of the treatment's variance explained by the covariates. This is the main reason the 0.064 of Section 1 is far below the truth: it is not that the forest learned badly, it is that the confounding in the treatment was never residualized out, and no matter how well the forest learns, this term is still there. The second term is the regularization bias layered on top: the estimation error $\ell_0 - \hat\ell$ decays at the machine-learning convergence rate, usually slower than $n^{-1/2}$ (the fate of high-dimensional nonparametric estimation, see the bias-variance tradeoff of Chapter 8; forests, kernels, and series estimators are around $n^{-1/4}$, while lasso can be faster under a sparsity assumption), and after multiplying by the un-residualized $\tilde D_i$, $\mathbb{E}[\tilde D_i (\ell_0 - \hat\ell)] \approx \mathbb{E}[(m_0(X) - \bar D)(\ell_0(X) - \hat\ell(X))]$ is generally nonzero, vanishing at about $n^{-1/4}$, slower than the sampling noise, and further amplifying the bias in finite samples. Both terms point to the same prescription: orthogonalize the treatment side as well. When Section 3.2 residualizes $D$ too, it is precisely to wipe out the third term at a stroke and to push the second term down to second order.

::: {.intuition}
Change the metaphor. You want to measure how much higher a table is than the floor, and the ruler in your hand (machine learning) itself has a systematic error, reading a bit too long. If your measurement method is clumsy, the reading of the table's height carries the ruler's bias, and the more the ruler is off, the more wrong your conclusion. Regularization bias is exactly this clumsy method of "inheriting the tool's bias directly into the conclusion." A clever method should make the final reading insensitive to the ruler's bias, so that even if the ruler is a bit off, the table's height is measured accurately. This property of being "insensitive to nuisance error" is exactly the Neyman orthogonality of the next subsection. It does not demand a more accurate ruler (we cannot change machine learning's convergence rate), it demands that the measurement method itself absorb the ruler's first-order error.
:::

### 3.2 Neyman Orthogonality: Making Estimation Immune to Nuisance Error

The first remedy for regularization bias is to build the estimator on a special moment condition, one that is first-order insensitive to small perturbations of the nuisance functions. This property is called Neyman orthogonality.

First look at the general structure. Suppose our estimator is defined by a moment condition $\mathbb{E}[\psi(W; \theta_0, \eta_0)] = 0$, where $W = (Y, D, X)$ is the data, $\theta_0$ is the target parameter, and $\eta_0 = (\ell_0, m_0)$ are the nuisance functions. Neyman orthogonality requires that the Gateaux derivative of this moment function with respect to the nuisances be zero at the truth:

$$\left.\frac{\partial}{\partial t}\, \mathbb{E}\big[\psi(W; \theta_0, \eta_0 + t(\eta - \eta_0))\big]\right|_{t=0} = 0, \quad \forall\, \eta.$$

Here $\eta$ ranges over all reasonable directions of the nuisance functions.

What this says is: push the nuisance a small step from the truth $\eta_0$ in any direction, and the expectation of the moment condition does not move. So when we substitute the biased $\hat\eta$ for $\eta_0$, the first-order effect disappears, and what remains is only second-order (squared) error. Second-order error decays much faster: if $\hat\ell$ and $\hat m$ each converge at $n^{-1/4}$, their product converges at $n^{-1/2}$, exactly the same order as the sampling noise, and no longer dominant. This is why DML allows the nuisances to be estimated by machine learning with slower convergence, and as long as each is faster than $n^{-1/4}$, orthogonality keeps the first-order damage of the slow rate outside the door.

For the PLR, the orthogonal moment function is that residual-on-residual one, but the key is to residualize both nuisances, not just one:

$$\psi(W; \theta, \ell, m) = \big(Y - \ell(X) - \theta (D - m(X))\big)\big(D - m(X)\big).$$

Compared with the naive FWL of Section 3.1, the difference is only in subtracting one more structural piece $\theta\,(D - m(X))$, and in using the residualized $D - m(X)$ rather than the raw $D$ in the denominator. It is precisely this change that makes the moment condition orthogonal to $\ell$ and $m$ at once. Intuitively, $Y$ residualized by $\hat\ell$ purges the confounding on the outcome side, and $D$ residualized by $\hat m$ purges the confounding on the treatment side, cleaning both sides at once, so that a biased nuisance estimate on either side is absorbed at first order by the residual orthogonality of the other. By contrast, Section 3.1 residualized only $Y$, the bias on the treatment side had no offset, and the first-order term leaked out. In Chernozhukov and coauthors' words, this is a nonlinear version of Frisch-Waugh that removes the original's reliance on linearity.

::: {.theorem}
(The debiasing property of DML, informal statement) Let $\hat\theta$ be defined by the sample moment $\frac{1}{n}\sum_i \psi(W_i; \hat\theta, \hat\ell, \hat m) = 0$, with $\psi$ satisfying Neyman orthogonality. If the nuisance estimates satisfy $\|\hat\ell - \ell_0\| \cdot \|\hat m - m_0\| = o_p(n^{-1/2})$ (the product-rate condition, for which each being faster than $n^{-1/4}$ suffices), and cross-fitting is used (see Section 3.3), then

$$\sqrt{n}\,(\hat\theta - \theta_0) \xrightarrow{d} \mathcal{N}(0, \sigma^2),$$

where $\sigma^2$ can be consistently estimated. That is, $\hat\theta$ attains $\sqrt{n}$ convergence and asymptotic normality, and its first-order asymptotic distribution is exactly the same as when the nuisances are known, so machine learning's regularization bias does not enter at first order.
:::

This theorem is the whole crux of DML. It does not say that machine learning estimates accurately (they converge slowly), it says that orthogonalization makes inference on the causal parameter immune to the nuisances' slow rate. The product-rate condition $\|\hat\ell - \ell_0\| \cdot \|\hat m - m_0\| = o_p(n^{-1/2})$ also brings a familiar benefit: it is a kind of double robustness, and as long as the product of the two nuisances' convergence rates is fast enough, either one being a bit worse can be made up for by the other, which is the same spirit as the double robustness of AIPW in Chapter 9, except that there it was about model specification being right or wrong, and here it is about convergence rates being fast or slow.

### 3.3 Cross-Fitting: Cutting the Self-Contamination of Overfitting

Orthogonalization cures regularization bias, but there is a more insidious disease left, called overfitting bias, and the second remedy, cross-fitting (also called sample splitting), is the specific cure for it.

The root of the disease is that if you use the same batch of data both to train the nuisance and to construct the moment condition, then $\hat\ell(X_i)$ has seen the $i$th observation's $Y_i$, and its prediction at $X_i$ secretly contains information about $u_i$. So the residual $Y_i - \hat\ell(X_i)$ is correlated with $\hat\ell$'s estimation error, the two sides that should be orthogonal are entangled again, and a bias term is produced that even Neyman orthogonality cannot remove. It comes from the data's self-reference: one observation serves as both a training sample and an evaluation sample. The 0.170 forest plug-in of Section 1 has exactly this disease: the treatment variable goes into the forest together with the covariates, and each observation's treatment-effect estimate borrows from its own outcome.

The cross-fitting fix is almost brutally simple: keep the data used to train the nuisance completely separate from the data used to construct the moment condition. Concretely, it is $K$-fold (usually $K = 5$):

First split the sample randomly into $K$ parts. For the $k$th part (the evaluation fold), use the remaining $K - 1$ parts (the training folds) to train $\hat\ell^{(-k)}$ and $\hat m^{(-k)}$, then use them to predict the nuisances for each observation in the $k$th part, and construct those observations' residuals and moment contributions. Cycle through all $K$ folds, so that each observation's nuisance is estimated from data that does not include itself, and the self-reference is cut. Finally combine the moment contributions from all $K$ folds and solve for $\hat\theta$.

Because each observation's $\hat\ell(X_i)$ comes from training folds that do not include $i$, $\hat\ell$'s error is independent of $u_i$, and the overfitting bias is zero by construction. The only cost is that each nuisance is trained on a fraction $(K-1)/K$ of the data, and the efficiency loss is tiny, negligible at $K = 5$. To eliminate the randomness introduced by the sample splitting itself, in practice one also repeats the whole cross-fitting several times (with different random splits) and takes the median. For clarity of exposition this chapter does only a single split and fixes the seed.

Put the two remedies together and the DML recipe is complete: a Neyman-orthogonal moment condition, paired with cross-fit nuisance estimates. The former blocks regularization bias (from biased nuisances), the latter blocks overfitting bias (from the entanglement of nuisance and residual). The theorem of Chernozhukov, Chetverikov, Demirer, Duflo, Hansen, Newey, and Robins (2018) proves that after these two steps, inference on $\hat\theta$ is as clean as when the nuisances are known. Section 6 quantifies the effect of each step on the Halcyon data: naive plug-in 0.064 and 0.170, orthogonal but not cross-fit about 0.140, orthogonal plus cross-fit 0.147, together with a confidence interval that covers the truth.

::: {.warning}
DML is not a master key. Its whole guarantee rests on the product-rate condition $\|\hat\ell - \ell_0\| \cdot \|\hat m - m_0\| = o_p(n^{-1/2})$, and this condition fails when the sample is not large enough, or when the nuisance functions are so complex that machine learning cannot learn them. When the sample is small, the forest learns the interaction-plus-kink confounding structure of Section 2 only crudely, both nuisances converge slowly, the product term has not yet been pushed below $n^{-1/2}$, and DML's point estimate still carries a sizable bias while its confidence interval undercovers. This chapter's Monte Carlo (Section 6) shows this honestly: in a moderate sample DML has already cut naive plug-in's bias by more than half, but for its interval to actually reach nominal coverage still requires a larger sample or better-tuned nuisances. "Using DML makes you unbiased" is a dangerous misreading: what it removes is the first-order bias from the two specific sources of regularization and overfitting, and it cannot remove the harm when the nuisances simply cannot be learned.
:::

### 3.4 From ATE to CATE: The Idea of the Causal Forest

DML gives us the average effect. But Halcyon's effect is heterogeneous, and what management really wants to know is whom the tool helps. Estimating the CATE $\tau(x)$ needs an estimator that can vary with $x$, and the causal forest (Wager and Athey 2018; the generalized random forest of Athey, Tibshirani, and Wager 2019) is the most mature answer so far. Its idea can be pushed all the way from the matching of Chapter 9.

To estimate the CATE at $x$, the most naive way is to find a bunch of merchants whose covariates are close to $x$ and compare the outcome difference between adopters and non-adopters among them. This is matching, and its Achilles' heel is that "close" has no good definition in high dimensions: in 20 dimensions the neighbors close to you are pitifully few, and which dimensions matter and how to weight them is not known in advance. The key step of the causal forest is to let the data itself learn how "close" should be defined.

It does this. Grow a stand of trees, but the objective of each tree's splits is not, as in an ordinary regression tree, to make the outcome $Y$ purer, but to maximize the difference in the treatment effect across child nodes, that is, to hunt specifically for the split directions where effect heterogeneity is strongest. Once the trees are grown, to estimate the CATE at $x$, look at which leaf $x$ falls into in each tree, gather all the training samples that land in the same leaf as $x$, and weight them by the frequency with which they share a leaf with $x$. This set of forest weights $\alpha_i(x)$ is a data-adaptive similarity: the neighbors the forest deems most relevant to the effect at $x$ get the highest weight. With the weights in hand, the CATE is a locally weighted treatment-effect estimate, formally like the solution to a local moment condition:

$$\hat\tau(x) = \frac{\sum_i \alpha_i(x)\, (D_i - \bar D_x)(Y_i - \bar Y_x)}{\sum_i \alpha_i(x)\, (D_i - \bar D_x)^2},$$

where $\bar D_x$ and $\bar Y_x$ are forest-weighted local means. You can see it is an adaptive nearest-neighbor version of residual-on-residual, of one lineage with DML, except that the weights are learned from data and the effect is allowed to vary with $x$.

::: {.assumption}
**honesty** Each tree uses one subsample to decide the split structure (where to cut) and another, non-overlapping subsample to fill in the effect estimates in the leaves (how much to estimate after cutting). The two subsamples do not overlap.
:::

Honesty is the key to the causal forest being able to do valid inference, and it deserves to be dissected on its own. An ordinary adaptive tree has a hidden hazard: it uses the data both to find the split points and to estimate the values in the leaves, so the splits chase the noise, cutting the observations that happen to be high into the same leaf, the leaf estimate is therefore overly optimistic, and the confidence interval is distorted. Honesty hands "finding the splits" and "estimating the effect" to two non-overlapping batches of data, so the split structure is exogenous to the within-leaf estimate, and the $\hat\tau$ in the leaf is therefore asymptotically unbiased and follows a normal distribution with computable variance. This is entirely of a piece with the spirit of cross-fitting in Section 3.3: both cut the self-reference through data splitting in exchange for honest inference, except that one is applied to nuisance estimation and the other to the splitting and filling of trees. Athey-Tibshirani-Wager prove that the honest causal forest's $\hat\tau(x)$ is pointwise asymptotically normal, so a confidence interval can be attached to each $x$.

The causal forest also comes with a by-product for estimating the ATE, and one that fits the estimand we want better than the PLR. Use the forest-estimated $\hat\tau(x)$, $\hat m(x)$, $\hat\ell(x)$ to construct the AIPW (augmented inverse propensity weighting) doubly robust score of Chapter 9, then average, and you get an estimate of the ATE. This AIPW version estimates the unweighted $\mathbb{E}[\tau(X)]$, the true ATE, rather than the PLR's overlap-weighted version. In Section 6 the causal forest's AIPW ATE is 0.148 (SE 0.007), and DML's IRM, taking the same AIPW route, gives 0.159, both aiming at the SATE of 0.137.

### 3.5 The Credibility of Heterogeneity: Calibration, GATES, and CLAN

The forest can spit out a $\hat\tau(x)$ curve, but a good-looking curve does not mean the heterogeneity is real. What machine learning is best at is fitting patterns that do not exist out of nothing, so after estimating the CATE, we must answer two questions: whether the curve is well calibrated overall, and whether the high- and low-effect groups it sorts out are really high and low. This subsection introduces three companion diagnostic tools that turn heterogeneity from "the forest says so" into "the data can verify it."

The first is the calibration test (Chernozhukov and coauthors). It regresses the outcome on two things: the forest's predicted average effect, and each observation's differential effect deviating from the average, $\hat\tau(X_i) - \bar{\hat\tau}$. Ideally the coefficient on the former (mean.forest.prediction) should be 1, showing the overall effect level is estimated right, and the coefficient on the latter (differential.forest.prediction) should also be 1, showing the forest captured both the direction and the magnitude of the heterogeneity. If the differential coefficient is significantly greater than zero, there is statistical evidence that the heterogeneity is real; if it is not significant, the ups and downs the forest draws are very likely noise. In Section 6 these two coefficients for Halcyon are 1.004 and 1.199, the former saying the level is almost perfectly calibrated, the latter significantly positive saying heterogeneity truly exists.

The second is GATES (group average treatment effects by heterogeneity score). Sort the merchants from low to high by the estimated $\hat\tau(x)$, cut them into several groups (say quintiles), and for each group compute a doubly robust group-average effect using AIPW scores. If the forest's ranking is meaningful, these group averages should increase monotonically, with the lowest group and the highest group pulling clearly apart. The beauty of GATES is that the definition of the groups uses the forest's $\hat\tau$, but the estimation of the group effects uses independent AIPW scores, and the two are not from the same source, so it is an out-of-sample validation of heterogeneity, not self-confirmation.

The third is CLAN (classification analysis). It asks the reverse question: what do the two groups of people the forest deems highest- and lowest-effect look like on the raw covariates. List and compare the covariate means of the highest-effect group and the lowest-effect group, and you can give heterogeneity an interpretable portrait. In Section 6 Halcyon's CLAN shows that the highest-effect group of merchants has a size ($X_1$) mean of $-1.31$ and a tenure ($X_2$) of $-1.06$, both clearly below the lowest-effect group's $1.10$ and $0.97$, which in plain language means: small and new merchants benefit most from the tool, while large and old merchants benefit least or even are hurt. This fits our mechanistic intuition about Growth Playbook: an automated business advisor is a lifeline to inexperienced small merchants, and superfluous to large merchants who have long since figured out the ropes. Only here does the heterogeneity stand up: endorsed by the calibration test, given the monotone staircase of GATES, and given an interpretable portrait by CLAN.

### 3.6 From CATE to Policy: Policy Learning

Estimating heterogeneity is only a means; the end is a decision: the tool is not free, and the platform should decide whom to give it to. This subsection turns the CATE into a deployable policy, and gives the language for judging whether a policy is good.

Let the cost of switching the tool on for one merchant be $c$ (converted into units of log GMV growth). A policy $\pi: \mathcal{X} \to \{0, 1\}$ is a rule that looks at a merchant's covariates $x$ and decides whether to switch it on. The value of a policy is the net growth it brings:

$$V(\pi) = \mathbb{E}\big[\pi(X)\,(\tau(X) - c)\big] + \mathbb{E}[Y(0)],$$

The second term has nothing to do with the policy (the baseline everyone gets whether or not switched on), so comparing policies looks only at the first term, the sum of net effects over the switched-on merchants. The optimal policy is obvious: switch on for, and only for, the merchants with $\tau(x) > c$, because only when the net effect is positive is it worth the investment. The problem is that $\tau(x)$ is unknown, we have only the forest-estimated $\hat\tau(x)$, and the policy $\hat\pi$ learned from it is necessarily imperfect. The quantity measuring how far it is from optimal is called regret:

$$\text{Regret}(\hat\pi) = V(\pi^\star) - V(\hat\pi) \geq 0,$$

where $\pi^\star$ is the oracle policy that knows the true $\tau$. The theoretical goal of policy learning is to prove that the learned policy's regret goes to zero at some rate.

There are two roads to learning a policy. One is plug-in: directly use the forest's $\hat\tau(x) > c$ as the rule. It is simple, but has two flaws: first, the noise in $\hat\tau$ near the decision boundary makes the rule flip back and forth, and second, the resulting rule is often a complex black box, impossible to explain to the business team and hard to deploy. The other road is to directly optimize a restricted, interpretable policy class, the approach of Athey and Wager (2021): rather than trusting $\hat\tau$ pointwise, use the AIPW doubly robust scores $\hat\Gamma_i$ of Chapter 9 (each an unbiased but noisy estimate of a merchant's treatment effect) as rewards, and pick out the one with the highest empirical value within a simple policy class. A common policy class is a shallow decision tree (policy tree), for instance a depth-2 tree, which uses only two or three covariates and a few thresholds, so the business team understands it at a glance and can write it straight into the targeting system. Athey-Wager prove that optimizing the doubly robust scores within such a restricted class, the learned policy's regret converges to the best in that class at nearly $n^{-1/2}$.

Section 6 learns a depth-2 policy tree on the Halcyon data. The rule it gives can be stated in one sentence: switch the tool on first for small merchants with short tenure, and for merchants with long tenure only when they are also small in size. Under this governance about 53% of merchants are switched on, delivering a net value of 0.073 per merchant (relative to baseline), whereas switching everyone on indiscriminately delivers only 0.037, and doing nothing is 0. More importantly, it approaches the oracle: the optimal rule that knows the true $\tau$ has a net value of 0.075, and the learned tree captures 97.5% of it, making the same decision as the oracle on 91.2% of merchants. Targeting nearly doubles the governance value, and this is the cash value of heterogeneity.

The points of this section, in brief: plugging machine-learning predictions straight into a causal estimate introduces regularization bias (first-order, fatal) and overfitting bias; Neyman orthogonality makes the moment condition immune to first-order error in the nuisances, cross-fitting cuts the self-reference, and the two together form DML, making causal inference immune to slow-rate machine learning, at the cost of needing the product of the nuisances' convergence rates to be faster than $n^{-1/2}$, so the guarantee is discounted when the sample is insufficient; applying the same splitting idea to trees (honesty) yields the causal forest, which estimates CATE pointwise with honest intervals, and the heterogeneity still has to stand up under the triple validation of calibration, GATES, and CLAN; finally policy learning turns the CATE into an interpretable targeting rule with a regret guarantee.

## 4 Estimation

The logic of identification and debiasing is clear, and this section covers implementation: how each estimator is computed concretely, how it is done in R, and in what settings they substitute for one another. The spine is still "where the previous method falls short gives rise to the next."

### 4.1 The Partially Linear Implementation of DML

The most basic DML is the cross-fit orthogonal implementation of the PLR of Section 3, in three steps. First, split the sample into $K$ folds, train two nuisance learners outside each fold, one regressing $Y$ on $X$ to get $\hat\ell$ and one regressing $D$ on $X$ to get $\hat m$, and use them to predict the fold's nuisances. Second, on the full sample construct the residuals $\tilde Y_i = Y_i - \hat\ell(X_i)$ and $\tilde D_i = D_i - \hat m(X_i)$. Third, solve the orthogonal moment condition, equivalent to running $\tilde Y$ on $\tilde D$ in a no-intercept regression, $\hat\theta = \sum_i \tilde D_i \tilde Y_i / \sum_i \tilde D_i^2$, with the standard error given by the sample variance of the influence function $\psi_i = \tilde D_i(\tilde Y_i - \hat\theta \tilde D_i) / \mathbb{E}[\tilde D^2]$. The nuisance learners can be random forests, gradient boosting, lasso, or ensembles of them, as long as each converges fast enough. The DoubleML package in R wraps up these three steps, and also allows you to write them by hand; Section 6 runs both for the reader, and the numbers are nearly identical.

### 4.2 IRM and AIPW: Aiming Directly at the ATE

The PLR's $\hat\theta$ estimates the overlap-weighted effect, and to get the unweighted ATE, you have to switch to the interactive regression model (IRM). It does not assume the treatment effect is constant, but models the two treatment states separately, using the AIPW doubly robust score of Chapter 9:

$$\hat\tau_0^{AIPW} = \frac{1}{n}\sum_i \Bigg[ \hat\mu_1(X_i) - \hat\mu_0(X_i) + \frac{D_i\,(Y_i - \hat\mu_1(X_i))}{\hat e(X_i)} - \frac{(1 - D_i)(Y_i - \hat\mu_0(X_i))}{1 - \hat e(X_i)} \Bigg],$$

where $\hat\mu_d(X) = \mathbb{E}[Y \mid X, D = d]$ are the two outcome regressions and $\hat e(X)$ is the propensity. This score is orthogonal to both $\mu_d$ and $e$ (which is exactly where the double robustness of Chapter 9 comes from), and paired with cross-fitting it is a DML estimator aimed at the ATE. Its difference from the PLR is that it splits the outcome regression into two, one for the treated group and one for the control, so it no longer mixes overlap weights into the estimand. The cost is that it explicitly puts the propensity in the denominator, so when overlap is poor ($\hat e$ close to 0 or 1) the variance blows up, a point Section 7 details. In Section 6 the IRM gives 0.159 (SE 0.009), slightly higher and slightly more variable than the PLR's 0.149: the PLR aims at the overlap-weighted effect (truth 0.144), the IRM aims at the ATE (truth 0.137), the difference coming from weighting or not, and secondarily from the extra noise that inverting the propensity adds to the IRM.

### 4.3 Causal Forest

To estimate the CATE, use a causal forest. In R the grf package's causal_forest grows one in a single line; internally it automatically uses honesty, automatically orthogonalizes $Y$ and $D$ once each (that is, first subtracting $\hat\ell$ and $\hat m$, of the same origin as DML), and then splits by effect heterogeneity. Once grown, predict gives each merchant's $\hat\tau(x)$ (by default out-of-bag predictions, a further guard against overfitting), average_treatment_effect gives the AIPW ATE, and test_calibration gives the calibration test of Section 3.5. The forest's tuning (number of trees, minimum leaf size, number of variables drawn at each split) is handed to cross-validation via the built-in tune.parameters, with no need to tune by hand. In Section 6 the forest-estimated $\hat\tau(x)$ has a correlation of 0.938 with the truth, showing it learned the ranking of the heterogeneity quite accurately.

### 4.4 GATES, Best Linear Projection, and RATE

Once heterogeneity is estimated, Section 4.4 covers how to compress it into a reportable, testable object. GATES was already covered in Section 3.5, and in implementation it is just grouping by $\hat\tau$ and computing an AIPW group mean for each group. The best linear projection (BLP) is another compression: project $\hat\tau(x)$ onto a few covariates of interest, to get a concise summary of "how the effect varies linearly with size and tenure," and grf's best_linear_projection gives the coefficients and robust standard errors directly. In Section 6 projecting the CATE onto size $X_1$ and tenure $X_2$ gives slopes of $-0.080$ and $-0.065$, consistent with the true $-0.101$ and $-0.061$ in the data-generating process, proving that the forest caught not only the existence of heterogeneity but also its direction and rough magnitude.

There is also a tool devoted to answering "is targeting actually worth it," RATE (rank-weighted average treatment effect). Along the path of "treating more people in order of $\hat\tau$ from high to low," it measures how much more you earn by prioritizing high-effect people over treating at random, and outputs a number called AUTOC or the area under the TOC curve, together with a standard error. If RATE is significantly positive, the forest's ranking really can be used to target for profit; if it is not significant, then however much heterogeneity "exists," it has no decision value. In Section 6 Halcyon's RATE is 0.134 (SE 0.010), significantly positive, so targeting is profitable, which paves the way for the policy learning of Section 4.5.

### 4.5 Policy Tree

Finally, turn heterogeneity into a rule. The policytree package's policy_tree takes two things, the covariates and a doubly robust reward matrix (each merchant's AIPW score under the two actions "switch on" and "do not switch on," with the switch-on column minus the cost $c$), and among all decision trees of a specified depth (say 2) it exhaustively finds the one with the highest empirical net value. It returns a shallow tree that can be drawn and written into the targeting system. Depth is a tradeoff between interpretability and approximation power: depth 1 can only cut on a single threshold, depth 2 can express interaction rules like "small and new," and any deeper is hard to interpret and prone to overfitting. The depth-2 tree learned in Section 6 uses only the two variables size and tenure and three thresholds, and captures 97.5% of the oracle policy's value.

The estimation ladder of this section, in brief: DML's partially linear implementation gives an overlap-weighted effect, the IRM switches to AIPW scores to aim directly at the ATE, and both rely on orthogonalization plus cross-fitting to be immune to machine learning's bias; for heterogeneity, use a causal forest, which pairs the same orthogonalization idea with honesty to estimate CATE pointwise; GATES, the best linear projection, and RATE compress this curve into an object that is reportable, testable, and able to judge the value of targeting; finally the policy tree lands heterogeneity into a deployable shallow rule. Each rung takes over where the previous one is not enough: when a constant effect is not enough, estimate a curve; when the curve itself is hard to believe, do the triple validation; only after validation passes do you talk about targeting.

## 5 Anchoring Papers

A method stands up only when it lands in real research. Three anchoring papers, one laying the theory of DML, one laying the causal forest, and one applying this toolkit to heterogeneous effects in a platform setting, each organized by the five elements of paper, method, data, results, and limitations, with the focus on how the assumptions are defended.

### 5.1 Chernozhukov et al. (2018)

::: {.case}
Paper and theoretical position: "Double/Debiased Machine Learning for Treatment and Structural Parameters," The Econometrics Journal. This paper integrates the ideas of orthogonalization and sample splitting, scattered across the semiparametric estimation literature, into a unified framework, providing a rigorous asymptotic theory for "estimating nuisances with machine learning and then doing causal inference," and is the original source of all the results in Section 3 of this chapter.

Method: the core is two things, a Neyman-orthogonal moment condition plus cross-fitting. The paper proves that as long as the nuisance learners satisfy the product-rate condition (each faster than $n^{-1/4}$), orthogonalization keeps machine learning's regularization bias out of first order, cross-fitting keeps overfitting bias out, and the target parameter attains $\sqrt{n}$ convergence and asymptotic normality, with asymptotic variance the same as when the nuisances are known. The framework applies to a large class of parameters: the PLR's treatment effect, the IRM's ATE, partially linear IV, LATE, and so on, of which this chapter uses only the first two.

Data: the paper gives an empirical demonstration on the effect of 401(k) eligibility on wealth, with the treatment being whether the employer offers 401(k) eligibility, the confounders a set of household financial characteristics, and the nuisances estimated separately by random forests, boosting, neural networks, and lasso, showing that different learners give stable and consistent effect estimates.

Results: the DML estimates are robust to the choice of nuisance learner, with different machine-learning methods giving effect estimates close to one another, whereas naive plug-in swings wildly with the learner. This robustness is itself an empirical footnote to the theory: orthogonalization makes the conclusion no longer depend on how accurately the nuisances are estimated.

Limitations: all guarantees rest on the product-rate condition, and when the sample is insufficient or the nuisance functions too complex, it fails and DML is no longer unbiased. The paper also concedes that DML handles the case where the nuisances are high-dimensional but the target parameter is low-dimensional, and problems where the target itself is high-dimensional (for instance the whole CATE curve) need extra structure of the causal-forest kind.
:::

This paper's defense strategy is one of pure theoretical honesty: it does not claim machine learning can estimate the nuisances accurately, but precisely delineates under how weak a condition, and at what rate, inference on the causal parameter can be made immune to the nuisances' inaccuracy, and lays the boundary of failure (the product-rate condition) plainly on the table.

### 5.2 Wager and Athey (2018)

::: {.case}
Paper: "Estimation and Inference of Heterogeneous Treatment Effects using Random Forests," Journal of the American Statistical Association. It reworks the random forest from a prediction tool into a causal estimator that can do valid statistical inference on the CATE, and is the founding work of the causal forest; the honesty and pointwise normality of Section 3.4 come from here and its follow-up, the generalized random forest (Athey, Tibshirani, and Wager 2019).

Method: there are two core innovations. The first is the honest tree: use one subsample to decide the splits and another, non-overlapping subsample to fill in the leaf estimates, so that the within-leaf treatment-effect estimate is exogenous to the split structure and can be asymptotically unbiased. The second is reinterpreting the forest as an adaptive kernel, with the forest weights $\alpha_i(x)$ giving a data-learned similarity, and the CATE being the local treatment effect under those weights. The paper proves that the honest causal forest's $\hat\tau(x)$ is pointwise asymptotically normal and can be paired with a consistently estimated variance, so a confidence interval can be given for each $x$, which earlier heterogeneous-effect methods could not do.

Data: the paper validates coverage and convergence on simulated systems, and demonstrates estimating CATE on observational data; the subsequent literature has widely used it for heterogeneous effects in medicine, labor, platforms, and other fields.

Results: the honest causal forest's confidence intervals reach nominal coverage, and the CATE estimates adaptively become fine in regions of strong heterogeneity and automatically smooth in weak regions. The slight increase in variance that honesty brings is far outweighed by the gain of valid inference.

Limitations: honesty sacrifices half the data to splitting and half to estimation, so precision suffers in small samples; the forest can also mistake noise for heterogeneity when the covariate dimension is very high and the heterogeneity signal very weak. Pointwise inference is valid, but joint inference on "the whole curve" still calls for caution.
:::

This paper's contribution is stitching machine learning's predictive power to statistics' inferential rigor, and this honesty step comes from the same insight as DML's cross-fitting: to get honest inference, you have to pay the cost of data splitting for "finding structure" and "estimating values."

### 5.3 Heterogeneous Treatment Effects in Platform Settings

::: {.case}
Paper and direction: the representative direction for applying this toolkit in platform settings is estimating the heterogeneous effect of an intervention (price, recommendation, subsidy, tool) across users or merchants, and targeting on that basis. A repeatedly cited methodological demonstration is Athey and Wager (2021, "Policy Learning with Observational Data," Econometrica), which connects the causal forest's CATE to policy learning, giving the full chain from observational data to a targeting rule with a regret guarantee, exactly the skeleton of Sections 3.6 and 4.5 of this chapter.

Method: first use a doubly robust score to estimate each individual's treatment effect unbiasedly (and noisily), then maximize the empirical policy value within a restricted, interpretable policy class (such as a shallow decision tree), and prove that the learned policy's regret relative to the best in that class converges at nearly $n^{-1/2}$. The key is that it does not require trusting the CATE pointwise, only that the doubly robust score be unbiased, which lets policy learning inherit DML's robustness.

Data: the method can be used for any observational or experimental data with treatment, outcome, and covariates; in platform settings a typical application is targeting promotions, features, or subsidies to the users or merchants whose CATE exceeds the cost.

Results: the targeting rule within the restricted policy class is both interpretable and near-optimal, converting heterogeneity into a sizable net-value gain, without deploying a black box that trusts the CATE pointwise.

Limitations: everything is built on unconfoundedness and overlap, and in observational data both can fail at any time; choosing the policy class too simple limits approximation power, while too complex makes it hard to interpret and prone to overfitting, so depth is a design choice to be traded off. The regret guarantee is relative to the best within the policy class, not relative to an omniscient oracle.
:::

Taken together, the meaning of the anchoring becomes clear: Chernozhukov and coauthors give the theoretical license for doing causal inference with machine learning, Wager and Athey extend it to the whole CATE curve while preserving inference, and Athey and Wager then connect heterogeneity to the decision. They share one methodological core: to make machine learning serve causality, you have to use orthogonalization and data splitting to reconcile "fitting well" and "inferring correctly," two things that would otherwise fight each other.

## 6 A Full Walkthrough on the Halcyon Data

Now run the earlier tools through once, in full, on the Halcyon panel. The code below uses R 4.5.3, fixing set.seed(1313) for reproducibility, and every number cited in the text comes from the actual run output of this code. The data are 10000 merchants and 20 covariates, adoption is driven by a nonlinear propensity, the true SATE is 0.137, and the effect decreases with size and tenure.

### 6.1 DGP

The design parameters are as follows: a persistent merchant quality is spread over 20 correlated covariates, and the adoption propensity and baseline outcome are driven by the same set of nonlinear terms (a size-by-momentum interaction and a quadratic kink in tenure). It is precisely this shared nonlinear structure that makes linear control fail while giving machine learning something to work with.

```r
set.seed(1313)
n <- 10000L; p <- 20L
Sigma <- outer(1:p, 1:p, function(i, j) 0.5^abs(i - j))
X <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)
X1 <- X[, 1]; X2 <- X[, 2]; X3 <- X[, 3]      # size, tenure, momentum

h_conf   <- 0.7 * X1 * X3 + 0.6 * (X2^2 - 1)  # shared nonlinear confounding
ps_index <- 0.40 * X1 + 0.32 * X3 + 0.50 * h_conf
e_true   <- pmin(pmax(plogis(ps_index), 0.05), 0.95)
D  <- rbinom(n, 1, e_true)
g0 <- 0.20 * X3 + 0.15 * X1 + 0.50 * h_conf + 0.10 * X[, 6]
tau_i <- 0.12 - 0.10 * X1 - 0.06 * X2 + 0.03 * X1 * X2   # CATE
Y  <- g0 + tau_i * D + rnorm(n, 0, 0.30)
```

The effect $\tau(x) = 0.12 - 0.10\,x_1 - 0.06\,x_2 + 0.03\,x_1 x_2$ decreases with size and tenure, small and new merchants benefit most, and a few large and old merchants have negative effects. The truths computed from the realized effects: the SATE is 0.137 and the ATT is 0.118. Note the ATT is below the SATE, because the platform gave the tool precisely to the batch of large merchants with small benefits, a mirror image of selection on gains, arising from the same logic as LATE below ATE in Chapter 12, that "the ones selected are not necessarily the ones who benefit most." The propensity was clipped to $[0.05, 0.95]$ at generation and its realized values land in $[0.10, 0.95]$, so overlap is ample.

### 6.2 Naive Estimates: Reproducing the Four Wrong Numbers

```r
dim <- mean(Y[D == 1]) - mean(Y[D == 0])                 # 0.517
ols <- coef(lm(Y ~ D + X))["D"]                          # 0.450
rf  <- ranger(Y ~ ., data.frame(Y, D, X), mtry = 8)      # plug-in
plugin <- mean(predict(rf, data.frame(D = 1, X))$pred -
               predict(rf, data.frame(D = 0, X))$pred)   # 0.170
```

Four numbers reproduce the opening: naive difference 0.517 (SE 0.015), linear control 0.450 (SE 0.014), forest plug-in 0.170, and the naive FWL of first residualizing $Y$ and then regressing $D$, 0.064. All deviate from the truth of 0.137, and the directions are not even consistent. Linear control does almost nothing, because the confounding hides in the $X_1 X_3$ interaction and the $X_2^2$ kink, which linear regression cannot absorb; the forest plug-in and the naive FWL each carry regularization bias, one overstating and one understating.

### 6.3 DML: Orthogonalization plus Cross-Fitting

```r
# 5-fold cross-fit, ranger nuisances for l(X)=E[Y|X], m(X)=E[D|X]
for (k in 1:5) {
  tr <- fold != k; te <- fold == k
  lhat[te] <- rf_fit(Y[tr], X[tr, ], X[te, ])
  mhat[te] <- rf_fit(D[tr], X[tr, ], X[te, ])
}
Yt <- Y - lhat; Dt <- D - mhat
theta_cf <- sum(Dt * Yt) / sum(Dt * Dt)                  # 0.147
```

First look at the separate contributions of orthogonalization and cross-fitting. Residualizing both sides but estimating the nuisances on the full sample (orthogonal, not cross-fit) gives 0.140; adding 5-fold cross-fitting gives 0.147 (SE 0.007). Both have already pulled that pile of wrong numbers from the opening back near the truth, showing that orthogonalization is the main workhorse here (pulling 0.064 and 0.170 to around 0.14 at a stroke) and cross-fitting is the finer second-order correction. The DoubleML package's PLR gives 0.149 (SE 0.007), almost identical to the hand-written implementation; its IRM (AIPW, aiming directly at the ATE) gives 0.159 (SE 0.009). Recall Section 2: the PLR estimates the overlap-weighted effect (truth 0.144) and the IRM estimates the ATE (truth 0.137), each estimator aiming at its own target, and the numerical difference comes mainly from the differing estimand, and secondarily from the variance that inverting the propensity adds to the IRM.

A single estimate is only one point. To see whether DML systematically debiases, run the whole procedure once each on 300 independent draws (4000 merchants each time, with a lighter forest to make repetition feasible), and compare the sampling distributions of the naive forest plug-in and the cross-fit DML.

![Distributions of the estimates from naive machine-learning plug-in (orange) and cross-fit DML (blue) over 300 draws, with the dashed line the true SATE. The naive plug-in systematically deviates from the truth, while DML's distribution is clearly closer to the truth and more concentrated.](assets/fig/fig_16_mc.svg)

The two distributions confirm the reasoning of Section 3. The naive plug-in is systematically too high, with a mean bias of 0.055; cross-fit DML cuts the mean bias to 0.034 and the distribution is more concentrated (standard deviation 0.011 versus 0.013), so the debiasing does work. But to be honest about it, on this moderate sample of only 4000 each time, DML's bias has not yet gone fully to zero, and its 95% confidence interval covers the truth only 0.16 of the time, far short of the nominal 0.95. This is exactly the discounting of the product-rate condition in insufficient samples that the warning box of Section 3.3 flags, made visible: the forest learns the interaction-plus-kink confounding of Section 2 not fast enough, and the second-order term has not been pushed below $n^{-1/2}$. The single estimate in the first half of this section uses the full sample of 10000, and only then does the point estimate land back near the truth (0.147 versus 0.137); for the interval to reach nominal coverage too, one needs a larger sample or faster-converging nuisance learners. Debiasing is not a guarantee of removing bias, it lowers the order of the bias, and the rest still has to be delivered by sample size.

### 6.4 Causal Forest: Estimating CATE and Validating Heterogeneity

```r
cf <- causal_forest(X, Y, D, num.trees = 4000, honesty = TRUE,
                    tune.parameters = "all")
tau_hat <- predict(cf)$predictions              # out-of-bag CATE
ate     <- average_treatment_effect(cf)         # AIPW ATE
cal     <- test_calibration(cf)                 # calibration test
```

The forest's AIPW ATE is 0.148 (SE 0.007), covering the truth. The estimated $\hat\tau(x)$ has a correlation of 0.938 with the true CATE, so the forest learned the ranking of the heterogeneity quite accurately. The calibration test gives a mean.forest.prediction coefficient of 1.004 and a differential.forest.prediction coefficient of 1.199 ($t = 16.8$), the former saying the overall level is almost perfectly calibrated, the latter significantly greater than zero saying heterogeneity truly exists and its direction is estimated right. Projecting the CATE onto size $X_1$ and tenure $X_2$ gives slopes of $-0.080$ and $-0.065$, consistent with the true generating coefficients $-0.101$ and $-0.061$.

GATES splits merchants into five groups by estimated CATE and computes an AIPW group effect for each group:

| Quintile (by estimated CATE) | AIPW group effect (SE) | True group mean |
|---|---|---|
| Q1 (lowest) | 0.025 (0.016) | -0.014 |
| Q2 | 0.050 (0.016) | 0.041 |
| Q3 | 0.090 (0.015) | 0.109 |
| Q4 | 0.194 (0.015) | 0.197 |
| Q5 (highest) | 0.381 (0.017) | 0.353 |

The group effects increase monotonically, climbing from nearly zero in the lowest group to 0.381 in the highest, each cell hugging the true group mean. The RATE (AUTOC) is 0.134 (SE 0.010), significantly positive, confirming that targeting by CATE really is profitable. CLAN paints the portrait of the heterogeneity: the highest-effect group has a size mean of $-1.31$ and tenure of $-1.06$, and the lowest-effect group has size $1.10$ and tenure $0.97$, which translated means small and new merchants benefit most and large and old merchants least, fitting the mechanism of Growth Playbook.

![Five groups formed by estimated CATE, with AIPW group effects (blue, with 95% confidence intervals) and true group means (orange triangles). The forest's ranking spreads merchants from low effect to high effect, and the group effects increase monotonically and hug the truth.](assets/fig/fig_16_cate.svg)

### 6.5 Policy Learning: Turning Heterogeneity into Targeting

```r
Gamma <- double_robust_scores(cf)               # AIPW reward matrix
Gamma[, 2] <- Gamma[, 2] - 0.10                 # subtract cost c = 0.10
ptree <- policy_tree(X[, c("X1", "X2")], Gamma, depth = 2)
```

Set the switch-on cost per merchant at $c = 0.10$. The learned depth-2 policy tree uses only the two variables size and tenure: for merchants with short tenure ($X_2 < -0.23$), switch on unless size is very large ($X_1 \geq 0.73$); for merchants with longer tenure, switch on only when size is also small ($X_1 < -0.26$). In one sentence, switch the tool on first for small and new merchants. Evaluating each policy's net value gain by the true $\tau$ (relative to switching no one on):

| Policy | Share switched on | Net value gain |
|---|---|---|
| Switch on no one | 0% | 0 |
| Switch on everyone | 100% | 0.037 |
| policy tree | 53% | 0.073 |
| oracle (knows true $\tau$) | 54% | 0.075 |

The learned tree governs only 53% of merchants yet captures a net gain of 0.073, nearly double the indiscriminate switch-on-everyone (0.037), reaching 97.5% of the oracle optimum (0.075), and making the same decision as the oracle on 91.2% of merchants. This is where the cash value of heterogeneity shows: identifying who benefits and targeting the tool precisely at them nearly doubles the value over a blanket policy.

![Along the path of "treating more merchants in order of estimated CATE from high to low," the net value gain as a function of the share treated (blue line). The curve peaks at about half, switching on everyone (orange dot) falls back, and the policy tree's operating point (green square) lands near the peak.](assets/fig/fig_16_policy.svg)

The walkthrough of this section, in brief: on the same data, with all identification assumptions holding, the naive procedures hand back the four mutually contradictory wrong numbers 0.517, 0.450, 0.170, 0.064, DML with orthogonalization plus cross-fitting hands back 0.147 to 0.159 and covers the truth; the causal forest learns the ranking of the CATE to a correlation of 0.938 with the truth, and the triple diagnostics confirm the heterogeneity is real; policy learning then nearly doubles the governance value and approaches the oracle. Machine learning did not replace causal inference; it is orthogonalization and data splitting that qualify it to serve causal inference.

## 7 Failure Modes and Robustness

In the simulation the identification assumptions are constructed to hold, but in real research they can fail at any time. This section sorts through the most common ways they fail and the actionable responses.

The most fundamental threat is unconfoundedness itself. Both DML and the causal forest take it as a given premise, and the whole apparatus of orthogonalization, cross-fitting, and honesty operates on top of this premise. They can fix regularization bias and overfitting bias, but they are powerless against omitted confounding. If Halcyon's adoption also depends on some factor not in $X$, say the ambition of the merchant's founder, and this factor directly affects growth, then even the most exquisite machine learning is only estimating a wrong estimand very precisely. Here we must flag a dangerous misconception: tossing dozens or hundreds of covariates into a forest does not guarantee unconfoundedness any more than tossing them into a linear regression. High dimension only makes the functional form more flexible, it does not conjure up a missing confounder out of thin air. The credibility of identification still depends only on the institutional background and the treatment-assignment mechanism, and has nothing to do with how complicated the machine learning is. The actionable response is the old road of Chapter 12: pin down exogenous variation to do IV, or use sensitivity analysis (such as Rosenbaum bounds, or the omitted-variable sensitivity framework tailored to DML) to quantify "how strong an omitted confounder would have to be to overturn the conclusion."

The second threat is insufficient overlap. Unconfoundedness guarantees identification, overlap guarantees estimability. When some class of merchants almost all adopt ($\hat e \to 1$) or almost all do not ($\hat e \to 0$), the denominator $\hat e$ or $1 - \hat e$ of the IRM's AIPW score approaches zero, the weights blow up, the variance runs out of control, and a few observations with extreme propensity can dominate the whole estimate. The diagnosis is direct: plot the distribution of $\hat e$ and see whether there is mass piling up near 0 or 1. There are several responses: trimming (dropping observations with extreme propensity, but this changes the estimand, so when reporting you must state that you are estimating the effect within the overlap region); switching to an overlap-weighted estimand (the PLR is naturally this, giving extreme propensities low weight and being more stable as a result); or simply conceding that in a region without overlap the effect is not identified, and confining the conclusion to merchants with support. Halcyon's propensity was clipped to $[0.05, 0.95]$ at generation and its realized values land in $[0.10, 0.95]$, so overlap is ample, but in real platform data overlap is often a hard constraint.

The third threat is the product-rate condition not being satisfied, a failure specific to DML that the warning box of Section 3.3 already pointed out. When the sample is not large enough, or the nuisance functions are so complex that machine learning cannot learn them, $\|\hat\ell - \ell_0\| \cdot \|\hat m - m_0\|$ has not fallen below $n^{-1/2}$, and beyond the first-order bias that orthogonalization blocks, the second-order term is still sizable, so DML's point estimate is biased and its confidence interval undercovers. This chapter's Monte Carlo shows this face on a moderate sample: DML has already cut most of naive plug-in's bias, but for the interval to actually reach nominal coverage still needs a larger sample or better-tuned, faster-converging nuisance learners. The response is pragmatic: use ensemble learning for the nuisances (stack forests, boosting, and lasso, taking the strengths of each to approach faster convergence), tune carefully by cross-validation, prioritize estimating the propensity accurately when the sample is tight (it goes in the denominator and is most sensitive), and keep reservations about small-sample coverage in the report.

The fourth class of problems surrounds the truth or falsity of heterogeneity. Machine learning can draw ups and downs that do not exist, so any CATE plot must first pass the triple validation of Section 3.5: is calibration's differential coefficient significant, is GATES monotone and do the high and low groups really pull apart, is RATE significantly positive. When none of the three passes, the honest conclusion is "there is no detectable heterogeneity, just report an ATE," rather than forcing noise to be read as "small merchants benefit more." Conversely, all three passing does not mean every detail the forest draws is credible: the pointwise CATE has large variance in covariate regions where data are sparse, so passing judgment on any specific merchant on that basis calls for caution, and the sound use is to interpret at the group level as in GATES, or on low-dimensional projections as in BLP.

The fifth is over-promising in policy learning. The policy tree's regret guarantee is relative to the best in the restricted class it lives in, not relative to an omniscient oracle. Choose the policy class too simple (depth 1) and the approximation power is limited, so the learned rule may systematically miss some people who should be treated; choose it too complex (too deep) and it is hard to interpret, prone to overfitting, and the constant in the regret worsens. Moreover, the learned policy faces a changing world at deployment: merchants react to the targeting itself (excluded merchants may churn), and the cost and effect of the intervention drift over time, so a rule optimal today may not be optimal tomorrow. Treating the policy tree as a decision aid that needs continual re-estimation, rather than a set-and-forget automaton, is the only sound attitude.

Stringing these failure modes together, the credibility of causal machine learning ultimately rests not on how advanced the algorithm is, but on two old things: whether the identification assumptions (unconfoundedness, overlap) are credible, and whether machine learning really learned the nuisances to a fast enough convergence. Orthogonalization and cross-fitting are remarkable tools, and they make causal inference immune to machine learning's regularization, but they cannot buy identification, nor make up for a nuisance that cannot be learned. No technique can substitute for substantive judgment about the treatment-assignment mechanism, a point that has not changed from Chapter 9 to this chapter.

## 8 Further Reading

::: {.readings}
Required reading, in suggested order:

- Chernozhukov, Chetverikov, Demirer, Duflo, Hansen, Newey and Robins (2018, The Econometrics Journal). The original framework of DML, the full version of Section 3 of this chapter; focus on the two sections on Neyman orthogonality and cross-fitting and on the product-rate condition.
- Wager and Athey (2018, Journal of the American Statistical Association). The founding work of the causal forest; focus on how honesty buys pointwise asymptotic normality.
- Athey, Tibshirani and Wager (2019, Annals of Statistics). The generalized random forest, placing the causal forest in a unified framework of local moment estimation, the theoretical basis of the grf package.
- Athey and Wager (2021, Econometrica). The regret theory of policy learning, the original text of Sections 3.6 and 4.5 of this chapter.
- Wager (2024, lecture notes/textbook Causal Inference: A Statistical Learning Approach). A modern synthesis integrating the above results into a course, good as an overview first read or a review.

Further reading:

- Robins and Rotnitzky (1995) and Bang and Robins (2005). The source of doubly robust and AIPW, the forerunner of DML's orthogonality, for understanding where the AIPW score of Section 4.2 comes from.
- Chernozhukov, Demirer, Duflo and Fernández-Val (2025, generic ML). The original formulation of GATES and CLAN, a general framework for inference on heterogeneity.
- Kennedy (2023, Electronic Journal of Statistics). The DR-learner, a modern method estimating CATE directly with doubly robust scores, complementary to the causal forest.
- Nie and Wager (2021, Biometrika). The R-learner, extending Robinson residualization to CATE estimation, connecting DML and heterogeneous effects.
- Chernozhukov, Cinelli, Newey, Sharma and Syrgkanis (2026). Omitted-variable sensitivity analysis for causal machine learning, a robustness tool for when the unconfoundedness of Section 7 fails.
- Hitsch, Misra and Zhang (2024, Quantitative Marketing and Economics). An application of heterogeneous treatment effects and targeting in marketing, landing this chapter's tools onto the empirics of platform pricing and promotion.
:::

::: {.apa-refs}
- Athey, S., Tibshirani, J., & Wager, S. (2019). Generalized random forests. *The Annals of Statistics, 47*(2), 1148-1178. https://doi.org/10.1214/18-AOS1709
- Athey, S., & Wager, S. (2021). Policy learning with observational data. *Econometrica, 89*(1), 133-161. https://doi.org/10.3982/ECTA15732
- Bang, H., & Robins, J. M. (2005). Doubly robust estimation in missing data and causal inference models. *Biometrics, 61*(4), 962-973. https://doi.org/10.1111/j.1541-0420.2005.00377.x
- Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal, 21*(1), C1-C68. https://doi.org/10.1111/ectj.12097
- Chernozhukov, V., Cinelli, C., Newey, W. K., Sharma, A., & Syrgkanis, V. (2026). Long story short: Omitted variable bias in causal machine learning. *The Review of Economics and Statistics*, 1-45. https://doi.org/10.1162/rest.a.1705
- Chernozhukov, V., Demirer, M., Duflo, E., & Fernández-Val, I. (2025). Fisher-Schultz lecture: Generic machine learning inference on heterogeneous treatment effects in randomized experiments, with an application to immunization in India. *Econometrica, 93*(4), 1121-1164. https://doi.org/10.3982/ECTA19303
- Hitsch, G. J., Misra, S., & Zhang, W. W. (2024). Heterogeneous treatment effects and optimal targeting policy evaluation. *Quantitative Marketing and Economics, 22*(2), 115-168. https://doi.org/10.1007/s11129-023-09278-5
- Kennedy, E. H. (2023). Towards optimal doubly robust estimation of heterogeneous causal effects. *Electronic Journal of Statistics, 17*(2), 3008-3049. https://doi.org/10.1214/23-EJS2157
- Nie, X., & Wager, S. (2021). Quasi-oracle estimation of heterogeneous treatment effects. *Biometrika, 108*(2), 299-319. https://doi.org/10.1093/biomet/asaa076
- Robins, J. M., & Rotnitzky, A. (1995). Semiparametric efficiency in multivariate regression models with missing data. *Journal of the American Statistical Association, 90*(429), 122-129. https://doi.org/10.1080/01621459.1995.10476494
- Robinson, P. M. (1988). Root-N-consistent semiparametric regression. *Econometrica, 56*(4), 931-954. https://doi.org/10.2307/1912705
- Wager, S. (2024). *Causal inference: A statistical learning approach* [Unpublished manuscript]. Stanford University.
- Wager, S., & Athey, S. (2018). Estimation and inference of heterogeneous treatment effects using random forests. *Journal of the American Statistical Association, 113*(523), 1228-1242. https://doi.org/10.1080/01621459.2017.1319839
:::
