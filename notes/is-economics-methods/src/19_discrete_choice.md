---
title: "Discrete Choice: From RUM to Mixed Logit"
subtitle: "Random Utility, IIA, and the Road to Random Coefficients"
seriesline: "Foundations of Information Systems Economics · Chapter 19"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 19 · Discrete Choice: From RUM to Mixed Logit"
---

## Introduction

If a food delivery platform raises the delivery fee on restaurant A by one dollar, where do the lost orders go? Perhaps they move to restaurant B down the same street with a similar menu, perhaps they switch to a different cuisine, perhaps the customer simply orders nothing at all. If the researcher only reports "how much did orders fall," the platform still does not know who A's true competitors are. And if the model treats every alternative as equally close, it may even deliver a very precise answer that happens to be very wrong. The difficulty in discrete choice research has never been merely to explain "which one the user picked," but to recover, from choice after choice, the substitution relationships among the alternatives.

Random utility maximization offers a clean language for this. Every option carries a part of its utility that the researcher can see and a part that the researcher cannot, and the consumer picks the one with the highest total utility. Assign a probability distribution to the unobservable part, and individual optimization turns into an estimable choice probability. This step looks like nothing more than adding a random term, but it actually completes the most crucial seam between consumer theory and data. The problems follow immediately: the level and scale of utility need to be normalized, price is often correlated with unobserved quality, and the distributional assumption on the error term directly dictates who substitutes for whom.

This chapter starts from the most computationally convenient model, the logit. Its closed form is a pleasure to work with, but the IIA restriction may force a similar restaurant and a completely different outside option to absorb the lost demand in the same proportion. Nested logit uses within-group correlation to admit that "close neighbors are more alike," and mixed logit lets consumers have different tastes. Price endogeneity is handled separately by instrumental variables or a control function. The Selva delivery-platform case that runs through the whole chapter is deliberately built so the truth is known, so that each model must not only fit the data but also account for exactly how far off it got the elasticities and diversion ratios.

What we really need to deliver, then, is not just a likelihood function that converges, but a substitution logic we can articulate: which normalizations the parameters depend on, whether IIA is credible here, and toward whom a price change pushes consumers. The differences among logit, nested logit, and mixed logit only reveal their economic content once you look at elasticities and diversion ratios. The next chapter will move the lens from individual choice to market shares: once "who bought whom" disappears, the Berry inversion takes over the job that the individual likelihood used to do.

## 1 A number that does not add up

::: {.case}
Selva is a fictional food delivery platform. We tell this case with simulated data, and the advantage is that the data generating process (DGP) is fully known, so every number an estimator hands back can be reconciled against the truth, a teaching condition that real data can never provide. The setup is as follows. The platform operates in 300 city-week markets, and in each market there are four restaurant brands to choose from, labeled A, B, C, D, plus an outside option of "not ordering." The four brands belong to two cuisine nests: A and B are Western comfort food (comfort), and C and D are Asian cuisine (asian). The data are at the order level: 1200 users, each placing 8 orders within the window, for 9600 choices in total, each time picking exactly one of the five options available in that market at that moment. The first thing the platform wants to know is basic: how sensitive are users to price, and if a brand raises its price by one dollar, how many orders does it lose, and to whom?
:::

Start with the most direct approach. Treat each choice as the dependent variable and run a standard conditional logit, putting brand fixed effects, star ratings, and price into utility, estimate the price coefficient by maximum likelihood, and then convert it into a demand elasticity. Brand A's own-price elasticity comes out to $-1.78$: a 1% price increase drops A's order volume by 1.78%. Taken in isolation, this looks like a conclusion you could use directly. The elasticity is less than 2 in absolute value, demand looks not especially sensitive, and written into a brief for the operations team the conclusion would read "users are not that price-conscious, A still has room to raise prices."

The problem is that we know the truth. These data were generated from a world with a distinctly higher price elasticity, and A's true own-price elasticity is $-2.78$. The naive logit understates the price sensitivity of demand by well over a third, and a recommendation to raise prices based on it could point in entirely the wrong direction. The regression itself is not miswritten, maximum likelihood converges cleanly, and the standard errors are small. So where exactly does the $-1.78$ go wrong?

The answer comes in two layers, and this chapter peels them off one at a time. The first layer is endogeneity: expensive brands are often also better along dimensions the researcher cannot see (faster preparation, more consistent flavor, a more prominent recommendation slot), and this unseen quality both raises price and attracts orders on its own, so the price coefficient is diluted by the illusion that "good things are naturally expensive and still sell well," and the elasticity is pulled toward zero. Section 3 will spell out where this bias comes from, and Section 4 uses a control function to push the price coefficient back from $0.141$ to $0.245$, with the elasticity moving accordingly from $-1.78$ to around $-3.00$. The second layer is more subtle, and is more the theme of this chapter: even with the price coefficient fixed, the logit still wrongly assumes that after A raises its price the lost orders are distributed in proportion to the other options' shares, whereas the truth is that most of them pour toward A's same-nest Western sibling B. Fixing this layer takes nested logit and mixed logit. By the time Section 6 has run all four estimators on Selva, A's own-price elasticity lands at $-2.75$, and the within-nest substitution is recovered as well.

## 2 The economic model and the estimand

Before touching estimation, pin down two things: what target quantity we want to estimate, and in what model language we write it down. The target quantity in discrete choice is very different from that in linear regression, and it is worth stating clearly first.

Demand estimation has two core deliverables, one being the own-price elasticity, and the other the substitution pattern, the latter often characterized by the diversion ratio. The own-price elasticity says how much a brand's own sales fall when it raises its price; the diversion ratio says where that lost fraction of orders goes, with $D_{A\to B}$ denoting the share of A's lost orders that B picks up. Together, these two answer the questions that operations and regulators actually care about: how much pricing power a brand has, whether a merger of two brands would markedly raise prices because they are close neighbors, and where delisting an option would push demand. Section 6 will show that both the own-price elasticity and the diversion ratio that the naive logit hands back are wrong, and that the direction of the error is set by the model specification, not by sampling noise. The estimand of this chapter is therefore not any single coefficient, but this whole set of elasticities and diversion ratios, all of which are functions of the utility parameters: estimate the parameters and substitute them into the closed-form formulas.

To write these quantities down we need the notation of random utility. Index the decision maker by $i$, let the choice set for one decision be $\{0,1,\dots,J\}$, where $0$ is the outside option (not ordering) and $1$ through $J$ are the $J$ inside products. Decision maker $i$ gets utility from option $j$ equal to

$$U_{ij} = V_{ij} + \varepsilon_{ij}$$

$V_{ij}$ is the observable part (also called the representative utility), a function of product characteristics and parameters; $\varepsilon_{ij}$ is the part the researcher cannot observe, which is given a probability distribution when it enters the model. The decision maker picks the option with the highest utility, so the probability of choosing $j$ is

$$s_{ij} = \Pr\big(U_{ij} > U_{ik}\ \ \forall k \neq j\big) = \Pr\big(\varepsilon_{ik} - \varepsilon_{ij} < V_{ij} - V_{ik}\ \ \forall k \neq j\big)$$

This expression stitches economics (pick the highest utility) to statistics (a probability over $\varepsilon$). The choice probability is a $J$-dimensional integral over the joint distribution of $\varepsilon$, and exactly what it looks like depends entirely on what distribution you assign to $\varepsilon$. Set it to multivariate normal and you get multinomial probit; set it to independent and identically distributed Type-I extreme value and you get the star of this chapter, the logit, whose integral has a closed form. Choosing a continuous, full-support distribution has two further advantages: the choice probability is smooth and differentiable in $V$ everywhere, and it always lands in the open interval $(0,1)$, so any observed choice pattern can be rationalized by the model.

The representative utility is usually set to be linear. In the Selva setup,

$$V_{ij} = \delta_j + \beta_r\, \mathrm{rating}_{j} - \alpha\, \mathrm{price}_{j} + \xi_{j}$$

$\delta_j$ is the brand fixed effect (the outside option is normalized to $\delta_0 = 0$), $\mathrm{rating}_j$ is the star rating, $\mathrm{price}_j$ is the price, $\alpha > 0$ is the disutility of price, and $\xi_j$ is a product quality that the researcher cannot observe but that both the decision maker and the platform know full well. $\xi_j$ is the source of the endogeneity problem that runs through the whole chapter: it both enters utility and is correlated with price, and Section 3 treats it in detail. For now just remember: that $\xi_j$ in $V_{ij}$ is precisely what the naive logit omits, and hence what biases its price coefficient.

::: {.intuition}
Think of $\varepsilon_{ij}$ as the decision maker's mood of the moment and every chance factor the researcher cannot measure: today they especially want something spicy, last time one place let them down, a friend just recommended some brand. The researcher sees none of this and can only treat it as a random variable. The cleverness of RUM is that it does not pretend to observe these moods; it admits it cannot, packs them into a distribution, and then asks "under this distribution, what is the probability of choosing each option." A deterministic model of individual optimization thereby becomes, in the researcher's eyes, a probability model, and probability is exactly the share we can count in the data. That is how theory and data get connected.
:::

To summarize: the estimand in discrete choice is the whole substitution structure of own-price elasticity and diversion ratio, the language that carries it is random utility, the presence of the outside option lets the "order or not" margin into the model too, and that unobservable $\xi_j$ in utility is the seed of endogeneity. This section has not yet touched what distribution to give $\varepsilon$, which is the first step in the next section's building up of the logit.

## 3 Identification

Identification comes logically before estimation. This section answers: given the RUM framework, which utility parameters and which substitution patterns can the choice shares in the data actually pin down, and which things are in principle not identified and must be supplied by normalization conventions or extra assumptions. This is the most finely analyzed section of the chapter, and as a depth benchmark it is split into numbered subsections, each of which gives intuition before formalization.

### 3.1 From random utility to logit: how the Gumbel distribution gives a closed form

Start with the most elementary question: if the unobservable utility $\varepsilon$ follows some distribution, can the probability of choosing an option be written as an expression we can actually compute. In general this is a high-dimensional integral with no closed form, and can only be simulated numerically. The whole magic of the logit lies in picking the right distribution, one that collapses this integral into a single line of algebra.

That distribution is the Type-I extreme value, also called the Gumbel distribution. Its cumulative distribution function and density are

$$F(\varepsilon) = \exp\big(-e^{-\varepsilon}\big), \qquad f(\varepsilon) = e^{-\varepsilon}\exp\big(-e^{-\varepsilon}\big)$$

The standard Gumbel has mean equal to the Euler-Mascheroni constant $\gamma \approx 0.577$ and variance $\pi^2/6$. It looks like a right-skewed bell, and the reason it is called an extreme value distribution is that it is the limiting distribution of the maximum of a large class of distributions, which makes it natural for characterizing "the most prominent among several unseen factors."

Now do that collapsing integral. Assume $\varepsilon_{i0},\dots,\varepsilon_{iJ}$ are independent and identically distributed standard Gumbel. Given $\varepsilon_{ij}$, choosing $j$ requires $\varepsilon_{ik} < \varepsilon_{ij} + V_{ij} - V_{ik}$ for all $k \neq j$, and by independence the conditional probability is the product of the individual CDFs,

$$\Pr(j \mid \varepsilon_{ij}) = \prod_{k \neq j} \exp\!\Big(-e^{-(\varepsilon_{ij} + V_{ij} - V_{ik})}\Big)$$

Then integrate over $\varepsilon_{ij}$ against its density. Folding the $k = j$ term into the product too (it equals $\exp(-e^{-\varepsilon_{ij}})$, which just completes a piece of the density), rearrange into

$$s_{ij} = \int_{-\infty}^{\infty} e^{-\varepsilon_{ij}} \exp\!\Big(-e^{-\varepsilon_{ij}} \underbrace{\textstyle\sum_{k} e^{-(V_{ij}-V_{ik})}}_{\equiv\, A}\Big)\, d\varepsilon_{ij}$$

Make the substitution $t = e^{-\varepsilon_{ij}}$, so $dt = -e^{-\varepsilon_{ij}}d\varepsilon_{ij}$, and the integral becomes $\int_0^\infty e^{-At}\,dt = 1/A$. And since $A = \sum_k e^{-(V_{ij}-V_{ik})} = e^{-V_{ij}}\sum_k e^{V_{ik}}$, substituting back gives

$$s_{ij} = \frac{e^{V_{ij}}}{\sum_{k} e^{V_{ik}}}$$

This is the logit choice probability, in multinomial logistic form. With the outside-option normalization $V_{i0} = 0$, that term in the denominator is just $e^0 = 1$, so

$$s_{ij} = \frac{e^{V_{ij}}}{1 + \sum_{k=1}^{J} e^{V_{ik}}}$$

This one line makes an important point clear: the probability of choosing an option depends only on the differences in representative utility across options, and varies with utility in a smooth, monotone way. The Gumbel assumption buys exactly this closed form, and the price it charges is the IIA that the next several subsections have to pay back.

::: {.intuition}
The quantity $\log\big(\sum_k e^{V_{ik}}\big)$ is called the inclusive value or log-sum, and it is in fact the expected maximum utility $\mathbb{E}[\max_k U_{ik}]$ available to the decision maker (up to a constant $\gamma$). It is a smooth approximation to the max: $\max(10, 20) = 20$, while $\log(e^{10}+e^{20}) \approx 20.00005$, almost identical, but differentiable everywhere. When we get to nested logit, this log-sum will reappear as "the attractiveness of a nest"; when we discuss welfare, it is the handle on the change in consumer surplus. Remember that in the logit world, "the total value of a set of options" is the log-sum-exp of their utilities.
:::

### 3.2 Only differences and scale are identified: two normalizations

RUM has two natural indeterminacies that cannot be pinned down without conventions, and they must be stated up front or the coefficients later will be uninterpretable.

The first is location (level). Choice depends only on the relative heights of utilities, so adding the same constant $a$ to the utilities of all options leaves the ranking untouched and the choice probabilities unchanged. Algebraically, replace each $V_{ij}$ with $V_{ij} + a$, and numerator and denominator both multiply by $e^a$ and cancel. The consequence is that the absolute level of utility is not identified, only differences are. The convention in practice is to set the outside option's utility to zero and read every other option's $\delta_j$ relative to it. By the same logic, anything that varies with the decision maker but not the option (for instance income $Y_i$ entering utility on its own) appears in the same amount across all options and cancels in the difference, and is uniformly not identified, unless it interacts with an option characteristic (for instance $\mathrm{price}_j / Y_i$), in which case the interaction does bring cross-option variation.

The second is scale. Multiply all utilities by the same positive number $\lambda$, and the ranking still does not change, but $\mathrm{Var}(\lambda\varepsilon) = \lambda^2\mathrm{Var}(\varepsilon)$ does. This says the scale of utility and the variance of $\varepsilon$ are one and the same thing, and pinning down one pins down the other. The logit convention fixes the variance of $\varepsilon$ at $\pi^2/6$ (the variance of the standard Gumbel). This convention has an easily overlooked corollary: every coefficient is really "the true coefficient divided by the standard deviation of $\varepsilon$." The larger the fluctuation in unseen utility (the stronger the noise), the more the estimated coefficient is compressed toward zero. So the logit coefficients from two data sets cannot be compared in magnitude directly, since their scales may differ; what can be compared is the ratio of coefficients, for instance $\beta_r / \alpha$, which is the marginal rate of substitution between rating and price, with the scale cancelling in numerator and denominator.

::: {.warning}
Interpreting the coefficients of a logit or any discrete choice model as "the absolute marginal utility" is wrong; they only have relative meaning. Comparing coefficient magnitudes across models or across samples, when the two error scales differ, may produce differences that are purely differences in scale. To compare, compare identified ratios (such as willingness-to-pay $= \beta_x/\alpha$) or elasticities, which are immune to the scale normalization.
:::

These two normalizations are not technical fastidiousness but one half of the price-elasticity puzzle from Section 1. The willingness-to-pay for a rating is $\beta_r / \alpha$, a dollar-denominated ratio; as long as the price coefficient $\alpha$ is biased by endogeneity, this ratio is wrong along with it, and the scale normalization by itself cannot save it, which takes the control function of Section 4.

### 3.3 IIA: the logit's strongest restriction

The logit's closed form is beautiful, but it quietly takes on a very strong behavioral restriction, and understanding it is the key to understanding the whole back half of the chapter. Take any two options $j$ and $k$; the ratio of their choice probabilities is

$$\frac{s_{ij}}{s_{ik}} = \frac{e^{V_{ij}}}{e^{V_{ik}}} = e^{V_{ij} - V_{ik}}$$

This ratio depends only on the utilities of $j$ and $k$ themselves, and has nothing whatever to do with what some third option $\ell$ is, whether it exists, or how good it is. This is the independence of irrelevant alternatives, IIA. It sounds like a harmless technicality, but it is actually a strong behavioral assertion.

The classic counterexample is red-bus/blue-bus. Imagine a commute with only two options, driving and taking a blue bus, each with half the share, so the ratio is 1. Now a red bus is introduced, identical to the blue bus except for color. IIA forces the ratio of "car to either bus" to remain 1, so the model predicts the three options split into thirds each, with driving's share dropping from 1/2 to 1/3. But common sense tells us that the red and blue buses are nearly the same thing to a commuter, that the red bus should mostly steal riders from the blue bus, that driving's share should stay near 1/2, and that the two buses should each take about 1/4. The logit goes wrong precisely by treating the red and blue buses as new options that are independent of each other just as "driving" is, without recognizing that the red and blue buses are close neighbors. The root cause is that the logit assumes the $\varepsilon$ are mutually independent, and does not allow "my preference for the blue bus and my preference for the red bus are highly positively correlated."

Translate IIA into elasticities and the problem is clearer. The logit's own-price and cross-price elasticities have closed forms:

$$\eta_{jj} = \frac{\partial s_{ij}}{\partial p_j}\frac{p_j}{s_{ij}} = -\alpha\, p_j\,(1 - s_{ij}), \qquad \eta_{jk} = \frac{\partial s_{ij}}{\partial p_k}\frac{p_k}{s_{ij}} = \alpha\, p_k\, s_{ik}$$

Fix your gaze on the cross elasticity $\eta_{jk}$: it depends only on the price and share of the option $k$ that raised its price, and has nothing at all to do with which option $j$ is being affected. In other words, after $k$ raises its price, the shares of all other products are pushed up by the same proportion, regardless of how similar they are to $k$. This is called proportional substitution. Carried down to the diversion ratio,

$$D_{j \to k} = \frac{s_{ik}}{1 - s_{ij}}$$

$j$'s lost demand is distributed entirely in proportion to the shares of the remaining options, and there is simply no place in the model for "who is closer to whom." This is exactly the deeper error from Section 1: when Selva's A raises its price, the truth is that orders flow mainly to B, which is also Western, but the logit merely apportions them by share across B, C, D, and the outside option, siphoning away the part that should have gone to B. The numbers in Section 6 will nail this down.

IIA is not all bad. It has one historically useful corollary: since the ratio of any two options is independent of a third, the $\beta$ estimated from any subset of options is consistent, and this was once used to subsample when the options were very numerous, to save computation. But for modern demand estimation, IIA is a restriction that must be relaxed, not a convenience to be exploited, because the substitution pattern is exactly the product we need to deliver.

### 3.4 Price endogeneity: another threat to identification

IIA is an ailment of the substitution structure, while price endogeneity is an ailment of the level parameters; the two are independent and must be treated separately. Return to that $\xi_j$ in utility. It is a product quality the researcher cannot observe but the platform and users know clearly: more reliable preparation, better packaging, a more prominent recommendation slot. The platform sees $\xi_j$ when it prices, good things are naturally priced higher, so $\mathrm{Cov}(\mathrm{price}_j, \xi_j) > 0$. The naive logit stuffs $\xi_j$ into the error, and since price is positively correlated with this error, the standard endogeneity bias arrives.

The direction of the bias can be worked out. High-$\xi$ products are both expensive and popular, and if you do not control for $\xi$, the model sees that "the expensive ones actually sell well" and pushes the disutility of price toward zero, so $\alpha$ is underestimated and demand looks less sensitive than it truly is. This is exactly where the $-1.78$ of Section 1 comes from: the true $\alpha = 0.30$, but the naive estimate is only $0.141$, less than half. Note that this bias is unrelated to sample size, and no amount of extra data will remove it; it is a problem at the level of identification. Note too that it is present just as much in individual choice data, and is not caused by "data that are not fine enough"; even if we observed every choice of every person, as long as price is correlated with $\xi$, $\alpha$ is biased.

Curing it requires an instrument: a variable that can move price but does not enter utility directly. In Selva this variable is each market's local input cost (courier hourly wages, food cost, and the like), which pushes up price but does not affect the user's preference for a given order, satisfying the exclusion restriction. With an instrument in hand, there are two roads to setting $\alpha$ right. One is the control function, which is especially handy on individual data and is detailed in Section 4; the other is to invert shares into a product-level $\delta_j$ and then do IV-GMM, which is the main battleground of market-level data, taken up later in this series when we discuss market-level demand estimation. Here we just set up the identification logic: the root of price endogeneity is that $\xi$ enters both utility and pricing, and the cure is to find exogenous variation that enters pricing only, not utility.

### 3.5 Relaxing IIA: what nested logit and random coefficients change

The root of IIA is that the $\varepsilon$ are independent, so relaxing it means allowing the $\varepsilon$ to be correlated, or equivalently, letting each person have their own tastes. What each of these two roads identifies is worth spelling out first.

Nested logit takes the road of "building a structure for the correlation." Sort the options into several nests, which in Selva are the Western nest $\{A, B\}$, the Asian nest $\{C, D\}$, and the outside option as a nest of its own. The unobservable utilities of same-nest products are allowed to be positively correlated, while across nests they remain independent. This is achieved through a generalized extreme value (GEV) joint distribution, and the choice probability it gives is

$$s_{ij} = \frac{e^{V_{ij}/\lambda_g}\Big(\sum_{k \in g} e^{V_{ik}/\lambda_g}\Big)^{\lambda_g - 1}}{\sum_{h}\Big(\sum_{k \in h} e^{V_{ik}/\lambda_h}\Big)^{\lambda_h}}$$

![The nest structure of the nested logit: one order is first allocated among the outside option and the two cuisine nests, and then among the products within a nest; the unobservable utilities of same-nest products are allowed to be positively correlated. The smaller the nest parameter $\lambda$, the stronger the within-nest correlation, and the value marked in the figure is the $\lambda = 0.875$ estimated on Selva in Section 4.3.](assets/fig/fig_19_nest.svg)

Here $g$ is the nest that $j$ belongs to (this chapter uses $g$ and $h$ for nests, unrelated to the book-wide convention of writing cohort as $g$), and $\lambda_g \in (0, 1]$ is the nest parameter. $\lambda_g = 1$ falls back to the ordinary logit (no extra within-nest correlation); $\lambda_g \to 0$ means same-nest products are nearly perfect substitutes (correlation tending to 1). IIA still holds within a nest but no longer holds across nests, so "A raises its price and flows mainly to B" can now be accommodated. What identifies $\lambda_g$? Variation of the kind "how changing the characteristic of some product in a nest affects the whole nest's share relative to other nests": how strong within-nest substitution is leaves an imprint in the response of the nest's total share to characteristics. In the log-sum notation of Section 3.1, you can read it as a two-level choice: first pick a nest by the nests' inclusive values, then pick a product within the nest by a logit, $s_{ij} = s_{ij \mid g}\cdot s_{ig}$, though the essence is a correlation structure of $\varepsilon$ rather than a genuine two-step process. The price is that the nests must be specified by the researcher before estimation, and if the nests are drawn wrong, the substitution pattern is wrong too.

Mixed logit (random coefficients logit) takes the road of "letting tastes differ across people," which is more flexible. Suppose each person's coefficients $\beta_i$ follow some distribution $F(\beta \mid \theta)$, for instance the price coefficient $\alpha_i = \bar\alpha + \sigma_\alpha \nu_i$ and the preference for the Western nest $b_i = \sigma_b \eta_i$ (the random coefficient $b_i$ here is unrelated to the nest parameter $\lambda_g$ of the nested logit, and is denoted separately to avoid confusion), where $\nu_i, \eta_i$ are the individual's taste draws. Given $\beta_i$, each person is internally still a logit; integrating the individual logit over the taste distribution gives the aggregate choice probability

$$s_{ij}(\theta) = \int \frac{e^{\,\beta\, x_{ij}}}{\sum_k e^{\,\beta\, x_{ik}}}\, dF(\beta \mid \theta)$$

This integral generally has no closed form and must be computed by simulation. Why does it relax IIA? Because substitution is no longer determined by shares but by "who is like whom in characteristics": two products that are both cheap are favored simultaneously by the same group of price-sensitive people, so when one raises its price, demand flows more toward the other cheap one, rather than being scattered by share. McFadden and Train (2000) proved that as long as the characteristics are chosen richly enough, mixed logit can approximate any random utility model, including probit, making it the most general tool in discrete choice.

The identification of mixed logit has one point that must be spelled out, and it also explains a deliberate design of this chapter's Selva data. The random coefficient distribution does not in principle require a panel: if the choice sets and product characteristics facing different decision makers have sufficiently rich variation satisfying the exclusion restriction, a cross section of single choices can identify the taste distribution too. The extra value of a panel is that the correlation of the same person's choices across different situations directly reveals persistent tastes, and helps separate this heterogeneity from the choice-by-choice independent utility shocks. Letting each Selva user place 8 orders uses precisely the internal consistency of these 8 choices (always picking the cheap one, always picking Western) to strengthen the identification of $\sigma_\alpha$ and $\sigma_b$. Under this chapter's particular DGP and limited cross-sectional variation, cutting each person's choices to 1 would make the variance-term estimates collapse toward zero; that is a symptom of insufficient information in this example, and cannot be extrapolated into a general claim that "a pure cross section cannot identify mixed logit."

The key points of this section can be summarized as follows: the logit's choice probability identifies only differences and ratios of utility, with the absolute level and scale supplied by normalizations; IIA locks substitution into proportional allocation and is a strong restriction that must be relaxed; price endogeneity is another threat, independent of IIA, arising because $\xi$ enters both utility and pricing, and is cured with an instrument; relaxing IIA has two roads, nested logit (build a structure for the correlation, identify the nest parameter) and mixed logit (let tastes vary across people, identify their distribution from rich choice-set variation, with a panel able to strengthen identification markedly). None of this has yet touched "how to estimate," which is the business of the next section.

## 4 Estimation: from plain logit to mixed logit

This section advances along the logic of "where the previous method fails gives rise to the next." Plain logit hands off two ailments, endogeneity bias and IIA; the control function cures endogeneity, and nested logit and mixed logit cure IIA.

### 4.1 Maximum likelihood for plain logit

With the choice probability in closed form, estimating plain logit is standard maximum likelihood. Having observed $N$ choices, with $d_{ij} = 1$ indicating that choice $i$ picked $j$, the log-likelihood is

$$\ell(\theta) = \sum_{i} \sum_{j} d_{ij}\, \log s_{ij}(\theta)$$

The logit log-likelihood is globally concave (log-sum-exp is a convex function), so there is no need to worry about separated local maxima. A finite optimum is unique only when the design matrix has full rank and the sample has no complete or quasi-complete separation; if separation occurs, the supremum of the likelihood is approached as the coefficients diverge, and the software reporting "converged" does not mean the MLE exists. Its score is also well structured,

$$\frac{\partial \ell}{\partial \beta} = \sum_i \sum_j \big(d_{ij} - s_{ij}\big)\, x_{ij}$$

The first-order condition is that "the difference between observed choice and predicted probability is orthogonal to the characteristics," which is the intuitive source of the moment condition. Running this plain logit on Selva, the price coefficient is $0.141$ (SE $0.007$), the rating coefficient $0.554$ (SE $0.089$), converting into a willingness-to-pay for a rating of $0.554 / 0.141 \approx \$3.93$ per star, and an own-price elasticity of $-1.78$. These are the set of numbers that did not add up in Section 1, ailing in two places: the price coefficient is biased by $\xi$, and substitution is locked by IIA. We treat them in turn below.

### 4.2 The control function cures endogeneity

The idea of the control function is plain: the root of endogeneity is the component of $\xi_j$ that is correlated with price, and if we can estimate this component, put it into the regression, and control it out, then the remaining price variation is clean. Two steps. First, regress price on the instrument and exogenous characteristics,

$$\mathrm{price}_{j} = f(\mathrm{rating}_j,\ w_j,\ \text{brand}) + e_{j}$$

$w_j$ is the local cost, the instrument, and take the residual $\hat e_j$. Second, put this residual as an extra control variable (that is, the control function) into utility,

$$V_{ij} = \delta_j + \beta_r\, \mathrm{rating}_j - \alpha\, \mathrm{price}_j + \rho\, \hat e_j$$

and run the logit again. The reasoning is that the component of $\xi_j$ correlated with price is absorbed by $\hat e_j$, so the price coefficient $\alpha$ is now identified on the "purified" variation and returns to consistency. In Selva's first step the coefficient on the instrument $w$ is $1.605$, with a partial-F as high as $2351$, a very strong instrument, so there is no need to worry about weak instruments. In the second step the price coefficient jumps from $0.141$ to $0.245$ (SE $0.008$), and the control function coefficient $\rho = 0.299$ is significantly positive, exactly as expected under "higher-quality products are priced higher." The own-price elasticity moves accordingly from $-1.78$ to $-3.00$. This layer of endogeneity is now cured.

::: {.warning}
The control function buys consistency but also comes with two disciplines. One is that it requires the first-step price equation to be essentially correctly specified, which is a stronger assumption than pure IV, because IV only needs the instrument to be valid, whereas the control function further depends on the specification of the relationship between $\xi$ and $e$. The other is that the second-step standard error must account for the first step being estimated, since $\hat e_j$ is not the true residual, and using the ordinary second-step SE directly will understate it. The reliable practice is to bootstrap the two steps together (resampling at the market level), or to use the Murphy-Topel two-step correction. This series' module on inference gives the variance formula for this kind of two-step estimator.
:::

The control function fixes the price coefficient but does not touch IIA: this logit still assumes proportional substitution. Section 6 will show that the diversion pattern of the control-function logit is almost identical to that of the naive version, with the diversion of A toward B still only around $0.25$, whereas the truth is $0.34$. To fix substitution, you have to change the model.

### 4.3 Nested logit cures within-nest substitution

Nested logit sorts A, B into the Western nest and C, D into the Asian nest, allowing same-nest products to substitute more strongly. Estimation is still maximum likelihood, with one extra parameter, the nest parameter $\lambda$ (here a common $\lambda$ is set for the two inside nests, and the outside option's nest is fixed at $\lambda_0 = 1$). The likelihood uses the GEV choice probability of Section 3.5, and to handle price endogeneity along the way, the utility continues to carry the control function residual $\hat e_j$. One point to watch: the nested logit likelihood is no longer globally concave and may have multiple local solutions, so in practice one must try several starting values and constrain $\lambda$ within $(0, 1]$ (this chapter uses the reparameterization $\lambda = \mathrm{logistic}(\theta)$ to guarantee this).

On Selva, $\hat\lambda = 0.875$ (SE $0.063$), significantly less than 1, indicating there is indeed extra within-nest correlation. The nest parameter can be translated into the correlation of same-nest unobservable utilities, about $1 - \lambda^2 = 0.235$. The price coefficient is $0.227$, the rating coefficient $0.644$. The most important change is in the substitution pattern: the diversion of A toward B upon a price rise climbs from the logit's $0.25$ to $0.32$, moving clearly toward the truth of $0.34$, because the model now recognizes that B is A's same-nest neighbor. The price of nested logit is that the nest structure is rigidly specified by the researcher, and drawing the nests wrong forfeits everything; moreover it assumes price sensitivity is the same for everyone, ignoring taste heterogeneity.

### 4.4 Mixed logit and simulated maximum likelihood

Mixed logit does not pre-draw nests, but lets each person have their own price sensitivity $\alpha_i = \bar\alpha + \sigma_\alpha\nu_i$ and Western preference $b_i = \sigma_b\eta_i$, and the substitution pattern emerges naturally from the data. Its choice probability is an integral with no closed form, and must be estimated by simulated maximum likelihood (SML): for each person draw $R$ sets of tastes $(\nu_i^{(r)}, \eta_i^{(r)})$, use them to compute $R$ individual logit probabilities, and average them as the simulated value of that person's choice probability.

The key in this example is that the likelihood preserves the panel structure. The same person $i$'s 8 choices share the same set of persistent tastes $\beta_i$, so one must first, given $\beta_i^{(r)}$, take the product of the probabilities of this person's 8 choices, and then average over $r$:

$$\hat P_i(\theta) = \frac{1}{R}\sum_{r=1}^{R} \prod_{t=1}^{T_i} \frac{e^{\,\beta_i^{(r)} x_{i c_{it} t}}}{\sum_k e^{\,\beta_i^{(r)} x_{i k t}}}, \qquad \ell(\theta) = \sum_i \log \hat P_i(\theta)$$

$c_{it}$ is the option chosen by person $i$ on the $t$-th occasion. The internal consistency of the 8 choices (whether this person always picks the cheap one, always picks Western) is the source of identification for $\sigma_\alpha, \sigma_b$. Price endogeneity is still handled with the control function residual.

Running $R = 100$ sets of simulation draws on Selva, mixed logit recovers the taste heterogeneity: the mean price coefficient $\bar\alpha = 0.293$ (SE $0.010$), the standard deviation $\sigma_\alpha = 0.116$ (SE $0.004$), and the Western-preference standard deviation $\sigma_b = 0.807$ (SE $0.043$). Compared with the true values $\bar\alpha = 0.30$, $\sigma_\alpha = 0.12$, $\sigma_b = 0.90$, nearly all on target. The substitution pattern is also the closest to the truth: the diversion of A toward B upon a price rise is $0.328$, and toward the outside option $0.250$, both ends closest to the true values $0.342$ and $0.238$; A's own-price elasticity is $-2.75$, almost exactly the true $-2.78$.

::: {.warning}
SML has a property one must know: when using the usual biased simulated log likelihood, a fixed number of simulation draws $R$ leaves a simulation bias that does not vanish with sample size, so it is generally inconsistent, and especially severe when the choice probability is small. In practice, let $R$ grow together with the sample size, or switch to the method of simulated moments (the simulation error in the moment conditions enters linearly and not through the $\log$, so the bias is smaller). In addition the mixed logit likelihood is non-convex, so one should supply an analytic gradient, use multiple starting values, and when necessary use quasi-Monte Carlo (such as Halton sequences) to reduce the simulation variance. This chapter's $R = 100$ already recovers the parameters stably in this simulated sample; if each person's choices are cut to 1, the variance terms in this example collapse toward zero, which shows the current cross-sectional variation is insufficient and cannot be used to claim that all cross-sectional mixed logits are unidentified.
:::

### 4.5 Elasticities and diversion: what demand estimation truly needs to deliver

Estimating the parameters is only half the product; what needs to be delivered are the elasticities and diversion ratios. All three models have computable forms for their elasticities. Logit and nested logit have closed-form elasticities (the former in Section 3.3, the latter obtained by differentiating the GEV probability, with the within-nest and cross-nest cross elasticities differing). The aggregate elasticity of mixed logit is a weighted average of the individual logit elasticities over the taste distribution,

$$\eta_{jk} = \frac{p_k}{s_j}\int \frac{\partial s_{ij}}{\partial p_k}\, dF_i, \qquad \frac{\partial s_{ij}}{\partial p_k} = \begin{cases} -\alpha_i\, s_{ij}(1 - s_{ij}) & k = j \\ +\alpha_i\, s_{ij} s_{ik} & k \neq j\end{cases}$$

It is precisely because different people are weighted differently (price-sensitive people carry more weight in the substitution among cheap products) that mixed logit breaks free of the shackle of proportional substitution. The diversion ratio is analogous, with $D_{j\to k}$ being a weighted average of the individual diversions $s_{ik}/(1 - s_{ij})$, the weights differing by situation (price rise, delisting, quality change). It is only when these quantities are computed that we get what operations and regulators want to see: who has pricing power, whether a merger would raise prices, and where a delisting pushes demand. The reconciliation table in Section 6 will place the elasticities and diversions of the four models side by side, against the truth.

A word on inference: the standard errors of logit and nested logit come from the Hessian of the maximum likelihood; mixed logit is the same but with attention to simulation error; with a control function, all must fold in the estimation uncertainty of the first step. If the same user's multiple choices are correlated (as a panel naturally is), the standard errors should be clustered at the user level, or they will be understated. The route of this section is now complete: plain logit reveals the ailments, the control function cures endogeneity, nested logit and mixed logit cure IIA, and elasticities and diversion are the final products.

## 5 Anchoring papers

A method only stands up once it lands in real research. Three anchoring papers: one is the methodological origin of random utility, one gives the standard solution for handling endogeneity on individual data, and one is a representative application of discrete choice in the context of information goods. Each is laid out along five elements, paper, method, data, results, limitations, with attention to how the assumptions are defended.

### 5.1 McFadden (1974)

::: {.case}
Paper and place in methodological history: "Conditional Logit Analysis of Qualitative Choice Behavior," in the volume Frontiers in Econometrics edited by Zarembka. This paper connected random utility maximization from the choice theory of psychology to econometrics, gave the estimation framework of conditional logit, and is the origin of the whole of discrete choice demand analysis, for which McFadden shared the 2000 Nobel Prize in Economics with Heckman.

Method: starting from RUM, set the unobservable utility to independent and identically distributed Type-I extreme value, derive the closed-form choice probability of conditional logit, and estimate by maximum likelihood. Unlike the earlier multinomial logit, conditional logit lets the choice probability depend on the option's own characteristics (such as travel time and cost), not just on the decision maker's characteristics, and it is this that makes demand analysis possible.

Data: the travel-mode choices of San Francisco Bay Area commuters, used to predict the ridership of the then-unbuilt BART rail transit. The options are driving, bus, rail, and so on, and the characteristics are the time and cost of each mode.

Results: the model's pre-BART prediction of the mode share was 6.3%, and the actual post-construction figure was about 6.2%, almost hitting the truth on the nose; by contrast, the official prediction of the time, based on aggregate methods, was as high as around 15%, nearly double the actual figure. This ability to "use a model to predict an option that does not yet exist" is exactly where the value of RUM over purely descriptive methods lies.

Limitations: conditional logit carries IIA, and the independence assumption of Type-I EV locks substitution into proportional allocation, unable to characterize red-bus/blue-bus close-neighbor substitution. McFadden himself was clear about this, and the later nested logit (GEV) and mixed logit were both created to relax it. The historical significance of this paper lies half in its being a pioneer, and half in its clearly identifying the very restriction that later generations would spend decades loosening.
:::

### 5.2 Petrin and Train (2010)

::: {.case}
Paper: "A Control Function Approach to Endogeneity in Consumer Choice Models," Journal of Marketing Research. It gives a practical and general solution for handling price endogeneity on individual choice data, and is the direct source of this chapter's Section 4.2.

Method: the control function. In the first step, regress the endogenous price on the instrument and exogenous characteristics and take the residual; in the second step, put the residual as an extra variable into utility and estimate the discrete choice model. Compared with BLP-style share inversion plus IV-GMM, the advantage of the control function is that it can be used directly on data that retain the complete individual choices, is simple to implement, and is naturally compatible with mixed logit.

Data: the paper demonstrates with consumer-goods data such as television-service choices, where price is correlated with unobservable product quality, constituting a typical endogeneity.

Results: the price coefficient and substitution elasticities after the control-function correction agree closely with those from BLP-style instrumental methods, showing that on individual data this simpler road can reach the same destination. This contrast of "two roads converging" is the source of its persuasiveness.

Limitations: the control function requires the first-step price equation to be correctly specified, an assumption stronger than pure IV; the standard errors of the two-step estimation need a special correction and cannot use the ordinary second-step SE directly. The authors account for both points, and this practice of putting the method's boundary of applicability out in the open is worth learning.
:::

### 5.3 Gentzkow (2007)

::: {.case}
Paper: "Valuing New Goods in a Model with Complementarity: Online Newspapers," American Economic Review. The question is how readers in the Washington area choose among the print Washington Post, its online edition, and reading neither, a representative work of discrete choice demand in the context of information goods, and it strikes directly at the core of this chapter: how unobservable taste heterogeneity distorts the substitution structure.

Method: an individual-level discrete choice demand model, whose specification allows the print and online editions to be either substitutes or complements (this is what complementarity in the title means, referring to what the model allows rather than what the conclusion asserts), with the key being to incorporate both observable and unobservable reader heterogeneity at once. Identification relies on the covariance structure of individual consumption bundles, not on the shares of single choices.

Data: individual-level reading microdata for Washington-area readers, including reader characteristics and consumption of the two media, print and online.

Results: the most memorable point of this paper is precisely a lesson about heterogeneity. A naive model without heterogeneity makes the online and print editions look like complements, as if launching the online edition would not cannibalize print; but once both observable and unobservable reader heterogeneity are put in, the two are in fact significant substitutes, and that layer of "complementarity" was a false appearance manufactured by unobservable heterogeneity. On the substitution estimate, the online edition pulled down print daily circulation by about 27000 copies, so the cannibalization is real. Ignoring heterogeneity not only underestimates substitution, it gets the sign of substitution backward, which is exactly the through-line of this chapter's mixed logit: flattening out heterogeneity systematically distorts the substitution pattern.

Limitations: identifying heterogeneity places high demands on model specification and data structure, and the individual-level preference correlation must be discernible from the data; the sample is newspaper readers of a specific period and region, and extrapolation to other information goods calls for caution. The value of this paper lies in demonstrating how the discrete choice framework can answer the core IS question of digital products, and why getting the heterogeneity right, and thereby getting the substitution structure right, is the key product.
:::

Putting the three together, the significance of anchoring is clear: McFadden set up RUM and logit and flagged IIA as the restriction to be loosened, Petrin-Train gave the standard solution for curing endogeneity on individual data, and Gentzkow showed how, on the home turf of information goods, getting reader heterogeneity right, and thereby getting even the sign of substitution right, changes the substantive conclusion. The advance of methods has always revolved around the two through-lines of "get heterogeneity right so as to loosen IIA, and get $\xi$ under control."

## 6 A full walkthrough on the Selva data

Now run all the tools of Section 4 through the full Selva panel. The code below uses R 4.5.3, with set.seed(15) fixed for reproducibility, and every number quoted in the text comes from the actual run output of this code.

### 6.1 DGP

The design parameters are as follows: 300 city-week markets, four brands A, B, C, D belonging to the Western nest $\{A, B\}$ and the Asian nest $\{C, D\}$, plus an outside option; 1200 users each placing 8 orders, for 9600 order-level choices in total. The true model is a random coefficients logit, where price sensitivity $\alpha_i$ and Western preference $b_i$ vary across people, and price is positively correlated with unobservable quality $\xi$.

```r
set.seed(15)
delta_true <- c(0.80, 0.50, 0.60, 0.30)   # brand fixed effects, outside = 0
beta_r     <- 0.90                          # rating taste (per star)
alpha_bar  <- 0.30; sigma_a <- 0.12         # price coefficient: mean and heterogeneity
sigma_b    <- 0.90                          # sd of Western-nest taste (within-nest correlation)
sigma_xi   <- 0.45                          # unobservable quality xi

# price inversion: rising with local cost w, quality xi, rating => price positively correlated with xi (endogenous)
grid[, price := price_base[j] + 1.6*w + 2.2*xi + 0.8*(rating-4) + rnorm(.N, 0, 0.6)]
# user-level random coefficients, fixed within their 8 orders
usr[, alpha_i := alpha_bar + sigma_a*rnorm(Nusr)]
usr[, b_i     := sigma_b*rnorm(Nusr)]
```

Simulate each choice from the true model: given the user's random coefficients, compute the utilities of the five options, add a Type-I EV draw, and take the argmax.

```r
choose_one <- function(mk, ai, bi){
  g <- grid[.(mk), on = "m"]
  V <- delta_true[g$j] + beta_r*g$rating - ai*g$price + bi*(g$nestj==1L) + g$xi
  V <- c(0, V)                                 # outside option utility 0
  which.max(V - log(-log(runif(J+1)))) - 1L    # take max after Type-I EV draw
}
```

In the generated data the outside-option share is $0.252$, and the four brands' shares are $0.186$, $0.209$, $0.160$, $0.193$. The correlation of price with $\xi$ is $0.355$, and with rating $0.567$, and price ranges from $6.6$ to $22.4$ dollars. The design intent of this DGP is: the true substitution pattern is determined by within-nest taste correlation, the true price sensitivity is $0.30$, and let us see which estimator can recover them.

### 6.2 Plain logit: reproducing the opening numbers

```r
negll_logit <- function(par){
  V  <- Xmat %*% par                       # brand dummies + rating + price
  Vm <- matrix(V, ncol = J+1L, byrow = TRUE)
  P  <- exp(Vm - apply(Vm,1,max)); P <- P / rowSums(P)
  -sum(log(rowSums(P * chosen_matrix)))
}
fit <- optim(start, negll_logit, method = "BFGS", hessian = TRUE)
```

The price coefficient is $0.141$ (SE $0.007$), the rating coefficient $0.554$ (SE $0.089$). The willingness-to-pay for a rating is $0.554 / 0.141 \approx \$3.93$ per star, and A's own-price elasticity is $-1.78$. These are the set of numbers that did not add up in Section 1, ailing because the price coefficient is biased by $\xi$ and substitution is locked by IIA.

### 6.3 Control function: setting the price coefficient right

```r
fs <- lm(price ~ w + rating + factor(j), data = grid)   # first step
grid[, cf := residuals(fs)]                              # control function residual
# second step: put cf as an extra variable into the logit utility and re-estimate
```

In the first step the coefficient on the instrument $w$ is $1.605$, with a partial-F as high as $2351$, a strong instrument. In the second step the price coefficient jumps from $0.141$ to $0.245$ (SE $0.008$), and the control function coefficient $\rho = 0.299$ is significantly positive, consistent with "higher-quality products are priced higher." The WTP for a rating becomes $\$2.82$, and A's own-price elasticity moves to $-3.00$. Endogeneity is cured, but the next subsection will show the substitution pattern has barely moved.

### 6.4 Nested logit: recovering within-nest substitution

```r
negll_nl <- function(par){
  lam <- plogis(par[8])                    # constrain lambda in (0,1)
  V  <- baseV(par) + par[6]*PRICE          # includes control-function residual
  Vs <- V; Vs[, inside] <- V[, inside]/lam
  eV <- exp(Vs)
  D1 <- rowSums(eV[, nest1]); D2 <- rowSums(eV[, nest2]); Dout <- eV[,1]
  den <- Dout + D1^lam + D2^lam
  # ... within-nest GEV probability, see the formula in Section 3.5
}
```

$\hat\lambda = 0.875$ (SE $0.063$), significantly less than 1, with the correlation of same-nest unobservable utilities about $1 - \lambda^2 = 0.235$. The price coefficient is $0.227$, the rating coefficient $0.644$. The diversion of A toward B upon a price rise climbs from the logit's $0.25$ to $0.32$, moving toward the truth of $0.34$.

### 6.5 Mixed logit: panel SML

```r
negll_mxl <- function(par){
  logLi_r <- matrix(0, Nuser, R)
  for (r in 1:R){                          # R = 100 sets of simulation draws
    alpha_r <- par["alpha_bar"] + par["sigma_a"]*NU[,r]   # user-level price coefficient
    b_r     <- par["sigma_b"]*ETA[,r]
    V <- baseV(par) - alpha_r*PRICE + b_r*NEST1
    P <- exp(V - apply(V,1,max)); P <- P/rowSums(P)
    # sum the log probabilities of a user's 8 orders first (panel structure), then average over r
    logLi_r[,r] <- tapply(log(rowSums(P*CHOSEN)), uid, sum)
  }
  -sum(log(rowMeans(exp(logLi_r))))        # stabilized numerically with log-sum-exp
}
```

With $R = 100$ sets of draws, the mean price coefficient $\bar\alpha = 0.293$ (SE $0.010$), the standard deviation $\sigma_\alpha = 0.116$ (SE $0.004$), and the Western-preference standard deviation $\sigma_b = 0.807$ (SE $0.043$), compared with the true values $0.30$, $0.12$, $0.90$, nearly all on target. A's own-price elasticity is $-2.75$, almost exactly the true $-2.78$.

![Where demand flows, by diversion ratio, after brand A raises its price by one dollar. Plain logit (red) apportions the lost orders roughly by share across the options, sending only about 0.25 to the same-nest sibling B and too much to the outside option; nested (blue) and mixed (olive) recover the high diversion to B, close to the truth (navy).](assets/fig/fig_19_substitution.svg)

### 6.6 The full reconciliation

Place the key products of the four estimators side by side with the truth. Price sensitivity and WTP:

| Estimator | Price coefficient $\alpha$ | Rating WTP (per star) | A own-price elasticity |
|---|---|---|---|
| Truth (DGP) | 0.300 | \$3.00 | -2.78 |
| plain logit (naive) | 0.141 | \$3.93 | -1.78 |
| logit + control function | 0.245 | \$2.82 | -3.00 |
| nested logit | 0.227 | \$2.84 | -3.04 |
| mixed logit (panel SML) | 0.293 | \$2.64 | -2.75 |

The diversion ratios (the proportion of lost orders flowing to each option) after A raises its price by one dollar:

| Estimator | $\to$ outside | $\to$ B (same nest) | $\to$ C | $\to$ D |
|---|---|---|---|---|
| Truth (DGP) | 0.238 | 0.342 | 0.188 | 0.231 |
| plain logit (naive) | 0.315 | 0.256 | 0.194 | 0.235 |
| logit + control function | 0.320 | 0.254 | 0.192 | 0.234 |
| nested logit | 0.291 | 0.321 | 0.175 | 0.213 |
| mixed logit (panel SML) | 0.250 | 0.328 | 0.191 | 0.231 |

The two tables together tell the whole through-line of the chapter. The first table looks at price sensitivity: the naive $\alpha = 0.141$ severely understates, the control function pushes it to $0.245$, and mixed logit further to $0.293$, almost hitting the truth of $0.30$. The second table looks at the substitution pattern, which is the deeper teaching point: the diversions of the naive and control-function versions of the logit are almost identical, both sending only around $0.25$ to the same-nest B and both sending too much, up to $0.32$, to the outside option, which shows that curing endogeneity does not touch IIA at all, and that the two problems are independent of each other. What truly pulls the diversion toward B from $0.25$ up to over $0.32$, close to the truth of $0.34$, is nested logit and mixed logit; and among them mixed logit gets even the $0.250$ going to the outside option closest to the truth of $0.238$.

![Mixed logit recovers the price-sensitivity heterogeneity that plain logit flattened out. The true distribution (navy solid line) N(0.30, 0.12) and the mixed logit estimated distribution (red dashed line) nearly coincide; plain logit gives only a single point compressed by endogeneity (gray), and the control function pushes this single point to the right (blue), but it is still a single point, blind to the spread within the population.](assets/fig/fig_19_heterogeneity.svg)

The reconciliation of this section can be summarized as follows: on the same data, plain logit is off by a third on price sensitivity and locked by IIA on the substitution pattern; the control function fixes the price coefficient but does not touch substitution; nested logit and mixed logit fix substitution, and among them mixed logit additionally recovers price sensitivity together with its population heterogeneity, with elasticities and diversions all close to the truth. The differences come not from sampling noise but entirely from the model specification.

## 7 Failure modes of identification and robustness

In the simulation the identification assumptions are constructed, but in real research they can fail at any moment. This section lays out the most common ways they fail and the actionable responses.

Misspecifying IIA is the most expensive error in discrete choice, because it is silent. Plain logit will always converge and always give a pretty price coefficient, you will not receive any error message, yet the substitution pattern it hands back is wrong, and the decisions of delisting, merging, and new-product introduction all ride on the substitution pattern. One diagnostic is the Hausman-McFadden test: re-estimate from a subset of the options, and if the $\hat\beta$ differs significantly from the full set, IIA is rejected. A more practical attitude is to treat IIA as a restriction to be relaxed by default, and unless there is strong reason to believe the options are genuinely independent of one another, go to nested or mixed logit and draw out the substitution pattern for the reader, rather than reporting only the average effect.

The nest structure of nested logit is a specification choice by the researcher, and drawing the nests wrong will steer the substitution pattern in a wrong direction. Along what dimension the nests should be split (cuisine, price tier, on- or off-platform) has no uniform answer, and depends on which dimension the close-neighbor relationships are strongest along. The robust practice is to try several reasonable nest structures and see whether the conclusion is stable, and if different partitions give substantively different substitution patterns, then the data themselves are insufficient to support nested logit, and one should turn to mixed logit, which needs no preset nests. Beware one common error: the within-nest share terms like $\log(s_{j\mid g})$ in nested logit are themselves endogenous, and forgetting to find a separate instrument for them will, through attenuation, push $\lambda$ toward 1, pretending that IIA holds.

The soft spots of mixed logit are at both the identification and computation ends. On identification, the variance parameters must be pinned by rich choice-set and characteristic variation; the repeated choices of a panel additionally provide the correlation of the same person across situations, which usually strengthens identification markedly. Only with single choices and choice sets lacking sufficient variation are the variance terms prone to weak identification and collapse toward zero. On specification, the shape of the taste distribution is chosen by the researcher, and setting the price coefficient to normal will allow a fraction of people to have a positive price coefficient (the more expensive the more they love to buy), which is unreasonable, so one often switches to lognormal to constrain it to be negative, but the tail of the lognormal may in turn give absurd elasticities, and there is no free lunch. On computation, the usual SML retains a simulation bias under a fixed number of draws, so let the number of draws grow with the sample, or use quasi-Monte Carlo to reduce variance, and confirm from multiple starting values that it is not a local solution. An honest paper reports the robustness of the conclusion to changing the distributional assumption and the number of draws, rather than reporting only one good-looking set of estimates.

Handling price endogeneity has its own failure modes. The control function depends on the first-step price equation being correctly specified, an assumption that cannot be tested directly, and if it is misspecified the second step is biased all the same. The validity of the instrument is likewise an assertion at the level of identification, not something a statistical test can prove, and the local cost must truly enter pricing only and not the user's preference for a given order; if a peak season both raises price and boosts the platform's overall popularity, it enters both sides at once and the exclusion restriction breaks. A weak instrument makes the variance of the corrected coefficient explode, and it may even be biased worse than without correction, and a first-step partial-F that is too low is a danger signal. Selva's instrument is extremely strong (partial-F 2351), and real data are often not so lucky, so report the first-step strength honestly.

There are two more easily overlooked checkpoints. One is the definition of the outside option, which determines the "market size," and hence the level of all shares and elasticities; define the market too large or too small, and the outside-option share is distorted accordingly, so the price elasticity is systematically biased. The other is the way income or budget enters, for if price enters utility through $\mathrm{price}/Y_i$ or $\log(Y_i - \mathrm{price})$, the implied substitution and welfare meaning are very different from linear entry, and choosing which requires a mechanistic reason.

Stringing these failure modes together, the credibility of discrete choice comes down, in the end, to two things: whether the substitution structure is specified right, and whether exogenous variation in price can be found. IIA tests, multiple nest structures, robustness of the taste distribution, and the strength and exclusion of the instrument are all evidence provided around these two, and none can be replaced by "the model converged and the coefficients are significant."

## 8 Further reading

::: {.readings}
Required reading, in suggested reading order:

- Train (2009, Discrete Choice Methods with Simulation, 2nd ed.). The standard textbook of discrete choice, covering logit, GEV, mixed logit, and simulation estimation in one volume, the authoritative source for all the methods in this chapter; read it first to establish the whole picture.
- McFadden (1974). The origin of conditional logit and RUM; read it to understand where IIA comes from and why it is the restriction to be loosened.
- McFadden and Train (2000, Journal of Applied Econometrics). The universality theorem that mixed logit can approximate any RUM; understand why random coefficients is the most general tool.
- Petrin and Train (2010, Journal of Marketing Research). The standard practice of handling endogeneity with a control function on individual data, the direct source of this chapter's Section 4.2.
- Gentzkow (2007, American Economic Review). An exemplar of discrete choice demand in the context of information goods, demonstrating why the substitution and complementarity structure is the key product.

Further reading:

- Berry (1994, RAND Journal of Economics). The origin of share inversion, connecting the individual choice model to data with only market-level shares, taken up later in this series when we discuss market-level demand estimation.
- Nevo (2001, Econometrica). A practical guide to random coefficients logit on market-level data, a mirror image of this chapter's individual version.
- Bhat (2001, Transportation Research Part B). Quasi-Monte Carlo (Halton sequence) simulation for mixed logit, a practical technique for reducing SML variance.
- Small and Rosen (1981, Econometrica). Welfare measurement for discrete choice models, the bridge from log-sum to changes in consumer surplus.
- Conlon and Mortimer (2021, RAND Journal of Economics). A modern treatment and aggregation of diversion ratios, a systematic exposition of the diversion ratio as a core product of demand.
:::

::: {.apa-refs}
- Berry, S. T. (1994). Estimating discrete-choice models of product differentiation. *The RAND Journal of Economics, 25*(2), 242-262. https://doi.org/10.2307/2555829
- Bhat, C. R. (2001). Quasi-random maximum simulated likelihood estimation of the mixed multinomial logit model. *Transportation Research Part B: Methodological, 35*(7), 677-693. https://doi.org/10.1016/S0191-2615(00)00014-X
- Conlon, C., & Mortimer, J. H. (2021). Empirical properties of diversion ratios. *The RAND Journal of Economics, 52*(4), 693-726. https://doi.org/10.1111/1756-2171.12388
- Gentzkow, M. (2007). Valuing new goods in a model with complementarity: Online newspapers. *American Economic Review, 97*(3), 713-744. https://doi.org/10.1257/aer.97.3.713
- McFadden, D. (1974). Conditional logit analysis of qualitative choice behavior. In P. Zarembka (Ed.), *Frontiers in econometrics* (pp. 105-142). Academic Press.
- McFadden, D., & Train, K. (2000). Mixed MNL models for discrete response. *Journal of Applied Econometrics, 15*(5), 447-470. <https://doi.org/10.1002/1099-1255(200009/10)15:5%3C447::AID-JAE570%3E3.0.CO;2-1>
- Nevo, A. (2001). Measuring market power in the ready-to-eat cereal industry. *Econometrica, 69*(2), 307-342. https://doi.org/10.1111/1468-0262.00194
- Petrin, A., & Train, K. (2010). A control function approach to endogeneity in consumer choice models. *Journal of Marketing Research, 47*(1), 3-13. https://doi.org/10.1509/jmkr.47.1.3
- Small, K. A., & Rosen, H. S. (1981). Applied welfare economics with discrete choice models. *Econometrica, 49*(1), 105-130. https://doi.org/10.2307/1911129
- Train, K. E. (2009). *Discrete choice methods with simulation* (2nd ed.). Cambridge University Press. https://doi.org/10.1017/CBO9780511805271
:::
