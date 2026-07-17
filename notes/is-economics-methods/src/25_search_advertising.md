---
title: "Consumer Search, Advertising & Information Goods"
subtitle: "Search Frictions, Rankings, and the Demand You Cannot See"
seriesline: "Foundations of Information Systems Economics · Chapter 25"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 25 · Consumer Search, Advertising & Information Goods"
---

## Introduction

Online shelves are nearly infinite and a price comparison takes only a few clicks, yet the same product can still carry very different prices; the items on the first page of search results capture most of the sales, while sellers further down are willing to pay for a spot near the top. If consumers really saw and compared every option, these facts would be hard to explain. A more ordinary fact is that people get tired: each additional product costs time and attention, and most searches stop once things are "already good enough." A very small search cost is enough to rewrite the conclusions of full-information competition wholesale.

For a platform, ranking, recommendation, and advertising all allocate scarce attention. For a researcher, this also creates a selection problem: a purchased product may fit preferences better, or it may simply have been seen first; a highly rated product may be higher quality, or it may just have attracted more reviews; a post-recommendation rise in sales may come from discovery, from persuasion, or from crowding out a purchase that would have happened anyway. Treating the final choice as an optimization over the complete choice set erases the search process from the data, and along with it misestimates price elasticities, substitution patterns, and the returns to advertising.

This chapter starts from Stigler's fixed-sample search, McCall-Weitzman sequential search, and the reservation-value rule, and explains the Diamond paradox and how search produces price dispersion. It then uses randomized rankings and platform experiments to discuss the identification of search costs, position effects, and advertising effects, and compares the likelihoods of fixed-sample and sequential search. Reviews and recommender systems, as extensions of the same attention mechanism, are folded in as well: when reviews are selectively generated, how reputation feeds back into later transactions, and whether recommendation expands discovery or raises sales concentration all require separating exposure, consideration, and choice.

Search costs, then, cannot be named merely by a well-fitting click model; they must be identified from clean variation. Whether advertising conveys information or shifts preferences, whether reviews reflect quality or selective expression, and whether recommendation expands discovery or merely reorders exposure all demand different experimental or quasi-experimental designs. Platform logs record what users clicked, but rarely record directly why they stopped; the work of search economics is precisely to recover the neglected choice process from the point where the user stopped.

## 1 A Number That Does Not Add Up

::: {.case}
Lumen is a fictional product-search platform: users enter a query and see a page listing several items. We tell this case with simulated data, whose advantage is that the data-generating process is fully known, so every estimator's behavior can be reconciled against the truth, a teaching condition that real data can never provide. The setup is as follows: 8000 consumers, each facing the same 8 products, but the platform displays them in a random order. Consumers inspect them one by one from the top down, paying a search cost for each item inspected, and only after inspecting one do they learn the utility it delivers (price and match are both revealed upon inspection); they stop and buy once they see something satisfactory, otherwise they keep scrolling or ultimately buy nothing. The platform's problem is to estimate demand: how price-sensitive consumers are, so that it can optimize pricing and ranking.

Because it is simulated, we know the answer. The true price sensitivity is $\alpha = 0.030$ (logit utility units, per dollar), the true search cost is $c = 0.40$ utility units, which converts to roughly 13 to 14 dollars per item inspected.
:::

The platform's data scientist runs the standard move: treat each buyer's choice as a discrete choice over all 8 products on the page, and estimate the price coefficient with a conditional logit. This is a full-information demand model, and it implicitly assumes the consumer has seen and compared all 8 options. The result is a price coefficient $\hat\alpha_{\text{naive}} = 0.0257$ (SE 0.0007). Taken in isolation, this is a decent demand estimate, ready to be used for elasticities, pricing, and ranking optimization.

The trouble lies in an implicit assumption of this model: it assumes the consumer has seen all 8 products. But the search log flatly contradicts this in black and white. In the data, consumers inspected on average only 1.717 of the 8 items, and 72.5% stopped searching after opening just the first one. In other words, the full-information model implicitly assumes zero search cost, while the true behavior corresponds to a search cost of about 13 dollars per inspection. A demand model that treats the search cost as zero, applied to a market where the search cost is plainly high, can you still trust the price coefficient it estimates?

The answer is no, at least not directly. $\hat\alpha_{\text{naive}} = 0.0257$ is about 14% below the truth of 0.030: price sensitivity is systematically underestimated. Section 3 will explain why: the full-information model treats the cheap, lower-ranked products the consumer never inspected as "seen but rejected," so the signal that price sends about choice is diluted and the elasticity is compressed. Based on this attenuated elasticity, the platform will set prices too high. Section 6 will use a structural search model to recover the true 0.030, and at the same time estimate the 13-dollar search cost that the full-information model erased. The through-line of this chapter is how to correctly reconstruct demand when you cannot see it in full.

## 2 The Economic Model and the Estimand

Before setting out to estimate, we need to understand how costly search shapes a market and exactly which quantities we want to estimate. This section first lays out the classic theory of search, which defines the objects of estimation and also explains the value of ranking and advertising.

### 2.1 Fixed-Sample Search and Sequential Search

Search theory begins with a plain question from Stigler (1961): a consumer faces a set of sellers whose posted prices are unknown, and each quote costs $c$ to obtain; how many should he ask? In fixed-sample search, the consumer decides in advance to draw $n$ quotes and then buys the cheapest among them. The benefit of asking one more is the chance of hitting a lower price; the cost is paying another $c$. The optimal $n^\ast$ falls where marginal benefit equals marginal cost: the amount by which the $n$-th search lowers the expected minimum price further is exactly $c$. As $c$ rises, $n^\ast$ falls, and people search less.

Fixed-sample search has an unnatural feature: it requires the consumer to fix how many sellers to search up front, drawing them all even if the first one already quotes an extremely low price. More intuitive is sequential search: look one at a time and, after each, decide whether to continue. McCall (1970) proved that the optimal strategy for sequential search is a reservation rule. There exists a reservation value $z$: anything better than it stops and is accepted, anything worse continues the search. $z$ is pinned down by an indifference condition: the expected gain from searching once more is exactly equal to the cost $c$. In the language of utility, if the utility a further search could deliver is denoted $u$, with distribution $F$, then

$$c = \mathbb{E}\big[\max(u - z, 0)\big] = \int_{z}^{\infty} (u - z)\, dF(u)$$

Read this equation as follows: the right side is the expected improvement from "inspecting one more item, switching to it if it beats the $z$ in hand and otherwise not," the left side is the cost of that one inspection, and where the two are equal the consumer is exactly indifferent about searching, so $z$ is the critical threshold. The higher the search cost $c$, the lower is $z$, and the sooner people are satisfied and the shallower they search. This reservation-value equation is the core mechanism of this chapter's running case: in Lumen, each consumer decides how many items to inspect before stopping precisely by this rule.

### 2.2 Weitzman's Pandora Problem and the Diamond Paradox

McCall's rule addresses homogeneous objects of search. Real search involves options that differ ex ante, some of which look more promising. Weitzman (1979)'s Pandora's boxes model handles this heterogeneous case and delivers a beautiful three-part optimal rule. Each box $j$ has a reservation value (also called an index) $z_j$, pinned down by its own cost equation:

$$c_j = \mathbb{E}\big[\max(u_j - z_j, 0)\big]$$

Given the boxes' $z_j$, the optimal strategy is: a selection rule, open boxes in decreasing order of $z_j$; a stopping rule, stop when the maximum realized utility among opened boxes exceeds the $z_j$ of all unopened boxes; and a choice rule, ultimately take the opened box with the highest utility. The depth of this rule is that it reduces a seemingly complex dynamic-programming problem to a set of static indices: each box's value of opening can be computed on its own, independently of the other boxes.

Costly search carries a striking equilibrium consequence, the Diamond (1971) paradox. Imagine many perfectly identical sellers selling the same thing, consumers searching sequentially, each quote costing $c > 0$. The claim is: the unique equilibrium price is the monopoly price, not the competitive price. The argument takes a single step. Suppose all sellers set some price $p$. Then no consumer is willing to search, because all sellers look alike and searching one more just wastes $c$. Since no one would leave over one seller raising its price slightly, any seller can raise its price a little (by less than $c$) without losing any customers. This logic pushes all the way up to the monopoly price before it stops. So even if the search cost is almost zero, as long as it is positive, Bertrand competition collapses entirely. The Diamond paradox is the intellectual origin of search economics: it shows the full weight of the single assumption that "search is costly," and it explains why price dispersion and markups are so common in real markets.

### 2.3 Equilibrium Price Dispersion

The Diamond paradox predicts a single monopoly price, yet what we see in reality is price dispersion: the same thing carries different prices from different sellers. Burdett and Judd (1983) show how price dispersion can arise as an equilibrium in a perfectly homogeneous, rational market. The key is that consumers are not uniform in search intensity: some are non-shoppers who see only one quote, others are shoppers who see two or more. Facing this mixed demand, sellers are caught in a dilemma: a high price earns more from non-shoppers, a low price wins the comparison among shoppers. In equilibrium there is no single optimal price; sellers use a mixed strategy, randomizing over a price interval, with every price in the interval yielding the seller equal expected profit. This indifference condition pins down the equilibrium price distribution $F(p)$. Price dispersion is therefore not a residue of irrationality or friction, but an equilibrium product of search costs together with heterogeneity in search intensity. This point is extremely important for empirical work, because it means the observed price distribution itself carries information about search costs, and this is exactly the footing for the Hong-Shum method in Section 4.

### 2.4 Advertising: Information or Persuasion

Advertising is a close neighbor of search, because both concern how consumers come to know products. Economics holds three views of advertising. The informative view treats advertising as information: it lowers search costs or directly puts the product into the consumer's consideration set, letting products that would not otherwise be seen be seen. The persuasive view treats advertising as a shift in preferences: it directly raises the consumer's willingness to pay for the product, independent of any information. The complementary view treats advertising itself as a consumption good that complements the product. The three views have sharply different policy and welfare implications: informative advertising improves matching and raises welfare, while persuasive advertising may merely transfer share or even distort choice. Distinguishing them is the core difficulty of advertising empirics; Ackerberg (2003)'s classic approach uses the contrast between first-time buyers and repeat customers to separate the informative effect from the persuasive effect (Yoplait 150 yogurt, see the further reading in Section 8).

In the platform setting, advertising and ranking are highly isomorphic. Placing a product at the top of a search results page is essentially a form of informative advertising: it raises the probability that the product is inspected, without necessarily changing the consumer's evaluation once she sees it. This distinction is crucial, and it is a key property of this chapter's running case: a randomized display position affects only which products get searched, not the utility a consumer assigns after inspection. Section 3 will show that this property is precisely the lever for identification.

### 2.5 Experience Goods and Bayesian Learning

There is a class of products whose quality is only learned through use, called experience goods. For these, the consumer chooses under uncertainty and continually updates her beliefs about quality through use and through advertising. Erdem and Keane (1996) characterize this process with a Bayesian learning model. The consumer is uncertain about the true quality $Q_j$ of product $j$, with a normal prior $\mathcal{N}(\mu_0, \sigma_0^2)$. Each use or exposure to advertising delivers a noisy signal $s = Q_j + \nu$, $\nu \sim \mathcal{N}(0, \sigma_\nu^2)$. Bayesian updating (the normal-normal case is the Kalman filter) gives the posterior, whose precision is the sum of the prior precision and the signal precision, and whose mean is the precision-weighted average of the two:

$$\frac{1}{\sigma_1^2} = \frac{1}{\sigma_0^2} + \frac{1}{\sigma_\nu^2}, \qquad \mu_1 = \sigma_1^2\left(\frac{\mu_0}{\sigma_0^2} + \frac{s}{\sigma_\nu^2}\right)$$

If the consumer is risk-averse, utility depends not only on the posterior expected quality $\mu$ but is also discounted by the posterior variance $\sigma^2$, so uncertainty itself lowers the choice probability. Advertising in this framework is a signal source: it affects choice by lowering the posterior variance, a purely informational mechanism, quite distinct from persuasion. Experience-good learning generalizes the static information acquisition of search to dynamic, intertemporal accumulation of information, and connects to the framework used later in this series on dynamic structural models; here it is enough to know that it formalizes "how a consumer's beliefs about quality evolve with information."

### 2.6 Reviews, Reputation, and Recommender Systems

Reviews are a mechanism for aggregating scattered use experiences, but they are by no means a random sample of quality. The average rating a product ends up displaying is filtered jointly by who bought it, who bothered to write, which reviews the platform chooses to place up front, and whether the merchant padded the count. Two kinds of selection are especially stubborn. One is reporting selection: people with especially good or especially bad experiences have more incentive to write, while a large mass of lukewarm middling experiences stays silent, so ratings are naturally biased toward the extremes. The other is reputation feedback: a few early good reviews push up sales, sales bring more reviews, and reviews further lift ranking and sales, snowballing. For this reason the positive correlation between rating and sales mixes three forces, quality, selection, and ranking, and cannot be read directly as "good reviews caused sales" in a causal sense.

Recommender systems, meanwhile, alter both the consumer's information set and demand at once. Recommending popular items may raise clicks in the short run, but through a feedback loop it concentrates exposure further onto a few products; personalized recommendation may improve match efficiency, or it may narrow the categories a consumer is exposed to more and more. Its welfare consequences must be decomposed into at least three parts: discovery, whether the consumer saw a product she could not have seen but that fits her better; persuasion, whether preferences were pushed to change by the interface itself; and displacement, whether the sales brought by recommendation were merely transferred from other products rather than growing the pie. Blend all three into one number, "recommendation raised sales," and the direction can be entirely different.

On identification, the cleanest design is to randomize review visibility, ranking, or recommendation eligibility, and to observe exposure, clicks, purchases, and later ratings at the same time. Comparing purchase rates only among recommended products introduces post-treatment selection; looking only at the mean star rating ignores review count, text content, and selective deletion. For platform governance, one should also treat fake reviews, retaliation, and professional reviewer incentives as mechanisms, rather than lumping them all under measurement noise.

### 2.7 The Estimand

Let us bring the theory down to the objects of estimation. This chapter has two core quantities to estimate: one is price sensitivity $\alpha$ (and, more generally, the preference parameters), the other is the search cost $c$. The two are entangled in observed data: an observed purchase is the product of the search process, and a consumer's failure to buy some cheap product may be because she dislikes it (low preference) or because she never searched it up at all (high search cost). Separating these two explanations is the goal of this chapter's identification. The search cost $c$ converted to money is $c / \alpha$, that is, how many dollars each additional item inspected is worth, an interpretable and comparable quantity. In addition, the value the platform cares about, that of a ranking slot or ad position, can be characterized as the causal effect on a product's sales of moving it from the bottom to the top, a quantity that can be cleanly identified under randomized ranking, as Section 3 develops.

To summarize this section: costly search, through the reservation-value rule, determines how deep a consumer looks and what she buys, and it can generate price dispersion and supra-competitive markups even in a homogeneous market; the objects of estimation in this chapter are price sensitivity $\alpha$ and search cost $c$, which are entangled in observed choices and must be separated; advertising and ranking, as information mechanisms, have value that can be characterized as causal effects on inspection probability and sales.

## 3 Identification: How Randomized Ranking Separates Search Cost from Preferences

Identification logically precedes estimation. The identification question here is: under what conditions do the observed search-and-purchase data separate the search cost $c$ from preferences $\alpha$? This is the most detailed section of the chapter, broken into numbered subsections that dissect the logic layer by layer.

### 3.1 The Gulf Between Limited Consideration and Full Information

Let us make the core of the problem clear first. The full-information demand model assumes each consumer chooses over the entire choice set, and that the observed purchase is the result of an all-in comparison. The world of costly search is not like this: the consumer inspected only a subset, the consideration set, and the purchase is a choice made within this subset, with products outside it never entering the decision.

The difference between these two worlds in choice data, applied to Lumen, is glaring. The full-information model assumes the consumer compared 8 products, yet the log shows she inspected on average only 1.717. This means that for a cheap product ranked lower that the consumer never opened, the full-information model reads "did not buy it" as "saw it, compared it, rejected it," when the truth is that this product was never seen. The former reading attributes the non-purchase to low preference (say, insensitivity to price, so the low price failed to attract), the latter attributes it to not having searched it up. Confusing the two is exactly the source of the naive estimate's error.

::: {.intuition}
Here is an analogy. You are in a large supermarket, and you browse only the two or three shelves near the entrance before grabbing a bottle of soy sauce and checking out. An analyst who sees only what you ultimately bought, but assumes you walked the whole store and compared every soy sauce, will conclude that "you are not interested in that cheaper bottle in the far corner." But the fact is you never walked to that corner. He has mistaken "you did not reach it" for "you reached it but did not like it." To correct this error, he needs to know which shelves you actually browsed, that is, your search path. The search log provides exactly this path, and the randomized shelf placement makes that path usable for identification.
:::

### 3.2 Why the Naive Full-Information Estimate Attenuates

Now let us explain where the 14% attenuation of Section 1 comes from. The naive conditional logit sets the choice set to all 8 products, including a large number that were never inspected. The price coefficient is identified from "the association between price and being chosen," and this association is diluted: those lower-ranked, cheaper products that should have strongly attracted price-sensitive consumers were not chosen because they were never searched up, so their low prices and non-purchases occur together, which in the full-information model becomes evidence that "even a low price failed to draw a purchase," pushing down the estimated price sensitivity.

The direction is definite: attenuation toward zero. The more unseen cheap products there are, the more heavily the price signal is diluted, and the more downward-biased is $\hat\alpha_{\text{naive}}$. In Lumen, consumers inspect only 1.717 products, a large mass of cheap options never enters the consideration set, so the attenuation is pronounced, 0.0257 against a truth of 0.030. This direction is exactly opposite to the upward-biased state dependence of Chapter 24, but the pathology is the same: both miscount the trace of an unmodeled mechanism (persistent heterogeneity there, limited consideration here) onto the target coefficient.

::: {.warning}
The bias of a full-information demand model in a search market cannot be rescued by a larger sample. However large the sample, as long as the model keeps assuming the consumer saw all options, the price coefficient keeps attenuating, because this is misspecification, not sampling noise. A data scientist who runs a full-information logit on a platform whose search log plainly shows shallow search will get systematically low elasticities, and will set systematically high prices on that basis. When you see that the average number inspected in the search log is far below the number available, the full-information model should be replaced.
:::

### 3.3 Randomized Ranking as the Lever for Identification

To separate the search cost from preferences, we need a lever that exogenously changes "who gets searched up." The platform's randomized ranking provides exactly that. When the platform displays the 8 products in a random order, a product's position has nothing to do with its own price or quality, being purely random. So which products are more likely to be inspected has an exogenous source that does not depend on product characteristics. For this lever to hold, it rests on two assumptions:

::: {.assumption}
**(Exogenous ranking + exclusion restriction)** The display position $r_{ij}$ is randomly assigned, independent of the product's price and quality and of the consumer's match value $\varepsilon_{ij}$; and position affects purchase only through affecting "whether product $j$ is inspected," and does not directly enter the post-inspection utility $u_{ij}$. The first condition makes the opportunity to be inspected exogenous, the second (the exclusion restriction) guarantees that position acts entirely through the single channel of search.
:::

Under these two assumptions, randomized ranking identifies two things.

Randomized ranking identifies two things. The first is the search cost. How deep a consumer searches directly reflects the level of the search cost: when the search cost is high, people inspect little and stop. Given the random display order, the observed average search depth (1.717 in Lumen) pins down $c$, because in the structural model it is $c$ that, through the reservation value $z$, determines on average how far people look before being satisfied. The second is preferences. Because which products enter a given consumer's consideration set is randomly determined (driven by random position), the composition of the consideration set is unrelated to product price, so within the consideration set the association between price and purchase is no longer contaminated by search selection and can be used to identify $\alpha$ cleanly. These two points together are why randomized ranking is the key to identification.

Ursu (2018) nails down this logic with a real randomized-ranking experiment, and also reveals a subtle property. She finds that ranking strongly affects which products consumers search, but conditional on being searched, position has no effect on purchase. In other words, products near the top sell more entirely because they are inspected by more people, not because being near the top itself changes how people evaluate them once seen. Position works through search, not through persuasion. This chapter's running case is built to carry this property: random position affects only which products are inspected and does not enter post-inspection utility, so conditional on being inspected, position is unrelated to purchase. This is both a credible behavioral assumption and a clear economic interpretation of the value of ranking: it sells exposure, not preference.

### 3.4 Why a Seemingly Reasonable Patch Overcorrects

Knowing the importance of the consideration set, a natural patch is: since the full-information model has the wrong choice set, replace the choice set with the products the consumer actually inspected, and run a conditional logit on the inspected set. This idea grabs the right direction, but it overcorrects.

The problem is that the inspected set is itself a product of the search process, and the search process depends endogenously on realized utility. A consumer inspects few products often because she ran into a cheap, high-match product early, satisfied the reservation value, and stopped; and a consumer who inspects many products usually did so because the earlier ones were all unappealing. So a correlation, manufactured by the stopping rule, arises between "the set of inspected products" and "price": a small consideration set is more likely to contain the cheap, good product that made the consumer stop. Running a logit on the inspected set counts this price-purchase correlation manufactured by the stopping rule into the price coefficient, leading to overestimation. In Lumen this inspected-set estimate gives 0.0344, overshooting the truth of 0.030, in the opposite direction from the naive estimate.

The lesson here is worth remembering: randomized ranking guarantees that the opportunity for which products are displayed, and therefore inspected, is exogenous, but when a consumer stops searching is endogenous. The former can be used for identification, the latter introduces selection if not modeled. The truly clean approach is not to run a regression on some observed subset, but to estimate the stopping rule itself explicitly as part of the model, which is the structural search estimation of the next section.

### 3.5 Sequential or Fixed-Sample

The identification logic above defaulted to consumers searching by a sequential rule. This assumption is itself testable, and the test may not turn out as expected. De los Santos, Hortaçsu and Wildenbeest (2012) test the sequential search model with data on web browsing plus purchasing, and find that the data reject it, favoring instead fixed-sample (non-sequential) search. Their evidence comes from a characteristic prediction of sequential search: if consumers search sequentially, they should stop upon seeing a good enough quote, so behavior like "searched many sellers, then went back and bought one seen earlier" should not appear. Yet such revisits and the absence of a "stop while ahead" pattern are common in the data, looking more like consumers deciding in advance how many to look at and then picking the best among them.

This finding has a direct implication for empirical work: sequential versus fixed-sample is not a technical detail to be chosen freely; they correspond to different behavioral assumptions and different likelihoods, and using the wrong one brings specification bias. This chapter's running case adopts sequential search, because it fits best with the narrative of the reservation-value rule and randomized ranking, but on real data the search protocol itself should be treated as an object to be tested rather than assumed. Section 7 returns to this robustness issue.

### 3.6 Identifying Advertising Effects

Identifying advertising effects has its own trap, and its core is selection. What products a platform gives advertising or a high rank to is usually not random; it goes to products already more likely to sell well, or ramps up during periods when demand is already strong. So the correlation "advertised products sell well" mixes in the reverse causality that "products that sell well anyway are the ones that get advertised." Without handling this layer of selection, any estimate of an advertising effect will be exaggerated.

Cleanly identifying advertising effects relies on exogenous variation. Randomized display position is one kind: it makes "being seen by more people" happen exogenously, thereby separating the causal effect of exposure from product quality. A field experiment is another: randomly advertise to a fraction of users and not to others, then compare the two groups. Sections 5 and 7 will show that when researchers actually run such randomized experiments, they often find the causal effect of advertising far smaller than a correlational analysis implies, because a large chunk of that correlation is in fact selection.

The main points of this section can be summarized as follows: in a search market, an observed non-purchase confounds low preference with not having been searched up, so the full-information model attenuates price sensitivity toward zero; randomized ranking provides exogenous variation, identifying the search cost from search depth and, by making the composition of the consideration set independent of price, identifying preferences; but stopping is endogenous, so simply regressing on the inspected subset overcorrects, and the correct approach is to model the search process explicitly; the distinction between sequential and fixed-sample is testable and may not turn out as expected, while identifying advertising effects relies on randomized experiments or randomized ranking to strip out selection.

## 4 Estimation: From Full Information to Structural Search

This section advances by the logic of "where the previous method fails gives rise to the next." The goal is to estimate price sensitivity $\alpha$ and search cost $c$ cleanly.

### 4.1 Naive Full-Information Conditional Logit

The easiest move is to treat each buyer's choice as a discrete choice over all 8 products on the page, and estimate with a conditional logit:

$$\Pr(\text{buy } j \mid \text{choice set} = \text{all } J) = \frac{\exp(-\alpha p_j + \xi_j)}{\sum_{k=1}^{J} \exp(-\alpha p_k + \xi_k)}$$

Its failure has already been laid out in Section 3: the choice set is stuffed with a large number of never-inspected products, the price signal is diluted, and $\hat\alpha$ attenuates toward zero. In Lumen it delivers 0.0257, 14% below the truth of 0.030. Its only use is as a benchmark for the attenuation, a hint that the full-information assumption cannot stand in this market.

### 4.2 Inspected-Set Conditional Logit

The second patch is to replace the choice set with the products the consumer actually inspected:

$$\Pr(\text{buy } j \mid \text{choice set} = \text{inspected set}_i)$$

It grabs the correct direction that "choice happens within the consideration set," but because of the endogenous stopping of Section 3.4 it overcorrects. In Lumen it gives 0.0344, overshooting the truth. Its significance is diagnostic: place the naive 0.0257 and the inspected-set 0.0344 side by side, and the truth is bracketed between them, which itself shows that the consideration set can be taken neither as everything (attenuation) nor simply as the inspected subset (overestimation), and the search process itself must be modeled.

### 4.3 Structural Sequential Search Estimation

The real solution is to estimate the entire search-and-purchase process as a structural model. Given parameters $(\gamma, \alpha, c)$, the Weitzman rule fully determines how deep each consumer looks and what she ultimately buys. Conversely, the observed distribution of search depth, the purchase rate, the prices of the products bought, and similar statistics carry information about these parameters. The estimation idea is the simulated method of moments: under candidate parameters, simulate the entire search process, compute a handful of simulated statistics, and tune the parameters until they match the same-named statistics in the data.

Concretely, pick a set of moments that identify the parameters: the average search depth (mainly determined by $c$, the higher the search cost the shallower the search), the purchase rate (determined by the products' attractiveness relative to the outside option, pinning down the utility level), the average price of the products bought (determined by $\alpha$, the more price-sensitive the cheaper the purchase), and the share who stop after inspecting just one. Under each candidate $(\gamma, \alpha, c)$, solve the reservation-value equation for $z$, simulate the search and purchase of $N$ consumers, compute these four moments, and minimize the weighted distance to the data moments. Because the search process is fully modeled, the endogeneity of stopping is handled correctly and the estimate is consistent. To keep simulation noise from making the objective function rugged, the simulation uses fixed random draws (common random numbers), and optimization starts from multiple points to avoid local minima.

In Lumen this structural estimate gives $\hat\alpha = 0.0298$, $\hat c = 0.409$, $\hat\gamma = 2.900$, all three landing almost dead center on the truth (0.030, 0.40, 2.90). The search cost converted to money is $\hat c / \hat\alpha = 13.71$ dollars, that is, each additional item inspected is worth about 13.7 dollars. The naive model treated this search cost as zero, the structural model recovers it, and at the same time corrects price sensitivity from the attenuated 0.0257 to 0.0298.

### 4.4 Estimating Search Costs from the Price Distribution: Hong-Shum

Structural sequential search estimation needs a search log, that is, knowing which products each consumer inspected. Such data are not always available. Hong and Shum (2006) give a method to estimate search costs using only the observed price distribution, requiring no quantity or search data, and its footing is exactly the equilibrium price dispersion of Section 2.3. The logic runs like this: in equilibrium, a seller's pricing over the interval must yield equal expected profit at every price, and this indifference condition links the price distribution $F(p)$ to the distribution of consumers' search costs. Sellers who set high prices earn from consumers who search little and compare little, those who set low prices earn from those who search a lot. Given the observed $F(p)$, inverting this relationship recovers the distribution of search costs, with the fixed-sample and sequential versions each having their own inversion formula. Using several graduate textbooks priced across different online bookstores, they estimate that the median search cost is modest in magnitude, and that the non-sequential and sequential versions differ substantially (about 2 dollars versus about 29 dollars), with a sizable share of consumers behaving as if they searched only once. The value of this method is its very low data requirement, needing only prices, at the cost of relying on the equilibrium pricing model being correctly specified.

### 4.5 Estimating Advertising Effects

Estimating advertising effects, as described in Section 3.6, succeeds or fails on whether selection can be stripped out. The correlational approach (comparing the sales of advertised and non-advertised products) almost inevitably overestimates, because placement itself is selective. Credible approaches come in two kinds. One is field experiments, randomly advertising to users or not, and comparing the two groups directly; Lewis and Reiley (2014) and Blake, Nosko and Tadelis (2015) are exemplars, the latter finding that heavily purchased branded keyword advertising has an incremental effect of nearly zero, because those users would have come anyway. The other is to use the platform's own randomization, such as random ranking, to identify the causal effect of exposure. In Lumen, randomized ranking lets us read off the causal value of position directly, which is part of the next section's walkthrough.

## 5 Anchoring Papers

A method holds up only when it lands in real research. Two anchoring papers, one establishing the method of estimating search costs from price data, one using a randomized-ranking experiment to nail down the relationship among search, position, and purchase, are each laid out along five elements, paper, method, data, results, and limitations, with the focus on how the identification assumptions are defended.

### 5.1 Hong and Shum (2006, RAND Journal of Economics)

::: {.case}
Paper and place in method history: "Using price distributions to estimate search costs," RAND Journal of Economics 37(2), 257-275. This paper's contribution is to prove that search costs can be estimated from price data alone, without observing quantities or search behavior, which is extremely valuable in data-scarce markets.

Method: use the equilibrium price dispersion model of Section 2.3 in reverse. Equilibrium requires sellers to be indifferent over the support of prices, and this condition maps the observed price distribution $F(p)$ to the distribution of consumer search costs. Given $F(p)$, inverting this mapping yields the search cost distribution. The paper gives separate inversions for fixed-sample search and sequential search, which correspond to different equilibrium indifference conditions and therefore give different search-cost estimates, and this contrast is itself a probe of the search protocol.

Data: the distribution of posted prices for several graduate textbooks (such as Billingsley, Stokey-Lucas, and others) across different online bookstores. The rationale for using price data is precisely the method's selling point: it does not require observing how many sellers a consumer searched or how much she bought, only cross-sectional price dispersion.

Results: the estimated median search cost is modest in magnitude, but the sequential and fixed-sample versions give substantially different numbers, the non-sequential version about 2 dollars, the sequential version about 29 dollars. This contrast itself reminds the user that the specification of the search protocol directly drives the estimate. A notable finding is that a sizable share of consumers behave as if they searched only once, consistent with later evidence of shallow search in many search markets.

Limitations: the method relies entirely on the equilibrium pricing model being correctly specified, and price dispersion must genuinely come from search frictions, not from product heterogeneity, quality differences, or cost differences. If bookstores differ in service, delivery, or trustworthiness, price differences need not mirror search costs, and the inverted search costs will be contaminated. Moreover the results are sensitive to the search protocol (sequential or fixed-sample), and the data themselves cannot fully determine which to use.
:::

This paper's defense strategy is to turn a hard-to-observe quantity (the search cost) into a function of an observable one (the price distribution), building the bridge between them out of an equilibrium model. Its persuasiveness depends on how solid that bridge (the equilibrium pricing model) is, which also reminds the user that the lower a method's data requirement, the heavier its dependence on model specification tends to be.

### 5.2 Ursu (2018, Marketing Science)

::: {.case}
Paper: "The power of rankings: Quantifying the effect of rankings on online consumer search and purchase decisions," Marketing Science 37(4), 530-552. The question is how exactly the ranking position of search results affects consumers, and whether it works through search or through directly changing preferences, which is precisely the real-world prototype of this chapter's running case.

Method: use a real randomized-ranking field experiment. The platform randomizes the display order of hotels, so position is unrelated to a hotel's own characteristics and is exogenous. On top of this exogenous variation, estimate a sequential search model to identify separately the effect of position on "whether it is searched" and on "whether it is purchased conditional on being searched." Randomization is the key to identification: it severs the usual endogenous correlation that "good hotels are ranked near the top."

Data: Expedia's hotel search data, with experimentally randomized rankings, recording which hotels a consumer inspected and which one she ultimately booked.

Results: ranking strongly affects which hotels consumers search, and those near the top are inspected far more. But a key finding is that, conditional on a hotel being inspected, its position has no effect on whether it is booked. Position works through search, not through persuasion or a position bias. This means the commercial value of a ranking slot is essentially the value of exposure: moving a product to the top raises sales purely because it is inspected by more people, not because being near the top makes people want to buy it more.

Limitations: the conclusion depends on the specification of the sequential search model, and if real search is closer to fixed-sample, the mechanism by which position affects search would differ. Randomization solves the endogeneity of ranking, but the model still makes parametric assumptions about how consumers search and how the outside option is specified. Moreover the external validity of a single platform and a single category (hotels) requires caution, and whether position also works purely through search in other categories is an empirical question.
:::

Put the two papers together, and the meaning of anchoring becomes clear: Hong-Shum show how, with only prices available, an equilibrium model estimates search costs, at the cost of heavy dependence on model specification; Ursu shows how a randomized experiment decomposes the causal role of ranking into "affects search" and "affects purchase" channels, and cleanly determines that it runs through the former. The latter especially provides the real-world grounding for this chapter's running case, and confirms the view that ranking should be understood as exposure rather than persuasion.

## 6 A Complete Walkthrough on the Lumen Data

Now let us run the tools of Section 4 through completely on the Lumen search log. The code below uses R 4.5.3, fixing set.seed(21) to ensure reproducibility, and every number cited in the text comes from the actual run output of this code.

### 6.1 DGP

The design parameters are as follows: 8000 consumers, 8 products, randomized display order. A consumer inspecting product $j$ gets utility $u_{ij} = \gamma - \alpha p_j + \varepsilon_{ij}$, $\varepsilon \sim \text{Gumbel}(0,1)$, with price and match value revealed together upon inspection. The truth is $\gamma = 2.90$, $\alpha = 0.030$, $c = 0.40$, plus a per-person outside option $u_{i0} = 0.60 + \text{Gumbel}(0,1)$.

```r
set.seed(21)
gamma <- 2.90; alpha <- 0.030; c_true <- 0.40; o0 <- 0.60
prices <- round(rnorm(8, 100, 18), 1)
mean_u <- gamma - alpha * prices
## reservation value z: c = E[max(u - z, 0)]
z_star <- uniroot(function(z) {
  uu <- gamma - alpha * sample(prices, 3e5, TRUE) - log(-log(runif(3e5)))
  mean(pmax(uu - z, 0)) - c_true }, c(-5, 20))$root
```

The reservation value solves to $z^\ast = 0.767$. Each consumer inspects one at a time in random order, maintaining the best utility among inspected products, and stops as soon as the best opened option (including the outside option) exceeds $z^\ast$, ultimately taking the larger of the best product and the outside option.

```r
for (r in 1:8) {
  j <- order_i[r]; u <- mean_u[j] + rgum(1)
  if (u > best_pu) { best_pu <- u; best_j <- j }
  if (max(u0, best_pu) >= z_star) break     # optimal stopping
}
buy <- if (best_pu > u0) best_j else 0       # buy only if best product beats outside
```

The key statistics from the simulation: purchase rate 0.521, average inspection of 1.717 products (out of 8), 72.5% of consumers stopping after inspecting just one. This shallow search is exactly what the full-information model misreads.

### 6.2 Naive and Inspected-Set Estimates

Naive full-information conditional logit, choice set being all 8 products:

```r
m_naive <- clogit(bought ~ price + strata(cons), data = buyers_all8)
```

It gives $\hat\alpha_{\text{naive}} = 0.0257$ (SE 0.0007), 14% below the truth of 0.030, which is the attenuated number from Section 1. Switching to the inspected set:

```r
m_insp <- clogit(bought ~ price + strata(cons), data = buyers_inspected)
```

It gives $\hat\alpha_{\text{insp}} = 0.0344$ (SE 0.0013), overshooting the truth. The truth of 0.030 is bracketed between the two estimates, one low and one high, confirming Section 3.4: taking the choice set as everything attenuates, taking the inspected subset overestimates because of endogenous stopping.

### 6.3 Structural Sequential Search SMM

Model the entire search process and recover the parameters with the simulated method of moments. The moments are the average search depth, the purchase rate, the average price of the products bought, and the share who inspect just one; under candidate $(\gamma, \alpha, c)$, solve for $z$, simulate search, compute moments, minimize the weighted distance to the data moments, with multi-start optimization:

```r
sim_moments <- function(g, a, cc) {
  z <- uniroot(...)$root
  ## simulate N consumers' search + purchase under Weitzman rule
  c(mean_depth, buy_rate, mean_purchased_price, share_depth1)
}
fit <- optim(start, obj_smm, method = "Nelder-Mead")   # from several starts
```

The structural estimate gives $\hat\alpha = 0.0298$, $\hat c = 0.409$, $\hat\gamma = 2.900$, all three landing almost dead center on the truth. The data moments (depth 1.717, purchase rate 0.521, average price 91.34, inspect just one 0.725) match the simulated moments (1.714, 0.520, 90.98, 0.722) tightly. The search cost converted to money is $\hat c / \hat\alpha = 13.71$ dollars, that is, each additional item inspected is worth about 13.7 dollars. The full-information model implicitly treated it as zero, and the structural model recovers it.

![Three estimates of price sensitivity $\alpha$ against the truth (dashed line, 0.030). The naive full-information logit attenuates toward zero, the inspected-set logit overshoots the truth, and the structural sequential search SMM lands almost dead center on the truth. Error bars are 95% confidence intervals.](assets/fig/fig_25_estimates.svg)

The figure above draws the three estimates together, with the dashed line the truth of 0.030. The naive estimate is on the left (attenuated), the inspected-set on the right (overestimated), and the structural estimate lands on the dashed line.

### 6.4 Search Depth and Limited Consideration

The structural model works because, at root, consumers search very shallowly and the consideration set is far smaller than the entire choice set. Plotting the distribution of search depth makes this clearest.

![The distribution of search depth: most consumers stop after inspecting just 1 to 2 items, averaging 1.717 (dashed line), while the full-information model assumes all 8 are seen. This gulf is exactly the source of the naive estimate's attenuation.](assets/fig/fig_25_depth.svg)

The figure above shows that 72.5% of consumers inspected just one item, the vast majority stop within two or three, averaging 1.717, far below the 8 available. The full-information model assumes the nonexistent world at the far right where all 8 are seen. Limited consideration is not a small bias; it is the basic fact of this market, and any demand estimate that ignores it gets the estimand wrong from the root.

### 6.5 The Value of Position: Advertising Effects Under Randomized Ranking

Randomized ranking lets the causal value of position be read off directly. Plotting the purchase probability against display position:

![Purchase probability by randomized display position (orange line), together with the probability of being inspected at each position (gray bars). A product ranked first has a purchase probability of 0.248, one ranked eighth only 0.008, purely because those near the top are inspected more. Position works through search and does not change the post-inspection evaluation.](assets/fig/fig_25_prominence.svg)

The figure above quantifies the value of position. Because position is randomly assigned, the relationship between position and purchase is a pure causal effect, containing no selection of "good products ranked near the top." A product ranked first is purchased with probability 0.248, one ranked eighth only 0.008, the former thirty times the latter. This huge gap comes entirely from exposure: products near the top are inspected by more people, and by construction position does not enter the consumer's post-inspection utility, so conditional on being inspected, position does not change the purchase. This is exactly Ursu (2018)'s finding, that position works through search rather than persuasion. For the platform, this gap of 0.248 versus 0.008 is why a ranking slot or ad position can command a high price: it sells attention.

### 6.6 The Grand Reconciliation

Line up the estimates and view them against the truth:

| Estimator | $\hat\alpha$ | Relative to truth 0.030 | Note |
|---|---|---|---|
| Naive full-info logit | 0.0257 | -14% | Unseen cheap items treated as "seen but not chosen," attenuated |
| Inspected-set logit | 0.0344 | +15% | Endogenous stopping, consideration-set selection biases high |
| Structural search SMM | 0.0298 | -1% | Models the search process, consistent |

The structural estimate also gives the search cost $\hat c = 0.409$ (truth 0.40) and $\hat\gamma = 2.900$ (truth 2.90), with a monetized search cost of 13.71 dollars per inspection.

This section's reconciliation can be summarized as follows: from the same search log, where consumers inspect on average only 1.717 products, the full-information logit attenuates price sensitivity by 14%, the inspected-set logit overestimates by 15% due to endogenous stopping, and only the SMM that models the search process structurally recovers the true price sensitivity, the 13.7-dollar search cost, and the utility level all at once; randomized ranking additionally makes the causal value of position plain, with the purchase probability at position 1 versus position 8 differing thirtyfold, driven purely by exposure.

## 7 Failure Modes of Identification and Robustness

In the simulation the search protocol and exogenous ranking are constructed, but in real research they can fail at any moment. This section lays out the most common failure modes and operational responses.

The specification of the search protocol is the top risk. This chapter's structural estimate assumes consumers search sequentially, and if real behavior is closer to fixed-sample (consumers deciding in advance how many to look at before picking the best), the likelihood is written wrong and the parameters are biased. De los Santos, Hortaçsu and Wildenbeest (2012) remind us this is no idle worry: their web data rejected sequential search precisely. The operational response is to estimate both the sequential and the fixed-sample versions and see whether conclusions are robust, and to use their characteristic predictions for a test: sequential search predicts "stop while ahead" and very few revisits to already-seen options, while fixed-sample does not, and the prevalence of revisits in the data is one criterion.

The exogeneity of ranking cannot be taken for granted. This chapter relies on the platform randomizing rankings, but real-world ranking algorithms are usually not random; they rank by relevance, conversion rate, or ad bids, with good products systematically near the top. Then position is entangled with product quality, and using position directly to identify the exposure effect credits quality to position. The operational response is to look for exogenous variation in the ranking, such as a ranking experiment the platform actually ran, or a quasi-random perturbation from some revision of the algorithm; absent such variation, the position effect can only be cautiously bounded rather than precisely estimated.

Heterogeneity in search costs gets averaged away. This chapter estimates a single representative search cost, but real consumers' search costs vary enormously, with some comparing prices at almost zero cost and others too lazy to look past the first click. If heterogeneity is large but the model estimates only one $c$, the estimated representative value may represent neither the active comparison-shoppers nor the buy-on-first-sight consumers. Hong-Shum's framework already allows estimating the entire distribution of search costs, and the structural model can also give $c$ a random coefficient, at the cost of higher data requirements and computation.

Forward-looking and learning behavior blur the boundary between search and experience. Search in this chapter is static information acquisition, but consumers of experience goods learn across periods, and this period's choice is partly to reduce next period's uncertainty (Section 2.5). Treating a forward-looking, learning consumer as a static searcher miscounts the learning motive into search cost or preferences. Whether a dynamic model is needed depends on the product's experience attributes and the consumer's decision horizon, and later parts of this series on dynamic structure will provide the corresponding tools.

Advertising effects are the easiest to be exaggerated by selection. As Section 3.6 discussed, placement is endogenous, and only products that sell well get advertised. The advertising effect from a correlational analysis is almost inevitably too high. The lesson from Blake, Nosko and Tadelis (2015) is especially striking: their field experiment found that heavily purchased branded keyword advertising has an incremental effect of nearly zero, because those users would have come anyway. The operational discipline is that a serious estimate of an advertising effect must rely on a randomized experiment or the platform's own randomization to net out the "would have happened anyway" part, and an advertising ROI derived from observed correlation alone should be heavily discounted.

Whether price is visible before search or learned only after search changes the whole model. This chapter sets prices to be revealed only upon inspection, suited to categories dominated by match value. But on many platforms the price shows on the listing page, and the consumer sees it before searching, in which case price affects the search order itself, and the model structure changes accordingly. Before building a model, first work out which attributes are pre-search visible and which are learned only post-search, a step that cannot be skipped; get it backwards and the entire logic of identification is built wrong.

Stringing these failure modes together, the credibility of demand estimation in a search market ultimately hangs on two things: whether the search protocol and the information structure are correctly specified, and whether there is credible exogenous variation (random ranking or a randomized experiment) to separate search, preferences, and advertising. The test of sequential versus fixed-sample, ranking experiments, the modeling of search-cost heterogeneity, and wariness toward advertising correlations are all evidence provided around these two questions, and none of them substitutes for a substantive judgment about the mechanism of "how consumers actually search and what they can see."

## 8 Further Reading

::: {.readings}
Required reading, in suggested order:

- Hong and Shum (2006, RAND Journal of Economics). The method benchmark for estimating search costs from the price distribution, the full version of Sections 4.4 and 5.1 of this chapter; focus on how the fixed-sample and sequential inversions correspond to different equilibrium conditions.
- Ursu (2018, Marketing Science). The exemplar of a randomized-ranking experiment, the real-world prototype of this chapter's running case; focus on why position works through search rather than persuasion.
- Weitzman (1979, Econometrica). The original text on Pandora's boxes and the reservation-value rule; focus on how the selection, stopping, and choice rules reduce a dynamic problem to a static index.
- Diamond (1971, Journal of Economic Theory). The origin of the search paradox, a one-page argument that shows the weight of a positive search cost.
- De los Santos, Hortaçsu and Wildenbeest (2012, American Economic Review). The test of sequential versus fixed-sample search; read how it uses behavioral features to distinguish the two protocols.

Further reading:

- Stigler (1961, Journal of Political Economy). The founding text of search theory, fixed-sample search and price dispersion.
- Burdett and Judd (1983, Econometrica). How equilibrium price dispersion arises in a homogeneous market, the theoretical basis of the Hong-Shum method.
- Kim, Albuquerque and Bronnenberg (2010, Marketing Science). An implementation embedding sequential search in an aggregate demand model; read how limited search changes the estimation of demand and substitution patterns.
- Ackerberg (2003, International Economic Review). Using the contrast between first-time and experienced buyers to distinguish informative from persuasive advertising, the classic Yoplait 150 analysis.
- Erdem and Keane (1996, Marketing Science). The Bayesian learning model for experience goods, the full version of Section 2.5 of this chapter.
- Blake, Nosko and Tadelis (2015, Econometrica). The field experiment finding paid search advertising has a near-zero incremental effect, a cautionary case for the causal identification of advertising.
- Lewis and Reiley (2014, Quantitative Marketing and Economics). A randomized experiment on the effect of online advertising on offline sales, how clicks severely understate the impact of advertising.
:::

::: {.apa-refs}
- Ackerberg, D. A. (2003). Advertising, learning, and consumer choice in experience good markets: An empirical examination. *International Economic Review, 44*(3), 1007-1040. https://doi.org/10.1111/1468-2354.t01-2-00098
- Anderson, S. P., & Renault, R. (2006). Advertising content. *American Economic Review, 96*(1), 93-113. https://doi.org/10.1257/000282806776157632
- Blake, T., Nosko, C., & Tadelis, S. (2015). Consumer heterogeneity and paid search effectiveness: A large-scale field experiment. *Econometrica, 83*(1), 155-174. https://doi.org/10.3982/ECTA12423
- Burdett, K., & Judd, K. L. (1983). Equilibrium price dispersion. *Econometrica, 51*(4), 955-969. https://doi.org/10.2307/1912045
- Crawford, G. S., & Shum, M. (2005). Uncertainty and learning in pharmaceutical demand. *Econometrica, 73*(4), 1137-1173. https://doi.org/10.1111/j.1468-0262.2005.00612.x
- De los Santos, B., Hortaçsu, A., & Wildenbeest, M. R. (2012). Testing models of consumer search using data on web browsing and purchasing behavior. *American Economic Review, 102*(6), 2955-2980. https://doi.org/10.1257/aer.102.6.2955
- Diamond, P. A. (1971). A model of price adjustment. *Journal of Economic Theory, 3*(2), 156-168. https://doi.org/10.1016/0022-0531(71)90013-5
- Erdem, T., & Keane, M. P. (1996). Decision-making under uncertainty: Capturing dynamic brand choice processes in turbulent consumer goods markets. *Marketing Science, 15*(1), 1-20. https://doi.org/10.1287/mksc.15.1.1
- Fleder, D., & Hosanagar, K. (2009). Blockbuster culture's next rise or fall: The impact of recommender systems on sales diversity. *Management Science, 55*(5), 697-712. https://doi.org/10.1287/mnsc.1080.0974
- Ghose, A., & Yang, S. (2009). An empirical analysis of search engine advertising: Sponsored search in electronic markets. *Management Science, 55*(10), 1605-1622. https://doi.org/10.1287/mnsc.1090.1054
- Goldfarb, A., & Tucker, C. (2011). Online display advertising: Targeting and obtrusiveness. *Marketing Science, 30*(3), 389-404. https://doi.org/10.1287/mksc.1100.0583
- Hong, H., & Shum, M. (2006). Using price distributions to estimate search costs. *The RAND Journal of Economics, 37*(2), 257-275. https://doi.org/10.1111/j.1756-2171.2006.tb00015.x
- Honka, E. (2014). Quantifying search and switching costs in the U.S. auto insurance industry. *The RAND Journal of Economics, 45*(4), 847-884. https://doi.org/10.1111/1756-2171.12073
- Kim, J. B., Albuquerque, P., & Bronnenberg, B. J. (2010). Online demand under limited consumer search. *Marketing Science, 29*(6), 1001-1023. https://doi.org/10.1287/mksc.1100.0574
- Lewis, R. A., & Reiley, D. H. (2014). Online ads and offline sales: Measuring the effect of retail advertising via a controlled experiment on Yahoo! *Quantitative Marketing and Economics, 12*(3), 235-266. https://doi.org/10.1007/s11129-014-9146-6
- Li, X., & Hitt, L. M. (2008). Self-selection and information role of online product reviews. *Information Systems Research, 19*(4), 456-474. https://doi.org/10.1287/isre.1070.0154
- McCall, J. J. (1970). Economics of information and job search. *The Quarterly Journal of Economics, 84*(1), 113-126. https://doi.org/10.2307/1879403
- Stigler, G. J. (1961). The economics of information. *Journal of Political Economy, 69*(3), 213-225. https://doi.org/10.1086/258464
- Ursu, R. M. (2018). The power of rankings: Quantifying the effect of rankings on online consumer search and purchase decisions. *Marketing Science, 37*(4), 530-552. https://doi.org/10.1287/mksc.2017.1072
- Weitzman, M. L. (1979). Optimal search for the best alternative. *Econometrica, 47*(3), 641-654. https://doi.org/10.2307/1910412
:::
