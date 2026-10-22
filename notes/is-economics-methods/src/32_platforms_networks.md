---
title: "Platforms, Network Effects & Switching Costs"
subtitle: "Two-Sided Pricing, Lock-In, and the Measurement of State Dependence"
seriesline: "Foundations of Information Systems Economics · Chapter 32"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 32 · Platforms, Network Effects & Switching Costs"
---

## Introduction

Why does a credit card pay cash back to its cardholders and yet charge fees to merchants? Why is an app store willing to subsidize developer tools while taking a cut from successful apps? In an ordinary market, pushing one side's price down to zero or below looks like leaving money on the table. In a platform market, it may be exactly the condition for making money on the other side. Participants attract one another across the two sides, so the platform cannot price each side in isolation, nor judge the strength of competition by looking at the profit from a single side.

The self-reinforcing nature of scale brings a thornier problem. Network effects can push a market past a critical size and tip it toward a handful of platforms, while multihoming lets users stand on several sides at once, turning one side into a competitive bottleneck. Compatibility, the installed base, preannouncements of new products, and even vaporware all reshape today's adoption, because what a consumer buys is not only the current functionality but also a bet on who else will still be in the ecosystem tomorrow. At the same time, the value of a platform comes from more than scale: P2P markets need liquidity and trust, and public collaborations like open-source communities and Wikipedia must additionally solve contribution incentives and governance.

Empirically the thing most easily misjudged is "staying." A merchant that has used the same payment platform for years may be bearing real switching costs, or may simply have preferred that platform all along. The first is state dependence, the second is persistent heterogeneity, and in a short panel the sequence of choices often makes the two nearly impossible to tell apart. If a naive logit counts all of preference-driven loyalty as lock-in, the switching cost comes out systematically too large. Network effects raise a similar reflection problem: a user influences peers, peers simultaneously influence the user, and correlation by itself cannot tell us which way the arrow points.

This chapter builds a map of platform-economics mechanisms from two-sided pricing, network effects, multihoming, and compatibility, and then uses a payment-platform case to focus on the identification of switching costs, initial conditions, and unobserved heterogeneity. You will learn how to interpret one side's price by putting it back into the platform's overall profit, how to separate genuine lock-in from choice persistence, and how to judge what evidence tipping, competitive bottlenecks, and community governance each require. The most dangerous shortcut in platform research is to see scale and call it a network effect, to see retention and call it a switching cost. The task of this chapter is to turn these attractive labels back into economic propositions that can be tested.

## 1 A Number That Does Not Sit Right

::: {.case}
Cadence is a fictional market for payment-processing platforms, in which two platforms A and B compete for the card-acquiring business of small merchants. We use simulated data to tell this story. The advantage is that the data-generating process is fully known, so the behavior of every estimator can be reconciled against the truth, a teaching condition that real data can never provide. The setup is as follows. The market contains 5000 small merchants observed over 10 months. Each month each merchant routes its card payments to either A or B. A merchant's choice moves with the fees the platforms charge each month, and is also pulled by two persistent forces: first, each merchant has a fixed, to us unobserved preference between A and B (which platform fits its business model better), and second, a real switching cost, in that having used one platform last month, staying with it this month saves a migration expense. The platform's question is direct: how deeply are merchants locked in, that is, how large is the switching cost?

Because this is a simulation, we know the answer. The true switching cost corresponds to a state-dependence parameter $\delta = 0.90$ (in logit utility units), which converts to roughly $1.64 per month in migration cost. The true fee sensitivity is $\beta = -0.55$, and the standard deviation of persistent preferences across merchants is $\sigma_\theta = 1.40$.
:::

Start with the most naive approach. The choices clearly show inertia: comparing each merchant's choice this month with last month's, across the 10 months 79.2% of observations stay with the same platform and only 20.8% switch. This high degree of persistence looks like exactly the evidence for lock-in. So a researcher runs a dynamic logit: regress whether A is chosen this month on this month's fee difference and last month's choice, and read the coefficient on last month's choice as state dependence, that is, the switching cost.

$$\Pr(\text{choose A}_{it}) = \Lambda\big(\theta_0 + \beta\, \Delta\mathrm{fee}_t + \delta\, L_{i,t-1}\big)$$

Here $\Delta\mathrm{fee}_t$ is the monthly fee difference between A and B, $L_{i,t-1} \in \{-1, +1\}$ marks which platform was chosen last month, and $\Lambda$ is the logistic function. The result is $\hat\delta_{\text{naive}} = 1.357$ (SE 0.012), with a fee coefficient $\hat\beta = -0.424$. Converting state dependence into money, $\hat\delta / |\hat\beta| = 3.20$, the naive estimate says merchants face a switching cost of about $3.20 per month. Written into a report to the platform, the conclusion would read "merchants are firmly locked in, switching is expensive, and the incumbency advantage is secure."

The problem is that we know the truth is $\delta = 0.90$, about $1.64 in money terms. The naive estimate overstates the switching cost by fully 51%, nearly doubling it in dollar terms. The regression is not written wrong, the logit is used correctly, and the fee is included, so what exactly is $\hat\delta_{\text{naive}}$ estimating? This is the main thread of the chapter. The answer comes in Section 3: when there is persistent, unobserved preference heterogeneity across merchants, a dynamic logit that does not control for that heterogeneity mistakes preference-generated persistence for lock-in, so the coefficient on last month's choice is systematically too high. In Section 6 we will take this 51% overstatement apart precisely, see where it comes from, and see how an appropriate estimator pushes it back down.

## 2 The Economic Model and the Estimand

Before estimating anything, we need to be clear about two things: why platform markets have dynamics like switching costs and network effects, and which of these quantities we actually want to estimate. The former determines the economic meaning of the object of estimation, the latter determines the target of identification. This section first lays out the theoretical skeleton of platform pricing and network effects, and then writes down precisely the object of estimation for switching costs.

### 2.1 Two-Sided Pricing: Why One Side Is Subsidized

The most counterintuitive thing about a platform is its pricing. In a one-sided market, a monopolist marks up over marginal cost according to demand elasticity, and price always lies above cost. A platform is not like this. It serves two interdependent groups of users, where the size of one side is the source of value for the other, so the optimal pricing can set one side's price below cost, or even at a loss, and recoup profit from the other side. To understand this we need a model that writes the cross-group externality into utility.

The membership model of Armstrong (2006) is the clearest starting point. Let the platform connect two groups of users $i \in \{1, 2\}$, where members of one group care about the size of the other group. The utility of a member of group $i$ is

$$u_i = \alpha_i\, n_j - p_i, \qquad j \neq i$$

where $n_j$ is the number of members of the other group $j$ who join the platform, $\alpha_i$ is the per-member external benefit that a group $i$ member obtains from each member of group $j$ (local notation for this section), and $p_i$ is the access price the platform charges group $i$. The number joining each group depends on that group's own utility, $n_i = \phi_i(\alpha_i n_j - p_i)$, and the two group sizes cause one another, forming a fulfilled-expectations fixed point. The platform incurs a connection cost $f_i$ for each member of group $i$, and monopoly profit is

$$\Pi = n_1 (p_1 - f_1) + n_2 (p_2 - f_2)$$

Optimizing over $p_1, p_2$, the profit-maximizing prices satisfy Armstrong's pricing formula:

$$p_i = f_i - \alpha_j\, n_j + \frac{n_i}{\partial n_i / \partial u_i}$$

Each of the three terms on the right has a meaning. $f_i$ is the own connection cost, the last term is the usual monopoly markup, and the term that truly distinguishes a platform from a one-sided market is the middle one, $-\alpha_j n_j$. It says that when pricing group $i$, the platform must subtract a piece from cost, and the piece it subtracts is exactly the external benefit that one more member of group $i$ brings to the other group: each additional group $i$ member benefits each of the $n_j$ members of group $j$ by $\alpha_j$, for a total of $\alpha_j n_j$. The platform internalizes this external benefit by lowering group $i$'s price to attract more group $i$ members, thereby enlarging the whole platform. This term is the direct effect: holding $n_j$ and the markup fixed, each unit increase in $\alpha_j$ lowers $p_i$ by $n_j$, so the larger the externality a group brings to the other, the lower the price that group is charged. In equilibrium, however, $n_j$ and the markup also move with $\alpha_j$, so the sign of the total effect is ambiguous. When $\alpha_j$ is large enough, $p_i$ can fall below cost, and that group becomes the subsidy side, while the other becomes the money side that pays.

::: {.case}
Let us work through this formula on Cadence. Suppose merchants (group 1) obtain an external benefit $\alpha_1 = 3.0$ from each cardholder (merchants value reach highly), cardholders (group 2) obtain an external benefit $\alpha_2 = 0.6$ from each merchant (one more store that accepts cards is a nice extra for a cardholder), and both groups have a unit connection cost of 1.0. Solving for the monopoly platform's optimal prices under linear membership demand gives $p_1^\ast = 1.384$ for merchants and $p_2^\ast = 0.684$ for cardholders, with group sizes $n_1^\ast = 0.263$ and $n_2^\ast = 0.237$. The cardholder side's price of 0.684 is below its connection cost of 1.0, so it is subsidized; the merchant side's price of 1.384 is above cost, so it is taxed. For comparison, if the platform ignored the cross-group externality and priced the two sides as isolated markets, it would obtain $p_1 = 1.100$ and $p_2 = 1.000$, with neither side priced below cost and neither side subsidized. Note, however, that the cardholder side's isolated price exactly equals its cost and its demand is zero: without the externality, this side would never get off the ground. The platform price differs from the isolated price through two channels: the $-\alpha_j n_j$ term, and the markup term $n_i / (\partial n_i / \partial u_i)$, which changes with scale. In the cardholder price, the externality discount $\alpha_1 n_1^\ast = 0.789$ exceeds the markup of 0.474, so the price falls below cost. In the merchant price, the externality discount is only $\alpha_2 n_2^\ast = 0.142$, while the markup rises from 0.100 in isolation to 0.526, so the price actually ends up above the isolated price. Which side is subsidized depends on which side's externality discount exceeds its markup: because merchants bring a small externality to cardholders and cardholders bring a large externality to merchants, the platform pushes the cardholder side below cost to build scale, and recoups the money from the merchants who value reach more. Real-world payment networks, freemium software, and the subsidized side of ride-hailing all follow this same logic.
:::

Rochet and Tirole (2003) reach the same insight from another angle, with a model that is usage-based rather than membership-based. The transaction volume of the two groups is the product of their respective willingness to participate, $V = D^1(p^1)\, D^2(p^2)$, because a transaction requires both sides to be willing, so how the total price $p = p^1 + p^2$ is split between the two sides is what the platform truly decides. Their central conclusion is that even for a monopoly platform, the problem cannot be reduced to setting price at marginal cost on each side; what matters is not the level of price but the price structure, that is, how a given total price $p$ is divided between the two sides. This structure is determined by the two sides' demand elasticities: at the optimum, the ratio of the two prices equals the ratio of the two sides' elasticities evaluated at the optimal prices, rather than the inverse ratio that an ordinary monopolist would use when pricing two separate markets. This also leads to the seesaw principle of the Rochet and Tirole (2006) survey: any factor that raises the optimal price on one side generally lowers the price on the other, so the two prices move like a seesaw, up on one side and down on the other, because the platform must continually reallocate the burden to keep both sides on board.

The distinction between these two models is worth remembering: Armstrong is about membership externalities, where you pay to join and your utility depends on how many people are on the other side; Rochet-Tirole is about usage externalities, where you pay per interaction and transaction volume is the product of the two sides' willingness. Real platforms often have both, charging an entry fee as well as a transaction commission. Rysman (2009) surveys the literature on two-sided markets, and Weyl (2010) develops a more general price theory of platforms, showing that the conclusions of the Rochet-Tirole model depend on the dimensions along which users are heterogeneous. For this chapter, the significance of two-sided pricing is that it explains why platforms care so much about attracting and retaining users: once one side's scale becomes the source of the other side's value, capturing and locking in users has a direct profit motive, and this is precisely why switching costs sit at the center of platform strategy.

### 2.2 Network Effects, Tipping, and Critical Mass

Network effects make platform markets naturally prone to concentration. When a platform's value rises with its user base, a user's willingness to join depends on their expectation of the eventual size, and the expectation is itself self-fulfilling. Katz and Shapiro (1985) characterize this with fulfilled-expectations demand. Let the willingness to pay of the marginal consumer at quantity position $q$ consist of two parts: first, a stand-alone intrinsic value $r(q)$ that decreases in $q$ (heterogeneous preferences), and second, a network value $v(y^e)$ that increases in the expected network size $y^e$. Equilibrium requires expectations to be correct, $y^e = q$, so the fulfilled-expectations demand curve is

$$P(q) = r(q) + v(q)$$

The key is that $P(q)$ can be non-monotonic. With $r(q)$ decreasing and $v(q)$ increasing, over the range where the network effect overpowers preference heterogeneity the demand curve rises with scale, forming a hump. Given a price $p$, a horizontal line can intersect this humped demand curve at multiple points, corresponding to multiple equilibria.

![The fulfilled-expectations demand curve under network effects (blue) and a given price (orange dashed line). The curve first rises then falls with scale; the price line intersects it at two points, the lower one being the unstable critical mass and the higher one the stable high-adoption equilibrium.](assets/fig/fig_32_tipping.svg)

The figure above draws out this mechanism. It writes adoption as $n$ and uses an illustrative product form for the demand curve, $P(n) = (1 - n)(0.15 + 1.35\, n)$: the marginal user's type $1 - n$ times a network benefit that rises with scale. This is not the additive form above, only a specification chosen to produce a hump. Under this specification, the price line $p = 0.22$ intersects the demand curve at two points: $n \approx 0.063$ is the unstable critical mass, and $n \approx 0.826$ is the stable high-adoption equilibrium, with a further equilibrium at adoption near zero. The meaning is this: if actual adoption exceeds the critical mass of 0.063, expectations spiral upward and the market rolls toward the high equilibrium; if it falls below, adoption collapses toward zero. This is tipping. When network effects are strong enough, only the low and the high equilibria are stable, one with essentially no adopters and one with high adoption (here 0.826); the critical mass in between is merely a watershed, so the market tips toward a single platform or standard. The condition for local stability is that the demand curve crosses the price line from above at that point, $P'(q^\ast) = r'(q^\ast) + v'(q^\ast) < 0$.

This logic of network effects has two implications for this chapter. First, it explains why platforms subsidize at a loss early on: crossing the critical mass is a matter of life and death, and a temporary low price buys the momentum to roll toward the high equilibrium, a motive that stacks on top of the two-sided subsidy motive of Section 2.1. Second, it reminds us that observed high concentration and a user's long stay on some platform need not come from switching costs; they can also come from the equilibrium selection of network effects. Conflating the two makes one read lock-in where there is none. This is precisely why defining the object of estimation clearly matters so much.

### 2.3 Multihoming, the Competitive Bottleneck, and Platform Competition

Competition among two-sided platforms cannot be assessed by asking only how many users each side has; one must also ask whether users multihome. If consumers join only one platform while merchants join several at once, the platform forms a competitive bottleneck on the consumer side: to reach the consumers that each platform holds exclusively, merchants must be present on multiple platforms. In this case the consumer side is fiercely contested and easily subsidized, while the merchant side may bear high fees. Food delivery is a ready example: most users keep a single delivery app on their phone (single-home), while restaurants list on nearly every platform (multihome), so each platform holds the sole gateway to "its own batch of exclusive users," showering subsidies on users to win installs on one side while charging high commissions to restaurants on the other. This is the typical face of a competitive bottleneck. Conversely, if both sides multihome widely, the value of a platform's exclusive access falls, and pricing and market power change accordingly.

This yields three empirical predictions. First, the scarcer the single-homing side and the more it attracts the other side, the stronger the subsidy. Second, interoperability or data-portability policies that lower multihoming cost should weaken bottleneck rents, but the effect depends on which side they change. Third, observing a price of zero on one side does not imply that side has no market power; a zero price may be exactly the price structure the platform chooses to win exclusive users.

A platform can also use a membership fee and a transaction fee at the same time. A two-part tariff in a two-sided market not only allocates surplus but also changes transaction intensity and the size of joining. When comparing platforms one must look at the price structure and not just compare total prices; the core criterion of Rochet-Tirole is exactly whether reallocating the two sides' prices while holding total charges fixed changes transaction volume.

### 2.4 Compatibility, Preannouncement, and Adoption Inertia

The installed base means that even a superior technology may not immediately replace an old standard. When users switch, they compare not only how good the product is on its own (stand-alone quality) but also weigh the network size after migrating, whether complements will keep up, and whether others will move together. Farrell and Saloner (1986) identify two inefficiencies that point in opposite directions. In the first, the market stays locked onto an inferior old standard and refuses to move; this is excess inertia, switching too late, and the coordination failure in which no one dares to move first and everyone is stuck in place is one form of it. In the second, a stampede onto a new standard prematurely abandons a still-valuable old one; this is excess momentum, switching too early. Standards wars over videotape formats or charging connectors are often used to illustrate this coordination dilemma: beyond how good a standard is to use, changing standards requires everyone to change at once, so either the matter drags on unresolved or a sudden shift in expectations ignites a collective switch. Katz and Shapiro (1994) survey these questions of adoption, competition between incompatible systems, and the choice of compatibility from the perspective of systems competition.

Product preannouncement changes this coordination game. A credible preannouncement can make users wait for a new platform and coordinate a future migration; incredible vaporware can also freeze the demand of competitors and delay effective adoption. So preannouncement is neither necessarily pro-competitive nor necessarily exclusionary, and its effect depends on credibility, the distance to the release date, the installed base, and compatibility. Empirically one should observe the adoption hazard, delays in purchasing the old product, and developer investment, rather than looking only at prices after the preannouncement.

### 2.5 P2P Markets, Trust, and Commons Governance

Sharing-economy platforms lower search, price-discovery, and settlement costs through digital matching, reputation, payment, and dispute resolution, enabling dispersed individuals to supply services. Einav, Farronato and Levin (2016) survey the design of such P2P markets. But a P2P market must solve liquidity, trust, and scale all at once, and these often pull against one another. The hardest part is the cold start. A newly launched short-term rental platform cannot attract guests without hosts and cannot retain hosts without guests; this is the chicken-and-egg problem of liquidity, exactly the problem that Caillaud and Jullien (2003) address in their study of competition among intermediation platforms. To make strangers willing to transact, the platform must also set up identity verification, reviews, and dispute handling; this is trust. Yet once the verification bar is raised, early supply is blocked at the door and scale cannot grow. Supply without demand creates no liquidity, scale without verifiable quality invites adverse selection, and strict vetting raises trust but may block entry. These three goals bear down at once, which is exactly what makes the early stage of a P2P platform so treacherous.

Another kind of digital platform is not a market but commons-based peer production, for example Wikipedia and open-source communities. Its core is not two-sided pricing but the provision of a public good, contributor motivation, content bias, and governance. Symbolic recognition may raise retention: in the large-scale natural field experiment of Gallus (2017) on Wikipedia, randomly assigned symbolic awards raised the retention of new contributors. Prefilled content or AI seeding, by contrast, may crowd out volunteer contribution. Thus "more auto-generated content" can raise quantity in the short run but weaken the human capital and community production base in the long run; platform governance should jointly evaluate contribution quantity, quality, diversity, and contributor retention.

### 2.6 The Estimand for Switching Costs

Now write down the object of estimation precisely. Individual $i$ (a merchant in Cadence) chooses among several platforms $j$ in period $t$, and the utility of choosing $j$ takes a random-utility form, with the previous period's choice entered explicitly:

$$u_{ijt} = x_{jt}'\beta + \alpha_{ij} + \delta\, \mathbf{1}\{y_{i,t-1} = j\} + \varepsilon_{ijt}$$

where $x_{jt}$ are the observable characteristics of platform $j$ in period $t$ (the fee in Cadence), $\alpha_{ij}$ is the individual's persistent, time-invariant preference fixed effect for platform $j$ (unobserved by the researcher; it follows the unified notation, which writes the unit effect as $\alpha_i$, and differs both from the $\alpha_i$ of Section 2.1, the cross-group external benefit, and from the $\alpha_i$ of Chapter 27, price sensitivity; in this chapter the price sensitivity to fees is folded into $\beta$), $\mathbf{1}\{y_{i,t-1} = j\}$ indicates whether $j$ was chosen last period, and $\varepsilon_{ijt}$ is an i.i.d. logit shock. The coefficient $\delta$ is the estimand of this chapter: structural state dependence, that is, the causal effect of the previous period's choice on this period's utility. $\delta > 0$ means that repeating last period's choice brings extra utility, equivalently that a switching cost is paid whenever $y_{it} \neq y_{i,t-1}$. It captures the real forces (learning, integration, accumulated data, staff habit) that make "staying put" easier, and is the econometric measure of platform lock-in.

Converting $\delta$ into money is direct: $\delta / |\beta|$ is how many dollars of fee the switching cost is equivalent to, because $\beta$ is the marginal response of utility to a dollar of fee. In Cadence the truth is $\delta = 0.90$, $\beta = -0.55$, so the switching cost is about $1.64 per month.

We must stress the division of labor between $\delta$ and $\alpha_{ij}$, which is the pivot of the whole chapter's identification. $\alpha_{ij}$ is persistent preference, making certain merchants always lean toward A regardless of history; $\delta$ is path dependence, making any merchant more willing to stay wherever it was last period, regardless of its intrinsic preference. Both produce choice persistence in the data, but they are two entirely different economic mechanisms with different meaning and different policy implications. In a market dominated by $\alpha_{ij}$, users stay because they have already found the platform that suits them, and there is no welfare loss; in a market dominated by $\delta$, users stay because migrating is costly, and the incumbent enjoys market power as a result, which can bring welfare loss. Estimating $\delta$ correctly means separating these two mechanisms.

This section can be summarized as follows: two-sided pricing and network effects explain why platforms have strong incentives to attract and lock in users, and remind us that persistence can come from several mechanisms; the estimand of this chapter is the state-dependence parameter $\delta$ in the utility model, the switching cost, which is meaningful only when separated from persistent preference heterogeneity $\alpha_{ij}$.

## 3 Identification: Separating Lock-In from Preference

Identification logically precedes estimation. The identification problem here is: under what assumptions can the observed sequence of choices separate structural state dependence $\delta$ from persistent heterogeneity $\alpha_{ij}$? The answer to this question depends only on the data structure and the assumptions, and is independent of whether MLE or Bayes is used afterward. This is the most detailed section of the chapter, broken into numbered subsections that dissect this classic difficulty layer by layer.

### 3.1 Two Sources of Persistence

Let us first make the core of the problem plain. Choice persistence, that is, "this period's choice is highly correlated with last period's," is an observable fact in the data. In Cadence, 79.2% of observations stay unchanged, and this number itself is not in dispute. What is in dispute is where it comes from. Two entirely different mechanisms can produce the same persistence.

The first is genuine state dependence, $\delta > 0$. Last period's choice causally changed this period's utility, and the user stays put because migrating is costly. This is the lock-in we want to measure.

The second is spurious state dependence, arising from persistent unobserved heterogeneity $\alpha_{ij}$. Imagine a set of merchants that innately prefer platform A (say A's interface fits their business better) and another set that innately prefer B. A merchant that prefers A chose A this period, chose A last period, and will choose A next period, not because of any lock-in but simply because it has always preferred A. Its choice sequence is highly persistent, but there is not a shred of switching cost in that persistence. This is spurious state dependence, and Heckman's (1981) characterization of it is the starting point for this whole literature.

The difficulty is that these two mechanisms look almost identical within a single choice sequence. For a merchant that is always at A, you cannot conclude it is locked in merely because it is always at A, since a merchant that simply prefers A looks exactly the same. To separate them requires more information, and also assumptions.

::: {.intuition}
Here is an analogy. You notice that a person goes to the same restaurant every year. Does this count as being "locked into" that restaurant? Two possibilities: one is that they find changing restaurants a hassle (a real switching cost), the other is that they simply love this restaurant's flavor the most (a persistent preference). Just seeing "they come every year" cannot tell the two apart. The way to tell them apart is to see what the person does after a temporary push has passed. Suppose the restaurant raises its prices sharply for a month and they go elsewhere. If they come straight back once prices return to normal, that is taste (preference) at work; if they stay at the other restaurant after prices return, that "not coming back" is the evidence of a switching cost. Identifying switching costs relies precisely on observing whether users return once a temporary exogenous push is gone, not on observing where they sit at a given moment.
:::

### 3.2 Why the Naive Estimate Is Systematically Too High

Now explain where the 51% overstatement of Section 1 comes from. The naive dynamic logit puts last period's choice in as an ordinary regressor while leaving out persistent heterogeneity $\alpha_{ij}$. The problem is that the omitted $\alpha_{ij}$ is positively correlated with the included lagged choice.

The reasoning is this: a merchant with high $\alpha_{ij}$ (strongly preferring $j$) is both more likely to have chosen $j$ last period and more likely to choose $j$ this period. So the regressor "last period's choice" is positively correlated with the favorable $\alpha_{ij}$ left in the error term. MLE cannot tell the two apart, so it attributes this positive correlation entirely to $\delta$. Formally, when $\alpha_{ij}$ is ignored, $\mathrm{plim}\,\hat\delta > \delta$, and the estimate is systematically too high. This is Heckman's (1981) classic conclusion: insufficient control of heterogeneity exaggerates state dependence.

::: {.warning}
The direction of this bias is definite: upward. As long as persistent heterogeneity is not adequately controlled, state dependence is overstated and almost never understated. The reason is that persistent preference and last period's choice are always positively correlated (a person who prefers $j$ was more likely to choose $j$ last period), and this positive correlation only pushes $\delta$ upward. So when you see a dynamic choice model that does not control for heterogeneity report a large switching cost, the correct reaction is not to believe it, but to ask how much preference it has counted as lock-in.
:::

The Cadence simulation displays this mechanism cleanly: the truth is $\delta = 0.90$, the naive logit estimates 1.357, 51% too high. The standard deviation of $\alpha_{ij}$ in the data-generating process is 1.40, a sizable amount, and it is precisely this ignored heterogeneity that is misread as lock-in.

### 3.3 The Initial Conditions Problem

Even if we decide to model persistent heterogeneity $\alpha_{ij}$ explicitly, there is a subtler hurdle, called the initial conditions problem, also from Heckman (1981). Section 2.1 of Chapter 11 distinguished state dependence from heterogeneity in linear dynamic panels, and Sections 3.7 and 3.8 of Chapter 12 covered the RE and CRE probit and the FE logit for binary panel models.

The problem is this: in a dynamic model, the first-period choice we observe, $y_{i0}$, is not exogenous. It is itself generated by the same $\alpha_{ij}$ process before the sample begins. A merchant that prefers $j$ already tended to choose $j$ before we began observing, so its initial choice $y_{i0}$ is correlated with its $\alpha_{ij}$. If we treat $y_{i0}$ as exogenously given during estimation, or simply condition it out, we are pretending that the initial state is independent of the individual effect, and this assumption is wrong. The consequence is that $\hat\delta$ remains inconsistent.

This is clearer in the language of likelihood. The joint likelihood of a sequence can be written as

$$\ell(y_{i1}, \dots, y_{iT} \mid y_{i0}, \alpha_i) = \prod_{t=1}^{T} \Pr\big(y_{it} \mid y_{i,t-1}, x_{it}, \alpha_i\big)$$

Each term conditions on $\alpha_i$. The difficulty is that to write the likelihood of the whole sequence with respect to the observed data, one must integrate over $\alpha_i$, and the distribution used in the integration should be the conditional distribution of $\alpha_i$ given the initial choice $y_{i0}$, not its unconditional distribution, because $y_{i0}$ carries information about $\alpha_i$. Use the wrong distribution and the estimate is biased. In the Cadence data-generating process, the first month's choice is determined precisely by the merchant's $\theta_i$ (that is, $\alpha_{ij}$) plus noise, and this correlation between the initial condition and the individual effect is deliberately built in, so that we can see the cost of ignoring it.

### 3.4 The Conditional Independence Assumption

To identify $\delta$ in a dynamic discrete-choice model, one usually cannot get around a key assumption, called conditional independence. It requires that, given last period's choice and observable characteristics, this period's utility shocks $\varepsilon_{ijt}$ are no longer correlated over time, so the remaining persistence in the sequence can only come from the explicitly modeled $\alpha_{ij}$ and state dependence $\delta$.

This assumption carries a lot of weight, and it is worth stating head-on what it constrains.

::: {.assumption}
**(conditional independence)** Given $(\alpha_i, y_{i,t-1}, x_{i1}, \dots, x_{iT})$, this period's idiosyncratic shock $\varepsilon_{ijt}$ is independent of the past, and the covariates are strictly exogenous given $\alpha_i$ (fees that respond to demand shocks would violate this). Equivalently, once persistent heterogeneity and last period's choice are controlled for, the choice sequence has no other source of persistence.
:::

Why it matters, and why it is uncomfortable: it routes all persistence in the sequence into two exits, either through $\alpha_{ij}$ (persistent preference) or through $\delta$ (state dependence). If reality has a third exit, say the idiosyncratic shock is itself serially correlated (a merchant does good business for several months in a row, or one of a platform's services is unstable for several months in a row), then this extra persistence has nowhere to go and gets squeezed into $\delta$, causing fresh overstatement. This is exactly why many people are reserved about the assumption: it is not a testable fact but a modeling decision that forces complex reality into two categories. Acknowledging this is more honest than pretending it must hold.

### 3.5 Isomorphism with the Reflection Problem

Pulling the lens back, the identification difficulty for switching costs is actually a special case of a broader class of problems, which in social economics is called the reflection problem, characterized precisely by Manski (1993). Seeing this isomorphism is especially useful for studying network and peer effects on platforms, because there one hits exactly the same identification wall.

Manski distinguishes three mechanisms that make individuals within a group behave similarly. The endogenous effect: an individual's outcome depends on the outcomes of others in the group (a genuine behavioral spillover, a social multiplier). The contextual effect: an individual's outcome depends on the exogenous characteristics of others in the group. The correlated effect: people in the group behave alike merely because they share an environment or shock, or self-selected into the same group, and there is no spillover at all. The classic linear-in-means model writes in the first two:

$$y = \alpha + \beta\, \mathbb{E}[y \mid g] + \mathbb{E}[x \mid g]'\gamma + x'\delta + u$$

Here the $\alpha, \beta, \gamma, \delta$ in Manski's model are local notation for this passage, unrelated to the $\delta$ denoting state dependence and the $\beta$ denoting price sensitivity earlier in the chapter.

Here $\beta$ is the endogenous spillover and $\gamma$ is the contextual effect. Manski's conclusion is that this model is not identified. Taking the group expectation solves for

$$\mathbb{E}[y \mid g] = \frac{\alpha}{1 - \beta} + \bar{x}_g'\,\frac{\gamma + \delta}{1 - \beta}$$

The group mean $\mathbb{E}[y \mid g]$ becomes an exact linear function of the group characteristic $\bar{x}_g$. Substituting it back into the original equation gives $y = \alpha/(1 - \beta) + \bar{x}_g'(\gamma + \beta\delta)/(1 - \beta) + x'\delta + u$. As long as $x$ varies within groups, the coefficient $\delta$ on one's own characteristics is still identified, but the structural parameters $(\alpha, \beta, \gamma)$ contain one more unknown than there are observable reduced-form coefficients, so they cannot be recovered from those coefficients. Intuitively, the "image in the mirror" (the peer mean) and the regressors that generate it move in lockstep, and you cannot tell whether group members behave alike because behavior is contagious among peers ($\beta$), because peer characteristics themselves matter ($\gamma$), or because everyone shares a common shock.

This problem and the switching-cost difficulty are two faces of the same thing. There we must separate "the causal effect of last period's choice" ($\delta$, state dependence) from "the correlation created by persistent preference" ($\alpha_{ij}$, a kind of correlated effect); here we must separate the endogenous spillover from the correlated effect. In both cases the obstacle is that an endogenous, self-reinforcing quantity (peer behavior, or one's own history) and a static source of confounding (peer characteristics plus common shocks, or persistent preference) move in lockstep in the observed data. The way out is also shared: break this lockstep. In the network setting, Bramoullé, Djebbari and Fortin (2009) point out that as long as reference groups do not overlap completely (my neighbors differ from yours, so the network has intransitivity), $Gy$ and $Gx$ are no longer collinear. Here $G$ is the row-normalized network adjacency matrix, and $Gy$ and $Gx$ are the vectors of each person's peer averages of outcomes and of characteristics. Then $\beta$ can be identified, and the characteristics of "friends of friends" can serve as instruments. In the switching-cost setting, the way out is to find temporary exogenous price variation that pushes users off their original platform, and then to see whether they stay on the new platform once the push is gone. This is the formal version of the analogy in Section 3.1, and the next subsection develops it.

### 3.6 What Can Pin Down $\delta$

Since observing where users sit at a given moment cannot separate lock-in from preference, identifying $\delta$ must rely on variation, specifically exogenous price variation that can push users to switch. Imagine a merchant who moves to B in a month when A's fee is temporarily high, after which A's fee returns to its usual level. If $\delta = 0$, the merchant's choice this month depends only on $\theta_i$ and this month's fee; having been on B last month leaves no trace, and the merchant returns to A with the same probability as before the fee rise. This round trip of "being pushed away by price and switching back" is fully explained by persistent preference plus price sensitivity. If $\delta > 0$, the merchant gains a retention bonus of $\delta$ on B: in a month with a zero fee gap, it returns to A with probability $\Lambda(\theta_i - \delta)$, below the $\Lambda(\theta_i)$ it would have without state dependence, so some merchants stay on B. What exposes state dependence is therefore the persistence after a temporary push, that is, a causal effect of past prices on today's choice, not the round trip itself. The degree of that persistence calibrates $\delta$, and $\beta$ converts it into $\delta / |\beta|$ dollars. This is why the fee difference in the Cadence data-generating process varies month to month; without this variation, $\delta$ cannot be identified.

Dubé, Hitsch and Rossi (2010) also use a permutation test to check whether their heterogeneity specification is adequate. Suppose each household's purchase sequence is randomly reshuffled. The persistence created by persistent heterogeneity does not depend on time order and survives the shuffle intact, whereas genuine state dependence does depend on time order and is destroyed. So when the same model is re-estimated on the reshuffled data, a model that has adequately controlled for heterogeneity should return a state-dependence estimate near zero; a clearly positive estimate means that heterogeneity is leaking into $\delta$. The test still re-estimates the same parametric model, and it presupposes the conditional independence of Section 3.4, because serially correlated shocks also depend on time order.

The key points of this section can be summarized as follows: the core obstacle to identifying switching costs is that genuine state dependence and persistent heterogeneity produce the same persistence in the choice sequence; the naive estimate systematically overstates $\delta$ by ignoring heterogeneity, and the initial conditions problem further aggravates the inconsistency; this difficulty is isomorphic to Manski's reflection problem, and the way out in both is to break the lockstep between the endogenous quantity and the confounding source using exogenous variation or network structure; and what can pin down $\delta$ is whether users remain on the new platform after a temporary exogenous price push, not whether they are seen staying put.

## 4 Estimation: From the Failed Naive Logit to Separating Heterogeneity

This section advances by "where the previous method fails gives rise to the next method." Every step answers the same question: how to peel persistent heterogeneity out of state dependence.

### 4.1 The Naive Pooled Logit

The easiest thing to do is to stack the whole panel and run one logit, treating last period's choice as an ordinary regressor and making no distinction between individuals:

$$\Pr(y_{it} = 1) = \Lambda\big(\theta_0 + \beta\, \Delta\mathrm{fee}_t + \delta\, L_{i,t-1}\big)$$

Its failure has already been laid out in Section 3: ignoring $\alpha_{ij}$, last period's choice absorbs persistent preference, and $\hat\delta$ is biased upward. In Cadence it delivers 1.357, 51% above the truth. This is not a bias that a larger sample can rescue; however large the sample, $\mathrm{plim}\,\hat\delta$ remains above the truth, because it is a specification error rather than sampling noise. Its only use is as a reference upper bound, telling you "how wrong you would be if you ignored heterogeneity entirely."

### 4.2 Random-Effects Logit: A Partial Fix

The natural first patch is to add a random intercept for each individual to absorb persistent preference:

$$\Pr(y_{it} = 1 \mid \theta_i) = \Lambda\big(\theta_i + \beta\, \Delta\mathrm{fee}_t + \delta\, L_{i,t-1}\big), \qquad \theta_i \sim \mathcal{N}(0, \sigma_\theta^2)$$

Estimation integrates over $\theta_i$ according to its normal distribution (numerically via a Laplace approximation, nAGQ=1), which is the standard random-effects logit. It gathers most of the persistent heterogeneity into $\theta_i$, so $\hat\delta$ falls substantially. But it does not handle the initial conditions problem of Section 3.3: the model assumes $\theta_i$ is independent of the initial choice $y_{i0}$, while in reality the two are correlated (a merchant that prefers A is initially more likely to be at A). This ignored correlation leaves $\hat\delta$ still slightly too high. In Cadence it estimates 0.978, about 9% above the truth of 0.90, already a huge improvement over the naive 1.357 but not yet on target, and $\hat\sigma_\theta = 1.257$ also understates the true 1.40, because part of the heterogeneity is still misassigned elsewhere.

### 4.3 The Wooldridge CRE Logit: Handling Initial Conditions

Wooldridge (2005b) gives a surprisingly simple solution to the initial conditions problem, and it is the workhorse of this chapter. Since $\theta_i$ is correlated with the initial choice, do not assume they are independent; instead model the distribution of $\theta_i$ explicitly as a conditional distribution depending on the initial choice and the history of covariates:

$$\theta_i \mid (y_{i0}, \bar{x}_i) \sim \mathcal{N}\big(\lambda_0 + \lambda_1 y_{i0} + \bar{x}_i'\lambda_2,\ \sigma_u^2\big)$$

where $y_{i0}$ is the initial choice and $\bar{x}_i = T^{-1}\sum_t x_{it}$ is the individual's time average of the covariates. This is the classic device of Chamberlain (1980) and Mundlak (1978), that is, the CRE specification of Section 3.7 of Chapter 12 (written there as $\alpha_i = \xi_0 + \bar X_i'\xi + a_i$), and the end of Section 3.8 of Chapter 12 also mentions this approach of Wooldridge (2005b). Here $\lambda_0$, $\lambda_1$ and $\lambda_2$ are local notation for this chapter (Chapter 12 writes $\xi$); they are coefficients, not time effects. The implementation is extremely light. Write $u_i = \theta_i - \lambda_0 - \lambda_1 y_{i0} - \bar{x}_i'\lambda_2$ for the deviation of $\theta_i$ from its conditional mean; then simply add the initial choice $y_{i0}$ (and the covariate means, if they vary across individuals; in Cadence the fee gap is the same for all merchants, so this term is absorbed into the intercept) as extra regressors in the random-effects logit, with $u_i$ as the random intercept.

$$\Pr(y_{it} = 1 \mid y_{i0}, u_i) = \Lambda\big(\lambda_0 + \lambda_1\, y_{i0} + \beta\, \Delta\mathrm{fee}_t + \delta\, L_{i,t-1} + u_i\big), \qquad u_i \sim \mathcal{N}(0, \sigma_u^2)$$

Entering $y_{i0}$ into the conditional distribution amounts to absorbing the information the initial choice carries about $\theta_i$, thereby sidestepping the initial conditions problem. Its cost is a parametric assumption: $\theta_i$ is normal given $(y_{i0}, \bar{x}_i)$. Under this assumption $\hat\delta$ is consistent. In Cadence it estimates 0.893 (SE 0.014), almost exactly the truth of 0.90, and in money $1.66 against the truth of $1.64, a negligible error. The standard deviation of the random intercept is now estimated at $\hat\sigma_u = 1.180$. This measures the conditional heterogeneity (conditional variance) given the initial choice, so it is naturally smaller than the unconditional RE version and should not be compared directly with the truth of 1.40; what really matters is that $\delta$ is now essentially unbiased. Taking the three steps in order, from the naive 1.357 to the RE 0.978 to the CRE 0.893, one can clearly see that each added layer of treatment for heterogeneity and initial conditions peels off one more layer of the spurious part of state dependence.

### 4.4 Honoré-Kyriazidou: No Distributional Assumption

Wooldridge's solution requires a parametric assumption on the distribution of $\theta_i$. If you are unwilling to assume even that, Honoré and Kyriazidou (2000) give a semiparametric dynamic-logit estimator that allows the fixed effect $\alpha_i$ to be entirely arbitrary. Its idea inherits Chamberlain's conditional logit: find events under which $\alpha_i$ is exactly canceled. For the binary dynamic logit, when each individual has at least four observations $y_{i1}, \dots, y_{i4}$ (the first of which serves as the initial condition, so that $T \geq 3$ in the notation of Section 3.3), one can construct a conditional likelihood on the subsample where the individual's first and last choices are fixed and the middle two periods switch, and its contribution takes the form

$$\Pr\big(y_{i2} = 1 \mid y_{i2} + y_{i3} = 1,\ y_{i1}, y_{i4}, x_{i3} = x_{i4}, x_i, \alpha_i\big) = \frac{1}{1 + \exp\big((x_{i3} - x_{i2})'\beta + 2\delta(y_{i4} - y_{i1})\big)}$$

On the event $x_{i3} = x_{i4}$, $\alpha_i$ drops out of this expression (the coefficient is $2\delta$ because this chapter codes last period's choice as $L_{i,t-1} = 2y_{i,t-1} - 1$). The cost is threefold: with continuous covariates, the event must be approximated by $x_{i3} \approx x_{i4}$ through kernel weighting $K((x_{i3} - x_{i4})/h_n)$, the convergence rate is slower than parametric methods ($\sqrt{n h_n}$ rather than $\sqrt{n}$), and it fails when the regressors contain time dummies or are highly persistent. It forms a trade-off with Wooldridge: Wooldridge places a parametric assumption on the distribution of $\alpha_i$ in exchange for fast convergence and simple implementation; Honoré-Kyriazidou makes no distributional assumption in exchange for robustness, at the cost of precision and range of applicability. Which to choose depends on which end you worry more about going wrong.

### 4.5 DHR's Mixture-of-Normals

Dubé, Hitsch and Rossi (2010) take the route of modeling heterogeneity extremely flexibly. Not content with a single normal random intercept, they let the entire parameter vector $\theta_i = [\alpha_i, \beta_i, \delta_i]$ (preference, price sensitivity, and state dependence all varying across individuals; in this section and Section 5.1, $\theta_i$ denotes this individual parameter vector, not the scalar random intercept of Section 4.2) follow a mixture of normals:

$$P(\theta_i \mid \pi, \mu, \Sigma) = \sum_{k} \pi_k\, \phi(\theta_i \mid \mu_k, \Sigma_k)$$

First draw a latent class $k$ with probability $\pi$, then draw $\theta_i$ from that class's multivariate normal, and four or five normal components can approximate a fairly complex heterogeneity distribution. Estimation uses hierarchical Bayes plus MCMC. The philosophy of this approach is: rather than worry that an underspecified heterogeneity overstates state dependence, build heterogeneity rich enough and let the data decide how much persistence remains attributable to genuine state dependence. After doing this on consumer-goods panels, they find that even after allowing for such rich heterogeneity, a substantial part of the inertia is still genuine state dependence, and they conclude that lock-in is not a pure artifact of heterogeneity.

### 4.6 The Supply Side: Harvesting or Investment

Having estimated the switching cost, it is natural to ask how it affects the platform's pricing, a step that takes us back to the two-sided pricing of Section 2.1. Switching costs turn the platform's pricing into a dynamic problem, and the switching-cost literature, surveyed by Farrell and Klemperer (2007), decomposes the dynamic optimal price into two opposing forces. Cabral (2011), who studies dynamic price competition with network effects, finds the same investment motive there: the entry price offered to win new users can fall below cost. The first-order condition of the platform's value function can be rearranged as

$$p_i - c_i = \underbrace{\frac{q_i}{-\,\partial q_i / \partial p_i}}_{\text{harvesting}} - \underbrace{\rho\, \frac{\partial \tilde{V}_i}{\partial q_i}}_{\text{investment}}$$

Here the index $i$ denotes a platform, unrelated to the group index of Section 2.1: $p_i$, $c_i$ and $q_i$ are platform $i$'s price, marginal cost and share, and $\tilde{V}_i$ is its continuation value (the discounted value of future profit), which moves with this period's share $q_i$. The first term is harvesting: switching costs make locked-in users insensitive to price, so $|\partial q_i / \partial p_i|$ falls with the switching cost, and the platform can therefore raise prices on existing users to harvest them. The second term is investment, where $\rho$ is the discount factor (local notation for this section, unrelated to the $\beta$ denoting fee price sensitivity earlier): winning one more user who will be locked in has value that is realized continually in the future, so the platform has an incentive to cut prices now to grab share, and this term vanishes when the switching cost is zero and grows as the switching cost rises. The two forces point in opposite directions, so the net effect of switching costs on prices and competition is ambiguous. The common intuition is that harvesting prevails and switching costs raise prices, but this is not a theorem; a very small switching cost may instead make the market more competitive, because the investment motive dominates during the fight for share. This ambiguity is itself an important conclusion: high switching costs do not necessarily mean high prices, and it depends on whether the market is in the phase of fighting for share or harvesting the installed base. It also explains why estimating $\delta$ correctly is so crucial: assessing a platform's market power and simulating the price effect of a platform merger both require knowing first how deeply users are locked in.

The comparison of the three main estimators is summarized in the table below.

| Dimension | Naive pooled logit | RE / Wooldridge CRE | Honoré-Kyriazidou |
|---|---|---|---|
| Heterogeneity $\alpha_i$ | Ignored entirely | Parametric distribution (normal) | Arbitrary fixed effect, no distribution |
| Initial conditions | Not handled | Wooldridge corrects with $y_{i0}$ | Canceled by conditioning |
| $\delta$ bias | Systematically upward (Cadence +51%) | Essentially unbiased after correction | Consistent, but low precision |
| Convergence rate | $\sqrt{n}$ | $\sqrt{n}$ | $\sqrt{n h_n}$, slower |
| Main cost | Overstates lock-in | Distributional assumption on $\alpha_i$ | Needs at least four observations per individual and covariate overlap |
| Implementation | glm | glmer plus initial choice | Semiparametric conditional likelihood |

The roadmap of this section is now complete: the naive logit systematically overstates lock-in by ignoring heterogeneity; random effects reclaim most of the heterogeneity but miss the initial conditions; Wooldridge patches this with the initial choice as an auxiliary regressor, a simple and reliable workhorse; when unwilling to assume a distribution use Honoré-Kyriazidou, at the cost of precision; when seeking flexibility use DHR's mixture-of-normals; and the estimated $\delta$ returns to platform pricing via the harvesting-investment decomposition.

## 5 Anchoring Papers

Methods prove themselves only in real research. This section anchors the chapter in two papers: one that establishes the methodology of switching-cost estimation and one that pushes it into a policy setting with a surprising welfare conclusion. Each is laid out along the five elements of paper, method, data, results, and limitations, with a focus on how the identifying assumptions are defended.

### 5.1 Dubé, Hitsch and Rossi (2010, RAND Journal of Economics)

::: {.case}
Paper and place in the history of the method: "State Dependence and Alternative Explanations for Consumer Inertia," RAND Journal of Economics 41(3), 417-445. This paper is the methodological benchmark for empirical switching costs, and it answers the main thread of this chapter head-on: of the observed brand inertia, how much is genuine state dependence and how much is merely an artifact of persistent heterogeneity?

Method: a random-coefficients logit with a lagged-choice dummy, modeling the heterogeneity of $\theta_i = [\alpha_i, \beta_i, \delta_i]$ extremely flexibly with a mixture of normals, estimated by hierarchical Bayes and MCMC (Section 4.5 of this chapter). Identification relies on whether consumers whom a price promotion temporarily pushed to another brand keep buying that brand after the promotion ends (Section 3.6 of this chapter). The paper also uses a permutation test to check the heterogeneity specification: once the time order of each household's purchase sequence is shuffled and the same model is re-estimated, the state-dependence estimate should be near zero if heterogeneity has been adequately controlled for.

Data: household purchase panels of consumer goods; the paper analyzes categories such as refrigerated orange juice and margarine, which are purchased frequently and heavily promoted, providing exactly the price variation identification requires.

Results: even after allowing for very rich heterogeneity, a substantial part of the inertia is still genuine structural state dependence, with $\delta > 0$ economically significant. In other words, inertia is not pure spurious heterogeneity, switching costs genuinely exist, and their magnitude relative to price is large enough to affect firms' dynamic pricing.

Limitations: the conclusion depends on the conditional independence assumption (Section 3.4 of this chapter), that idiosyncratic shocks are no longer serially correlated after controlling for rich heterogeneity, and this assumption is untestable. The price variation needed for identification must be sufficiently exogenous; if promotions themselves respond to demand, price is no longer a clean push. The heterogeneity is approximated by a finite number of normal components, and the choice of the number of components is a modeling decision. The anchoring paper of Section 5.3 of Chapter 12 offers a cautionary counterpoint: on a household scanner panel of bottled water, the standard random-coefficients logit finds clearly positive state dependence, while a quasi-experiment that exploits hurricane stock-outs detects none. A positive $\delta$ from a lagged-choice model is therefore a conclusion conditional on the heterogeneity specification.
:::

This paper's defense strategy is worth emulating. It does not simply assert that heterogeneity has been fully controlled for; it tests the claim with the permutation test, in which a clear state-dependence estimate from the same model on time-shuffled data would show that heterogeneity has not been controlled for. Turning an identification claim into a robustness test that can actually be carried out is far more forceful than piling up assumptions in the text.

### 5.2 Handel (2013, American Economic Review)

::: {.case}
Paper: "Adverse Selection and Inertia in Health Insurance Markets: When Nudging Hurts," American Economic Review 103(7), 2643-2682. The question is how large the inertia in health insurance is, and whether reducing that inertia is good or bad for welfare, with a surprising answer.

Method: write switching costs into a CARA expected-utility model of insurance choice, in which employees choose among several plans, and utility contains a switching-cost parameter $\eta$ that is paid only when "not re-enrolling in last period's plan," estimated by simulated MLE of a random-coefficients probit. The key to identification comes from an institutional detail: in one year the employer forced everyone to re-choose a plan, so there was no default and the switching cost was zero, and this period's choices cleanly reveal preferences free of inertia; in later years the default is to stay in the current plan and the switching cost reappears. Using the inertia-free year as a baseline separates inertia from preference.

Data: the insurance-choice and claims panel of employees at a large self-insured employer, where the plans differ in financial terms but share the same medical network, which makes the argument that "choice differences come from inertia rather than network differences" cleaner.

Results: the estimated switching cost is large, about $1729 for single employees and about $2480 per year for those with dependents. More important is the counterfactual: cutting the switching cost by three-quarters, with premiums fixed, appears to improve matching and add about $114 of consumer surplus; but once premiums are allowed to adjust endogenously with the enrollment mix, welfare instead falls by about 7.7% (roughly $115 per person per year), because the adverse selection that inertia had been holding down is released, healthy people rush to cheap plans and push generous plans into a death spiral, so switchers gain while stayers lose.

Limitations: the CARA and claims-distribution specifications affect the estimated risk premium. Identification leans on the exogeneity of the forced re-enrollment year; if employees already anticipated the future default rule in that year and chose strategically, the baseline is no longer clean. The policy implication of the conclusion is highly sensitive to how premiums respond endogenously, which is both where the paper's insight lies and where it is fragile.
:::

Taking the two papers together, the meaning of anchoring is clear: Dubé-Hitsch-Rossi show how to use flexible heterogeneity to separate genuine lock-in from artifacts, and a permutation test to check the heterogeneity specification; Handel shows how a clean exogenous event (forced re-enrollment) calibrates the switching cost and reveals a counterintuitive welfare lesson, namely that reducing switching costs need not improve welfare because it loosens the constraint on adverse selection. The latter especially reminds platform and policy researchers that measuring lock-in is itself a means, and what really must be evaluated is its role in the whole market equilibrium.

## 6 A Full Walkthrough on the Cadence Data

Now run the tools of Section 4 in full on the Cadence panel. The code below uses R 4.6.1 and fixes set.seed(20) for reproducibility, and every number cited in the text comes from the actual run output of this code.

### 6.1 The DGP

The design parameters are as follows: 5000 merchants, 10 months, the first month's choice as the initial condition, and months 2 through 10 (9 periods) used for dynamic estimation. The latent utility of a merchant choosing A relative to B is

$$U^\ast_{it} = \theta_i + \beta\, \Delta\mathrm{fee}_t + \delta\, L_{i,t-1} + \nu_{it}, \qquad \nu_{it} \sim \text{Logistic}(0, 1)$$

$\Delta\mathrm{fee}_t$ is the monthly fee difference between A and B (in dollars, an exogenous AR(1) time series), $L_{i,t-1} \in \{-1, +1\}$ marks which platform was chosen last month, and $\theta_i \sim \mathcal{N}(0, \sigma_\theta^2)$ is persistent preference. A is chosen when $U^\ast_{it} > 0$.

```r
set.seed(20)
N_merch  <- 5000; T_months <- 10
beta_true <- -0.55; delta_true <- 0.90; sigma_theta <- 1.40

theta <- rnorm(N_merch, 0, sigma_theta)
dfee  <- numeric(T_months); dfee[1] <- rnorm(1, 0, 1.2)
for (t in 2:T_months) dfee[t] <- 0.5 * dfee[t - 1] + rnorm(1, 0, 1.2)

choiceA <- matrix(NA_integer_, N_merch, T_months)
u0 <- theta + rlogis(N_merch) + beta_true * dfee[1]   # initial condition
choiceA[, 1] <- as.integer(u0 > 0)
for (t in 2:T_months) {
  L_prev <- 2L * choiceA[, t - 1] - 1L
  ustar  <- theta + beta_true * dfee[t] + delta_true * L_prev + rlogis(N_merch)
  choiceA[, t] <- as.integer(ustar > 0)
}
```

There are two key design decisions. First, the initial choice is generated by $\theta_i$ plus noise, correlated with persistent preference, which deliberately creates the initial conditions problem of Section 3.3. Second, the fee difference varies from month to month, which supplies what identification needs: exogenous variation that temporarily pushes merchants to switch. The truth is recorded as follows: $\delta = 0.90$, a monetized switching cost of about $1.64 per month; and the overall share of merchants choosing A over months 2 through 10 is 0.466.

### 6.2 The Observed Inertia

First look at the persistence in the data. Compare each merchant's choice this month with last month's:

```r
est <- dt[month >= 2]
stay_rate <- est[, mean(choiceA == lagA)]     # 0.792
```

79.2% of observations stay with the same platform, a switch rate of 20.8%. Viewed in isolation, this high inertia looks like evidence of deep lock-in. Section 6.5 uses counterfactual simulations to see how much of it each of the two channels can support.

### 6.3 The Naive and Random-Effects Estimates

The naive pooled logit, making no distinction between individuals:

```r
m_naive <- glm(choiceA ~ dfee + Lprev, family = binomial, data = est)
```

gives $\hat\delta_{\text{naive}} = 1.357$ (SE 0.012), $\hat\beta = -0.424$, a monetized switching cost of $3.20. This is the number that did not sit right in Section 1, 51% above the truth. Add a random intercept to absorb persistent preference, but do not handle initial conditions:

```r
m_re <- glmer(choiceA ~ dfee + Lprev + (1 | merchant),
              family = binomial, data = est, nAGQ = 1)
```

$\hat\delta_{\text{RE}} = 0.978$ (SE 0.014), $\hat\sigma_\theta = 1.257$, monetized $1.87. A substantial drop from the naive 1.357, but still about 9% above the truth of 0.90, because the correlation between the initial choice and $\theta_i$ is not handled, and this correlation keeps pushing $\hat\delta$ upward.

### 6.4 The Wooldridge CRE

Add the initial choice initA as an auxiliary regressor to handle initial conditions:

```r
m_cre <- glmer(choiceA ~ dfee + Lprev + initA + (1 | merchant),
               family = binomial, data = est, nAGQ = 1)
```

$\hat\delta_{\text{CRE}} = 0.893$ (SE 0.014), $\hat\beta = -0.537$, $\hat\sigma_u = 1.180$, a monetized switching cost of $1.66. These are almost exactly the true values of 0.90 and $1.64. The three-step progression is plain to see: naive 1.357, RE 0.978, CRE 0.893, and each layer of confounding handled peels off one more layer of spurious state dependence.

### 6.5 Decomposing the Inertia: What Each Channel Holds Up

How much of the 79.2% inertia actually comes from genuine switching costs? A simulation cannot give a unique split, but it can show how much each channel supports on its own: turn off one channel at a time, hold the rest of the structure fixed, re-simulate, and see where the inertia drops to.

```r
sim_stay <- function(delta_c, sigma_c, seed) { ... }  # re-simulate stay rate
stay_full      <- sim_stay(delta_true, sigma_theta, 20)  # 0.790  both
stay_no_sd     <- sim_stay(0,          sigma_theta, 20)  # 0.628  heterogeneity only
stay_no_hetero <- sim_stay(delta_true, 0,           20)  # 0.708  state dep only
stay_neither   <- sim_stay(0,          0,           20)  # 0.509  fees + noise only
```

![Decomposing the inertia: the stay rate in four counterfactual worlds. With only fees and noise the baseline is about 0.509; turning on heterogeneity alone raises it to 0.628, turning on state dependence alone raises it to 0.708, and both together reach 0.790. The observed high inertia is supported jointly by the two forces, and dropping either one causes a marked decline.](assets/fig/fig_32_inertia.svg)

The figure above draws the four counterfactual worlds together. With only fees and noise, the baseline stay rate is 0.509, the repetition produced by random fluctuation alone in a world with neither persistent preference nor lock-in. Turning on persistent heterogeneity alone raises it to 0.628; turning on state dependence alone raises it to 0.708; both together reach 0.790, agreeing with the observed inertia of 0.792 computed directly from the panel in Section 6.2 (the two come from two simulation paths of the same data-generating process with slightly different random-number sequences, so they differ slightly in the third decimal). The meaning is clear: the observed high inertia is not the work of any single mechanism. Heterogeneity and state dependence each contribute a sizable piece, and dropping either one makes inertia fall markedly. The two pieces do not simply add up, however, and the size of each depends on which channel is switched on first; what the decomposition shows is therefore that no single channel can support 0.790 on its own, not what exact share is lock-in. This makes the question Dubé-Hitsch-Rossi set out to answer concrete on Cadence, and it explains why a 79% inertia number, viewed statically, can be read neither directly as lock-in nor directly as preference.

### 6.6 The Overall Reconciliation

Line up the estimators and view them against the truth:

| Estimator | $\hat\delta$ | SE | Monetized (USD/month) | Relative to truth |
|---|---|---|---|---|
| True (DGP) | 0.900 | | 1.64 | 100% |
| Naive pooled logit | 1.357 | 0.012 | 3.20 | +51% |
| RE logit (no IC fix) | 0.978 | 0.014 | 1.87 | +9% |
| Wooldridge CRE logit | 0.893 | 0.014 | 1.66 | -1% |

![The three estimates of state dependence $\delta$ against the truth (dashed line, 0.90). The naive pooled logit lies far above the truth, random effects fall back significantly but remain too high, and the Wooldridge CRE is almost exactly on target. Error bars are 95% confidence intervals.](assets/fig/fig_32_estimates.svg)

The figure above draws the three estimates together with their 95% confidence intervals, with the dashed line at the truth of 0.90. The naive estimate's point and interval both lie far to the right of the truth, random effects pull it back a good deal but the interval still does not cover the truth, and the Wooldridge CRE's point falls almost on the dashed line with its interval firmly covering the truth.

The reconciliation of this section can be summarized as follows: on the same data, with the price variation needed for identification available, the naive logit overstates the switching cost by 51% by ignoring heterogeneity, random effects correct most of it but still run 9% high by missing initial conditions, and the Wooldridge CRE is essentially unbiased once it patches initial conditions; the observed 79% inertia, through counterfactual decomposition, is seen to be supported jointly by persistent heterogeneity and genuine state dependence, and separating the two is exactly the purpose of the whole chapter.

## 7 Failure Modes of Identification and Robustness

In the simulation, persistent heterogeneity and state dependence are constructed, but in real research their separation can fail at any moment. This section lays out the most common ways it fails and the actionable responses.

Underspecified heterogeneity is the number-one threat and also the most hidden. Section 3 proved that insufficient control of heterogeneity biases upward, but there is no hard standard for "enough." A single normal random intercept may not capture the true preference distribution, and residual heterogeneity keeps seeping into $\hat\delta$. There are two actionable responses. First, like Dubé-Hitsch-Rossi, model heterogeneity more flexibly (a mixture of normals) and see whether $\hat\delta$ keeps falling as heterogeneity grows richer; if it falls, the previous layer was underspecified. Second, run a permutation test: shuffle the time order of each individual's sequence and re-estimate the same model. The shuffle destroys genuine state dependence, so the estimate should be near zero; a clearly positive estimate means the model has not absorbed persistent differences across individuals, so the original $\hat\delta$ also contains heterogeneity.

The handling of initial conditions cannot be skipped. This is especially true in short panels, since the shorter the panel the greater the weight of the initial period, and the larger the bias from ignoring its correlation with $\theta_i$. Cadence has only 10 periods, and the gap between RE and CRE (0.978 versus 0.893) comes mainly from this initial-conditions term. In practice, as long as the panel is not long, one should default to Wooldridge's approach and put the initial choice in, at almost zero cost yet eliminating a systematic bias.

Once conditional independence fails, $\hat\delta$ becomes too high again. If the idiosyncratic shock is itself serially correlated (a merchant does especially good business for several months in a row, or one of a platform's services is unstable for several months in a row), this extra persistence has nowhere to go and gets counted into state dependence. The diagnostic idea is to see whether the residuals or the choice sequence still have persistence beyond first-order Markov after the model is controlled for; if so, consider adding a richer dynamic structure or explicitly modeling serially correlated shocks, rather than quietly handing it to $\delta$.

Insufficient exogenous price variation leaves identification without a fulcrum. As Section 3.6 explained, $\delta$ is calibrated with price variation that temporarily pushes users to switch, by seeing whether they stay once the push has passed. If prices hardly move, or price changes are themselves endogenous responses to demand, identification hangs in the air. The diagnostic is to check how many price-driven switches there actually are in the data, and whether the source of price variation is credibly exogenous (a platform's pricing experiment, a regulatory shock, or cost pass-through are relatively clean sources). When variation is too scarce, the honest thing is to admit that the switching cost can only be bounded within a wide interval, rather than to report a precise but fragile point estimate.

The reflection problem shadows research on network and peer effects on platforms. If the research question is not the history dependence of a single user but the spillovers among users (whether one merchant's adoption of a feature drives its neighbors to adopt), then the Manski dilemma of Section 3.5 enters directly: the endogenous spillover and the correlated effect move in lockstep at the group-mean level and cannot be separated. The response is to exploit network structure itself, since when reference groups partially overlap (my neighbors differ from yours) collinearity can be broken, using the exogenous characteristics of "friends of friends" as instruments, or relying on exogenous group assignment. Treating this problem as a simple regression usually yields spillovers that are artifacts of common shocks.

Forward-looking behavior blurs the boundary between switching costs and expectations. The state dependence in this chapter is backward-looking, with last period's choice affecting this period. But if users are forward-looking, they will choose cautiously now because they anticipate future lock-in, which entangles with switching costs. Identifying forward-looking behavior requires a dynamic structural model (the single-agent dynamic framework of Chapter 31 in this series), and a static state-dependence specification will fold part of the forward-looking behavior into $\delta$. Shcherbakov (2016), who estimates switching costs in the pay-television market with a dynamic model of consumer behavior that includes switching costs, is one example of this route. Whether you need a dynamic model depends on whether the switching cost is large enough relative to the user's decision horizon to be worth looking ahead over.

Stringing these failure modes together, the credibility of switching-cost estimation ultimately rests not on the sophistication of the estimator but on two things: whether heterogeneity is modeled fully enough that the residual no longer contaminates $\delta$, and whether there is credibly exogenous price variation to calibrate switching. The permutation test, flexible heterogeneity, the initial-conditions correction, and the argument for an exogenous price source are all evidence provided around these two questions, and none of them can substitute for a substantive judgment about the mechanism of "why users stay."

## 8 Further Reading

::: {.readings}
Required reading, in suggested order:

- Dubé, Hitsch and Rossi (2010, RAND Journal of Economics). The methodological benchmark for empirical switching costs, the full version of Sections 3.6 and 4.5 of this chapter; focus on how the permutation test checks the heterogeneity specification.
- Handel (2013, American Economic Review). A model for the policy consequences of measuring inertia; focus on how forced re-enrollment calibrates the switching cost, and why reducing inertia instead worsens adverse selection.
- Armstrong (2006, RAND Journal of Economics). The core of two-sided pricing; focus on membership externalities and the competitive bottleneck, understanding why one side is subsidized.
- Rochet and Tirole (2003, Journal of the European Economic Association). The founding work on two-sided markets; focus on the argument that the price structure rather than the price level is the platform's decision.
- Manski (1993, Review of Economic Studies). The origin of the reflection problem; understand why the endogenous, contextual, and correlated effects are inseparable in the linear-in-means model.

Further reading:

- Katz and Shapiro (1985, American Economic Review). The classic source of fulfilled-expectations demand and tipping.
- Rysman (2004, Review of Economic Studies). An exemplar of structural estimation of two-sided platforms, with indirect network effects and welfare in the Yellow Pages market.
- Wooldridge (2005b, Journal of Applied Econometrics). The simple solution to the initial conditions problem, the theoretical basis for the workhorse of Section 4.3 of this chapter.
- Honoré and Kyriazidou (2000, Econometrica). The dynamic FE logit with no distributional assumption; read how the conditional likelihood cancels the fixed effect.
- Farrell and Klemperer (2007, Handbook of Industrial Organization). The authoritative survey of switching costs and network effects, with a full discussion of harvesting and investment.
- Bramoullé, Djebbari and Fortin (2009, Journal of Econometrics). The identification result that cracks the reflection problem using network structure.
- Dubé, Hitsch and Chintagunta (2010, Marketing Science). A structural model and counterfactual of tipping and concentration under indirect network effects.
:::

::: {.apa-refs}
- Armstrong, M. (2006). Competition in two-sided markets. *The RAND Journal of Economics, 37*(3), 668-691. https://doi.org/10.1111/j.1756-2171.2006.tb00037.x
- Bramoullé, Y., Djebbari, H., & Fortin, B. (2009). Identification of peer effects through social networks. *Journal of Econometrics, 150*(1), 41-55. https://doi.org/10.1016/j.jeconom.2008.12.021
- Cabral, L. (2011). Dynamic price competition with network effects. *The Review of Economic Studies, 78*(1), 83-111. https://doi.org/10.1093/restud/rdq007
- Caillaud, B., & Jullien, B. (2003). Chicken & egg: Competition among intermediation service providers. *The RAND Journal of Economics, 34*(2), 309-328. https://doi.org/10.2307/1593720
- Chamberlain, G. (1980). Analysis of covariance with qualitative data. *The Review of Economic Studies, 47*(1), 225-238. https://doi.org/10.2307/2297110
- Dubé, J.-P., Hitsch, G. J., & Chintagunta, P. K. (2010). Tipping and concentration in markets with indirect network effects. *Marketing Science, 29*(2), 216-249. https://doi.org/10.1287/mksc.1090.0541
- Dubé, J.-P., Hitsch, G. J., & Rossi, P. E. (2010). State dependence and alternative explanations for consumer inertia. *The RAND Journal of Economics, 41*(3), 417-445. https://doi.org/10.1111/j.1756-2171.2010.00106.x
- Einav, L., Farronato, C., & Levin, J. (2016). Peer-to-peer markets. *Annual Review of Economics, 8*(1), 615-635. https://doi.org/10.1146/annurev-economics-080315-015334
- Farrell, J., & Klemperer, P. (2007). Coordination and lock-in: Competition with switching costs and network effects. In M. Armstrong & R. Porter (Eds.), *Handbook of industrial organization* (Vol. 3, pp. 1967-2072). Elsevier.
- Farrell, J., & Saloner, G. (1986). Installed base and compatibility: Innovation, product preannouncements, and predation. *American Economic Review, 76*(5), 940-955.
- Gallus, J. (2017). Fostering public good contributions with symbolic awards: A large-scale natural field experiment at Wikipedia. *Management Science, 63*(12), 3999-4015. https://doi.org/10.1287/mnsc.2016.2540
- Handel, B. R. (2013). Adverse selection and inertia in health insurance markets: When nudging hurts. *American Economic Review, 103*(7), 2643-2682. https://doi.org/10.1257/aer.103.7.2643
- Heckman, J. J. (1981). The incidental parameters problem and the problem of initial conditions in estimating a discrete time-discrete data stochastic process. In C. F. Manski & D. McFadden (Eds.), *Structural analysis of discrete data with econometric applications* (pp. 179-195). MIT Press.
- Honoré, B. E., & Kyriazidou, E. (2000). Panel data discrete choice models with lagged dependent variables. *Econometrica, 68*(4), 839-874. https://doi.org/10.1111/1468-0262.00139
- Katz, M. L., & Shapiro, C. (1985). Network externalities, competition, and compatibility. *American Economic Review, 75*(3), 424-440.
- Katz, M. L., & Shapiro, C. (1994). Systems competition and network effects. *Journal of Economic Perspectives, 8*(2), 93-115. https://doi.org/10.1257/jep.8.2.93
- Manski, C. F. (1993). Identification of endogenous social effects: The reflection problem. *The Review of Economic Studies, 60*(3), 531-542. https://doi.org/10.2307/2298123
- Mundlak, Y. (1978). On the pooling of time series and cross section data. *Econometrica, 46*(1), 69-85. https://doi.org/10.2307/1913646
- Rochet, J.-C., & Tirole, J. (2003). Platform competition in two-sided markets. *Journal of the European Economic Association, 1*(4), 990-1029. https://doi.org/10.1162/154247603322493212
- Rochet, J.-C., & Tirole, J. (2006). Two-sided markets: A progress report. *The RAND Journal of Economics, 37*(3), 645-667. https://doi.org/10.1111/j.1756-2171.2006.tb00036.x
- Rysman, M. (2004). Competition between networks: A study of the market for Yellow Pages. *The Review of Economic Studies, 71*(2), 483-512. https://doi.org/10.1111/0034-6527.00512
- Rysman, M. (2009). The economics of two-sided markets. *Journal of Economic Perspectives, 23*(3), 125-143. https://doi.org/10.1257/jep.23.3.125
- Shcherbakov, O. (2016). Measuring consumer switching costs in the television industry. *The RAND Journal of Economics, 47*(2), 366-393. https://doi.org/10.1111/1756-2171.12131
- Weyl, E. G. (2010). A price theory of multi-sided platforms. *American Economic Review, 100*(4), 1642-1672. https://doi.org/10.1257/aer.100.4.1642
- Wooldridge, J. M. (2005b). Simple solutions to the initial conditions problem in dynamic, nonlinear panel data models with unobserved heterogeneity. *Journal of Applied Econometrics, 20*(1), 39-54. https://doi.org/10.1002/jae.770
:::
