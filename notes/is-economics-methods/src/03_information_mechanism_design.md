---
title: "Information Economics, Contracts & Mechanism Design"
subtitle: "Private Information, Incentives, and the Design of Rules"
seriesline: "Foundations of Information Systems Economics · Chapter 3"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 3 · Information Economics, Contracts & Mechanism Design"
---

## Introduction

A platform asks advertisers to report how much a click is worth to them, then hands the best ad slot to the highest bidder. The rule looks fair, except that it skips one detail: the advertiser knows its own value and the platform does not, and when the bid decides both the slot and the payment, telling the truth need not pay. An institution that asks participants to volunteer their private information, without first answering whether lying pays, tends to collect strategy rather than information.

Asymmetric information wears two faces. What is hidden before contracting is a type: low-quality sellers, high-risk borrowers, and high-value advertisers may all try to disguise themselves. What is hidden after contracting is an action: effort, maintenance, and content moderation cannot be fully monitored. Signaling, screening, and contracts each handle one of these difficulties, while mechanism design turns the problem around: instead of taking the rules as given and predicting an equilibrium, it starts from the outcome we want and looks for rules under which self-interested participants cooperate of their own accord. Incentive compatibility and individual rationality are therefore not side conditions but the boundary within which an institution can work at all.

The chapter proceeds as follows. Section 1 opens with a miscalculated number in Vireo's ad slot auction. Section 2 writes down the general form of a mechanism and the target quantities. Section 3 takes up, in turn, adverse selection, signaling, screening, moral hazard, the revelation principle and the Myerson lemma, VCG, information design, and cheap talk. Section 4 turns to auctions, bilateral trade, and matching. Section 5 reads three anchoring papers. Section 6 settles the accounts on Vireo's five data blocks: V0, a hand-computed block; V1, ten thousand ad exchange auctions; V2, a SaaS usage menu; V3, a linear contract for ad reviewers; and V4, a 200 × 200 talent market. Section 7 discusses failure modes. The question running through the chapter is this: the platform cannot order participants to hand over their true information, but can it arrange payments, menus, or signals so that telling the truth happens to be in their own interest?

## 1 An Auction That Cannot Extract the Truth

::: {.case}
The search platform Vireo wants to sell the ad slots on a page of search results to advertisers. To shrink the problem to a size that can be done by hand, take a concrete setup: there are two ad slots, the top slot brings 100 clicks a month and the second slot 40 clicks (the higher the position, the more clicks; this gap is called the position effect). Three advertisers have private per-click values of USD 10, 6, and 4, and each value is known only to the advertiser itself. Vireo wants to give the slots to the advertisers who value them most, and to collect as much money as possible.
:::

The product manager proposes a straightforward scheme: each of the three advertisers reports the price it is willing to pay per click, the highest report gets the top slot and the next-highest the second slot, and each pays the price it reported. This is the pay-your-bid, first-price rule. He estimates the revenue on the spot: the advertisers will report their true values 10, 6, and 4, so the top slot pays $10 \times 100 = 1000$ dollars and the second slot $6 \times 40 = 240$ dollars, for USD 1240 a month.

Buried in this estimate is a fatal assumption: that the advertisers bid truthfully. When each pays its own bid, an advertiser with value 10 that really bid 10 would win the top slot but hand over all of its surplus, for a net payoff of zero; so it shades its bid down, to the point where it still just wins while keeping as much surplus as possible. How far it shades depends on how it guesses its rivals will bid, and the rivals are running the same calculation, so bidding becomes a game of mutual second-guessing in which nobody's bid equals anybody's value. The USD 1240 can neither be computed accurately nor reveal what the product manager really wants to know: who values these slots most.

What Vireo needs is a set of rules under which "report the true value" is itself each advertiser's best choice, whatever the others report. Does such a rule exist? How much money does it collect, and is it the one that collects the most? Section 3.6 gives the VCG mechanism, which makes truthful bidding a dominant strategy, and Section 6.1 applies it to this example: it should collect only USD 680, far less than USD 1240, but every cent of it is collected on truthful bids. Vireo's ad exchange also runs a large number of auctions every day, and there the platform sees only the bids. Section 4 answers the questions that arise there: whether switching formats collects more money, how high the reserve price should be, and whether values can be recovered from bids.

## 2 The Economic Model and the Target Quantities

To analyze private information, first write the mechanism down in full. There are $n$ participants; participant $i$'s private type $\theta_i$ takes values in a type space $\Theta_i$ and encodes information only it knows: a value, a quality, a cost; the joint distribution of types is common knowledge. Throughout this chapter $\theta$ denotes a private type and $\Theta_i$ a type space; in the later chapters of the book $\theta$ and $\Theta$ mean something else. A mechanism specifies the messages each participant can send, that is, a strategy set $\mathcal{S}_i$ (in the notation of Chapter 2), and the outcome under each message profile $s = (s_1, \dots, s_n)$: the chosen alternative $a(s)$ (who gets what) and the money $t_i(s)$ each participant pays. The mechanism is written $(\mathcal{S}_1, \dots, \mathcal{S}_n; a, t)$; participants choose messages according to their types, and the outcome of the mechanism is the equilibrium outcome of this game, with the types as its inputs. A mechanism with $\mathcal{S}_i = \Theta_i$ is called a direct mechanism: participants only report types, the profile of reports is written $\theta' = (\theta_1', \dots, \theta_n')$, and the prime marks a report.

Two sets of constraints run through all of mechanism design. The first is incentive compatibility (IC): the mechanism leads each participant, acting out of self-interest, to make the choice the designer wants; in a direct mechanism, this means that truthful reporting is no worse than any misreport. When this holds whatever the others report, the mechanism is dominant-strategy incentive compatible (DSIC); when it holds in expectation, and only given that the others report truthfully, it is Bayesian incentive compatible (BIC). The second is individual rationality (IR, also called the participation constraint): participating must be no worse than staying out, or the participant walks away. Almost every mechanism design problem optimizes some objective subject to IC and IR.

Most results of the chapter need only a narrower environment: each participant's private information is one number, and money enters payoffs linearly. The assumptions below fix this environment.

::: {.assumption}
**I.1 (independent private values; MWG 23.D, p. 887; Proposition 23.D.3, p. 890)** Participants' types are mutually independent; $\theta_i$ takes values in an interval $[\underline\theta_i, \bar\theta_i]$, $\underline\theta_i < \bar\theta_i$, with a density that is positive everywhere on the interval; participant $i$'s payoff depends only on its own type.

**I.1' (symmetry; MWG Example 23.B.5, p. 865; Guerre, Perrigne and Vuong 2000)** On top of I.1, participants' types are identically distributed, with distribution function $F$ and density $f$.

**I.2 (single-parameter quasilinearity; MWG 23.D, p. 887; Example 23.B.4, p. 862)** Participant $i$'s payoff is $u_i = \theta_i q_i - t_i$, where $q_i$ is the quantity $i$ receives or its probability of receiving the object, and $t_i$ is the money it pays.
:::

I.1 restricts the information structure: independence keeps a participant's beliefs about the others' types from varying with its own type, so "how much a given report gets and how much it pays" depends only on what was reported, and the Myerson lemma of Section 3.5 rests on exactly this; the interval support and the positive density make integration by parts valid. The identical distribution of I.1' is used only in the symmetric formulas (first-price bids, the GPV inversion, a uniform reserve price). The typical failures of I.1 are common values and correlated values (Section 7). I.2 has two layers: money enters payoffs linearly, that is, risk neutrality and no budget constraint; and the type multiplies quantity linearly, that is, the private information is one number (single-parameter). Quasilinearity alone is not enough; the type must also enter linearly. Risk-averse bidders, advertisers with budget caps, and customers with private information along two dimensions all violate I.2. In a single-parameter environment the type is the per-unit value, so from Section 3.5 on we write $\theta_i = v_i$, with support $[\underline v_i, \bar v_i]$; a direct mechanism is written as an allocation rule $q_i(v)$ and a payment rule $t_i(v)$, $v = (v_1, \dots, v_n)$.

The target quantities come in three kinds: efficiency, that is, giving things to those who value them most and maximizing total surplus $\sum_i v_i q_i$, with VCG as the leading example; revenue, with the optimal auction as the leading example, which often sacrifices a little efficiency; and quantities tied to the specific problem, such as the probability of trade, the probability that persuasion succeeds, and the stability of a matching. The primitives of moral hazard differ slightly: what is private is not a type but an action, usually effort $e$ with cost $c(e)$. Effort shifts the distribution of performance $x$ (density $f(x \mid e)$) but is itself unseen, so the principal can only pay by performance, $w(x)$. The objective is the principal's expected profit, and IC becomes "make the agent choose, of its own accord, the effort the designer wants." As in Chapter 1, fix the target quantity first and discuss implementation second: "is this mechanism any good" must first be translated into "efficiency or revenue, for which type distribution, under which set of constraints."

## 3 Core Structure: From Adverse Selection to Mechanism Design

This section is written as in Chapter 1: first Vireo's numbers show what is to be proved, then come numbered assumptions and theorems, proof sketches that name the assumption each step uses, and an "In words" paragraph that answers five questions. Assumptions carry the chapter's prefixes I. (private information), K. (contracts), and MT. (matching); theorems are numbered consecutively within the chapter, with MWG's section, proposition, and page numbers in parentheses, and results that are not in the book cite the original papers. Sections 3.1 to 3.4 are ordered by the kind of information asymmetry and by who moves first, Sections 3.5 and 3.6 are the two general results of mechanism design, and Sections 3.7 and 3.8 treat information itself as a lever.

### 3.1 Adverse Selection and Market Unraveling

Vireo runs a used-equipment channel. The seller knows the quality $\theta$ of the equipment, while the buyer knows only that $\theta$ is uniformly distributed on $[0, 1]$; the seller values equipment of quality $\theta$ at $\theta$, and the buyer, who puts it to better use, values it at $1.5\theta$. Under full information every piece of equipment should trade. But the buyer cannot see quality: at an asking price $p$ the only sellers willing to sell are those with $\theta \leq p$, whose average quality is $p/2$; the buyer rationally anticipates this and is willing to pay only $1.5 \times (p/2) = 0.75p$, which is below the asking price at every positive price.

::: {.theorem}
**Theorem 1 (market unraveling under adverse selection; MWG Definition 13.B.1, p. 439; Example 13.B.1, p. 442; Akerlof 1970)** Let quality $\theta$ be uniformly distributed on $[0, 1]$ and known only to the seller; the seller's reservation value is $\theta$; buyers are risk neutral and compete with one another, and they value quality $\theta$ at $\lambda\theta$, $1 < \lambda < 2$ ($\lambda$ is used only in this section). A competitive equilibrium is a price $p$ and a set of sellers willing to trade, $\{\theta : \theta \leq p\}$, such that $p = \lambda\,\mathbb{E}[\theta \mid \theta \leq p]$. Then the unique equilibrium price is $p = 0$: only the lowest-quality equipment trades, even though every piece of equipment is worth more in the buyer's hands.
:::

The proof sketch has two steps. Step one (sellers decide whether to participate on the basis of their private information): for $p \leq 1$ the sellers willing to sell are those with $\theta \leq p$, and the uniform distribution gives an average quality of $p/2$; for $p > 1$ everyone sells, and the average quality is $1/2$. Step two (competition among buyers pushes the price to the expected value): $p > 1$ requires $p = \lambda/2 < 1$, a contradiction; $p \leq 1$ requires $p = \lambda p/2$, which for $\lambda < 2$ leaves only $p = 0$. MWG writes the same example with the reservation value equal to $\alpha$ times productivity, and $\alpha = 1/\lambda > 1/2$ is exactly the unraveling condition; Akerlof's original example has $\alpha = 2/3$ (MWG footnote 10), that is, this section's $\lambda = 1.5$.

In words: in a market with asymmetric information about quality, where price is the only contractible term, which types trade is decided by a fixed point: the average value of the sellers a price attracts must live up to that price. The uniform distribution with $\lambda < 2$ is the most extreme case, in which the whole market unravels. With a differently shaped distribution or larger gains from trade, the usual result is partial unraveling: the set of trades is cut off at the high-quality end, the goods observed to trade are of lower average quality than the population, and there may be more than one equilibrium (MWG Figure 13.B.3, p. 442). The result is determined by two primitives, the dispersion of quality and the multiple of the gains from trade, $\lambda$: the smaller the gains relative to the dispersion of quality, the more complete the unraveling. What the data can test are selection effects such as "at the same price, traded goods are worse than untraded ones," and identifying them takes exogenous variation in buyers' information, such as the introduction of disclosure, certification, or mandatory warranties; Lewis (2011), in Section 5.3, studies exactly how seller disclosure changes prices on eBay Motors. Vireo's numbers: at every price step from 0.1 to 1, buyers are willing to pay only $0.75p$ for the sellers who stay, which is 0.75 at $p = 1$; for the platform to do business in this channel, it must invest in ratings, certification, and guarantees, the information mechanisms that fight adverse selection.

### 3.2 Signaling: The Informed Party Speaks Up

There are two directions for fighting adverse selection, depending on who breaks the information asymmetry. The first is for the party holding the private information to do something to prove itself; this is signaling, and the education model of Spence (1973) is the prototype. The version in the used-equipment channel is a warranty: the seller of good equipment will almost never be asked to pay out, so a long warranty is cheap for that seller, while it is expensive for the seller of bad equipment; buyers cannot see quality, but they can see the length of the warranty. Below, the type $\theta$ of the seller (the sender) is the equipment's value to the buyer, the signal $s \geq 0$ is the warranty length in years, the signal costs $c(s, \theta)$, and the buyer (the receiver) bids $w(s)$ on the basis of the signal; in the model the signal itself creates no value and serves only to tell types apart.

::: {.assumption}
**I.3 (single crossing; MWG 13.C, pp. 450, 453)** The signal cost satisfies $c(0, \theta) = 0$, is increasing and strictly convex in $s$ and decreasing in $\theta$, and has $\partial^2 c(s, \theta)/\partial s\,\partial\theta < 0$: the higher the type, the lower the marginal cost of one more unit of signal.
:::

Geometrically, I.3 says that in the (signal, bid) plane the indifference curves of two types cross at most once, the high type's being flatter, so there is a range of signals the high type is willing to pay for and the low type is not. It rules out signals that are equally expensive for everyone, such as a badge anyone can pick up at zero cost. The payoff $\theta_i q_i - t_i$ of I.2 has the same property (MWG footnote 32, p. 888), and the screening of Section 3.3 and the Myerson lemma of Section 3.5 use the same single crossing.

The equilibrium concept for signaling games is the weak PBE of Theorem 7 (ii) in Section 3.6 of Chapter 2: strategies are sequentially rational given beliefs, and beliefs on the path, $\mu(\theta \mid s)$, follow from Bayes' rule. For the education model with two firms bidding, MWG adds one more condition: the two receivers hold the same belief after every signal (off the path included) and bid as in a Nash equilibrium given that belief, so competition pushes the bid to the expected value $w(s) = \mathbb{E}_\mu[\theta \mid s]$ (MWG 13.C, p. 452). Below, PBE means a weak PBE with this condition added.

::: {.theorem}
**Theorem 2 (Spence signaling: separating equilibria, least-cost separation, and the intuitive criterion; MWG Lemmas 13.C.1 and 13.C.2, pp. 453-454; pp. 455, 458; equation (13.AA.4), p. 471; Cho and Kreps 1987)** There are two types $\theta_H > \theta_L > 0$, each with positive prior probability; the signal creates no value, I.3 holds, and only pure strategies are considered; two receivers bid under a common belief, $w(s) = \mathbb{E}_\mu[\theta \mid s]$. (i) In every separating PBE, $w(s(\theta_L)) = \theta_L$, $w(s(\theta_H)) = \theta_H$, and $s(\theta_L) = 0$; $s(\theta_H)$ lies in $[\underline s, \bar s]$, where $\theta_H - c(\underline s, \theta_L) = \theta_L$ and $\theta_H - c(\bar s, \theta_H) = \theta_L$, and every value in the interval can be supported as a separating PBE by suitable off-path beliefs. (ii) The least-cost separating equilibrium, with $s(\theta_H) = \underline s$, Pareto dominates the other separating equilibria. (iii) The only PBE outcome that satisfies the intuitive criterion is least-cost separation.
:::

The proof sketch has four steps. Step one (competition among receivers and Bayes' rule): under separation the signals on the path reveal the type, and competition pushes the bid to that type's value. Step two ($c$ increasing in $s$): the low type gets only $\theta_L$ from any signal on the path, so it sends $s = 0$. Step three (I.3): for the low type not to mimic requires $\theta_H - c(s(\theta_H), \theta_L) \leq \theta_L$, that is, $s(\theta_H) \geq \underline s$; for the high type not to switch to $0$ requires $\theta_H - c(s(\theta_H), \theta_H) \geq \theta_L$, that is, $s(\theta_H) \leq \bar s$; single crossing makes the same signal cheaper for the high type, so $\underline s < \bar s$; with off-path beliefs that always treat the sender as the low type, nobody deviates. Lowering $s(\theta_H)$ only saves the high type's cost and leaves everyone else's payoff unchanged, which is (ii). Step four (I.3, (iii)): for signals between $\underline s$ and $s(\theta_H)$, even if the low type got the highest bid $\theta_H$, it would do worse than its equilibrium payoff $\theta_L$, so these signals are equilibrium-dominated for it. The intuitive criterion then requires beliefs to put weight only on the high type, and the high type strictly gains by switching to a signal slightly above $\underline s$, so only $\underline s$ survives; a pooling outcome is broken in the same way by a signal that only the high type is willing to send (MWG p. 458). MWG states the result only for two types; with more than two types the intuitive criterion need not leave a single outcome.

In words: under I.3, the signaling equilibria form a set: the separating outcomes trace out a whole segment indexed by $s(\theta_H) \in [\underline s, \bar s]$, and there are pooling outcomes besides; the intuitive criterion keeps only least-cost separation, the outcome usually called the Riley outcome. In welfare terms the signal is pure waste: to avoid being mimicked, the high type pays $c(\underline s, \theta_H)$, which creates no value at all. When single crossing fails (a badge everyone can afford), separation is no longer possible and only pooling remains. Where in the equilibrium set the outcome falls is decided by off-path beliefs, a refinement choice of the analyst; $\underline s$ is determined by the low type's cost of mimicking. In the data, identifying the effect of a signal requires exogenously changing its cost or its visibility, which on a platform usually means a randomized experiment (Chapter 25). A hand-computed example: with $\theta_L = 1$, $\theta_H = 2$, and $c(s, \theta) = s^2/\theta$, we get $\underline s = 1$ and $\bar s = \sqrt 2$; in least-cost separation the good seller offers a 1-year warranty at a cost of $1/2$, the bad seller offers none, and this $1/2$ is the price of the information asymmetry.

### 3.3 Screening: The Uninformed Party Designs a Menu

The second direction runs the other way: the party without the information designs a self-selection menu that induces each type to reveal itself; this is screening. It differs from signaling only in who moves first. In MWG, "screening" refers to two different models: 13.D has several firms competing in menus (competitive screening, where an equilibrium may fail to exist, MWG pp. 465-466), while 14.C and 23.F have a single monopolist designer offering a menu (monopolistic screening). This section and V2 are both about the latter.

Vireo sells a SaaS API product, tiered by monthly usage. Start with two kinds of customers: a low-usage customer values a thousand calls at $\theta_L = 4$ dollars and a high-usage customer at $\theta_H = 7$ dollars, and the share of high types is $\mu = 0.4$ (in this section $\mu$ denotes only the share of high types); a customer who buys $q$ thousand calls and pays $t$ dollars gets $\theta q - t$, and serving $q$ thousand calls costs Vireo $q^2/2$ dollars. If it could see types, Vireo would give each type the efficient usage $q = \theta$ and take all the surplus: $(4, 16)$ for the low type and $(7, 49)$ for the high type, earning USD 14.6 per customer on average. When types are unseen this menu fails: the high type would switch to the low type's tier and get $7 \times 4 - 16 = 12$ dollars. The optimal menu makes the low type's IR and the high type's IC bind exactly; the high type gets the efficient usage of 7, and the low type's usage is pushed down to

$$q_L = \theta_L - \frac{\mu}{1 - \mu}\,(\theta_H - \theta_L) = 4 - \frac{0.4}{0.6} \times 3 = 2$$

The low type pays USD 8 and keeps no surplus; the high type pays USD 43 and keeps an information rent of $(\theta_H - \theta_L)\,q_L = 6$ dollars, which it takes away by virtue of its private information; average profit falls to USD 11. To hold down the information rent owed to the high type, Vireo deliberately distorts the low end: one more unit of usage for the low type gives the high type an extra $\theta_H - \theta_L$ when it mimics, and Vireo would have to concede that much more information rent. This structure of "distort the bottom, extract at the top" returns when Chapter 4 treats second-degree price discrimination.

![The V2 two-type menu: the participation lines of the low and the high type, $t = 4q$ and $t = 7q$, and the high type's indifference line through the two bundles, $t = 7q - 6$. The low type gets the pushed-down bundle B ($q = 2$, $t = 8$) and sits on its own participation line; the high type gets bundle A with the efficient usage ($q = 7$, $t = 43$), USD 6 below its own participation line at $q = 7$, and this vertical gap is the information rent; the figure also marks the bundles $(4, 16)$ and $(7, 49)$ that would be offered if types were visible.](assets/fig/fig_03_screening.svg)

Two types are just two grid points. The continuous-type result uses the virtual value $\Psi(\theta) = \theta - \big(1 - F(\theta)\big)/f(\theta)$: the marginal revenue the seller can collect from type $\theta$, which falls short of $\theta$ by an information-rent correction. The marginal cost is $q$, so the $q_L$ above is exactly the low type's virtual value under the two-point distribution.

::: {.theorem}
**Theorem 3 (monopolistic screening with continuous types: monotonicity plus the envelope condition, virtual surplus, no distortion at the top, and ironing; Mussa and Rosen 1978; Maskin and Riley 1984; MWG Proposition 23.D.2, p. 888; Example 23.F.1, pp. 900-903; Exercise 23.F.2, p. 923)** A customer's type $\theta$ lies in $[\underline\theta, \bar\theta]$, with distribution function $F$ and a density $f$ that is positive everywhere (I.1); the payoff is $\theta q - t$ (I.2), and zero from not buying; the seller offers a menu $\{(q(\theta), t(\theta))\}$, with a strictly convex serving cost $C(q)$, $C(0) = C'(0) = 0$ ($C$ is used only in this section and Section 6.4). Write $U(\theta) = \theta q(\theta) - t(\theta)$. (i) The menu is incentive compatible if and only if $q(\theta)$ is nondecreasing and $U(\theta) = U(\underline\theta) + \int_{\underline\theta}^{\theta} q(z)\,dz$. (ii) The seller's expected profit equals $\mathbb{E}\big[\Psi(\theta)\,q(\theta) - C(q(\theta))\big] - U(\underline\theta)$; if $\Psi$ is increasing (regular), the optimal menu sets $U(\underline\theta) = 0$ and, for each $\theta$, chooses $q(\theta)$ to maximize the virtual surplus $\Psi(\theta)\,q - C(q)$: types with $\Psi(\theta) \leq 0$ are excluded, and the rest satisfy $C'(q(\theta)) = \Psi(\theta)$. (iii) No distortion at the top: $\Psi(\bar\theta) = \bar\theta$, so the highest type gets the efficient quantity; every other type served gets less than the efficient quantity. (iv) When $\Psi$ is not monotone, replace $\Psi$ in (ii) by the ironed virtual value $\bar\Psi$: the slope, in quantile space, of the convex hull of $\int \Psi\,dF$ (the largest convex function that lies below it); on an interval where $\bar\Psi$ is flat, all types get the same quantity (bunching).
:::

The proof sketch has four steps. Step one (I.2, the type multiplies quantity): type $\theta$ mimicking $\theta'$ gets $U(\theta') + (\theta - \theta')\,q(\theta')$; swapping the two types gives the bounds $(\theta - \theta')\,q(\theta') \leq U(\theta) - U(\theta') \leq (\theta - \theta')\,q(\theta)$, so $q$ is nondecreasing; letting $\theta' \to \theta$ gives $U' = q$, and integrating gives the envelope formula; conversely, when $q$ is nondecreasing, $U(\theta) - U(\theta') = \int_{\theta'}^{\theta} q \geq (\theta - \theta')\,q(\theta')$ and IC holds, which is (i). Step two (I.1, integration by parts): substituting $t = \theta q - U$ and the envelope formula into expected profit and integrating by parts gives $\mathbb{E}[U(\theta)] = U(\underline\theta) + \mathbb{E}\big[q(\theta)\big(1 - F(\theta)\big)/f(\theta)\big]$, which is the profit expression of (ii); $U(\underline\theta)$ only lowers profit, and IR lets it be zero. Step three ($\Psi$ increasing): the $q$ obtained by maximizing virtual surplus pointwise increases with $\Psi$, hence is nondecreasing and satisfies (i) automatically, so the solution found while ignoring the monotonicity constraint is the optimum; $F(\bar\theta) = 1$ gives (iii). Step four (ironing): when $\Psi$ is not monotone, the pointwise solution is no longer monotone and the monotonicity constraint binds; replacing $\int \Psi\,dF$ by its convex hull amounts to averaging $\Psi$ out over the intervals where monotonicity fails; for the rigorous argument, see the treatment of auctions in Myerson (1981); this book only states it.

In words: under I.1, I.2, and regularity, the theorem characterizes the profit-maximizing menu: who is excluded (types with $\Psi < 0$), how much each type gets, and how much information rent each keeps; the information rent $U(\theta) = \int q$ is uniquely determined by the allocation, and it is what the designer must give up to make high types tell the truth. Distortion runs only downward: every type served, except at the very top, gets less than the efficient quantity, by an amount set by $(1 - F)/f$. When regularity fails, ironing is required and a range of types is bunched on the same tier; a menu without ironing lets some types gain by buying another type's tier; with multidimensional types there is no comparably clean characterization (Section 7). The result is determined jointly by the designer's objective and the type distribution: with the same distribution, a designer who maximizes total surplus does not distort. The data show the shape of the menu, the share choosing each tier, and bunching; a customer's choice among a few tiers is a discrete choice problem, and the logit of Chapter 27 estimates preferences from such data; identification needs exogenous variation in prices or tiers. The V2 numbers (Section 6.4): types range from USD 2 to 10 per thousand calls; under the regular distribution, types below 4.632 are excluded, the top gets $q = 10$, exactly the efficient quantity, and profit is USD 10.29, only a little more than half of the efficient total surplus of USD 19.55, and the information rent is USD 5.54. The bimodal distribution needs ironing, the types between 3.78 and 6.94 are bunched on a single tier at $q = 1.84$, and without ironing some type could gain USD 5.29 by buying another type's tier.

### 3.4 Moral Hazard: Incentivizing an Unseen Action

Everything above is hidden information before contracting. Moral hazard is a hidden action after contracting: a principal (such as the platform) hires an agent (such as a review team) to do a job; the agent's effort shifts the distribution of performance but is itself unseen, and the principal can only pay by performance. Start with an example that can be done by hand. The team can choose high or low effort, with high effort costing 5 and low effort 0; under high effort the probability of high performance (worth 100) is 0.8, under low effort only 0.2, and low performance is worth 0. The platform pays by performance, $w_H$ for high performance and $w_L$ for low, subject to limited liability ($w \geq 0$: it cannot charge the team). To induce high effort, high effort must be no worse for the team than low effort:

$$0.8\, w_H + 0.2\, w_L - 5 \;\geq\; 0.2\, w_H + 0.8\, w_L$$

Rearranging gives $w_H - w_L \geq 25/3$. Under limited liability the platform pushes $w_L$ down to 0, so $w_H = 25/3 \approx 8.33$, and the team's expected payoff is $0.8 \times 25/3 - 5 = 5/3 \approx 1.67$, strictly positive: this is a rent forced out by limited liability. When the platform induces high effort, its expected profit is $0.8 \times (100 - 25/3) = 220/3 \approx 73.33$; letting effort stay low yields only 20, and with observable effort it would be 75; the agency cost $5/3$ exactly equals that rent. In the example the team is risk neutral, so only the limited-liability rent appears; when the agent is risk averse there is also the risk-incentive tradeoff: performance is laced with luck, and the stronger the incentive, the more risk is loaded onto the agent.

Now let effort be continuous: effort $e$ with an increasing, convex cost $c(e)$, performance $x$ with density $f(x \mid e)$, pay $w(x)$, and the agent's payoff $u(w) - c(e)$. Tackling the inner optimization "which $e$ will the agent choose" head-on is hard; the standard move is the first-order approach: replace the whole incentive constraint by the agent's first-order condition $\int u(w(x))\,\partial_e f(x \mid e)\,dx = c'(e)$. The conditions under which the first-order approach is reliable depend on the distribution of performance.

::: {.assumption}
**K.1 (monotone likelihood ratio, MLRP; MWG 14.B, p. 485; for the form with continuous effort, see p. 504)** The likelihood ratio $\partial_e f(x \mid e)/f(x \mid e)$ is increasing in performance $x$: the higher the performance, the more it looks like the result of high effort rather than of luck.

**K.2 (convexity of the distribution function in effort, CDFC; Rogerson 1985)** The distribution function of performance, $F(x \mid e)$, is convex in $e$.
:::

K.1 says that high performance is evidence of effort; first-order stochastic dominance does not imply it (MWG Figure 14.B.1). K.2 says that effort raises the probability of high performance at a decreasing rate; this condition is strong, and MWG also notes that the known sufficient conditions are all demanding (p. 504); Jewitt (1988) gives another set of sufficient conditions.

::: {.theorem}
**Theorem 4 (the first-order approach and the informativeness principle; Rogerson 1985; Holmström 1979, Proposition 3, p. 84; MWG 14.B, pp. 487-488)** The agent's payoff is $u(w) - c(e)$, with $u$ strictly increasing and strictly concave, and effort is one-dimensional. (i) Under K.1 and K.2, the principal can replace the incentive constraint by the agent's first-order condition, and the contract found this way solves the original problem. (ii) Suppose the first-order approach is valid and the optimal contract among those that depend only on $x$ induces a unique, interior effort; let there be another contractible signal $y$, and let $(x, y)$ have joint density $f(x, y \mid e)$. There is a contract depending on $(x, y)$ that strictly Pareto dominates it if and only if $y$ is informative about effort, that is, $f(x, y \mid e)$ cannot be written as $g(x, y)\,h(x, e)$: $x$ is not a sufficient statistic for $(x, y)$ with respect to $e$.
:::

The proof sketch: (i) is only stated, and the proof is in Rogerson (1985). K.1 makes the optimal pay of the relaxed problem increase with performance, and K.2 then makes the agent's expected payoff concave in effort, so the first-order condition finds the global optimum; the rigorous argument has to settle whether the solutions of the relaxed and the original problems coincide, which takes considerable space. (ii) has two steps, following Holmström's proof. Step one ($y$ not informative, using the strict concavity of $u$): given any $w(x, y)$, for each $x$ take a certain amount whose utility equals the average of the utility of $w(x, y)$ under $g(x, y)$; the factorization gives the agent the same expected utility under the two contracts, so it chooses the same effort, and Jensen's inequality makes the certain amounts cheaper on average, so $y$ is worthless. Step two ($y$ informative, using the first-order approach): the optimality condition requires $1/u'(w)$ to be an increasing affine function of the likelihood ratio $\partial_e f/f$, and here the likelihood ratio depends on $y$, so a contract that looks only at $x$ violates this condition and can be improved upon (Holmström 1979, §5).

In words: when the first-order approach is valid, the theorem answers which performance measures belong in pay: any measure that still carries information about effort after controlling for $x$. A measure that creates no value should be written into the contract as long as it helps the principal separate luck from effort; a pure-noise measure only adds risk for the agent. When K.1 or K.2 fails, the "optimal contract" may not actually be optimal, and one must turn to the analysis of Grossman and Hart (1983), which does not rely on the first-order approach. The result is determined by the joint distribution of the signals: a measure unrelated to effort but correlated with the noise in performance is also informative, and with a positive correlation it enters the contract with a negative weight, which is where relative performance evaluation comes from. In the data, pay-performance sensitivity should change after a monitoring technology goes live; variation in the rollout date can be handled with the event studies of Chapter 19. The V3 numbers (Section 6.5): when only the report rate is added, its weight is $-0.5$; when only the peer benchmark is added, its weight is $-0.2632$, exactly $-0.6$ times the weight of 0.4386 that this contract puts on the reviewer's own performance; the pure-noise measure gets a weight of 0.

Linear contracts are the most widely used class, and in one particular setting they have a closed-form solution.

::: {.assumption}
**K.3 (the LEN setting; Holmström and Milgrom 1987)** Performance is $x = e + \varepsilon$, $\varepsilon \sim \mathcal{N}(0, \sigma^2)$; the agent has CARA utility over $w - c(e)$, with coefficient of absolute risk aversion $A$ (as in Chapter 1); the principal is risk neutral; contracts are restricted to be linear, $w(x) = w_0 + \beta x$, where $\beta$ is the contract slope; $c(e) = c'' e^2/2$, with $c''$ a constant.
:::

L, E, and N stand for the linear contract, exponential (CARA) utility, and normal noise. Linearity is a restriction, not a conclusion: Holmström and Milgrom (1987) prove that linear contracts are optimal in a dynamic model in which performance accumulates over time and the agent adjusts effort continuously; in a static model, linearity serves only to obtain a closed form. $c''$ is constant only under quadratic cost; with a general cost, the optimal slope solves a fixed point involving $c''(e^\ast)$.

::: {.theorem}
**Theorem 5 (closed-form solution of the linear contract; the LEN setting of Holmström and Milgrom 1987; (i) and (ii) follow by direct computation in the three steps below)** Under K.3, the agent's certainty equivalent must be no lower than zero. (i) The optimal slope is $\beta^\ast = 1/(1 + A\sigma^2 c'')$, and the agent chooses $e^\ast = \beta^\ast/c''$. (ii) If there is another signal $y = k e + \eta$ ($k$ is the effect of effort on $y$ and $\eta$ is noise; both are used only in this section and Section 6.5), $(\varepsilon, \eta)$ is jointly normal with covariance matrix $\Sigma$, and the contract is $w = w_0 + \beta_x x + \beta_y y$, then, writing $\boldsymbol{k} = (1, k)'$ and $\tilde\sigma^2 = 1/(\boldsymbol{k}'\Sigma^{-1}\boldsymbol{k})$, the optimal incentive power $\beta = \beta_x + k\beta_y$ is $\beta^\ast = 1/(1 + A\tilde\sigma^2 c'')$, with weights $(\beta_x, \beta_y)' = \beta^\ast\tilde\sigma^2\,\Sigma^{-1}\boldsymbol{k}$.
:::

The proof sketch has three steps. Step one (CARA and normality in K.3): with normally distributed pay, the agent's certainty equivalent is $w_0 + \beta e - c(e) - A\beta^2\sigma^2/2$. Step two ($c$ quadratic): the agent's first-order condition is $\beta = c'' e$, so incentives work only through $\beta$. Step three (the principal is risk neutral): the principal uses $w_0$ to push the certainty equivalent down to zero, its own expected profit becomes $e - c(e) - A\beta^2\sigma^2/2$, and substituting $e = \beta/c''$ and differentiating with respect to $\beta$ gives (i). (ii) is analogous: given the incentive power $\beta$, the risk the agent bears is a quadratic form in the weights; minimizing it under the linear constraint $\beta_x + k\beta_y = \beta$ is a minimum-variance weighting of the two signals, and its effect is to replace $\sigma^2$ by $\tilde\sigma^2$.

In words: under K.3, the optimal slope writes incentives and insurance into a single formula: the noisier the performance, the more risk averse the agent, and the faster the marginal cost of effort rises, the weaker the incentive; with no noise or a risk-neutral agent, $\beta^\ast = 1$ and the agent takes the whole marginal product. When normality or CARA fails, the formula no longer holds, and restricting attention to linear contracts may also lower the principal's profit. The result is determined by $A$, $\sigma^2$, $c''$, and which measures the principal writes into the contract; one more informative measure amounts to cutting the noise in performance from $\sigma^2$ to $\tilde\sigma^2$. In the data, the testable prediction is that, among comparable jobs, the noisier the performance, the less pay responds to it; but noisy jobs often differ in other ways as well, so identification needs exogenous variation in the noise itself. The V3 numbers: Vireo's ad reviewers have $\sigma^2 = 4$, $A = 0.5$, and $c'' = 1$ (amounts in thousands of dollars); paying by performance alone gives $\beta^\ast = 1/3$, and each reviewer brings the platform USD 166.67 of profit a month. Once the user report rate is written into the contract, $\tilde\sigma^2$ falls from 4 to 2, the incentive power rises from $1/3$ to $1/2$, yet the weight on performance itself falls from $1/3$ to $1/4$, as part of the pay moves onto the report rate, and profit rises to USD 250, an increase of USD 83.33, or 50%.

### 3.5 The Revelation Principle and the Myerson Lemma

The space of mechanisms is frighteningly large. V1's ad exchange runs first-price auctions, and with two bidders the bids average only 62.7% of values, with six bidders 88.4%. Could the platform simply have advertisers report their values directly and then place the equilibrium bids on their behalf? The result below says yes, and that nothing is lost by doing so.

::: {.theorem}
**Theorem 6 (the revelation principle; MWG Proposition 23.C.1, p. 871; Proposition 23.D.1, p. 884; Myerson 1979)** Given the mechanism $(\mathcal{S}_1, \dots, \mathcal{S}_n; a, t)$ and a strategy profile $s^\ast(\theta) = (s_1^\ast(\theta_1), \dots, s_n^\ast(\theta_n))$ of it, construct a direct mechanism from the composite rules $a(s^\ast(\cdot))$ and $t(s^\ast(\cdot))$. (i) If $s^\ast$ is a dominant-strategy equilibrium, then truthful reporting is a dominant strategy for every participant in this direct mechanism; (ii) if $s^\ast$ is a BNE, then truthful reporting is a BNE of the direct mechanism. In both cases, the outcome of the direct mechanism under truthful reporting is the same as that of the original mechanism.
:::

The proof sketch has two steps. Step one (the construction): reporting $\theta_i'$ in the direct mechanism is the same as having $s_i^\ast(\theta_i')$ carried out on one's own behalf in the original mechanism. Step two ($s^\ast$ is an equilibrium): if some type gained by misreporting in the direct mechanism, it would gain just as much by switching to $s_i^\ast(\theta_i')$ in the original mechanism, contradicting the fact that $s^\ast$ is an equilibrium; the dominant-strategy version holds whatever the others report, and the Bayesian version holds in expectation when the others report truthfully.

In words: any outcome a mechanism can implement in some equilibrium can be implemented by a direct mechanism with truthful reporting, so the search for the optimal mechanism can be confined to direct mechanisms $(q, t)$ that satisfy IC and IR; the reduction is lossless. It has two boundaries: the direct mechanism may also have untruthful equilibria, and the theorem guarantees only that truthful reporting is one of them; and the designer must be able to commit to carrying out the outcome rule, so if the platform would exploit the reported information after the fact, the reduction fails. It is a tool of analysis, not advice to actually run direct mechanisms: in practice, platforms pick first price or the generalized second-price (GSP) auction for simplicity and transparency. Its implication for empirical work is that the outcome is ultimately a function of types, so what is to be recovered is the type distribution, and the GPV method of Section 4.4 does exactly that. In V1, if the platform placed the equilibrium bids on behalf of the advertisers in two-bidder auctions, the bids would still be about 62.7% of values, with allocation and payments unchanged to the cent.

In a single-parameter environment the reduction goes one step further: payments are determined by the allocation.

::: {.theorem}
**Theorem 7 (the Myerson lemma, Bayesian version; MWG Proposition 23.D.2, p. 888; Myerson 1981)** Under I.1 and I.2, write the interim allocation and payment of a direct mechanism $(q, t)$ as $\bar q_i(v_i) = \mathbb{E}[q_i(v) \mid v_i]$ and $\bar t_i(v_i) = \mathbb{E}[t_i(v) \mid v_i]$, and the information rent as $U_i(v_i) = v_i\bar q_i(v_i) - \bar t_i(v_i)$. The mechanism is BIC if and only if, for every $i$: (i) $\bar q_i$ is nondecreasing; (ii) $U_i(v_i) = U_i(\underline v_i) + \int_{\underline v_i}^{v_i} \bar q_i(z)\,dz$.
:::

The proof sketch has four steps. Step one (the independence in I.1): the distribution of the others' types does not vary with $v_i$, so when type $v_i$ reports $v_i'$, its probability of winning and its expected payment are $\bar q_i(v_i')$ and $\bar t_i(v_i')$, which depend only on the report. Step two (the linearity in I.2): the interim payoff from reporting $v_i'$ is $U_i(v_i') + (v_i - v_i')\,\bar q_i(v_i')$, and BIC says that $U_i(v_i)$ is no lower than this. Step three (squeezing the incentive inequalities): swapping the two types gives $(v_i - v_i')\,\bar q_i(v_i') \leq U_i(v_i) - U_i(v_i') \leq (v_i - v_i')\,\bar q_i(v_i)$, so $\bar q_i$ is nondecreasing; letting $v_i' \to v_i$, the derivative of $U_i$ equals $\bar q_i$ almost everywhere, and integrating gives (ii). This is the envelope formula in integral form, which MWG derives directly from the incentive inequalities, without the differentiability of the value function that the envelope theorem of Chapter 1 assumes. Step four (sufficiency): when (i) and (ii) hold, $U_i(v_i) - U_i(v_i') = \int_{v_i'}^{v_i} \bar q_i \geq (v_i - v_i')\,\bar q_i(v_i')$, which is exactly BIC.

In words: in a single-parameter environment with independent private values, once the allocation rule is fixed, every type's expected payment is pinned down up to a single constant, the information rent of the lowest type: $\bar t_i(v_i) = v_i\bar q_i(v_i) - U_i(\underline v_i) - \int_{\underline v_i}^{v_i}\bar q_i(z)\,dz$. All the designer can really choose is the allocation and the treatment of the lowest type, and revenue equivalence, the optimal auction, and Myerson-Satterthwaite all start from this sentence. With correlated types, payments are no longer pinned down, and the designer can exploit the correlation to extract more (Section 7); with multidimensional types, monotonicity plus the envelope condition no longer characterizes incentive compatibility. The Myerson lemma as the market design literature usually states it is the dominant-strategy version: for every profile of the others' reports, $q_i$ is nondecreasing in one's own report, and payments are given pointwise by the same integral formula; this chapter uses the weaker Bayesian version, which suffices for the results of Section 4. For the data, it supplies the logic for recovering values from bids (Section 4.4). The V1 numbers: first price and second price both give the object to the bidder with the highest value, so $\bar q_i(v) = F(v)^{n-1}$ is the same, and the lowest type's rent is zero in both, so the theorem predicts that each type's expected payment is the same under the two formats; Section 4.2 tests this statement on ten thousand auctions.

### 3.6 VCG: Making Truth a Dominant Strategy

Now return to the question of Section 1: is there a mechanism that makes truthful bidding a dominant strategy while giving the allocation to those who value it most? There is: the Vickrey-Clarke-Groves (VCG) mechanism. The construction takes two sentences. Choose the alternative $a^\ast$ that maximizes the sum of reported values. Each participant pays the externality it imposes on everyone else, namely "the largest total value the others could have obtained without it" minus "the total value the others actually obtain under the chosen alternative." This way of pricing is called the Clarke pivot: you pay only when your presence changes the others' allocation, and what you pay is exactly the value you crowd out from the others.

::: {.theorem}
**Theorem 8 (VCG and Groves mechanisms; MWG Proposition 23.C.4, p. 877; the Clarke mechanism, p. 878; Propositions 23.C.5 and 23.C.6, pp. 879-880; Green and Laffont 1977)** A quasilinear private-values environment: participant $i$'s value for alternative $a$ is $v_i(a)$, and its reported value function is written $v_i'$. (i) If the alternative $a^\ast$ maximizes $\sum_j v_j'(a)$ and participant $i$ pays $t_i = h_i(v_{-i}') - \sum_{j \neq i} v_j'(a^\ast)$, where $h_i$ does not depend on $i$'s report (a Groves mechanism), then truthful reporting is a (weakly) dominant strategy for every participant; taking $h_i = \max_a \sum_{j \neq i} v_j'(a)$ gives the VCG mechanism, whose payments are nonnegative. (ii) If each participant's value can be any function on the set of alternatives, then the only efficient and dominant-strategy incentive compatible mechanisms are Groves mechanisms, and none of them is budget balanced at every type profile.
:::

The proof sketch: the proof of (i) is short and is given in full. Let participant $i$'s true value be $v_i(\cdot)$; the VCG payment is

$$t_i = \max_{a}\sum_{j \neq i} v_j'(a) \;-\; \sum_{j \neq i} v_j'(a^*)$$

Measured at true values, $i$'s net utility is

$$u_i = v_i(a^*) - t_i = \Big[v_i(a^*) + \sum_{j \neq i} v_j'(a^*)\Big] - \max_{a}\sum_{j \neq i} v_j'(a)$$

The second term on the right does not involve $i$'s report and is a constant beyond its control, so $i$ need only maximize the bracketed term $v_i(a^\ast) + \sum_{j \neq i} v_j'(a^\ast)$. And $a^\ast$ is the alternative the mechanism picks from the reports to maximize $v_i'(a) + \sum_{j \neq i} v_j'(a)$; under truthful reporting the mechanism's choice maximizes exactly the sum $i$ truly cares about, and a misreport can only lead the mechanism to pick an alternative that is worse for the reporter, or at best as good, so truthful reporting is dominant; the same holds for any $h_i$. Payments are nonnegative because the largest total value of the others without $i$ is no lower than their total value under any alternative. For the proof of the first half of (ii), see MWG pp. 879-880, and MWG does not prove the second half; this book only states it, because the proof has to rule out non-Groves payments one by one on a rich domain of values.

In words: with quasilinear private values, VCG prices by externality, turning each participant's net utility into total social welfare minus a constant beyond its control, so individual self-interest and efficiency are perfectly aligned. The second-price auction (Vickrey 1961) is its single-object special case: the winner pays the second-highest price, which is exactly the externality it imposes on the runner-up. How much is collected is determined by $h_i$, and the Clarke pivot makes payments nonnegative; (ii) says that dominant strategies plus efficiency leave only Groves mechanisms, at the price of a budget that is generally unbalanced, and in the bilateral trade of Section 4.5 VCG's imbalance becomes the platform's subsidy. When quasilinearity fails (budget constraints, for example), truth is no longer dominant; VCG is also vulnerable to shill bidding (one bidder bidding under several false identities) and to collusion (Section 7). For the data, dominant strategies mean that bids are values and need no inversion, provided participants understand and believe this. The V0 numbers (Section 6.1): the top slot pays USD 520, or USD 5.2 per click, and the second slot pays USD 160, or USD 4 per click, for a total revenue of USD 680; the two winners' net utilities, USD 480 and USD 80, are both nonnegative, so IR holds.

### 3.7 Information Design: Designing Information Rather Than Allocation

The mechanisms so far all design "who gets what and who pays how much," but the platform holds one more lever: it has information, and it can decide how information about quality or fit is revealed to the other side. Designing the information itself so as to steer the other side's choices is called information design, and its foundational framework is the Bayesian persuasion of Kamenica and Gentzkow (2011). Its key is commitment: the sender (such as the platform) first publicly commits to an experiment that maps the true state into a signal; the receiver (such as a buyer) sees the signal, updates by Bayes' rule, and then acts in its own interest. The sender cannot lie, but it can choose how much to reveal and how.

Vireo has a product recommendation slot. A product is in one of two states, good or bad; the user's prior is $\Pr(\text{good}) = 0.3$, the user buys only when the posterior is at least 0.5, and Vireo only wants the user to buy. With no information revealed, the posterior stays at 0.3 and the user does not buy; with full disclosure only good products are bought, with a purchase probability of 0.3. A better scheme: always recommend a good product, and recommend a bad product too with probability $x$ (in this section $x$ is not performance), pushing the posterior after a recommendation, $0.3/(0.3 + 0.7x)$, down to exactly 0.5; this gives $x = 3/7$, and the purchase probability rises to $0.3 + 0.7 \times 3/7 = 0.6$.

::: {.theorem}
**Theorem 9 (Bayes plausibility and concavification; Kamenica and Gentzkow 2011)** The state $\omega$ takes finitely many values, with prior $\mu_0$; the sender first commits to a signal structure (the conditional distribution of the signal given the state); the receiver sees the signal, forms the posterior $\mu$, and chooses the action that is best for itself, breaking ties in the sender's favor; the sender's payoff depends on the posterior only through this action and is written $\mathcal{V}(\mu)$, and $\mathrm{cav}\,\mathcal{V}$ is the smallest concave function no lower than $\mathcal{V}$. (i) A distribution of posteriors can be generated by some signal structure if and only if its mean equals $\mu_0$ (Bayes plausibility). (ii) The sender's optimal expected payoff (the supremum) is $\mathrm{cav}\,\mathcal{V}(\mu_0)$; the sender gains from disclosure if and only if $\mathrm{cav}\,\mathcal{V}(\mu_0) > \mathcal{V}(\mu_0)$.
:::

The proof sketch has three steps. Step one (Bayes' rule): the expected posterior equals the prior, which is necessity; conversely, given a distribution of posteriors with mean $\mu_0$, in state $\omega$ let the probability of sending the signal "the posterior is $\mu$" equal the weight of that posterior times $\mu(\omega)/\mu_0(\omega)$; Bayes' rule then yields exactly the posterior $\mu$, which is sufficiency. Step two (the definition of $\mathcal{V}$): the sender's expected payoff is the mean of $\mathcal{V}$ under the distribution of posteriors, so the problem becomes maximizing it over distributions of posteriors with mean $\mu_0$. Step three (the definition of the concave closure): over all ways of writing $\mu_0$ as a convex combination of posteriors, the supremum of the weighted average of $\mathcal{V}$ is $\mathrm{cav}\,\mathcal{V}(\mu_0)$; for the supremum to be attained, a regularity condition such as the upper semicontinuity of $\mathcal{V}$ is needed; see the original paper.

In words: provided the sender can commit and the receiver updates by Bayes' rule, optimal disclosure is neither "tell everything" nor "tell nothing" but a split of the prior into several posteriors, such that at the posteriors favorable to the sender the receiver is just willing to act; if $\mathrm{cav}\,\mathcal{V}(\mu_0) = \mathcal{V}(\mu_0)$, no disclosure beats disclosing nothing. Without commitment we are back to the cheap talk of Section 3.8; when the receiver does not update by Bayes' rule or has other sources of information, "the mean posterior equals the prior" no longer characterizes the feasible set. In the unified framework of Bergemann and Morris (2019), signals can be taken to be action recommendations that the receiver is willing to obey. The result is determined by the sender's objective and the receiver's action threshold, here the step at 0.5. The data must record both the disclosure policy and the actions; the predictions are in the table of Section 4.7, and testing them takes exogenous variation in the prior or the threshold. Vireo's numbers: $\mathcal{V}(\mu)$ jumps from 0 to 1 at 0.5, and the concave closure on $[0, 0.5]$ is the straight line $2\mu$, which takes the value 0.6 at the prior 0.3, matching the hand-computed purchase probability exactly; the probability that a bad product is slipped into the recommendations is $3/7 \approx 0.429$, and optimal disclosure doubles the purchase probability of full disclosure.

### 3.8 Cheap Talk: Soft Information Without Commitment

The power of Bayesian persuasion comes from commitment. A great deal of real-world information transmission has no commitment: a seller's description of a product and an analyst's recommendation are soft information in the sense of cheap talk, where messages are costless, unverifiable, and binding on no one after the fact. How much information can be transmitted depends on how aligned the two sides' interests are.

::: {.theorem}
**Theorem 10 (partition equilibria of cheap talk; Crawford and Sobel 1982)** The state $\omega$ is uniformly distributed on $[0, 1]$ and seen only by the sender; the receiver chooses a real-valued action; both sides have quadratic losses, the receiver's ideal action is $\omega$ and the sender's is $\omega + b$, $b > 0$ (in this subsection $b$ denotes only the bias and is unrelated to the bid function $b(v)$ of Section 4.1); messages are costless and unverifiable, and the sender cannot commit. Then every equilibrium is outcome-equivalent to a partition equilibrium: $[0, 1]$ is divided by $0 = \omega_0 < \omega_1 < \dots < \omega_N = 1$ into $N$ intervals, the sender reveals only which interval $\omega$ falls in, the receiver chooses the midpoint of that interval, and the endpoints satisfy $\omega_{i+1} - \omega_i = \omega_i - \omega_{i-1} + 4b$. An equilibrium with $N$ intervals exists if and only if $N(N - 1) < 1/(2b)$.
:::

The proof sketch has three steps. Step one (costless messages, no commitment): the sender always sends the message it likes best; payoffs satisfy single crossing, so the states that induce the same action form an interval, which gives the partition structure. Step two (the receiver acts by Bayes' rule): the posterior within an interval is uniform, and under quadratic loss the optimal action is the interval's midpoint. Step three (indifference at the endpoints): the sender at an endpoint state $\omega_i$ is indifferent between the two neighboring midpoints, $\omega_i + b = (\omega_{i-1} + 2\omega_i + \omega_{i+1})/4$, so each interval is longer than the previous one by $4b$; the recursion from $\omega_0 = 0$ gives $\omega_i = i/N + 2b\,i(i - N)$, and the first interval has positive length if and only if $N(N - 1) < 1/(2b)$.

In words: without commitment and with costless messages, how fine the transmitted information can be is determined by the bias $b$. With perfectly aligned interests, full truth-telling is possible; once there is a systematic bias, the sender always wants to exaggerate, and only the segment the state falls in can be revealed. Since the sender wants to push the action up, intervals must get wider toward the top for the sender at each endpoint to be indifferent. The babbling equilibrium ($N = 1$) always exists; equilibrium selection is the analyst's call, usually the one with the most intervals. With commitment we are back in Section 3.7, where the sender can actively shape beliefs; when messages can be verified, information is often forced out in full. In the data, the theory predicts coarser recommendations the less aligned the interests are, and testing it requires exogenous variation in the bias itself, such as a change in the commission structure. The V0 numbers: for biases of 0.01, 0.05, 0.1, 0.2, and 0.25, the largest numbers of intervals that can be sustained are 7, 3, 2, 2, and 1; with a bias of 0.05, the endpoints of the three intervals are 0, 2/15, 7/15, and 1, and each interval is longer than the one before by $4b = 3/15$.

## 4 Auctions, Bilateral Trade, and Matching

Auctions are the application of mechanism design closest to platform business. This section finishes the story of efficiency and revenue in auctions, then shows how values can be recovered from bids, and ends with bilateral trade and matching.

### 4.1 Four Basic Auctions and the First-Price Equilibrium

Single-object auctions come in four textbook forms: first-price sealed-bid (the highest bidder wins and pays its own bid), second-price sealed-bid (the Vickrey auction: the highest bidder wins and pays the second-highest bid), English (an open ascending auction in which the price rises until only one bidder remains), and Dutch (a descending auction in which the price falls from a high level until someone accepts). With private values, the Dutch auction is strategically equivalent to the first-price auction, and the English auction is outcome-equivalent to the second-price auction. In the second-price auction truthful bidding is dominant (Theorem 8); in the first-price auction, bidding the true value means zero surplus, so bidders must shade, and the equilibrium below says by how much.

::: {.theorem}
**Theorem 11 (the symmetric equilibrium of the first-price auction; Krishna 2010; for the special case $n = 2$ with uniform values, see MWG Example 23.B.5, p. 865)** Under I.1, I.1', and I.2, there are $n$ bidders and no reserve price; write $v_{-i,(1)} = \max_{j \neq i} v_j$ for the highest value among the other bidders. The symmetric, strictly increasing BNE bid is

$$b(v) = \mathbb{E}\big[v_{-i,(1)} \mid v_{-i,(1)} < v\big] = v - \frac{\int_{\underline v}^{v} F(z)^{n-1}\,dz}{F(v)^{n-1}}$$

When values are uniform on $[0, 1]$, $b(v) = \frac{n-1}{n}\, v$.
:::

The proof sketch has three steps. Step one (I.1, I.1'): when all other bidders use the increasing $b$, a bid of $b(z)$ wins with probability $F(z)^{n-1}$, and a bidder with value $v$ has expected payoff $\big(v - b(z)\big)F(z)^{n-1}$. Step two (equilibrium requires $z = v$ to be optimal): the first-order condition gives $\frac{d}{dv}\big[b(v)F(v)^{n-1}\big] = v\,\frac{d}{dv}F(v)^{n-1}$; integrating with $b(\underline v) = \underline v$ and then integrating by parts gives the formula above. Step three (single crossing from I.2): the payoff has increasing differences in the bid and the value, so a local optimum is a global optimum; this step can also be verified directly with the monotonicity plus envelope condition of Theorem 7.

In words: with symmetric independent private values, each bidder bids "the conditional expectation of the highest rival value given that it wins," which is exactly what it would pay on winning a second-price auction; the more bidders, the closer this conditional expectation is to its own value, and the less it shades. With asymmetric bidders there is generally no closed form, and the bidder with the highest value need not win; risk-averse bidders bid more. The bid is determined by $F$ and $n$, so $n$ must be exogenous and observable; when the number of participants is endogenous to the auction, entry has to be modeled as well, and Chapter 30 covers entry models. The numbers: uniform on $[0, 1]$ with two bidders, $b(v) = v/2$; in V1 the average shading ratio $(v - b)/v$ falls from 0.373 with two bidders to 0.116 with six.

### 4.2 Revenue Equivalence

Consider two auctions, one that cannot extract the truth and forces shading and one that can. Will their revenues be worlds apart? With values uniform on $[0, 1]$, the second-price auction's expected revenue is the expected second-highest value, $\frac{n-1}{n+1}$; the first-price auction's is the expected highest bid, $\frac{n-1}{n}\,\mathbb{E}[v_{(1)}] = \frac{n-1}{n} \cdot \frac{n}{n+1}$; the two are equal: for $n = 2$ both are $1/3$, and for $n = 3$ both are $1/2$. This is no coincidence.

::: {.theorem}
**Theorem 12 (revenue equivalence, a corollary of the Myerson lemma; MWG Proposition 23.D.3, p. 890)** Under I.1 and I.2 (bidders are risk neutral, and values are independent, each drawn from an interval with a density that is positive everywhere), if the BNEs of two auctions give each bidder the same probability of winning at every profile of values, and give each bidder's lowest type the same expected utility, then the two auctions bring the seller the same expected revenue.
:::

The proof sketch has three steps. Step one (Theorem 6): the equilibrium of each auction corresponds to a BIC direct mechanism. Step two (Theorem 7): each type's expected payment is pinned down by $\bar q_i$ and $U_i(\underline v_i)$. Step three (the interval and positive density in I.1): taking the expectation of the expected payment over $v_i$ and integrating by parts gives MWG equation (23.D.14):

$$\mathbb{E}\big[\bar t_i(v_i)\big] = \mathbb{E}\big[\Psi_i(v_i)\,\bar q_i(v_i)\big] - U_i(\underline v_i), \qquad \Psi_i(v) = v - \frac{1 - F_i(v)}{f_i(v)}$$

The right side contains only the allocation and the treatment of the lowest type; if both are the same, so is revenue.

In words: the seller's expected revenue is determined only by the allocation rule and the treatment of the lowest type, not by the surface rules of the auction, so if the platform changes the auction format without changing the allocation, expected revenue should not change. When the data show a difference, first ask which of I.1 and I.2 has broken: with risk-averse bidders the first-price auction collects more; with affiliated values the second-price and English auctions collect more (Milgrom and Weber 1982); with asymmetric bidders there is no general ranking. The choice of auction format should therefore turn on other dimensions: resistance to collusion, simplicity, resistance to manipulation. Testing the result requires randomizing formats across comparable auctions, while watching out for interference among the same set of bidders (Section 5.3 of Chapter 25). The V1 numbers: on the same values of ten thousand auctions, the first-price auction collects USD 1.4173 per auction on average and the second-price auction USD 1.4222, a difference of $-0.0049$, about 1.6 times its standard error; the true value is 1.4172 for both. Across 1,000 replications of the simulation, the difference has a mean of 0.00007 and a 2.5th to 97.5th percentile range of $[-0.0060, 0.0063]$, which covers zero.

### 4.3 Optimal Auctions and Reserve Prices

Revenue equivalence says only that the auctions that "give the object to the highest-value bidder" collect the same revenue; it does not say that they collect the most. Set a reserve price $r$ and sell to no one whose value falls below it: with values uniform on $[0, 1]$ and two bidders, taking $r = 1/2$ raises expected revenue from $1/3$ to $5/12 \approx 0.417$, at the cost of a $1/4$ probability that both values fall below $1/2$ and the object goes unsold.

::: {.theorem}
**Theorem 13 (optimal auctions and reserve prices; MWG Example 23.F.2, pp. 903-905, equation (23.F.9); Exercise 23.F.3, p. 923; Myerson 1981; Riley and Samuelson 1981)** Under I.1 and I.2, the seller values the object at zero, and every $\Psi_i$ is increasing (regular). The revenue-maximizing mechanism gives the object to the bidder whose $\Psi_i(v_i)$ is largest and nonnegative (and does not sell if all are negative), and leaves every lowest type an information rent of zero. Adding I.1', this is equivalent to: the bidder with the highest value gets the object, provided that value is no lower than the reserve price $r^\ast$, $\Psi(r^\ast) = 0$, that is, $r^\ast = \big(1 - F(r^\ast)\big)/f(r^\ast)$ (no reserve price when $\Psi(\underline v) \geq 0$); a second-price or a first-price auction with reserve price $r^\ast$ implements it, and $r^\ast$ does not depend on the number of bidders.
:::

The proof sketch has four steps. Step one (the revenue expression in the proof of Theorem 12): expected revenue equals $\mathbb{E}\big[\sum_i \Psi_i(v_i)\,q_i(v)\big] - \sum_i U_i(\underline v_i)$. Step two (IR): maximize $\sum_i \Psi_i q_i$ pointwise, giving the object to the bidder with the largest nonnegative $\Psi$; the lowest information rent that IR allows for the lowest type is zero. Step three (regularity): an increasing $\Psi_i$ makes $\bar q_i$ nondecreasing under this rule, so the monotonicity condition of Theorem 7 holds automatically; without regularity, ironing is needed (Theorem 3 (iv)). Step four (I.1'): when the $\Psi$ are the same and increasing, comparing $\Psi$ is comparing values, and $\Psi \geq 0$ means $v \geq r^\ast$; a second-price auction with $r^\ast$ implements this allocation with zero rent for the lowest type, and by Theorem 12 a first-price auction with the same reserve price collects the same revenue.

In words: under I.1, I.1', I.2, and regularity, the revenue-optimal auction is a standard auction with a reserve price equal to the monopoly price: facing the highest-value bidder, the seller is like a monopolist; raising the threshold loses the sale with some probability but collects more when the sale goes through, and $r^\ast$ balances the two. It sacrifices efficiency: bidders whose values are positive but below $r^\ast$ are shut out. Without regularity, ironing comes first; with asymmetric bidders, the optimal mechanism no longer gives the object to the highest-value bidder; and if the number of participants falls in response to the reserve price, part of the extra revenue is offset. The result is determined by what the seller knows about $F$ and by the revenue objective: $r^\ast$ depends only on $F$, not on $n$, but the extra money a reserve price can collect shrinks quickly as $n$ grows. Testing it requires first estimating $F$ (Section 4.4) and then randomizing the reserve price; the field experiment of Ostrovsky and Schwarz (2023), in Section 5.2, set reserve prices according to this theory. The V1 numbers: $r^\ast = 1.0615$ dollars per click; the reserve price raises the second-price auction's revenue by USD 0.0228 per auction on average (standard error 0.0019), against a true value of 0.0253, with a 2.5th to 97.5th percentile range of $[0.0216, 0.0291]$ across 1,000 replications of the simulation, and 3.4% of auctions go unsold. At true values, two-bidder auctions collect 0.0991 more, close to a tenth of their revenue, while auctions with three to six bidders collect only 0.0206, 0.0050, 0.0013, and 0.0004 more.

![V1 revenue: the left panel shows average revenue per auction against the number of bidders $n$ for four formats, first price, second price, and each with the reserve price $r^\ast$; the dots are sample means over the 2,000 auctions at each $n$, the vertical bars are the 2.5th to 97.5th percentiles across 1,000 replications of the simulation, and the curves are expected revenue with and without the reserve price; the right panel shows the revenue gain from the reserve price, about USD 0.10 with two bidders, about 0.02 with three, and almost zero with four or more.](assets/fig/fig_03_revenue.svg)

### 4.4 Recovering Values from Bids: GPV

Vireo sees bids, not values. The first-price bid is a strictly increasing function of the value; can this be run in reverse, reading the distribution of values off the distribution of bids?

::: {.theorem}
**Theorem 14 (GPV identification; Guerre, Perrigne and Vuong 2000; the quantile form is derived in this chapter)** Under I.1, I.1', and I.2, $n$ is known, bids come from the symmetric equilibrium of Theorem 11, and $b(\cdot)$ is strictly increasing and differentiable; write $F_b$ and $f_b$ for the distribution function and density of a single bid, with $f_b > 0$ in the interior of the support. Then the value corresponding to each equilibrium bid $b$ is

$$v = b + \frac{F_b(b)}{(n-1)\, f_b(b)}$$

Equivalently, for $\tau \in (0, 1)$, the $\tau$ quantile of values satisfies $F^{-1}(\tau) = F_b^{-1}(\tau) + \dfrac{\tau}{n-1}\,\dfrac{d}{d\tau} F_b^{-1}(\tau)$. The value distribution $F$ is therefore identified from the bid distribution.
:::

The proof sketch has three steps. Step one (the first-order condition of Theorem 11): $b'(v) = (n-1)\big(v - b(v)\big)f(v)/F(v)$. Step two (bids strictly increasing): $F_b\big(b(v)\big) = F(v)$, and differentiating gives $f_b\big(b(v)\big)\,b'(v) = f(v)$; substituting into step one solves for $v - b = F_b/\big((n-1)f_b\big)$. Step three (quantiles): $F_b^{-1}(\tau) = b\big(F^{-1}(\tau)\big)$ and $\frac{d}{d\tau}F_b^{-1}(\tau) = 1/f_b\big(F_b^{-1}(\tau)\big)$, so $F_b/f_b$ at the $\tau$ quantile equals $\tau\,\frac{d}{d\tau}F_b^{-1}(\tau)$. Near the extreme quantiles $f_b$ is very small or the derivative is undefined, so the result holds only at interior quantiles.

In words: under symmetric independent private values, risk neutrality, and equilibrium bidding, bid data identify the whole value distribution: each value is the bid plus a shading term that uses only the bid distribution and $n$. This is the empirical version of the statement in Section 3.5: the outcome is a function of types, and when the function is known, the types can be recovered. When the assumptions fail, the inversion recovers something else: with risk-averse bidders the true shading is smaller, so the inversion overstates values; with affiliated values, bids mix in inferences about others' signals; and if an agent such as an autobidding tool bids on advertisers' behalf on the platform, the inversion recovers the agent's "values," not the advertisers'. With the wrong $n$, the inversion is also systematically biased. For estimation, this chapter uses only the empirical quantiles of bids and their differences; the original GPV paper estimates $f_b$ by kernel density estimation, and kernel smoothing and bandwidth choice are in Chapter 14. One testable restriction is that the recovered values must increase with the bids. The V1 numbers: over the quantiles from 5% to 95%, the median relative error of the inversion is 0.83% pooled, 2.5% with two bidders (at most 8.6%), and 0.38% with six; across 1,000 replications of the simulation, the 2.5th and 97.5th percentiles of the pooled median error are 0.66% and 0.98%. The error is large with two bidders partly because $1/(n-1)$ amplifies the noise in the quantile differences.

![The V1 GPV inversion: three panels for $n = 2$, 4, and 6, with the quantile level $\tau$ on the horizontal axis; the dashed line is the true value quantile, the lowest line is the observed quantile of first-price bids, the solid line is the value recovered from this sample, and the shaded band is the 2.5th to 97.5th percentile band across 1,000 replications of the simulation, widest with two bidders.](assets/fig/fig_03_gpv.svg)

### 4.5 Myerson-Satterthwaite: The Impossibility of Efficient Bilateral Trade

Mechanism design also has bad news that draws boundaries. Take the simplest brokering problem: a seller privately knows its cost $c$, a buyer privately knows its value $v$, and efficiency requires trade when $v \geq c$; the platform wants every mutually beneficial trade to happen, without paying out of its own pocket.

::: {.theorem}
**Theorem 15 (Myerson-Satterthwaite; MWG Proposition 23.E.1, p. 895; Myerson and Satterthwaite 1983)** There are one seller and one buyer, both risk neutral; the cost $c$ and the value $v$ are independent, each drawn from an interval with a density that is positive everywhere, and the interiors of the two intervals intersect. Then no mechanism satisfies all three of the following: BIC; ex post efficiency, that is, trade always occurs when $v > c$ and never when $v < c$, and money changes hands only between buyer and seller (budget balance); and interim individual rationality, that is, every type's expected payoff from participating is nonnegative.
:::

The proof sketch has three steps. Step one (Theorem 7 applied to both sides): the buyer's information rent is determined by the probability of trade through the envelope formula; the higher the seller's type, the less it wants to sell, so its envelope formula runs in the opposite direction. Step two (integration by parts): write the buyer's virtual value as $\Psi(v)$, the distribution and density of the cost as $F_c$ and $f_c$, and the seller's virtual cost as $\Psi_c(c) = c + F_c(c)/f_c(c)$; under any BIC and interim individually rational mechanism, the money collected from the buyer minus the money paid to the seller has an expectation no greater than $\mathbb{E}\big[\big(\Psi(v) - \Psi_c(c)\big)\,q(v, c)\big]$, where $q$ is the probability of trade. Step three (ex post efficiency): under the ex post efficient rule $q = 1\{v \geq c\}$ this bound is negative, so the expected deficit is positive, contradicting budget balance; when $v$ and $c$ are both uniform on $[0, 1]$, the bound is $\mathbb{E}\big[(2v - 1 - 2c)\,1\{v \geq c\}\big] = -1/6$.

In words: a platform that brokers bilateral trade under private values cannot make every mutually beneficial trade happen, keep participation voluntary, and do without outside subsidies, all at once. When $v$ and $c$ are both uniform on $[0, 1]$, the second-best mechanism, which maximizes the gains from trade subject to BIC, IR, and budget balance, trades only when $v \geq c + 1/4$ (Myerson and Satterthwaite 1983), and the linear equilibrium of the double auction that sets the price at the midpoint of the two bids implements it exactly. The probability of trade drops from $1/2$ to $9/32 \approx 0.281$, and the gains from trade fall from $1/6$ to $9/64$, which is 84.4% of the efficient level, while every trade whose gains lie between 0 and $1/4$ is given up. With outside subsidies, efficiency is attainable: the pivot mechanism in which the buyer pays $c$ and the seller receives $v$ is efficient, with an expected deficit of $1/6$, equal to the entire gains from trade; by step three (MWG Exercise 23.E.7 (a), pp. 922-923), in the uniform case no efficient, BIC, and IR mechanism has a smaller deficit. With discrete types the theorem can be circumvented (MWG Exercise 23.E.3, p. 922). The result is determined by the overlap of the two distributions and by "no money from outside"; real platforms therefore make do with posted prices, commissions, and thickening the market. In the data, values and costs must first be recovered from bids (the logic of Section 4.4) before the forgone mutually beneficial trades become visible.

### 4.6 Matching and Market Design

Auctions allocate objects through prices, but platforms also face a large class of allocation problems that bypass prices and pair up two sides: labor platforms match workers with firms, and admissions systems match students with schools; the question is how to match so that the matching is stable.

::: {.assumption}
**MT.1 (two-sided one-to-one matching, strict preferences; Gale and Shapley 1962; Roth 1982, §2)** There are two finite, disjoint sets of participants (proposers and receivers, such as workers and firms); each participant $i$ has complete, transitive, strict preferences $\succ_i$ over every member of the other side and over "staying single"; a matching $m$ pairs each participant with at most one member of the other side, and $m(i)$ is $i$'s partner.
:::

MT.1 rules out ties and complementarities: with ties in preferences, a tie-breaking rule must be fixed first, and complementary preferences, such as a couple who want to be in the same city, can make stable matchings cease to exist. A matching $m$ is stable if nobody would rather be single and there is no blocking pair: $i$ and $j$ on opposite sides with $j \succ_i m(i)$ and $i \succ_j m(j)$. An unstable matching unravels on its own, since a blocking pair always has a motive to pair off privately.

::: {.theorem}
**Theorem 16 (a stable matching exists, and DA gives the proposer-optimal stable matching; Gale and Shapley 1962; Roth 1982, Theorems 1 and 2, p. 620)** Under MT.1: (i) a stable matching exists; (ii) deferred acceptance (DA) stops in finitely many steps with a stable result: proposers propose to receivers in order of preference, each receiver tentatively holds the best proposal in hand and rejects the rest, and those rejected move on to their next choice; (iii) the DA outcome is proposer-optimal: every proposer is at least as well off under it as under any other stable matching.
:::

The proof sketch has three steps. Step one (finite sets): each proposer proposes to each receiver at most once, so the algorithm stops in finitely many steps. Step two (the strict preferences of MT.1): if $i$ prefers receiver $j$ to its match, $i$ must have proposed to $j$ earlier and been rejected; $j$ rejected it because it held someone better, and what a receiver holds only improves over time, so $j$ does not prefer $i$, and there is no blocking pair; this gives (i) and (ii). Step three (induction): call $j$ possible for $i$ if some stable matching pairs the two. Suppose that up to some step no proposer has been rejected by a possible receiver, and that $i$ is now rejected by $j$; then in any stable matching that pairs $i$ with $j$, the proposer whom $j$ kept can only be matched to a possible partner it likes less than $j$ (by the induction hypothesis), so that proposer and $j$ form a blocking pair, a contradiction. Hence $j$ is not possible for $i$, so nobody is ever rejected by a possible partner, which is (iii) (Roth 1982, pp. 620-621).

In words: under MT.1 a stable matching always exists, DA is an algorithm that constructs one, and whoever proposes comes out ahead: the proposing side gets its best partner among all stable matchings and the receiving side its worst (Roth and Sotomayor 1990). With complementarities in preferences, a stable matching may not exist. Which side proposes is the designer's decision, and it decides at which end of the set of stable matchings the outcome lands; the set itself is determined by preferences, and the more preferences share a common component, the smaller the set. The data show the matched pairs and (in a centralized clearinghouse) the submitted preference lists, and estimating preferences from them is a separate identification problem. The V0 numbers: with 3 workers and 3 firms, workers 1 and 2 both get their first choice when workers propose, and only their second choice when firms propose; in V4's 200 × 200 market, DA is stable in both directions and only 16 workers have different partners. The average rank of a worker's partner on the worker's own list is 36.61 when workers propose and 39.27 when firms propose, and for firms the two figures are 43.16 and 41.11, in the same order, so the proposer's advantage is real but small.

::: {.theorem}
**Theorem 17 (stability and incentives; Roth 1982, Theorem 3, p. 622; Theorem 5, p. 623; Dubins and Freedman 1981)** Under MT.1: (i) no stable matching mechanism makes truthful reporting of preferences a dominant strategy for all participants; (ii) in the mechanism that always selects the proposer-optimal stable matching (DA run from the proposing side), truthful reporting is a dominant strategy for every proposer.
:::

The proof sketch: (i) uses Roth's counterexample. Construct a 3-by-3 preference profile with exactly two stable matchings, one optimal for the men and one optimal for the women; if the mechanism picks the former, some woman swaps the last two entries on her list, which makes the latter the unique stable matching, and she gains; if the mechanism picks the latter, some man gains by the same maneuver; so whichever the mechanism picks, someone gains by misreporting. The proof of (ii) is longer; see Roth (1982), §5, and Dubins and Freedman (1981). What it says is that a proposer's misreport can only change the set of stable matchings; it cannot get the proposer anyone better than its best stable partner under its true preferences.

In words: stability and universal truth-telling cannot both be had. DA protects only the proposing side, and a receiver may gain by truncating its preference list: if it cuts the list so that the last entry is exactly its partner in the receiver-optimal stable matching, DA gives it that partner. Choosing who proposes is choosing whose truth-telling to protect. The room for gain depends on how large the set of stable matchings is; when the set is small, receivers have almost nothing to gain. When MT.1 fails because of complementarities, a stable matching may not exist, and neither result applies. In the data, if a centralized clearinghouse keeps the submitted lists, systematically short lists on the receiving side may be a sign of truncation. The V4 numbers: none of 2,000 misreports by workers (1,000 truncations and 1,000 swaps of two choices) gained; 216 of them changed the outcome, all to the misreporter's detriment. On the firm side, 14 of 1,000 random truncations gained, and with targeted truncation as described above, the firms that gained were exactly the 16 whose partners differ between the two stable matchings, moving up 25.56 places on their own lists on average, and the 30 profitable reports all came from these 16 firms.

### 4.7 From Incentive Constraints to Testable Predictions

The empirical content of information models comes from how the incentive constraints move with the primitives. The table below writes several of the chapter's comparative statics as predictions that can be taken to data; the last column lists the identification traps that signs alone cannot separate.

| Mechanism | Comparative static | Observable prediction | Identification warning |
|---|---|---|---|
| Screening | Distortion of the low type rises as the information rent gets more expensive | Quality falls at the low end of the menu, exclusion rises | Type distribution and cost are often endogenous together |
| Signaling | Signal cost falls for the low type | $\underline s$ rises, signal-quality gradient flattens | The quality gap between signal holders and non-holders may be only selection |
| Moral hazard | Performance signal becomes more precise | Pay-performance sensitivity or contract cost changes | Monitoring may directly raise productivity |
| Auction competition | Bidder count rises | First-price shading shrinks, revenue rises | Bidder entry is endogenous, needs an entry shifter |
| Persuasion | Common prior or receiver threshold changes | Action rate moves continuously with the prior; only the no-information benchmark jumps at the threshold | Requires observing the information policy and beliefs |

IC, IR, or Bayes plausibility provide only necessary constraints; data that satisfy them do not single out any one mechanism. Identifying a mechanism takes institutional changes under which competing models give different joint responses to the same shock.

## 5 Anchoring Papers

Three studies bring the chapter's language of mechanisms down to real platforms: the auction rules of search advertising, reserve prices set by the theory, and disclosure in a used-car market.

### 5.1 Edelman, Ostrovsky and Schwarz (2007)

::: {.case}
Paper: "Internet Advertising and the Generalized Second-Price Auction: Selling Billions of Dollars Worth of Keywords," *American Economic Review*. It brought the search ad slot auction into the analysis of mechanism design and is the prototype of this chapter's case; Varian (2007) independently gave a related analysis.

Method: what search engines use is not VCG but the GSP auction: advertisers bid per click, higher bids get higher positions, and the winner of each position pays, per click, the bid of the next one down. Despite the "second price" in its name, GSP is not a generalization of VCG: truthful bidding in it is not a dominant strategy. The authors characterize the equilibria of GSP and propose the refinement of locally envy-free equilibrium: no advertiser wants to swap positions and prices with the advertiser just above it.

Results: GSP has multiple equilibria, and the locally envy-free equilibrium with the lowest revenue reproduces VCG exactly in allocation and payments, which explains why GSP, with its simpler rules, can still approximate the properties of VCG. Section 6.1 verifies this on Vireo's two ad slots.

Limitations: the model is static, with a single keyword and exogenous click-through rates; real ad auctions have quality scores, budget constraints, competition across keywords, and dynamic learning. The causal size of the position effect itself has to be identified by randomizing the ranking, which is exactly what Section 6.5 of Chapter 33 does.
:::

### 5.2 Ostrovsky and Schwarz (2023)

::: {.case}
Paper: "Reserve Prices in Internet Advertising Auctions: A Field Experiment," *Journal of Political Economy*. It takes Theorem 13 from the blackboard to an online advertising platform.

Method: a large-scale field experiment that sets reserve prices for online advertising auctions according to optimal auction theory, adapted to the setting of sponsored search. By Theorem 13, the reserve price depends on the value distribution, which must be estimated first, and this is exactly what an inversion of the kind in Section 4.4 is for; for how the paper estimates it, see the original.

Results: consistent with the theory, revenue rose significantly after the new reserve prices were introduced. The platform, the size of the increase, and the exact algorithm for the reserve prices are left to the original paper; this book cites no numbers from it.

Limitations: the gain from a reserve price hinges on how accurately the value distribution is estimated, and it can be eroded by bidders' responses: participation may fall and bids may adjust; if treated and control auctions share the same advertisers, they may interfere with each other (Chapter 26).
:::

### 5.3 Lewis (2011)

::: {.case}
Paper: "Asymmetric Information, Adverse Selection and Online Disclosure: The Case of eBay Motors," *American Economic Review*. Since Akerlof it has been known that asymmetric information causes adverse selection in markets for used goods, which makes the rapid growth of online used-car auctions a puzzle.

Method and data: used-car auctions on eBay Motors. Sellers voluntarily disclose private information on the auction page, and the disclosure forms a precise contract: the seller agrees to deliver the car shown on the page at the transaction price, which protects buyers from adverse selection.

Results: online disclosure is an important determinant of prices; the cost of disclosure affects how much sellers disclose and also the transaction prices.

Limitations: for the size of the effects and the details of identification, see the original paper; disclosure is itself the seller's choice, and sellers of good cars are more willing to disclose, which is exactly the selection problem flagged in the table of Section 4.7.
:::

## 6 Running Case: Settling Vireo's Accounts

The code below runs under R 4.6.1. Only the key statements are shown, with several statements sometimes merged onto one line; every number in the text comes from the output of an actual run, and the result lines in the code that begin with `#` are copied verbatim from what the run printed. Functions such as `vcg_position`, `fpa_bid`, `gpv_quantile`, `optimal_menu`, `len_contract`, and `deferred_acceptance` implement the formulas of Sections 3 and 4; their names say what they do.

### 6.1 V0: Three Accounts for Two Ad Slots

Efficiency requires that the advertiser with value 10 get the top slot, the one with value 6 the second slot, and the one with value 4 nothing. VCG payments are computed as externalities: without advertiser 1, advertiser 2 moves up to the top slot and advertiser 3 fills the second slot, so the others get $6 \times 100 + 4 \times 40 = 760$ dollars in total, while with it present the others get only $6 \times 40 = 240$ dollars, so advertiser 1 pays USD 520; by the same logic advertiser 2 pays $1160 - 1000 = 160$ dollars.

```r
v <- c(10, 6, 4); ctr <- c(100, 40, 0)                        # USD per click (sorted); clicks a month; 0 = no slot
pm_naive <- v[1] * ctr[1] + v[2] * ctr[2]                      # pay-your-bid, if bids were the true values
pay1 <- (v[2] * ctr[1] + v[3] * ctr[2]) - (v[2] * ctr[2] + v[3] * ctr[3])   # others' best without 1 - others' total with 1
pay2 <- (v[1] * ctr[1] + v[3] * ctr[2]) - (v[1] * ctr[1] + v[3] * ctr[3])
gsp_truth <- v[2] * ctr[1] + v[3] * ctr[2]                     # GSP, truthful bids: each slot pays the next bid
p1_low <- (v[2] * (ctr[1] - ctr[2]) + v[3] * ctr[2]) / ctr[1]   # advertiser 2 indifferent between the two slots
c(pay_your_bid = pm_naive, vcg_pay1 = pay1, vcg_pay2 = pay2, per_click_1 = pay1 / ctr[1], per_click_2 = pay2 / ctr[2],
  vcg_revenue = pay1 + pay2, gsp_truthful = gsp_truth, gsp_lowest_lef = p1_low * ctr[1] + v[3] * ctr[2])
#   pay_your_bid       vcg_pay1       vcg_pay2    per_click_1    per_click_2
#         1240.0          520.0          160.0            5.2            4.0
#    vcg_revenue   gsp_truthful gsp_lowest_lef
#          680.0          760.0          680.0
```

Under VCG the two winners net $10 \times 100 - 520 = 480$ and $6 \times 40 - 160 = 80$ dollars. Under GSP, if all three bid truthfully, the top slot pays the second-ranked bid of USD 6 per click and the second slot the third-ranked bid of USD 4, for revenue of USD 760, even more than VCG. Truthful bidding is not dominant in GSP, but in this example it happens to be a locally envy-free equilibrium: advertiser 2, moving to the top slot and paying USD 6, would get only 0, less than its current USD 80; advertiser 3, moving to the second slot and paying USD 4, would also get only 0; advertiser 1, dropping to the second slot, would get $40 \times (10 - 4) = 240$ dollars, less than its current USD 400. The locally envy-free equilibrium with the lowest revenue makes advertiser 2 exactly indifferent between the two slots, with a top-slot price of USD 5.2 per click, and reproduces VCG's USD 680. So USD 680 is the lower bound on the revenue of these equilibria, and USD 760 is the revenue of just one of them.

![V0's two ad slots: the left panel shows each advertiser's per-click value (bars), the VCG per-click prices (USD 5.2 and 4, payments of USD 520 and 160), and the prices under truthful GSP bidding (USD 6 and 4); the right panel shows monthly revenue under three rules: USD 1240 for pay-your-bid assuming truthful bids (not an equilibrium), USD 760 for GSP with truthful bids, and USD 680 for the lowest-revenue locally envy-free GSP equilibrium and for VCG.](assets/fig/fig_03_position.svg)

### 6.2 V0: The Remaining Small Accounts

The remaining hand-computed numbers of Sections 3 and 4 come from the lines below: revenue equivalence and the reserve price with uniform values, Akerlof, the two effort levels, persuasion at the recommendation slot, bilateral trade, the number of cheap-talk intervals, and the 3 × 3 DA.

```r
n <- c(2, 3); r <- 0.5                                         # values uniform on [0, 1]
rbind(second_price = (n - 1) / (n + 1), first_price = (n - 1) / n * n / (n + 1))
c(reserve_rev_n2 = 2 * r * (1 - r) * r + (1 - r)^2 * (r + (1 - r) / 3), unsold = r^2)
p <- seq(0.1, 1, by = 0.1); all(1.5 * p / 2 < p)               # Akerlof: buyers' value of the sellers who stay
gap <- 5 / (0.8 - 0.2); c(w_H = gap, rent = 0.8 * gap - 5, principal = 0.8 * (100 - gap), first_best = 0.8 * 100 - 5)
x <- persuasion_x(0.3, 0.5); c(x = x, P_buy = 0.3 + 0.7 * x, full_disclosure = 0.3)
wd <- 1 / 4                                                    # bilateral trade: trade iff v >= c + 1/4
c(p_trade = (1 - wd)^2 / 2, gains = ((1 - wd^3) / 3 - wd^2 * (1 - wd)) / 2, gains_first_best = 1 / 6)
b <- c(0.01, 0.05, 0.1, 0.2, 0.25); setNames(cheap_talk_N(b), b)  # cheap talk: most intervals N for bias b
wp <- list(c(1, 2, 3), c(2, 1, 3), c(1, 2, 3)); fp <- list(c(2, 1, 3), c(1, 2, 3), c(1, 2, 3))
rbind(workers_propose = deferred_acceptance(wp, fp), firms_propose = invert_match(deferred_acceptance(fp, wp), 3))
#                   [,1] [,2]
# second_price 0.3333333  0.5
# first_price  0.3333333  0.5
# reserve_rev_n2         unsold
#      0.4166667      0.2500000
# [1] TRUE
#        w_H       rent  principal first_best
#   8.333333   1.666667  73.333333  75.000000
#               x           P_buy full_disclosure
#       0.4285714       0.6000000       0.3000000
#          p_trade            gains gains_first_best
#        0.2812500        0.1406250        0.1666667
# 0.01 0.05  0.1  0.2 0.25
#    7    3    2    2    1
#                 [,1] [,2] [,3]
# workers_propose    1    2    3
# firms_propose      2    1    3
```

### 6.3 V1: Ten Thousand Ad Exchange Auctions

V1 has 10,000 auctions, 2,000 for each number of bidders $n$ from 2 to 6; values are i.i.d. draws from a normal distribution truncated to $[0, 3]$ dollars per click (mean parameter 1.3, standard deviation parameter 0.6), whose density is log-concave, hence regular, and positive everywhere (I.1). First-price equilibrium bids are computed by numerical integration of the formula in Theorem 11 (with the reserve price, the lower limit becomes $r^\ast$), second-price bids are truthful, and each format is also run with the reserve price $r^\ast$; the four formats are run on the same values, the revenue difference is the mean of the auction-by-auction differences, and its standard error is the standard deviation of the differences divided by the square root of the number of auctions.

```r
P1 <- list(F = list(family = "truncnorm", mean = 1.3, sd = 0.6, lower = 0, upper = 3),
           n_set = 2:6, n_auc_per_n = 2000L, taus = seq(0.05, 0.95, by = 0.01), h = 0.02)   # h: half-width of the quantile difference
D1 <- value_dist(P1$F); rs <- optimal_reserve(D1)              # Psi(r*) = 0
set.seed(3101); V1 <- gen_V1(P1)                               # values; exact equilibrium bids, with and without r*
auc <- V1[, {
  o <- order(v_true, decreasing = TRUE); v1 <- v_true[o[1]]; v2 <- v_true[o[2]]
  list(n = n[1], rev_fp = max(b_fp), rev_sp = v2,               # first price: top bid; second price: second value
       rev_fp_res = if (v1 < rs) 0 else max(b_fp_res, na.rm = TRUE), rev_sp_res = if (v1 < rs) 0 else max(v2, rs))
}, by = auction]
d <- auc$rev_fp - auc$rev_sp; g <- auc$rev_sp_res - auc$rev_sp  # paired: the same value draws
round(c(r_star = rs, fp = mean(auc$rev_fp), sp = mean(auc$rev_sp), diff = mean(d), se_diff = sd(d) / sqrt(length(d)),
        fp_res = mean(auc$rev_fp_res), sp_res = mean(auc$rev_sp_res), gain = mean(g), se_gain = sd(g) / sqrt(length(g))), 4)
gpv <- rbindlist(lapply(P1$n_set, function(k) {
  q <- gpv_quantile(V1[n == k]$b_fp, k, P1$taus, h = P1$h)       # v_hat = b + tau b'(tau) / (n - 1)
  q[, .(n = k, relerr = abs(v_hat - D1$quantile(tau)) / D1$quantile(tau))]
}))
gpv[, .(median_relerr = round(median(relerr), 4), max_relerr = round(max(relerr), 4)), by = n]
#  r_star      fp      sp    diff se_diff  fp_res  sp_res    gain se_gain
#  1.0615  1.4173  1.4222 -0.0049  0.0031  1.4420  1.4449  0.0228  0.0019
#        n median_relerr max_relerr
#    <int>         <num>      <num>
# 1:     2        0.0247     0.0862
# 2:     3        0.0095     0.0636
# 3:     4        0.0074     0.0392
# 4:     5        0.0069     0.0368
# 5:     6        0.0038     0.0253
```

Sampling uncertainty is reported by simulation: repeating the same design 1,000 times (seed 3102), the gain from the reserve price is positive every time; in terms of the 2.5th to 97.5th percentiles, the first-price minus second-price difference is $[-0.0060, 0.0063]$ and the pooled GPV median error is 0.66% to 0.98%.

### 6.4 V2: The SaaS Usage Menu

The type $\theta$ is the value per thousand API calls, on a grid of 2,001 points between USD 2 and 10; the regular distribution is a truncated normal $\mathcal{N}(6, 2^2)$, and the bimodal distribution is the truncation of $0.7\,\mathcal{N}(4, 0.7^2) + 0.3\,\mathcal{N}(8, 0.7^2)$. `optimal_menu` computes the menu by Theorem 3: it runs a monotone regression of $\Psi$ in quantile space to get $\bar\Psi$, sets $q = \max(0, \bar\Psi)$, and then computes the information rents and prices from the envelope formula.

```r
th <- seq(2, 10, length.out = 2001)                            # types: USD per thousand API calls
dens <- list(regular = list(family = "truncnorm", mean = 6, sd = 2, lower = 2, upper = 10),
             bimodal = list(family = "truncnorm_mix", weights = c(0.7, 0.3), means = c(4, 8), sds = c(0.7, 0.7),
                            lower = 2, upper = 10))
t(sapply(dens, function(par) {
  m <- optimal_menu(th, type_density(par, th), cost2 = 1)      # q = max(0, ironed Psi); U = int q; t = theta q - U
  m0 <- optimal_menu(th, type_density(par, th), cost2 = 1, iron = FALSE)
  s <- menu_summary(m, 1)
  round(c(excluded_below = s$theta_exclusion, q_top = s$q_top, profit = s$profit, rent = s$rent,
          surplus_efficient = s$surplus_efficient, bunch_lo = s$bunch_lo, bunch_hi = s$bunch_hi, bunch_q = s$bunch_q,
          ic_gain_without_ironing = menu_ic_gain(m0)), 4)
}))
unlist(screening_two_type(theta_L = 4, theta_H = 7, mu = 0.4, cost2 = 1)[c("q_L", "q_H", "t_L", "t_H", "rent_H",
                                                                          "profit", "q_L_fb", "profit_fb")])
#         excluded_below q_top  profit   rent surplus_efficient bunch_lo bunch_hi
# regular          4.632    10 10.2911 5.5385           19.5475       NA       NA
# bimodal          3.356    10  7.9376 4.4662           15.4422     3.78     6.94
#         bunch_q ic_gain_without_ironing
# regular  0.0000                  0.0000
# bimodal  1.8398                  5.2912
#       q_L       q_H       t_L       t_H    rent_H    profit    q_L_fb profit_fb
#       2.0       7.0       8.0      43.0       6.0      11.0       4.0      14.6
```

![The V2 optimal menu: the left panel is the regular distribution and the right panel the bimodal one; the solid line is the optimal menu (by the ironed virtual value), the diagonal dashed line is the efficient usage $q = \theta$, and the dotted line is the menu without ironing (which is not incentive compatible). Under the bimodal distribution, the types between USD 3.78 and 6.94 are bunched on a single tier at $q = 1.84$.](assets/fig/fig_03_menu.svg)

### 6.5 V3: The Reviewers' Linear Contract

Each of Vireo's ad reviewers creates $x = e + \varepsilon$ thousand dollars of value for the platform a month, with $\sigma^2 = 4$; three other contractible measures are available: the rate at which users report approved ads, $y = -0.5e + \eta$, where $\eta$ has variance 1 and is independent of $\varepsilon$; the output of another team in the same group, unrelated to effort, whose noise also has variance 4 and a correlation of 0.6 with $\varepsilon$; and a measure unrelated to both effort and noise.

```r
A <- 0.5; c2 <- 1; s2 <- 4                                     # CARA per thousand USD; c(e) = e^2 / 2; Var(eps) = 4
two <- function(k, var, rho) {                                 # x = e + eps and a second signal y = k e + eta
  cv <- rho * sqrt(s2 * var); len_contract(A, c2, matrix(c(s2, cv, cv, var), 2), c(1, k))
}
K <- list(output_only = len_contract(A, c2, matrix(s2), 1), report_rate = two(-0.5, 1, 0),
          peer_output = two(0, 4, 0.6), pure_noise = two(0, 1, 0))
t(sapply(K, function(k) round(c(sigma2_tilde = k$sigma_s2, beta = k$beta, weight_x = k$b[1],
                                weight_y = if (length(k$b) > 1) k$b[2] else 0, effort = k$effort,
                                risk_premium = k$risk_premium, profit_usd = 1000 * k$profit), 4)))
#             sigma2_tilde   beta weight_x weight_y effort risk_premium
# output_only         4.00 0.3333   0.3333   0.0000 0.3333       0.1111
# report_rate         2.00 0.5000   0.2500  -0.5000 0.5000       0.1250
# peer_output         2.56 0.4386   0.4386  -0.2632 0.4386       0.1231
# pure_noise          4.00 0.3333   0.3333   0.0000 0.3333       0.1111
#             profit_usd
# output_only   166.6667
# report_rate   250.0000
# peer_output   219.2982
# pure_noise    166.6667
```

Both the report rate and the peer benchmark make $\tilde\sigma^2$ smaller, and the pure-noise measure changes nothing. A small simulation (seed 3103) compares the two contracts over 500 replications of 1,000 reviewer-months each: adding the report rate raises profit by USD 84.07 on average, with a 2.5th to 97.5th percentile range of USD 51.00 to 118.34, positive in all 500 replications, consistent with the closed-form USD 83.33.

### 6.6 V4: A 200 × 200 Talent Market

There are 200 workers and 200 firms; a worker's utility from firm $j$ is a common component of the firm plus an idiosyncratic component, firms' utilities over workers are built the same way, the two components have equal variances, and preference lists are complete and strict (MT.1).

```r
set.seed(3104)
g4 <- gen_V4(list(n_workers = 200L, n_firms = 200L, sd_common = 1, sd_idio = 1))   # common + idiosyncratic utility
pp <- g4$prefs$prop; pr <- g4$prefs$recv                       # workers' lists of firms; firms' lists of workers
mP <- deferred_acceptance(pp, pr)                              # workers propose: worker -> firm
mR <- invert_match(deferred_acceptance(pr, pp), 200)           # firms propose, written worker -> firm
fP <- invert_match(mP, 200); fR <- invert_match(mR, 200)       # firm -> worker
c(stable_P = is_stable(mP, pp, pr), stable_R = is_stable(mR, pp, pr)); sum(mP != mR)   # workers with different partners
round(c(worker_rank_P = mean(partner_rank(mP, pp)), worker_rank_R = mean(partner_rank(mR, pp)),
        firm_rank_P = mean(partner_rank(fP, pr)), firm_rank_R = mean(partner_rank(fR, pr))), 3)
gain <- sapply(1:200, function(j) {                            # firm j truncates its list just below its partner in mR
  prj <- pr; prj[[j]] <- truncate_pref(pr[[j]], match(fR[j], pr[[j]]))
  f2 <- invert_match(deferred_acceptance(pp, prj), 200)
  partner_rank(fP, pr)[j] - partner_rank(f2, pr)[j]            # places gained on the true list
})
c(firms_gaining = sum(gain > 0), mean_gain = mean(gain[gain > 0]))
# stable_P stable_R
#     TRUE     TRUE
# [1] 16
# worker_rank_P worker_rank_R   firm_rank_P   firm_rank_R
#        36.610        39.265        43.155        41.110
# firms_gaining     mean_gain
#       16.0000       25.5625
```

The misreporting checks are run on the same preferences (seed 3105); the results are in Section 4.6.

## 7 Failure Modes and Robustness

The first trap is mistaking VCG's theoretical halo for an engineering solution. Its revenue can fall as bidders are added, and it is easily manipulated by shill bidding and collusion; its payments require solving for the counterfactual allocation in each participant's absence, so advertisers struggle to understand why they pay what they pay. Search engines chose GSP, trading dominant strategies for transparency.

The second is private values and independence. In many platform auctions, values have a common-value component: the commercial value of a keyword is roughly the same for all advertisers, and each has only a noisy estimate of it. The winner tends to be the one who overestimates most, so winning is itself bad news, and a rational bidder must shade more in advance; this is the winner's curse. Revenue equivalence then fails: in the English auction, exits are public and reveal information step by step, and revenue is often higher. With correlated values the seller can even exploit the correlation to extract the entire surplus (Crémer and McLean 1988), and standard auctions are no longer equivalent. In these cases the GPV inversion does not recover values.

The third is multidimensional types. The clean solution of Theorem 3 holds only for one-dimensional types; real customers are heterogeneous at once in usage, time of day, and how much they care about reliability, and multidimensional screening generally has no simple characterization: the optimal menu can be extremely complex and may even involve randomization.

The fourth is the sensitivity of moral hazard to the specification. The first-order approach is reliable only under conditions such as K.1 and K.2; performance measures can often be gamed, and strong incentives induce behavior that makes the measure look good rather than behavior that gets the job done; with multiple tasks, strong incentives on the measurable tasks crowd out tasks that cannot be measured but matter just as much (Holmström and Milgrom 1991). In designing incentives, measurability and gameability deserve thought before mathematical optimality.

Last is commitment. The revelation principle, information design, and the optimal auction all assume that the designer can commit to the rules; if the platform would exploit the reported information after the fact, participants will anticipate this, and all the conclusions above have to be rewritten. The sharpest theorems of mechanism design all rest on private values, independence, one-dimensionality, and commitment; robust design lies in recognizing which of these is violated and in which direction the violation pushes the conclusion, and often in preferring a second-best but transparent, manipulation-resistant mechanism.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Mas-Colell, Whinston and Green (1995, Microeconomic Theory), Chapter 13 (adverse selection, signaling, screening), Chapter 14 (principal-agent and moral hazard), and Chapter 23 (mechanism design and auctions). The theoretical skeleton of this chapter comes from these three chapters.
- Krishna (2010, Auction Theory). The standard textbook on auction theory, more detailed than MWG.
- Myerson (1981). The complete framework of the optimal auction and revenue equivalence; read it for the virtual value and the envelope argument.

Further reading:

- Akerlof (1970), Spence (1973), and Cho and Kreps (1987). The original papers on adverse selection, signaling, and the intuitive criterion.
- Mussa and Rosen (1978) and Maskin and Riley (1984). The original papers on screening with continuous types, the source of Theorem 3.
- Holmström (1979), Grossman and Hart (1983), Rogerson (1985), and Jewitt (1988). The informativeness principle and the conditions under which the first-order approach is valid; Holmström and Milgrom (1987, 1991) cover linear contracts and multitasking.
- Myerson (1979), Green and Laffont (1977), and Vickrey (1961). The revelation principle, the characterization of Groves mechanisms, and the origin of the second-price auction.
- Riley and Samuelson (1981), Milgrom and Weber (1982), and Myerson and Satterthwaite (1983). Reserve prices, affiliated values, and the impossibility of efficient bilateral trade.
- Guerre, Perrigne and Vuong (2000) and Athey and Haile (2002). The founding works on identifying value distributions from bids, the starting point of structural auction econometrics.
- Edelman, Ostrovsky and Schwarz (2007), Varian (2007), and Ostrovsky and Schwarz (2023). The theory of search ad auctions and a field experiment on them, required reading for information systems research; Einav, Farronato, Levin and Sundaresan (2018) discuss why online markets moved from auctions to posted prices.
- Kamenica and Gentzkow (2011), Bergemann and Morris (2019), and Crawford and Sobel (1982). Bayesian persuasion, the unified survey of information design, and the original paper on cheap talk.
- Gale and Shapley (1962), Dubins and Freedman (1981), Roth (1982, 1984), and Roth and Sotomayor (1990). Stable matching, incentives, and the empirics of the medical residency match.
:::

::: {.apa-refs}
- Akerlof, G. A. (1970). The market for "lemons": Quality uncertainty and the market mechanism. *The Quarterly Journal of Economics, 84*(3), 488-500. https://doi.org/10.2307/1879431
- Athey, S., & Haile, P. A. (2002). Identification of standard auction models. *Econometrica, 70*(6), 2107-2140. https://doi.org/10.1111/j.1468-0262.2002.00435.x
- Bergemann, D., & Morris, S. (2019). Information design: A unified perspective. *Journal of Economic Literature, 57*(1), 44-95. https://doi.org/10.1257/jel.20181489
- Cho, I.-K., & Kreps, D. M. (1987). Signaling games and stable equilibria. *The Quarterly Journal of Economics, 102*(2), 179-221. https://doi.org/10.2307/1885060
- Crawford, V. P., & Sobel, J. (1982). Strategic information transmission. *Econometrica, 50*(6), 1431-1451. https://doi.org/10.2307/1913390
- Crémer, J., & McLean, R. P. (1988). Full extraction of the surplus in Bayesian and dominant strategy auctions. *Econometrica, 56*(6), 1247-1257. https://doi.org/10.2307/1913096
- Dubins, L. E., & Freedman, D. A. (1981). Machiavelli and the Gale-Shapley algorithm. *The American Mathematical Monthly, 88*(7), 485-494. https://doi.org/10.1080/00029890.1981.11995301
- Edelman, B., Ostrovsky, M., & Schwarz, M. (2007). Internet advertising and the generalized second-price auction: Selling billions of dollars worth of keywords. *American Economic Review, 97*(1), 242-259. https://doi.org/10.1257/aer.97.1.242
- Einav, L., Farronato, C., Levin, J., & Sundaresan, N. (2018). Auctions versus posted prices in online markets. *Journal of Political Economy, 126*(1), 178-215. https://doi.org/10.1086/695529
- Gale, D., & Shapley, L. S. (1962). College admissions and the stability of marriage. *The American Mathematical Monthly, 69*(1), 9-15. https://doi.org/10.1080/00029890.1962.11989827
- Green, J., & Laffont, J.-J. (1977). Characterization of satisfactory mechanisms for the revelation of preferences for public goods. *Econometrica, 45*(2), 427-438. https://doi.org/10.2307/1911219
- Grossman, S. J., & Hart, O. D. (1983). An analysis of the principal-agent problem. *Econometrica, 51*(1), 7-45. https://doi.org/10.2307/1912246
- Guerre, E., Perrigne, I., & Vuong, Q. (2000). Optimal nonparametric estimation of first-price auctions. *Econometrica, 68*(3), 525-574. https://doi.org/10.1111/1468-0262.00123
- Holmström, B. (1979). Moral hazard and observability. *The Bell Journal of Economics, 10*(1), 74-91. https://doi.org/10.2307/3003320
- Holmström, B., & Milgrom, P. (1987). Aggregation and linearity in the provision of intertemporal incentives. *Econometrica, 55*(2), 303-328. https://doi.org/10.2307/1913238
- Holmström, B., & Milgrom, P. (1991). Multitask principal-agent analyses: Incentive contracts, asset ownership, and job design. *The Journal of Law, Economics, and Organization, 7*(special issue), 24-52. https://doi.org/10.1093/jleo/7.special_issue.24
- Jewitt, I. (1988). Justifying the first-order approach to principal-agent problems. *Econometrica, 56*(5), 1177-1190. https://doi.org/10.2307/1911363
- Kamenica, E., & Gentzkow, M. (2011). Bayesian persuasion. *American Economic Review, 101*(6), 2590-2615. https://doi.org/10.1257/aer.101.6.2590
- Krishna, V. (2010). *Auction theory* (2nd ed.). Academic Press.
- Lewis, G. (2011). Asymmetric information, adverse selection and online disclosure: The case of eBay Motors. *American Economic Review, 101*(4), 1535-1546. https://doi.org/10.1257/aer.101.4.1535
- Mas-Colell, A., Whinston, M. D., & Green, J. R. (1995). *Microeconomic theory*. Oxford University Press.
- Maskin, E., & Riley, J. (1984). Monopoly with incomplete information. *The RAND Journal of Economics, 15*(2), 171-196. https://doi.org/10.2307/2555674
- Milgrom, P. R., & Weber, R. J. (1982). A theory of auctions and competitive bidding. *Econometrica, 50*(5), 1089-1122. https://doi.org/10.2307/1911865
- Mussa, M., & Rosen, S. (1978). Monopoly and product quality. *Journal of Economic Theory, 18*(2), 301-317. https://doi.org/10.1016/0022-0531(78)90085-6
- Myerson, R. B. (1979). Incentive compatibility and the bargaining problem. *Econometrica, 47*(1), 61-73. https://doi.org/10.2307/1912346
- Myerson, R. B. (1981). Optimal auction design. *Mathematics of Operations Research, 6*(1), 58-73. https://doi.org/10.1287/moor.6.1.58
- Myerson, R. B., & Satterthwaite, M. A. (1983). Efficient mechanisms for bilateral trading. *Journal of Economic Theory, 29*(2), 265-281. https://doi.org/10.1016/0022-0531(83)90048-0
- Ostrovsky, M., & Schwarz, M. (2023). Reserve prices in internet advertising auctions: A field experiment. *Journal of Political Economy, 131*(12), 3352-3376. https://doi.org/10.1086/725702
- Riley, J. G., & Samuelson, W. F. (1981). Optimal auctions. *American Economic Review, 71*(3), 381-392.
- Rogerson, W. P. (1985). The first-order approach to principal-agent problems. *Econometrica, 53*(6), 1357-1367. https://doi.org/10.2307/1913212
- Roth, A. E. (1982). The economics of matching: Stability and incentives. *Mathematics of Operations Research, 7*(4), 617-628. https://doi.org/10.1287/moor.7.4.617
- Roth, A. E. (1984). The evolution of the labor market for medical interns and residents: A case study in game theory. *Journal of Political Economy, 92*(6), 991-1016. https://doi.org/10.1086/261272
- Roth, A. E., & Sotomayor, M. A. O. (1990). *Two-sided matching: A study in game-theoretic modeling and analysis* (Econometric Society Monographs No. 18). Cambridge University Press.
- Spence, M. (1973). Job market signaling. *The Quarterly Journal of Economics, 87*(3), 355-374. https://doi.org/10.2307/1882010
- Varian, H. R. (2007). Position auctions. *International Journal of Industrial Organization, 25*(6), 1163-1178. https://doi.org/10.1016/j.ijindorg.2006.10.002
- Vickrey, W. (1961). Counterspeculation, auctions, and competitive sealed tenders. *The Journal of Finance, 16*(1), 8-37. https://doi.org/10.1111/j.1540-6261.1961.tb02789.x
:::
