---
title: "IV: LATE, MTE & Control Functions"
subtitle: "From Encouragement to the Marginal Treatment Effect"
seriesline: "Foundations of Information Systems Economics · Chapter 12"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 12 · IV: LATE, MTE & Control Functions"
---

## Introduction

Meridian's data hand us a pleasant answer: merchants who adopt the new tool do far better than those who do not, and the gap stays large even after we throw in a long list of observable characteristics. What is genuinely unsettling is that the most capable merchants are also the most eager to adopt. If managerial ability goes unrecorded, no list of controls, however long, can coax it out of the error term. OLS does its arithmetic diligently, yet it very likely mistakes "who adopts" for "what adoption does."

Instrumental variables switch the question. Rather than compare adopters with non-adopters directly, we look for an exogenous nudge: something that changes the probability of adoption but that should not change the outcome through any other channel. Randomized encouragement supplies exactly this kind of variation. The change in the outcome per unit change in the encouragement, divided by the change in the adoption rate per unit change in the encouragement, gives us not some magically cleansed population average effect, but the effect on those whom this particular encouragement pushes toward adoption.

This distinction is the crux of the whole chapter. IV trades endogeneity away for less but cleaner variation, usually at the price of a larger standard error; and the LATE it identifies belongs to the compliers, not automatically to every merchant. Who gets moved by the instrument is who enters the answer. MTE will go further and show that different instruments, OLS, ATE, ATT, and LATE are all just weighting different marginal populations. So every time you see an IV coefficient, you should ask twice: does this instrument really act only through the treatment? And even if the answer is yes, whose effect exactly is it reporting?

## 1 An Absurdly Large Number

::: {.case}
Meridian is a fictional B2B SaaS platform that sells an analytics dashboard to a large base of client firms. We tell the story with simulated data, and the payoff is that the data-generating process is fully known, so every estimator can be reconciled against the truth, a teaching luxury that real data can never offer. The outcome $Y$ is a client firm's monthly engagement (a 0-to-100 activity index). The treatment $D$ is whether the firm has adopted the analytics dashboard. Adoption is voluntary: the sharper, faster-growing firms are the ones that adopt on their own, and "sharpness" is precisely the thing the platform cannot observe. Management's question is direct: how much more engaged does adopting the dashboard actually make a client?

Start with the most natural calculation: subtract the average engagement of firms that did not adopt from the average engagement of firms that did. That difference is 12.104 (SE 0.102). Laid out on a growth-review slide, it says that clients who adopt the dashboard are 12 points more engaged, a number juicy enough for the product team to take credit for.

The platform's data scientist did not take it at face value. She worried that adopters were simply a more engaged bunch to begin with, so she threw every observable firm characteristic she had into the regression, the most important being firm size $W$ (larger firms are both sharper and keener to try new things). After the adjustment, the coefficient on $D$ fell from 12.104 to 5.285 (SE 0.063). A big drop, but still high.
:::

Neither number is credible, and Chapter 9 has already told us why. The naive comparison, 12.104, equals the true treatment effect plus a selection bias: firms that adopt the dashboard would have been more engaged even without it. Controlling for firm size shaves off the part of that bias proxied by size, which is why the number falls from 12.104 to 5.285. The trouble is that the "sharpness" driving adoption is unobserved, and firm size is only an imperfect shadow of it. Once size is held fixed, the part of sharpness that size fails to proxy stays in the residual and keeps propping up adopters' engagement. This is exactly the ineradicable remainder from Chapter 9: selection lands on the unobservables, and no amount of clever covariate adjustment leaves no tail behind.

So what is the truth? Because this is a simulation, we know it. The average effect of adopting the dashboard on firms that did adopt (the ATT) is 3.056, and the average effect on a randomly drawn firm from the whole population (the ATE) is 2.018. The naive comparison 12.104 minus the ATT 3.056 leaves 9.048, all of it selection bias, nearly three-quarters of the reported number. The adjusted 5.285 still exceeds the ATT by 2.230, and that 2.230 is the residual bias that size cannot proxy away, attached to the unobserved sharpness.

::: {.intuition}
Meridian's predicament fits in one sentence: adopters are highly engaged half because of the dashboard and half because they were sharp to begin with. Covariate adjustment is like taking a not-quite-clean rag to wipe off this illusion. It removes what the rag can reach (the sharpness that size proxies) and cannot remove what it cannot reach (the sharpness that size fails to proxy). The lesson of Chapter 9 is that as long as selection still clings to an unobserved factor, this rag will never wipe the surface clean. To cure it, we need a completely different idea: not to wipe, but to go around.
:::

Going around is something Meridian happened to do once. The platform's customer success team ran a round of randomized encouragement: from the client base they drew a random batch and gave them an extra push, an onboarding invitation plus a set of getting-started tutorials; the rest of the clients carried on as usual. Call this randomly assigned encouragement $Z$, with $Z=1$ for firms that got the push and $Z=0$ for those that did not. The encouragement was decided by a coin flip, so it has nothing to do with whether a firm is sharp or not, a fact guaranteed by design and requiring no modeling whatsoever.

Encouragement is not the same as adoption. Only some of the encouraged firms actually went on to adopt the dashboard (this batch, which moves only when given a push, are the compliers); some firms would have adopted on their own even without encouragement (always-takers, the die-hard adopters); and some firms would not adopt no matter how hard they were pushed (never-takers, the immovable). This imperfect compliance is exactly where IV earns its keep, and it is also the root of why IV can see only the compliers. What the rest of this chapter sets out to make clear is how this random push pulls the causal effect of adoption out of the contamination of selection, and whose effect it is that gets pulled out.

## 2 The Economic Model and the Estimand

Before doing anything, get two things straight: which causal quantity we actually want, and what notation writes down endogeneity and the instrument precisely.

### 2.1 Potential Outcomes Plus Potential Treatments

The potential outcomes from Chapter 9 carry over unchanged. Each firm has two potential engagement levels, $Y_i(1)$ being the engagement it would realize if it adopted the dashboard and $Y_i(0)$ the engagement it would realize if it did not, with the individual effect $\tau_i = Y_i(1) - Y_i(0)$ varying from firm to firm (heterogeneous effects). Observed engagement is given by the switching equation, $Y_i = Y_i(0) + \tau_i D_i$. So far this is all old notation.

What is new in IV is that we give the treatment itself a set of potential outcomes too. Introduce potential treatments $D_i(z)$: whether firm $i$ would adopt when the encouragement takes value $z$. Here $D_i(1)$ is its adoption status when encouraged and $D_i(0)$ its adoption status when not encouraged. Just as with the outcome, we only ever observe one of them, because a firm is either encouraged or not, and the realized adoption status is

$$D_i = D_i(0) + \big(D_i(1) - D_i(0)\big)\, Z_i.$$

Here $D_i(1) - D_i(0)$ is the individual causal effect of encouragement on this firm's adoption. It takes three values: equal to 1 means this firm adopts only when pushed (complier), equal to 0 means adoption is unaffected by encouragement (always-taker or never-taker), and equal to $-1$ means it perversely does the opposite of what the encouragement suggests (defier). Because the pair $\big(D_i(0), D_i(1)\big)$ is never seen in full, which compliance type a firm belongs to is not identified at the individual level, and this is the starting point for understanding all of IV's subtleties.

### 2.2 Endogeneity: Why OLS Is Neither ATE nor ATT

Write the outcome equation in the decomposition form from Chapter 9, $Y_i = \alpha + \tau_i D_i + \eta_i$, where $\eta_i = Y_i(0) - \alpha$ collects everything other than adoption that determines baseline engagement (that is, the deviation of $Y_i(0)$ from its mean). The precise meaning of endogeneity is $\mathbb{E}[\eta_i \mid D_i] \neq 0$: adopters and non-adopters already differ in ways unrelated to adoption. At Meridian, sharp firms both adopt more readily and are more engaged to begin with (higher $\eta$), so adoption is positively correlated with the disturbance, and this is the selection problem.

Once there is endogeneity, the OLS estimate $\hat\beta$ converges to neither the ATE nor the ATT. That skeleton decomposition from Chapter 9 still holds here: the naive comparison equals the ATT plus selection bias, $\mathbb{E}[Y\mid D=1] - \mathbb{E}[Y\mid D=0] = \text{ATT} + \big(\mathbb{E}[Y(0)\mid D=1] - \mathbb{E}[Y(0)\mid D=0]\big)$. When the selection-bias term is nonzero, OLS departs from the ATT; and even if the selection bias were miraculously zero, heterogeneous effects would still make the ATT depart from the ATE. Meridian's 12.104 is this decomposition made flesh: the ATT of 3.056 plus selection bias of 9.048, exact to the decimal. What OLS converges to is a number inflated by selection, corresponding to none of the causal quantities we care about.

::: {.intuition}
This endogeneity problem has an ancient incarnation, and it is the one that forced IV into being. Economics first met it in demand estimation (Working 1927). You have a cloud of price-and-quantity data points and want to estimate the slope of the demand curve, but price never falls from the sky; it is set jointly by supply and demand, and $P$ carries both demand shocks and supply shocks. Regress $Q$ on $P$ directly and the "statistical demand curve" you estimate is not the demand curve in any economic sense, because price and the demand disturbance are tangled together. To identify the demand slope, you need something that moves supply only and does not touch demand (weather affecting the catch, say), and you use it to manufacture exogenous variation in price. That supply-shifting-only thing is the earliest instrument. Meridian's randomized encouragement pushes adoption only and does not touch the engagement disturbance; the logic is one and the same.
:::

And so the suspense is set: if OLS estimates neither the ATE nor the ATT, then what exactly does IV with randomized encouragement estimate? The answer, in the next section, is the LATE, an effect on the compliers. To see why it is that, and whether it is worth trusting, we have to lay the four assumptions IV rests on out on the table one at a time.

## 3 Identification

Section 1's randomized encouragement gave us a design that goes around selection, but going around is not the same as identifying. We will see immediately that dividing the two differences (encouragement on adoption, encouragement on engagement) yields a clean and tidy number, but by what right is it the causal effect of adoption, and whose causal effect is it? Answering these two questions rests not on arithmetic but on a set of assumptions about counterfactual worlds. Identification logically precedes estimation. Identification asks whether, under a set of such assumptions, the causal quantity we want can be written as a function of the observed data distribution; that answer depends only on the design and the institutional background, and has nothing to do with which estimator or how large a sample we use afterward. This is the most finely worked section of the whole chapter. We start from the intuition of a ratio, to make plain what IV is actually computing, then dissect one by one the four assumptions it rests on, then prove the LATE theorem, discuss why LATE is not the ATE, and finally set up MTE as the master object that unifies everything, giving the control function as a parametric equivalent route.

### 3.1 A Ratio: First, What IV Is Computing

Lay out the arithmetic of IV first, because it is simpler than you would think. Randomized encouragement is a coin flip, so the encouraged group and the non-encouraged group have no systematic difference of any kind before the encouragement. Comparing two things across the two groups is enough.

The first thing is the effect of encouragement on adoption, called the first stage: the adoption rate of the encouraged group minus the adoption rate of the non-encouraged group. At Meridian, about 42.2% of firms adopt on their own when not encouraged, and the adoption rate rises to 88.8% once encouraged, so the first stage is this rise, about 0.466. The second thing is the effect of encouragement on the outcome, called the reduced form or intent-to-treat (ITT, the intention effect): the average engagement of the encouraged group minus the average engagement of the non-encouraged group. Note that ITT never once mentions adoption; it looks only at the change in engagement brought about by the act of "giving a push," and that change is 0.369.

The crucial step is to divide the two. Encouragement lifted the adoption rate by 0.466 and lifted engagement by 0.369, so for every firm that adopts because of the encouragement, how much engagement improvement does it bring on average?

$$\hat\beta_{IV} = \frac{\text{ITT}}{\text{first stage}} = \frac{0.369}{0.466} = 0.792.$$

This ratio is the Wald estimator. Its intuition is to undo the dilution of the ITT: encouragement makes only 46.6% of people actually change their adoption behavior, and the rest either adopt anyway or refuse no matter what, so for them encouragement has no effect on engagement working through adoption. The ITT spreads everyone's engagement change flat; dividing by the first stage amounts to spreading it only across those who were genuinely moved, giving their effect per unit change in adoption.

The value 0.792 is worth stopping to examine. It is far below the OLS 12.104, below the covariate-adjusted 5.285, and even below the true ATT of 3.056 and ATE of 2.018. Why it is credible, and whose effect it is, are two questions answered respectively by the assumptions and the theorem in the next three subsections. For now, remember the shape of this ratio: a clean numerator (the ITT) divided by a clean denominator (the first stage), and that is the entire machinery of IV.

### 3.2 Four Identifying Assumptions

For the Wald ratio to really equal a meaningful causal quantity, four assumptions are needed. Their division of labor is clean; we state each in turn, explain each with Meridian, and say clearly for each which can be tested and which can only be argued.

::: {.assumption}
**A1 (relevance)** The instrument really affects the treatment: $\text{Cov}(D_i, Z_i) \neq 0$, equivalently $\mathbb{E}[D_i \mid Z_i = 1] - \mathbb{E}[D_i \mid Z_i = 0] \neq 0$.
:::

Relevance is just the first stage being nonzero. Encouragement has to actually push some firms into adopting, otherwise the denominator is zero and the ratio is meaningless. This is the only one of the four assumptions that can be tested directly, because both $D$ and $Z$ are observed, and one just regresses $D$ on $Z$ and checks whether the coefficient is significant. Meridian's first stage is 0.466, encouragement lifted the adoption rate by nearly half, and relevance is beyond question. But how weak relevance has to get before things go wrong is the theme of the weak-instrument discussion in Section 4, where we will see that a seemingly significant first stage can still be nowhere near enough.

::: {.assumption}
**A2 (random assignment / as-good-as-randomly-assigned)** The instrument is independent of all potential outcomes and potential treatments. To let random assignment and the exclusion below each mind its own business, we index the potential outcomes here doubly by the encouragement value $z$ and the adoption value $d$:

$$Z_i \perp \big(\{Y_i(z,d)\}_{z,d \in \{0,1\}},\, D_i(0), D_i(1)\big).$$
:::

A2 says the instrument is as good as randomly assigned. It requires that the assignment of encouragement be unrelated to any of a firm's potential outcomes or any compliance type: the platform cannot have singled out firms that were going to take off anyway for the push. At Meridian this is designed in, since encouragement rests on a coin flip and is independent by construction. In real research, instruments are rarely handed out by a coin flip, and this condition then has to be argued from institutional detail; and it is itself untestable, because what it constrains is the relationship between the instrument and potential outcomes that are never observed.

::: {.assumption}
**A3 (exclusion restriction)** The instrument affects the outcome only through the treatment: for all $z, z', d$, $Y_i(z, d) = Y_i(z', d)$, so that the outcome can be written as a function of the treatment alone, $Y_i(d)$.
:::

Exclusion is the soul of IV, and also the one most easily talked past. It requires that the entire effect of encouragement on engagement run through the single channel of "bringing about adoption," with no other channel at all. At Meridian the threat is concrete: could that onboarding invitation and the tutorials themselves make a firm run its business more attentively, even if it ends up not adopting the dashboard? Could the extra few calls the customer success team placed themselves improve the client relationship and engagement? As long as encouragement retains such a direct channel that bypasses adoption, exclusion is broken, the ITT gets mixed with effects not running through adoption, and the Wald ratio is no longer the effect of adoption.

Here it is worth spelling out the division of labor between A2 (random assignment) and A3 (exclusion), because many textbooks blur them into a single word, "exogeneity," when in fact they are two different things. Random assignment is about whether the instrument is cleanly tossed. Exclusion is about the channel: whether the instrument has only the one path, adoption, to the outcome. An instrument can perfectly well be randomly assigned (A2 holds) yet still have a side door straight to the outcome (A3 fails). The most typical vulnerability of randomized encouragement is exactly this: the coin flip is beyond reproach, yet the act of encouragement itself may carry a placebo-like direct effect. A2 and A3 together are what we usually call instrument exogeneity, and it is exactly the kind of independence the LATE proof later actually uses. Exclusion, like random assignment, is untestable; it is an economic claim about mechanism, one you can only reason about, never let the data prove for you.

::: {.assumption}
**A4 (monotonicity)** The instrument works in a consistent direction for everyone; there are no defiers: $D_i(1) \geq D_i(0)$ holds for all $i$.
:::

Monotonicity requires that the push of encouragement point the same way for every firm, and that nobody perversely refuses to adopt because they were encouraged. What it rules out are the defiers. This condition too is untestable, because it requires knowing at once a firm's adoption status both when encouraged and when not, and we see only one of the two. It is usually secured by a latent-index model: a firm adopts when its adoption propensity crosses a threshold, $D_i(z) = \mathbf{1}\{\pi_0 + \pi_{1,i} z + \varepsilon_i > 0\}$, and as long as the sign of $\pi_{1,i}$ does not vary from person to person (encouragement pushes everyone toward adoption, or at most fails to move them), monotonicity holds. At Meridian this is natural: one more onboarding invitation will not make any firm dig in and refuse out of spite. But it can fail in other settings, and Section 7 gives a counterexample.

::: {.warning}
Of the four assumptions, only relevance can be tested directly. Random assignment, exclusion, and monotonicity are all untestable; what they constrain are counterfactuals that are never observed. This means the credibility of an IV paper rests almost entirely on the argument for these three, and not on any statistical test. Exclusion especially so: it is often waved through with a "of course this instrument works only through the treatment," and that very sentence is the one most in need of being broken down and argued out. In Angrist and Pischke's words, the reason the instrument is excluded "is just a story," and you have to tell the story convincingly enough.
:::

### 3.3 The Four Compliance Types and How Monotonicity Erases the Defiers

With potential treatments in hand, we can sort firms into four types by $\big(D_i(0), D_i(1)\big)$.

| | $D_i(0) = 0$ | $D_i(0) = 1$ |
|---|---|---|
| $D_i(1) = 0$ | never-taker (does not adopt even when pushed) | **defier** (does the opposite) |
| $D_i(1) = 1$ | **complier** (adopts only when pushed) | always-taker (adopts even without a push) |

A never-taker is $(0,0)$, not adopting whether pushed or not. An always-taker is $(1,1)$, adopting whether pushed or not. A complier is $(0,1)$, not adopting without a push and adopting once pushed, the only group the instrument genuinely moves. A defier is $(1,0)$, going out of its way to do the opposite of the encouragement.

The role of monotonicity (A4) is to strike the defier cell from this table. Why it must be struck is something the proof in the next subsection makes very clear: the defiers' effect enters the Wald ratio with a negative sign and cancels against the compliers' effect, and once you allow them to exist, the ratio becomes an uninterpretable soup. After striking the defiers, the two observed groups each mix only two types: those who adopt in the encouraged group ($D=1, Z=1$) are always-takers plus the moved compliers, and those who adopt in the non-encouraged group ($D=1, Z=0$) are pure always-takers. This "only two types left per group" structure is the key that lets the LATE theorem go through.

Here you can already see IV's fundamental limitation, and the line this chapter returns to again and again: IV can learn only about the compliers. The adoption status of always-takers and never-takers does not vary with the instrument at all, so the instrument manufactures no contrast on them, and IV knows nothing of their effects. Whichever compliers this instrument of yours pushes, that is the group you get to learn about. Switch to another instrument, it pushes another group, and what you learn is another effect.

### 3.4 The LATE Theorem

Now connect the ratio of 3.1 with the type classification of 3.3, and prove that the Wald ratio equals precisely the average treatment effect on the compliers. This is the LATE theorem of Imbens and Angrist (1994) and Angrist, Imbens, and Rubin (1996).

::: {.theorem}
**(LATE theorem)** Under the four assumptions of relevance, random assignment (independence), exclusion, and monotonicity,

$$\beta_{IV} = \frac{\mathbb{E}[Y_i \mid Z_i = 1] - \mathbb{E}[Y_i \mid Z_i = 0]}{\mathbb{E}[D_i \mid Z_i = 1] - \mathbb{E}[D_i \mid Z_i = 0]} = \mathbb{E}\big[Y_i(1) - Y_i(0) \mid \text{complier}\big].$$

The numerator is the ITT (reduced form), the denominator is the first stage, and the ratio equals the average treatment effect on the compliers, that is, the LATE.
:::

The derivation goes step by step, and each step uses up one assumption. Start from the numerator. Exclusion lets us write the outcome as $Y_i = Y_i(0) + \tau_i D_i$ (depending only on adoption, not directly on encouragement). Substituting into the ITT leaves an extra term $\mathbb{E}[Y_i(0)\mid Z_i=1]-\mathbb{E}[Y_i(0)\mid Z_i=0]$, and random assignment guarantees that encouragement is unrelated to $Y_i(0)$, so this term is zero, leaving only:

$$\mathbb{E}[Y_i \mid Z_i = 1] - \mathbb{E}[Y_i \mid Z_i = 0] = \mathbb{E}[\tau_i D_i \mid Z_i = 1] - \mathbb{E}[\tau_i D_i \mid Z_i = 0].$$

Independence guarantees that the joint distribution of $\tau_i$ with the potential treatments does not vary with $Z$, so we can replace $D_i$ with the corresponding potential treatment, and the right-hand side becomes $\mathbb{E}\big[\tau_i \big(D_i(1) - D_i(0)\big)\big]$. Now split by the three values of $D_i(1) - D_i(0)$:

$$\mathbb{E}\big[\tau_i (D_i(1) - D_i(0))\big] = \underbrace{\Pr(\cdot = 1)\, \mathbb{E}[\tau_i \mid \cdot = 1]}_{\text{compliers}} + \underbrace{\Pr(\cdot = 0)\cdot 0}_{\text{always/never}} + \underbrace{\Pr(\cdot = -1)\, \mathbb{E}[-\tau_i \mid \cdot = -1]}_{\text{defiers}}.$$

The middle term is mechanically zero, because for always-takers and never-takers $D_i(1) - D_i(0) = 0$, and the instrument manufactures no adoption change on them. The third term is the defiers, carrying a negative sign. Enter monotonicity: it sets $\Pr(D_i(1) - D_i(0) = -1) = 0$, and the defier term vanishes entirely. This is why the defiers must be ruled out, otherwise they would contaminate the numerator with a negative sign. What is left is only the compliers:

$$\mathbb{E}[Y_i \mid Z_i = 1] - \mathbb{E}[Y_i \mid Z_i = 0] = \Pr(\text{complier})\, \mathbb{E}[\tau_i \mid \text{complier}].$$

Finally look at the denominator. The first stage equals exactly the share of compliers, $\mathbb{E}[D_i \mid Z_i = 1] - \mathbb{E}[D_i \mid Z_i = 0] = \Pr(\text{complier})$, because the extra people who adopt in the encouraged group relative to the non-encouraged group are precisely the compliers who adopt only when pushed. Dividing numerator by denominator, the complier share cancels, leaving $\mathbb{E}[\tau_i \mid \text{complier}]$, the LATE. QED.

::: {.intuition}
The intuition of the LATE theorem can be remembered this way. The ITT is "how much engagement one push brings," the first stage is "how much adoption one push brings," and both are clean contrasts at the level of the full sample, because encouragement is random. But in the full sample only the compliers' adoption is genuinely moved; the always-takers' and never-takers' adoption does not budge, and their contribution to both numerator and denominator is zero. So this ratio automatically, and only on the compliers, does its division, computing "one push brings about one adoption, and this adoption brings how much engagement," precisely the compliers' average effect. The instrument is like a spotlight, illuminating only the group it can move, and those it cannot reach have no bearing on the estimate.
:::

### 3.5 LATE Is Not the ATE

LATE is the effect on the compliers, not the effect on everyone. This is no wordplay; it has precise algebra. Denote the individual first-stage strength $\Delta_i = D_i(1) - D_i(0)$ (which under monotonicity takes 0 or 1). One can show that the probability limit of 2SLS is an effect average weighted by first-stage strength:

$$\hat\beta_{IV} \xrightarrow{p} \frac{\mathbb{E}[\tau_i\, \Delta_i]}{\mathbb{E}[\Delta_i]} = \text{ATE} + \frac{\text{Cov}(\tau_i, \Delta_i)}{\mathbb{E}[\Delta_i]}.$$

The first equality says LATE weights each person's effect by "how much the instrument moves them," with people whose adoption propensity is untouched by the instrument ($\Delta_i = 0$, the always/never-takers) getting weight zero and dropping straight out of the estimate. The second equality writes it as the ATE plus a covariance term, making it plain at a glance when LATE equals the ATE: if and only if the individual effect $\tau_i$ is uncorrelated with the individual instrument responsiveness $\Delta_i$, the covariance term is zero and LATE equals the ATE.

At Meridian these two numbers are far apart: LATE is 0.784, ATE is 2.018. The covariance term is $0.784 - 2.018 = -1.234$, significantly negative. Why negative? Because in this DGP, the sharpest firms with the lowest resistance to adoption are both die-hard always-takers and the group that gains the most from the dashboard; whereas the compliers genuinely moved by randomized encouragement are exactly the firms with moderate-to-high resistance and modest gains to begin with. The instrument moves the low-gain people, so the compliers' average effect is below the population average. This is not an estimation error; it is the setting speaking: what randomized encouragement can pry loose is by nature not the die-hards who should adopt most and gain most (they adopt without a push), but the fence-sitters who were undecided in the first place.

This leads straight to a consequence that is deadly in practice: two valid instruments will trace out different complier populations, give different LATEs, and both are right. The classic example is raising college enrollment: use a need-based grant as the instrument, and you pry loose poor students who would not have enrolled only because their pockets were empty (the compliers are the liquidity-constrained, whose returns are not necessarily the highest); use a merit scholarship as the instrument, and you pry loose top students who were just shy of qualifying (the compliers are the high-ability, whose returns may be higher). Both instruments are clean, yet they give different LATEs because they move different people. This is also why we will later see that an over-identification test can reject even when both instruments are perfectly valid: what it rejects is not validity but "the two LATEs are equal," an assumption that under heterogeneous effects should not hold in the first place.

### 3.6 MTE: The Master Object That Unifies Everything

LATE is the effect on the compliers, but who the compliers are depends in turn on the instrument. Is there a more fundamental object of which ATE, ATT, LATE, and even OLS are all different weightings? There is: the marginal treatment effect (MTE), first proposed by Björklund and Moffitt (1987) in a parametric selection model, then systematized through a line of work by Heckman and Vytlacil as its core.

Write the latent-index selection model of 3.2 in tidy form. A firm adopts when its adoption propensity crosses a threshold, which after standardization reads

$$D_i = \mathbf{1}\{U_{D,i} \leq P(Z_i)\}, \qquad P(Z_i) = \Pr(D_i = 1 \mid Z_i),$$

where $U_{D,i}$ is an unobserved adoption resistance standardized to a uniform distribution on $[0,1]$: firms with smaller $U_D$ are more easily pushed into adopting (the sharpest, the keenest to try new things), and larger $U_D$ marks the immovable die-hard hold-outs. Here $P(Z)$ is the propensity score, the adoption probability given the instrument value. The beauty of this formulation is that it organizes "who adopts" into a clean line: everyone with resistance below the threshold $P(Z)$ adopts, and the instrument's role is to move this threshold. Encouragement pushes the threshold from $P(0) = 0.421$ up to $P(1) = 0.885$, and the firms swept over by this stretch of threshold movement ($0.421 \leq U_D < 0.885$) are the compliers.

The MTE is defined as the marginal effect taken along this resistance dimension: the average treatment effect of the group of firms standing exactly on the edge of adopting or not ($U_D = u$, that is, indifferent exactly when $P(Z) = u$),

$$\Delta^{MTE}(u) = \mathbb{E}\big[Y_i(1) - Y_i(0) \mid U_{D,i} = u\big].$$

It is a function of resistance $u$, drawing the treatment effect as a curve along the spectrum "from most willing to adopt to least willing to adopt." In Meridian's DGP this curve is $\Delta^{MTE}(u) = 6 - 8u$: the least-resistant firms ($u \to 0$) gain 6 points from adopting the dashboard, the effect hits zero at resistance $u = 0.75$, and beyond that (the most stubborn firms) adoption is even slightly harmful. That the MTE slopes downward is itself the signature of selection on gains: the more a firm resists adoption, the smaller its return to adopting, which is why it resists in the first place.

The MTE has three equivalent characterizations, each giving one kind of intuition. The first is as a limit of LATE (for the formal result see Heckman and Vytlacil 1999, 2005, and also Heckman 1997): do a Wald between two infinitesimally close instrument values, $\Delta^{MTE}(P(z)) = \lim_{z' \to z} \text{Wald}(z, z')$, so that at the threshold $u = P(z)$, MTE is a LATE that is as local as it can possibly get, prying loose the thin layer of compliers at the threshold. The second is local IV: $\Delta^{MTE}(p) = \partial\, \mathbb{E}[Y \mid P(Z) = p] / \partial p$, differentiating the conditional mean of the outcome in the propensity. The third comes from direct derivation. Split the adoption gain into $Y_i(1) - Y_i(0) = \text{ATE} + (U_{1i} - U_{0i})$, where $U_{1i} - U_{0i}$ is the unobserved heterogeneous part of the gain (the amount by which the individual gain deviates from the ATE). Inside $\mathbb{E}[Y \mid P(Z) = p]$ sits a term $\mathbb{E}[D\,(U_1 - U_0) \mid P(Z) = p]$ (the average unobserved gain of adopters at the threshold), and differentiating it in $p$ gives exactly $\Delta^{MTE}(p) = \text{ATE} + \mathbb{E}[U_1 - U_0 \mid U_D = p]$, the latter term being the deviation-from-average unobserved gain of the group at the threshold.

Vytlacil's (2002) equivalence theorem endorses this framework: the four LATE assumptions (essentially independence plus monotonicity) are completely equivalent to the constraints that this nonparametric latent-index selection model imposes on counterfactual data. That is, whether you use the language of LATE or the language of the selection model, the identifying power is the same, and LATE is no more general than a Heckman-style selection model. The difference is only in the framing.

The most elegant conclusion is that every treatment-effect parameter we know is a weighted integral of this one MTE curve, only the weights differing:

$$\text{ATE} = \int_0^1 \Delta^{MTE}(u)\, du,$$

$$\text{LATE}(z, z') = \frac{1}{P(z) - P(z')} \int_{P(z')}^{P(z)} \Delta^{MTE}(u)\, du,$$

$$\text{ATT} = \int_0^1 \Delta^{MTE}(u)\, \omega_{ATT}(u)\, du, \qquad \omega_{ATT}(u) = \frac{\Pr(P(Z) > u)}{\mathbb{E}[P(Z)]}.$$

The ATE averages the whole curve evenhandedly. The LATE averages uniformly only over the threshold interval $[P(z'), P(z)]$ of the compliers. The ATT gives more weight to the low-resistance end, because adopters are concentrated at the easy-to-adopt end to begin with. OLS is a weighted integral too, only that its weights carry an extra selection term $\mathbb{E}[U_1 \mid U_D = u]$ and $\mathbb{E}[U_0 \mid U_D = u]$, and this is exactly the algebraic root of why OLS, under selection on gains, is neither the ATE nor the ATT. A picture is worth a thousand words: ATE, ATT, and LATE are three ways of weighing the same MTE curve, measuring out three different numbers. At Meridian these three numbers are ATE 2.018, ATT 3.056, LATE 0.784, and Section 6 will reduce them all to weighted integrals of this one curve.

To identify the entire MTE curve, the price is not small: it requires the propensity $P(Z)$ to sweep the whole $(0,1)$ interval (full support), that is, an instrument strong enough to push the adoption threshold all the way from near 0 to near 1. This is far stronger than the relevance of ordinary IV. In real data instruments are rarely this powerful, and what can be identified is often only a segment of the curve, which is also why MTE estimation is data-hungry by birth.

### 3.7 The Control Function and the Forbidden Regression

That residual term inside the propensity score of the MTE machinery is in fact the idea of the control function, a parametric equivalent route to IV. The thought is to grab the first-stage residual and stuff it into the outcome equation as a control variable. Written out it is two steps:

$$D_i = \pi_0 + \pi_1' Z_i + v_i, \qquad Y_i = \beta_0 + \beta_1 D_i + \rho\, \hat v_i + e_i,$$

where $\hat v_i$ is the first-stage residual. The intuition is this: adoption is endogenous because the part of $D_i$ not explained by the instrument ($v_i$) is correlated with the outcome disturbance; put $v_i$ explicitly into the outcome equation and control it, and the remaining variation in $D_i$ is clean. Under joint normality or linear-projection assumptions, conditioning on $\hat v_i$ washes out the endogeneity and $\hat\beta_1$ is consistent. Moreover the coefficient $\rho$ is itself an endogeneity test: $\rho = 0$ says adoption is in fact exogenous and OLS suffices; $\rho \neq 0$ says there is indeed endogeneity, and this is exactly the same intuition behind the Hausman test and the Heckman selection correction. In the parametric MTE model, the control function and the MTE machinery are one and the same: that term $\mathbb{E}[D\,(U_1 - U_0) \mid P(Z) = p]$ inside $\mathbb{E}[Y \mid P(Z) = p]$ is a control function written into the propensity, and differentiating it gives the MTE. The virtue of the control function is that it generalizes cleanly to nonlinear, discrete-dependent-variable outcomes that 2SLS cannot handle, where there is no ready two-stage least squares.

::: {.warning}
Buried near the control function is a classic trap, the forbidden regression. When the first stage is nonlinear (say a probit fit to a binary adoption), some people will stuff the probit's fitted value $\hat D_i$ straight into the second stage, or worse, take a nonlinear transform of the fitted value $g(\hat D_i)$ in place of $\widehat{g(D_i)}$. This is inconsistent, because projection does not pass through a nonlinear function, $\mathbb{E}[g(D) \mid Z] \neq g(\mathbb{E}[D \mid Z])$. There are only two safe roads: either honestly use the fitted value of a linear first stage in linear 2SLS, or use an explicit control function, but never substitute a nonlinear fitted value. The same trap has a smaller version: even if the first stage is linear, do not manually run two regressions yourself, because then the second stage's standard errors are wrong (they do not account for the uncertainty of the first-stage estimate). The correct practice is to use an integrated IV command (such as `feols(Y ~ 1 | D ~ Z)`), letting the software construct the residual from the true $D$ and correct the variance.
:::

The main points of this section summarize as follows: the Wald ratio is the ITT divided by the first stage; under relevance, random assignment, exclusion, and monotonicity it equals the compliers' LATE, of which only relevance is testable; LATE is generally not the ATE, but the ATE plus a covariance term between the effect and instrument responsiveness, and different instruments trace different compliers; MTE is the marginal effect along the unobserved-resistance dimension, and ATE, ATT, LATE, and OLS are all its weighted integrals; the control function is IV's parametric equivalent route, to be used while steering clear of the forbidden regression.

## 4 Estimation

Identification made clear what IV estimates; this section covers how to estimate it, and under what circumstances it cannot be estimated at all. The through-line remains "where the last method fails is what gives birth to the next."

### 4.1 From Wald to 2SLS

With only a single binary instrument, estimation is just the Wald ratio of 3.1, using sample means in place of population means. When the instrument is continuous, or there are multiple instruments, or there are exogenous covariates along for the ride, we need the more general machine of two-stage least squares (2SLS). It is just two steps: the first stage regresses the endogenous adoption $D$ on the instrument $Z$ (together with the exogenous covariates) and takes the fitted value $\hat D$; the second stage regresses the outcome $Y$ on $\hat D$, and the coefficient is the IV estimate. In matrix form, $\mathbf{b}_{2SLS} = (\hat{\mathbf{X}}'\hat{\mathbf{X}})^{-1}\hat{\mathbf{X}}'\mathbf{y}$, where $\hat{\mathbf{X}} = \mathbf{Z}(\mathbf{Z}'\mathbf{Z})^{-1}\mathbf{Z}'\mathbf{X}$ is the fitted value of the endogenous variable projected onto the instrument space. When just-identified (the number of instruments equals the number of endogenous variables), 2SLS collapses to $\mathbf{b}_{IV} = (\mathbf{Z}'\mathbf{X})^{-1}\mathbf{Z}'\mathbf{y}$, and with a single binary instrument it is the Wald ratio, numerically identical.

The intuition for using the fitted value is to pull out the clean piece of the endogenous variable on its own. $D$ contains both the exogenous variation pushed by the instrument and the endogenous variation tangled with the disturbance, and $\hat D$ keeps only the former, so using it to explain $Y$ goes around the endogeneity. This also explains why the second stage's standard error cannot use what OLS spits out directly: $\hat D$ is estimated, and ordinary OLS standard errors treat it as known data and get it wrong. The IV residuals satisfy $\mathbf{Z}'\mathbf{e}_{IV} = 0$ (the moment condition $\mathbb{E}[z_i e_i] = 0$, which is the GMM viewpoint, with OLS corresponding to $\mathbb{E}[x_i e_i] = 0$), and an integrated IV command uses it to get the variance right.

### 4.2 The Price of IV: Precision

IV buys consistency and sells precision. Under homoskedasticity the two variances compare directly:

$$\text{Var}(\hat\beta_{2SLS}) = \frac{1}{N}\frac{\sigma_\varepsilon^2}{\sigma_X^2\, R^2_{XZ}}, \qquad \text{Var}(\hat\beta_{OLS}) = \frac{1}{N}\frac{\sigma_\varepsilon^2}{\sigma_X^2}.$$

The two differ only by an $R^2_{XZ}$, the fraction of the endogenous variable explained by the instrument. Because $R^2_{XZ} \leq 1$, the variance of IV is never smaller than that of OLS, and the weaker the instrument ($R^2_{XZ} \to 0$) the more the variance explodes. The reason echoes 3.5: IV uses only the piece of treatment variation pushed by the instrument and throws away the rest, so with less information, larger standard errors follow naturally. At Meridian the instrument is extremely strong, and the SE of 2SLS (0.247) is still in the acceptable range; let the instrument weaken and the price shows itself at once, which is the theme of the next subsection.

### 4.3 Over-identification and Optimal Weighting

When the number of instruments exceeds the number of endogenous variables (over-identified), there are multiple moment conditions available, and combining them with optimal weighting is GMM (generalized method of moments), of which 2SLS is the special case under homoskedasticity. The extra instruments allow an over-identification test, the Sargan $J$ statistic: if all instruments are valid, the coefficients estimated using different subsets of instruments should agree, and the $J$ test checks whether their differences are large enough to be more than sampling noise.

But the Sargan $J$ has already fallen out of favor under heterogeneous effects, for exactly the reason discussed in 3.5: different instruments trace different compliers and give different LATEs, and $J$ will reject because two LATEs are unequal, even when both instruments are perfectly valid. So a $J$ rejection cannot be read simply as "some instrument is invalid"; it may only be "two instruments illuminate different populations." By the same token, the Durbin-Wu-Hausman test used to check endogeneity (comparing whether OLS and 2SLS agree) will also reject because OLS estimates some other quantity and 2SLS estimates the LATE, and what it rejects is "the two are equal," not "endogeneity exists." Both of these tests must be read with care.

### 4.4 Weak Instruments: Modern Practice

Weak instruments are IV's most treacherous failure mode, and it is not as simple as "estimated imprecisely," but rather "the estimate cannot be inferred from by conventional methods." When the first stage is near zero, three bad things happen, and they must be told apart. The first is a finite-sample bias toward OLS: even if the instrument is perfectly clean, 2SLS in finite samples biases toward OLS, with the bias approximately equal to the OLS bias times $1/(F+1)$ ($F$ being the first-stage F statistic, see Angrist-Pischke), and the smaller $F$, the more it biases. The second is amplification of any uncleanliness: if the instrument is even a touch unclean ($\text{Cov}(Z,\varepsilon)$ slightly nonzero), the ever-smaller denominator in the probability limit $\hat\beta^{IV} \xrightarrow{p} \beta_1 + \frac{\text{Cov}(Z, \varepsilon)}{\text{Cov}(Z, D)}$ blows this bit of contamination up out of control, even flipping the sign, and this is an asymptotic malady, different from the previous one. The third is that $\hat\beta^{IV}$ becomes the ratio of two correlated normals, the normal approximation of the central limit theorem fails, and the conventional $t$ test and confidence interval are badly distorted.

Meridian builds a weak-instrument variant to act out this malady: turn the push of encouragement on adoption down to almost zero, and the first stage is left at just 0.011. A single estimate of 2SLS is 2.817, with an SE as high as 9.480 and a confidence interval so wide as to carry no information. More telling is repeated sampling: rerun this weak design 1000 times, and the median of 2SLS is pushed from the true value 0.784 to 3.382, so the direction does bias toward OLS (consistent with weak-instrument theory), but what is truly lethal is not this bit of shift, rather that the mean reaches 29.133 with tails so heavy as to have almost no finite moments, and the whole sampling distribution has already rotted. The coverage of the conventional 95% confidence interval for the true value is only 0.931, short of the nominal 0.95; at times like this one should switch to weak-instrument-robust inference such as Anderson-Rubin. This is a weak instrument: a single estimate looks respectable, yet the distribution is already unusable. It has a finite-sample side (bias toward OLS, easing as the sample grows) and an asymptotic side (amplifying the instrument's uncleanliness), and "how strong is strong enough" depends on $F$, not on whether the first stage is significant.

Diagnosis and remedy have a set of modern tools, far more discerning than the antiquated "F greater than 10."

| Tool | What it does | When to use |
|---|---|---|
| first-stage F | Report the first-stage F of the excluded instruments (after removing exogenous regressors) | The basic checkup for relevance; with multiple endogenous variables each first stage must pass |
| effective F (Montiel-Olea-Pflueger 2013) | A robust substitute when the first-stage F fails under heteroskedasticity/clustering, compared against MOP critical values | When errors are non-homoskedastic, replacing the mechanical "10" |
| tF correction (Lee et al. 2022) | Inflate the critical value of the $t$ test by the observed first-stage F | When you want an honest confidence interval without setting a hard threshold |
| Anderson-Rubin | Test $H_0: \beta = \beta_0$ by regressing $Y - \beta_0 D$ on $Z$ and checking whether the instruments are jointly significant | When the instrument may be weak, for size-correct robust inference |

Traditionally, for a single endogenous variable one checks whether the first-stage F exceeds 10 (the Staiger-Stock rule of thumb). This threshold is now considered "large enough in practice, but probably not enough in theory": Lee, McCrary, Moreira, and Porter point out that for the conventional $|t| > 1.96$ to truly correspond to a 5% test level, the first-stage F must be far greater than 10 (roughly above 100). Their tF procedure writes the critical value of $t$ as a smooth function of F, with a smaller F giving a larger critical value and a wider interval, so that an honest confidence interval is delivered without setting a hard threshold. Under heteroskedasticity or clustering, Staiger-Stock's F does not itself hold, and one should switch to Montiel-Olea and Pflueger's effective F, comparing against their critical values computed for a target worst-case bias. The most robust move is the Anderson-Rubin test: it tests a specific value $\beta_0$, its size is correct whether the instrument is strong or weak, and inverting it gives a weak-instrument-robust confidence set (which may be unbounded or even empty, and that itself is honest information). In Meridian's weak variant, the conventional interval is $[-15.764, 21.398]$, while the Anderson-Rubin interval is $[-\infty, 1317.636]$; the latter honestly admits with an unbounded interval that "this instrument simply cannot pin down $\beta$," while the former hands you a range that looks finite but deceives.

There is a neighboring class of problem, many-instrument bias: once the number of instruments is large relative to the sample size, the fitted value $\hat D$ of 2SLS overfits the disturbance, and the estimate biases toward OLS. Here LIML (limited-information maximum likelihood) is approximately median-unbiased with better weak-instrument properties (at the cost of fatter tails), while JIVE (jackknife IV) removes the overfitting of a unit's own observation via a leave-one-out approach. These few will come up again in Section 7.

### 4.5 The Control Function Two-Step

The control function of 3.7 lands on estimation as two steps: first run the first stage and store the residual $\hat v$, then put $\hat v$ into the outcome equation as an extra regressor. Because the second step uses the residual estimated in the first, using the second step's OLS standard errors directly understates the uncertainty, and the correct practice is to bootstrap the whole two-step procedure, or use an analytic variance that folds in the estimation uncertainty of both steps. Its extra benefit relative to 2SLS is that it can handle nonlinear outcomes and test endogeneity along the way with $\rho$, at the cost of having to posit one more functional-form assumption about how $v$ enters the outcome, which 2SLS does not need.

## 5 Anchoring Papers

A method only stands firm when it lands in real research. Three anchoring papers: one a classic vehicle for LATE, one a cautionary cause célèbre for weak instruments, and one a representative work in the IS field that treats randomized encouragement as an instrument. Each is combed through along the five elements of paper, method, data, results, and limitations, with a focus on how the assumptions are defended.

### 5.1 Angrist (1990)

::: {.case}
Paper: "Lifetime Earnings and the Vietnam Era Draft Lottery: Evidence from Social Security Administrative Records," American Economic Review. It is the most classic empirical vehicle for the LATE theorem, turning instrumental variables from an old trick of demand estimation into a tool for identifying treatment effects.

Method: to estimate the causal effect of military service on later earnings, but who serves is not random, and whether one enlists is tangled with family background, ability, and employment prospects, so comparing veterans with non-veterans directly is contaminated by selection. Angrist uses the Vietnam-era draft lottery as the instrument: lottery numbers are randomly assigned by birth date, and lower numbers are more likely to be drafted. The lottery number is $Z$, service is the endogenous $D$, and earnings are $Y$. The exclusion story is that the lottery number can affect earnings only by affecting whether one serves, and the number itself does not enter the wage equation; monotonicity is that a lower number only increases, never decreases, the probability of service.

Data: administrative earnings records from the Social Security Administration, matched to lottery numbers defined by birth date.

Results: the lottery number has a strong first stage on service and a measurable reduced form on earnings, and the Wald ratio shows that white veterans earned about 15% less per year than non-veterans in the early 1980s. The compliers here are those who "served only because they drew a low number and would not have served otherwise," and what is estimated is the LATE on them, not the effect on all veterans, and certainly not the effect on voluntary enlistees.

Limitations: LATE's external validity is the core limitation. The compliers the lottery pries loose need not represent all veterans, and voluntary enlistees (always-takers) and those who would never serve no matter what (never-takers) are not in this number. Exclusion also has debatable points, for instance a low number might make people stay in school longer to dodge the draft, thereby affecting earnings while bypassing service.
:::

This paper's defense strategy is the model for IV: the randomness of the instrument comes from an institution (the lottery), exclusion is argued from mechanism rather than statistical test, and the estimate's scope is honestly confined to the compliers, without overreaching to speak of everyone.

### 5.2 Angrist and Krueger (1991) and the Weak-Instrument Cause Célèbre

::: {.case}
Paper: "Does Compulsory School Attendance Affect Schooling and Earnings?", Quarterly Journal of Economics. It is both an ingenious natural experiment and the number-one cautionary case for the weak-instrument problem.

Method: to estimate the return to one more year of schooling, and schooling is endogenous (ability affects both years of schooling and earnings). The authors use quarter of birth (QOB) as the instrument: compulsory schooling laws permit dropping out by age, while the age of school entry is drawn by birth date, so children born early in the year have already read a few months more than those born late by the time they reach the legal dropout age, and are forced to attend a bit more school. Quarter of birth thus becomes an instrument for schooling, the rationale being that it is nearly random yet affects earnings only through years of schooling.

Data: a large sample from the US Census, with quarter of birth, years of schooling, and earnings.

Results: the main specification estimates a return to schooling of about 7% to 8%, close to OLS. The authors also interact quarter of birth with year of birth and state of birth, creating hundreds of instruments to improve precision.

Limitations: it is precisely those hundreds of instruments that laid the mine. Bound, Jaeger, and Baker (1995) point out that quarter of birth is a very weak instrument, with a feeble first stage, and combined with the many-instrument bias of over-identification, 2SLS biases toward OLS, so a seemingly precise estimate is in fact unreliable. They ran a beautiful and cruel demonstration: replace the real quarter of birth with computer-generated random fake instruments (with no causal effect on schooling whatsoever), and one can still "reproduce" comparable estimates and significance. This shows that much of the original estimate's apparent solidity comes from a mechanical illusion of weak plus many instruments, not from genuine identification.
:::

This cause célèbre is the real-world version of Section 4.4 of this chapter: a first stage that looks significant does not mean the instrument is strong enough, many instruments aggravate the bias toward OLS, and random fake instruments reproducing the result is the most glaring diagnostic of weak instruments. It also pushed the whole field toward tools like the effective F and Anderson-Rubin.

### 5.3 Bapna and Umyarov (2015)

::: {.case}
Paper: "Do Your Online Friends Make You Pay? A Randomized Field Experiment on Peer Influence in Online Social Networks," Management Science 61(8):1902-1920. It is a representative work in the IS field of using randomized gifting as exogenous variation to identify causality, sharing a core with Meridian: relying on a randomly assigned encouragement to render the endogenous paying status exogenous.

Method: to estimate peer influence in online social networks, that is, whether a friend becoming a paying member really makes one upgrade to paying oneself, rather than mere correlation. The difficulty is that a friend's paying is highly positively correlated with one's own paying, but this could be genuine influence or just birds of a feather (homophily) and common shocks, and observational comparison cannot tell them apart. The authors worked with the music community Last.fm, randomly selecting about 1000 people from roughly 3.8 million users to be gifted a stretch of premium membership (the manipulated group), with another roughly 1000 as controls. The random gift makes this batch of users exogenously become premium (as good as randomly assigned), so one can cleanly see whether their friends are more likely to upgrade along with them. The random gift is the encouragement that renders the friend's paying status exogenous.

Data: a panel of about 3.8 million Last.fm users, including the random-gift assignment, the friendship network, paying status, and social activity.

Results: at the observational scope, a friend's paying is highly positively correlated with one's own paying, but the causal peer influence given by the randomized experiment is much smaller and still significantly positive: the friends of those randomly gifted do have a higher probability of upgrading to paying. This confirms that the naive correlation mixes in a great deal of homophily, and is not pure peer influence.

Limitations: what is estimated is the effect on the portion of people pried loose by the random gift, which need not extrapolate to everyone; exclusion requires arguing that receiving the gift does not directly affect friends' behavior through some other channel (say being specially favored by the platform); the external validity of a single platform is the usual reservation.
:::

Put the three together and the meaning of anchoring becomes clear: Angrist shows how a lottery gives a clean LATE and honestly confines its scope, Angrist-Krueger and its critique show how weak instruments make a seemingly precise estimate collapse, and Bapna and Umyarov show how a platform's random gift pulls peer influence out of homophily. All three share one core: find an as-good-as-randomly-assigned piece of variation and render the endogenous treatment status exogenous. Meridian puts this core to the most direct question of all: randomly encourage adoption, then look at the effect of adoption on a firm itself.

## 6 A Complete Walkthrough on the Meridian Data

Now run the earlier tools on Meridian from start to finish. The code below uses R 4.5.3 with the random seed fixed at 404, and every number cited in the text comes from the actual run output of this code. IV estimation uses the IV syntax of `fixest::feols`.

### 6.1 The Data-Generating Process

The design parameters are as follows: a Heckman-Vytlacil threshold-crossing model with $n = 30000$ firms. One unobserved adoption resistance $U_D \sim \text{Uniform}(0,1)$ drives adoption, baseline engagement, and the adoption gain all at once, and this is the "sharpness" that observable covariates cannot fix. The smaller $U_D$, the sharper and the quicker to adopt when pushed; the larger, the more of a stubborn hold-out.

```r
set.seed(404); n <- 30000
U_D <- runif(n)                                  # unobserved adoption resistance
W   <- -0.7 * qnorm(U_D) + sqrt(0.51) * rnorm(n) # firm size, correlated with sophistication
Y0  <- 60 + 5 * W - 12 * U_D + rnorm(n, 0, 3)     # baseline engagement
tau <- 6 - 8 * U_D + rnorm(n, 0, 1.5)             # adoption gain, declining in resistance

Z   <- rbinom(n, 1, 0.5)                          # randomized encouragement (binary instrument)
a0 <- -0.20; a1 <- 1.40
D   <- as.integer(pnorm(a0 + a1 * Z) > U_D)       # adoption: threshold crossing, monotone in Z
Y   <- Y0 + D * tau                               # observed engagement
```

A few design intentions need explaining. The gain $\tau = 6 - 8 U_D$ declines in resistance, and the baseline $Y_0$ also declines in resistance, so adopters (low $U_D$) look good for two reasons: high true gain and a high starting point, which is why OLS is severely biased upward. The adoption thresholds $P(0) = \Phi(a_0) = 0.421$ and $P(1) = \Phi(a_0 + a_1) = 0.885$ mean that the medium-to-high-resistance firms swept over by this stretch of threshold are the compliers, and their gains are modest to begin with, so the ATT exceeds the ATE and the ATE exceeds the LATE, forming a stark ordering. Firm size $W = -0.7\,\Phi^{-1}(U_D) + \sqrt{0.51}\,\varepsilon$ is correlated with sharpness, so covariate adjustment can shave off part of the bias but cannot shave it clean, because $W$ is only an imperfect shadow of $U_D$. The threshold is monotone in $Z$ by construction, and there are no defiers.

The true values computed from the individual gains: the ATE is 2.018, the ATT is 3.056, the LATE on the compliers is 0.784, and the MTE curve is $6.000 - 8.000\, u$. These few numbers are the target every estimator that follows is reconciled against.

### 6.2 OLS's Two Bad Answers

```r
b_ols <- coef(lm(Y ~ D))["D"]                    # naive
b_adj <- coef(lm(Y ~ D + W))["D"]                # control for firm size
```

Naive OLS gives 12.104 (SE 0.102), and controlling for firm size drops it to 5.285 (SE 0.063). The true ATT is 3.056, the selection bias is $12.104 - 3.056 = 9.048$, and the absurdly large number of Section 1 is confirmed here. Controlling for size shaves off the part of the bias size can proxy, yet it still exceeds the ATT by 2.230, and that 2.230 is the residual attached to the unobserved resistance $U_D$ that size cannot reach. This is exactly the suspense of Chapter 9: selection lands on the unobservables, and covariate adjustment cannot save it.

### 6.3 IV: First Stage, ITT, Wald, and 2SLS

```r
fs  <- feols(D ~ Z)                              # first stage
itt <- feols(Y ~ Z)                             # reduced form / ITT
iv  <- feols(Y ~ 1 | D ~ Z)                     # just-identified 2SLS
F_iid <- fitstat(fs, "wald", vcov = "iid")$wald$stat
F_eff <- fitstat(fs, "wald", vcov = "hetero")$wald$stat
```

The first stage is 0.466 (SE 0.005), and encouragement lifts the adoption rate from 42.2% to 88.8%. The first-stage F reaches about 9483 (effective F about 9486), and the instrument is strong beyond dispute. The ITT is 0.369 (SE 0.117), that is, a random push brings on average a 0.369-point lift in engagement. The Wald ratio, ITT divided by first stage, is $0.369 / 0.466 = 0.792$ (SE 0.252), and the integrated 2SLS gives 0.792 (SE 0.247), the two agreeing and both hugging the true LATE of 0.784 (off by 0.008).

Worth savoring is the drop across this string of numbers: OLS says 12.104, controlling for size says 5.285, and clean IV says 0.792. IV pulls the adoption effect out of the mire of selection, and what it pulls out is a far smaller number, because the compliers randomized encouragement pries loose are exactly the least-gaining fence-sitters.

### 6.4 Who the Compliers Are

The first stage of 0.466 is the estimate of the complier share: about 42.2% adopt on their own before encouragement (always-takers), encouragement lifts the adoption rate to about 88.8%, and the roughly 46.6 percentage points lifted are the compliers, with the remaining roughly 11.2% being never-takers. To describe what the compliers look like, one can use Abadie's $\kappa$ weighting, taking the mean of any observable characteristic over the complier subpopulation.

```r
pz <- mean(Z)
kappa <- 1 - D * (1 - Z) / (1 - pz) - (1 - D) * Z / pz
w_complier <- sum(kappa * W) / sum(kappa)        # E[W | complier]
```

The estimate is $\mathbb{E}[W \mid \text{complier}] = -0.299$ (oracle true value $-0.296$), while the mean $W$ over all firms is about 0.006. Compliers are smaller-than-average firms, which matches intuition: large firms tend to be always-takers who adopt without a push, and what randomized encouragement pries loose are the small firms that were originally hesitant. Being able to characterize the compliers' features is a practical way to ease the extrapolation anxiety of LATE: at least one can say clearly whom this local effect is local to.

![Left: randomized encouragement lifts the adoption rate from about 42.2% (Z=0) to about 88.8% (Z=1), and the lifted 0.466 is the first stage, also the estimate of the complier share. Right: cut the line of unobserved adoption resistance U_D into three segments by adoption propensity, with always-takers about 42.2%, compliers about 46.6%, and never-takers about 11.2%, and IV sees only the middle segment of compliers.](assets/fig/fig_12_compliers.svg)

### 6.5 The Weak-Instrument Variant

```r
a1_weak <- 0.028                                 # turn the push down to almost zero
D_weak  <- as.integer(pnorm(a0 + a1_weak * Z) > U_D)
iv_weak <- feols(Y ~ 1 | D_weak ~ Z)
```

After turning the push weak, the first stage is left at just 0.011, and the first-stage F is 3.77 (effective F 3.77), far below any reasonable threshold. A single 2SLS is 2.817 with an SE as high as 9.480 and a conventional 95% interval of $[-15.764, 21.398]$; the Anderson-Rubin interval is $[-\infty, 1317.636]$, honestly admitting with unboundedness that it cannot be pinned down. Repeated sampling 1000 times gives a 2SLS median of 3.382 and mean of 29.133 (tails so heavy as to have almost no finite moments), with the median biased toward OLS (about 15.9) rather than the true value 0.784, and the conventional interval's coverage of the true value only 0.931. This confirms Section 4.4: under weak instruments a single estimate looks respectable, but the sampling distribution has in fact rotted, and robust inference must rest on Anderson-Rubin rather than the conventional interval.

### 6.6 MTE: One Curve, Three Weightings

To draw the entire MTE curve, one needs an instrument that can sweep the full range of adoption propensity, so we switch to a continuous instrument (think of it as a randomly assigned tier of onboarding-push strength), whose propensity sweeps from near 0 to near 1, achieving full support.

```r
Zc  <- runif(n)                                            # continuous instrument
D_c <- as.integer(pnorm(-2.40 + 4.80 * Zc) > U_D)         # propensity covers (0,1)
P   <- pnorm(-2.40 + 4.80 * Zc)                            # estimated propensity
# local IV: differentiate the conditional mean of Y in P (slope of a quadratic)
liv <- lm(Y ~ poly(P, 2), subset = (D_c == 1 | D_c == 0))
# control function two-step: residual into the outcome equation, recover MTE intercept and slope
```

The MTE estimated by local IV is $5.661 - 7.098\, u$ (slope SE 1.346), right in direction and magnitude but noisy, because the propensity piles up near 0 and 1 and is sparse in the middle, and nonparametric differentiation is data-hungry by birth (at $n = 3{,}000{,}000$ the same estimator's slope converges to $-8.03$, proving it consistent, only with large variance at this sample size). The control function two-step gives a more precise curve, $6.135 - 8.184\,u$ (implied ATE 2.043) (bootstrap SE 0.079), hugging the true value.

The most elegant step is to integrate this curve under three weightings, recovering three familiar quantities: uniform weight (the full interval) gives ATE 2.037 (the result of numerical grid integration of the line above, with the small difference from its analytic value 2.043 coming only from grid discretization; true value 2.018); averaging uniformly only over the complier threshold interval $[0.421, 0.885]$ gives LATE 0.794 (2SLS 0.792, true value 0.784); and using the low-resistance-weighted ATT weights gives 2.954 (the corresponding continuous-design oracle ATT is 2.922). The same MTE, weighed differently, measures out three different numbers, and this is the literal redemption of "everything is a weighted MTE."

![The MTE curve drawn along unobserved adoption resistance U_D (true value 6-8u, the control-function estimate, and the noisy local-IV dashed line). The low-resistance sharp firms gain the most, and the stubborn firms with resistance above about 0.75 are harmed by adopting. The ATE (uniform over the full interval, 2.04), LATE (averaged over the complier threshold window [0.42,0.88], 0.79), and ATT (weighted toward the low-resistance end, 2.95) are all weighted integrals of this same curve.](assets/fig/fig_12_mte.svg)

### 6.7 The Full Reconciliation

Put all the estimators together, with the target being the true LATE 0.784 (and, as reference, the ATE 2.018 and ATT 3.056):

| Estimator | Estimate | SE | Target |
|---|---|---|---|
| Naive OLS | 12.104 | 0.102 | none (biased by 9.048) |
| OLS controlling for firm size | 5.285 | 0.063 | none (still exceeds ATT by 2.230) |
| ITT (reduced form) | 0.369 | 0.117 | a diluted version of complier LATE |
| Wald = ITT / first stage | 0.792 | 0.252 | LATE |
| 2SLS (feols IV) | 0.792 | 0.247 | LATE (true value 0.784) |
| Weak-instrument 2SLS (single) | 2.817 | 9.480 | broken, cannot infer |

The reconciliation of this section summarizes as follows: naive OLS is propped up by selection to 12.104, controlling for the observable firm size can only press it down to 5.285, still unable to remove the residual attached to the unobserved resistance; IV with randomized encouragement, first stage 0.466, ITT 0.369, divided gives 0.792, hugging the true complier-LATE 0.784; this number is far below the ATE and ATT, because what encouragement pries loose are the least-gaining fence-sitters; switch the instrument to a continuous version that can sweep full support, and the entire MTE curve can be drawn, with ATE, ATT, and LATE all becoming its weighted integrals; and turn the instrument weak, and the whole inference machinery breaks at once, with only Anderson-Rubin still standing.

## 7 Failure Modes and Robustness

In a simulation the assumptions are built in, and in real research they can fail at any moment. This section combs through the most common ways of failing and the actionable responses.

Weak instruments are the number-one threat, and their destructive power has already been quantified above. There are three points of practical discipline. First, report the strength of the first stage, and report it right: the first-stage F of the excluded instruments (after removing the exogenous regressors), switched to the effective F under heteroskedasticity or clustering, and stop worshiping the "F greater than 10" old threshold that is theoretically not up to standard. Second, when the instrument may be weak, base the main inference on weak-instrument-robust methods like Anderson-Rubin rather than the conventional $t$ interval; the interval Anderson-Rubin gives, even if unbounded, is more honest than a conventional interval that looks finite but is in fact distorted. Third, avoid manufacturing a big pile of instruments for the sake of precision, since many instruments aggravate the bias toward OLS, and this is exactly the lesson of the Angrist-Krueger cause célèbre: better one strong instrument than a hundred weak ones.

Exclusion violation is the most insidious threat, because it is untestable and can only be argued. At Meridian the threat is concrete: that onboarding invitation and the tutorials themselves may raise engagement, even if the firm ends up not adopting the dashboard; the extra calls the customer success team placed themselves improve the client relationship. As long as encouragement retains such a direct channel bypassing adoption, the ITT gets mixed with effects not running through adoption, and the Wald ratio is no longer the effect of adoption. The response can only be design and argument: design the encouragement to be as "pure" as possible (doing only actions directly related to adoption, carrying no other kind of care alongside), argue institutionally that no side channel exists, and if necessary find a cleaner instrument with no suspicion of a direct effect. No statistical test can prove exclusion for you, and this point must be faced head-on in the paper, with a whole paragraph making clear by what right this instrument is excluded.

Monotonicity failure means defiers exist, and the Wald ratio loses its LATE reading. It is reasonable in many settings (one more push will not make anyone dig in and refuse out of spite), but it can also fail (de Chaisemartin 2017 discusses the "tolerating defiers" case). The classic hazard appears when the instrument itself has a heterogeneous direction: an encouragement design that is an incentive for some clients but an annoyance for others (frequently pestered clients simply churn and refuse to adopt) manufactures defiers. Judging monotonicity rests on a substantive understanding of the instrument's mechanism, not on a statistical test.

LATE's extrapolation limitation is the property most in need of an honest reckoning, and it is the through-line of this chapter. IV estimates the effect on the compliers, and the compliers are in turn determined by the instrument. At Meridian LATE 0.784 is far below the ATE 2.018, because what randomized encouragement pries loose are the least-gaining fence-sitters; if management took this 0.784 to argue "the dashboard is useless and not worth promoting," it would be badly wrong, because the firms that should adopt most and gain most (the always-takers) are not in this number at all. Conversely, if the policy question is exactly "should we give the hesitant fence-sitters one more push," then LATE is precisely the most fitting answer. Before using LATE, ask clearly: is the population your policy targets the same batch as the compliers this instrument pries loose? Characterizing the compliers' observable features (like the Abadie $\kappa$ of Section 6.4) and seeing, within the MTE framework, which segment of the resistance spectrum LATE falls on, are both means of easing this extrapolation anxiety.

The forbidden regression is a purely technical but extremely common trap, warned about earlier: the fitted value of a nonlinear first stage cannot be stuffed straight into the second stage, because projection does not pass through a nonlinear function. The only safe roads are linear 2SLS or an explicit control function. Even with a linear first stage, do not manually run two regressions, because then the second stage's standard errors are wrong, and one should use an integrated IV command.

The choice of estimator under heteroskedasticity and many instruments is also worth mentioning. Under heteroskedasticity 2SLS is still consistent but not optimal, and the effective F is the right ruler for judging instrument strength here; when many-instrument bias is severe, LIML is approximately median-unbiased with better weak-instrument properties (fatter tails being the cost), the Fuller correction can bound its moments, and JIVE removes a unit's own overfitting via leave-one-out, all being more robust choices than 2SLS.

Stringing these failure modes together, the credibility of IV rests ultimately not on the mechanics of 2SLS but on two things: whether the instrument is really both strong and clean (relevance plus independence plus exclusion), and whether you have thought clearly about who the compliers it pries loose are, and whether this local effect is on point for your policy question. When the instrument does not stand, the right reaction is not to switch to a fancier estimator, but to admit that this design cannot answer this question, or to look elsewhere for identification (say using a panel structure to difference out time-invariant unobserved heterogeneity, which this series will take up later when discussing Difference-in-Differences). IV's honest boundary is to lay these two things out clearly, rather than papering over them with a pretty first-stage F.

## 8 Further Reading

::: {.readings}
Required reading, in suggested reading order:

- Angrist and Pischke (2009, *Mostly Harmless Econometrics*), Chapter 4. The most readable introduction to IV and LATE, the source of the language of complier, as-good-as-randomly-assigned, and exclusion restriction, and the tone of this chapter is modeled on it.
- Imbens and Angrist (1994, *Econometrica*). The original paper on the LATE theorem, read for how the four assumptions jointly lock the Wald ratio onto the complier effect.
- Angrist, Imbens and Rubin (1996, *JASA*). Recasts IV in the language of potential outcomes, the source of the never-taker / always-taker / complier / defier classification, to be read together with the discussion of exclusion and monotonicity.
- Heckman and Vytlacil (2005, *Econometrica*). The magnum opus establishing MTE as the unifying master object, and the weighting table for ATE/ATT/LATE/OLS is right here, the full version of Section 3.6 of this chapter.

Further reading:

- Vytlacil (2002, *Econometrica*). Proves that the LATE assumptions are equivalent to the nonparametric latent-index selection model, and reading it makes clear that LATE and the Heckman selection model are two sides of one coin.
- Montiel Olea and Pflueger (2013, *JBES*). The effective F and robust weak-instrument testing, the modern practice that replaces "F greater than 10" under heteroskedasticity/clustering.
- Andrews, Stock and Sun (2019, *Annual Review of Economics*). A practical survey of weak instruments, covering the first-stage F, Anderson-Rubin, and conditional likelihood ratio in one sweep, a self-check list before writing an IV paper.
- Angrist and Krueger (1991, *QJE*) and Bound, Jaeger and Baker (1995, *JASA*). The quarter-of-birth instrument and its critique, the most classic cause célèbre of weak and many-instrument bias, and the two must be read together.
- Carneiro, Heckman and Vytlacil (2011, *AER*). The benchmark empirical work on MTE, estimating the marginal return to college education with the NLSY79, demonstrating how to draw the entire curve.
:::

::: {.apa-refs}
- Abadie, A. (2003). Semiparametric instrumental variable estimation of treatment response models. *Journal of Econometrics, 113*(2), 231-263. https://doi.org/10.1016/S0304-4076(02)00201-4
- Anderson, T. W., & Rubin, H. (1949). Estimation of the parameters of a single equation in a complete system of stochastic equations. *The Annals of Mathematical Statistics, 20*(1), 46-63. https://doi.org/10.1214/aoms/1177730090
- Andrews, I., Stock, J. H., & Sun, L. (2019). Weak instruments in instrumental variables regression: Theory and practice. *Annual Review of Economics, 11*(1), 727-753. https://doi.org/10.1146/annurev-economics-080218-025643
- Angrist, J. D. (1990). Lifetime earnings and the Vietnam era draft lottery: Evidence from Social Security administrative records. *American Economic Review, 80*(3), 313-336.
- Angrist, J. D., Imbens, G. W., & Krueger, A. B. (1999). Jackknife instrumental variables estimation. *Journal of Applied Econometrics, 14*(1), 57-67. <https://doi.org/10.1002/(SICI)1099-1255(199901/02)14:1%3C57::AID-JAE501%3E3.0.CO;2-G>
- Angrist, J. D., Imbens, G. W., & Rubin, D. B. (1996). Identification of causal effects using instrumental variables. *Journal of the American Statistical Association, 91*(434), 444-455. https://doi.org/10.1080/01621459.1996.10476902
- Angrist, J. D., & Krueger, A. B. (1991). Does compulsory school attendance affect schooling and earnings? *The Quarterly Journal of Economics, 106*(4), 979-1014. https://doi.org/10.2307/2937954
- Angrist, J. D., & Pischke, J.-S. (2009). *Mostly harmless econometrics: An empiricist's companion*. Princeton University Press. https://doi.org/10.1515/9781400829828
- Bapna, R., & Umyarov, A. (2015). Do your online friends make you pay? A randomized field experiment on peer influence in online social networks. *Management Science, 61*(8), 1902-1920. https://doi.org/10.1287/mnsc.2014.2081
- Björklund, A., & Moffitt, R. (1987). The estimation of wage gains and welfare gains in self-selection models. *The Review of Economics and Statistics, 69*(1), 42-49. https://doi.org/10.2307/1937899
- Bound, J., Jaeger, D. A., & Baker, R. M. (1995). Problems with instrumental variables estimation when the correlation between the instruments and the endogenous explanatory variable is weak. *Journal of the American Statistical Association, 90*(430), 443-450. https://doi.org/10.1080/01621459.1995.10476536
- Carneiro, P., Heckman, J. J., & Vytlacil, E. J. (2011). Estimating marginal returns to education. *American Economic Review, 101*(6), 2754-2781. https://doi.org/10.1257/aer.101.6.2754
- de Chaisemartin, C. (2017). Tolerating defiance? Local average treatment effects without monotonicity. *Quantitative Economics, 8*(2), 367-396. https://doi.org/10.3982/QE601
- Fuller, W. A. (1977). Some properties of a modification of the limited information estimator. *Econometrica, 45*(4), 939-953. https://doi.org/10.2307/1912683
- Hausman, J. A. (1978). Specification tests in econometrics. *Econometrica, 46*(6), 1251-1271. https://doi.org/10.2307/1913827
- Heckman, J. J. (1997). Instrumental variables: A study of implicit behavioral assumptions used in making program evaluations. *The Journal of Human Resources, 32*(3), 441-462. https://doi.org/10.2307/146178
- Heckman, J. J., & Vytlacil, E. (1999). Local instrumental variables and latent variable models for identifying and bounding treatment effects. *Proceedings of the National Academy of Sciences, 96*(8), 4730-4734. https://doi.org/10.1073/pnas.96.8.4730
- Heckman, J. J., & Vytlacil, E. (2005). Structural equations, treatment effects, and econometric policy evaluation. *Econometrica, 73*(3), 669-738. https://doi.org/10.1111/j.1468-0262.2005.00594.x
- Imbens, G. W., & Angrist, J. D. (1994). Identification and estimation of local average treatment effects. *Econometrica, 62*(2), 467-475. https://doi.org/10.2307/2951620
- Lee, D. S., McCrary, J., Moreira, M. J., & Porter, J. (2022). Valid t-ratio inference for IV. *American Economic Review, 112*(10), 3260-3290. https://doi.org/10.1257/aer.20211063
- Montiel Olea, J. L., & Pflueger, C. (2013). A robust test for weak instruments. *Journal of Business & Economic Statistics, 31*(3), 358-369. https://doi.org/10.1080/00401706.2013.806694
- Sargan, J. D. (1958). The estimation of economic relationships using instrumental variables. *Econometrica, 26*(3), 393-415. https://doi.org/10.2307/1907619
- Staiger, D., & Stock, J. H. (1997). Instrumental variables regression with weak instruments. *Econometrica, 65*(3), 557-586. https://doi.org/10.2307/2171753
- Vytlacil, E. (2002). Independence, monotonicity, and latent index models: An equivalence result. *Econometrica, 70*(1), 331-341. https://doi.org/10.1111/1468-0262.00277
- Wald, A. (1940). The fitting of straight lines if both variables are subject to error. *The Annals of Mathematical Statistics, 11*(3), 284-300. https://doi.org/10.1214/aoms/1177731868
- Working, E. J. (1927). What do statistical "demand curves" show? *The Quarterly Journal of Economics, 41*(2), 212-235. https://doi.org/10.2307/1883501
:::
