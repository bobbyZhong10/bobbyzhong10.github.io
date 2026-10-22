---
title: "AI, Algorithms, Technology Adoption & Data Economics"
subtitle: "Algorithms as Agents, Machine Predictions as Data, and the Economics of Data"
seriesline: "Foundations of Information Systems Economics · Chapter 34"
series: "Foundations of Information Systems Economics: Empirical Methods for Information Systems Economics"
sidebartitle: "Foundations of Information Systems Economics"
sidebarsub: "Theory · Methods · Applications"
chapterlabel: "Chapter 34 · AI, Algorithms, Technology Adoption & Data Economics"
---

## Introduction

A call center gives all of its staff access to generative AI, and three months later average handling speed has risen. The result makes it sound as if technology adoption is complete, but it leaves at least four questions open: who actually used the tool, who redesigned their workflow, whether the speed gain came at the cost of error rates or long-run skill, and whether colleagues who were not granted access also benefited from shared prompts. AI exposure is not AI use, use is not workflow redesign, and a task-level improvement is even further from the automation of an entire occupation. If all these layers are lumped into a single "AI treatment," the regression may look clean, but its economic meaning is murky.

AI enters the research scene in two other ways as well. Firms hand pricing, recommendation, and bidding over to learning algorithms, and the repeated interaction among those algorithms can produce equilibrium consequences the designers never explicitly ordered. Researchers, for their part, use large language models to turn text, images, and logs into variables such as sentiment, quality, or intent. The latter looks like it merely saves annotation cost, but it brings prediction error into the estimation: a machine label carries different consequences depending on whether it sits in the dependent variable or the regressors, and a pretty out-of-sample accuracy is no guarantee that the downstream parameters are unbiased. Meanwhile, training data has non-rivalry, externalities, and increasing returns to scale, so you cannot simply apply the intuitions built for traditional capital goods.

This chapter places these problems inside a single economic picture. The technology-adoption part distinguishes access, actual use, and organizational complementary investment, uses factorial assignment to identify the complementarity between AI and training or process redesign, and discusses the Productivity J-curve, task reallocation, augmentation, delegation, and the possibility that immediate performance and long-run deskilling may point in opposite directions. The algorithm part examines learning and potential coordination in repeated pricing, but insists on keeping collusion in simulation separate from evidence in the real world. The measurement part picks up from Chapter 17, handling machine-generated variables with representative gold-standard annotation, error correction, and prediction-powered inference.

What this means is that "used artificial intelligence" is far from a sufficient definition of a treatment, and single-shot assisted performance is far from a sufficient outcome variable. Credible research needs to state, at the same time, on which margin adoption happened and what complementary changes the organization made, and to put assisted-time performance, retention after the tool is taken away, error identification, and skill transfer all side by side; and if the key variable is model-generated, it must keep asking how measurement error enters the estimation and the inference. What AI research is least short of is surprising demonstrations; what is genuinely scarce is a design that can say which technology, through which organizational mechanism, for whom, and over how long, produced what economic consequence.

## 1 A number that does not add up

::: {.case}
Imagine a researcher studying the effect of product-review sentiment on sales. The true sentiment $x^\ast$ is a latent variable that cannot be observed directly, so the researcher uses a large language model to score the sentiment of each product's reviews, obtains a measure $\hat{x}$, then regresses sales $y$ on $\hat{x}$ and reads off the coefficient. We use simulated data to tell this story, and the advantage is that the data-generating process is fully known, so every estimator's performance can be reconciled against the truth, a teaching condition that real data can never give you. The setup is as follows: 4000 products, and the true model is $y = 1 + 2\,x^\ast + \varepsilon$, so the true effect of sentiment on sales is $\beta = 2.0$. The LLM's score is noisy, $\hat{x} = x^\ast + u$, where $u$ is the LLM's measurement error, whose variance makes the reliability of this measure exactly 0.5, meaning that half of the variance of $\hat{x}$ is true signal and half is noise.

Because this is a simulation, we know the answer is $\beta = 2.0$.
:::

The researcher runs the most direct regression: OLS of sales on the LLM sentiment score. The result is $\hat\beta_{\text{naive}} = 0.998$ (SE 0.023). Taken in isolation, this is a pretty result, a precise coefficient, a tiny standard error, highly significant, ready to go into the paper, with the conclusion that "a one-unit increase in sentiment raises sales by about 1 unit."

The problem is that we know the truth is 2.0. The naive regression underestimates the effect by a full half: $\hat\beta_{\text{naive}} = 0.998$ is exactly the truth of 2.0 times the reliability of 0.5. This is no coincidence but the classic consequence of measurement error: when a noisy measure is used as a regressor, the coefficient attenuates toward zero, scaled by a factor that is exactly the reliability. Worse still is that pretty standard error of 0.023. It is not miscomputed: it measures the sampling variation of the naive slope around its own target (the truth of 2.0 times the reliability of 0.5). But the researcher wants 2.0, and a narrow interval centered on 0.998 almost never covers 2.0, so the researcher ends up falsely confident in a coefficient that is off by half.

A regression that takes an AI measure as a precise measurement and plugs it straight in falls into two pitfalls at once: the coefficient is biased, and the inference fails. The main thread of this chapter is to explain where these two pitfalls come from, and how to use a small batch of human-verified gold-standard annotations to recover the truth of 2.0 and bring the standard error to an honest level. Section 6 will show that a proper correction pulls the coefficient from 0.998 back to about 2.07, while the standard error changes from the naive 0.023 to 0.125, a figure that also accounts for the uncertainty of the correction step itself.

## 2 The economic model and the estimand

This chapter has four threads, and this section lays out the theory for each and defines the estimand for the empirical main thread. We first discuss the strategic consequences of algorithms as economic agents, then the economics of data as a factor, next AI adoption and organizational complements, and finally arrive at the estimand for machine-generated variables.

### 2.1 Algorithmic pricing and collusion

Handing pricing over to algorithms changes the nature of competition. A recurring worry is that a few independently operating pricing algorithms, each merely maximizing its own profit, may over long interaction learn on their own to sustain high prices, reaching a kind of collusion with no explicit agreement whatsoever. This worry is not fanciful; Calvano, Calzolari, Denicolò and Pastorello (2020) give powerful evidence for it through simulation.

The setup is a repeated Bertrand oligopoly game with logit demand. There are $n$ firms, firm $i$ sets price $p_i$, and demand is

$$q_i = \frac{\exp((a_i - p_i)/\mu)}{\sum_{j=1}^{n} \exp((a_j - p_j)/\mu) + \exp(a_0/\mu)}$$

where $a_i$ is the product quality index, $a_0$ is the outside option, and $\mu$ is the degree of horizontal differentiation. Per-period profit is $\pi_i = (p_i - c_i)\, q_i$. Given this demand, you can numerically compute two benchmarks: the one-shot Bertrand-Nash price $p^N$ (the competitive benchmark), and the joint-profit-maximizing monopoly price $p^M$ (the full-collusion benchmark).

Each firm is a Q-learning algorithm. It takes the previous period's prices of all firms as its state $s$, maintains a Q-table $Q(s, a)$ recording the long-run value of setting price $a$ in state $s$, and updates according to the following rule:

$$Q(s_t, a_t) \leftarrow (1 - \alpha)\, Q(s_t, a_t) + \alpha\Big[ r_t + \delta\, \max_{a'} Q(s_{t+1}, a') \Big]$$

where $\alpha$ is the learning rate, $\delta$ is the discount factor, $r_t$ is the realized profit in the current period, and $\max_{a'} Q(s_{t+1}, a')$ is the bootstrapped estimate of future value. The algorithm explores with epsilon-greedy: with probability $1 - \epsilon_t$ it picks the price with the current highest Q-value, and with probability $\epsilon_t$ it randomly tries a price, with $\epsilon_t$ decaying over time, shifting gradually from more exploration toward more exploitation.

The key finding is this: two such algorithms, which do not communicate, do not share a goal, and are not written with any collusion instruction, nonetheless converge after long interaction to prices clearly above the competitive level. The degree of collusion is captured by an index,

$$\Delta = \frac{\bar{p} - p^N}{p^M - p^N}$$

where $\bar{p}$ is the average price after convergence, $\Delta = 0$ corresponds to competition, and $\Delta = 1$ corresponds to full monopoly collusion. This index can be defined on prices as above, or on profits (replacing $\bar{p}, p^N, p^M$ with the corresponding profits); Calvano and coauthors use the profit version, while the simulation in Section 6 of this chapter uses the price version, and the two are of comparable magnitude but not exactly equal in value. Calvano and coauthors find that under the duopoly benchmark, the profit-based collusion index captures roughly 70% to 90% of the way from competition to monopoly. The mechanism supporting this result is a reward-punishment strategy the algorithms learn on their own: once one side cuts price and deviates, the other enters a limited-duration price-war punishment phase, then gradually returns to high prices. This is a version, self-emerging from reinforcement learning, of the trigger strategy in the repeated-game folk theorem, reached with no communication at all, purely through each side's own trial-and-error learning. Section 6 will reproduce this phenomenon with a simplified simulation.

The policy puzzle raised by algorithmic collusion is profound. Antitrust law (the first section of the US Sherman Act, Article 101 in the EU) usually requires an "agreement" or "concerted practice," whereas the algorithms reach a collusive outcome without any meeting of the minds, so their conduct may fall outside the reach of existing prohibitions. Klein (2021) shows in simulations that under sequential pricing, algorithms converge to supra-competitive price cycles. Assad, Clark, Ershov and Xu (2024) find field evidence in the German retail gasoline market: the adoption of algorithmic pricing software raised markups, and noticeably so only when all stations in a market adopted it, consistent with the interpretation that the algorithms soften competition. Others caution against overgeneralizing: Miklós-Thal and Tucker (2019) point out that better AI demand forecasting may instead unravel collusion, because it strengthens the temptation to undercut secretly and grab business in high-demand periods. Whether algorithmic pricing helps or harms competition depends on whether the algorithms learn collusion sustained by reward and punishment or a sharper pursuit of each firm's own profit, and that is itself an empirical question.

### 2.2 The economics of data

AI's capability is built on data, and data as a factor of production has three features that differ from traditional inputs, and together they define the economics of data.

The first is data network effects, which must be distinguished from classic network effects. A classic network effect is when a user's utility rises directly with the number of other users (Chapter 32). A data network effect runs an indirect loop: more users generate more data, more data trains a better product, and a better product in turn attracts more users. Value does not come directly from headcount but is mediated through product quality. Hagiu and Wright (2023) characterize this loop and point out that whether the resulting first-mover data advantage is durable depends on the shape of the learning function by which quality improves with data, the asymmetry in learning ability across firms, the speed of data accumulation and depreciation, and other factors, so it is not automatically winner-take-all.

The second is non-rivalry. Data can be used by countless firms at the same time without being used up, which is entirely unlike a machine or a plot of land that can only be used in one place. Jones and Tonetti (2020) infer from this that, from a social standpoint, widespread use of data is valuable, because non-rivalry implies increasing returns to scale; but firms, out of fear of being disrupted, tend to hoard data, which leads to data being used too narrowly. One of their counterintuitive conclusions is that assigning property rights over data to consumers is often closer to optimal than assigning them to firms, because consumers internalize the privacy cost and therefore decide more reasonably whether to share or withhold.

The third is the privacy externality. One person's data can often predict things about other people related to them, so data creates a cross-individual externality. Acemoglu, Makhdoumi, Malekian and Ozdaglar (2022) show that precisely because your data can be inferred from others' data, once everyone else has shared, protecting your own data becomes nearly useless, so the price of your data is driven down, and the market therefore ends up with over-sharing, that is, "too much data, too cheap." This forms one side of a coin whose other side is the Jones-Tonetti worry about "too little use" from hoarding, and together they characterize how the allocation of data can deviate from the optimum in both directions.

The economics of data is both the background and a constraint for the other three threads of this chapter: the capabilities of algorithmic pricing, of the AI tools that firms adopt, and of AI measurement all grow with data, and these three features of data determine how that capability is distributed, who holds it, and with what efficiency it is allocated.

### 2.3 Technology adoption, organizational complements, and the Productivity J-curve

AI is not just a software license you can buy. To turn it into productivity, a firm usually also has to rewrite the workflow, reallocate decision rights, train employees, clean up the knowledge base, and build new quality control. Most of these complementary investments are intangible, counted as cost in the short run, yet not necessarily recorded in a timely way by the capital account or by standard productivity measures.

Express this with a minimal model. The firm chooses AI capital $A$ and organizational capital $O$, and output is

$$Y=F(K,L,A,O),\qquad \frac{\partial^2 F}{\partial A\,\partial O}>0.$$

A positive cross partial means complementarity: the more thorough the workflow redesign, the higher the marginal product of AI; and the stronger the AI, the higher the return on organizational capital. This is exactly the increasing differences property of Section 3.10 of Chapter 2, and the factorial contrast $\tau_{AO}$ of Section 4.6 corresponds to the $\Delta_{as}$ there. The treatment effect of deploying AI alone is not a stable technological constant but depends on the organizational state $O$. This explains why the same tool may have completely different effects for novice workers, experienced workers, different teams, or different management regimes.

Complementarity also changes adoption. A firm adopts when the expected present value of the profit gain from adoption, $\Delta\pi_{is}$, exceeds the fixed and adjustment costs:

$$D_{it}=\mathbf 1\Big\{\mathbb E_t\sum_{s\geq t}\delta^{s-t}\Delta\pi_{is}(A,O)-C_{it}^{adopt}-C_{it}^{reorg}\geq0\Big\}.$$

Therefore early adopters are not random firms. They may have better managers, more compatible data foundations, and lower reorganization costs, or they may be facing a stronger negative shock. A simple comparison of adopters and non-adopters mixes the treatment effect, selection, anticipation, and dynamic option value.

The Productivity J-curve of Brynjolfsson, Rock and Syverson (2021) adds the time dimension. Early in adoption, the firm invests heavily in unmeasured intangible capital, and measured productivity may first fall; only once the complements are built and the returns begin to be harvested does productivity rise. A short-run negative coefficient in an event study therefore does not necessarily mean the technology is ineffective, and a long-run positive coefficient is not necessarily all a contemporaneous treatment effect. The research design needs to distinguish, in advance, the installation period, the organizational investment period, and the harvesting period.

Technology diffusion also involves both an extensive and an intensive margin. Comin and Hobijn (2010) distinguish the lag until a technology is first adopted from the penetration speed after adoption. A country or firm can own a technology very early yet use it only in a very small scope; using only a binary adoption date will miss intensity, quality, and use. For AI, "buying a copilot," "employees actually invoking it," "writing the suggestions into the production process," and "the organization restructuring around AI" are four different margins.

::: {.theorem data-label="Organizational complementarity"}
If $F$ is twice continuously differentiable and $F_{AO}>0$, the marginal product of AI, $F_A$, is strictly increasing in organizational capital $O$ (proof: $\partial F_A/\partial O = F_{AO}>0$). So the AI effect is not a technological constant, and estimating only the average adoption coefficient will mask this joint production structure; the result applies to workflow redesign, training, data readiness, and decentralized decision rights only insofar as they are components of $O$.
:::

Empirical predictions must stay restrained. If the interaction of AI and training is positive, it supports complementarity, but it may also come from differential compliance, measurement, or the selection of training toward high-potential workers. The most powerful design is factorial randomization, which randomizes AI access and the organizational intervention separately; next best is the interaction of an exogenous rollout with a predetermined training schedule. Grouping after the fact by observed usage usually reintroduces selection.

### 2.4 Exposure is not automation: task, job, and system

The distinction most likely to cause a serious misreading in the AI-and-labor literature is this: exposure measures whether a technology, at its current capability, can significantly change how a task is performed (the task-level exposure concept of Eloundou et al. 2023); it does not directly say the task will be automated, still less that the entire job will disappear. A job is a bundle of interdependent tasks, and a firm is a system built out of tasks, workflows, responsibilities, and risk controls. Even if many tasks are exposed, the organization may choose augmentation, recombine tasks, or create new services rather than replace labor.

Therefore, going from task-level exposure to wages, employment, or productivity requires at least four additional objects: the cost of adoption and the speed of diffusion, the complementarity relations among tasks, the elasticity of product demand, and labor supply and general-equilibrium adjustment. High exposure can raise the marginal productivity of a class of skills, or it can lower the demand for it; the direction cannot be decided by the exposure score alone. Empirical research should measure exposure, actual adoption, task reallocation, and labor-market outcomes separately.

### 2.5 Deskilling, learning, and dynamic human capital

The short-run productivity effect and the long-run skill effect of AI adoption can point in opposite directions. Suppose employee skill accumulates according to

$$h_{t+1}=(1-\rho)h_t+\phi e_t,$$

where $\rho\in(0,1)$ is the depreciation rate of skill, $\phi>0$ is the efficiency with which practice turns into skill, and $e_t$ is the effort that requires active diagnosis, practice, and feedback. When AI does the hard part on the worker's behalf, current output rises but $e_t$ falls, so future skill $h_{t+1}$ declines; if the AI instead provides explanation, feedback, and moderate scaffolding, it may raise current performance and learning at the same time. The direction is determined by work design, not by the binary variable "whether AI is used."

This yields a testable dynamic prediction. Pure delegation should lower performance under the no-AI condition, on a delayed test, in error-detection ability, or on novel-task performance; instructional augmentation should improve transfer. Heterogeneity also matters: senior workers may spend the saved cognitive resources on higher-level tasks, while novices lose the practice that forms a basic mental model. A high-standard experiment should jointly report assisted performance, unassisted retention, error detection, task breadth, and long-run occupational mobility, not just immediate throughput.

### 2.6 Text as data and the estimand for generated variables

Now we settle onto the empirical main thread. As a measurement tool, AI's most typical scenario is text-as-data, turning unstructured text into analyzable variables. Gentzkow, Kelly and Taddy (2019) organize this process into a clear pipeline: first represent the documents as a document-term matrix, that is, the count of each document on each word or phrase (often with tf-idf weighting), an object of extremely high dimension whose vocabulary size often exceeds the number of documents; then use one of three classes of methods to map the high-dimensional text features into the economic variable you want. Dictionary methods weight with a preset word list, supervised methods learn a prediction model from labeled data, and unsupervised methods (such as topic models) extract latent structure from the text. The LLM belongs to a powerful version of the supervised class: it directly outputs a prediction of the target concept.

Whichever class of method you use, the variable $\hat{x}$ that comes out is a generated variable, that is, a noisy quantity produced by a prior estimation or prediction step, not the truth. Once you put it into a downstream regression, you are back to this chapter's estimand problem.

This chapter's estimand is the true structural coefficient $\beta$, the causal or structural effect of the latent true variable $x^\ast$ (true sentiment, true quality, true intent) on the outcome $y$. We cannot observe $x^\ast$, only its machine-generated measure $\hat{x} = x^\ast + u$. The goal is, from the noisy $\hat{x}$ and the outcome $y$, to consistently estimate the true coefficient $\beta$ of $x^\ast$ and to give honest inference. Between this estimand and its measure lies a layer of measurement error, and how to get through this layer is the core of identification in Section 3. It must be stressed that the challenge here is not causal identification (we assume the relationship between $x^\ast$ and $\varepsilon$ is already clean), but a purely measurement and inference challenge, which exists independently of the causal problem, and which happens all the same even in a perfect randomized experiment, as long as the key variable is machine-generated.

To summarize this section: algorithms as economic agents may learn, in repeated interaction, collusion without agreement; data as a factor has a distinctive allocation economics due to network effects, non-rivalry, and privacy externalities; the return to AI adoption depends on the complementarity of workflow, skills, and organizational capital; and AI as a measurement tool produces noisy generated variables. The four threads correspond respectively to market equilibrium, factor allocation, firm production, and empirical measurement, and cannot be replaced by one sweeping "AI effect" coefficient.

## 3 Identification: the two pitfalls of machine-generated variables

Identification logically precedes estimation. The question here is: when the key variable is machine-predicted, under what conditions, and with what information, can we consistently estimate the true coefficient $\beta$ and make valid inference? This section is the most detailed in the chapter, broken into numbered subsections.

### 3.1 Attenuation: measurement error as a regressor

Start with the first pitfall, bias. The true model is $y = \beta x^\ast + \varepsilon$, and we replace $x^\ast$ with the noisy measure $\hat{x} = x^\ast + u$. To get a clean attenuation result, we need an assumption about the measurement error:

::: {.assumption}
**(classical measurement error)** The error $u$ of the machine measure satisfies $\hat{x} = x^\ast + u$, where $\mathbb{E}[u] = 0$, and $u$ is independent of both the true variable $x^\ast$ and the structural error $\varepsilon$. That is, the LLM's score is the truth plus a layer of pure noise unrelated to the truth and to the outcome.
:::

Under this assumption, running OLS of $y$ on $\hat{x}$ has probability limit

$$\mathrm{plim}\; \hat\beta_{\text{OLS}} = \beta \cdot \frac{\sigma_{x^\ast}^2}{\sigma_{x^\ast}^2 + \sigma_u^2} = \beta \cdot \lambda$$

where $\lambda = \sigma_{x^\ast}^2 / (\sigma_{x^\ast}^2 + \sigma_u^2) \in (0, 1)$ is the reliability. The coefficient is multiplied by a factor smaller than 1 and so attenuates toward zero; this factor is exactly the share of the measure's variance that is true signal. The noisier the measure, the smaller $\lambda$, and the more severe the attenuation.

::: {.intuition}
Why is the direction always toward zero? The intuition is this: the noise $u$ mixed into $\hat{x}$ has nothing to do with $y$; it merely adds, on top of the true signal, a layer of jitter unrelated to the outcome. The regression slope measures, when $\hat{x}$ moves by one unit, how much $y$ moves on average. Part of the variation in $\hat{x}$ is pure noise that does not move $y$ at all, so the average movement in $y$ per unit of movement in $\hat{x}$ is diluted. The larger the noise share, the heavier the dilution. This is why running a regression with a fairly noisy LLM score gives an estimated effect that is noticeably too small, not too large. The direction is certain: toward zero.
:::

In this case the reliability is exactly 0.5, so the naive coefficient 0.998 is precisely half of the truth 2.0. Within the above setting of a single regression, classical additive error, and positive reliability, the naive coefficient is attenuated in absolute value; once you add other mismeasured covariates, or the error is differential, or the prediction error is correlated with the truth and the structural error, this directional conclusion no longer holds. Therefore, do not treat "an AI measure pushes the coefficient toward zero" as a general law for machine-generated variables.

There is one important contrast to make clear right away, to avoid misuse. The attenuation above happens when the machine-generated variable is a regressor (on the right-hand side). If the machine generates the outcome variable (the left-hand side), the situation is different: if $\hat{y} = y^\ast + v$ and $v$ is classical error, the OLS slope is actually consistent, only the residual variance grows larger and the standard error widens. But the error of a machine prediction is often not classical; it is correlated with the features (the prediction is systematically off at certain values), and then even on the left-hand side it brings bias. So a practical discipline is: first ask which side of the regression this machine-generated variable is on; on the right-hand side, attenuation is the default; on the left-hand side, look at whether the error is correlated with the covariates.

### 3.2 Overconfidence: inference failure as a generated variable

The second pitfall is inference. The naive standard error 0.023 is not miscomputed: it correctly measures the sampling variation of the naive slope around $\beta\lambda$. The naive interval fails because it is centered on $\beta\lambda$ rather than $\beta$; that is bias, not a standard error that is too small (Section 3.6 of Chapter 5). The inference problem arises after the correction: the correction uses $\hat\lambda$ or a calibration equation estimated on the gold sample, so the corrected regressor becomes a generated regressor in the sense of Pagan (1984). The conventional standard error of the second-step regression pretends that the first-step estimate is exact, leaves out its sampling error, and is therefore generally invalid.

In this case the corrected estimator's standard error is 0.125, more than five times the naive 0.023. The increase has two sources: dividing by $\hat\lambda$ (0.482) roughly doubles the standard error, and the rest comes from the sampling error of $\hat\lambda$, which is estimated from only 400 labels; this second piece is exactly what the conventional second-step standard error leaves out. So 0.125 and 0.023 are standard errors of two different estimators, and it is wrong to say that the naive standard error is five times too small. Correct inference must add the uncertainty of the first step back in. Section 4 gives the specific method, and Section 4.4 of Chapter 8 gives the general correction for two-step standard errors.

### 3.3 Why prediction error does not vanish as the sample grows

Some will wonder whether these problems automatically vanish when the sample is large enough. They do not. This is the point to be most wary of with machine-generated variables. The attenuation from measurement error is a specification problem, not sampling noise; no matter how large the sample, $\mathrm{plim}\,\hat\beta$ is still $\beta\lambda$ rather than $\beta$. Increasing the LLM corpus from ten thousand to a million texts does not shrink the measurement-error variance $\sigma_u^2$ of $\hat{x}$, because it comes from the model's inherent uncertainty in predicting each individual and is unrelated to the sample size. In fact, the larger the sample, the narrower the naive interval, yet it remains centered on $\beta\lambda$: you grow ever more confident in a biased coefficient, so the problem is amplified rather than diluted. The bias and inference failure of machine-generated variables can only be solved by correction, not by big data.

### 3.4 How a small batch of gold standard rescues the truth

Since sample size will not save you, where does the rescue come from? The answer is a small batch of gold-standard annotations, that is, on a small subsample, using human verification to obtain the true variable $x^\ast$ (or the true label), so that we know exactly where and by how much the machine measure is wrong. This small batch of gold data is the key to identification, because it makes the unobservable measurement-error structure observable.

How it works depends on the correction method, but the underlying logic is the same. On the gold subsample we observe both the truth $x^\ast$ and the machine measure $\hat{x}$, so we can estimate their relationship, such as the reliability $\lambda$ or the systematic bias of the machine measure relative to the truth. With this relationship in hand, we can apply it to the entire large sample and calibrate the machine measure back to the truth. The large sample provides efficiency (machine measures in vast numbers are almost free), and the small gold sample provides validity (it corrects the bias and gives honest uncertainty). This is exactly the core idea of modern prediction-powered methods and the main thread of this chapter's estimation section; Sections 3.1 to 3.3 of Chapter 17 set out the same logic from the side of the outcome variable.

::: {.warning}
The gold-standard subsample must be randomly drawn from the same population, and the true labels must genuinely be reliable. If the batch sent for human verification is itself selective (say only the easy-to-judge cases were picked), or the so-called gold labels themselves carry systematic error, the correction will swap one bias for another. The upper bound on the credibility of the correction is set by the quality and representativeness of that small batch of gold data, and this step cannot be skimped on or faked.
:::

### 3.5 Identifying AI adoption and complementarity

The first threat in studying AI adoption is endogenous timing. Managers adopt when they anticipate rising demand, worsening costs, or the start of organizational reform, so pre-treatment trends may already be shifting. Staggered adoption can use the group-time estimands of Chapter 19, but parallel trends must be argued with respect to "how the treated cohort would have evolved had it not yet received access," and insignificant coefficients in the few periods before adoption are no substitute for that argument.

The second threat is confounding among access, use, and redesign. Randomized access identifies the ITT; actual use is selected by skills, workload, and beliefs, so directly comparing users with non-users is no longer random. If access shifts use significantly, you can use access as an instrument to identify the use effect for compliers, but exclusion requires that access not directly affect the outcome through experimentation, manager attention, or coworker spillover.

The third threat is the interaction interpretation of complementarity. In the outcome model

$$Y=\tau_0+\tau_A A+\tau_O O+\tau_{AO}AO+\varepsilon,$$

$\tau_{AO}>0$ identifies average complementarity only when the joint assignment of $(A,O)$ is exogenous, the variable scales have economic meaning, and the outcome model is adequate. A factorial experiment is the most direct; if $O$ is a workflow redesign that the firm chooses on its own after treatment, the interaction is post-treatment selection, not a clean estimate of organizational complementarity.

The fourth threat is interference. AI may let treated workers take over more tasks, leave harder cases to untreated coworkers, or improve a shared knowledge base. A worker-level difference then also contains task reallocation. This needs team-level randomization, a task-level outcome decomposition, or the exposure mapping of Chapter 26.

The fifth threat is productivity measurement. Issues per hour, quality, customer satisfaction, and worker learning may point in different directions. When AI raises speed but lowers long-run skill accumulation, short-run flow productivity does not equal the long-run human-capital effect. At a minimum, jointly report quantity, quality, escalation, retention, and learning outcomes.

The main points of this section can be summarized as follows: under classical measurement error, a machine-generated variable on the right-hand side attenuates the coefficient toward zero, and the bias does not vanish as the sample grows; the gold standard makes the measurement-error structure estimable. Causal identification of AI adoption must additionally deal with endogenous timing, the access-use gap, post-treatment organizational choice, interference, and multidimensional productivity. Measurement correction and adoption identification are two separate hurdles.

## 4 Estimation: from naive regression to prediction-powered correction

This section proceeds by "where the previous method fails gives rise to the next." The goal is to consistently estimate the true coefficient $\beta$ and give an honest standard error.

### 4.1 Naive OLS

The most direct approach is to take the machine measure as the truth and run OLS of $y$ on $\hat{x}$. Its failure was made clear in Section 3: the coefficient attenuates toward $\beta\lambda$, and the interval is centered on the wrong value. In this case it returns 0.998 (SE 0.023), half the truth of 2.0, with an interval far from 2.0. Its only use is as a lower-bound reference for the attenuation.

### 4.2 Gold standard only

The second idea is to simply run the regression on the gold subsample that has the truth $x^\ast$, avoiding measurement error entirely. This is correct on bias, since OLS of $y$ on the true $x^\ast$ is consistent. The cost is precision: the gold sample is small, so the estimation is noisy. In this case, using 400 gold annotations gives 1.984 (SE 0.080), nearly unbiased, but the standard error is several times larger than with the full sample, wastefully throwing away the information in the other several thousand machine measures. This exposes the core tension: machine measures are many and dirty, gold annotations are accurate and few, and the ideal estimator should have both.

### 4.3 Reliability correction (errors-in-variables)

If you can estimate the reliability $\lambda$, the attenuation can be inverted directly. Estimate $\lambda$ from the gold subsample (where both $x^\ast$ and $\hat{x}$ are present, the slope of regressing $x^\ast$ on $\hat{x}$ is one estimate of the reliability), then divide the full-sample naive coefficient by it:

$$\hat\beta_{\text{corrected}} = \frac{\hat\beta_{\text{OLS}}}{\hat\lambda}$$

In this case the gold sample gives $\hat\lambda = 0.482$ (truth 0.500), and the corrected coefficient is about 2.068 (computed from the unrounded value), almost exactly the truth of 2.0. The standard error requires a bootstrap that also captures the uncertainty of the generating step; it comes to 0.125, honestly wider than the naive 0.023. But it is also wider than the gold-only 0.080: $\hat\lambda$ is estimated from only 400 labels, its error passes through the division into the coefficient, and the ratio estimator does not use the $y$ in the gold sample. So in this design the reliability correction is not more precise than the gold sample alone; its usefulness lies in settings where the gold sample verifies $x^\ast$ but contains no $y$. This correction is equivalent to first calibrating the machine measure back to the truth with the gold sample and then running the regression, a procedure known as regression calibration; the two routes give the same number.

### 4.4 SIMEX and modern prediction-powered methods

The reliability correction requires the error to be classical and $\sigma_u^2$ to be estimable. When the attenuation has no analytic inverse (in a nonlinear model, say) or the error structure is more complex, there are two classes of more general tools.

SIMEX (simulation-extrapolation) does not solve for the inverse analytically but relies on simulation. It deliberately adds noise of varying intensity to the measure, observes how the coefficient attenuates further as noise grows, fits this attenuation curve, then extrapolates to the "zero noise" point to obtain the corrected coefficient. Its advantage is generality; the cost is that it still assumes the error is additive and unrelated to the truth and the outcome (the assumption of Section 3.1), and its variance must be known or estimated.

Better suited to the AI setting is the class of methods known as prediction-powered inference. They do not require the machine model to be correctly specified or unbiased, but they still impose requirements on the target population, the sampling mechanism of the gold sample, label quality, and the relevant moment conditions. The data structure is a small batch of gold annotations plus a large batch of machine predictions. The PPI of Angelopoulos, Bates, Fannjiang, Jordan and Zrnic (2023) is the leading example, and Section 3.3 of Chapter 17 has already set out its construction and variance from the side of the outcome variable. The idea is first to compute an estimator on the vast set of machine predictions and then to correct its bias using the "difference between machine prediction and truth" on the gold sample; this correction term is called the rectifier. For the mean, for example, with $f$ the machine prediction model, a gold sample $(X_i, Y_i)$ of size $n$ in which $Y_i$ is the truth, and a sample $\tilde{X}_j$ of size $N$ that has machine predictions only,

$$\widehat\theta_{PPI} = \underbrace{\frac{1}{N}\sum_{j} f(\tilde{X}_j)}_{\text{estimate on predictions}} - \underbrace{\frac{1}{n}\sum_{i} \big( f(X_i) - Y_i \big)}_{\text{rectifier}}$$

The first term is the estimate computed on the vast set of machine predictions, where the large sample provides efficiency. The second term is the rectifier: it uses the small gold sample to estimate the systematic bias of the machine predictions and subtracts it, which yields an unbiased estimate and a valid confidence interval without trusting the machine model. The same logic underlies the post-prediction inference of Wang, McCormick and Leek (2020) for machine-predicted outcome variables and the design-based supervised learning (DSL) of Egami, Hinck, Stewart and Wei (2023) for LLM annotation; both use a small batch of expert labels to correct the bias of a large batch of machine labels. Battaglia, Christensen, Hansen and Sacher (2024) go further and combine the upstream variable generation and the downstream regression into one joint estimation step, handling this two-step problem at its root. These methods share one structure: the large prediction set gives efficiency, and the small gold set gives validity. The reliability correction of Section 4.3 borrows only this idea and is not a special case of PPI: it builds no rectifier, uses the gold sample only to estimate $\lambda$, needs the classical-error assumption for validity, and need not be more precise than the gold sample alone.

### 4.5 Estimating algorithmic collusion

The second thread of this chapter, algorithmic collusion, has two forms of "estimation." One is simulation, as in Section 2.1, directly dropping Q-learning algorithms into a specified game and running it, observing the prices and collusion index they converge to; this is Calvano and coauthors' approach and what Section 6 will reproduce. It answers "given such algorithms and such a market, will collusion emerge on its own," and is a computational experiment rather than an estimation on real data. The other is field identification, as in Assad and coauthors, estimating the causal effect of adopting algorithmic pricing in real market data, where the difficulty is that the timing of adoption is unobservable and endogenous; their response is to use structural breaks in pricing behavior to identify adopters, then use headquarters-level adoption as an instrument to handle the endogeneity. The two forms answer different questions: the simulation probes the possibility of the mechanism, and field identification measures the actual consequences in reality.

### 4.6 Estimating adoption, diffusion, and organizational complementarity

The estimator for the adoption effect depends on the assignment. Randomized access first reports the ITT, and retains the assignment probability and the rollout strata. If actual use is of interest, report the first stage, the complier population, and IV/LATE, rather than plugging the fitted use probability straight into a nonlinear outcome model.

Diffusion reports at least two margins. The extensive margin is the first-adoption hazard, which can use a duration model or an event-history design; the intensive margin is active users, tasks covered, API calls, or workflow share. Binary adoption treats a firm that trials once and a firm that fully restructures its process as the same, systematically flattening heterogeneity.

Organizational complementarity is best suited to a $2\times2$ factorial design: randomize AI access $A$ and workflow training $O$, giving the four cells $(0,0),(1,0),(0,1),(1,1)$. The difference-in-differences in experimental cells

$$\tau_{AO}=\big(\bar Y_{11}-\bar Y_{01}\big)-\big(\bar Y_{10}-\bar Y_{00}\big)$$

directly estimates how much the AI marginal effect increases due to training. If team interference exists, both randomization and inference move up to the team level. If the outcome is on a nonlinear scale, the economic interpretation of the interaction depends on the scale, and you should report both cell means and policy-relevant contrasts.

The J-curve needs dynamic estimands. Defining event time zero as access, installation, or effective use will change the path. A reasonable approach is to record the rollout, training, workflow migration, and active-use dates at the same time, predefine the investment window, and report the cumulative effect, rather than cherry-picking one post coefficient.

The three main generated-variable correction estimators are compared and summarized in the table below.

| Dimension | Naive OLS | reliability / EIV correction | prediction-powered (PPI/DSL) |
|---|---|---|---|
| Bias | attenuates toward zero $\beta\lambda$ | consistent after correction | consistent after correction |
| Standard error | valid only for $\beta\lambda$ (overconfident) | bootstrap-corrected | valid inference built in |
| What is needed | machine measure only | gold sample to estimate $\lambda$, error must be classical | gold sample, any error structure |
| Assumption on the machine model | treated as truth | classical additive error | no assumption |
| Result in this case | 0.998 (SE 0.023) | 2.068 (SE 0.125) | same class of idea, see text |

The roadmap of this section is now complete: naive OLS attenuates and is overconfident because of measurement error; using the gold sample only is unbiased but wastes information; the reliability correction inverts the attenuation under classical error and bootstrap-corrects the standard error; SIMEX corrects the attenuation by simulation and extrapolation without an analytic inverse, but still needs the error variance to be known or estimated; the PPI class of methods makes no assumption about the error structure and secures validity from the small gold set, but basic PPI is more precise than the gold sample alone only when the predictions are accurate enough, and guaranteeing that it does no worse than the gold sample alone requires PPI++ (Section 3.3 of Chapter 17); and algorithmic collusion is studied in two estimation forms, simulation and field identification.

## 5 Anchor papers

Methods prove themselves only in real research. This section anchors the chapter in four papers: one that is the simulation benchmark for algorithmic collusion, one that is the methodological benchmark for correcting inference from machine predictions, one that is a field study of generative AI's effect on real knowledge-work productivity, and one that is the mechanism explanation of how organizational complementary investment makes returns follow a Productivity J-curve. Each is laid out along the five elements of paper, method, data, results, and limitations, with a focus on how the assumptions are defended.

### 5.1 Calvano, Calzolari, Denicolò and Pastorello (2020, American Economic Review)

::: {.case}
Paper and place in the history of methods: "Artificial Intelligence, Algorithmic Pricing, and Collusion," American Economic Review 110(10), 3267-3297. This paper turned "will algorithms collude on their own" from a conjecture into reproducible computational evidence, and is the benchmark of the algorithmic-pricing literature.

Method: let two (or more) independent Q-learning pricing algorithms interact repeatedly in a repeated Bertrand-logit game, each algorithm knowing only its own profit, updating by the Q-learning rule, exploring with epsilon-greedy, not communicating, not sharing a goal, and containing no collusion instruction whatsoever. After the algorithms converge, the collusion index $\Delta = (\bar\pi - \pi^N)/(\pi^M - \pi^N)$ measures where the outcome lands between competition and monopoly, and the response pattern of the algorithms upon deviation is examined.

Data: not real market data, but large-scale simulation. The authors run repeatedly over a grid of learning rates and exploration rates, recording the converged prices, the collusion index, and the dynamic path after a deviation, to confirm the robustness of the result rather than a single stroke of luck.

Result: the two algorithms consistently converge to supra-competitive prices, and under the duopoly benchmark the collusion index is roughly 0.7 to 0.9, capturing most of the way from competition to monopoly. More crucial is the mechanism. The authors' simulations show that the algorithms learn a reward-punishment strategy on their own: a deviating price cut by one side triggers a limited-duration price-war punishment, followed by a gradual return to high prices. This is a classic device for sustaining collusion in repeated games, yet here it emerges spontaneously, with no one designing it. Collusion weakens as the number of firms grows but remains significant.

Limitation: the conclusions come from a specific algorithm (tabular Q-learning) and a specific game setup, the Q-table expands exponentially with the state space, and real-world pricing algorithms are not necessarily this kind. Convergence requires an extremely long learning period, and real markets do not necessarily give algorithms so much room for trial and error. Most fundamental is external validity: the simulation proves collusion can arise in algorithmic interaction, but this does not mean algorithmic pricing in real markets is definitely colluding, which field evidence must answer.
:::

This paper's defense strategy is to place a phenomenon hard to observe cleanly in reality into a fully controllable simulation, use a large number of repeated runs to rule out chance, and use the examination of the deviation dynamics to nail down the mechanism. Its persuasiveness lies in the clarity of the mechanism, and its limit lies exactly there too: the simulation can prove possibility, but it cannot replace measurement of a real market.

### 5.2 Angelopoulos, Bates, Fannjiang, Jordan and Zrnic (2023, Science)

::: {.case}
Paper: "Prediction-powered inference," Science 382(6671), 669-674. The question is how, as research increasingly relies on machine-learning predictions as data, one can still obtain valid estimates and confidence intervals without trusting the prediction model; this is exactly the problem that this chapter's empirical main thread sets out to solve. Section 5.3 of Chapter 17 uses the same paper as an anchor from the side where the machine label is the outcome variable.

Method: PPI. The data set consists of two parts: a small gold sample with both true labels and machine predictions, and a large sample with machine predictions only. The core is correction with a rectifier: first compute the statistic you want (mean, quantile, regression coefficient) on the vast set of predictions, then subtract the average bias, "prediction minus truth," measured on the gold sample, and finally construct the confidence interval from the sum of two independent variances. The whole method makes no assumption about the machine model.

Data: the method itself is general, and the paper demonstrates it in several fields, including scenarios where machine-predicted protein structures, galaxy properties, and gene expression serve as data, showing how the same correction gives valid inference on different machine predictions.

Result: the confidence intervals from PPI have valid coverage and, when the predictions are accurate enough, are narrower than the classical interval based on the gold sample alone; without sacrificing validity, PPI turns the vast machine predictions into precision. It moves the answer to the question "can machine predictions be used as data?" from "using them directly fails" to "using them with this correction is valid."

Limitation: validity depends on the gold sample being randomly drawn from the same population, with reliable true labels. If the gold sample is too small, the correction term is itself noisy and the precision gain is limited. If the machine predictions are nearly unrelated to the truth, basic PPI does not fall back to the gold-sample-only estimate: the extra prediction term only adds variance, and the PPI interval can be wider than the classical one. Guaranteeing that PPI does no worse than the classical interval requires PPI++ with power tuning (Section 3.3 of Chapter 17). PPI addresses the validity of inference, not causal identification; if the downstream regression itself suffers from endogeneity, PPI cannot fix that for you.
:::

### 5.3 Brynjolfsson, Li and Raymond (2025, Quarterly Journal of Economics)

::: {.case}
Paper: "Generative AI at Work," Quarterly Journal of Economics 140(2), 889-942. It is the representative field study of generative AI's effect on real knowledge-work productivity.

Method: studies a Fortune 500 business-process software firm deploying a conversational assistant in waves, using a staggered event-study design with agent and time fixed effects, and analyzes heterogeneity by experience, skill, and stage of use.

Data: high-frequency work records of 5,179 customer-support agents, with outcomes including issues resolved per hour, handle time, quality, and worker retention.

Result: AI access raises productivity (issues resolved per hour) by about $14\%$ on average, with the gains highly concentrated among novice and lower-skill workers (about $34\%$); the model helps them adopt the communication patterns of high-performing workers, exhibiting knowledge diffusion and skill compression.

Limitation: the rollout was not fully randomized for academic research, and identification depends on the rollout timing and untreated trends; the study comes from one firm and one support task, and cannot be automatically extrapolated to all occupations. The access effect, active use, and long-run skill accumulation still need to be distinguished.
:::

### 5.4 Brynjolfsson, Rock and Syverson (2021, AEJ: Macroeconomics)

::: {.case}
Paper: "The Productivity J-Curve: How Intangibles Complement General Purpose Technologies," American Economic Journal: Macroeconomics 13(1), 333-372.

Method: builds a dynamic model of general-purpose technology and intangible complementary capital, showing how the investment period and the harvesting period make measured productivity first low then high.

Data: the paper calibrates against, and illustrates the mechanism with, the productivity patterns of historical GPTs and digital-technology investment, and is not a single-firm treatment-effect study.

Result: the large organizational investment early in a new technology is treated as a current cost yet not counted as capital, leading to productivity being understated; later, when these complements are harvested, productivity may be overstated. The J-curve explains the timing mismatch between adoption and measured returns.

Limitation: the J-curve is one structural explanation, not uniquely identified simply from seeing a fall-then-rise. Demand recovery, selection, survivorship, and concurrent reorganization can all produce a similar path, so additional variation is needed to distinguish them.
:::

Taken together, the four papers make the meaning of anchoring more complete: Calvano and coauthors show AI as an economic agent changing the competitive equilibrium, Angelopoulos and coauthors show how AI as a measurement tool is corrected, Brynjolfsson, Li and Raymond show the heterogeneous effect of AI access on real-work productivity, and Brynjolfsson, Rock and Syverson explain why organizational complementary investment makes returns follow a J-curve over time. Market strategy, measurement, adoption, and organizational transformation are four different research objects, and cannot be replaced by one sweeping "AI works" verdict.

## 6 A full walkthrough

Now we run three of these threads (machine-generated variables, algorithmic collusion, and AI adoption) in full, each on its own simulated data set. The following code uses R 4.6.1 and fixes the random seeds to ensure reproducibility (Section 6.1 uses set.seed(22); the sessions of Section 6.2 use seeds 1 to 6; Section 6.3 uses set.seed(2604)), and every number cited in the text comes from the actual run output of this code.

### 6.1 Generated regressor: from 0.998 to 2.07

The design parameters are as follows: 4000 products, the true model $y = 1 + 2\,x^\ast + \varepsilon$, the LLM measure $\hat{x} = x^\ast + u$, with measurement error that makes the reliability 0.5; 400 of the products are human-verified, which yields the gold standard $x^\ast$.

```r
set.seed(22)
N <- 4000; n_L <- 400; b1 <- 2.0
x_star <- rnorm(N, 0, 1)
y      <- 1.0 + b1 * x_star + rnorm(N, 0, 1.5)
x_hat  <- x_star + rnorm(N, 0, 1.0)      # LLM score = truth + noise (reliability 0.5)
```

Naive OLS, taking the LLM measure as the truth:

```r
m_naive <- lm(y ~ x_hat)
```

gives $\hat\beta_{\text{naive}} = 0.998$ (SE 0.023), exactly the truth 2.0 times the reliability 0.5, which is the attenuated number from Section 1. Using only the 400 gold annotations:

```r
lab <- sample(N, n_L)
m_gold <- lm(y[lab] ~ x_star[lab])
```

gives 1.984 (SE 0.080), nearly unbiased but with poor precision. The reliability correction, first estimating $\lambda$ from the gold sample, then inverting:

```r
reliab_hat <- coef(lm(x_star[lab] ~ x_hat[lab]))[2]          # 0.482
b_eiv <- coef(m_naive)["x_hat"] / reliab_hat                 # 2.068
```

The gold sample gives $\hat\lambda = 0.482$ (truth 0.500), and the corrected coefficient is 2.068, almost exactly on the truth. The standard error uses a bootstrap to fold in the uncertainty of the generating step, giving 0.125.

![The three estimates of the sentiment effect $\beta$ against the truth (dashed line, 2.0). Naive OLS takes the LLM measure as the truth, so the coefficient attenuates to about half, with an extremely narrow standard error; the gold-sample-only estimate is unbiased, with a wider interval than the naive one; the reliability correction is almost exactly on the truth and has the widest interval, because it also accounts for the estimation error in $\hat\lambda$.](assets/fig/fig_34_genreg.svg)

The figure above plots the three estimates together with their confidence intervals. The naive estimate's point is far to the left of the truth with a suspiciously narrow interval, the gold-sample estimate is on target with a wider interval, and the corrected estimate is on target with the widest interval. Together the three show this chapter's two pitfalls and its one remedy at a glance: attenuation, overconfidence, and calibration with gold annotations.

### 6.2 Algorithmic collusion: two Q-learners learn to raise prices

The second thread drops two Q-learning pricing algorithms into a Bertrand-logit game. The demand parameters are $a = 2$, $a_0 = 0$, $\mu = 0.25$, cost $c = 1$, and we first numerically compute the competitive and monopoly benchmarks:

```r
pNash <- 1.473          # Bertrand-Nash price
pMon  <- 1.925          # monopoly price
```

The price is discretized into 11 grid points, covering from slightly below the competitive price to slightly above the monopoly price. Each algorithm maintains a Q-table, updates by the Q-learning rule, with the epsilon-greedy exploration rate decaying exponentially over time, and runs for one million periods:

```r
Q1[s1, s2, a1] <- (1 - lr) * Q1[s1, s2, a1] +
                  lr * (r1 + delta * max(Q1[a1, a2, ]))   # Q-learning update
```

To rule out a single stroke of luck, six sessions are run independently (with different random seeds). The result: the average prices after convergence in the six sessions are 1.775, 1.791, 1.764, 1.657, 1.764, and 1.829, all clearly above the competitive price 1.473. The price-based collusion index $\Delta = (\bar{p} - p^N)/(p^M - p^N)$ is 0.668, 0.703, 0.645, 0.407, 0.643, and 0.787 across the six sessions, with a mean of 0.642, meaning it captures on average about 64% of the way from competition to monopoly, and all six sessions are supra-competitive.

![The average price of two Q-learning algorithms evolving over the learning period (blue line, a moving average), with the green dashed line below being the Bertrand-Nash competitive price and the orange dashed line above being the monopoly price. Starting from random exploration, the algorithms gradually learn to raise the price to a level clearly above competition, approaching monopoly, with no communication or collusion instruction throughout.](assets/fig/fig_34_collusion.svg)

The figure above plots the price path of one of the sessions. At the start the algorithms explore over a wide range and the price rises and falls; as exploration decays and the Q-values converge, the two algorithms steadily raise the price to between the competitive price and the monopoly price. This path climbing from competition toward collusion is walked purely by each side's own profit-seeking trial-and-error learning, with no one told to "collude."

The limits of this exercise must be stated honestly. This simulation is a simplified version of Calvano and coauthors' setup, with fewer price grid points and a learning period fixed at one million periods, whereas the original runs until the strategies converge. The price-based collusion index it yields (mean 0.64) is measured on a different basis from the profit-based 0.7 to 0.9 of the original, so the gap does not show that the simplified version is more conservative (see Section 2.1); the two agree in direction and magnitude. The collusion index varies quite a bit across sessions, from 0.41 to 0.79, and this itself is a real feature: the convergence point of reinforcement learning is random, and you cannot expect it to land in the same place every time. The correct reading treats this simulation as a demonstration of the mechanism rather than a precise quantitative prediction.

### 6.3 The complementarity of AI access and workflow training

The third walkthrough uses a $2\times2$ factorial design, randomizing AI access and workflow training separately. In the true DGP the AI standalone effect is $0.05$, the training standalone effect is $0.12$, and the interaction of the two is $0.18$:

```r
log_productivity <- 2 + 0.05*ai_access + 0.12*workflow_training +
  0.18*ai_access*workflow_training + 0.25*baseline_skill + error
fit <- lm(log_productivity ~ ai_access*workflow_training + baseline_skill)
```

In the actual estimation, the AI effect is $0.073$ without workflow training and $0.237$ with it; the difference between the two, that is, the interaction estimate, is $0.164$ (SE $0.028$). Random factorial assignment gives this interaction a clear policy interpretation: training raises the marginal effect of AI access by about $0.164$ on the log scale.

![The factorial cell means of AI access and workflow training. Without training the slope between the two points is smaller, and after training the productivity gain from AI access becomes clearly larger; the difference in the slopes of the two lines is the organizational complementarity.](assets/fig/fig_34_complements.svg)

This walkthrough also shows why an observational usage interaction is not enough. If training is chosen by managers after the fact, or only the workers who feel AI is useful keep using it, a positive interaction may come from selection. The reason this example can give causal complementarity is that the two margins, access and training, are both independently randomized.

### 6.4 Where the three threads converge

Putting the three walkthroughs together, the theme of this chapter becomes clear. In the generated-regressor thread, AI is the researcher's tool, and its output carries measurement error; in the algorithmic-collusion thread, AI is a market agent, and its pricing carries learned strategy; in the factorial-adoption thread, AI is a production technology, and its return depends on organizational complements. Ignore any one layer, and you will mis-write a conditional result as a stable "AI effect."

The reconciliation of this section can be summarized as follows: in the generated-regressor walkthrough, naive OLS attenuates the truth 2.0 to 0.998, and the reliability correction pulls it back to 2.068; in the algorithmic-collusion walkthrough, two non-communicating Q-learners all converge to supra-competitive prices across six sessions, with a price-based collusion index averaging 0.64; in the adoption walkthrough, workflow training raises the AI effect from $0.073$ to $0.237$, with an interaction of $0.164$.

## 7 Failure modes and robustness

In the simulation the measurement error is classical and the algorithm setup is given, but in real research these premises can fail to hold at any moment. This section lays out the most common failure modes and their operational responses.

Non-classical prediction error is the number-one threat. The attenuation formula in Section 3 assumes the measurement error is independent of the truth and additive, but machine prediction error often is not like this; it may be systematically high at some values and low at others, or correlated with other covariates. In this case the simple reliability correction no longer applies, and the direction of the bias is no longer guaranteed to be toward zero. The operational response is to prefer methods that make no assumption about the error structure, that is, the class of tools like PPI and DSL that correct directly using the gold sample, which are robust to the form of the prediction error; and at the same time, check on the gold sample directly the relationship of the prediction error with the truth and with the covariates, to judge whether the classical assumption is wildly off.

The representativeness of the gold sample cannot be taken for granted. The validity of all corrections rests on the one thing that "the gold sample is randomly drawn from the same population, with reliable true labels." If the easy-to-judge samples were picked during human verification, or the annotators themselves have a systematic bias, the correction swaps one bias for another. The response is to truly randomize the drawing of the gold sample, apply quality control and inter-annotator agreement checks to the annotation, and honestly state in the report the size and drawing method of the gold sample.

Outcome variable or regressor, be sure to tell them apart. Section 3.1 stressed that the bias form differs depending on whether the machine-generated variable is on the left-hand side or the right-hand side. Taking a machine-predicted outcome variable as a clean outcome and regressing directly, if the prediction error is correlated with the regressor, biases the coefficient, and the direction requires specific analysis. Before building the model, first make clear which variable is machine-generated and which side of the equation it is on, an indispensable step.

The conclusions of the algorithmic-collusion simulation must be extrapolated with caution. The simulation in Section 6.2 proves collusion can arise on its own, but it depends on a specific algorithm, a specific game, and an extremely long learning period. Real-world pricing algorithms are not necessarily tabular Q-learning, the market does not necessarily give enough time for trial and error, and demand and cost shocks are far richer than in the simulation. Reading "algorithms collude in the simulation" directly as "algorithms are colluding in real markets" does not hold, and the latter needs field evidence like Assad and coauthors'. Conversely, dismissing the simulation result as entirely unrealistic is also wrong; the mechanism it reveals (agreement-free collusion sustained by reward and punishment) is a real theoretical possibility that regulators should take seriously.

The legal definition of collusion remains unsettled. Whether algorithms reaching a collusive outcome with no meeting of the minds constitutes an "agreement" under antitrust law is currently undecided. Whether multiple firms sharing the same third-party pricing software constitute a hub-and-spoke collusion is also under debate. This is not a question econometrics can answer, but an empirical researcher, when interpreting evidence on algorithmic pricing, should be clear that what they measure is a market consequence, and that the assignment of legal liability is a separate layer of judgment.

The externalities of data will distort data-based measurement and policy. The privacy externality discussed in Section 2.2 means that individuals' data-sharing decisions are not independent, and one person's data leak reveals information about others. Any measurement or intervention that relies on user data may be distorted in welfare evaluation by this externality, and simply summing the effects on individuals will miss the cross-individual spillover. Research involving data policy should explicitly incorporate this externality rather than treat it as a simple aggregation of independent individuals.

Stringing these failure modes together, the credibility of empirical work in the AI era ultimately hinges on one clear-eyed recognition: the output of AI is not the truth. As a measurement tool, its predictions carry error and need gold annotation to calibrate and a skepticism about the error structure; as an economic agent, its behavior carries strategy, and simulation can reveal the mechanism but field evidence is needed to ground it; as a factor, data carries externalities, and both allocation and welfare evaluation must count in non-rivalry and privacy externalities. No technical device can replace substantive judgment about "how this number was machine-generated and where it is wrong."

## 8 Further reading

::: {.readings}
Required reading, in the suggested reading order:

- Angelopoulos, Bates, Fannjiang, Jordan and Zrnic (2023, Science). The benchmark for valid inference using machine predictions as data, the full version of Sections 4.4 and 5.2 of this chapter; focus on how the rectifier corrects bias without trusting the model.
- Calvano, Calzolari, Denicolò and Pastorello (2020, American Economic Review). The simulation benchmark for algorithmic collusion, the full version of Sections 2.1 and 5.1 of this chapter; focus on how the reward-punishment mechanism emerges on its own.
- Gentzkow, Kelly and Taddy (2019, Journal of Economic Literature). The framework survey of text-as-data, understanding how text becomes a generated variable and the inference problems that follow.
- Pagan (1984, International Economic Review). The classic result on generated regressors, understanding why the standard error of a two-step estimation must be corrected.
- Jones and Tonetti (2020, American Economic Review). The economics of data non-rivalry, understanding why the allocation of data may be too narrow due to hoarding.
- Brynjolfsson, Li and Raymond (2025, Quarterly Journal of Economics). Productivity, heterogeneity, and knowledge diffusion of generative AI in real customer-support work.
- Brynjolfsson, Rock and Syverson (2021, AEJ: Macroeconomics). GPTs, organizational complements, and the Productivity J-curve.

Further reading:

- Wang, McCormick and Leek (2020, PNAS). Post-prediction inference for machine-predicted outcome variables.
- Egami, Hinck, Stewart and Wei (2023, NeurIPS). Design-based supervised learning for LLM annotation, using expert labels to correct the bias of LLM annotation.
- Battaglia, Christensen, Hansen and Sacher (2024). Joint estimation combining variable generation and the downstream regression into one step.
- Klein (2021, RAND Journal of Economics). Price cycles of algorithmic collusion under sequential pricing.
- Assad, Clark, Ershov and Xu (2024, Journal of Political Economy). Field evidence of algorithmic pricing in the German gasoline market.
- Miklós-Thal and Tucker (2019, Management Science). Why better demand prediction may instead unravel collusion.
- Hagiu and Wright (2023, RAND Journal of Economics). Data network effects and when a data advantage is durable.
- Acemoglu, Makhdoumi, Malekian and Ozdaglar (2022, AEJ: Microeconomics). The privacy externality of data and "too much data."
- Ludwig and Mullainathan (2024, Quarterly Journal of Economics). Treating machine learning as a tool for hypothesis generation rather than mere prediction.
- Comin and Hobijn (2010, American Economic Review). Cross-technology evidence on first adoption and intensive diffusion.
:::

::: {.apa-refs}
- Acemoglu, D., Makhdoumi, A., Malekian, A., & Ozdaglar, A. (2022). Too much data: Prices and inefficiencies in data markets. *American Economic Journal: Microeconomics, 14*(4), 218-256. https://doi.org/10.1257/mic.20200200
- Angelopoulos, A. N., Bates, S., Fannjiang, C., Jordan, M. I., & Zrnic, T. (2023). Prediction-powered inference. *Science, 382*(6671), 669-674. https://doi.org/10.1126/science.adi6000
- Assad, S., Clark, R., Ershov, D., & Xu, L. (2024). Algorithmic pricing and competition: Empirical evidence from the German retail gasoline market. *Journal of Political Economy, 132*(3), 723-771. https://doi.org/10.1086/726906
- Battaglia, L., Christensen, T., Hansen, S., & Sacher, S. (2024). *Inference for regression with variables generated by AI or machine learning* (Cowles Foundation Discussion Paper No. 2421). Yale University. https://arxiv.org/abs/2402.15585
- Brynjolfsson, E., Li, D., & Raymond, L. R. (2025). Generative AI at work. *The Quarterly Journal of Economics, 140*(2), 889-942. https://doi.org/10.1093/qje/qjae044
- Brynjolfsson, E., Rock, D., & Syverson, C. (2021). The productivity J-curve: How intangibles complement general purpose technologies. *American Economic Journal: Macroeconomics, 13*(1), 333-372. https://doi.org/10.1257/mac.20180386
- Calvano, E., Calzolari, G., Denicolò, V., & Pastorello, S. (2020). Artificial intelligence, algorithmic pricing, and collusion. *American Economic Review, 110*(10), 3267-3297. https://doi.org/10.1257/aer.20190623
- Comin, D., & Hobijn, B. (2010). An exploration of technology diffusion. *American Economic Review, 100*(5), 2031-2059. https://doi.org/10.1257/aer.100.5.2031
- Egami, N., Hinck, M., Stewart, B. M., & Wei, H. (2023). Using imperfect surrogates for downstream inference: Design-based supervised learning for social science applications of large language models. In *Advances in Neural Information Processing Systems 36*. https://arxiv.org/abs/2306.04746
- Eloundou, T., Manning, S., Mishkin, P., & Rock, D. (2023). *GPTs are GPTs: An early look at the labor market impact potential of large language models*. https://arxiv.org/abs/2303.10130
- Gentzkow, M., Kelly, B., & Taddy, M. (2019). Text as data. *Journal of Economic Literature, 57*(3), 535-574. https://doi.org/10.1257/jel.20181020
- Hagiu, A., & Wright, J. (2023). Data-enabled learning, network effects, and competitive advantage. *The RAND Journal of Economics, 54*(4), 638-667. https://doi.org/10.1111/1756-2171.12453
- Jones, C. I., & Tonetti, C. (2020). Nonrivalry and the economics of data. *American Economic Review, 110*(9), 2819-2858. https://doi.org/10.1257/aer.20191330
- Klein, T. (2021). Autonomous algorithmic collusion: Q-learning under sequential pricing. *The RAND Journal of Economics, 52*(3), 538-558. https://doi.org/10.1111/1756-2171.12383
- Ludwig, J., & Mullainathan, S. (2024). Machine learning as a tool for hypothesis generation. *The Quarterly Journal of Economics, 139*(2), 751-827. https://doi.org/10.1093/qje/qjad055
- Miklós-Thal, J., & Tucker, C. (2019). Collusion by algorithm: Does better demand prediction facilitate coordination between sellers? *Management Science, 65*(4), 1552-1561. https://doi.org/10.1287/mnsc.2019.3287
- Pagan, A. (1984). Econometric issues in the analysis of regressions with generated regressors. *International Economic Review, 25*(1), 221-247. https://doi.org/10.2307/2648877
- Wang, S., McCormick, T. H., & Leek, J. T. (2020). Methods for correcting inference based on outcomes predicted by machine learning. *Proceedings of the National Academy of Sciences, 117*(48), 30266-30275. https://doi.org/10.1073/pnas.2001238117
:::
