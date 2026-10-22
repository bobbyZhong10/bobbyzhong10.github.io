set.seed(1111)
N <- 20000
x <- rnorm(N)
treatment <- rbinom(N, 1, plogis(-0.1 + 0.5 * x))
true_satisfaction <- 3 + 0.40 * treatment + 0.55 * x + rnorm(N, sd = 0.80)

# The LLM is systematically more generous on treatment conversations and
# its error changes with customer complexity x: nonclassical label error.
llm_satisfaction <- true_satisfaction + 0.30 * treatment - 0.20 * x +
  rnorm(N, sd = 0.55)

# A platform logging rule preferentially retains escalated and treated chats.
logged <- rbinom(N, 1, plogis(-0.3 + 0.5 * treatment - 0.35 * true_satisfaction))

# Gold labels are audited by simple random sampling from the target population.
gold <- integer(N)
gold[sample.int(N, 1200)] <- 1L

lyra <- data.frame(x, treatment, true_satisfaction, llm_satisfaction, logged, gold)
saveRDS(lyra, "lyra_labels.rds")
cat(sprintf("N = %d | gold = %d | logged share = %.3f\n", N, sum(gold), mean(logged)))
cat(sprintf("mean label error: control %.3f | treatment %.3f\n",
            mean(llm_satisfaction[treatment == 0] - true_satisfaction[treatment == 0]),
            mean(llm_satisfaction[treatment == 1] - true_satisfaction[treatment == 1])))
