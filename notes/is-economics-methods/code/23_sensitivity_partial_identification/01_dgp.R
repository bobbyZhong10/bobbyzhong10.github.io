set.seed(1515)
N <- 8000
x <- rnorm(N)
u_management <- rnorm(N)
adopt <- rbinom(N, 1, plogis(-0.2 + 0.7*x + 0.9*u_management))
tau <- 0.20
log_sales <- 2 + tau*adopt + 0.45*x + 0.65*u_management + rnorm(N, sd=0.75)
pre_policy_complaints <- 0.55*u_management + 0.20*x + rnorm(N, sd=0.80)
nova <- data.frame(x, adopt, log_sales, pre_policy_complaints, u_management)
saveRDS(list(data=nova, truth=tau), "nova_sensitivity.rds")
cat(sprintf("N = %d | adoption = %.3f | true effect = %.3f\n", N, mean(adopt), tau))
