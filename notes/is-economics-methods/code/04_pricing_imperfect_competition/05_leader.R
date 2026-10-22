# ------------------------------------------------------------------
# 05_leader.R -- block C4 (deterministic): price leadership in the Cirrus-Garnet logit duopoly at the C1
# baseline. Simultaneous pricing (Bertrand-Nash) against sequential pricing, with Cirrus leading (price_leader)
# and with Garnet leading (price_leader on the reversed duopoly); joint maximisation for reference. Only the
# comparison "with a price leader both prices exceed the simultaneous prices" is the chapter's claim (the card's
# Brown and MacKay 2023 item); nothing here models pricing frequencies or algorithms.
# Reads truth04.rds. Writes res_leader.rds.
# ------------------------------------------------------------------
local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1L) setwd(dirname(normalizePath(gsub("~+~", " ", f, fixed = TRUE))))
})
source("pricing_fns.R")
P4 <- readRDS("truth04.rds")$params$C4
d <- unname(P4$delta); a <- P4$alpha; cc <- P4$c

p_sim <- bertrand_logit(d, a, cc)
p_seq <- price_leader(d, a, cc)                                       # Cirrus (product 1) leads
p_seq_garnet <- rev(price_leader(rev(d), a, rev(cc)))                 # Garnet leads; back in (Cirrus, Garnet) order
p_joint <- joint_max_logit(d, a, cc)
prof <- function(p) logit_profit(p, d, a, cc)
shr <- function(p) logit_share(p, d, a)

res <- list(p_sim = p_sim, p_seq = p_seq, p_seq_garnet = p_seq_garnet, p_joint = p_joint,
            profit_sim = prof(p_sim), profit_seq = prof(p_seq), profit_seq_garnet = prof(p_seq_garnet), profit_joint = prof(p_joint),
            s_sim = shr(p_sim), s_seq = shr(p_seq), s_seq_garnet = shr(p_seq_garnet), s_joint = shr(p_joint),
            cs_sim = log(1 + sum(exp(d - a * p_sim))) / a, cs_seq = log(1 + sum(exp(d - a * p_seq))) / a,
            cs_seq_garnet = log(1 + sum(exp(d - a * p_seq_garnet))) / a, cs_joint = log(1 + sum(exp(d - a * p_joint))) / a)
saveRDS(res, "res_leader.rds")
tab <- rbind(simultaneous = c(p_sim, prof(p_sim)), cirrus_leads = c(p_seq, prof(p_seq)),
             garnet_leads = c(p_seq_garnet, prof(p_seq_garnet)), joint = c(p_joint, prof(p_joint)))
colnames(tab) <- c("p_Cirrus", "p_Garnet", "profit_Cirrus", "profit_Garnet")
print(round(tab, 4))
cat("Saved res_leader.rds\n")
