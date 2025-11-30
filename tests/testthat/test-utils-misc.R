library(testthat)
library(muscal)

# Tests for %||% operator
expect_equal(NULL %||% 5, 5)
expect_equal(3 %||% 4, 3)

# Tests for adam_update_block
set.seed(1)
V <- matrix(rnorm(4),2,2)
G <- matrix(rnorm(4),2,2)
M <- matrix(0,2,2)
V2 <- matrix(0,2,2)
res <- muscal:::adam_update_block(V,G,M,V2,step_count=1,beta1=0.9,beta2=0.999,adam_epsilon=1e-8,learning_rate=0.01)

# manual computation
exp_M <- 0.1*G
exp_V2 <- 0.001*(G*G)
M_hat <- exp_M / (1 - 0.9)
V_hat <- exp_V2 / (1 - 0.999)
step <- 0.01 * M_hat / (sqrt(V_hat) + 1e-8)
exp_V_new <- V - step

expect_equal(res$M, exp_M)
expect_equal(res$V2, exp_V2)
expect_equal(res$V, exp_V_new)
