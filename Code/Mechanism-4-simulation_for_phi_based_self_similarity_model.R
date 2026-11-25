


#####################################################################
# END
#####################################################################

### simulation for phi-based self-similarity model 
simulate_system <- function(phi, T = 200, sigma = 1) {
  X <- numeric(T)
  eps <- rnorm(T, mean = 0, sd = sigma)
  
  X[1] <- eps[1]
  
  for (t in 2:T) {
    #feedback <- sum(phi^(2:t) * X[(t-1):1])
	feedback <- sum(phi^(2:t) * X[(t-1):1])
    X[t] <- eps[t] - feedback
  }
  
  data.frame(
    time = 1:T,
    X = X,
    cumulative = cumsum(X),
    phi = phi
  )
}

#

#set.seed(123)

res <- do.call(rbind,lapply(c(0.01,(sqrt(5)-1)/2,0.99),function(x){simulate_system(x,T=5000,sigma=0.2)}))

#
library(ggplot2)


ggplot(res, aes(time, cumulative, color = factor(phi))) +
  geom_line(alpha = 0.9) +
  geom_hline(yintercept=0,lty=2)+
  labs(
    x = "Developmental time",
    y = "Cumulative velocity deviation",
    color = expression(Phi),
    title = "Golden-ratio–specific stabilization of stochastic velocity noise"
  ) +
  theme_classic()


(1-phi^2)*phi/(1-phi)=(1+phi)*phi

ggplot(res, aes(x = X, color = factor(phi))) +
  geom_density(alpha = 0.9) +
  geom_hline(yintercept=0,lty=2)+
  labs(
    x = "Developmental time",
    y = "Cumulative velocity deviation",
    color = expression(Phi),
    title = "Golden-ratio–specific stabilization of stochastic velocity noise"
  ) +
  theme_classic()





