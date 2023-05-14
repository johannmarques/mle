library(tidyverse)

# Simulating data generating processes

n = 2000 # Sample size

# Setting random seed
set.seed(1981)

# Exercise 2
# Generating x
x = rnorm(n = n, mean = 5, sd = sqrt(1))

# Defining betas

beta = c(2.35, 0.75)

# Simulating model 1
# Since \mu_1, \mu_2, \sigma_1 and \sigma_2 were not specified, I'm considering the following
mu_1 = 2
mu_2 = 1
sigma_1 = sqrt(1.5)
sigma_2 = sqrt(3)

e1 = rnorm(n = n, mean = mu_1, sd = sigma_1)
e2 = rnorm(n = n, mean = mu_2, sd = sigma_2)

# Defining an auxiliary random variable P ~ Bernoulli(1/2)

P = rbinom(n = n, size = 1, prob = 1/2)

e = P * e1 + (1-P) * e2

y_1 = beta[1] + beta[2] * x + e

tibble(e1, e2) %>%
  pivot_longer(cols = c('e1', 'e2')) %>%
  ggplot(aes(x = value, group = name)) +
    geom_histogram(aes(y = after_stat(density), fill = name),
                   alpha = .3, position = 'identity', binwidth = 0.25, color = 'black')

ggplot() +
  geom_point(aes(y = y_1, x = x))

# Simulating model 2
# Since a and b were not specified, I'm considering the following

a = 1.75
b = .8

u = rgamma(n = n, shape = a, rate = b)

y_2 = beta[1] + beta[2] * x + u

# Simulating model 3

v = rt(n = 500, df = 1)

y_3 = beta[1] + beta[2] * x + v

rbind(
  tibble(y = y_1, x = x, group = '1'),
  tibble(y = y_2, x = x, group = '2'),
  tibble(y = y_3, x = x, group = '3')
) %>%
  ggplot(aes(x = y, color = group)) +
    geom_histogram(alpha = .3) + facet_wrap(~group, scales = 'free')


# Exercise 3.2
# I'm assuming y - Xb ~ N(0, sigma^2)
# For the remaining exercises, I'll be using OLS estimates as initial conditions

# Objective funciton
l_norm0 <- function(theta, y, x){
  -length(y)/2*log(2*pi*theta[3]) + sum(-1/(2*theta[3])*(y - theta[1] - theta[2]*x)^2)
}

l_norm0 <- function(theta, y, x){
  (1/sqrt(2*pi*theta[3]))^length(n)*exp(-1/(2*theta[3])*sum((y - theta[1] - theta[2]*x)^2))
}

# Model 1

# Initial condition
init_norm = c(coef(
  lm(y_1 ~ x)
), sigma2 = ols_mod$sigma^2) # Adding OLS estimate for \sigma^2

optim(par = init_norm,
      fn = function(theta){-l_norm0(theta = theta, y = y_1, x=x)})

# Model 2

optim(par = init_norm,
      fn = function(theta){-l_norm0(theta = theta, y = y_2, x=x)})


# Model 3

optim(par = init_norm,
      fn = function(theta){-l_norm0(theta = theta, y = y_3, x=x)})

# Exercise 4

# Model 1

# Objective funciton
l_normmix <- function(theta, mu1, mu2, sigma2_1, sigma2_2, y, x){
  (length(y)/2)*log((2/(pi*(sigma2_1 + sigma2_2)))) + sum(-2/(sigma2_1+sigma2_2)*(y - theta[1] - theta[2]*x - (mu1+mu2)/2)^2)
}

# Initial condition previously defined #

optim(par = init_norm[-3],
      fn = function(theta){-l_normmix(theta = theta, mu1 = 2, mu2 = 1, sigma2_1 = 1.5, sigma2_2 = 3, y = y_1, x=x)})

# Model 2

# Objective funciton
l_gamma <- function(theta, a, b, y, x){
  if(prod(y - theta[1] - theta[2]*x > 0)){ # Checking whether u_i > 0 \forall i = 1, ..., n
    return(length(y)*log(b^a/gamma(a)) + (a-1)*sum(log(y - theta[1] - theta[2]*x)) - b*sum(y - theta[1] - theta[2] * x))
  }else{
    return(-Inf)
    }
}

# Initial condition

init_gamma = coef(
  lm(y_2 ~ x)
)

# The objective function evaluated at the initial condition returns an error.
# So I'll be subtracting 1

optim(par = init_gamma-1,
      fn = function(theta){-l_gamma(theta = theta, a = 1.75, b = .8, y = y_2, x = x)})

# Model 3

# Objective funciton
l_t <- function(theta, y, x){
  length(y)*(log(gamma((theta[3]+1)/2)) - log(sqrt(theta[3]*pi)) - log(gamma(theta[3]/2))) -
    (theta[3]+1)/2*sum(log(1 + ((y - theta[1] - theta[2]*x)^2)/theta[3]))
}

# Initial condition

init_t = c(coef(
  lm(y_3 ~ x)
), df = 0.0001) # Setting a pretty small degree of freedom as an initial guess

optim(par = init_t,
      fn = function(theta){-l_t(theta = theta, y = y_3, x=x)})

