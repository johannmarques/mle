library(tidyverse)

theme_set(theme_bw())

purple_rain_colors <- c(
    "#533a70",
    "#a489c2",
    "#cbbcdc",
    "#e5ddee")

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
# Since \mu_1, \mu_2, \sigma_1 and \sigma_2 were not specified,
# I'm considering the following
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
                   alpha = .8, position = 'identity', binwidth = 0.25, color = 'black') +
  labs(fill = '', x = '', y = '') + scale_fill_manual(values = purple_rain_colors[c(1,2)])

ggplot() +
  geom_histogram(aes(y = after_stat(density), x = e),
                 alpha = .8, position = 'identity', binwidth = 0.25, color = 'black', fill = purple_rain_colors[3]) +
  labs(fill = '', x = '', y = '')

tibble(e, e1, e2)

# Simulating model 2
# Since a and b were not specified, I'm considering the following

a = 1.75
b = .8

u = rgamma(n = n, shape = a, rate = b)

y_2 = beta[1] + beta[2] * x + u

# Simulating model 3

v = rt(n = n, df = 1)

y_3 = beta[1] + beta[2] * x + v

rbind(
  tibble(y = y_1, x = x, group = 'Modelo 1'),
  tibble(y = y_2, x = x, group = 'Modelo 2'),
  tibble(y = y_3, x = x, group = 'Modelo 3')
) %>%
  ggplot(aes(x = x, y = y, color = group)) +
    geom_point(alpha = .75) + facet_wrap(~group, scales = 'free') + scale_y_log10() +
  scale_color_manual(values = as.character(purple_rain_colors[1:3])) + labs(color = '')

summary(data.frame(x, y_1, y_2, y_3))

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
init_norm = lm(y_1 ~ x) %>%
  {c(coef(.), sigma2 = summary(.)$sigma^2)} # Adding OLS estimate for \sigma^2

mod1_norm <- optim(par = init_norm,
      fn = function(theta){-l_norm0(theta = theta, y = y_1, x=x)})

# Model 2

# Initial condition

init_gamma = lm(y_2 ~ x) %>%
  {c(coef(.), sigma2 = summary(.)$sigma^2)}

mod2_norm <- optim(par = init_gamma,
      fn = function(theta){-l_norm0(theta = theta, y = y_2, x=x)})

# Model 3

# Initial condition

init_t = lm(y_3 ~ x) %>%
  {c(coef(.), sigma2 = summary(.)$sigma^2)}

mod3_norm <- optim(par = init_t,
      fn = function(theta){-l_norm0(theta = theta, y = y_3, x=x)})

# Exercise 4

# Model 1

# Objective funciton
l_normmix <- function(theta, mu1, mu2, sigma2_1, sigma2_2, y, x){
  -length(y)*log(2)+ sum(log(1/sqrt(2*pi*sigma2_1)*
                               exp(-1/(2*sigma2_1)*
                                     (y - theta[1] - theta[2] * x - mu1)^2) +
                           1/sqrt(2*pi*sigma2_2)*
                           exp(-1/(2*sigma2_2) *
                                 (y - theta[1] - theta[2] * x - mu2)^2)))
}

# Initial condition previously defined #

mod1_true <- optim(par = init_norm[-3],
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

# Initial condition previously defined #

# The objective function evaluated at the initial condition returns an error.
# So I'll be subtracting 1

mod2_true <- optim(par = init_gamma[-3]-1,
      fn = function(theta){-l_gamma(theta = theta, a = 1.75, b = .8, y = y_2, x = x)})

# Model 3

# Objective funciton
l_t <- function(theta, y, x){
  length(y)*(log(gamma((theta[3]+1)/2)) - log(sqrt(theta[3]*pi)) - log(gamma(theta[3]/2))) -
    (theta[3]+1)/2*sum(log(1 + ((y - theta[1] - theta[2]*x)^2)/theta[3]))
}

init_t2 = c(coef(
  lm(y_3 ~ x)
), df = 0.0001) # Setting a pretty small degree of freedom as an initial guess

mod3_true <- optim(par = init_t2,
      fn = function(theta){-l_t(theta = theta, y = y_3, x=x)})

plot_curves <- function(...){
  ggplot() +
    geom_point(aes(y =eval(parse(text = paste0('y_', ...))), x=x), alpha = .5, color = purple_rain_colors[1]) +
    geom_abline(slope = beta[2],
                intercept = beta[1],
                linetype = 'F1',
                color = 'grey',
                linewidth = 1) +
    geom_abline(slope = eval(parse(text = paste0('mod', ..., '_true$par[2]'))),
                intercept = eval(parse(text =paste0('mod', ..., '_true$par[1]'))),
                linetype = 'dashed',
                color = purple_rain_colors[1],
                linewidth = 1,
                alpha = .75) +
    geom_abline(slope = eval(parse(text = paste0('mod', ..., '_norm$par[2]'))),
                intercept = eval(parse(text = paste0('mod', ..., '_norm$par[1]'))),
                linetype = 'solid',
                color = purple_rain_colors[1],
                linewidth = 1,
                alpha = .75) +
    labs(x = 'x', y = 'y')
}

invoke_map(.f = plot_curves, .x = 1:3)

eval(parse(text = paste0('y_', 1)))
