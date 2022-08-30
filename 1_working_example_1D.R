library(dplyr)
library(tmg)
library(kergp)
library(ggplot2)
library(viridis)
library(rstan)
library(ggpubr)
library(ggExtra)
library(MVN)
library(mvtnorm)

options(mc.cores = parallel::detectCores())

# pdf_1 <- function(x){
#   ifelse(x<=0.2, 0.75, 
#          ifelse(x<=0.7, 0.25, 0.75))
# }
pdf_1 <- function(x){
  res <- ifelse(x<=0.1, 0, 
                ifelse(x<=0.2, 4*(x-0.1),
                       ifelse(x<=0.3, 0.4-2*(x-0.2),
                              ifelse(x<=0.4, 0.2,
                                     ifelse(x<=0.6, 0.2+4*(x-0.4),
                                            ifelse(x<=0.7, 1-4*(x-0.6),
                                                   ifelse(x<=0.8, 0.6+2*(x-0.7),
                                                          ifelse(x<=0.9, 0.8, 0.8-2*(x-0.9)))))))))
  res <- 1/(1+exp(-6*res+3))
  return(res)
}
class_1 <- function(x){
  y <- pdf_1(x)
  x <- c(x)
  r <- runif(length(x))
  return(1*(r<=y))
}
setwd("C:/Users/athen/Research/SurelyLipschitzGaussianClassifiers")
source("./0_source_all_material.R")

u <- seq(0, 1,, 101)
plot(u, pdf_1(u), "l", ylim=c(0, 1))
set.seed(1)
u <- seq(0, 1,, 21)

n_samp <- 1000

data <- data.frame(x1=c(sample(u, n_samp, replace=TRUE)))
data$class <- class_1(data$x1)

name_index <- c("x1")
kernel_choice <- "Mat52"
param_lengthscale_start <- c(0.1)

domain_bounds <- matrix(c(0,  1), nrow=1)
dim <- 1
n_nodes <- 11

C_lip <- 60

## OPTIM

prepared_optim <- prepare_optim(data=data,
                                name_index=name_index,
                                domain_bounds = domain_bounds,
                                dim=dim,
                                n_nodes=n_nodes,
                                C_lip = C_lip, 
                                kernel_choice = "Mat52",
                                param_lengthscale = param_lengthscale_start,
                                mean_val = 0)
weight_start=rep(0, nrow(prepared_optim$coord_nodes))
lengthscale_start=param_lengthscale_start
cleaned_data <- prepared_optim$data_clean

sd_logtheta<-1.5

list_len <- seq(0.01, 1, 0.01)
neg_log_like<-sapply(seq_along(list_len), function(i){
  res_opt <- run_optim_fixed_len(prepared_optim, sd_logtheta=sd_logtheta,
                                 lengthscale_candidates=list_len[i])$value
})

lengthscale <-list_len[which.min(neg_log_like)]

ggplot(data.frame(theta=list_len, y=neg_log_like), aes(x=theta, y=y))+
  geom_line()+
  theme_bw()+
  ylab("Negative log posterior")+
  ggtitle("Fixed theta, optimised epsilon")+
  geom_vline(xintercept=lengthscale, lty=2, col="red")


res_optim <- run_optim_fixed_len(prepared_optim, sd_logtheta=sd_logtheta,
                                 lengthscale_candidates=lengthscale, return_Hess = TRUE)
new_coord <- data.frame(x1=seq(0, 1,, 101))

res <- field_predict(new_coord, weights=res_optim$weights, name_index, domain_bounds, dim,
                     n_nodes)
ggplot(res, aes(x=x1, y=prob, group=real))+
  geom_line(col="black", lty=2)+
  geom_line(data=data.frame(x1=u, prob=pdf_1(u), real="true"), col="red")+
  theme_bw()+
  geom_point(prepare_data(data, name_index), inherit.aes = FALSE,
             mapping=aes(x=x1, y=success/multiplicity, size=multiplicity),
             pch=21, col="blue", fill="white", alpha=0.5)+
  geom_text(prepare_data(data, name_index), inherit.aes = FALSE,
            mapping=aes(x=x1, y=success/multiplicity, label=paste0("nx=", multiplicity)),
            col="blue", size=1.5)+
  scale_size_area()+
  theme(legend.position = "none")+
  ggtitle("Ground-truth (red), Lipschitz-GP classifiers (MAP, in black) and\nempirical probabilities in the labelled trainset (blue circles)")+
  ylab("Probability to be in class 1 at x")+
  xlab("x")+
  labs(caption=paste0("Size of trainset: ", nrow(data)))

# MCMC 
prep_rstan <- prepare_rstan(data=data, name_index = name_index, domain_bounds = domain_bounds,
                            dim=dim, n_nodes=n_nodes,
                            C_lip=C_lip, kernel_choice = kernel_choice,
                            param_lengthscale = c(lengthscale), mean_val = 0)

set.seed(1)
fit <- stan(model_code = prep_rstan$scode, model_name = "My model", 
            data=prep_rstan$data, init=res_optim$weights, 
            iter = 100000, verbose = FALSE, thin=10, warmup=50000)
print(fit)
samp <- rstan::extract(fit)$x
samp_thin <- samp[seq(1, 10000, 100), ]

# Laplace approximation
n_simu <- 100
invCov <- res_optim$Hess
eiginvCov <- eigen(invCov)
sqrtCov <- (eiginvCov$vectors) %*% diag(sqrt(1/eiginvCov$values)) %*% t(eiginvCov$vectors)
Cov <- (eiginvCov$vectors) %*% diag(1/eiginvCov$values) %*% t(eiginvCov$vectors)

candidate_weights <- t(res_optim$weights+
                         sqrtCov%*%matrix(rnorm(n_simu*length(res_optim$weights)), ncol=n_simu))

res_laplace <- field_predict(new_coord, weights=candidate_weights, name_index, domain_bounds, dim,
                             n_nodes)
ggplot(res_laplace, aes(x=x1, y=prob, group=real))+
  geom_line(col="grey", alpha=0.5)+
  geom_line(data=res, col="black", lty=2)+
  geom_line(data=data.frame(x1=u, prob=pdf_1(u), real="true"), col="red")+
  theme_bw()+
  geom_point(prepare_data(data, name_index), inherit.aes = FALSE,
             mapping=aes(x=x1, y=success/multiplicity, size=multiplicity),
             pch=21, col="blue", fill="white", alpha=0.5)+
  geom_text(prepare_data(data, name_index), inherit.aes = FALSE,
            mapping=aes(x=x1, y=success/multiplicity, label=paste0("nx=", multiplicity)),
            col="blue", size=1.5)+
  scale_size_area()+
  theme(legend.position = "none")+
  ggtitle("Ground-truth (red), Lipschitz-GP classifiers (MAP, in black) and\nempirical probabilities in the labelled trainset (blue circles)")+
  ylab("Probability to be in class 1 at x")+
  xlab("x")+
  labs(caption=paste0("Size of trainset: ", nrow(data)))

indPlot <- c(1, 5)
joint_density_Laplace <- data.frame(t(t(as.matrix(expand.grid(seq(-1, 1,, 201), 
                                                   seq(-1, 1,, 201))))*2*
                             sqrt(diag(Cov))[indPlot]+
                             res_optim$weights[indPlot]))

colnames(joint_density_Laplace) <- c("x", "y")

joint_density_Laplace$density <- apply(joint_density_Laplace, 1, function(x){
  dmvnorm(x, mean=res_optim$weights[indPlot], sigma=Cov[indPlot, indPlot])
})
# plot of the joint density
plot1 <- ggplot(joint_density_Laplace ) +
  geom_raster(mapping=aes(x=x, y=y, fill=density))+
  # geom_point(data.frame(x=candidate_weights[,indPlot[1]],
              #           y=candidate_weights[,indPlot[2]]),
              # mapping=aes(x=x, y=y), col="grey", pch="x")+
  coord_cartesian(ylim=range(joint_density_Laplace$y),
                  xlim=range(joint_density_Laplace$x))+
  scale_fill_viridis(
    option = "magma", 
    direction = 1)+
  theme(legend.position = "none")+
  ylab(paste0("Weight ", indPlot[2], " (at x=", prepared_optim$coord_nodes[indPlot[2]],")" ))+
  xlab(paste0("Weight ",indPlot[1], " (at x=", prepared_optim$coord_nodes[indPlot[1]],")" ))+
  ggtitle("Laplace approximation")+
  theme_bw()


# plot of the joint density
plot2 <- ggplot(data.frame(x=samp[, indPlot[1]],
                  y=samp[, indPlot[2]]),
       mapping=aes(x=x, y=y) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
  geom_point(data=data.frame(x=samp[, indPlot[1]],
                             y=samp[, indPlot[2]]), col="grey", pch=".", alpha=0.4)+
  coord_cartesian(ylim=range(joint_density_Laplace$y),
                  xlim=range(joint_density_Laplace$x))+
  scale_fill_viridis(
    option = "magma", 
    direction = 1)+
  theme_bw()+
  theme(legend.position = "none")+
  ylab(paste0("Weight ", indPlot[2], " (at x=", prepared_optim$coord_nodes[indPlot[2]],")" ))+
  xlab(paste0("Weight ",indPlot[1], " (at x=", prepared_optim$coord_nodes[indPlot[1]],")" ))+
  ggtitle("MCMC")

ggarrange(plot1, plot2, ncol=2, legend = "bottom", common.legend = TRUE)

ggMarginal(plot2, type = "histogram")
mvn(samp)
