getwd()
setwd('Documents/R programming/Econ 613/assignment3/')

library(bayesm)
data(margarine)
choice <- margarine$choicePrice
demo <- margarine$demos

## Exercise 1
# Average and dispersion in product characteristics
mean_price <- as.matrix(apply(choice[,3:12],2,mean))
sd_price <- as.matrix(apply(choice[,3:12],2,sd))
stat <- cbind(mean_price,sd_price)
colnames(stat) <- c('Avg Price','Price Sd')

# market share count
count <- matrix(0,nrow=10,ncol = 1)
for (i in 1:10){
  count[i,1] <- nrow(choice[choice$choice== i,])
}
rownames(count) <- colnames(choice)[3:12]
colnames(count) <- 'Amount'

## by product
count <- as.data.frame(count)
count$byProduct_share <- count$Amount/(sum(count$Amount))
count$byProduct_share

## by brand
by_brand_share <- count[1:7,2]
by_brand_share[1] <- by_brand_share[1] +count[8,2]
by_brand_share[3] <- by_brand_share[3] +count[9,2]
by_brand_share[4] <- by_brand_share[4] +count[10,2]
by_brand_share <- as.matrix(by_brand_share)
rownames(by_brand_share) <- c("PPk","PBB","PFl","PHse","PGen","PImp","PSS")

## by tub/stk
tub_v_stk <- matrix(0,ncol = 1,nrow = 2)
rownames(tub_v_stk) <- c('tub','stk')
colnames(tub_v_stk) <- 'Market Share'
tub_v_stk[1,1] <- sum(count[7:10,2]) # for tub
tub_v_stk[2,1] <- sum(count[1:6,2])  # for stk

## mapping
merged <- merge(choice,demo,by = 'hhid',all.x = T)
mapping <- table(merged$choice,merged$Income)
col_sum <- apply(mapping,2,sum)
for (i in 1:14){
  mapping[,i] <- mapping[,i]/col_sum[i]
}
round(mapping,2)

## Exercise 2
# The product characteristic does change and we should use conditional logit model
d <- matrix(0,ncol=10,nrow = nrow(choice))
colnames(d) <- colnames(choice)[3:12]
for (i in 1:nrow(d)){
  k <- choice[i,2]
  d[i,k] <- 1
}
price <- choice[,3:12]

prob_clogit <- function(beta,X){
  l <- length(beta)
  b <- beta[1]
  alpha <- beta[2:l]
  mat <- X*b
  for (i in 1:length(alpha)){
    mat[,i+1] <- mat[,i+1]+alpha[i]
  }
  exp_mat <- exp(mat)
  row_sum <- apply(exp_mat,1,sum) 
  exp_mat <- exp_mat/row_sum
  return(exp_mat)
}

ll_clogit <- function(d,X,beta){
  prob <- prob_clogit(beta,X)
  mat <- d*log(prob)
  summation <- -sum(mat)
  return(summation)
}
start <- rep(0,10)

op_clogit <- optim(par=start,ll_clogit,X=price,d=d,method = 'BFGS')
clogit_beta <- op_clogit$par
clogit_beta # beta:-6.6566340 alpha: -0.9543259  1.2969965 -1.7173298 -2.9040264 -1.5153021  0.2517927  1.4648942  2.3575437 -3.8966267
## Interpretation: since beta is negative, we can expect that an increase in price of one alternatives
#                  decreases the probability of choosing this alternative and increase the probability
#                  of choosing other alternatives

## Exercise 3
# The product characteristic does not change, we should use multinomial logit
income <- as.matrix(merged$Income,ncol=1)

prob_mlogit <- function(X,coef){
  alpha <- coef[1:9]
  beta <- coef[10:18]
  l <- length(alpha)
  mat <- matrix(0,ncol = l,nrow = nrow(X))
  for (i in 1:l){
    mat[,i] <- beta[i]*X+alpha[i]
  }
  exp_mat <- exp(mat)
  e <- matrix(1,ncol = 1,nrow = nrow(X))
  exp_mat <- cbind(e,exp_mat)
  row_sum <- apply(exp_mat,1,sum) 
  exp_mat <- exp_mat/row_sum
  return(exp_mat)
}

ll_mlogit <- function(d,X,coef){
  prob <- prob_mlogit(X,coef)
  mat <- d*log(prob)
  summation <- -sum(mat)
  return(summation)
}

start1 <- rep(0,18)
op_mlogit <- optim(par=start1,ll_mlogit,X=income,d=d,method='BFGS')
mlogit_beta <- op_mlogit$par
mlogit_beta
## alpha: -0.843545649 -2.397656003 -1.199428121 -1.688616844 -4.137055731 -1.529169108 -2.846055103 -2.573291074 -4.279712750
## beta:  -0.003156338  0.014507166  0.003980338 -0.001328126  0.030527384 -0.007002723  0.022807121  0.017661767  0.010698254
# Interpretation: if beta[i] is positive, it means that an increase in income will leads to greater probability of choosing
#                 this alternative relative to the base case. if beta[i] is negative, it means that an increase in income
#                 will leads to smaller probability of choosing this alternative relative to the base case

## Exercise 4
library(ramify)
marginal_clogit <- function(X,beta,d){
  prob_ik <- prob_clogit(beta,X)
  b <- beta[1]
  mean_ME <- matrix(0,ncol = 10,nrow = 10)
  for (i in 1:ncol(prob_ik)){
    prob_ij <- resize(rep(prob_ik[,i],10),ncol = ncol(prob_ik),
                      nrow = nrow(prob_ik))
    margin <- prob_ij*(d-prob_ik)*b
    mean_ME[i,] <- colMeans(margin)
  }
  return(mean_ME)
}
ME_clogit <- marginal_clogit(price,clogit_beta,d)
round(ME_clogit,4)
## Intepretation: the mean marginal effect measures the effect of change in probability of choosing choice J 
#                 by changing in one unit the value of the regressor for kth alternative in average level

marginal_mlogit <- function(X,coef){
  prob <- prob_mlogit(X,coef)
  beta <- c(0,coef[10:18])
  beta_bar <- matrix(0,ncol = ncol(prob),nrow = nrow(prob))
  for (i in 1:nrow(beta_bar)){
    sums <- sum(prob[i,]*beta)
    beta_bar[i,] <- sums
  }
  mat <- matrix(0,ncol = ncol(prob),nrow = nrow(prob))
  for (j in 1:ncol(mat)){
    mat[,j] <- beta[j]
  }
  marginal <- prob*(mat-beta_bar)
  mean_ME <- colMeans(marginal)
  return(mean_ME)
}

ME_mlogit <- marginal_mlogit(income,mlogit_beta)
ME_mlogit # in average level
# -0.0010504137 -0.0009016311  0.0006266867  0.0001660472 -0.0002794477  0.0004431356 -0.0006821378  0.0008861440  0.0007338590  0.0000577577
## Interpretation: the mean ME for mlogit measures the average effect on the probability of choosing J 
#                  by one unit a regressor that takes the same value across all alternatives 

# Exercise 5
prob_mixlogit <- function(X_clogit,X_mlogit,beta){
  l <- (length(beta)-1)/2
  alpha <- beta[1:l]
  beta_mlogit <- beta[(l+1):(2*l)]
  beta_clogit <- beta[length(beta)]
  mat <- beta_clogit*X_clogit
  for (i in 1:l){
    mat[,i+1] <- alpha[i]+beta_mlogit[i]*X_mlogit+mat[,i+1]
  }
  exp_mat <- exp(mat)
  row_sum <- apply(exp_mat,1,sum) 
  exp_mat <- exp_mat/row_sum
}

ll_mixlogit <- function(X_clogit,X_mlogit,d,beta){
  prob <- prob_mixlogit(X_clogit,X_mlogit,beta)
  mat <- d*log(prob)
  summation <- -sum(mat)
  return(summation)
}
start2 <- rep(0,19)

op_mixed <- optim(par=start2,ll_mixlogit,X_clogit=price,X_mlogit=income,d=d,method = 'BFGS')
op_mixed$par # first 9 values are alphas, 10:18 are for coef of X_mlogit
# alphas: -0.838705945  0.891148169 -1.826370582 -2.871247434 -2.454001559  0.498968897  0.805453868  1.866785193 -4.140083624
# betas:  -0.004333800  0.014258958  0.004025557 -0.001264787  0.029710007 -0.009327126  0.021914644  0.016902350  0.008674428
# betas: extra for clogit part -6.659699884

# remove choice 10 data
merged2 <- merged[merged$choice!=10,]
price2 <- merged2[,3:11]
income2 <- as.matrix(merged2$Income,ncol=1)
start3 <- rep(0,17)
d2 <- d[d[,10]!=1,]
d2 <- d2[,-10]

op_mixed2 <- optim(par=start3,ll_mixlogit,X_clogit=price2,X_mlogit=income2,d=d2,method = 'BFGS')
op_mixed2$par # first 8 values are alphas, 9:16 are for coef of X_mlogit
# alphas: -0.843787287  0.892575161 -1.825287978 -2.871731853 -2.452713672  0.497321521  0.809064852  1.870859459
# betas:  -0.004210987  0.014186227  0.003967073 -0.001313038  0.029612711 -0.009289775  0.021785560  0.016771027
# betas: extra for clogit part -6.662222287

coeff <- op_mixed$par[c(-9,-18)]
LL_bf <- ll_mixlogit(price2,income2,d2,coeff)
LL_br <- op_mixed2$value

MTT <- 2*(LL_bf-LL_br)
MTT 
## the critical value of chi-square test with 95% confidence & df = 9 is 16.919, 
## can't reject the null hypothesis. The market share for choice 10 is very small
## only 33 out of 4470 choose it. It might not cause significant effect on our 
## estimation when we exclude it from our data