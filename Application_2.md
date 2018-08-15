Zero-Adjusted Birnbaum-Saunders regression model
================
Vera Tomazella, JuvÃªncio S. Nobre, Gustavo, H.A. Pereira and Manoel
Santos-Neto
August 14, 2018

## BS nonlinear regression model

In this application, the data set considered is the biaxial fatigue data
analyzed by Rieck and Nedelman (1991), Lemonte and Cordeiro (2011),
Lemonte and Patriota (n.d.), Vanegas, Rondon, and Cysneiros (2012) and
Lemonte, Cordeiro, and Moreno-Arenas (2016). This data set has 46
observations and consists of lifes of a metal piece subjected to cyclic
stretching and compressing, where the response variable \(N\) denotes
the number of cycles to failure of the metal specimen and the
explanatory variable \(W\) is the work per cycle (\(mj/m^3\)). Here, we
use the response variable on a reduced scale to facilitate the
estimation procedure, i.e., \(Y \times 100 \equiv N\). Therefore, we fit
a RBS nonlinear regression model \[
\mu_i = \beta_1 \textrm{e}^{\frac{\beta_2}{w_i}}, \quad i=1,\ldots,46,
\] where \(Y_i \sim \text{BS}(\mu_i,\sigma)\). The adjustment was
performed using the function **nlgamlss()** of the package **gamlss.nl**
whose objective is to allow nonlinear fitting within a GAMLSS model. The
estimates and standard errors obtained were:

``` r
library(gamlss.nl)
```

    ## Warning: package 'gamlss.nl' was built under R version 3.4.4

    ## Loading required package: gamlss

    ## Loading required package: splines

    ## Loading required package: gamlss.data

    ## Loading required package: gamlss.dist

    ## Loading required package: MASS

    ## Loading required package: nlme

    ## Loading required package: parallel

    ##  **********   GAMLSS Version 5.0-6  **********

    ## For more on GAMLSS look at http://www.gamlss.org/

    ## Type gamlssNews() to see new features/changes/bug fixes.

    ## Loading required package: survival

``` r
library(RBS)
```

    ## Loading required package: gmm

    ## Warning: package 'gmm' was built under R version 3.4.4

    ## Loading required package: sandwich

    ## Warning: package 'sandwich' was built under R version 3.4.4

    ## Loading required package: pracma

    ## 
    ## Attaching package: 'RBS'

    ## The following object is masked from 'package:gamlss.data':
    ## 
    ##     oil

    ## The following object is masked from 'package:stats':
    ## 
    ##     residuals

``` r
library(ggplot2)
data("Biaxial", package="ssym")
N <- Biaxial$Life
N <- N/100;N
```

    ##  [1] 32.80 50.46 15.63 47.07  9.77 28.34 22.66 22.08 10.40  7.00 15.83
    ## [12]  4.82  8.04 10.93 11.25  8.84 13.00  8.52  5.80 10.66 11.14  3.86
    ## [23]  7.45  7.36  7.50  3.16  4.56  5.52  3.55  2.42  1.90  1.27  1.85
    ## [34]  2.55  1.95  2.83  2.12  3.27  3.73  1.25  1.87  1.35  2.45  1.37
    ## [45]  2.00  1.90

``` r
W <- Biaxial$Work;W
```

    ##  [1]  11.5  13.0  14.3  15.6  16.0  17.3  19.3  21.1  21.5  22.6  22.6
    ## [12]  24.0  24.0  24.6  25.2  25.5  26.3  27.9  28.3  28.4  28.6  30.9
    ## [23]  31.9  34.5  40.1  40.1  43.0  44.1  46.5  47.3  48.7  52.9  56.6
    ## [34]  59.9  60.2  60.3  60.5  62.1  62.8  66.5  67.0  67.1  67.9  68.8
    ## [45]  75.4 100.5

``` r
plot(N~W)
```

![](Application_2_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
mod1 <- nlgamlss(y=N,mu.fo=~p1*exp(p2/W),sigma.formula = ~1,mu.start =c(1,40),sigma.start =1,family = RBS(mu.link = "identity",sigma.link = "identity"))
mod1$converged
```

    ## [1] 1

``` r
summary(mod1)
```

    ## *******************************************************************
    ## Family:  c("RBS", "BirnbaumSaunders") 
    ## 
    ## Call:  nlgamlss(y = N, mu.fo = ~p1 * exp(p2/W), sigma.formula = ~1,  
    ##     mu.start = c(1, 40), sigma.start = 1, family = RBS(mu.link = "identity",  
    ##         sigma.link = "identity")) 
    ## 
    ## Fitting method: "JL()" 
    ## 
    ## -------------------------------------------------------------------
    ## Mu link function:  identity
    ## Mu Coefficients:
    ##     Estimate  Std. Error  t-value    p-value
    ## p1     1.276      0.1719    7.423  1.148e-13
    ## p2    47.957      3.4411   13.936  3.812e-44
    ## 
    ## -------------------------------------------------------------------
    ## Sigma link function:  identity
    ## Migma Coefficients:
    ##              Estimate  Std. Error  t-value    p-value
    ## (Intercept)     9.814       2.048    4.792  1.651e-06
    ## 
    ## -------------------------------------------------------------------
    ## No. of observations in the fit:  46 
    ## Degrees of Freedom for the fit:  3
    ##       Residual Deg. of Freedom:  43 
    ##                       at cycle:  22 
    ##  
    ## Global Deviance:     214.695 
    ##             AIC:     220.695 
    ##             SBC:     226.1809 
    ## *******************************************************************

``` r
df <- data.frame(N=N,W=W,fv=mod1$mu.fv)
ggplot(df,aes(W,N)) + geom_point() + geom_line(aes(W,fv),col=2)
```

![](Application_2_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
envelope.nlrbs <- function(model,k=19, color = "grey50", xlabel = "Theorical Quantile",ylabel = "Empirical Quantile",font="serif")
{
 
  n <- model$N
  td  <- model$residuals
  sigma <- model$sigma.fv[1]
  mu <- model$mu.fv
  re <- matrix(0,n,k)

  for(i in 1:k)
  { #inicio do for do residuos
    y1 <- mapply(rRBS,n=1,mu=mu,sigma=sigma)
    model1 <- nlgamlss(y=y1,mu.fo=model$mu.formula,sigma.formula = ~1,mu.start =c(1.276,40.954),sigma.start =9.817,family = RBS(mu.link = "identity",sigma.link = "identity"))
    rdf=model1$residuals
    re[,i]=sort(rdf)
  }
  e10 <- numeric(n)
  e20 <- numeric(n)
  e11 <- numeric(n)
  e21 <- numeric(n)
  e12 <- numeric(n)
  e22 <- numeric(n)
  for(l in 1:n)
  {#inicio do Bootstrap
    eo = sort(re[l,])
    e10[l] = eo[ceiling(k*0.01)]
    e20[l] = eo[ceiling(k*(1 - 0.01))]
    e11[l] = eo[ceiling(k*0.05)]
    e21[l] = eo[ceiling(k*(1 - 0.05))]
    e12[l] = eo[ceiling(k*0.1)]
    e22[l] = eo[ceiling(k*(1 - 0.1))]
  }#fim do Bootstrap
  
  a <- qqnorm(e10, plot.it = FALSE)$x
  r <- qqnorm(td, plot.it = FALSE)$x
  xb = apply(re, 1, mean)
  rxb <- qqnorm(xb, plot.it = FALSE)$x
  
  df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
  ggplot(df,aes(r,td))+geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel)+theme(text=element_text(size=14,family=font)) 
}
envelope.nlrbs(mod1)
```

![](Application_2_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->

``` r
df<- data.frame(res=mod1$residuals,ind=1:mod1$N,fv=mod1$mu.fv)
ggplot(df,aes(ind,res)) + geom_point() + geom_hline(yintercept = c(-3,3),lty=2)
```

![](Application_2_files/figure-gfm/unnamed-chunk-1-4.png)<!-- -->

``` r
ggplot(df,aes(fv,res)) + geom_point() + geom_hline(yintercept = c(-3,3),lty=2)
```

![](Application_2_files/figure-gfm/unnamed-chunk-1-5.png)<!-- -->

``` r
model <- mod1
x1 <- exp(47.954/W)
x2 <- (1.276/W)*exp(47.954/W)
x <- cbind(x1,x2)
z <- matrix(1,nrow=model$N)
y <- model$y
p <- ncol(x)
q <- ncol(z)
linkstr <- "identity"
linkobj <- make.link(linkstr)
linkfun <- linkobj$linkfun
linkinv <- linkobj$linkinv
mu.eta <- linkobj$mu.eta
sigma_linkstr <- "identity"
sigma_linkobj <- make.link(sigma_linkstr)
sigma_linkfun <- sigma_linkobj$linkfun
sigma_linkinv <- sigma_linkobj$linkinv
sigma_mu.eta <- sigma_linkobj$mu.eta
B = function(Delta, I, M) {
  B = (t(Delta) %*% (I - M) %*% Delta)
  return(B)
}
loglik <- function(vP) {
  betab = vP[1:p]
  alpha = vP[-(1:p)]
  eta = as.vector(x %*% betab)
  tau = as.vector(z %*% alpha)
  mu = linkinv(eta)
  sigma = sigma_linkinv(tau)
  f <- 0.5 * sigma - 0.5 * log((sigma + 1)) - 0.5 * log(mu) - 
    1.5 * log(y) + log((sigma * y) + y + (sigma * mu)) - 
    y * (sigma + 1)/(4 * mu) - (sigma * sigma * mu)/(4 * 
                                                       y * (sigma + 1)) - 0.5 * log(16 * pi)
  return(sum(f))
}
muest <- model$mu.coefficients
sigmaest <- model$sigma.coefficients
x0 <- c(muest, sigmaest)
h0 <- hessian(loglik, x0)
Ldelta = h0[(p + 1):(p + q), (p + 1):(p + q)]
Lbeta = h0[1:p, 1:p]
b11 = cbind(matrix(0, p, p), matrix(0, p, q))
b12 = cbind(matrix(0, q, p), solve(Ldelta))
B1 = rbind(b11, b12)
b211 = cbind(solve(Lbeta), matrix(0, p, q))
b212 = cbind(matrix(0, q, p), matrix(0, q, q))
B2 = rbind(b211, b212)
b311 = cbind(matrix(0, p, p), matrix(0, p, q))
b312 = cbind(matrix(0, q, p), matrix(0, q, q))
B3 = rbind(b311, b312)
mu <- model$mu.fv
sigma <- model$sigma.fv
eta <- linkfun(mu)
ai <- mu.eta(eta)
dmu <- (-1/(2 * mu)) + sigma/((y * sigma) + y + (sigma*mu)) + ((sigma + 1) * y)/(4 * (mu^2)) - (sigma^2)/(4 * 
                                                                                                          y * (sigma + 1))
  Deltamu <- crossprod(x, diag(ai * dmu))
  tau <- sigma_linkfun(sigma)
  bi <- sigma_mu.eta(tau)
  dsigma <- (y + mu)/((sigma * y) + y + (sigma * mu)) - 
    y/(4 * mu) - (sigma * (sigma + 2) * mu)/(4 * (sigma + 
                                                    1) * (sigma + 1) * y) + sigma/(2 * (sigma + 1))
  Deltasigma <- crossprod(z, diag(bi * dsigma))
  Delta <- rbind(Deltamu, Deltasigma)
  BT <- B(Delta, solve(h0), B3)
  autovmaxthetaPC <- eigen(BT, symmetric = TRUE)$val[1]
  vetorpcthetaPC <- eigen(BT, symmetric = TRUE)$vec[, 1]
  dmaxG.theta <- abs(vetorpcthetaPC)
  vCithetaPC <- 2 * abs(diag(BT))
  Cb0 <- vCithetaPC
  Cb.theta <- Cb0/sum(Cb0)
  BM <- B(Delta, solve(h0), B1)
  autovmaxbetaPC <- eigen(BM, symmetric = TRUE)$val[1]
  vetorpcbetaPC <- eigen(BM, symmetric = TRUE)$vec[, 1]
  dmaxG.beta <- abs(vetorpcbetaPC)
  vCibetaPC <- 2 * abs(diag(BM))
  Cb1 <- vCibetaPC
  Cb.beta <- Cb1/sum(Cb1)
  BD <- B(Delta, solve(h0), B2)
  autovmaxdeltaPC <- eigen(BD, symmetric = TRUE)$val[1]
  vetordeltaPC <- eigen(BD, symmetric = TRUE)$vec[, 1]
  dmaxG.alpha = abs(vetordeltaPC)
  vCideltaPC = 2 * abs(diag(BD))
  Cb2 = vCideltaPC
  Cb.alpha = Cb2/sum(Cb2)
  local.infl <- list(dmax.beta = dmaxG.beta, dmax.alpha = dmaxG.alpha, 
                 dmax.theta = dmaxG.theta, Ci.beta = Cb.beta, Ci.alpha = Cb.alpha, 
                 Ci.theta = Cb.theta)
 
  Cibeta <<- local.infl$Ci.beta
  Cialpha <<- local.infl$Ci.alpha
  Cigamma <<- local.infl$Ci.gamma
  Citheta <<- local.infl$Ci.theta
  
  dmbeta <<- local.infl$dmax.beta
  dmalpha <<- local.infl$dmax.alpha
  dmtheta <<- local.infl$dmax.theta 
  
  index <<- 1:length(Cibeta)
  CbxI <- data.frame(index,Cibeta,Cialpha,Citheta)
  CbxI. <- CbxI[c(1),]
  p1 <- ggplot(CbxI,aes(index,Cibeta)) + geom_point(pch=4) + geom_text(data=CbxI.,aes(index,Cibeta,label=index),hjust=2) 
  p2 <- ggplot(CbxI,aes(index,Cialpha)) + geom_point(pch=4) + geom_text(data=CbxI.,aes(index,Cialpha,label=index),hjust=2)
  p4 <- ggplot(CbxI,aes(index,Citheta)) + geom_point(pch=4) + geom_text(data=CbxI.,aes(index,Citheta,label=index),hjust=2)
  p1 + ylim(0,1)+xlab("I")
```

![](Application_2_files/figure-gfm/unnamed-chunk-1-6.png)<!-- -->

``` r
  p2 + ylim(0,1)+xlab("I")
```

![](Application_2_files/figure-gfm/unnamed-chunk-1-7.png)<!-- -->

``` r
  p4 + ylim(0,1)+xlab("I")
```

![](Application_2_files/figure-gfm/unnamed-chunk-1-8.png)<!-- -->

``` r
  I <- 1:length(dmbeta)
  dxI <- data.frame(I,dmbeta,dmalpha,dmtheta)
  dxI... <- dxI[c(1),]
  p1 <- ggplot(CbxI,aes(I,dmbeta)) + geom_point(pch=4) 
  p2 <- ggplot(CbxI,aes(I,dmalpha)) + geom_point(pch=4) + geom_text(data=dxI...,aes(I,dmalpha,label=I),hjust=2)
  p4 <- ggplot(CbxI,aes(I,dmtheta)) + geom_point(pch=4) 
  p1 + ylim(0,1)
```

![](Application_2_files/figure-gfm/unnamed-chunk-1-9.png)<!-- -->

``` r
  p2 + ylim(0,1)
```

![](Application_2_files/figure-gfm/unnamed-chunk-1-10.png)<!-- -->

``` r
  p4 + ylim(0,1)
```

![](Application_2_files/figure-gfm/unnamed-chunk-1-11.png)<!-- -->

We have that the mean of the number of cycles (\(\times 100\)) to
failure of the metal specimen can be described by
\(\hat \mu(w_i) = 1.276\times \textrm{e}^{\frac{47.954}{w_i}}\).

## References

<div id="refs" class="references">

<div id="ref-lc:09">

Lemonte, A.J., and G.M. Cordeiro. 2011. “Birnbaum-Saunders nonlinear
regression models.” *Journal of Statistical Computation and Simulation*
53: 4441–52.

</div>

<div id="ref-Lemonte:2016aa">

Lemonte, Artur J., Gauss M. Cordeiro, and Germán Moreno-Arenas. 2016.
“Improved Likelihood-Based Inference in Birnbaum–Saunders Nonlinear
Regression Models.” *Applied Mathematical Modelling* 40 (19): 8185–8200.

</div>

<div id="ref-Lemonte:2011aa">

Lemonte, Artur J., and Alexandre G. Patriota. n.d. “Influence
Diagnostics in Birnbaum–Saunders Nonlinear Regression Models.” *Journal
of Applied Statistics* 38 (5). Taylor & Francis: 871–84.

</div>

<div id="ref-rn:91">

Rieck, J.R., and J.R. Nedelman. 1991. “A log-linear model for the
Birnbaum-Saunders distribution.” *Technometrics* 3: 51–60.

</div>

<div id="ref-vrc:12">

Vanegas, L.H., L.M. Rondon, and F.J.A. Cysneiros. 2012. “Diagnostic
procedures in Birnbaum-Saunders nonlinear regression models.”
*Computational Statistics and Data Analysis* 56: 1662–80.

</div>

</div>
