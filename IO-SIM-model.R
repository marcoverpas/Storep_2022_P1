#This R code reproduces some of the experiments discussed at STOREP Conference,
#May 26 2022, Session 3: "Economic modelling". Presentation title:
#Circular economy in a simplified input-output stock-flow consistent model

#Created on May 26 2022 by Marco Veronese Passarella

#UPLOAD LIBRARIES
library(expm)

#CLEAR
rm(list=ls(all=TRUE))

#PERIODS (i= 1 to 65)
nPeriods=65  

#PARAMETERS
alpha1 = 0.6  #Propensity to consume out of income
alpha2 = 0.4  #Propensity to consume out of wealth
theta = 0.2   #Tax rate
mu10 = 0.1    #Autonomous component of mark-up on price of product 1
mu11 = 0.75   #Output-gap elasticity of price of product 1
mu20 = 0.1    #Autonomous component of mark-up on price of product 2
mu21 = 0.25   #Output-gap elasticity of price of product 2
r = 0.02      #Interest rate
lambda0 = 0.1 #Fixed proportion of government bills to total wealth
lambda1 = 0   #Sensitivity of portfolio choices to interest rate
lambda2 = 0   #Sensitivity of portfolio choices to transactions demand for money

#VARIABLES
#Gross product
y=matrix(data=0,nrow=1,ncol=nPeriods)
#Net product
yn=matrix(data=0,nrow=1,ncol=nPeriods)
#Total consumption
c=matrix(data=0,nrow=1,ncol=nPeriods)
#Government expenditures
g=matrix(data=0,nrow=1,ncol=nPeriods) 
#Taxes
t=matrix(data=0,nrow=1,ncol=nPeriods)
#Disposable income
yd=matrix(data=0,nrow=1,ncol=nPeriods)
#Cash demand
h_h=matrix(data=0,nrow=1,ncol=nPeriods)
#Cash supply
h_s=matrix(data=0,nrow=1,ncol=nPeriods)
#Supply of bills
b_s=matrix(data=0,nrow=1,ncol=nPeriods)
#Private demand for bills
b_h=matrix(data=0,nrow=1,ncol=nPeriods)
#CB holdings of bills
b_cb=matrix(data=0,nrow=1,ncol=nPeriods)
#Net wealth of households
v=matrix(data=0,nrow=1,ncol=nPeriods)
#Labour 
n=matrix(data=0,nrow=1,ncol=nPeriods) 
#Product per worker
pr=matrix(data=c(1.2,0.8),nrow=2,ncol=nPeriods)
#Additional matrix to calculate labour coefficients
I2=matrix(data=c(1,1),nrow=2,ncol=nPeriods)
#Wage rate
w=matrix(data=0.86,nrow=1,ncol=nPeriods)
#Vector of mark-up
mu=matrix(data=c(mu10,mu20),nrow=2,ncol=nPeriods)
#Matrix of coefficients
A=matrix(data=c(0.1,0.1,0.1,0.1),nrow=2,ncol=2)
#Auxiliary matrix of coefficients
B=matrix(data=c(0,0,0,0),nrow=2,ncol=2)
#Identity matrix
I <- diag(2)
#Production vector
x=matrix(data=0,nrow=2,ncol=nPeriods)
#Fully-adjusted production vector
x_star=matrix(data=0,nrow=2,ncol=nPeriods)
#Final demand vector
d=matrix(data=0,nrow=2,ncol=nPeriods)
#Price vector
p=matrix(data=c(0.9,1.1),nrow=2,ncol=nPeriods)
#Share of government consumption
sigma=matrix(data=c(0.4,0.6),nrow=2,ncol=nPeriods)
#Share of households consumption
beta=matrix(data=c(0.6,0.4),nrow=2,ncol=nPeriods)
#General price level
pa=matrix(data=0,nrow=1,ncol=nPeriods)
#Government spending deflator
pg=matrix(data=0,nrow=1,ncol=nPeriods)

#MODEL

#Define time loop
for (i in 2:nPeriods){

  #Define iterations loop (to force convergence to simultaneous solution)  
  for (iterations in 1:20){
    
    #Introduce shock to government spending
    if (i>=15){g[1,i]=20} #Government expenditures passes from 0 to 20 after 15 periods    
    
    #Create the model (system of difference equations)
    
    #Total consumption
    if (i<=2) {c[,i] = alpha1*yd[,i-1] + alpha2*v[,i-1]}
    else{c[,i] = alpha1*yd[,i-1]/pa[,i-1] + alpha2*v[,i-1]/pa[,i-1]}
    
    #Final demand/consumption vector
    d[,i] = beta[,i]*c[,i] + sigma[,i]*g[,i] 
    
    #Gross production vector (full adjustment)
    x_star[,i] = solve(I-A) %*% d[,i]
    
    #Gross product in nominal terms
    y[,i] = t(p[,i]) %*% (I+B) %*% d[,i] #t(p[,i]) %*% x[,i]
    
    #Net product in nominal terms
    yn[,i] = t(p[,i]) %*% d[,i]
    
    #Endogenous price vector (using Hadamard product)
    p[1,i]  =  (p[2,i]*A[3]*(1+mu[1,i]) + w[,i]/pr[1,i]) /(1 - A[1]*(1+mu[1,i]) )
    p[2,i]  =  (p[1,i]*A[2]*(1+mu[2,i]) + w[,i]/pr[2,i]) /(1 - A[4]*(1+mu[2,i]) )
    
    #Endogenous mark-up
    mu[1,i] = mu10 + mu11*(x_star[1,i-1]-x[1,i-1])
    mu[2,i] = mu20 + mu21*(x_star[2,i-1]-x[2,i-1])
    
    #General price level faced by households and firms
    pa[,i] = t(p[,i]) %*% beta[,i]
    
    #Deflator of government spending 
    pg[,i] = t(p[,i]) %*% sigma[,i]
    
    #Tax payments
    t[,i] = theta*(yn[,i] + r * b_h[,i-1])
    
    #Disposable income
    yd[,i] = yn[,i] + r * b_h[,i-1] - t[,i]
    
    #Supply of government bills
    b_s[,i] = b_s[,i-1] + g[,i]*pg[,i] - t[,i] + r * b_h[,i-1]
    
    #CB holdings of bills
    b_cb[,i] = b_s[,i] - b_h[,i]
    
    #Supply of cash money
    h_s[,i] = h_s[,i-1] + (b_cb[,i] - b_cb[,i-1])
    
    #Net wealth of households
    v[,i] = v[,i-1] + yd[,i] - c[,i]*pa[,i]              
    
    #Private demand for bills
    b_h[,i] = lambda0 * v[,i] + lambda1 * r * v[,i] - lambda2 * yd[,i]

    #Cash held by households
    h_h[,i] = v[,i] - b_h[,i]
    
  }
  
  #OUT-OF-ITERATION CALCULATIONS
  
  #Gross production vector (adjustment over time)
  if(i>=15){B = B + A%^%(i-14)}
  else{B = c(0,0,0,0) }
  x[,i] = (I+B) %*% d[,i]
  
  #Employment 
  n[,i] = (I2[,i]/pr[,i]) %*% x[,i]
  
}


#FIGURE 2
layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))

#Consisency check
plot(h_h[1,2:65]-h_s[1,2:65],type="l",col=2,lwd=2,lty=1,font.main=1.5,cex.main=1.5, main="a) Consistency check",ylab = '£',xlab = '',cex.axis=1.5,cex.lab=1.5,ylim=range(-0.1,0.1))

#Output components: nominal
plot(yn[,2:65],type="l",col=1,lwd=2,lty=1,font.main=1.5,cex.main=1.5, main="b) Output components (nominal)",ylab = '£',xlab = '',cex.axis=1.5,cex.lab=1.5,ylim=range(0,170))
lines(y[,2:65],type="l",col=2,lwd=2,lty=2)
lines(c[,2:65]*pa[,2:65],type="l",col=3,lwd=2,lty=3)
lines(g[,2:65]*pa[,2:65],type="l",col=4,lwd=2,lty=4)
legend("right",c("Net output","Gross output","Households consumption","Government consumption"),  bty = 'n', cex=1.5, lty=c(1,2,3,4), lwd=c(2,2,2,2), col = c(1,2,3,4), box.lty=0)

#Output components: real
plot(yn[,2:65]/pa[,2:65],type="l",col=1,lwd=2,lty=1,font.main=1.5,cex.main=1.5, main="c) Output components (real)",ylab = '(real)',xlab = '',cex.axis=1.5,cex.lab=1.5,ylim=range(0,150))
lines(y[,2:65]/pa[,2:65],type="l",col=2,lwd=2,lty=2)
lines(c[,2:65],type="l",col=3,lwd=2,lty=3)
lines(g[,2:65],type="l",col=4,lwd=2,lty=4)
legend("right",c("Net output","Gross output","Households consumption","Government consumption"),  bty = 'n', cex=1.5, lty=c(1,2,3,4), lwd=c(2,2,2,2), col = c(1,2,3,4), box.lty=0)

#Consumption types
plot(beta[1,2:65]*c[,2:65],type="l",col=3,lwd=2,lty=1,font.main=1.5,cex.main=1.5, main="d) Final demand composition (real)",ylab = '(real)',xlab = '',cex.axis=1.5,cex.lab=1.5,ylim=range(0,75))
lines(beta[2,2:65]*c[,2:65],type="l",col=3,lwd=2,lty=2)
lines(sigma[1,2:65]*g[,2:65],type="l",col=4,lwd=2,lty=1)
lines(sigma[2,2:65]*g[,2:65],type="l",col=4,lwd=2,lty=2)
legend("right",c("Good 1: households","Good 2: households","Good 1: government","Good 2: government"),  bty = 'n', cex=1.5, lty=c(1,2,1,2), lwd=c(2,2,2,2), col = c(3,3,4,4), box.lty=0)

#Disposable income and consumption
plot(yd[,2:65],type="l",col="purple",lwd=2,lty=1,font.main=1.5,cex.main=1.5, main="e) Disposable income (nominal)",ylab = '£',xlab = '',cex.axis=1.5,cex.lab=1.5,ylim=range(0,120))
lines(pa[,2:65]*c[,2:65],type="l",col=3,lwd=2,lty=2)
legend("right",c("Disposable income","Household consumption"),  bty = 'n', cex=1.5, lty=c(1,2), lwd=c(2,2), col = c("purple",3), box.lty=0)

#Wealth and saving
plot(v[,2:65],type="l",lwd=2,lty=1,col="dodgerblue3",font.main=1.5,cex.main=1.5,main="f) Wealth level and saving (nominal)",ylab = '£',xlab = '',cex.axis=1.5,cex.lab=1.5)
par(new="TRUE")
plot(diff(v[,2:65]),type="l",lwd=2,lty=1,col="orange",xlab = '',ylab = '',xaxt='n',yaxt='n')
axis(side=4)
par(xpd=TRUE)
legend("right",c("Stock of financial assets","Household saving (right axis)"),  bty = 'n', cex = 1.5, lty=c(1,1), lwd=c(2,2), col = c("dodgerblue3","orange"), box.lwd=0)


#FIGURE 3
layout(matrix(c(1,2), 1, 2, byrow = TRUE))

#Output gap
plot(x_star[1,13:25]-x[1,13:25],type="l",col=1,lwd=2,lty=1,font.main=1,cex.main=1, main="a) Output gap",ylab = '(real)',xlab = '',cex.axis=1,cex.lab=1,ylim=range(-0.6,0.6))
lines(-x_star[2,13:25]+x[2,13:25],type="l",col=2,lwd=2,lty=1)
legend("topright",c("Product 1","Product 2 (reversed sign)"),  bty = 'n', cex=1, lty=c(1,1), lwd=c(2,2), col = c(1,2), box.lty=0)

#Prices
plot(p[1,13:25]/p[1,13],type="l",col=1,lwd=2,lty=1,font.main=1,cex.main=1, main="b) Unit prices",ylab = '£ (normalised)',xlab = '',cex.axis=1,cex.lab=1,ylim=range(0.95,1.15))
lines(p[2,13:25]/p[2,13],type="l",col=2,lwd=2,lty=1)
lines(pa[,13:25]/pa[,13],type="l",col=3,lwd=2,lty=2)
legend("right",c("Product 1","Product 2","Price level"),  bty = 'n', cex=1, lty=c(1,1,2), lwd=c(2,2,2), col = c(1,2,3), box.lty=0)