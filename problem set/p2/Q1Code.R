rm(list=ls())

# 1 (c)
diceDP=function(n=1,pi1=0.5,pi2=0.25,pi3=1/8,pi4=1/16,pi5=1/32,pi6=1/32){
    if (n == 1){
        return(c(pi1,pi2,pi3,pi4,pi5,pi6))
    }
    
    # n = 1, initialize table to store result
    table=matrix(rep(0,n*6*n), nrow=n)
    table[1,1:6]=c(pi1,pi2,pi3,pi4,pi5,pi6)
    
    # using i to indicate how many times to roll the dice, j indicate the result each roll
    for (i in seq(2,n,1)){
        for (y in seq(i,6*n,1)){
            for (k in 1:(min(7,y)-1)){
                table[i,y] = table[i,y] + table[i-1,y-k] * table[1,k]
            }
         
        }
    }
    vector_re = table[n,seq(n,6*n,1)]
    names(vector_re)=seq(n,6*n,1)
    return(vector_re)
}

# compute f_Y50
Y50 = diceDP(n=50)

# 10 highest and lowest probability outcomes
sortY50 = sort(Y50)
# 10 lowest 
cbind(names(sortY50[1:10]),sortY50[1:10])
# 10 highest
cbind(names(sortY50[length(sortY50):(length(sortY50)-9)]),sortY50[length(sortY50):(length(sortY50)-9)])


####################################
# 1 (d)
pi1=0.5;pi2=0.25;pi3=1/8;pi4=1/16;pi5=1/32;pi6=1/32

Ex = 1*pi1 + 2*pi2 + 3*pi3 + 4*pi4 + 5*pi5 + 6*pi6
Ex2= 1^2*pi1 + 2^2*pi2 + 3^2*pi3 + 4^2*pi4 + 5^2*pi5 + 6^2*pi6
Vx = Ex2 - Ex^2

# negative binomial approximation
diceNB = function(n,phi){
    vector_re = NULL
    for (i in seq(n,6*n,1)){
        pro = choose((i-1),(n-1))*(1-phi)^(i-n)*phi^n
        vector_re = c(vector_re,pro)
    }
    names(vector_re)=seq(n,6*n,1)
    return(vector_re)
}

# normal approximation
diceNormal = function(n,Ex,Vx){
    E = n*Ex
    V = n*Vx
    vector_re = NULL
    for (i in seq(n,6*n,1)){
        pro = pnorm(i, mean=E, sd=sqrt(V)) - pnorm(i-1, mean=E, sd=sqrt(V))
        vector_re = c(vector_re,pro)
    }
    names(vector_re)=seq(n,6*n,1)
    return(vector_re)
}

DP50=diceDP(n=50)
NB50=diceNB(50,phi=1/Ex)
Normal50=diceNormal(50,Ex,Vx)

plot(as.numeric(names(DP50)),DP50,pch=20,main="Problem1 (d)",ylab="Probability",xlab="y")
points(as.numeric(names(DP50)),NB50,col="red",pch=3)
points(as.numeric(names(DP50)),Normal50,col="green",pch=4)
legend('topright',c("Dynamice Programming","Negative Binomial","Normal"),pch=c(20,3,4),col=c("black","red","green"))

# variance
# variance of negative binomial approxiamation
n=50; phi=1/Ex
V_NB=n*(1-phi)/phi^2
V_NB

# variance of normal approximation
V_Normal=n*Vx
V_Normal

# variance of dynamice programming
Ex_DP = sum(seq(n,6*n,1) * DP50)
Ex_DP_square = sum(seq(n,6*n,1)^2 * DP50)
V_DP = Ex_DP_square - Ex_DP^2
V_DP

#######################
# 1 (e)
DP50=diceDP(n=50)
# p value of dynamic programming
DP50['300']

# exact number of p value is 
pi6^50

# p value of normal approximation
n=50
pnorm(300,mean=n*Ex,sd=sqrt(n*Vx),lower.tail = F)

# p value of binomial approximation
pnbinom(6*n-n,n,1/Ex,lower.tail = F)
