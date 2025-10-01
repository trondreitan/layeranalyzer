source("numeric_stoch_diff.R")

n=200

mu=3


# Generate timeseries realizations based on 2 processes, one cause and one effect
t1=10
t2=20
s1=2
s2=3
sigma1=s1*sqrt(2/t1)
sigma2=s2*sqrt(2/t2)
rho=0.9

der2=function(x)
{
  ret=c(-(1/t1)*(x[1]-mu),-(1/t2)*(x[2]-mu))
  return(ret)
}

diff2=function(x)
{
  ret=matrix(c(sigma1^2,rho*sigma1*sigma2,rho*sigma1*sigma2,sigma2^2),nrow=2)
  return(ret)
}

x02=rep(3,2)

traj2=explicit.stoch(x02, der2, diff2,100,0.01)
#plot(traj2$t, traj2$x[,1], type="l")
#lines(traj2$t, traj2$x[,2],col="red")

index1=sort(sample(1:length(traj2$t),n,replace=F))
T1=traj2$t[index1]
x1=traj2$x[index1,1]+rnorm(length(T1),0,0.1)


index2=sort(sample(1:length(traj2$t),n,replace=F))
T2=traj2$t[index1]
x2=traj2$x[index1,2]+rnorm(length(T2),0,0.1)


# plot(T1,x1, type="b")
# lines(T2,x2,type="b",col="red")


write(t(cbind(T1,x1)),"test_corr1.txt",ncol=2)
write(t(cbind(T2,x2)),"test_corr2.txt",ncol=2)

