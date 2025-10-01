source("numeric_stoch_diff.R")

n=200

mu=3


# Generate timeseries realizations based on 2 processes, one cause and one effect
t1=3
t2=20
s1=1
s2=4
sigma1=s1*sqrt(2/t1)
sigma2=s2*sqrt(2/t2)

der2=function(x)
{
  ret=c(-(1/t1)*(x[1]-x[2]),-(1/t2)*(x[2]-mu))
  return(ret)
}

diff2=function(x)
{
  ret=diag(c(sigma1^2,sigma2^2))
  return(ret)
}

x02=rnorm(2,mu,s1)

traj2=explicit.stoch(x02, der2, diff2,100,0.01)
#plot(traj2$t, traj2$x[,1], type="l")
#lines(traj2$t, traj2$x[,2],col="red")

index1=sort(sample(1:length(traj2$t),n,replace=F))
t1=traj2$t[index1]
x1=traj2$x[index1,1]+rnorm(length(t1),0,0.1)


index2=sort(sample(1:length(traj2$t),n,replace=F))
t2=traj2$t[index1]
x2=traj2$x[index1,2]+rnorm(length(t2),0,0.1)


# plot(t1,x1, type="b")
# lines(t2,x2,type="b",col="red")

write(t(cbind(t1,x1)),"test_effect.txt",ncol=2)
write(t(cbind(t2,x2)),"test_cause.txt",ncol=2)

