source("numeric_stoch_diff.R")

n=400

mu=3

# Generate timeseries realizations based on OU
t1=10
s1=2
sigma1=s1*sqrt(2/t1)

der1=function(x)
{
  ret=-(1/t1)*(x-mu)
  return(ret)
}

diff1=function(x)
{
  ret=as.matrix(sigma1^2)
  return(ret)
}

x01=rnorm(1,mu,s1)

traj1=explicit.stoch(x01, der1, diff1,100,0.01)
#plot(traj1$t, traj1$x[,1], type="l")

t=traj1$t
x1=traj1$x[,1]
y=x1+rnorm(length(x1),0,0.1)

#plot(t,y, type="l")

index=sort(sample(1:length(t),n,replace=F))

write(t(cbind(t[index],y[index])),"test_1layer.txt",ncol=2)
write(t(cbind(t,x1)),"test_1layer_orig1.txt",ncol=2)




# Generate timeseries realizations based on 2 layers
t1=3
t2=20
s1=0.5
s2=2
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

t=traj2$t
x2=traj2$x[,1]
x22=traj2$x[,2]
y=x2+rnorm(length(x2),0,0.1)

#plot(t,y, type="l")

index=sort(sample(1:length(t),n,replace=F))

write(t(cbind(t[index],y[index])),"test_2layer.txt",ncol=2)
write(t(cbind(t,x2)),"test_2layer_orig1.txt",ncol=2)
write(t(cbind(t,x22)),"test_2layer_orig2.txt",ncol=2)


png("2layer.png",height=800,width=1000)
plot(t,x2,type="l",ylim=c(-14,7),lwd=2)
lines(t,x22-10,lwd=2)
dev.off()


png("poisson_2layer.png",height=800,width=1000)
par(cex=1.5)
counts=rpois(length(index),exp(y[index]-3))
plot(t[index],counts,log="y",xlab="time",ylab="counts",lwd=3)
dev.off()




# Generate timeseries realizations based on 3 layers
t1=2
t2=10
t3=50
s1=1
s2=3
s3=5
sigma1=s1*sqrt(2/t1)
sigma2=s2*sqrt(2/t2)
sigma3=s3*sqrt(2/t3)

der3=function(x)
{
  ret=c(-(1/t1)*(x[1]-x[2]),-(1/t2)*(x[2]-x[3]),-(1/t3)*(x[3]-mu))
  return(ret)
}

diff3=function(x)
{
  ret=diag(c(sigma1^2,sigma2^2,sigma3^2))
  return(ret)
}

x03=rnorm(3,mu,s1)

traj3=explicit.stoch(x03, der3, diff3,100,0.01)
#plot(traj3$t, traj3$x[,1], type="l")
#lines(traj3$t, traj3$x[,2],col="red")
#lines(traj3$t, traj3$x[,3],col="blue")

t=traj3$t
x3=traj3$x[,1]
x32=traj3$x[,2]
x33=traj3$x[,3]
y=x3+rnorm(length(x3),0,0.1)

#plot(t,y, type="l")

index=sort(sample(1:length(t),n,replace=F))

write(t(cbind(t[index],y[index])),"test_3layer.txt",ncol=2)
write(t(cbind(t,x3)),"test_3layer_orig1.txt",ncol=2)
write(t(cbind(t,x32)),"test_3layer_orig2.txt",ncol=2)
write(t(cbind(t,x33)),"test_3layer_orig3.txt",ncol=2)




