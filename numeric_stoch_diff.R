
# Code based on description in Kloeden & Platen (1992):
# Numerical Solution of Stochastic Differential Equations.
# Springer Verlag, Berlin.


explicit.stoch=function(x0, # initial state
  state.der, # vector function of the state, deterministic time derivative
  diffusion, # diffusion matrix function of state
  T, # end of time interval
  dt) # time steps
  {
    t=seq(0,T,dt) # time vector
    N=length(t) # length of time interval
    n=length(x0) # size of state vector
    x=array(0,c(N,n)) # state matrix (time and state size)
    one=rep(1,n)
    
    x[1,]=x0 # start the simulation
    
    for(i in 2:N)
      {
        W=rnorm(n)*sqrt(dt)
        x.old=x[i-1,]

        x1=x.old+state.der(x.old)*dt+diffusion(x.old)%*%W

        R.plus=(x.old+state.der(x.old)*dt)%*%t(one)+diffusion(x.old)*sqrt(dt)
        R.minus=(x.old+state.der(x.old)*dt)%*%t(one)-diffusion(x.old)*sqrt(dt)
        U.plus=x.old%*%t(one)+diffusion(x.old)*sqrt(dt)
        U.minus=x.old%*%t(one)-diffusion(x.old)*sqrt(dt)

        second.term=0*x0
        for(j in 1:n)
          {
            for(r in 1:n)
              {
                if(r==j)
                  {
                    second.term=second.term+
                      (diffusion(R.plus[,j])[,j]+
                       diffusion(R.minus[,j])[,j]+
                       2*diffusion(x.old)[,j])*W[j]
                  }
                else
                  {
                    second.term=second.term+
                      (diffusion(U.plus[,r])[,j]+
                       diffusion(U.minus[,r])[,j]-
                       2*diffusion(x.old)[,j])*W[j]
                  }
              }
          }
        second.term=second.term*0.25

        third.term=0
        V=array(0,c(n,n))
        for(j1 in 1:n)
          {
            V[j1,j1]=dt
            if(j1<n)
              {
                for(j2 in (j1+1):n)
                  {
                    if(runif(1)<0.5)
                      {
                        V[j1,j2]=dt
                        V[j2,j1]=-dt
                      }
                    else
                      {
                        V[j1,j2]=-dt
                        V[j2,j1]=dt
                      }
                  }
              }
          }
        for(j in 1:n)
          {
            for(r in 1:n)
              {
                if(r==j)
                  {
                    third.term=third.term+
                      (diffusion(R.plus[,j])[,j]-
                       diffusion(R.minus[,j])[,j])*
                         (W[j]^2-dt)/sqrt(dt)
                  }
                else
                  {
                    third.term=third.term+
                      (diffusion(U.plus[,r])[,j]-
                       diffusion(U.minus[,r])[,j])*
                         (W[j]*W[r]+V[r,j])/sqrt(dt)
                  }
              }
          }
        third.term=third.term*0.25
        
        x[i,]=x.old+0.5*(state.der(x.old)+state.der(x1))*dt+
          second.term+third.term
      }

    ret=list(t=t, x=x)

    return(ret)
  }

