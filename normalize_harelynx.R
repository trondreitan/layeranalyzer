# Normalizes log(hare) and log(lynx) in files "hare.txt" and "lynx.txt":

# Read data:
hare.orig=read.table("http://folk.uio.no/trondr/layered/hare.txt",sep=" ")
lynx.orig=read.table("http://folk.uio.no/trondr/layered/lynx.txt",sep=" ")
names(hare.orig)=c("year","val")
names(lynx.orig)=c("year","val")

hare.norm=hare.orig
lynx.norm=lynx.orig

# Log-transform:
hare.norm$val=log(hare.orig$val)
lynx.norm$val=log(lynx.orig$val)

# Test normality
shapiro.test(hare.norm$val)
#no
shapiro.test(lynx.norm$val)
#maybe?


# Transform hare data:
y=hare.norm$val

# Make cumulative density approximation:
# First, determine accuracy from the range of the data.
scale=10^ceiling(log10(max(y)-min(y)))
# Density approximation from R's kernel density estimation method:
dd.hare=density(y, n=100001, from=min(y)-scale, to=max(y)+scale, adjust=1) 
dd.hare$F=cumsum(dd.hare$y)*(dd.hare$x[2]-dd.hare$x[1])
dd.hare$F=dd.hare$F/max(dd.hare$F)

# Transform values to something standard normalized-like:
trans.hare=function(numhares)
  qnorm(dd.hare$F[max(which(dd.hare$x<numhares))])

# Transform values back again:
invtrans.hare=function(normhares)
  dd.hare$x[min(which(dd.hare$F>pnorm(normhares)))]

# Transform data:
x=0*y
for(i in 1:length(x))
  x[i]=trans.hare(y[i])
# Put into table:
hare.norm[,2]=x


# Try transforming back:
z=0*y
for(i in 1:length(x))
  z[i]=invtrans.hare(x[i])
# Comparison with x goes well...

# Write transformed hare data to file:
write(t(as.matrix(hare.norm)),file="hare_norm.txt",ncolumns=2)



# Transform lynx data:
y=lynx.norm$val

# Make cumulative density approximation:
# First, determine accuracy from the range of the data.
scale=10^ceiling(log10(max(y)-min(y)))
# Density approximation from R's kernel density estimation method:
dd.lynx=density(y, n=100001, from=min(y)-scale, to=max(y)+scale, adjust=1) 
dd.lynx$F=cumsum(dd.lynx$y)*(dd.lynx$x[2]-dd.lynx$x[1])
dd.lynx$F=dd.lynx$F/max(dd.lynx$F)

# Transform values to something standard normalized-like:
trans.lynx=function(numlynx)
  qnorm(dd.lynx$F[max(which(dd.lynx$x<numlynx))])

# Transform values back again:
invtrans.lynx=function(normlynx)
  dd.lynx$x[min(which(dd.lynx$F>pnorm(normlynx)))]

# Transform values to something standard normalized-like:
x=0*y
for(i in 1:length(x))
  x[i]=trans.lynx(y[i])
# Put into table:
lynx.norm[,2]=x

# Write transformed lynx data to file:
write(t(as.matrix(lynx.norm)),file="lynx_norm.txt",ncolumns=2)

