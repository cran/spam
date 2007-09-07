# This script should not be run... 

# plotting
if (F) { # test ploting is a bad idea...


n <- 10
m <- 5

set.seed(124)
tt <- matrix(rnorm(m*n),n,m)
tt[tt<0] <- 0
ss <- as.spam(tt)



par(mfcol=c(1,2))
.Spam$imagesize=10
display(ss)
.Spam$imagesize=n*m+1
display(ss)



par(mfcol=c(1,2))
plot(tt)
plot(ss)


plot(      tt[,1])
plot(spam( tt[,1]))

plot(      t( tt[,2]))
plot(spam( t( tt[,2]),nrow=1))



nl <- length(ss)  #ok
ss@entries <- 1:nl
z <- ss
br <- c(seq(0.1,max(z)/2,l=nl),max(z))
par(mfcol=c(1,2))
.Spam$imagesize=1000
image(z, breaks=br,col=tim.colors(nl))
.Spam$imagesize=10
image(z, breaks=br,col=tim.colors(nl))


nl <- length(ss)
ss@entries <- 1:nl
par(mfcol=c(1,2))
.Spam$imagesize=1000
image(ss, col=tim.colors(nl))
.Spam$imagesize=10
image(ss, col=tim.colors(nl))


# very large sample
nz <- 128
ln <- nz^2
ss <- spam(0,ln,ln)
for (i in 1:nz) ss[sample(ln,1),sample(ln,1)] <- i
image(ss, col=tim.colors(nl),cex=10)




}
