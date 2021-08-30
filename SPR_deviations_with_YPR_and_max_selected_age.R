# calculate differences in SPR

##########################
#### INPUT PARAMETERS ####
##########################
# set working directory and read in the data
setwd('your/path/here')
species.df <- read.csv('species_parameters.csv')

# input species name, and code will pull relevant parameters
# for convenience, here are species names:
# Thunnus albacares, Scomber scombrus, Gadus morhua, Clupea harengus, Sebastes mystinus
sp.name <- "Sebastes mystinus"
sp.index <- which(species.df$Species == sp.name); sp.index

# populate parameters from data
M    <- species.df[sp.index, "M"]; M
vB.Linf <- species.df[sp.index, "Linf"]; vB.Linf
vB.k    <- species.df[sp.index, "k"]; vB.k
vB.t0   <- 0 # can replace with actual values later

# fecundity exponents
b.mass    <- species.df[sp.index, "b.mass"]; b.mass
b.length  <- species.df[sp.index, "b.length"]; b.length # not using this at the moment

# use Bev/Holt equation for L.mat50, could use Marc's if M depends on L
L.mat50 <- with(species.df[sp.index,], Linf * b.mass * k / (M + b.mass * k)); L.mat50

# set length at selectivity equal to a fraction (sel.mult) of L.mat50
sel.mult <- 1 # set equal to 1 for selectivity equal to maturity
L.sel50 <- sel.mult * L.mat50; L.sel50

# weight-length (g and cm), W=aL^b
# hard-wired for now, can include in input file later
WL.a <- 0.01
WL.b <- 3
a.mass <- 1 # coef. of fecundity-mass relationship

# slope parameters for maturity ogive and selectivity
# large values (e.g. 5+) approximate full selection beginning with length in data file (L.mat50 and L.sel50)
mat.slp    <- 10
sel.slp    <- 10

# range of SPR targets as proportion of unfished egg production
SPR.target <- seq(0.95, 0.15, by=-0.05); SPR.target # grid of target SPR targets (F_X%)

# plausible age range for SPR calculations
A.min  <- 0
A.max  <- round(5/M, 0) # a ripe old age
A.plus <- 3*A.max       # approximation of "plus group"

# age after which selectivity is set to zero 
A.sel.max <- A.plus # default; no decline in selectivity with age
A.sel.max <- A.max-20 # for example, ages greater than Amax-10 have selectivity = 0

#################
### FUNCTIONS ###
#################
# von Bert growth
L <- function(A, Linf, k, t0) { Linf * (1-exp(-k*(A-t0))) }

# weight-length function
W <- function(L, W.a, W.b) { W.a * L^W.b }

# maturity at length function
Mat <- function(L, slope, L.50) { 1/(1+exp(-slope*(L-L.50))) }

# selectivity at length function
Sel <- function(L, slope, L.50) { 1/(1+exp(-slope*(L-L.50))) }

# fecundity-weight function (eggs per gram)
Fec <- function(W,a,b) { a * W^b }

# SPR function takes scalar F & M, plus vectors of sel.at.age, mat.at.age, and fec.at.age
SPR <- function(F, M, age.vec, sel.at.age, mat.at.age, fec.at.age)
{
   unfished <- sum( exp(-M*age.vec) * mat.at.age * fec.at.age )
   Z <- numeric(length(age.vec)) # Z is cumulative mortality by age a
   for (i in 2:length(age.vec))
   {
      Z[i] <- Z[i-1] + M + F*sel.at.age[i-1]
   }
   fished <- sum( exp(-Z) * mat.at.age * fec.at.age )
   spr <- fished/unfished
   return(spr)
}

# objective function to find F associated with SPR target
obj <- function(F.tmp, SPR.targ, M.tmp, a.vec, sel.vec, mat.vec, fec.vec) {
   spr.tmp <- SPR(F = F.tmp, M = M.tmp, age.vec=a.vec, sel.at.age=sel.vec, mat.at.age=mat.vec, fec.at.age=fec.vec)
   sq.error <- (spr.tmp - SPR.targ)^2
   return(sq.error)
}

YPR <- function(F, M, age.vec, sel.at.age, w.at.age)
{
   F.at.age <- F*sel.at.age
   Z.at.age <- M + F.at.age
   cumulative.Z <- numeric(length(age.vec))
   for (i in 2:length(age.vec))
   {
      cumulative.Z[i] <- sum( M + F.at.age[1:i] )
   }
   ypr <- sum( exp(-cumulative.Z) * (F.at.age/Z.at.age) *
              (1 - exp(-cumulative.Z)) * w.at.age )
   return(ypr)
}

#########################
### PRELIMINARY CALCS ###
#########################
a.vec	<- A.min:A.plus	                              # age vector
l.vec	<- L(A=a.vec, Linf=vB.Linf, k=vB.k, t0=vB.t0) # length-at-age vector
w.vec	<- W(L=l.vec, W.a=WL.a, W.b=WL.b)             # weight-at-age vector
fec.vec <- Fec(w.vec, a=a.mass, b=b.mass)              # fecundity-at-age vector
mat.vec <- Mat(L=l.vec, slope=mat.slp, L.50=L.mat50)	# maturity-at-age vector
sel.vec <- Sel(L=l.vec, slope=sel.slp, L.50=L.sel50)	# selectivity-at-age vector

# make selectivity at age = 0 for ages greater than A.sel.max
index <- which(a.vec == A.sel.max)
if (A.sel.max < A.plus)
{
   sel.vec[(index+1):length(sel.vec)] <- 0
}

# CHECK: calculate SPR with F=0 (should be 1)
SPR(0, M=M, age.vec=a.vec, sel.at.age=sel.vec, mat.at.age=mat.vec, fec.at.age=fec.vec)

# write csv file with age-specific quantities
age.table <- cbind.data.frame(Age=a.vec, Length=l.vec, Weight=w.vec, Fecundity=fec.vec, Maturity=mat.vec, Selectivity=sel.vec)
head(age.table)
write.csv(age.table, file=paste(sp.name,' age table.csv',sep=''), row.names=FALSE, quote=FALSE)

####################################
### CALCULATE DIFFERENCES IN SPR ###
####################################

# create object to hold results
out <- matrix(NA, nrow=length(SPR.target), ncol=5)

# loop over target SPR values
# search for F between 0 and F.max
F.max <- 5

for (i in 1:length(SPR.target)) {
   # record the actual SPR target
   SPR.targ <- SPR.target[i]; SPR.targ
   # find F for current target SPR given b.mass=1
   soln <- optimize(obj, M.tmp=M, SPR.targ=SPR.target[i], fec.vec=Fec(w.vec, a=a.mass, b=1), 
                    a.vec=a.vec, sel.vec=sel.vec, mat.vec=mat.vec, interval=c(0,F.max))
   soln
   F.SPR <- soln[[1]]; F.SPR
   # flag solutions where obj function is greater than (0.0005)^2
   if(soln[[2]] > 2.5e-07) { print(paste("Solution to iteration",i,"has squared error greater than 2.5e-07")) }
   # plug F.SPR back into the SPR formula to be sure that it gives you the right answer
   SPR.c <- SPR(F.SPR, M=M, age.vec=a.vec, sel.at.age=sel.vec, mat.at.age=mat.vec, fec.at.age=Fec(w.vec, a=a.mass, b=1))
   SPR.c
   # plug F.SPR into the SPR formula with fecundity at age based on b.mass != 1
   SPR.real <- SPR(F.SPR, M=M, age.vec=a.vec, sel.at.age=sel.vec, mat.at.age=mat.vec, fec.at.age=fec.vec)
   SPR.real
   # calculate difference in SPR
   SPR.Diff <- SPR.c - SPR.real; SPR.Diff
   out[i,] <- c(SPR.targ, SPR.c, F.SPR, SPR.real, SPR.Diff)
}
out.df <- as.data.frame(out)
names(out.df) <- c('SPR.target', 'SPR.est', 'F.SPR.b1', 'SPR.true', 'SPR.diff')
out.df

# solve for F that reduces SPR to target based on 'real' b.mass (not b.mass=1, like above)
# create object to hold results
out2 <- matrix(NA, nrow=length(SPR.target), ncol=3)

for (i in 1:length(SPR.target)) {
   # record the actual SPR target
   SPR.targ <- SPR.target[i]; SPR.targ
   # find F for current target SPR given b.mass=1
   soln <- optimize(obj, M.tmp=M, SPR.targ=SPR.target[i], fec.vec=fec.vec, # fec.vec is based on b.mass from input species
                    a.vec=a.vec, sel.vec=sel.vec, mat.vec=mat.vec, interval=c(0,F.max))
   soln
   F.SPR <- soln[[1]]; F.SPR
   # flag solutions where obj function is greater than (0.0005)^2
   if(soln[[2]] > 2.5e-07) { print(paste("Solution to iteration",i,"has squared error greater than 2.5e-07")) }
   # plug F.SPR back into the SPR formula to be sure that it gives you the right answer
   SPR.c <- SPR(F.SPR, M=M, age.vec=a.vec, sel.at.age=sel.vec, mat.at.age=mat.vec, fec.at.age=fec.vec)
   SPR.c
   out2[i,] <- c(SPR.targ, SPR.c, F.SPR)
}
out2 # columns 1 and 2 should be very close (+/- 0.001)

# append the F that achieves target SPR with b.mass > 1
out.df$F.SPR.bmass <- out2[,3]
out.df

# create objects to hold YPR with b=1 and b.mass for input species
YPR.F.b1 <- YPR.F.bmass <- numeric(length(SPR.target))
for (i in 1:length(SPR.target)) {
   YPR.F.b1[i]    <- YPR(F=out.df[i,"F.SPR.b1"], M=M, age.vec=a.vec, sel.at.age=sel.vec, w.at.age=w.vec)
   YPR.F.bmass[i] <- YPR(F=out.df[i,"F.SPR.bmass"], M=M, age.vec=a.vec, sel.at.age=sel.vec, w.at.age=w.vec)
}

# add YPR to output, including YPR(SPR.target, b=1) - YPR(SPR.target, b>1)
out.df <- cbind(out.df, YPR.F.b1, YPR.F.bmass, YPR.diff=YPR.F.b1-YPR.F.bmass)
out.df

# plot
windows()
par(mfrow=c(2,2))
# SPR vs. F
with(out.df, plot(F.SPR.b1, SPR.target, type='o', xlim=c(0,1.2*max(out.df[,c(3,6)])), ylim=c(0,1), ylab="SPR"))
# YPR vs. F
with(out.df, plot(F.SPR.b1, YPR.F.b1, type='o', xlim=c(0,1.2*max(out.df[,c(3,6)])), ylim=c(0,1.2*max(out.df[,c(7,9)])), ylab="YPR"))
# plot YPR as a function of the SPR targets
# first, plot YPR assuming b.mass=1 (in black)
with(out.df, plot(SPR.target, YPR.F.b1, type='o', xlim=c(0,1), ylim=c(0,1.2*max(out.df[,c(7,9)])), ylab="YPR"))
# now plot YPR using b.mass for input species
with(out.df, points(SPR.target, YPR.F.bmass, pch=20, col='red'))
with(out.df, lines(SPR.target, YPR.F.bmass, col='red'))
legend("topleft", c("F(b=1)","F(b)"), col=1:2, bty="n", text.col=1:2)
plot(a.vec[1:(A.max+1)], sel.vec[1:(A.max+1)], type='l', ylab='selectivity', xlab='age', lwd=2)
lines(a.vec[1:(A.max+1)], mat.vec[1:(A.max+1)], col='red')
text(0.75*A.max, 0.55, "Selectivity")
text(0.75*A.max, 0.45, "Maturity", col='red')

# Since YPR has a maximum, YPR(b=1) produces higher yields at low F (high SPR.target)
# and lower yields for high F (low SPR.target)
#### However, this is PER RECRUIT, and does not account for reductions in recruitment resulting from high F

write.csv(out.df, file=paste(sp.name,' SPR_differences_and_YPR.csv',sep=''), row.names=FALSE, quote=FALSE)
