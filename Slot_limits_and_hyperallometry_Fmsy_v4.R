# Evaluating the performance of slot limits assuming hyperallometric fecundity

# clear and set working directory
rm(list=ls(all=TRUE))
# NEED TO UPDATE THIS LINE WITH YOUR WORKING DIRECTORY
setwd('your/path/here')

# options
Charnov.M <- TRUE ###NEW### replaces constant M assumption with M[t] = k*(L[t]/Linf)^(-1.5)
plot.M <- TRUE    ###NEW### plot Charnov.M versus constant M value from input file?
do.biology.plots <- TRUE # 1 png file for each scenario; takes a minute or two, worth looking at once
do.Fmsy <- TRUE # replace input F with Fmsy based on maximizing yield
harvest.type <- 'biomass' # choose either 'biomass' (default) or 'numbers'
do.yield.vs.F.plot <- TRUE # plots yield over a grid of F values to help debug Fmsy search (adds run time)
do.egg.contour.plots <- TRUE
if(do.egg.contour.plots) { library(lattice) }

# read parameter file from working directory
parms.df <- read.csv('parameter_file_global_stocks_v3.csv')
# drop runs with 'exclude' not equal to zero
parms.df <- droplevels(subset(parms.df, exclude==0))
parms.df

# initialize list object for results, with one element per run; result will be a list of lists
results.list <- as.list(1:nrow(parms.df))
names(results.list) <- with(parms.df, paste(scenario,"_",species,sep=""))

# function to apply to each scenario
# units of individual fish length and weight are cm and kg, respectively
# Barneche et al. used weight in grams, so intercepts have been transformed to kg
model.fun <- function(i) {
# specify parameters
description <- with(parms.df[i,], paste(scenario,"_",species,sep=""))
R0 <- parms.df[i,"R0"]
A  <- parms.df[i,"A"]
M  <- parms.df[i,"M"]
h <- parms.df[i,"steepness"]
Linf.mean <- parms.df[i,"Linf.mean"]
Linf.min <- parms.df[i,"Linf.min"]
Linf.max <- parms.df[i,"Linf.max"]
k  <- parms.df[i,"k"]
t0 <- parms.df[i,"t0"]
Lmat.over.Linf <- parms.df[i,"Lmat.over.Linf"]
a  <- parms.df[i,"a"] # W=a*L^b
b  <- parms.df[i,"b"]
Lmin <- parms.df[i,"Lmin"]
Lmax <- parms.df[i,"Lmax"]
F  <- parms.df[i,"F"]
G <- parms.df[i,"G"]
CV <- parms.df[i,"CV"]
a.fec <- parms.df[i,"a.fec"] # F=(a.fec)*W^(b.fec)
b.fec <- parms.df[i,"b.fec"]

################################
### Preliminary Calculations ###
################################
# convert steepness to Goodyear's Compensation Ratio (see Martell et al. 2008, CJFAS, Appendix B)
CR <- 4*h / (1-h) # CR as function of h

# SPR.crash per Brooks et al. 2009
SPR.crash <- (1-h)/(4*h) # inverse of CR

### Growth and Fecundity ###
# given a CV for Linf (0.1 is default), compute range for Linf.min, Linf.max
qnorm(0.999, mean=1000, sd=0.1*1000)/1000 # if CV=0.1 and mean=1000, this is the multiplier for Linf.max

# vector of values for Linf
if (G==1) {
   Linf <- Linf.mean
} else {
   Linf <- seq(Linf.min, Linf.max, length=G); Linf
   # calculate proportion of recruits in each trajectory based on normal dist.
   half.bin <- diff(Linf)[1]/2 # used to center the length bins around the mean when G>1
}
# function to integrate normal dist over bin width (get proportions by bin)
prop.fun <- function(x, h) {
   out <- pnorm(x+h, mean=Linf.mean, sd=Linf.mean*CV) -
          pnorm(x-h, mean=Linf.mean, sd=Linf.mean*CV)
   return(out)
}
# if number of trajectories = 1, p.g = 1, else, use prop.fun
if(G==1)
{
   p.g <- 1
} else {
   p.g <- sapply(Linf, prop.fun, h=half.bin); p.g
   # normalize to sum to 1
   p.g <- p.g / sum(p.g)
}
p.g <- matrix(p.g,nrow=A,ncol=length(p.g),byrow=TRUE)
dimnames(p.g)[[1]] <- 1:A
dimnames(p.g)[[2]] <- Linf

# calculate growth trajectories
# length at age
age.vec <- 1:A
vB.fun <- function(age, LINF, K, T0) { LINF*(1-exp(-K*(age-T0))) }
L.at.age <- sapply(Linf, vB.fun, age=age.vec, K=k, T0=t0)
dimnames(L.at.age) <- list(age.vec, 1:G)
# weight at age
W.at.age <- a*L.at.age^b

# maturity by age and growth trajectory (Lmat = Linf * Lmat.over.Linf)
Lmat <- Linf * Lmat.over.Linf
Lmat <- matrix(Lmat,nrow=A,ncol=length(Lmat),byrow=TRUE)
p.mat <- (Lmat < L.at.age) + 0 # adding zero converts logical to numeric
p.mat
rm(Lmat)

# replace constant M (scalar from input file) with Charnov estimator (matrix; function of L, k, and Linf)
# if M remains constant (Charnov.M=FALSE), still need to specify M as a matrix
if( Charnov.M ) {
   M.constant <- M # save the value from the input file
   L.over.Linf <- sweep(L.at.age,2,Linf,"/") # this is the same for all growth trajectories when k is the same
   M <- k * L.over.Linf^(-1.5)
} else {
   M.constant <- M # save the value from the input file
   M <- matrix(M, nrow=nrow(L.at.age), ncol=ncol(L.at.age))
}

########################
### Population Model ###
########################
# Vulnerability by age and growth trajectory (V=1 for Lmin < L < Lmax, else 0)
V.at.age <- matrix(0, nrow=nrow(L.at.age), ncol=ncol(L.at.age))
dimnames(V.at.age) <- dimnames(L.at.age)
V.at.age[L.at.age > Lmin & L.at.age < Lmax] <- 1

# function that returns yield as a function of F
# user chooses to return catch (in numbers) or yield (in biomass)
Y.given.F <- function(F, type='biomass', M, V, A, G, R0, CR, W, a.fec, b.fec, p.mat, p.g) {
   Z.at.age  <- M + F*V
   S.at.age <- matrix(0, nrow=nrow(W), ncol=ncol(W))
   dimnames(S.at.age) <- dimnames(W)
   S.at.age[1,] <- 1
   S0.at.age <- S.at.age
   for (g in 1:G) {
      for (a in 2:A) { # skip age 1 (S=1)
         S0.at.age[a,g] <- S0.at.age[a-1,g] * exp(-M[a-1,g]) # only natural mortality
      }
   }
   for (g in 1:G) {
      for (a in 2:A) { # skip age 1 (S=1)
         S.at.age[a,g] <- S.at.age[a-1,g] * exp(-Z.at.age[a-1,g]) # natural and fishing mortality (Z)
      }
   }
   Fec.at.age <- a.fec*W.at.age^b.fec
   phi.unfished <- sum(S0.at.age * p.g * p.mat * Fec.at.age) # equation 7 (F=0)
   phi.fished   <- sum(S.at.age * p.g * p.mat * Fec.at.age)  # equation 7 (F>0)
   # NB: allow Req to go negative during search for Fmsy... otherwise optimize() can fail.
   #     Req definition in main script (below) has minimum of zero, as it should
   Req <- R0 * (CR - phi.unfished/phi.fished) / (CR - 1)
   N.at.age  <- Req * S.at.age * p.g
   C.at.age <- N.at.age * ( F*V.at.age / Z.at.age ) * ( 1 - exp(-Z.at.age) )
   catch <- sum(C.at.age)
   Y.at.age <- C.at.age * W.at.age
   yield <- sum(Y.at.age)
   if(type=='numbers') { return(catch) } else { return(yield) }
}

# use production model approximation to define range for Fmsy search and plotting
# since Charnov M estimator is a function of size, use the constant M estimate to define range
Fmsy.ballpark <- M.constant * (sqrt(4*h/(1-h))-1) # from Mangel et al.

if(do.Fmsy) { # if TRUE, solve for F that maximizes yield (in biomass or numbers), otherwise use input value of F
   # Find Fmsy; define search interval as F=0.0001 to F=3*Fmsy.ballpark
   F.opt <- optimize(Y.given.F, lower=0.0001, upper=3*Fmsy.ballpark, maximum=TRUE, type=harvest.type,
                     M=M, V=V.at.age, A=A, G=G, R0=R0, CR=CR, W=W.at.age, a.fec=a.fec, b.fec=b.fec, p.mat=p.mat, p.g=p.g)
   # replace input value of F with new value of Fmsy
   F <- F.opt[[1]]
   parms.df[i,"F"] <<- F # change to dataframe object outside of function requires global assignment ("<<-")
}

# calculate yield over a grid of 100 F values (slows things down a bit, but useful for comparisons and debugging
F.grid <- seq(0.0001, 3*Fmsy.ballpark, length=100)
Y.grid <- sapply(F.grid, FUN=Y.given.F, type=harvest.type, M=M, V=V.at.age, A=A, G=G, R0=R0,
                 CR=CR, W=W.at.age, a.fec=a.fec, b.fec=b.fec, p.mat=p.mat, p.g=p.g)
Y.grid[Y.grid < 0] <- 0 # Y.given.F function allows for negative yield (smooth function for optimize)

# Total Mortality by age and growth trajectory
Z.at.age  <- M + F*V.at.age

# Survival (S[a,g] is the same as 'l[a,g]' in Gwinn et al.)
S.at.age <- matrix(0, nrow=nrow(L.at.age), ncol=ncol(L.at.age)) # initialize
dimnames(S.at.age) <- dimnames(L.at.age)
S.at.age[1,] <- 1

# survival for unfished equilibrium state (F=0)
S0.at.age <- S.at.age # initialize object
for (g in 1:G) {
   for (a in 2:A) { # skip age 1 (S=1)
      S0.at.age[a,g] <- S0.at.age[a-1,g] * exp(-M[a-1,g])
   }
}

# survival for fished equilibrium state (F>0)
for (g in 1:G) {
   for (a in 2:A) { # S(1)=1 (see above)
      S.at.age[a,g] <- S.at.age[a-1,g] * exp(-Z.at.age[a-1,g])
   }
}

# fecundity at age
Fec.at.age <- a.fec*W.at.age^b.fec

# egg production at age
eggs.unfished.at.age <- S0.at.age * p.g * p.mat * Fec.at.age
eggs.fished.at.age <- S.at.age * p.g * p.mat * Fec.at.age

# total egg production = sum(survival * proportion in each growth trajectory * fecundity)
phi.unfished <- sum(S0.at.age * p.g * p.mat * Fec.at.age) # equation 7 (F=0)
phi.fished   <- sum(S.at.age * p.g * p.mat * Fec.at.age)  # equation 7 (F>0)

# Fished equilibrium number of age-1 recruits (equation 6 in Gwinn et al.)
Req <- max(0, R0 * (CR - phi.unfished/phi.fished) / (CR - 1))

# Equilibrium numbers at age by growth trajectory
N0.at.age <- R0 * S0.at.age * p.g # unfished equilibrium numbers at age

# exploited equilibrium numbers at age
N.at.age  <- Req * S.at.age * p.g

# C.at.age = catch in numbers
# Y.at.age = catch in biomass
# Baranov Catch Equation for catch in numbers at age
C.at.age <- N.at.age * ( F*V.at.age / Z.at.age ) * ( 1 - exp(-Z.at.age) )

# yield at age (catch at age times weight at age)
Y.at.age <- C.at.age * W.at.age

# Gwinn et al. define "J" as the proportion of fecundity produced by the older half of age classes
age.mid <- floor(A/2)+1 # could adjust for odd values of "A"; keep "A" even for best accuracy
# calculate J for unfished population (J0), for comparison to J
# sum across growth trajectories for ages "age.mid" to A, then divide by total fecundity
J0 <- sum((S0.at.age * p.g * p.mat * Fec.at.age)[age.mid:A,]) / phi.unfished
# calculate J for fished population (same as J in Gwinn et al.)
J <- sum((S.at.age * p.g * p.mat * Fec.at.age)[age.mid:A,]) / phi.fished

# create list of model outputs for run i
model.list <- list(parameters = parms.df[i,],
                   CR = CR,
				   Linf.vector = Linf,
				   prop.recruits.at.length = p.g,
				   length.at.age = L.at.age,
				   weight.at.age = W.at.age,
				   prop.mature.at.age = p.mat,
				   fecundity.at.age = Fec.at.age * p.mat, # mature fish only
				   vulnerability.at.age = V.at.age,
				   natural.mortality.at.age = M,
				   total.mortality.at.age = Z.at.age,
				   survival.unfished.at.age = S0.at.age,
				   survival.fished.at.age = S.at.age,
				   eggs.unfished.at.age = eggs.unfished.at.age,
				   eggs.fished.at.age = eggs.fished.at.age,
				   phi.unfished = phi.unfished,
				   phi.fished = phi.fished,
				   SPR = phi.fished / phi.unfished,
				   SPR.crash = SPR.crash,
				   fished.equilibrium.recruitment = Req,
				   numbers.unfished.at.age = N0.at.age,
				   N_0 = sum(N0.at.age),
				   numbers.fished.at.age = N.at.age,
				   N_F = sum(N.at.age),
				   B_0 = sum(N0.at.age * W.at.age),
				   B_F = sum(N.at.age * W.at.age),
				   catch.at.age = C.at.age,
				   catch = sum(C.at.age),
				   yield.at.age = Y.at.age,
				   yield = sum(Y.at.age),
				   harvest.type = harvest.type,
				   yield.vs.F = cbind.data.frame(F=F.grid, yield=Y.grid),
				   J0 = J0, # 'unfished' version of J
				   J = J # per footnote to Table 3 in Gwinn et al.
				   )
print(paste("Finished scenario",i))
return(model.list)
}

print(paste("Harvest units in", harvest.type))

# loop over runs (could do this in parallel if it gets too slow)
for (i in 1:nrow(parms.df)) { results.list[[i]] <- model.fun(i) }
print("Population model results are in 'results.list'")

#####################
### Model Outputs ###
#####################

# uncomment lines for more output
out.df <- cbind.data.frame(parms.df, # append parameter inputs; F replaced with Fmsy if do.Fmsy=TRUE
                           #compensation.ratio=sapply(results.list, with, CR),
						   SPR.crash=sapply(results.list, with, SPR.crash),
                           Req=sapply(results.list, with, fished.equilibrium.recruitment),
						   #N_0=sapply(results.list, with, N_0),
						   #B_0=sapply(results.list, with, B_0),
						   #N_F=sapply(results.list, with, N_F),
						   #B_F=sapply(results.list, with, B_F),
						   phi.unfished=sapply(results.list, with, phi.unfished),
						   phi.fished=sapply(results.list, with, phi.fished),
						   SPR=sapply(results.list, with, SPR),
						   #J0=sapply(results.list, with, J0),
						   #J=sapply(results.list, with, J), # definition depends on max. age; use with caution
                           catch=sapply(results.list, with, catch),
						   yield=sapply(results.list, with, yield))
						   #rel.catch=sapply(results.list, with, catch)/max(sapply(results.list, with, catch)),
						   #rel.yield=sapply(results.list, with, yield)/max(sapply(results.list, with, yield)))
out.df
write.csv(out.df, "results.csv", row.names=FALSE)

#########
# PLOTS #
#########
if(do.biology.plots) {
   bio.plot <- function(x) {
      description <- paste(x$parameters$scenario,"_",x$parameters$species,sep="")
	  png(filename = paste("biology_",description,".png",sep=""), width=7, height=7, units='in', res=600, antialias = "cleartype")
      par(mfrow=c(2,2))
      plot(x$Linf.vector, x$prop.recruits.at.length[1,], type='h', ylim=c(0,max(x$prop.recruits.at.length[1,])),
		   ylab="prop. of recruits",main=description, xlim=c(x$parameters$Linf.min, x$parameters$Linf.max),
		   bty='l', cex.axis=0.8, xlab="Linf", yaxs="i")
      legend('topleft', paste("G =",x$parameters$G), bty='n')
      plot(1:x$parameters$A, x$length.at.age[,x$parameters$G], xlab="age", ylab="length", type='l', lty=1, bty='l',
           ylim=c(0,1.05*x$parameters$Linf.max), cex.axis=0.9, xlim=c(0,x$parameters$A), xaxs="i", yaxs="i")
      if(x$parameters$G>1) { for (j in 1:(x$parameters$G-1)) { lines(1:x$parameters$A, x$length.at.age[,j]) } }
      abline(h=x$parameters$Lmin, lty=3); text(0.9*x$parameters$A, x$parameters$Lmin+0.05*x$parameters$Lmin, "Lmin", cex=0.7)
      abline(h=x$parameters$Lmax, lty=3); text(0.1*x$parameters$A, x$parameters$Lmax+0.05*x$parameters$Lmax, "Lmax", cex=0.7)
      plot(1:x$parameters$A, x$weight.at.age[,x$parameters$G], xlab="age", ylab="weight", xaxs="i", yaxs="i",
           type='l', lty=1, bty='l', cex.axis=0.9, xlim=c(0,x$parameters$A), ylim=c(0,max(x$weight.at.age[,x$parameters$G])))
      if(x$parameters$G>1) { for (j in 1:(x$parameters$G-1)) { lines(1:x$parameters$A, x$weight.at.age[,j]) } }
      plot(1:x$parameters$A, x$fecundity.at.age[,x$parameters$G], xlab="age", ylab="fecundity", type='l', xaxs="i", yaxs="i",
           lty=1, bty='l', cex.axis=0.9, xlim=c(0,x$parameters$A), ylim=c(0,max(x$fecundity.at.age[,x$parameters$G])))
      if(x$parameters$G>1) { for (j in 1:(x$parameters$G-1)) { lines(1:x$parameters$A, x$fecundity.at.age[,j]) } }
      dev.off()
   }
   lapply(results.list, bio.plot)
}

if (plot.M) {
   plot.M.fun <- function(x) {
      description <- paste(x$parameters$scenario,"_",x$parameters$species,sep="")
      png(filename = paste("M_",description,".png",sep=""), width=7, height=7, units='in', res=600, antialias = "cleartype")
      plot(1:x$parameters$A, rep(x$parameters$M, x$parameters$A), type='l', col='red', ylim=c(0,max(x$natural.mortality.at.age)),
           xlab="age", ylab="Natural Mortality Rate", main=description)
      lines(1:x$parameters$A, x$natural.mortality.at.age[,1])
      legend("topright",c("constant M","Charnov M"),lty=c(1,1), col=c(2,1))
      dev.off()
   }
   lapply(results.list, plot.M.fun)
}


# one plot of yield vs. F per species, comparing scenarios
if(do.yield.vs.F.plot) {
   spp <- as.character(unique(parms.df$species))
   H.type <- ifelse(harvest.type=="numbers", "Catch (numbers)", "Yield (biomass)")
   for(i in 1:length(spp)) { # one plot per species, comparing all scenarios
      tmp.parms <- subset(parms.df, species==spp[i])
	  tmp.list <- results.list[parms.df$species==spp[i]]
      max.yield <- max(unlist(lapply(tmp.list, function(x) max(x$yield.vs.F$yield))))
      max.F <- max(unlist(lapply(tmp.list, function(x) max(x$yield.vs.F$F))))
	  if(harvest.type=='numbers') {
         png(filename = paste("Catch_",spp[i],".png",sep=""), width=7, height=5, units='in', res=600, antialias = "cleartype")
	  } else {
         png(filename = paste("Yield_",spp[i],".png",sep=""), width=7, height=5, units='in', res=600, antialias = "cleartype")
	  }
	  plot(tmp.list[[1]][['yield.vs.F']]$F, tmp.list[[1]][['yield.vs.F']]$yield, main=spp[i],
	       type='l', xlab="F", ylab=H.type, bty='l', xlim=c(0,max.F), ylim=c(0,max.yield))
      abline(v=tmp.parms[1,"F"], lty=3)
	  if(length(tmp.list)>1) {
	     for(j in 2:length(tmp.list)) {
	        lines(tmp.list[[j]][['yield.vs.F']]$F, tmp.list[[j]][['yield.vs.F']]$yield,lty=j,col=j)
		    abline(v=tmp.parms[j,"F"],lty=3,col=j)
         }
      legend("topright",legend=tmp.parms[,"scenario"],lty=1:j, col=1:j)
      }
   dev.off()
   }
   rm(spp,H.type,tmp.parms,tmp.list,max.yield,max.F,i,j)
}

if(do.egg.contour.plots) {
   plot.eggs <- function(x) {
	  description <- paste(x$parameters$scenario,"_",x$parameters$species,sep="")
	  unfished.egg.production <- with(x, eggs.unfished.at.age/max(eggs.unfished.at.age))
	  png(filename = paste("Unfished_Eggs_",description,".png",sep=""), width=7, height=7, units='in', res=600, antialias = "cleartype")
	  # need to use print() with lattice plots (lattice objects not printed by default)
      print(contourplot(t(unfished.egg.production), xlab="growth trajectory", ylab="age", at=seq(0,1,by=0.1), labels=list(cex=0.7),
	                    aspect="xy", main=paste("unfished egg production, ",description, ", b=", x$parameters$b.fec, sep="")))
      dev.off()
   }
   lapply(results.list, plot.eggs)
}
