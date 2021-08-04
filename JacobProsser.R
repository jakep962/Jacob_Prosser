# Submitted by Jacob Prosser

# The following code is a subset from my MSc theisis. The main dataframe, 
# paramater values, and section had been modified. This may result in some
# funny or unrealistic population estimates. 

# The code is designed to 1) take environmental data such as temperature and 
# salinity and produce separate probability distributions that allows the 
# model to introduce random/stochastic aspects into a population model for 
# a four-staged marine parasite.

# I've selected to submit this code as it contains several moving parts (pull 
# datasets, produces probability distrubutions, calculates life-stages), uses
# several common and uncommon packages and shows how I organize my Rcode. 


#  Start: Stochastic temperature and salinity code
rm(list=ls())
setwd("/Users/jakeprosser/Desktop/Shared-with-Amy/thesis_work/second_chapter/data/Prosser_profiles_09_12/")

# Packages needed 
packs.needed <-c("tidyverse", "stringr", "dplyr", "base", 
                 "ggplot2", "patchwork", "bbmle", "chron", "cowplot")

lapply(packs.needed, require, character.only = TRUE)

# Loading pre-fitted functions data
#baydespoir_prof <- read.csv("../../../first_chapter/winter_2020/model_code/data_files/baydespoir/baydespoir_prof.csv")
load("../../../first_chapter/winter_2020/model_code/data_files/amy_regional/ParamsFuncs.RData")
load("../../../first_chapter/winter_2020/model_code/data_files/amy_regional/farmNL.RData")

# Importing the enviromental data
water.profile.2017.18 <- read.csv("../Prosser_profiles_09_12/a_MEDS_profile_prof_17_18.csv", header = T)
water.profile.2019.20 <- read.csv("../Prosser_profiles_09_12/a_MEDS_profile_prof_19_20.csv", header = T)

# Filtering to only the required columns and combining  
water.profile.2017.18.f <- water.profile.2017.18 %>% 
  select(DATA_TYPE:Q_POS, PSAL, TEMP)
water.profile.2019.20.f <- water.profile.2019.20 %>% 
  select(DATA_TYPE:Q_POS, PSAL, TEMP)

water.profile.df <- rbind(water.profile.2017.18.f, water.profile.2019.20.f) # The combined and filtered df

# The temperature data showed a seasonal pattern that was not observed 
# within the salinity column. Because of this, I produced a Juliean sequence 
# (first observation = 1 and following count up) to allow for a sin/cos function 
# to be fit. 
monthly.dates <- dates(paste(water.profile.df$OBS_MONTH, 
                             water.profile.df$OBS_DAY, 
                             water.profile.df$OBS_YEAR,
                             sep="/"))

# The following section formate the data to be used for the Tideverse/ggplot and Chron packages. 
water.profile.df$ggplot.date <- as.Date(paste(water.profile.df$OBS_YEAR, 
                                        water.profile.df$OBS_MONTH, 
                                        water.profile.df$OBS_DAY,
                                        sep="-")) 

water.profile.df$chron_date <-chron(dates=monthly.dates, 
                              origin. = c(month = 01,day = 01,
                                          year = 2009)) # setting into chron format 
# and setting the first observation

water.profile.df <- arrange(water.profile.df, chron_date) # ordering dates 

o <- as.Date("2009-01-01") # The origin or first observation date

water.profile.df$days_since_origin <- (as.numeric(as.Date(water.profile.df$chron_date) - o)) # Producing a seq from first to last observation

# Formating the data to only show the geological area wanted and removing empty 
# observations (only want data that has recorded temperature)
water.profile.df <- filter(water.profile.df, !is.na(TEMP)) 
profile.df.sub <- filter(water.profile.df, LATITUDE...N. > 46.1 & LATITUDE...N. < 48.04)
profile.df.sub <- filter(water.profile.df, LONGITUDE...E. < -54 & LONGITUDE...E. > -59.5)



# Showing the change in temperature and time and double checking that temperature 
# was subset and recorded correctly. 
ggplot(profile.df.sub, aes(x=days_since_origin, y=TEMP)) +
  geom_point() +
  geom_jitter() +
  xlab("Days Since Origin") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=60, hjust=1))

# Now that the Julian sequence has been recorded and subset to the specific geological region needed, 
# the seasonal function can be created. 

#Two methods were looked at, both used produced similar results but method two 
# was used do it's simplicity. Method 1 has been # out but can be quickly un-# 
# through the shortcut command-shift-c on Mac or by deleting each. 

# Method 1 to fit the data (found online)
# Time <- profile.df.sub$days_since_origin
# temperature <- profile.df.sub$TEMP
# 
# xc<-cos(2*pi*Time/366) # becomes b2
# xs<-sin(2*pi*Time/366) # becomes b1
# 
# fit.lm <- lm(temperature~xc+xs)
# 
# # access the fitted series (for plotting)
# fit <- fitted(fit.lm)
# 
# # find predictions for original time series
# pred <- predict(fit.lm, newdata=data.frame(Time=Time))
# plot(TEMP ~ days_since_origin, data= profile.df.sub)
# lines(fit, col="red")
# lines(Time, pred, col="blue")
# summary(fit.lm)

# Method 2
SSTlm <- lm(profile.df.sub$TEMP ~ sin((2*pi*profile.df.sub$days_since_origin)/365)
            + cos((2*pi*profile.df.sub$days_since_origin)/365)
)
summary(SSTlm)

# Model Parameters ---------------------------------------------------------------------------------------------------------
# Building the stochastic temperature and deterimistic temperature.
# Parameters are here to placed here to speed up code when jumping back into it. 
a <- mean(profile.df$TEMP, na.rm=T) # Mean temperature
b1 <- -6.959576 # Two parameters for the sin-cos function. 
b2 <- -3.567477 #Two parameters for the sin-cos function. 

# The theta parameter acts as the stochastic daily variation parameter for the seasonal temperature 
# function. The level of randomness/stochasticsty is controlled through the sd = level and mean= skewness. 
theta=rnorm(length(profile.df$days_since_origin), 0, 1)

# The temperature function with daily variation
Temp = function(t){
  Temp = a + b1*sin(2 * pi * t/365) + b2*cos(2 * pi * t/365) + theta
  return(Temp)
}

# The temperature function wihtout daily variaiton. 
det.Temp = function(t){
  Temp = a + b1*sin(2 * pi * t/365) + b2*cos(2 * pi * t/365)
  return(Temp)
}
# Salinity Asymetric Laplace Distrubution -----------------------------------------------
# Changing the bin and x-vectors
dx <- 0.5 # The scale that time will increase in
sub.edge <- seq(15, 40, dx) # Subset of values that salinity could fall in (observed within the data)
sub.mid <- seq(15+dx/2, 40-dx/2, dx) 
sub.x <- sub.mid
sub.bay <- subset(profile.df.sub,PSAL <= 35 & PSAL>= 15) # Subsetting to only take values that are greater then 15 (removed outliers)

# Plotting the subset distrubution 
sub.hist <- hist(sub.bay$PSAL, sub.edge) # Plotting the histogram when subsetted from 15-40
sub.freq <- sub.hist$counts
sub.freq <- sub.freq/sum(sub.freq) # Converting the count data into frequency
sum(sub.freq) # Should sum to 1

# Checking to make sure that they're the same length 
length(sub.freq)
length(sub.x)


# From frequency distrubution of the salinity we saw that it would best be represented by an asymetric Laplace distbrution 
# that contains 3 parameter values. The m, location parameter, pho and k. The location parameter can easily be calculated as
# it's the most frequently observed parameter. The other can be calculated through maxiumn likelihood. 

# Finding the m value 
sub.df <- as.data.frame(cbind(sub.x, sub.hist$counts))
sub.df[which.max(sub.df$V2),]  # The maxium value for m therefore is 32.2 
sub.m <- sub.df[which.max(sub.df$V2),][1,1] # Pulls only the m value

# The asymetric Laplace distrubution code within R. 
sub.asy.laplace <- function(sub.m, pho, k){
  d <- dx*(pho/(k+(1/k)))*exp(-((sub.x-sub.m)*pho*(sign(sub.x-sub.m))*k^(sign(sub.x-sub.m))))
} # Same distrubution but x was changed to sub.x. Sub.x is the x-vector for the subset data

# The following section caclulates the pho and k parameters through maxiumn likelhood.
sigma=1
sub.negLL = function(pho, k){
  -sum(dnorm(sub.freq, mean=sub.asy.laplace(sub.m, pho, k), sd=sigma, log=TRUE))
} # Only want to fit to pho and k so only include them within the code

sub.mle.asy.laplace.Ex <- mle2(sub.negLL, 
                               start=list(pho=1, k=3), ) # Similary, only want to fit to pho and k

# removing to only the parameter values needed. 
a.est <- unname(coef(sub.mle.asy.laplace.Ex))[-4]

# Fitted asymmetric Laplace distribution parameters.
print(a.est)

# Fitted parameter values 
sub.m <- 31.75
pho0 <- 2.3396181
k0 <-0.6618817

# The function that are being fitted 
sub.asy.laplace <- function(sub.m, pho, k){
  d <- (pho/(k+(1/k)))*exp(-((sub.x-sub.m)*pho*(sign(sub.x-sub.m))*k^(sign(sub.x-sub.m))))
}

# Life-history functions -------------------------------------------------------------------
epsilon <- function(T){
  b0 = 1.86
  b1 = -1.30
  T[T<0.000001]=0.000001
  A = exp(b0+b1*log(T/10))
  # Do not extrapolate beyond the longest egg string production time
  A[A>45.1]=45.1
  # The rate is the reciprocal of the time to produce the egg string
  epsilon = 1/A
  return(epsilon)
}

eta <- function(Temp){
  Temp[Temp <0.000001]=0.000001
  exp(5.6 - 0.43*log(Temp/10) - 0.78*log(Temp/10)^2)
}

gammaP <- function(T){
  b0 = 1.216
  b1 = -1.38
  T[T<0.000001]=0.000001
  A = exp(b0+b1*log(T/10))
  A[A>23.5]=23.5
  gammaP = 1/A
  return(gammaP)
}

V <- function(T,S){
  V=NULL
  Tvec = T
  Svec=S
  # The function needs to be capable of taking vector arguments for
  # use in other files
  for(i in seq(1,length(Tvec))){
    T = Tvec[i]
    S = Svec[i]
    T[T<0.000001]=0.000001
    # The fitted coefficients
    V = rbind(V,-2.2579 + 0.4557*log(T)+0.6373*log(S))
  }
  V[V<0]<-0
  V[V>1]<-1
  return(V)
}

# Life-history parameter -------------------------------------------------------
# Start of the stochastic temperature and salinity df

# Life stage initial conditions 
P0 <- 1
C0 <- 1
A0 <- 1 
I0 <- 1
P = P0
C = C0
A = A0
I = I0

# Parameter conditions 
iota <- iota1 # Attachment rate of lice
f <- 500000 # Number of fish within each pen
N <- 50 # The number of simulations 
TT <-365*3 # The total time
dtime = 0.5 # Time steps size
time <-  seq(0, TT, dtime) # The sequence of time
time0 <-time[1] # Time 0
temp.sal.simulation.df <- NULL # The df that the observations will go to

# The distribution is seperated into four loop sections. They are seperated on 
# different parameter value sections

for(sd in seq(0, 5, 1)){  # Daily variation temperature loop
  
  rmew <- 0 # Fitted mean value
  sd <- sd # standard deviation
  theta=rnorm(length(time), rmew, sd)
  
  Temp = function(t){
    Temp = a + b1*sin(2 * pi * t/365) + b2*cos(2 * pi * t/365) + sample(theta, 1, replace=TRUE)
    return(Temp)
  }
  
  Temp0 <- Temp(1) # First observed temp in the prediction.
  
  for(pho in c(0.2913, 0.4826, 0.6739, 1.0565, 
               1.4391, 1.8217, 2.339618, 2.5,
               3, 3.5)){ # Daily variation salinity loop
    
    # Initial value
    sub.m <- 31.75 # Should stay 32.25
    pho <- pho     # Should follow the vector above
    k <- 0.6618817 # Should stay 1.122308
    
    sub.sal.df <- as.data.frame(cbind(sub.x, sub.asy.laplace(sub.m, pho, k),
                                      cumsum(sub.asy.laplace(sub.m, pho, k))))
    names(sub.sal.df) <- c("Salinity Values", "PDF", "CDF")
    
    # Salinity Value
    salinity0 = sample(sub.sal.df$`Salinity Values`, 1, replace = TRUE, prob=sub.sal.df$PDF) # Need to define salinity0 within the loop
    Temp0 <- Temp0
    
    for(v in seq(1, N)){ # The number of simulations that will run 
      
      # Setting initial conditions again
      S <- salinity0[1]
      P = P0
      C = C0
      A = A0
      I = I0
      
      temp.sal.simulation.df <- rbind(temp.sal.simulation.df,c(sd, pho, v, time0, Temp0, S, P, C, A, I))
      
      for(t in time[2:length(time)]){ # Populates the matrix for each time step (each day a new estimate)
        
        # Random variable where the loop reruns v 10 times 
        S = sample(sub.sal.df$`Salinity Values`, 1, replace = TRUE, prob=sub.sal.df$PDF)
        
        # Stage
        P = max((P + (eta(Temp(t))*epsilon(Temp(t))*A*V(Temp(t), S) - mup(S)*P - gammaP(Temp(t))*P)*dtime), 0)
        I = max((I + (gammaP(Temp(t))*P -iota*f*I - mui(S)*I)*dtime), 0)
        C = max((C + (iota*f*I - gammaC(Temp(t))*C - muc(S)*C)*dtime), 0)
        A = max((A + (gammaC(Temp(t))*C - mua(S)*A)*dtime), 0)
        
        # Row binding everything into a single data frame
        temp.sal.simulation.df <- rbind(temp.sal.simulation.df,c(sd, pho, v, t, Temp(t), S, P, C, A, I))
        
      } # End of the time loop
    } # End of the v loop
    print(c("pho=", pho))
    
  } # End of the j salinity loop
  print(c("sd=", sd))
  print(Sys.time())   
} # End of the temperature loop



# Calculating the growth rate of the population ----------------------------------

Floq.sd.pho.data <- NULL # The dataframe that will hold the calculated Floquet exponents
sim.data <- ts.sd.total
for(p in seq(0, 10, 1)){
  test.sd <- subset(sim.data, sim.data$sd==p) # separates sim.data matrix into smaller systems
  
  for(s in c(0.2913, 0.4826, 0.6739, 1.0565, 
             1.4391, 1.8217, 2.339618, 2.5,
             3, 3.5, 4)){ 
    
    test.pho <- subset(test.sd, test.sd$pho ==s)
    
    for (z in seq(1, 50)){
      test.sub <- subset(test.pho, test.pho$Simulation.Number==z) # seperates and calculates the
      # Floquet theory for each of the n=100 simulations 
      
      SpectralBond = function(lambda, g){
        index <- NULL
        pastvals <- matrix(c(rep(A0, TT+1), rep(C0, TT+1), 
                             rep(I0, TT+1), rep(P0, TT+1)),
                           ncol=4, nrow=(TT+1)) # Produces a vector of A0 of length = TotalTime (TT)
        
        for(i in seq(1,50)){
          # Calculates 
          index[i] = max(c(max(test.sub$Adult.Females), 
                           max(test.sub$Chalimi.and.pre.adults), 
                           max(test.sub$Copepodids), 
                           max(test.sub$Nauplii)))# Finds the max abundance of adult females
          
          A0 <- tail(test.sub$Adult.Females, n=1)/index[i]
          C0 <- tail(test.sub$Chalimi.and.pre.adults, n=1)/index[i]
          I0 <- tail(test.sub$Copepodids, n=1)/index[i]
          P0 <- tail(test.sub$Nauplii, n=1)/index[i]
          
          pastvals=matrix(c(test.sub$Adult.Females/index[i],
                            test.sub$Chalimi.and.pre.adults/index[i],
                            test.sub$Copepodids/index[i], 
                            test.sub$Nauplii/index[i],
                            ncol = 4, nrow=(TT+1)
          ))
          
          # if(i>1 && (abs(index[i]-index[i-1])<1e-5)){
          #   break
          # }
          
          
          r=index[i]-1
          return(r)
        }
      }
      # Now calculating the Floquet
      FEst = log(SpectralBond(1,1)+1)/365
      Floq.sd.pho.data <-rbind(Floq.sd.pho.data, c(FEst, z, p, s)) # combines 
      # the calculated Floquet exponents, 
      # simulation number, and the scale parameter
      
    }
  }
}
Floq.sd.pho.data <- as.data.frame(Floq.sd.pho.data)
names(Floq.sd.pho.data) <- c("floquet", "sim", "sd", "pho")
head(Floq.sd.pho.data)

# After calculating the Floquet values for each sd and pho value, I know format
# the data frame so that the mean can be calculated.
Floq.sd.data.mean <- Floq.sd.pho.data %>% 
  split(list(.$sd, .$pho), lex.order = F) %>%
  map(~mean(.$floquet))

Floq.sd.data.mean <- t(as.data.frame(Floq.sd.data.mean))
Floq.sd.data.mean <- as.data.frame(format(
  round(Floq.sd.data.mean, 5),
  nsmall = 3))

Floq.sd.data.mean$V2 <- row.names(Floq.sd.data.mean)
Floq.sd.data.mean <- Floq.sd.data.mean %>% 
  separate(V2, c("sd", "zero", "pho"))

Floq.sd.data.mean$pho[is.na(Floq.sd.data.mean$pho)] <- 0
Floq.sd.data.mean$pho_combined <- paste(Floq.sd.data.mean$zero,
                                        Floq.sd.data.mean$pho, sep=".")

Floq.sd.data.mean$sd <- gsub("X","", 
                             as.character(Floq.sd.data.mean$sd))

Floq.sd.pho.data.og <- Floq.sd.pho.data
Floq.sd.pho.data <- subset(Floq.sd.pho.data, Floq.sd.pho.data$pho <3.9)

# Formating the data 
temp.sal.simulation.df

ts.sd.total.sub <- filter(temp.sal.simulation.df, 
                          temp.sal.simulation.df$`Simulation Number` ==1 &
                            temp.sal.simulation.df$pho < 3.9)

ts.sd.total.abund <- filter(ts.sd.total.sub, 
                            ts.sd.total.sub$sd == c(0, 1, 5, 10))

ts.sd.total.abund <- filter(ts.sd.total.abund, 
                            ts.sd.total.abund$Time > 2919.5)
# Producing a name column 
ts.sd.total.abund$pho[ts.sd.total.abund$pho == 2.3396181] <- 2.339618

ts.sd.total.abund$test <- ts.sd.total.abund$pho
ts.sd.total.abund$test[ts.sd.total.abund$pho == c(0.2913, 0.4826, 0.6739,
                                                  1.0565, 1.4391, 1.8217, 
                                                  2.339618, 2.5, 3, 3.5)] <- c(expression(lambda~ "=0.2913"),
                                                                               expression(lambda~"=0.4826)"),
                                                                               expression(lambda~"=0.6739"),
                                                                               expression(lambda~"=1.0565"),
                                                                               expression(lambda~"=1.4391"),
                                                                               expression(lambda~"=1.8217"),
                                                                               expression(lambda~"=2.3396"),
                                                                               expression(lambda~"=2.5"),
                                                                               expression(lambda~"=3.0"),
                                                                               expression(lambda~"=3.5"))


ts.sd.total.abund$test[ts.sd.total.abund$pho == 0.2913] <- expression(lambda~ "=0.2913")

# Looking at the Floquet exponents 
ggplot(Floq.sd.pho.data, aes(x=factor(sd), y=floquet)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(x="Theta Standard Deviation Value", y="") +
  facet_wrap(~ pho, ncol = 5) +
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(x=expression(paste("Standard Deviation Value of Theta "(theta))), 
       y=expression("Floquet Exponents  "(phi))) +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_x_discrete(limits = levels(Floq.sd.pho.data$sd)) + 
  annotate("text", x=1:11, y = 0.131, label = floq.mean, angle = 45) 


ts.sd.total.abund$pho <- as.character(ts.sd.total.abund$pho)

ggplot(ts.sd.total.abund, aes(x=factor(Time), y=log(`Adult Females`), 
                              group=factor(sd))) + 
  geom_line(aes(colour=factor(sd))) +
  facet_wrap(~ pho, ncol = 5) +
  scale_color_manual(name = expression(paste("Standard \nDeviation \nValue of Theta "(theta))), 
                     breaks = c("10", "5", "1", "0"),
                     values=rev(c("#E69F00", "#56B4E9", "#009E73", "#D55E00"))) + 
  scale_x_discrete(breaks=c("2920", "3285", "3650"),
                   labels=c("2920", "3285", "")) +
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(x="Time (days)", 
       y= "Log Abundance of Adult Females")

ggplot(ts.sd.total.abund, aes(x=Time, y= log(Adult.Females))) + 
  geom_line(colour=sd) + 
  labs(x="Theta Standard Deviation Value", y="") +
  facet_wrap(~ pho) +
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(x=expression(paste("Standard Deviation Value of Theta "(theta))), 
       y=expression("Floquet Exponents  "(phi))) +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  
  
  # Abundance Plot 
  # Abundance, temperature, and ratio plot
  adult.abundance.sub <- filter(sd.total, sd.total$sd==c(0, 1, 5, 10))

ggplot(sd.total, aes(x= Time,  y=log(`Adult Females`), 
                     colour=factor(sd))) +
  geom_line(data = subset(sd.total, sd == 10), size = 1) +
  geom_line(data = subset(sd.total, sd == 5), size = 1) +
  geom_line(data = subset(sd.total, sd == 1), size = 1) +
  geom_line(data = subset(sd.total, sd == 0), size = 1) +
  labs(x="Time (days)", y="Log Abundance of Adult Females") +
  theme_classic() + 
  scale_color_manual(name = "Theta \nValue", 
                     breaks = c("10", "5", "1", "0"),
                     values=c("#000000", "#FF9900", "#006600", "#0000FF")) 


ggplot(sd.total, aes(x= Time,  y=Temperature, 
                     colour=factor(sd))) +
  geom_line(data = subset(sd.total, sd == 10), size = 1) +
  geom_line(data = subset(sd.total, sd == 5), size = 1) +
  geom_line(data = subset(sd.total, sd == 1), size = 1) +
  geom_line(data = subset(sd.total, sd == 0), size = 1) +
  labs(x="Time (days)", y="Temperature (celsius)") +
  theme_classic() + 
  scale_color_manual(name = "Theta \nValue", 
                     breaks = c("10", "5", "1", "0"),
                     values=c("#000000", "#FF9900", "#006600", "#0000FF")) 

ggplot(sd.total, aes(x=Time,  y= Salinity)) +
  geom_line(data = subset(sd.total, sd == 10), size = 1) +
  labs(x="Time (days)", y="Salinity (psu)") +
  theme_classic() +  
  scale_color_manual(breaks = c("1.8217", "0.6739"),
                     values=c("#FF0000", "#00FFFF")) + 
  theme(legend.position="NA")  


plot_grid(ratio.plot, ratio.plot.log, labels = c("", "log"))

plot_grid(adult.salinity, adult.abundance, ratio.plot, labels = c("(a)", "(b)", "(c)"), 
  