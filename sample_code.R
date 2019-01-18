#This code will take takes the win odds from a horse race and defines probability
#density functions for each horse using a Gaussian curve. This is done by finding
#parameter estimates for the mean of each independent distribution via maximum
#likelihood mehtods. After this is done generic values are chosen for the win odds
#and parameter estimates are found. Then to races are simulated and a the winning
#horse is recorded and plotted via a bar plot. The plot is animated as each race
#occurs. Finally the pdf's and bar plots are plotted together showing the means
#for each horse's pdf and below it the final proportion of times each horse won
#the race when simulated.


#This function determines the parameter estimates with inputs being the win odds
#and the error present in the maximum likelihood estimates. The output is a matrix
#containing the parameter estimates, probabilities and scale parameters.
normal.first<-function(winOdds = "", epsilon=""){
    n  <-  length(winOdds)

    #odds_conversion is a function that converts odds to probabilities that I wrote
    oD<-winOdds[order(winOdds)]
    probs  <-  (1/(1+oD))/sum(1/(1+oD))
    c  <-  1
    model.probs <- seq_len(n)

    #protects program from entering infinite loop by updating scale parameter if
    #program takes too long to find parameter estimates
    while ( sum(abs(probs-model.probs)) > epsilon ) {
        model.probs <- gen.vec <- seq_len(n)

        #reasonable initial parameter estimates
        new.shift <- probs ^ (-1)/sum(probs ^ (-1))
        jacob <- matrix(nrow = n,ncol = n)

        #starting time for first while loop
        s <- Sys.time()

        #updating scale parameter
        norm.scale <- c * rep(1, n)
        #continues estimated parameter search until error is less than epsilon
        while ( sum(abs(probs-model.probs)) > epsilon ) {
            shift <- new.shift
            #calculating parameter values for all n horse contestants
            for ( i in gen.vec ) {
                #function that calculates probabilities given current parameter states
                deriv.int <- function(shift){
                    center.int <- function(x){
                        p.center <- dnorm(x, shift[i], norm.scale[i])
                        for ( k in gen.vec[-i] ) {
                            p.center <- p.center * (1 - pnorm(x, shift[k], norm.scale[k]))
                        }
                        p.center
                    }
                    #Vectorize allows for an array to passed into a single input function
                    g.center <- Vectorize(center.int)
                    as.numeric(integrate(g.center, -Inf, Inf)[1])
                }
                #creating an environment that allows for numerical derivatives to be
                #found over the integral based maximum likelihood functions in "deriv.int"
                myenv <- new.env()
                assign("shift", as.numeric(shift), envir = myenv)
                prob.deriv <- numericDeriv(quote(deriv.int(shift)), "shift", myenv)
                jacob[, i] <- attr(prob.deriv, "gradient")
                model.probs[i] <- prob.deriv[1]
            }
            #updating the parameter estimates
            new.shift <- c(shift + (probs - model.probs) %*% ginv(jacob))

            #breaking the second while loop if taking too long
            if ( abs(s-Sys.time()) > 5 ) { break() }
        }
        #updates the scale parameter value
        c <- c + 1
    }
    #returning output matrix with final parameter estimates
    return(matrix(c(new.shift, model.probs, norm.scale), nrow = n, ncol = 3))	
}


#loading the pachage MASS to find generalized inverses of matrices
install.packages("MASS")
library(MASS)

#calling function "odds_conversion.R" to calculate odds probabilities
source("odds_conversion.R")

#calculating parameter values
winOdds <- seq(1, 10, 2)
(parameters <- normal.first(winOdds, 0.0001))
n <- dim(parameters)[1]

#Drawing s values from each pdf for each horse
s <- 10000
rDraws <- matrix(nrow = s, ncol = n)
for ( i in 1:n ) {
    rDraws[,i] <- rnorm(s, parameters[i, 1], parameters[i, 3])
}

#plotting and calculating simulated races to compare with empirical values
#need to install "animation"
#NOTE: depending on system compatabilities further installation may need
#to be done for the package "magick" and "curl".
#for ubuntu:
#sudo apt-get install libcurl4-openssl-dev
#sudo apt-get install libmagick++-dev
#install.packages(c("curl","magick"));library(curl);library(magick)
install.packages("animation")
library(animation)

#creating animated plot
par(bg = "white")
ani.record(reset = TRUE)
Hwins <- integer(n)
for ( i in 1:s ) {
    tempWin<-which(rDraws[i,] == min(rDraws[i,]))
    Hwins[tempWin] <- Hwins[tempWin] + 1

    #plotting only every 100th frame for visual asthetics.
    if ( i %% 100 == 0 ) {
        tempCol<-rep("grey",n)
        tempCol[tempWin]<-"light blue"
        barplot(Hwins/s, 10, col = tempCol, main = "Simulated Proportions of Wins",
                xlab = "Horses", ylab = "Proportion Won",
                names.arg = as.character(ceiling(Hwins / s * 100) / 100))
        ani.record()        
    }
}

#playing animated plot with 0.25 second intervals between frames
oopts = ani.options(interval = 0.25)
ani.replay()

#running this code will play the simulation online
saveHTML(ani.replay(), img.name = "Simulated_BarPlot")

#Plotting the model derived pdf's and bar plot
par(mfrow = c(2, 1))

#pdf curves
curve(dnorm(x, parameters[1, 1], parameters[1, 3]), parameters[1, 1] - 3, parameters[n, 1] + 3,
      main = "Gaussian Curves for each Horse", xlab = "", ylab = "", lwd = 2)
abline(v = parameters[1,1], lwd = 2)
for ( i in 2:n ) {
    curve(dnorm(x, parameters[i, 1], parameters[i, 3]), parameters[1, 1] - 3, parameters[n, 1] + 3,
         add=TRUE, lwd = 2)
    abline(v = parameters[i, 1], lwd = 2)
}

#bar plot
barplot(Hwins/s, n, main = "Simulated Proportions of Wins",
        xlab = "Horses", ylab = "Proportion Won",
        names.arg = as.character(ceiling(Hwins / s * 100) / 100))
