#The detections and nondetections of eDNA must be saved in a data frame that includes a column containing the survey locations, a column containing integer-valued sample numbers (i.e., 1, 2, \ldots) of each survey location, and one or more columns (one for each subsample or replicate of the sample) containing  binary-valued indicators of whether eDNA was detected (1) or not detected (0).  Thus, each row of this data frame corresponds to the detections and nondetections of eDNA in a single sample.  The function \code{occData()} uses this data frame to compute \code{detectionMats}, the list of data matrices required by \code{occModel()}.

#If multi-sca#le occupancy models are to be fitted using covariates of eDNA occurrence or detection, these covariates must be saved as columns of a data frame.  This data frame also must include a column of survey locations.  If the values of covariates do not differ among samples of a location, the data frame will contain a single row for each survey location.  This kind of data frame is passed using the \code{siteData} argument of \code{occModel()}.    If the values of one or more covariates differ among samples of a survey location, the data frame must include an additional column for the integer-valued sample numbers of each survey location.  This kind of data frame is passed using the \code{siteAndSampleData} argument of \code{occModel()}.  Only one data frame of covariates can be used in a single model, that is, either the \code{siteData} or the \code{siteAndSampleData} argument is used, not both arguments.





library(eDNAoccupancy)
data(fungusDetectionData)
data(fungusSurveyData)

fungusDetections = occData(fungusDetectionData, siteColName = 'site',
                            sampleColName = 'sample')
 #
 ## number of detections per sample
 head(fungusDetections$y)
 ## number of PCR replicates per sample
 head(fungusDetections$K)

#We fit a multi-scale occupancy model without covariates and print a summary of the parameter estimates using the following code.

 set.seed(0157)
 fit = occModel(detectionMats=fungusDetections, niter=11000,
                niterInterval=5000)
 posteriorSummary(fit, burnin=1000, mcError=TRUE)


## Center and scale numeric-valued covariate measurements
fungusSurveyData.sc = scaleData(fungusSurveyData)

set.seed(0157)
fit = occModel(formulaSite          = ~ 1,
               formulaSiteAndSample = ~ frogs,
               formulaReplicate     = ~ 1,
               detectionMats        = fungusDetections,
               siteData             = fungusSurveyData.sc,
               niter                = 6000,
               niterInterval        = 2000,
               siteColName = 'site'
               )
posteriorSummary(fit, burnin=1000, mcError=TRUE)


#If we want to assess whether the Markov chain used to compute these estimates appears to have converged, trace plots of the parameters may be created as follows (Fig.~\ref{fig:TracePlotFungusAnalysis}).
plotTrace(fit, c('beta.(Intercept)', 'alpha.(Intercept)', 'alpha.frogs',
            'delta.(Intercept)'),  burnin=1000)

#Autocorrelation plots of the parameters are created similarly
 plotACF(fit, c('beta.(Intercept)', 'alpha.(Intercept)', 'alpha.frogs',
             'delta.(Intercept)'),  burnin=1000)

#After inspection of these plots, suppose we decide that the MCMC algorithm needs to be run longer, either to eliminate the transient portion of the Markov chain or to reduce Monte Carlo errors in the parameter estimates.  We can resume the MCMC algorithm for the currently fitted model as follows.
 fit = updateOccModel(fit, niter=5000, niterInterval=2000)
 posteriorSummary(fit, burnin=1000,  mcError=TRUE)

 #These estimates of the parameters are computed using the updated Markov chain containing \rinline{fit$niterations} iterations.  The Monte Carlo errors in these parameter estimates are slightly lower, but the estimates are otherwise similar to those computed with only \rinline{fit$niterations-5000} iterations.

#In addition to estimating posterior summaries of the model's formal parameters, we also may be interested in estimating posterior summaries of derived parameters.  For example, in the second model fitted to the \emph{Bd} data, the probability of eDNA occurrence in ponds was assumed to be constant ($\psi$), the conditional probability of eDNA occurrence in samples was assumed to be a function of the frog density index \code{frogs}, and the conditional probability of eDNA detection was assumed to be constant ($p$).  The posterior medians of these derived parameters are estimated as follows.

psi = posteriorSummaryOfSiteOccupancy(fit, burnin=1000)
theta = posteriorSummaryOfSampleOccupancy(fit, burnin=1000)
p = posteriorSummaryOfDetection(fit, burnin=1000)

 ## output estimates of posterior medians
 cbind(psi=psi$median, theta=theta$median[,1], p=p$median[,1])

 frogs = fungusSurveyData[, 'frogs']
 plot(frogs, theta$median[,1], ylim=c(0,1), xlim=c(0,0.8), cex=2)
 segments(frogs, theta$lower[,1], frogs, theta$upper[,1], lwd=2)
 
#One way to assess the relative importance of such estimated relationships is to compare competing models using model-selection criteria.  For example, we compute the PPLC and WAIC criteria for the previously fitted model as follows.
 posteriorPredictiveLoss(fit, burnin=1000)
 WAIC(fit, burnin=1000)
