# fit crawl and draw random imputations from a fitted track
crawl_fit_and_imp_examp <- function(track_num){

library(crawl)

load("~/Desktop/NMML_Phase_2/R_scripts/NAP_Hex_Grid_and_Covariate_Operations/an355.RData")

initial = list(
  a=c(coordinates(td)[1,1],0,0,coordinates(td)[1,2],0,0),
  P=diag(c(10000^2,5400^2,5400^2,10000^2,5400^2,5400^2))
)

fixPar = c(log(250), log(500), log(1500), rep(NA,6))
displayPar( mov.model=~1, err.model=list(x=~loc_class-1), drift = TRUE,
            data=td,fixPar=fixPar)
constr=list(
  lower=c(rep(log(1500),2), rep(-Inf,4)),
  upper=rep(Inf,6)
)

fit1 <- crwMLE(
  mov.model=~1, err.model=list(x=~loc_class-1), drift=TRUE,
  data=td, Time.name="Time", 
  #theta=theta,
  initial.state=initial, fixPar=fixPar, constr=constr,
  control=list(trace=1, REPORT=1)
)

## Add regular time grid for simulation
low_time = min(grid_times[grid_times>=min(fit1$data$Time)])
hi_time = max(grid_times[grid_times<=max(fit1$data$Time)])
predTime <- seq(low_time, hi_time, 600)
mu <- crwPredict(object.crwFit=fit1, predTime, flat=TRUE)

# Create track simulator
simObj <- crwSimulator(fit1, predTime, parIS=0)

##
print("Generating imputations...")
reps = 30
samp = vector("list",reps)
for(i in 1:reps){
  samp[[i]] <- crwPostIS(simObj, fullPost = FALSE)$alpha.sim[simObj$locType=="p",c("mu.x","mu.y")]
  # lines(samp[[i]])
}
output <- list(imputations=samp,time=predTime)
return(output)
}
