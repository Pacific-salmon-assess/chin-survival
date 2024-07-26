# Dipper example script
# From: https://www.afsc.noaa.gov/Publications/ProcRpt/PR2013-01.pdf

library(RMark)
library(tidyverse)
data(dipper)

# Process data specifying data, model (CJS is recaptures only), grouping
# variable, and initial time
dipper.proc <- process.data(dipper, model="CJS", groups="sex", begin.time=1981)
glimpse(dipper.proc)

# Generate design data list
dipper.ddl <- make.design.data(dipper.proc)

# Add dummy variable specifying a flood effect on phi (survival) in two years
dipper.ddl$Phi$flood <- ifelse(dipper.ddl$Phi$time %in% c(1982:1983), 1, 0)

dipper.analysis=function()
{
 # Create specifications for Phi and p
 Phi.1=list(formula=~time)
 Phi.2=list(formula=~-1+time,link="sin")
 Phi.3=list(formula=~sex+weight)
 Phi.4=list(formula=~flood)
 p.1=list(formula=~1)
 p.2=list(formula=~time)
 p.3=list(formula=~Time)
 p.4=list(formula=~sex)
 # Create a list of combinations of parameter specifications;
 # the argument "CJS" tells it to look for CJS parameters
 cml=create.model.list("CJS")
 # Call mark.wrapper; the arguments are cml and then like with
 # the mark function the processed data and ddl are specified,
 # but they must be named data= and ddl= because they are passed
 # through to the mark function; I've also set output=FALSE
 # so it doesn't show a summary for the 16 models but it will
 # show the model being fitted and any notations or errors.
 mark.wrapper(cml, data=dipper.proc, ddl=dipper.ddl, output=FALSE)
}

dipper.results = dipper.analysis()
