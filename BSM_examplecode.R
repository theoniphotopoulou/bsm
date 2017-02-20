###################################################################################
# NOTES ON USING FUNCTIONS: get.BSMdive(), find.Xzone(), xzone.wrapper()
# These two functions work on dive profiles abstracted using the broken-stick model, or detailed time-depth dive data
#
# Theoni Photopoulou 20140109
###################################################################################

# Set your workspace to where you save these files and load the data
# setwd("/Users/user/Desktop")
setwd("/Users/theoniphotopoulou/Documents/PUBLICATIONS/02BSM/MEE/S1")
source("BSM_functions.R")
load('BSM_sampledata.RData')

################ get.BSMdive
# This function does two things, depending on the kind of dive data you give it.
# a) If you give it a vector of depths sampled at regular intervals, it carries out abstraction using the broken-stick model
# b) If you give it a dive that has already been abstracted on-board an SRDL, it works out the order in which the points were added to the abstracted profile 

## ARGUMENTS TO FUNCTION - use args(get.BSMdive) to see the defaults

# BS.dive: is a logical argument indicating whether the dive data have already been abstracted by the broken-stick model
# detailed.depth: is a vector of depths and this argument should be provided when you have high resolution data
# res: takes an integer and is the resolution of the high resolution data in seconds, its default is 4 sec
# divetable: is a dataframe and should be provided when you have already abstracted dive data. this should have column titles the same as the SMRU divetables 
# breakpoints: is an integer and is the number of breakpoints you want in your abstracted dive. the default is 4
# dn: is an integer indicates the dive number, essentially the row index in the divetable for the dive you want to abstract
# plot: is a logical argument and indicates whether you want your dive to be plotted
# dev.new: is a logical argument and indicates whether you want your dive to be plotted in a new window
# plot.truth: is a logical argument and indicates whether you want the detailed dive to be plotted. this only applies when you have detailed dive data.
# plot.bsm: is a logical argument and indicates whether you want the broken-stick abstracted dive to be plotted
# cex.axis: is a real number and indicates the size of the numbers in the axes
# cex.lab: is a real number and indicates the size of the axis labels
# draw.numbers: is a logical argument and indicates whether you want numbers indicating the order of the points selected by the broken-stick be plotted
# draw.lines: is a logical argument and indicates whether you want lines to be plotted, one horizontal line at max.depth and one vertical one at the time of max.depth

################
# CHECK HOW IT WORKS WITH DETAILED DIVES
################

dtest1 <- get.BSMdive(BS.dive=FALSE, detailed.depth=detailed.dive.data[[1]]$DEPTH, res=4, divetable, breakpoints=4, dev.new=T, plot.truth=TRUE, plot.bsm=TRUE, cex.axis=1.3, cex.lab=1.3, draw.numbers=TRUE, draw.lines=TRUE)
dtest2 <- get.BSMdive(BS.dive=FALSE, detailed.depth=detailed.dive.data[[2]]$DEPTH, res=4, divetable, breakpoints=4, dev.new=T, plot.truth=TRUE, plot.bsm=TRUE, cex.axis=1.3, cex.lab=1.3, draw.numbers=TRUE, draw.lines=TRUE)

################
# CHECK WITH ABSTRACTED DIVES
################

dn <- 2
bsmtest <- get.BSMdive(BS.dive=TRUE, res=4, divetable=abstracted.dive.data, breakpoints=4, dn=dn, plot=TRUE, dev.new=FALSE, plot.truth=FALSE, plot.bsm=TRUE, cex.axis=1.3, cex.lab=1.3, draw.numbers=TRUE, draw.lines=TRUE)

################ find.Xzone
# This function calculates the dive zone for a dive that has been abstracted using the broken-stick algorithm 

## ARGUMENTS TO FUNCTION - use args(find.Xzone) to see the defaults

# output: the object holding the output from get.BSMdive 
# BS.dive: is a logical argument indicating whether the dive data have already been abstracted by the broken-stick model
# plot: is a logical argument and indicates whether you want your dive to be plotted
# plot.legend: is a logical argument and indicates whether you want ta legend to be plotted
# cex.numbers: is a real number and indicates the size of the numbers indicating the BSM points
# cex.axis: is a real number and indicates the size of the numbers in the axes
# cex.lab: is a real number and indicates the size of the axis labels 
# lwd: is a real number and indicates the width of the line of the dive zone  
# dev.new: is a logical argument and indicates whether you want your dive to be plotted in a new window 
# plot.xzone: is a logical argument and indicates whether you want the dive zone to be plotted
# plot.truth: is a logical argument and indicates whether you want the detailed dive to be plotted. this only applies when you have detailed dive data. 
# plot.bsm: is a logical argument and indicates whether you want the broken-stick abstracted dive to be plotted
# plot.title: is a logical argument and indicates whether you want a title to be plotted
# draw.numbers: is a logical argument and indicates whether you want numbers indicating the order of the points selected by the broken-stick be plotted

################
# CHECK HOW IT WORKS WITH DETAILED DIVES
################

d.xzone <- find.Xzone(output=dtest1, BS.dive=F, plot.truth=TRUE)

################
# CHECK WITH ABSTRACTED DIVES
################

bsm.xzone <- find.Xzone(output=bsmtest, BS.dive=T, plot.truth=FALSE)

################ xzone.wrapper
# This wrapper function lets you run the abstraction and/or the dive zone function(s) on whole datasets of dives, in "SMRU divetable" form

## ARGUMENTS TO FUNCTION - use args(xzone.wrapper) to see the defaults

# BS.dive: the object holding the output from get.BSMdive 
# divetable: is a logical argument indicating whether the dive data have already been abstracted by the broken-stick model
# start.dive: is a logical argument and indicates whether you want your dive to be plotted
# end.dive: is a logical argument and indicates whether you want ta legend to be plotted
# breakpoints: is an integer and is the number of breakpoints you want in your abstracted dive. the default is 4
# plot.bsm: is a logical argument and indicates whether you want the broken-stick abstracted dive to be plotted
# plot.xzone: is a logical argument and indicates whether you want the dive zone to be plotted
# cex.axis: is a real number and indicates the size of the numbers in the axes
# cex.lab: is a real number and indicates the size of the axis labels 

################
# CHECK HOW IT WORKS WITH DETAILED DIVES
################

# d.w <- xzone.wrapper(BS.dive=FALSE, divetable=divetable, start.dive=1, end.dive=3, breakpoints=4, plot.bsm=F, plot.xzone=T, cex.axis=1.3, cex.lab=1.3)

################
# CHECK WITH ABSTRACTED DIVES
################
divetable <- abstracted.dive.data
bsm.w <- xzone.wrapper(BS.dive=TRUE, divetable=divetable, start.dive=1, end.dive=nrow(divetable), breakpoints=4, plot.bsm=F, plot.xzone=T, cex.axis=1.3, cex.lab=1.3)
str(bsm.w)


