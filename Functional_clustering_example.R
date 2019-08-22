# MA902 FINAL REPORT - R EXAMPLE - REFERENECE NUMBER: 1700264 
# Package used: "funFEM" package
# Data used: "Canadian Temperature" data (Ramsay & Silverman)
library(fda)
library(funFEM)
weather_data <- CanadianWeather
str(weather_data)
print((weather_data$dailyAv[,2,"Temperature.C"]), col = red)
plot(weather_data$dailyAv[ ,35,"Temperature.C"], col = "red")

## CREATE THE B-SPLINE-BASIS 
b_spline_basis <- create.bspline.basis(c(0, 365), nbasis = 21, norder = 4)
# c(0,365) = Interval of the functional data
# nbasis = 21: number of basis functions
# norder = 4: order of polynomilas + 1 (4 corresponds to cubic)

## USE THE B-SPLINE-BASIS TO APPROXIMATE THE DATA
func_data <- smooth.basis(day.5, weather_data$dailyAv[ , ,"Temperature.C"], b_spline_basis,
                      fdnames=list("Days", "Station", "Degrees in C"))$fd
# 1st argument: 
# 2nd argument: discrete data
# 3rd argument: the basis that will be used
# $fd: functional data object containing a smooth of the data

## CLUSTER THE DATA WITH THE funFEM PACKAGE
clusters <- funFEM(func_data, K = 4)
# func_data: functional data to be clustered
# K: number of clusters

# VISUALIZE THE CLUSTERED FUNCTIONS
plot(func_data, col = clusters$cls, lwd = 2, lty = 1, main = "Clustered functions of size 4")
