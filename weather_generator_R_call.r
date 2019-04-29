
###
## Function to extractmet data from global field ECMWF data
###

daily_mean <-function(var, interval, missing_allowed) {
      # work out how many intervals fit
      # i.e. number of days possible
      nos_days = ceiling(length(var)/interval)
      output = array(NaN, dim=c(nos_days))
      b = 0
      for (i in seq(1, nos_days)) {
	         if (length(which(is.na(var[b:(b+interval)]))) < missing_allowed) {
	             output[i] = mean(var[b:(b+interval)], na.rm=T)
	         } else {
               output[i] = NaN
           }
	         b = b+interval
      }
      return(output)
}

daily_min <-function(var, interval, missing_allowed) {
      # work out how many intervals fit
      # i.e. number of days possible
      nos_days = ceiling(length(var)/interval)
      output = array(NaN, dim=c(nos_days))
      b = 0
      for (i in seq(1, nos_days)) {
	         if (length(which(is.na(var[b:(b+interval)]))) < missing_allowed) {
	             output[i] = min(var[b:(b+interval)], na.rm=T)
	         } else {
               output[i] = NaN
           }
	         b = b+interval
   }
   return(output)
}

daily_max <-function(var, interval, missing_allowed) {
      # work out how many intervals fit
      # i.e. number of days possible
      nos_days = ceiling(length(var)/interval)
      output = array(NaN, dim=c(nos_days))
      b = 0
      for (i in seq(1, nos_days)) {
	         if (length(which(is.na(var[b:(b+interval)]))) < missing_allowed) {
	             output[i] = max(var[b:(b+interval)], na.rm=T)
	         } else {
               output[i] = NaN
           }
	         b = b+interval
      }
      return(output)
}

calc_photoperiod_sec<-function(lat,days){

	     # function calculates the day length in hours based on day of year and latitude (degrees).
	     # the output is daylength converted to seconds

       declin    = - asin ( sin ( 23.45 * ( pi / 180 ) ) * cos ( 2. * pi * ( days + 10. ) / 365. ) )
       sinld     = sin ( lat*(pi/180.) ) * sin ( declin )
       cosld     = cos ( lat*(pi/180.) ) * cos ( declin )
       aob       = sinld / cosld
       daylength = 12.0 * ( 1. + 2. * asin ( aob ) / pi )
	     # convert hours to seconds
	     daylength=daylength*3600
	     # now return
	     return(daylength)

} # end function calc_photoperiod_sec

sp_humidity_to_vpd<-function(sp_moist,atmos_press,air_temperature) {

       # Converts specific humidity (kg/kg) into VPD (Pa)
       # Determine vapour pressure (Pa); based on specific humidity, air
       # pressure (Pa input)
       # (hPa) and ratio of gas constants (p95 McIlveen 1986, Basic Meteorology -
       # a physical outline)

       # calculate vapour pressure of the air
       vpair = ((sp_moist*(atmos_press*1.0e-2))/0.62197)*1.0e2
       # Saturation vapour pressure (Pa) calculation from Jones p110; uses
       # absolute air temperature (oC)
       vpsat = (0.061375*exp((17.502*air_temperature)/(240.97+air_temperature)))*1.0e4
       # Difference between the vapour pressure of saturation and air, i.e. the
       # VPD (Pa)
       vpd_pa = vpsat-vpair
       # return to user
       return(vpd_pa)

} # end function sp_humidity_to_vpd

vpd_to_rh<-function(vpd_in,air_temperature) {

      # Converts VPD (Pa) to rel humidity (frac)
      # Determine vapour pressure (Pa); based on specific humidity, air
      # pressure (Pa input)
      # (hPa) and ratio of gas constants (p95 McIlveen 1986, Basic Meteorology -
      # a physical outline)

      # Saturation vapour pressure (Pa) calculation from Jones p110; uses
      # absolute air temperature (oC)
      vpsat = (0.061375*exp((17.502*air_temperature)/(240.97+air_temperature)))*1.0e4
      # Difference between the vapour pressure of saturation and air, i.e. the
      # VPD (Pa)
      rh = vpd_in/vpsat
      # return to user
      return(rh)

} # end function vpd_to_rh

rh_to_vpd<-function(rh_in,air_temperature) {

      # Converts rel humidity (frac) into VPD (Pa)
      # Determine vapour pressure (Pa); based on specific humidity, air
      # pressure (Pa input)
      # (hPa) and ratio of gas constants (p95 McIlveen 1986, Basic Meteorology -
      # a physical outline)

      # Saturation vapour pressure (Pa) calculation from Jones p110; uses
      # absolute air temperature (oC)
      vpsat = (0.061375*exp((17.502*air_temperature)/(240.97+air_temperature)))*1.0e4
      # Difference between the vapour pressure of saturation and air, i.e. the
      # VPD (Pa)
      vpair = vpsat*rh_in
      vpd_pa = vpsat-vpair
      # return to user
      return(vpd_pa)

} # end function rh_to_vpd

###
## Script to test call the weather generator
## for the purpose downscaling daily meteorological drivers to hourly
###

setwd("/home/lsmallma/WORK/R/Scripts/weather_generator")

library(ncdf4)
duke = nc_open("/home/lsmallma/WORK/GREENHOUSE/observations/site_level/US-Dk3/AMF_USDk3_2005_L2_WG_V003.nc")
names(duke$var)

# read raw
time = ncvar_get(duke, "HRMIN")
time = ncvar_get(duke, "DOY")
airt = ncvar_get(duke, "TA")
#lat = ncvar_get(duke, "Latitude")
swrad = ncvar_get(duke, "Rg") # W.m-2
windsp = ncvar_get(duke, "WS") # m.s-1
precip = ncvar_get(duke, "PREC") # mm.t-1
vpd = ncvar_get(duke, "VPD") # kPa
rh = ncvar_get(duke, "RH") # hPa
#pressure = ncvar_get(duke, "Pressure") # kPa

par(mfrow=c(3,4))
plot(as.vector(airt))
plot(as.vector(rh))
plot(as.vector(vpd))
plot(as.vector(swrad))
#plot(as.vector(pressure))
plot(as.vector(precip))
plot(as.vector(windsp))

# transform some drivers
#rh=vpd_to_rh(vpd,airt)
rh = rh*1e-2

# hack
swrad[which(is.na(swrad))] = 0
swrad[which(swrad == -9999)] = 0

airt_hr = array(NA, dim=(length(airt)*0.5))
precip_hr = array(NA, dim=(length(airt)*0.5))
swrad_hr = array(NA, dim=(length(airt)*0.5))
rh_hr = array(NA, dim=(length(airt)*0.5))
vpd_hr = array(NA, dim=(length(airt)*0.5))
windsp_hr = array(NA, dim=(length(airt)*0.5))
time_hr = array(NA, dim=(length(airt)*0.5))
b = 1
for (i in seq(2,length(airt),2)){
    time_hr[b] = time[i]
    airt_hr[b] = mean(airt[(i-1):i],na.rm=T)
    precip_hr[b] = mean(precip[(i-1):i]/(30*60),na.rm=T)
    swrad_hr[b] = mean(swrad[(i-1):i],na.rm=T)
    rh_hr[b] = mean(rh[(i-1):i],na.rm=T)
    vpd_hr[b] = mean(vpd[(i-1):i],na.rm=T)
    windsp_hr[b] = mean(windsp[(i-1):i],na.rm=T)
    b = b+1
}

# generate max, min, mean values
time = daily_mean(time_hr, 24, 24)
sat_avg_in = daily_mean(airt_hr, 24, 24)
rh_avg_in = daily_mean(rh_hr, 24, 24)
sat_max_in = daily_max(airt_hr, 24, 24)
rh_max_in = daily_max(rh_hr, 24, 24)
sat_min_in = daily_min(airt_hr, 24, 24)
rh_min_in = daily_min(rh_hr, 24, 24)
swrad_in = daily_mean(swrad_hr, 24, 24)
ppt_in = daily_mean(precip_hr, 24, 24)
wind_in = daily_mean(windsp_hr, 24, 24)
coa_in = rep(360,length.out=length(wind_in))

par(mfrow=c(3,4))
plot(sat_avg_in)
plot(rh_avg_in)
plot(sat_max_in)
plot(rh_max_in)
plot(sat_min_in)
plot(rh_min_in)
plot(swrad_in)
plot(ppt_in)
plot(wind_in)
plot(coa_in)

# make up some data
latitude = 35.97
days_in_yr = 365.25

# create the shared object
system('R CMD SHLIB -o weather_generator.so weather_generator.f90 wg_interface.f90')
# load the shared object
dyn.load("weather_generator.so")
# call function
tmp=.Fortran("weathergeneratorinterface",latitude=as.single(latitude),nos_days=as.integer(length(time)),days_in_yr=as.single(days_in_yr)
					,time=as.integer(time)
					,sat_avg_in=as.single(sat_avg_in),sat_max_in=as.single(sat_max_in),sat_min_in=as.single(sat_min_in)
					,ppt_in=as.single(ppt_in),swrad_in=as.single(swrad_in),coa_in=as.single(coa_in)
					,rh_avg_in=as.single(rh_avg_in),rh_max_in=as.single(rh_max_in),rh_min_in=as.single(rh_min_in),wind_in=as.single(wind_in)
          ,sat_out=as.single(array(0,dim=c(length(sat_avg_in)*24))),ppt_out=as.single(array(0,dim=c(length(sat_avg_in)*24)))
					,swrad_out=as.single(array(0,dim=c(length(sat_avg_in)*24))),coa_out=as.single(array(0,dim=c(length(sat_avg_in)*24)))
					,rh_out=as.single(array(0,dim=c(length(sat_avg_in)*24))),wind_out=as.single(array(0,dim=c(length(sat_avg_in)*24))) )

# extract wanted output
sat_out=tmp$sat_out ; ppt_out=tmp$ppt_out
coa_out=tmp$coa_out ; swrad_out=tmp$swrad_out
rh_out=tmp$rh_out   ; wind_out=tmp$wind_out
# unload the object
dyn.unload("weather_generator.so")
# and clean up
rm(tmp) ; gc()

par(mfrow=c(2,3))
plot(sat_out,col="red", type="l") ; lines(airt_hr,col="blue")
plot(ppt_out,col="red", type="l") ; lines((precip_hr),col="blue")
plot(coa_out,col="red", type="l")
plot(swrad_out,col="red", type="l") ; lines(swrad_hr,col="blue")
plot(rh_out,col="red", type="l") ; lines(rh_hr,col="blue")
plot(wind_out,col="red", type="l") ; lines(windsp_hr,col="blue")

summary(sat_out) ; summary(airt_hr)
summary(ppt_out) ; summary(precip_hr)
summary(swrad_out) ; summary(swrad_hr)
summary(rh_out) ; summary(rh_hr)

((sum(sat_out)-sum(airt_hr))/sum(airt_hr))*100
((sum(ppt_out)-sum(precip_hr))/sum(precip_hr))*100
((sum(swrad_out)-sum(swrad_hr))/sum(swrad_hr))*100
((sum(rh_out)-sum(rh_hr))/sum(rh_hr))*100

# 80 with sunset; 85 with out sunset ; 85 with 23.94 # 86.99 # interp with sunset
summary(lm(swrad_hr[1:length(swrad_out)]~swrad_out))

par(mfrow=c(1,1))
plot(swrad_out[1:48],col="red", type="l") ; lines(swrad_hr[1:48],col="blue")
lines(swrad_out[1:48],col="yellow")
