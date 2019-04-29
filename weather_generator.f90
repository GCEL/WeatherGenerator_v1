module weather_generator

  implicit none

  private ! assume all private unless specifically stated

  public :: hours_per_day, daily_to_hourly

  real,parameter :: seconds_per_hour = 3600.0
  real,parameter :: hours_per_day    = 24.0
  real,parameter :: lag_in_max_T     = 1.0  ! Nos of hrs after Noon that max temp occurs.
  real,parameter :: precip_period    = 1.0  ! Nos of hrs overwhich rainfall occurs.
  real,parameter :: pi = 3.14159265
  real,parameter :: earth_axis_tilt  = 23.4 ! tilt of earth's axis degrees

  ! variables that we want saved b/wn calls, e.g. to ensure matching b/wn days
  real           :: sat_at_midnight = -999. ! A way to make sure each day's T-curves match..
  real           :: rh_at_midnight  = -999. ! As above for RH

  save

contains
  !
  !----------------------------------------------------------------------------
  !
  subroutine daily_to_hourly( latitude_degs, day_number, days_in_year, & ! intent(in)
                                    co2_day,   tavg_day,     tmax_day, & ! intent(in)
                                   tmin_day,    ppt_day,       sw_day, & ! intent(in)
                                     rh_day, rh_max_day,   rh_min_day, & ! intent(in)
                                   wind_day,                           & ! intent(in)
                                     co2_hr,     sat_hr,       ppt_hr, & ! intent(out)
                                      sw_hr,      rh_hr,      wind_hr  ) ! intent(out)

    !> description <!

    implicit none

    ! arguments..
    integer,intent(in) :: day_number
    real,   intent(in) :: days_in_year, latitude_degs, co2_day, ppt_day, rh_day, &
         rh_min_day, rh_max_day, sw_day,     &
         Tavg_day, Tmax_day, Tmin_day, wind_day
    real,dimension(nint(hours_per_day)),&
         intent(out) :: co2_hr, ppt_hr, sat_hr, sw_hr, rh_hr, wind_hr

    ! local variable..
    real :: declination, & !
              daylength, & ! day length in hours
               latitude, & ! site latitude in radians
                sunrise    ! hour of sunrise

    ! Convert latitude into radians..
    latitude = latitude_degs * ( pi / 180.0 )
    ! Start with establishing the daylength and sunrise, as we use this throughout..
    call sunrise_daylength( latitude , day_number , days_in_year , declination , sunrise , daylength )

    ! Check whether in polar winter (ie sun never rises)..
    if ( sunrise .ge. (0.5 * hours_per_day) ) then

       ! No sunlight, so...
       sw_hr  = 0.0
       sat_hr = Tavg_day
       rh_hr  = rh_day

    else

      ! the solar daily cycle..
      call calc_sw( day_number , days_in_year , latitude , sw_day , declination , &
                                                       sunrise , daylength , sw_hr )

      ! then calculate the daily temperature cycle...
      call calc_temp( Tmax_day , Tmin_day , Tavg_day , sunrise , sat_hr )

      ! Then I think relative humidity came into the picture? (from which we work out vpd..)
      ! Ideally we want the daily-max/min of RH, in order to create a curve.
      ! Assign max to time of Tmin (ie dawn), and min to time of Tmax (ie noon+T-lag)
      call calc_rh( rh_max_day , rh_min_day , rh_day , sunrise , rh_hr )

    endif

    ! calculate daily cycle in rain - assume all falls just after dawn
    ! mean kg.m-2.s-1
    ppt_hr  = calc_ppt( ppt_day , sunrise )

    ! calculate daily cycle in CO2 - assumed constant at moment.
    co2_hr  = calc_co2( co2_day )

    ! calculate daily cycle in wind - assumed constant at moment.
    wind_hr = calc_wind( wind_day )


  end subroutine daily_to_hourly
  !
  !----------------------------------------------------------------------------
  !
  subroutine sunrise_daylength( latitude , day_number , days_in_year , declination , sunrise , daylength )

    ! Calculate the solar declination (in radians), and !
    ! then the daylength and from that sunrise.         !

    implicit none

    ! arguments..
    integer,intent(in) :: day_number
    real,   intent(in) :: days_in_year, & ! number of days in year
                              latitude    ! site latitude in radians
    real,  intent(out) :: daylength, & ! day length in hours
                        declination, & !
                            sunrise    ! hour of sunrise

    ! local variables..
    real :: calcs, day_angle
    integer :: steps_per_day

    steps_per_day = nint(hours_per_day)

    ! how far round a sun-earth orbit (equation A.5
    ! in Hartman's Global Physical Climatology)..
    day_angle = ( 2.0 * pi * day_number ) / days_in_year

    ! Solar declination..
    declination = calc_declination( day_angle )

    ! Daylength in hours..
    calcs = -tan(latitude) * tan(declination)
    if ( ( calcs .lt. -1.0  .or.  calcs .gt. 1.0 ) .and. declination < 0.0) then
      ! Days that have no sunlight at all...
      daylength = 0.0
    else
      daylength = min(hours_per_day,hours_per_day * acos( calcs ) / pi) ! 23.94
    endif

    ! Time of sunrise..
    !( 0hr=midnight, 12=noon, 'day' is evenly spread either side of noon)
    sunrise =  (hours_per_day * 0.5) - ( 0.5 * daylength )

  end subroutine sunrise_daylength
  !
  !----------------------------------------------------------------------------
  !
  function calc_declination( day_angle )

    ! Method taken from Appendix A of Hartmann's !
    ! Global Physical Climatology (1994). Uses a !
    ! truncated Fourier series approximation.    !

    implicit none

    ! arguments..
    real,intent(in) :: day_angle

    ! function result..
    real :: calc_declination

    ! local variables..
    real :: A(4), B(4)

    ! Fourier series coefficients..
    A(1) =  0.006918       ;   B(1) = 0.0
    A(2) = -0.399912       ;   B(2) = 0.070257
    A(3) = -0.006758       ;   B(3) = 0.000907
    A(4) = -0.002697       ;   B(4) = 0.001480

    ! actual calculation (equation A.6)..
    calc_declination = A(1) &
                     + A(2) * cos( 1.0 * day_angle ) &
                      + B(2) * sin( 1.0 * day_angle ) &
                     + A(3) * cos( 2.0 * day_angle ) &
                      + B(3) * sin( 2.0 * day_angle ) &
                     + A(4) * cos( 3.0 * day_angle ) &
                      + B(4) * sin( 3.0 * day_angle )

  end function calc_declination
  !
  !----------------------------------------------------------------------------
  !
  subroutine calc_sw( day_number, days_in_year , latitude , sw_day , declination , sunrise , daylength , sw_hr )

    ! This routine estimates the daily solar cycle.                  !
    ! Input (sw_day) is the daily-mean downward sw at the surface.   !
    ! The method is to follow Harmann (pg 29, Global Physical        !
    ! Climatology, 1994).  We calculate how the sunlight would be    !
    ! spread across the day (ie dependent upon latitude and !
    ! declination that would have been received under clear !
    ! skies is first calculated, according to..                      !
    ! q = cos( solar-zenith-angle ) !
    ! and then we rescale to match the daily average, which we know. !

    implicit none

    ! arguments..
    integer,intent(in)  :: day_number
    real,   intent(in)  :: daylength, & ! day length in hours
                         declination, & !
                            latitude, & !
                             sunrise, & ! hour of sunrise
                              sw_day, & !
                        days_in_year
    real,   intent(out) :: sw_hr(nint(hours_per_day))

    ! local variables..
    integer :: hr,minhr,maxhr
    real    :: cos_solar_zenith_angle(nint(hours_per_day)), day_angle, hour_angle

    ! initialise..
    sw_hr = 0.0
    cos_solar_zenith_angle = 0.0

    ! how far round a sun-earth orbit..
    day_angle = ( 2.0 * pi * day_number ) / days_in_year

    ! Now calculate the solar-zenith-angle, as per..
    !  sin(latitude)sin(declinaion) + cos(latitude)cos(declination)cos(hour-angle)
    ! where hour-angle = ( hr - 12 ) * 15. * ( pi / 180.d0 )
    ! and from that the solar flux.

    do hr = 1 , nint(hours_per_day)

      ! only calculate when sun is up..
      if ( ( real(hr) .gt. sunrise ) .and. ( real(hr) .lt. sunrise+daylength ) ) then

         hour_angle = real( hr - nint(hours_per_day * 0.5))  * 15. * ( pi / 180.0 )

         ! eqn 2.15 in Hartmann's Global Physical Climatology...
         cos_solar_zenith_angle(hr) = sin(latitude) * sin(declination) &
                   + cos(latitude) * cos(declination) * cos(hour_angle)
         ! sanity check..
         if ( cos_solar_zenith_angle(hr) .lt. 0.0 ) then
           print*,"The cosine of the solar zenith angle is less than 0 (it shouldn't be!)."
           print*,"At step ",hr," of ",hours_per_day,"on day ",day_number," at latitude ",latitude*(180/pi), &
                      "the cosine of the solar zenith angle is ",cos_solar_zenith_angle(hr)
!           stop "Stopping!"
         endif

       endif
     enddo

     ! This is what we would have received if there had 1Wm/2 of incoming solar,
     ! but actually we know that the number was much higher, and that it averaged
     ! to sw_day, so scale..
     ! sw_hr = hr_ratio * scale_factor, where scale_factor = 24*sw-day / sum(hr-ratio)
     minhr = max(1,ceiling(sunrise))
     maxhr = min(real(hours_per_day),floor(sunrise+daylength))
     sw_hr( minhr:maxhr ) = max( ( ( hours_per_day * sw_day ) / sum(cos_solar_zenith_angle( minhr:maxhr )) ) * &
                                 cos_solar_zenith_angle( minhr:maxhr ) , 0.0 )

     ! catch the error if cloudy day is sunnier that clear-sky day!
     if ( abs( (hours_per_day*sw_day) - sum(sw_hr)) .gt. 1.0 ) then
         print*,"ERROR! More sunlight on cloudy day than in our clear-sky calculations!"
         print*,"Clear sky SW radiation ",24.0*sw_day,"Cloudy (input) SW radiation ",sum(sw_hr)
         print*,"Sunrise ",sunrise,"day length ",daylength,"Latitude ",latitude*(180.0/pi),"Day number ",day_number
!         stop "Stopping"
     endif

  end subroutine calc_sw
  !
  !----------------------------------------------------------------------------
  !
  subroutine calc_temp( tmax_day , tmin_day , tavg_day , sunrise , t_hr )

    ! Create a curve of hourly temperatures from daily data. !
    ! Minimum is assumed to occur at dawn.                   !
    ! Maximum is assumed to occur at noon + T-lag.           !
    ! The midnight-temperature from each day is remembered,  !
    !  to ensure the curves match from one day to the next.  !
    ! During polar winter (ie no sunrise all day) then we    !
    !  assume the temperature is equal to tavg_day all day.  !

    implicit none

    ! arguments..
    real,intent(in)  :: tmax_day, tmin_day, tavg_day, & ! daily max / min / mean temperature (oC)
                        sunrise ! hour of sunrise
    real,intent(out) :: t_hr(nint(hours_per_day))

    ! local variables..
    real    :: hottest_time, sum_one, sum_two
    integer :: index_hot_time, hr

    ! if very first day, set the midnight-temp to equal the daily minimum..
    if ( sat_at_midnight .eq. -999. ) sat_at_midnight = 0.5*(tavg_day + tmin_day)
!TLS    sat_at_midnight = 0.5 *(tavg_day+tmin_day)

    ! We break the daily curve into three parts:
    !  midnight--dawn, dawn--maximum, maximum--midnight.

    do hr = 1 , int(sunrise)
      ! straight line from yesterday's midnight to today's t-min (at sunrise).
      t_hr( hr ) =  ( 1 / ( sunrise - 0.0 ) ) * &
               ( tmin_day * ( hr - 0.0 ) + sat_at_midnight * ( sunrise - hr ) )
    enddo

    hottest_time = (hours_per_day*0.5) + lag_in_max_T
    index_hot_time = int(hottest_time)
    do hr = int(sunrise)+1 , index_hot_time
      ! straight line from Tmin to Tmax
      t_hr( hr ) = ( 1 / ( hottest_time - sunrise ) ) * &
               ( tmax_day * ( hr - sunrise ) + tmin_day * ( hottest_time - hr ) )
    enddo

    ! calculate the new midnight-temperature, based on maintaining the daily-avg temp..
    sum_one = 0.0 ; sum_two = 0.0
    do hr = index_hot_time+1 , nint(hours_per_day)
      sum_one = sum_one + ( hr - hottest_time )
      sum_two = sum_two + ( hours_per_day - hr )
    enddo
    sat_at_midnight = ( 1 / sum_one) * ( &
            ( hours_per_day - hottest_time ) * ( hours_per_day * Tavg_day - sum(t_hr(1:index_hot_time) ) ) &
            - tmax_day * sum_two )

    ! Now calculate those post-PeakT hours..
    do hr = int(hottest_time)+1, hours_per_day
      ! straight line from Tmax downwards, in such a way as to maintain Tavg.
      t_hr( hr ) = ( 1 / ( hours_per_day - hottest_time ) ) * &
               ( sat_at_midnight * ( hr - hottest_time ) + tmax_day * ( hours_per_day - hr ) )
    enddo

  end subroutine calc_temp
  !
  !----------------------------------------------------------------------------
  !
  subroutine calc_rh( rh_max , rh_min , rh_avg , sunrise , rh_hr )

    ! We can get hold of specific humidity, from which we should  !
    ! be able to calculate the actual amounts of water correctly, !
    ! and then regain relative humidity.                          !
    ! But for now, we simply mirror the T-curve..                 !
    !(note polar-winter => constant RH all day)                   !

    implicit none

    ! arguments..
    real,intent(in)  :: rh_min, rh_max, rh_avg, & ! min / max / mean daily relative humidity (0-1)
                       sunrise ! hour of sunrise
    real,intent(out) :: rh_hr(nint(hours_per_day))

    ! local variables..
    integer :: index_hot_time, hr
    real    :: hottest_time, sum_one, sum_two

    ! if very first day, set the midnight-RH to equal the daily maximum..
    if ( rh_at_midnight .eq. -999. ) rh_at_midnight = rh_max
!TLS    rh_at_midnight = rh_max

    ! We break the daily curve into three parts:
    !  midnight--dawn, dawn--maximum, maximum--midnight.

    do hr = 1 , int(sunrise)
      ! straight line from yesterday's midnight to today's RH-max (at sunrise).
      rh_hr( hr ) =  ( 1 / ( sunrise - 0.0 ) ) * &
                        ( rh_max * ( hr - 0.0 ) + rh_at_midnight * ( sunrise - hr ) )
      ! sanity check however
      rh_hr( hr ) = min(1.0,max(0.0,rh_hr( hr )))
    enddo

    hottest_time = (hours_per_day*0.5) + lag_in_max_T
    index_hot_time = int(hottest_time)
    do hr = int(sunrise)+1 , index_hot_time
      ! straight line from Tmin to Tmax
      rh_hr( hr ) = ( 1 / ( hottest_time - sunrise ) ) * &
                       ( rh_min * ( hr - sunrise ) + rh_max * ( hottest_time - hr ) )
      ! sanity check however
      rh_hr( hr ) = min(1.0,max(0.0,rh_hr( hr )))
    enddo

    ! calculate the new midnight-RH, based on maintaining the daily-avg RH..
    sum_one = 0.0 ; sum_two = 0.0
    do hr = index_hot_time+1 , nint(hours_per_day)
      sum_one = sum_one + ( hr - hottest_time )
      sum_two = sum_two + ( hours_per_day - hr )
    enddo
    rh_at_midnight = ( 1 / sum_one) * ( &
                   ( hours_per_day - hottest_time ) * ( hours_per_day * rh_avg - sum(rh_hr(1:index_hot_time) ) ) &
                   - rh_min * sum_two )

    ! Now calculate those post-PeakT hours..
    do hr = int(hottest_time)+1, nint(hours_per_day)
      ! straight line from RH_min upwards, in such a way as to maintain RH-avg.
      rh_hr( hr ) = ( 1 / ( hours_per_day - hottest_time ) ) * &
                       ( rh_at_midnight * ( hr - hottest_time ) + rh_min * ( hours_per_day - hr ) )
      ! sanity check however
      rh_hr( hr ) = min(1.0,max(0.0,rh_hr( hr )))
    enddo

  end subroutine calc_rh
  !
  !----------------------------------------------------------------------------
  !
  function calc_ppt( ppt_day, sunrise )

    ! Calculate the daily-cycle in precipitation.                       !
    ! Generally Mat assumes it falls in the first few hours after dawn. !
    ! There is a possibility that knowing whether we are looking at     !
    ! stratiform vs convective might tell us something more.            !

    implicit none

    ! argument..
    real,intent(in) :: ppt_day, & !
                       sunrise    ! hour of sunrise

    ! function declaration..
    real,dimension(nint(hours_per_day)) :: calc_ppt

    ! local variable..
    integer :: hr, start, finish

    calc_ppt = 0.0

    ! ppt_day is the daily-avg (kg.m-2-s-1).  If we only rain for 5hrs (e.g. b/wn 5am and 9am) then
    ! we need hourly-avg = 24/5 * daily-avg to conserve the daily

    ! Only rain between 5am and 9am...
!    do hr = 6 , 10
!      calc_ppt( hr ) = ( 24.0 / 5.0 ) * (ppt_day*3600)
!    enddo

    ! estimate start and finish points for precipitation
    start = nint(sunrise) ; start = max(1, start)
    finish = start + nint(precip_period) ; finish = min(nint(hours_per_day), finish)

    ! Only rain between 5am and 9am...
    do hr = start , finish
       calc_ppt( hr ) = ( hours_per_day / real(finish-start+1) ) * (ppt_day*seconds_per_hour)
    enddo

  end function calc_ppt
  !
  !----------------------------------------------------------------------------
  !
  function calc_co2( co2_day )

    ! Calculate the daily-cycle in CO2 !

    implicit none

    ! argument..
    real,intent(in) :: co2_day

    ! function declaration..
    real,dimension(nint(hours_per_day)) :: calc_co2

    ! local variable..
    integer :: hr

    do hr = 1 , nint(hours_per_day)
      calc_co2( hr ) = co2_day
    enddo

  end function calc_co2
  !
  !----------------------------------------------------------------------------
  !
  function calc_wind( wind_day )

    ! Calculate the daily-cycle in wind speed !

    implicit none

    ! argument..
    real,intent(in) :: wind_day

    ! function declaration..
    real,dimension(nint(hours_per_day)) :: calc_wind

    ! local variable..
    integer :: hr

    do hr = 1 , nint(hours_per_day)
      calc_wind( hr ) = wind_day
    enddo

  end function calc_wind
  !
  !----------------------------------------------------------------------------
  !
end module weather_generator
