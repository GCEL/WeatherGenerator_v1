
  !! This provides an interface b/wn the wrapper and whatever !!
  !! weather generator we wish to call.                       !!

  !****************************************

  subroutine weathergeneratorinterface( latitude   , nos_days   , & ! in
                                        days_in_yr , time       , & ! in
                                        sat_avg_in , sat_max_in , & ! in
                                        sat_min_in , ppt_in     , & ! in
                                        swrad_in   , coa_in     , & ! in
                                        rh_avg_in  , rh_max_in  , & ! in
                                        rh_min_in  , wind_in    , & ! in
                                        sat_out    , ppt_out    , & ! out
                                        swrad_out  , coa_out    , & ! out
                                        rh_out     , wind_out    )  ! out

    use weather_generator

    implicit none

    ! arguments..
    ! location and timing information
    integer,intent(in)              :: nos_days       !
    real,dimension(nos_days),intent(in) :: latitude & ! degrees
                                          ,time       ! day of time series
    real,intent(in)                 :: days_in_yr

    ! declare daily input variables
    real,dimension(nos_days),intent(in)    :: sat_avg_in & ! daily avg surface air temperature (oC)
                                             ,sat_max_in & ! daily max surface air temperature (oC)
                                             ,sat_min_in & ! daily min surface air temperature (oC)
                                             ,ppt_in     & ! mean (precip) (kg.m-2.s-1)
                                             ,swrad_in   & ! mean short wave in (W.m-2)
                                             ,coa_in     & ! mean CO2 (ppm)
                                             ,rh_avg_in  & ! daily avg rel humidity (frac)
                                             ,rh_max_in  & ! daily max rel humidity (frac)
                                             ,rh_min_in  & ! daily min rel humidity (frac)
                                             ,wind_in      ! mean wind (ms.-1)

    real,dimension(nos_days*nint(hours_per_day)),intent(inout) :: sat_out   & ! hourly surface air temperature (oC)
                                                ,ppt_out   & ! hourly precip (kg.m-2.s-1)
                                                ,swrad_out & ! hourly sw radiation (W.m-2)
                                                ,coa_out   & ! hourly CO2 (ppm)
                                                ,rh_out    & ! hourly rel humidity (frac)
                                                ,wind_out    !

    ! local variables..
    integer             :: i, day_number, a, b
    real, dimension(nint(hours_per_day)) :: coa_hrly, ppt_hrly, rh_hrly, sat_hrly, sw_hrly, wind_hrly

    ! calculate needed timing information
    a =  nint(hours_per_day + 1.0) ; b = nint(hours_per_day)

    do i = 1 , nos_days
       ! calculate day-of-year..
       day_number = max(1,nint(mod( time(i) , days_in_yr )))

       ! Call the weather-generator with the daily-means..
       call daily_to_hourly( latitude(i),   day_number,    days_in_yr,  &  ! in
                             coa_in(i),     sat_avg_in(i),       &         ! in
                             sat_max_in(i), sat_min_in(i),       &         ! in
                             ppt_in(i),     swrad_in(i),         &         ! in
                             rh_avg_in(i),  rh_max_in(i),        &         ! in
                             rh_min_in(i),  wind_in(i),          &         ! in
                             coa_hrly,      sat_hrly,            &         ! out
                             ppt_hrly,      sw_hrly,             &         ! out
                             rh_hrly,       wind_hrly            )         ! out

       ! Update the output variables..
         coa_out( (i-1)*a : i*b ) =  coa_hrly
         ppt_out( (i-1)*a : i*b ) =  ppt_hrly
          rh_out( (i-1)*a : i*b ) =   rh_hrly
         sat_out( (i-1)*a : i*b ) =  sat_hrly
       swrad_out( (i-1)*a : i*b ) =   sw_hrly
        wind_out( (i-1)*a : i*b ) = wind_hrly

    enddo

  end subroutine weathergeneratorinterface
