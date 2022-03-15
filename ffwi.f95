!
! Name: Lauren Lobo
! ID: 1087364
! Date: 2022-02-04
! Modern Fortran program of 1976 Forest Fire Weather Index program using Fortran f95 standards
!

program main
use FFWIndices
    implicit none

    ! Fine Fuel Moisture Code variables
    real :: today_final_moisture,correct_rain_route,temp_c,start_moisture,prev_ffmc,emc_drying
    real :: emc_wetting,logdry_rate,inter_x,ffmc
    
    ! Duff Moisture Code variables
    real :: drying_factor,dmc_afterrain,moisture_content,rain_effect,effective_rain,dc
    
    ! Drought Code Variables
    real ::yest_drought,drying_factor_drought,moisture_equiv_drought,drought_afterrain
    
    ! Indexes Variables
    real :: todays_ffmc,ffm_function,bui_less_dmc,dmc_bui_less_dmc,inter_fwi,log_fwi
    
    ! Extra variables 
    real :: wmr,ffm,dmc,si,bui,fwi,r,after_rain,ra,start_dmc,rain, tx,humidity,wind_speed
   
    character (len = 50) :: input_filename
    character (len = 50) :: output_filename
    integer :: l,j,ffm_moisture_after,ndays,nn,idays,i,idc,iffm,idmc,isi,ibui,ifwi,ih,iw
    integer , dimension(12) :: month_length
    real,dimension(12) :: daylength, month_daylength

        !FORMATS FOR PRINTING
1004    format(2X(/),'PROGRAM NO.: F-4 Forest Fire Weather Index')
1005    format(1X,' FFMC: Fine Fuel Moisture Code, DMC: Duff Moisture Code, DC: Drought Code')
1006    format(1X,' ISI: Initial Spread Index, BUI: Build Up Index, FWI: Fire Weather Index')
100     format(I2,F4.1,F4.1)
101     format(F4.1,I4,I4,F4.1)
102     format(F4.1,F4.1,F5.1,I2,I2)
1002    format(1(/) 1X,'  DATE  TEMP  RH   WIND  RAIN   FFMC   DMC   DC     ISI   BUI   FWI')
1003    format( 1X,'  mm/dd (C)   (%) (km/h) (mm)'/)

    ! File input and output names
    call file_input(input_filename)
    open(1,file=input_filename,status='old') 
    call file_output(output_filename)
    open(2,file=output_filename,status='unknown') 

    ! print out title and legend for table
    write(2,1004)
    write(2,1005)
    write(2,1006)

    ! Loop to read initial values of months and day lengths for arrays
    j=1
    do while(j <= 12)
        read(1,100) month_length(j), daylength(j),month_daylength(j)
        j = j + 1
    end do
    
    !    read in the previous days rain, prvious days duff moisture, yesterdays drought, fine fuel moisture after drying, and number of days
    read(1,102) prev_ffmc,start_dmc,yest_drought,ffm_moisture_after,ndays
    
    do j = ffm_moisture_after, 12
        nn = month_length(j)
        if (j == ffm_moisture_after ) then
            idays = month_length(j) - ndays + 1
        else
            idays = 1
        end if
        
        ! read weather data
        l=0 
        i = idays

        do i =idays, nn
            l = l + 1
            read(1,101,END=2000) temp_c,ih,iw,r
            if (l == 1) then
                ! titles and units of table
                write(2,1002)  
                write(2,1003)  
            end if
                ! collecting data from input of each days weather tempurature (C), humidity(%), wind speed(km/h), and rain(mm)
                tx = temp_c
                humidity = ih
                wind_speed = iw
                rain = r
            
            !   FINE FUEL MOISTURE CODE
            call FFMC_func(today_final_moisture,correct_rain_route,temp_c,start_moisture,prev_ffmc,emc_drying,logdry_rate, &
            inter_x,ffmc,r,humidity,wind_speed,after_rain,emc_wetting,ffm,ra)
            
            !     DUFF MOISTURE CODE 
            call DUFF_func(wmr,humidity,r,start_dmc,temp_c,drying_factor,dmc_afterrain,moisture_content,rain_effect, &
            effective_rain,dmc,daylength,j,ra)

            !        DROUGHT CODE 
            call DROUGHT_func(j,month_daylength,dc,r,ra,temp_c,yest_drought,drying_factor_drought,moisture_equiv_drought, &
            drought_afterrain,effective_rain)

            !  INITIAL SPREAD INDEX, BUILDUP INDEX, FIRE WEATHER INDEX
            call ISI_BUI_FWI_func(todays_ffmc,ffm_function,bui_less_dmc,dmc_bui_less_dmc,inter_fwi,log_fwi,si,bui,fwi,idc,iffm,&
            idmc,isi,ibui,ifwi,ffm,wind_speed,dc,dmc)
            
            ! write to file
            write(2,1001) j,i,tx,ih,iw,rain,iffm,idmc,idc,isi,ibui,ifwi
1001        format(1X,2I3,F6.1,I4,I6,F7.1,6I6)

            ! set variables for previous day before loop continues to next day
            prev_ffmc = ffm
            start_dmc = dmc
            yest_drought = dc

        end do
    end do

2000    stop 
end program main