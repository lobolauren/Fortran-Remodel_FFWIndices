
!
! Name: Lauren Lobo
! ID: 1087364
! Date: 2022-02-04
! file and calculation functions for ffwi program based on equations in Pub.1333
!
module FFWIndices
    implicit none
    contains
    subroutine file_input(input_filename)
            character (len = 50), intent(inout) :: input_filename
103     format(1X,'Enter Input Filename: ')
104     format(10A)

            write(*,103)
            read(*,104) input_filename
    end subroutine file_input

    subroutine file_output(output_filename)
            character (len = 50), intent(inout) :: output_filename
104     format(10A)
105     format(1X,'Enter Output Filename: ')

            write(*,105)
            read(*,104) output_filename
    end subroutine file_output

    subroutine FFMC_func(today_final_moisture,correct_rain_route,temp_c,start_moisture,prev_ffmc,emc_drying,logdry_rate, &
        inter_x,ffmc,r,humidity,wind_speed,after_rain,emc_wetting,ffm,ra)
        real, intent(inout) :: today_final_moisture,correct_rain_route,temp_c,start_moisture,prev_ffmc,emc_drying
        real, intent(inout) :: logdry_rate,inter_x,ffmc,after_rain,emc_wetting,ffm,r,humidity,wind_speed,ra
        
        if (r <= 0.5) then
            r = 0.0
            after_rain = prev_ffmc   !set previous days rain to be after_rain 
        else
            ra = r 
            if (ra > 1.45) then
                if (ra-5.75 < 0) then
                    ffmc=57.87-(18.2*alog(ra-1.016)) !equation 10a
                end if
                if (ra-5.75 == 0) then
                    ffmc=57.87-(18.2*alog(ra-1.016)) !equation 10b
                end if
                if (ra-5.75 > 0) then
                    ffmc=40.69-(8.25*alog(ra-1.905)) !equation 10c
                end if
            else 
                ffmc = 123.85-(55.6*alog(ra+1.016))
            end if
            correct_rain_route = 8.73*exp(-0.1117*prev_ffmc) ! C by equation 11
            after_rain = (prev_ffmc/100.)*ffmc+(1.0-correct_rain_route) !ffmv after rain by equation 9
            if (after_rain < 0) then
                after_rain=0.0
            end if
        end if

        start_moisture = 101. -after_rain 
        emc_drying = 0.942*(humidity**0.679)+(11.*exp((humidity-100.)/10.))+0.18*(21.1-temp_c) &
        *1*(1.-1./exp(0.115*humidity)) ! equation 11

        if (start_moisture - emc_drying < 0) then
            emc_wetting = 0.618*(humidity**0.753)+(10.*exp((humidity-100.)/10.))+0.18*(21.1-temp_c) &
            *1*(1.-1./exp(0.115*humidity)) ! equation 2a
            if (start_moisture > emc_wetting) then
                today_final_moisture = start_moisture
            else
                today_final_moisture = emc_wetting -(emc_wetting-start_moisture)/1.9953 !equation 2b
            end if                
        else if (start_moisture - emc_drying == 0) then
            today_final_moisture = start_moisture
        else 
            inter_x = 0.424*(1.-(humidity/100.)**1.7)+(0.0694*(wind_speed**0.5))*(1.-(humidity/100.)**8)
            logdry_rate = inter_x*(0.463*(exp(0.0365*temp_c)))
            today_final_moisture = emc_drying+(start_moisture-emc_drying)/10.**logdry_rate
        end if
        ffm=101.- today_final_moisture ! todays FFMC
        if (ffm <= 101) then
            if (ffm < 0 ) then
                ffm = 0.0
            end if
        else
            ffm = 101
        end if
        end subroutine FFMC_func


    subroutine DUFF_func(wmr,humidity,r,start_dmc,temp_c,drying_factor,dmc_afterrain,moisture_content,rain_effect, &
        effective_rain,dmc,daylength,j,ra)
        real, intent(inout) ::wmr,humidity,r,start_dmc, temp_c,drying_factor,dmc_afterrain,moisture_content,rain_effect, &
        effective_rain,dmc,ra
        integer,intent(in) :: j
        real,dimension(12),intent(inout) :: daylength
        
        if (temp_c + 1.1 < 0) then
            temp_c =-1.1
        end if
        drying_factor = 1.894 * (temp_c + 1.1) * (100. - humidity) * (daylength(j) * 0.0001)
        if (r > 1.5) then
            ra = r
            effective_rain = 0.92 * ra - 1.27 !equation 13
            moisture_content = 20.0+280./exp(0.023*start_dmc) !equation 12
            if(start_dmc <=33) then
                rain_effect = 100./(0.5+0.3*start_dmc) !equation 15 a
            else if (start_dmc - 65 <= 0) then
                rain_effect = 14.-1.3*alog(start_dmc) !equation 15 b
            else
                rain_effect = 6.2*alog(start_dmc)-17.2 !equation 15 c
            end if
            wmr = moisture_content+(1000.*effective_rain)/(48.77+rain_effect*effective_rain) ! equation 14 to find todays misture
            dmc_afterrain = 43.43*(5.6348-alog(wmr-20.))
        else
            dmc_afterrain = start_dmc
        end if 
        if (dmc_afterrain < 0) then
            dmc_afterrain = 0.0
        end if
        dmc = dmc_afterrain + drying_factor !todays duff moisture
    end subroutine DUFF_func


    subroutine DROUGHT_func(j,month_daylength,dc,r,ra,temp_c,yest_drought,drying_factor_drought, &
        moisture_equiv_drought,drought_afterrain,effective_rain)
        real, intent(inout) :: dc,r,ra, temp_c,yest_drought,drying_factor_drought, &
        moisture_equiv_drought,drought_afterrain,effective_rain
        integer,intent(in) :: j
        real,dimension(12),intent(inout) :: month_daylength
        
        if (temp_c + 2.8 < 0) then
            temp_c =-2.8
        end if
        drying_factor_drought = (.36*(temp_c+2.8)+month_daylength(j))/2. ! set yesterdays drying factor as previous
        if (r > 2.8) then
            ra = r
            effective_rain = 0.83*ra-1.27
            moisture_equiv_drought = 800.* exp(-yest_drought/400.) !equation 19
            drought_afterrain = yest_drought-400.* alog(1.+((3.937*effective_rain)/moisture_equiv_drought))
            if (drought_afterrain <= 0) then
                drought_afterrain=0.0
            end if
            dc = drought_afterrain + drying_factor_drought !todays drought
            if(dc < 0) then
                dc=0.0
            end if 
        else
            drought_afterrain = yest_drought
            dc = drought_afterrain + drying_factor_drought !todays drought 
            if (dc <= 0) then
                dc = 0.0
            end if 
        end if
    end subroutine DROUGHT_func

    subroutine ISI_BUI_FWI_func( todays_ffmc,ffm_function,bui_less_dmc,dmc_bui_less_dmc,inter_fwi,log_fwi,si, &
        bui,fwi,idc,iffm,idmc,isi,ibui,ifwi,ffm,wind_speed,dc,dmc)
        real, intent(inout) ::todays_ffmc,ffm_function,bui_less_dmc,dmc_bui_less_dmc,inter_fwi,log_fwi,si,bui, &
        fwi,ffm,wind_speed,dc,dmc
        integer, intent(inout) :: idc,iffm,idmc,isi,ibui,ifwi
        
        ! Initial Spread Index code
        todays_ffmc = 101.- ffm
        ffm_function = 19.1152*exp(-0.1386*todays_ffmc)*(1.+todays_ffmc**4.65/7950000.)
        si = ffm_function*exp(0.05039*wind_speed)

        !buildup index
        bui = (0.8*dc*dmc)/(dmc+0.4*dc)
        if (bui < dmc) then
            bui_less_dmc = (dmc-bui)/dmc
            dmc_bui_less_dmc = 0.92+(0.0114*dmc)**1.7
            bui = dmc-(dmc_bui_less_dmc*bui_less_dmc)
            if (bui < 0) then
                bui=0
            end if
        end if

        ! fire weather index
        if (bui <= 80) then
            inter_fwi = 0.1*si*(0.626*bui**0.809+2.)
        else
            inter_fwi = 0.1*si*(1000./(25.+108.64/exp(0.023*bui)))
        end if
        if( inter_fwi-1.0 > 0) then
            log_fwi = 2.72*(0.43*alog(inter_fwi))**0.647
            fwi = exp(log_fwi)
        else
            fwi = inter_fwi
        end if
        
        !cast real to integers
        idc = int(dc+0.5)
        iffm = int(ffm+0.5)
        idmc = int(dmc+0.5)
        isi = int(si+0.5)
        ibui = int(bui+0.5)
        ifwi = int(fwi+0.5)
                
    end subroutine ISI_BUI_FWI_func
end module FFWIndices