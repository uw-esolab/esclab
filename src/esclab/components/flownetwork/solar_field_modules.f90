module SF_props
    implicit none

        contains
        double precision function Density_SF(fnumd,T,P)
            !This function accepts as inputs temperature [K] and pressure [Pa]
            !This function outputs in units of [kg/m^3]
            double precision::xlo,xhi, Dens_fluid, Td,HTFPropsav
            double precision::T,P,fnumd
            !!double precision,dimension(size(fprop(1,:)))::dxx,dyy !Create dummy arrays
            integer::fnum,lb,ub,dum,t_warn
            !Density_SF=1.
            fnum=nint(fnumd)
            Td=T-273.15             !Convert from K to C
                
            select case(fnum)
            case(1)   !    1.) Air
            Density_SF = P/(287.*T)
            case(2)   !    2.) Stainless_AISI316
                Density_SF=8349.38 - 0.341708*T - 0.0000865128*T*T  !EES
            case(3)   !    3.) Water (liquid)
                Density_SF = 1000 
            case(4)   !    4.) Steam
                continue
            case(5)   !    5.) CO2
                continue
            case(6)   !    6.) Salt (68% KCl, 32% MgCl2)
            Density_SF = 1E-10*T*T*T - 3E-07*T*T - 0.4739*T + 2384.2
            case(7)   !    7.) Salt (8% NaF, 92% NaBF4)
            Density_SF = 8E-09*T*T*T - 2E-05*T*T - 0.6867*T + 2438.5
            case(8)   !    8.) Salt (25% KF, 75% KBF4)
            Density_SF = 2E-08*T*T*T - 6E-05*T*T - 0.7701*T + 2466.1
            case(9)   !    9.) Salt (31% RbF, 69% RbBF4)
            Density_SF = -1E-08*T*T*T + 4E-05*T*T - 1.0836*T + 3242.6
            case(10)   !    10.) Salt (46.5% LiF, 11.5%NaF, 42%KF)
            Density_SF =  -2E-09*T*T*T + 1E-05*T*T - 0.7427*T + 2734.7
            case(11)   !    11.) Salt (49% LiF, 29% NaF, 29% ZrF4)
            Density_SF = -2E-11*T*T*T + 1E-07*T*T - 0.5172*T + 3674.3
            case(12)   !    12.) Salt (58% KF, 42% ZrF4)
            Density_SF =  -6E-10*T*T*T + 4E-06*T*T - 0.8931*T + 3661.3
            case(13)   !    13.) Salt (58% LiCl, 42% RbCl)
            Density_SF = -8E-10*T*T*T + 1E-06*T*T - 0.689*T + 2929.5
            case(14)   !    14.) Salt (58% NaCl, 42% MgCl2)
            Density_SF = -5E-09*T*T*T + 2E-05*T*T - 0.5298*T + 2444.1
            case(15)   !    15.) Salt (59.5% LiCl, 40.5% KCl)
            Density_SF = 1E-09*T*T*T - 5E-06*T*T - 0.864*T + 2112.6
            case(16)   !    16.) Salt (59.5% NaF, 40.5% ZrF4)
            Density_SF =  -5E-09*T*T*T + 2E-05*T*T - 0.9144*T + 3837.
            case(17)   !    17.) Salt (60% NaNO3, 40% KNO3)
            Density_SF = dmax1(-1E-07*T*T*T + 0.0002*T*T - 0.7875*T + 2299.4,1000.d0)
            case(18)
            !Density_SF of Nitrate Salt, [kg/m3]
            Density_SF = dmax1(2090 - 0.636 * (T-273.15),1000.d0)
            case(19)
            !Density_SF of Caloria HT 43 [kg/m3]
            Density_SF = dmax1(885 - 0.6617 * Td - 0.0001265 * Td*Td,100.d0)
            case(20)
            !Density_SF of HITEC XL Nitrate Salt, [kg/m3]
            Density_SF = dmax1(2240 - 0.8266 * Td,800.d0)
            case(21)
            !Density_SF of Therminol Oil [kg/m3]
            Density_SF = dmax1(1074.0 - 0.6367 * Td - 0.0007762 * Td*Td,400.d0)
            case(22)
            !Density_SF of HITEC Salt, [kg/m3]
            Density_SF = dmax1(2080 - 0.733 * Td,1000.d0)
            case(23)
            !Density_SF of Dowtherm Q [kg/m3]
            Density_SF = dmax1(-0.757332 * Td + 980.787,100.d0)                               ! Russ 10-2-03
            case(24)
            !Density_SF of Dowtherm RP [kg/m3]
            Density_SF = dmax1(-0.000186495 * Td*Td - 0.668337 * Td + 1042.11,200.d0)  !Russ 10-2-03
            case(25)
            !Density_SF of HITEC XL Nitrate Salt, [kg/m^3]
            Density_SF = dmax1(2240 - 0.8266 * Td,800.d0)
            case(26) !Argon
            Density_SF = dmax1(P/(208.13*T),1.e-10)
            case(27) !Hydrogen
            Density_SF = dmax1(P/(4124.*T),1.e-10)
            case(28)    !T-91 Steel: "Thermo hydraulic optimisation of the EURISOL DS target" - Paul Scherrer Institut
            Density_SF = -0.3289*Td + 7742.5
            case(29)    !Therminol 66: Reference: Therminol Reference Disk by Solutia: http://www.therminol.com/pages/tools/toolscd.asp
            Density_SF = -0.7146*Td + 1024.8
            case(30)    !Therminol 59: Reference: Therminol Reference Disk by Solutia: http://www.therminol.com/pages/tools/toolscd.asp
            Density_SF = -0.0003*Td*Td - 0.6963*Td + 988.44
            case(31:35) 
            continue !no informaion
            !!case(36:) !Any integer greater than 35
            !!!Call the user-defined property table
            !!lb=fl_bounds(fnum-35)
            !!ub=fl_bounds(fnum-35+1)-1
            !!if(ub.lt.lb) ub=size(fprop(1,:))
            !!dxx(:)=fprop(1,lb:ub)
            !!dyy(:)=fprop(3,lb:ub)
            !!call interp(Td,size(dxx),dxx,dyy,Gjsav,Density_SF)
            !!if((Gjsav.eq.ub).or.(Gjsav.eq.lb)) dum=t_warn(Td,dxx(lb),dxx(ub),"User-specified fluid")
            ! case(36) !36-User defined SF HTF
            ! call NR_LINEAR_INTERPOLATION_00(Td,size(T_SF_HTF),T_SF_HTF,den_NN_SF_HTF,HTFPropsav,Density_SF) !Density_SF, Td is in [C], Density_SF in [kg/m3]
            ! case(37) !37-User defined TES HTF
            ! call NR_LINEAR_INTERPOLATION_00(Td,size(T_TES_HTF),T_TES_HTF,den_NN_TES_HTF,HTFPropsav,Density_SF) !Density_SF, Td is in [C], Density_SF in [kg/m3]
            case(40)
            !Density_SF of Dowtherm A [kg/m3]
            Density_SF = dmax1(1063.61 - 0.605235*Td - 0.000860877*Td*Td,400.d0)!Density_SF, Td is in [C], Density_SF in [kg/m3]
            end select

        end function


        !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        double precision function Viscosity_SF(fnumd,T,P)
        implicit none
        !This function accepts as inputs temperature [K] and pressure [Pa]
        !This function outputs in units of [Pa-s]
        double precision::Tx,xlo,xhi, Td,HTFPropsav
        double precision,intent(in)::T,P,fnumd
        !!double precision,dimension(size(fprop(1,:)))::dxx,dyy !Create dummy arrays
        integer::fnum,lb,ub,dum,t_warn
        Viscosity_SF=1.
        fnum=nint(fnumd)
        Td = T-273.15

        select case(fnum)
        case(1)   !    1.) Air
        Viscosity_SF = dmax1(0.0000010765 + 7.15173E-08*T - 5.03525E-11*T*T + 2.02799E-14*T*T*T,1.e-6)
        case(2)   !    2.) Stainless_AISI316
        continue
        case(3)   !    3.) Water (liquid)
        continue 
        case(4)   !    4.) Steam
        continue
        case(5)   !    5.) CO2
        continue
        case(6)   !    6.) Salt (68% KCl, 32% MgCl2)
        Viscosity_SF = .0146*exp(2230./T)*0.001 !convert cP to kg/m-s
        case(7)   !    7.) Salt (8% NaF, 92% NaBF4)
        Viscosity_SF = .0877*exp(2240./T)*0.001 !convert cP to kg/m-s
        case(8)   !    8.) Salt (25% KF, 75% KBF4)
        Viscosity_SF = .0431*exp(3060./T)*0.001 !convert cP to kg/m-s
        case(9)   !    9.) Salt (31% RbF, 69% RbBF4)
        Viscosity_SF = .0009
        case(10)   !    10.) Salt (46.5% LiF, 11.5%NaF, 42%KF)
        Viscosity_SF = .0400*exp(4170./T)*0.001 !convert cP to kg/m-s
        case(11)   !    11.) Salt (49% LiF, 29% NaF, 29% ZrF4)
        Viscosity_SF = .0069
        case(12)   !    12.) Salt (58% KF, 42% ZrF4)
        Viscosity_SF = .0159*exp(3179./T)*0.001 !convert cP to kg/m-s
        case(13)   !    13.) Salt (58% LiCl, 42% RbCl)
        Viscosity_SF = .0861*exp(2517./T)*0.001 !convert cP to kg/m-s          !
        case(14)   !    14.) Salt (58% NaCl, 42% MgCl2)
        Viscosity_SF = .0286*exp(1441./T)*0.001 !convert cP to kg/m-s
        case(15)   !    15.) Salt (59.5% LiCl, 40.5% KCl)
        Viscosity_SF = .0861*exp(2517./T)*0.001 !convert cP to kg/m-s          !
        case(16)   !    16.) Salt (59.5% NaF, 40.5% ZrF4)
        Viscosity_SF = .0767*exp(3977./T)*0.001 !convert cP to kg/m-s
        case(17)   !    17.) Salt (60% NaNO3, 40% KNO3)
        Tx=T-273.15  !This particular equation is in terms of degrees celsius
        Viscosity_SF = dmax1(-1.473302E-10*Tx**3 + 2.279989E-07*Tx**2 - 1.199514E-04*Tx + 2.270616E-02,.0001d0)
        case(18)
        !Absolute Viscosity_SF of Nitrate Salt, [Pa s]
        Viscosity_SF = dmax1((22.714 - 0.12 * Td + 0.0002281 * Td *Td - 0.0000001474 * Td*Td*Td) / 1000,1.e-6)
        !case(19)
        !Absolute Viscosity_SF of Caloria HT 43 [m2/s]
        !Viscosity_SF = (0.040439268 * max(10.d0,Td)**-1.946401872) * density(19.d0, T, 0.d0)
        case(20)  
        !Absolute Viscosity_SF of HITEC XL Nitrate Salt, [Pa s]
        Viscosity_SF = 1372000 * Td**-3.364
        case(21)
        !Absoute Viscosity_SF of Therminol Oil [Pa s]
        Viscosity_SF = 0.001 * (10**0.8703 * dmax1(Td,20.)**(0.2877 + Log10(dmax1(Td,20.)**-0.3638)))
        case(22)
        !Absolute Viscosity_SF of HITEC Salt, [Pa s]
        Viscosity_SF = dmax1(0.00622 - 0.0000102 * Td,1.e-6)
        case(23)
        !Absoute Viscosity_SF of Dowtherm Q [Pa s]
        Viscosity_SF = 1 / (132.40658 + 4.36107 * Td + 0.0781417 * Td*Td - 0.00011035416 * Td*Td*Td) !Hank 10-2-03
        case(24)
        !Absoute Viscosity_SF of Dowtherm RP [Pa s]
        Viscosity_SF = 1 / (4.523003 + 0.39156855 * Td + 0.028604206 * Td*Td)  !Hank 10-2-03
        case(25)
        !Absolute Viscosity_SF of HITEC XL Nitrate Salt, [Pa s]
        Viscosity_SF = 1372000 * Td**-3.364
        case(26)   !Argon 
        Viscosity_SF = 4.4997e-6 + 6.38920E-08*T - 1.24550E-11*T*T
        case(27)  !Hydrogen
        Viscosity_SF=0.00000231 + 2.37842E-08*T - 5.73624E-12*T*T
        case(28)
        continue
        case(29)    !Therminol 66: Reference: Therminol Reference Disk by Solutia: http://www.therminol.com/pages/tools/toolscd.asp
            IF(Td < 80.)THEN
                Viscosity_SF = 1.31959963 - 0.171204729*Td + 0.0100351594*Td**2 - 0.000313556341*Td**3 + 0.0000053430666*Td**4 - 4.66597650E-08*Td**5 + 1.63046296E-10*Td**6
            ELSE
                Viscosity_SF = 0.0490075884 - 0.00120478233*Td + 0.0000130162082*Td**2 - 7.58913847E-08*Td**3 + 2.47856063E-10*Td**4 - 4.26872345E-13*Td**5 + 3.01949160E-16*Td**6
            ENDIF
        case(30)    !Therminol 59: Reference: Therminol Reference Disk by Solutia: http://www.therminol.com/pages/tools/toolscd.asp
            IF (Td < 25.)THEN
                Viscosity_SF = 0.0137267822 - 0.000218740224*Td + 0.0000759248815*Td**2 - 0.00000473464744*Td**3 - 1.97083667E-07*Td**4 + 4.35487179E-09*Td**5 + 2.40243056E-10*Td**6
            ELSE
                Viscosity_SF = 0.0114608807 - 0.000313431056*Td + 0.00000416778121*Td**2 - 3.04668508E-08*Td**3 + 1.23719006E-10*Td**4 - 2.60834697E-13*Td**5 + 2.22227675E-16*Td**6
            ENDIF
        case(31:35)
        continue  !no information
        !!case(36:) !Any integer greater than 35
        !!!Call the user-defined property table
        !!lb=fl_bounds(fnum-35)
        !!ub=fl_bounds(fnum-35+1)-1
        !!if(ub.lt.lb) ub=size(fprop(1,:))
        !!dxx(:)=fprop(1,lb:ub)
        !!dyy(:)=fprop(4,lb:ub)
        !!call interp(Td,size(dxx),dxx,dyy,Gjsav,Viscosity_SF)
        !!if((Gjsav.eq.ub).or.(Gjsav.eq.lb)) dum=t_warn(Td,dxx(lb),dxx(ub),"User-specified fluid")

        ! case(36) !36-User defined SF HTF
        !     call NR_LINEAR_INTERPOLATION_00(Td,size(T_SF_HTF),T_SF_HTF,viscosity_SF_HTF,HTFPropsav,Viscosity_SF) !HTF Viscosity_SF, Td is in [C], Viscosity_SF in [Pa-s]
        ! case(37) !37-User defined TES HTF
        !     call NR_LINEAR_INTERPOLATION_00(Td,size(T_TES_HTF),T_TES_HTF,viscosity_TES_HTF,HTFPropsav,Viscosity_SF) !TES Viscosity_SF, Td is in [C], Viscosity_SF in [Pa-s]
        case(40)
        !Absoute Viscosity_SF of Dowtherm A [Pa-s]
        Viscosity_SF = 0.786512*dmax1(Td,20.)**-1.44263 !HTF Viscosity_SF, Td is in [C], Viscosity_SF in [Pa-s]
        end select

        continue
        end function

        double precision function spec_SF(fnumd,T,P)
            !This function accepts as inputs temperature [K] and pressure [Pa]
            !This function outputs in units of [kJ/kg-K]
            double precision::xlo,xhi, Td,HTFPropsav
            double precision,intent(in)::T,P,fnumd
            !!double precision,dimension(size(fprop(1,:)))::dxx,dyy !Create dummy arrays
            integer::fnum,lb,ub,dum,t_warn

            spec_SF=1.
            fnum=nint(fnumd)
            Td = T - 273.15
            select case(fnum)
            case(1)   !    1.) Air
                spec_SF = 1.03749 - 0.000305497*T + 7.49335E-07*T*T - 3.39363E-10*T*T*T
            !spec_SF = 1.03749 - 0.000305497*T + 7.49335E-07*T*T - 3.39363E-10*T*T*T
            case(2)   !    2.) Stainless_AISI316
                spec_SF = 0.368455 + 0.000399548*T - 1.70558E-07*T*T !EES
            case(3)   !    3.) Water (liquid)
                spec_SF = 4.181d0  !mjw 8.1.11 
            case(4)   !    4.) Steam
                continue
            case(5)   !    5.) CO2
                continue
            case(6)   !    6.) Salt (68% KCl, 32% MgCl2)
                spec_SF = 1.156
            case(7)   !    7.) Salt (8% NaF, 92% NaBF4)
                spec_SF = 1.507
            case(8)   !    8.) Salt (25% KF, 75% KBF4)
                spec_SF = 1.306
            case(9)   !    9.) Salt (31% RbF, 69% RbBF4)
                spec_SF = 9.127
            case(10)   !    10.) Salt (46.5% LiF, 11.5%NaF, 42%KF)
                spec_SF = 2.010
            case(11)   !    11.) Salt (49% LiF, 29% NaF, 29% ZrF4)
                spec_SF = 1.239
            case(12)   !    12.) Salt (58% KF, 42% ZrF4)
                spec_SF = 1.051
            case(13)   !    13.) Salt (58% LiCl, 42% RbCl)
                spec_SF = 8.918
            case(14)   !    14.) Salt (58% NaCl, 42% MgCl2)
                spec_SF = 1.080
            case(15)   !    15.) Salt (59.5% LiCl, 40.5% KCl)
                spec_SF = 1.202
            case(16)   !    16.) Salt (59.5% NaF, 40.5% ZrF4)
                spec_SF = 1.172
            case(17)   !    17.) Salt (60% NaNO3, 40% KNO3)
                spec_SF = -1E-10*T*T*T + 2E-07*T*T + 5E-06*T + 1.4387
            case(18) !Heat Capacity of Nitrate Salt, [J/kg/K]
                spec_SF = (1443. + 0.172 * (T-273.15))/1000.d0
            case(19)
            !Specific Heat of Caloria HT 43 [J/kgC]
                spec_SF = (3.88 * (T-273.15) + 1606.0)/1000.
            case(20)
            !Heat Capacity of HITEC XL Nitrate Salt, [J/kg/K]
                spec_SF = dmax1(1536 - 0.2624 * Td - 0.0001139 * Td * Td,1000.d0)/1000.
            case(21)
            !Specific Heat of Therminol Oil, kJ/kg/K
                spec_SF = (1.509 + 0.002496 * Td + 0.0000007888 * Td*Td)
            case(22)
            !Heat Capacity of HITEC Salt, [J/kg/K]
                spec_SF = (1560 - 0.0 * Td)/1000.
            case(23)
            !Specific Heat of Dowtherm Q, J/kg/K
                spec_SF = (-0.00053943 * Td*Td + 3.2028 * Td + 1589.2)/1000.               ! Russ 10-2-03
            case(24)
            !Specific Heat of Dowtherm RP, J/kg/K
                spec_SF = (-0.0000031915 * Td**2 + 2.977 * Td + 1560.8)/1000.       !Russ 10-2-03
            case(25)
            !Heat Capacity of HITEC XL Nitrate Salt, [J/kg/K]
                spec_SF = dmax1(1536 - 0.2624 * Td - 0.0001139 * Td * Td,1000.d0)/1000.
            case(26)    ! Argon
                spec_SF = 0.5203 !Cp only, Cv is different
            case(27)    ! Hydrogen
                spec_SF = dmin1(dmax1(-45.4022 + 0.690156*T - 0.00327354*T*T + 0.00000817326*T*T*T - 1.13234E-08*T*T*T*T + 8.24995E-12*T*T*T*T*T - 2.46804E-15*T*T*T*T*T*T,11.3d0),14.7d0)
            case(28)    !T-91 Steel: "Thermo hydraulic optimisation of the EURISOL DS target" - Paul Scherrer Institut
                spec_SF = 0.0004*Td*Td + 0.2473*Td + 450.08
            case(29)    !Therminol 66: Reference: Therminol Reference Disk by Solutia: http://www.therminol.com/pages/tools/toolscd.asp
                spec_SF = 0.0036*Td + 1.4801   
            case(30)    !Therminol 59: Reference: Therminol Reference Disk by Solutia: http://www.therminol.com/pages/tools/toolscd.asp
                spec_SF = 0.0033*Td + 1.6132
            case(31:35)	
            continue
            !!case(36:) !Any integer greater than 35
            !!!Call the user-defined property table
            !!lb=fl_bounds(fnum-35)
            !!ub=fl_bounds(fnum-35+1)-1
            !!if(ub.lt.lb) ub=size(fprop(1,:))
            !!dxx(:)=fprop(1,lb:ub)
            !!dyy(:)=fprop(2,lb:ub)
            !!call interp(Td,size(dxx),dxx,dyy,Gjsav,spec_SF)
            !!        if((Gjsav.eq.ub).or.(Gjsav.eq.lb)) dum=t_warn(Td,dxx(lb),dxx(ub),"User-specified fluid")
            ! case(36) !36-User defined SF HTF
            !     call NR_LINEAR_INTERPOLATION_00(Td,size(T_SF_HTF),T_SF_HTF,spec_SF_SF_HTF,HTFPropsav,spec_SF) !Specific heat, Td is in [C], cp in [kJ/kg/K]
            ! case(37) !37-User defined TES HTF
            !     call NR_LINEAR_INTERPOLATION_00(Td,size(T_TES_HTF),T_TES_HTF,spec_SF_TES_HTF,HTFPropsav,spec_SF) !Specific heat, Td is in [C], cp in [kJ/kg/K]
            case(40)
            !Specific Heat of Dowtherm A, kJ/kg/K
                spec_SF = 1.47524 + 0.00368606*(Td-273) - 0.00000516458*(Td-273)**2 + 8.99399E-09*(Td-273)**3 !Specific heat, Td is in [K], cp in [kJ/kg/K]
            end select

        end function

        !
        !*************** Dowtherm A **************************
        !
        !Enthalpy of Dowtherm A [J/kg]
        Double Precision Function H_Dowtherm_A(T) !T [K]
         ! chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.appliedthermalfluids.com/wp-content/uploads/2018/02/Dowtherm-A-heat-transfer-fluid-TDS.pdf
        implicit none
        Double Precision T, Td
        Td = T - 273.15 ! [C]
        !H_Dowtherm_A = (-19.8113 + 1.50647*T + 0.00144152*T**2) * 1000      !
        !H_Dowtherm_A = (-38.0792 + 1.50904*T + 0.00142671*T**2) * 1000
        H_Dowtherm_A = (-12.7078 + 1.481714*Td + 0.0014292857*Td**2) * 1000 ! [J/kg]      
        End Function


end module SF_props

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module SF_piping_functions
implicit none
    
    contains
        !***************************************************************************************************
        ! Friction factor (taken from Piping loss model)
        !***************************************************************************************************
        ! Uses an iterative method to solve the implicit friction factor function.
        ! For more on this method, refer to Fox, et al., 2006 Introduction to Fluid Mechanics.
        function FricFactor_piping(Rough, Reynold)result(f)

        implicit none

        double precision,intent(in)::Rough,Reynold
        double precision:: Test, TestOld, X, Xold, Slope, f
        double precision,parameter:: Acc = .01 !0.0001
        integer:: NumTries

        X = 33.33333  !1. / 0.03
        TestOld = X + 2. * Log10(Rough / 3.7 + 2.51 * X / Reynold)
        Xold = X
        X = 28.5714  !1. / (0.03 + 0.005)
        NumTries = 0

        do
            NumTries = NumTries + 1
            Test = X + 2 * Log10(Rough / 3.7 + 2.51 * X / Reynold)
            If (Abs(Test - TestOld) <= Acc) Then
                f = 1. / (X * X)
                exit 
            End If

            If (NumTries > 20) Then
                !call Messages(-1," Could not find friction factor solution",'Warning',0,250) 
                return
            End If

            Slope = (Test - TestOld) / (X - Xold)
            Xold = X
            TestOld = Test
            X = dmax1((Slope * X - Test) / Slope,1.e-5)
        enddo
        End Function
        !***************************************************************************************************
        ! Trough system piping loss model
        !***************************************************************************************************
        !
        ! This piping loss model is derived from the pressure drop calculations presented in the 
        ! following document:
        !
        !   Parabolic Trough Solar System Piping Model
        !   Final Report May 13, 2002 � December 31, 2004
        !
        !   B. Kelly
        !   Nexant, Inc. San Francisco, California
        !
        !   D. Kearney
        !   Kearney & Associates
        !   Vashon, Washington
        !
        !   Subcontract Report
        !   NREL/SR-550-40165
        !   July 2006
        !
        ! ----------------------------
        ! Note on use of this function
        ! ----------------------------
        ! The function returns the pressure drop across a given length of pipe, and also accounts for 
        ! a variety of possible pressure-loss components. This function should be called multiple times -
        ! once for each section under consideration.  For example, separate calls should be made for the
        ! HCE pressure drop, the pressure drop in each section of the header in which flow/geometrical 
        ! conditions vary, the section of pipe leading to the header, and so on.
        !
        ! ----------------------------
        ! Inputs
        ! ----------------------------
        ! No | Name         | Description                           | Units     |  Type
        !===================================================================================
        !  1 | Fluid        | Number associated with fluid type     | none      | dble
        !  2 | m_dot        | Mass flow rate of the fluid           | kg/s      | dble
        !  3 | T            | Fluid temperature                     | K         | dble
        !  4 | P            | Fluid pressure                        | Pa        | dble
        !  5 | D            | Diameter of the contact surface       | m         | dble
        !  6 | Rough        | Pipe roughness                        | m         | dble
        !  7 | L_pipe       | Length of pipe for pressure drop      | m         | dble
        !  8 | Nexp         | Number of expansions                  | none      | dble
        !  9 | Ncon         | Number of contractions                | none      | dble
        ! 10 | Nels         | Number of standard elbows             | none      | dble
        ! 11 | Nelm         | Number of medium elbows               | none      | dble
        ! 12 | Nell         | Number of long elbows                 | none      | dble
        ! 13 | Ngav         | Number of gate valves                 | none      | dble
        ! 14 | Nglv         | Number of globe valves                | none      | dble
        ! 15 | Nchv         | Number of check valves                | none      | dble
        ! 16 | Nlw          | Number of loop weldolets              | none      | dble
        ! 17 | Nlcv         | Number of loop control valves         | none      | dble
        ! 18 | Nbja         | Number of ball joint assemblies       | none      | dble
        !===================================================================================
        ! ----------------------------
        ! Outputs
        ! ----------------------------
        ! 1. PressureDrop  (Pa)

        double precision function PressureDrop(Fluid,m_dot,T,P,D,Rough,L_pipe,&
                                            Nexp,Ncon,Nels,Nelm,Nell,Ngav,Nglv,Nchv,Nlw,Nlcv,Nbja)

        use SF_props

        implicit none                                       
                                            
        real(8),intent(in):: Fluid,m_dot,T,P,D,Rough,L_pipe,Nexp,Ncon,Nels,Nelm,Nell,Ngav,Nglv,Nchv,&
                            Nlw,Nlcv,Nbja
        real(8):: rho, v_dot, mu, nu, u_fluid, Re, f, DP_pipe, DP_exp,DP_con,DP_els,DP_elm,DP_ell,DP_gav,&
                DP_glv,DP_chv,DP_lw,DP_lcv,DP_bja, FricFactor,&
                HL_pm,pi,g
        pi=3.1415928; g=9.80665

        !Calculate fluid properties and characteristics
        rho = density_SF(fluid,T,P)
        mu = viscosity_SF(fluid,T,P)
        nu = mu/rho
        v_dot = m_dot/rho   !fluid volumetric flow rate
        u_fluid = v_dot/(pi*(D/2.)*(D/2.))  !Fluid mean velocity

        !Dimensionless numbers
        Re = u_fluid*D/nu
        if(Re<2300.) then
            f = 64./dmax1(Re,1.d0)
        else
            f = FricFactor_piping(Rough/D,Re)
        endif

        !Calculation of pressure loss from pipe length
        HL_pm = f*u_fluid*u_fluid/(2.*D*g)
        DP_pipe = HL_pm*rho*g*L_pipe

        !Calculation of pressure loss from Fittings
        DP_exp = 0.25*rho*u_fluid*u_fluid*Nexp
        DP_con = 0.25*rho*u_fluid*u_fluid*Ncon
        DP_els = 0.9 * D / f * HL_pm * rho * g * Nels
        DP_elm = 0.75 * D / f * HL_pm * rho * g * Nelm
        DP_ell = 0.6 * D / f * HL_pm * rho * g * Nell
        DP_gav = 0.19 * D / f * HL_pm * rho * g * Ngav
        DP_glv = 10.0 * D / f * HL_pm * rho * g * Nglv
        DP_chv = 2.5 * D / f * HL_pm * rho * g * Nchv
        DP_lw = 1.8 * D / f * HL_pm * rho * g * Nlw
        DP_lcv = 10.0 * D / f * HL_pm * rho * g * Nlcv
        DP_bja = 8.69 * D / f * HL_pm * rho * g * Nbja

        PressureDrop = sum((/DP_pipe, DP_exp,DP_con,DP_els,DP_elm,DP_ell,DP_gav,&
                        DP_glv,DP_chv,DP_lw,DP_lcv,DP_bja/))

        end 


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! C_v Computation for various valves
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        double precision function CV_data(Valve_Type, D_in, Valve_position)result(C_v)
            !!!!!!!!!!!!!!!!!!
            ! INPUTS
            !! Valve_Type: Select Valve Type [Integer]
            !! D_in: Valve Diameter [m]
            !! Valve_position: Position of Valve [0-1] with 1 being fully open
            !!!!!!!!!!!!!!!!!!
            ! Valve Types
            !! 1: Concentric Butterfly
            !! 2: Triple Offset Butterfly
            

            implicit none    
            
            INTEGER, intent(in) :: Valve_Type
            DOUBLE PRECISION, intent(in) ::  Valve_position, D_in

            Double precision, dimension(3,2) :: test_mat
            Double precision, dimension(7, 10) :: CV_concentric, CV_triple
            Double Precision, Dimension(7) :: D_concentric, D_triple
            Double Precision, Dimension(10) :: Pos_concentric, Pos_triple
            Integer :: Found, N, Match_D, Match_Pos, Ind_D, Ind_Pos, Ind_D_low, Ind_D_high, Ind_Pos_low, Ind_Pos_High
            Double Precision :: x0, x1, y0, y1, yL0, yL1, yH0, yH1, D, C_V_min
            
            D = D_in*39.3701 ! [m -> inches]
            C_v_min = 0.0001d0 !minimum CV allowed


            if(Valve_position == 0.d0)then
                C_v = C_V_min ! Subject to change, need a small amount of flow so hydraulic model won't crash
                return
            endif
            
            test_mat = Transpose(reshape((/1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0/), (/2,3/)))

            select case(Valve_Type)

            case(1) ! Concentric Butterfly Valve: https://www.valvesonline.com.au/references/flow-rates/butterfly-valves/
                D_concentric = (/ 4.d0, 8.d0, 12.d0, 16.d0, 20.d0, 24.d0, 32.d0 /)
                Pos_concentric = (/ 0.d0, 1.d0/9.d0, 2.d0/9.d0, 3.d0/9.d0, 4.d0/9.d0, 5.d0/9.d0, 6.d0/9.d0, 7.d0/9.d0, 8.d0/9.d0, 1.d0 /)
                CV_concentric = TRANSPOSE(RESHAPE(&
                                (/&
                                0.d0, 0.5d0, 17.d0, 36.d0, 78.d0, 139.d0, 230.d0, 364.d0, 546.d0, 600.d0, &
                                0.d0, 3.d0, 89.d0, 188.d0, 408.d0, 727.d0, 1202.d0, 1903.d0, 2854.d0, 3136.d0,&
                                0.d0, 5.d0, 234.d0, 495.d0, 1072.d0, 1911.d0, 3162.d0, 5005.d0, 7507.d0, 8250.d0,&
                                0.d0, 8.d0, 464.d0, 983.d0, 2130.d0, 3797.d0, 6282.d0, 9942.d0, 14913.d0, 16388.d0,&
                                0.d0, 14.d0, 791.d0, 1674.d0, 3628.d0, 6465.d0, 10698.d0, 16931.d0, 25396.d0, 27908.d0,&
                                0.d0, 22.d0, 1222.d0, 2587.d0, 5605.d0, 9989.d0, 16528.d0, 26157.d0, 39236.d0, 43116.d0,&
                                0.d0, 45.d0, 2387.d0, 4791.d0, 8736.d0, 13788.d0, 20613.d0, 31395.d0, 48117.d0, 68250.d0&
                                /)&
                                , (/10, 7/) )  )

                ! Interpolate to get C_v value
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                ! Determine indices for Diameter interpolation
                found = 0
                match_D = 0
                do n = 1, Size(D_concentric)
                    if(ABS(D-D_concentric(n)) < 0.05d0)then
                        match_D = 1
                        ind_D = n
                    elseif(D_concentric(n)>D .and. found == 0)then
                        ind_D_low = n-1
                        ind_D_high = N
                        found = 1
                    endif
                end do

                ! Determine indices for Position interpolation
                found = 0 
                match_Pos = 0
                do n = 1, Size(Pos_concentric)
                    if(Pos_concentric(n) == Valve_Position)then
                        match_Pos = 1
                        ind_Pos = n
                    elseif(Pos_concentric(n)>Valve_Position .and. found == 0)then
                        ind_Pos_low = n-1
                        ind_Pos_high = N
                        found = 1
                    endif
                end do

                ! Check if the indices match the lookup table
                if(match_Pos == 1 .and. match_D == 1)then
                    C_v = CV_concentric(ind_D, ind_Pos)
                    C_v = max(C_v, C_v_min)
                ! Check if the diameter is the same as the lookup table
                elseif(match_D == 1)then
                    x0 = Pos_concentric(ind_Pos_low)
                    x1 = Pos_concentric(ind_Pos_high)
                    y0 = CV_concentric(ind_D, ind_Pos_low)
                    y1 = CV_concentric(ind_D, ind_Pos_high)
                    C_v = y0 + (Valve_Position-x0)*(y1-y0)/(x1-x0)
                    C_v = max(C_v, C_v_min)

                ! Check if the position is the same as the lookup table
                elseif(match_Pos == 1)then
                    x0 = D_concentric(ind_D_low)
                    x1 = D_concentric(ind_D_high)
                    y0 = CV_concentric(ind_D_low, ind_Pos)
                    y1 = CV_concentric(ind_D_high, ind_Pos)
                    C_v = y0 + (D-x0)*(y1-y0)/(x1-x0)
                    C_v = max(C_v, C_v_min)
                ! else interpolate between diameter and position
                else
                    ! First interpolate position for each diameter
                    x0 = Pos_concentric(ind_Pos_low)
                    x1 = Pos_concentric(ind_Pos_high)
                    yL0 = CV_concentric(ind_D_low, ind_Pos_low)
                    yL1 = CV_concentric(ind_D_low, ind_Pos_high)
                    yH0 = CV_concentric(ind_D_high, ind_Pos_low)
                    yH1 = CV_concentric(ind_D_high, ind_Pos_high)

                    y0 = yL0 + (Valve_Position-x0)*(yL1-yL0)/(x1-x0)
                    y1 = yH0 + (Valve_Position-x0)*(yH1-YH0)/(x1-x0)

                    ! Lastly interpolate between diameters
                    x0 = D_concentric(ind_D_low)
                    x1 = D_concentric(ind_D_high)
                    C_v = y0 + (D-x0)*(y1 - y0)/(x1-x0)
                    C_v = max(C_v, C_v_min)

                endif
                

            case(2) ! Triple Offset Butterfly: https://www.valvesonline.com.au/references/flow-rates/butterfly-valves/
                D_triple = (/ 4.d0, 8.d0, 12.d0, 16.d0, 20.d0, 24.d0, 32.d0 /)
                Pos_triple = (/ 0.d0, 1.d0/9.d0, 2.d0/9.d0, 3.d0/9.d0, 4.d0/9.d0, 5.d0/9.d0, 6.d0/9.d0, 7.d0/9.d0, 8.d0/9.d0, 1.d0 /)
                CV_triple = TRANSPOSE(RESHAPE(&
                                (/&
                                0.d0, 8.4d0, 29.3d0, 58.6d0, 92.d0, 140.d0, 200.d0, 330.d0, 370.d0, 420.d0, &
                                0.d0, 38.2d0, 140.d0, 270.d0, 420.d0, 640.d0, 900.d0, 1500.d0, 1690.d0, 1920.d0,&
                                0.d0, 88.4d0, 310.d0, 620.d0, 980.d0, 1460.d0, 2080.d0, 3450.d0, 3890.d0, 4420.d0,&
                                0.d0, 150.d0, 530.d0, 1060.d0, 1660.d0, 2490.d0, 3540.d0, 5870.d0, 6620.d0, 7520.d0,&
                                0.d0, 270.d0, 930.d0, 1850.d0, 2900.d0, 4350.d0, 6190.d0, 10300.d0, 11600.d0, 13200.d0,&
                                0.d0, 420.d0, 1450.d0, 2890.d0, 4530.d0, 6800.d0, 9680.d0, 16100.d0, 18100.d0, 20600.d0,&
                                0.d0, 810.d0, 2820.d0, 5630.d0, 8840.d0, 13300.d0, 18900.d0, 31400.d0, 35400.d0, 40200.d0&
                                /)&
                                , (/10, 7/) ) )
                
                ! Interpolate to get C_v value
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                ! Determine indices for Diameter interpolation
                found = 0
                match_D = 0
                do n = 1, Size(D_triple)
                    if(ABS(D-D_triple(n)) < 0.05d0)then
                        match_D = 1
                        ind_D = n
                    elseif(D_triple(n)>D .and. found == 0)then
                        ind_D_low = n-1
                        ind_D_high = N
                        found = 1
                    endif
                end do

                ! Determine indices for Position interpolation
                found = 0 
                match_Pos = 0
                do n = 1, Size(Pos_triple)
                    if(Pos_triple(n) == Valve_Position)then
                        match_Pos = 1
                        ind_Pos = n
                    elseif(Pos_triple(n)>Valve_Position .and. found == 0)then
                        ind_Pos_low = n-1
                        ind_Pos_high = N
                        found = 1
                    endif
                end do

                ! Check if the indices match the lookup table
                if(match_Pos == 1 .and. match_D == 1)then
                    C_v = CV_triple(ind_D, ind_Pos)
                    C_v = max(C_v, C_v_min)

                ! Check if the diameter is the same as the lookup table
                elseif(match_D == 1)then
                    x0 = Pos_triple(ind_Pos_low)
                    x1 = Pos_triple(ind_Pos_high)
                    y0 = CV_triple(ind_D, ind_Pos_low)
                    y1 = CV_triple(ind_D, ind_Pos_high)
                    C_v = y0 + (Valve_Position-x0)*(y1-y0)/(x1-x0)
                    C_v = max(C_v, C_v_min)

                ! Check if the position is the same as the lookup table
                elseif(match_Pos == 1)then
                    x0 = D_triple(ind_D_low)
                    x1 = D_triple(ind_D_high)
                    y0 = CV_triple(ind_D_low, ind_Pos)
                    y1 = CV_triple(ind_D_high, ind_Pos)
                    C_v = y0 + (D-x0)*(y1-y0)/(x1-x0)
                    C_v = max(C_v, C_v_min)

                ! else interpolate between diameter and position
                else
                    ! First interpolate position for each diameter
                    x0 = Pos_triple(ind_Pos_low)
                    x1 = Pos_triple(ind_Pos_high)
                    yL0 = CV_triple(ind_D_low, ind_Pos_low)
                    yL1 = CV_triple(ind_D_low, ind_Pos_high)
                    yH0 = CV_triple(ind_D_high, ind_Pos_low)
                    yH1 = CV_triple(ind_D_high, ind_Pos_high)

                    y0 = yL0 + (Valve_Position-x0)*(yL1-yL0)/(x1-x0)
                    y1 = yH0 + (Valve_Position-x0)*(yH1-YH0)/(x1-x0)

                    ! Lastly interpolate between diameters
                    x0 = D_triple(ind_D_low)
                    x1 = D_triple(ind_D_high)
                    C_v = y0 + (D-x0)*(y1 - y0)/(x1-x0)
                    C_v = max(C_v, C_v_min)

                endif 
            end select


        end function CV_DATA

        DOUBLE PRECISION FUNCTION ControlValve_deltaP(Fluid, Valve_Type, m_dot, D_in, Valve_Position, T)result(dP)
            !!!!!!!!!!!!!!!!!!
            ! INPUTS
            !! Fluid: Fluid label for property look up [double precision]
            !! Valve_Type: Select Valve Type [Integer]
            !! m_dot: Mass flow of fluid [kg/s]
            !! D_in: Valve Diameter [m]
            !! Valve_position: Position of Valve [0-1] with 1 being fully open
            !! T: Temperature of fluid [K]
            !!!!!!!!!!!!!!!!!!
            ! OUTPUTS
            !! dP: Pressure drop across valve [Pa]
        
            use SF_props
        
            implicit none    
            
            INTEGER, intent(in) :: Valve_Type
            DOUBLE PRECISION, intent(in) ::  Valve_position, D_in, Fluid, m_dot, T

            DOUBLE PRECISION :: Q, SG, Cv, rho_fluid

            ! Compute Volumetric Flow Rate
            rho_fluid = Density_SF(Fluid, T, 0.d0)
            Q = m_dot/rho_fluid
            ! Convert flowrate to gpm
            Q = Q*15850.323140625002 
            ! Compute specific gravity of fluid 
            SG = rho_fluid/1000.d0
            ! Compute Cv of valve
            Cv = CV_data(Valve_Type, D_in, Valve_position)
            ! Compute pressure drop [Pa]
            dP = SG*Q**2/(Cv**2) * 6894.76
        end function ControlValve_deltaP

        Function pipe_dTdt(n_nodes, T, Vol, mass_flow, mc_mult, Fluid_ID, heat_loss, L_cv)result(dt)
            !!!!!!!!!!!!!!!!!!!!!!!
            ! INPUTS
            !! n_nodes: number of nodes [integer]
            !! T: array of temperatures corresponding to each node [K]
            !! Vol_norm: The volume of a typical control volume [m^3]
            !! mass_flow: Mass flow through pipe [kg/s]
            !! mc_mult: thermal capacitance multiplier to account for piping and supports [-]
            !! Fluid_ID: Fluid identifer for function calls
            !!!!!!!!!!!!!!!!!!!!!!
            ! OUTPUTS
            !! dtdt: array of temperature derivatives corresponding to each node
            use SF_props

            implicit none

            integer, intent(in) :: n_nodes
            Double precision, dimension(n_nodes), intent(in) :: T
            Double precision, intent(in) :: Vol, mass_flow, Fluid_ID, mc_mult, heat_loss, L_cv
            Double precision, dimension(n_nodes) :: dt
            Double precision, dimension(n_nodes-1) :: dtdt_bar
            Double precision :: T_ave, rho, c, q_out
            Integer :: n


            ! Loop Through Control Volumes to compute CV temperature rate
            do n = 1,n_nodes-1
                ! Heat loss
                q_out = heat_loss*L_cv
                ! Compute CV average temp and properties
                T_ave = (T(n)+T(n+1))/2.d0
                rho = Density_SF(Fluid_ID, T_ave, 0.d0)
                c = spec_SF(Fluid_ID, T_ave, 0.d0)*1000.d0 ! [J/kg-K]
                ! Compute CV temperature rate
                dtdt_bar(n) = 1/Vol/rho/c/mc_mult*(mass_flow*c*(T(n)-T(n+1)) - q_out)
            end do

            ! Compute Nodal Temperature Rates
            dt(1) = 0.0d0 ! First node is set by input to type
            dt(n_nodes) = dtdt_bar(n_nodes-1) ! Cannot average last node
            do n = 2,n_nodes-1 ! Remaining middle nodes
                dt(n) = (dtdt_bar(n-1) + dtdt_bar(n))/2.d0
            end do
            
        end function pipe_dTdt

    
end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module Header_functions
implicit none

    contains

        function diams_inlet(n_cv, n_loop, L_row, L_exp, diam_file)result(Diam_cv)
            integer :: n, cc, loop_count
            integer, intent(in) :: n_cv, n_loop
            double precision, intent(in) :: L_row, L_exp
            double precision :: pi, holder, D_curr
            double precision, dimension(n_cv) :: Diam_cv
            double precision, dimension(n_loop/2) :: Diam
            character(len=100) :: diam_file

            ! Get specified geometries
            OPEN(2, File = diam_file)
            do cc = 1,n_loop/2
                READ(2,*), Diam(cc), holder
            end do
            CLOSE(2)
            
            cc = 1
            loop_count = 1
            ! Loop through all control volumes
            do n = 1,n_cv
                D_curr = Diam(loop_count)
                ! Check if control volume is part of an expansion loop
                if(cc > 1)then
                    Diam_cv(n) = D_curr
                    cc = cc + 1
                    if(cc == 4)then
                        cc = 1
                        loop_count = loop_count + 1
                    end if
                ! Else control volume is the segment of header between sf loop inlets
                else
                    Diam_cv(n) = D_curr
                    cc = cc +1
                    loop_count = loop_count + 1
                end if
            end do

        end function diams_inlet

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! This function creates an array with the volumes of each control volume in an inlet header
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function vols_inlet(n_cv, n_loop, L_row, L_exp, diam_file)result(vols)
            integer :: n, cc, loop_count
            integer, intent(in) :: n_cv, n_loop
            double precision, intent(in) :: L_row, L_exp
            double precision :: D_curr, pi, holder, L_tot
            double precision, dimension(n_cv) :: vols
            double precision, dimension(n_loop/2) :: Diam
            character(len=60), intent(in) :: diam_file

            ! Get specified geometries
            OPEN(2, File = diam_file)
            do cc = 1,n_loop/2
                READ(2,*), Diam(cc), holder
            end do
            CLOSE(2)

            pi = 3.141529
            cc = 1
            loop_count = 1
            L_tot = 0.d0
            ! Loop through all control volumes
            do n = 1,n_cv
                D_curr = Diam(loop_count)
                ! Check if control volume is part of an expansion loop
                if(cc > 1)then
                    vols(n) = (D_curr/2)**2*pi*(L_row*2 + L_exp*2)/2
                    L_tot = L_tot + (L_row*2 + L_exp*2)/2
                    cc = cc + 1
                    if(cc == 4)then
                        cc = 1
                        loop_count = loop_count + 1
                    end if
                ! Else control volume is the pipe between two sf loop returns
                else
                    vols(n) = (D_curr/2)**2*pi*L_row*2
                    cc = cc +1
                    loop_count = loop_count + 1
                    L_tot = L_tot + L_row*2
                end if
            end do
            holder = 1
            
        end function vols_inlet

        function diams_return(n_cv, n_loop, L_row, L_exp, diam_file)result(Diam_cv)
            integer :: n, cc, loop_count
            integer, intent(in) :: n_cv, n_loop
            double precision, intent(in) :: L_row, L_exp
            double precision :: pi, holder, D_curr
            double precision, dimension(n_cv) :: Diam_cv
            double precision, dimension(n_loop/2) :: Diam
            character(len=60) :: diam_file

            ! Get specified geometries
            OPEN(2, File = diam_file)
            do cc = 1,n_loop/2
                READ(2,*), holder, Diam(cc)
            end do
            CLOSE(2)
            
            cc = 1
            loop_count = 1
            ! Loop through all control volumes
            do n = 1,n_cv
                D_curr = Diam(loop_count)
                ! Check if control volume is part of an expansion loop
                if(cc > 1)then
                    Diam_cv(n) = D_curr
                    cc = cc + 1
                    if(cc == 4)then
                        cc = 1
                        loop_count = loop_count + 1
                    end if
                ! Else control volume is the segment of header between sf loop inlets
                else
                    Diam_cv(n) = D_curr
                    cc = cc +1
                    loop_count = loop_count + 1
                end if
            end do
            holder = 1
        end function diams_return


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! This function creates an array with the volumes of each control volume in a return header
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function vols_return(n_cv, n_loop, L_row, L_exp, diam_file)result(vols)
            integer :: n, cc, loop_count
            integer, intent(in) :: n_cv, n_loop
            double precision, intent(in) :: L_row, L_exp
            double precision :: D_curr, pi, holder
            double precision, dimension(n_cv) :: vols
            double precision, dimension(n_loop/2) :: Diam
            character(len=60) :: diam_file

            ! Get specified geometries
            OPEN(2, File = diam_file)
            do cc = 1,n_loop/2
                READ(2,*), holder, Diam(cc)
            end do
            CLOSE(2)


            pi = 3.141529
            cc = 1
            loop_count = 1
            do n = 1,n_cv
                D_curr = Diam(loop_count)
                
                if(cc > 1)then
                    vols(n) = (D_curr/2)**2*pi*(L_row*2 + L_exp*2)/2
                    cc = cc + 1
                    if(cc == 4)then
                        cc = 1
                        loop_count = loop_count + 1
                    end if
                else
                    vols(n) = (D_curr/2)**2*pi*L_row*2
                    cc = cc +1
                    loop_count = loop_count+1
                end if
                    
                    
            end do
        end function vols_return

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! This function computes the nodal temperature rate of change for the inlet header
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function dT_dt_inlet(m_dots, t, Vols, L_cv, mc, t_bar, n_nodes, fluid, heat_loss)result(dT)
            integer, intent(in) :: n_nodes
            double precision, dimension(n_nodes-1), intent(in) :: m_dots, t_bar, Vols, L_cv
            double precision, dimension(n_nodes), intent(in) :: T
            double precision, intent(in) :: fluid, mc, heat_loss
            integer :: n
            double precision, dimension(n_nodes-1) :: dT_bar
            double precision, dimension(n_nodes) :: dT
            double precision :: c, rho, h1, h2


            ! Compute control volume temperature rate of change
            do n = 1, (n_nodes-1)
                c = 1000.d0*spec_NN(fluid, t_bar(n), 0.d0)
                rho = den_NN(fluid, t_bar(n), 0.d0)
                h1 = H_Dowtherm_A(t(n))
                h2 = H_Dowtherm_A(t(n+1))
                ! dT_bar(n) = 1/mc/(Vols(n)*rho*c)*(m_dots(n)*(h1 - h2) - 3500.d0*L_cv(n))
                dT_bar(n) = 1/mc/(Vols(n)*rho*c)*(m_dots(n)*(h1 - h2) - heat_loss*L_cv(n))
            end do

            ! Compute nodal temperature rate of change
            dT(1) = 0.d0
            do n = 2,(n_nodes-1)
                dT(n) = 1.0/2.0*(dT_bar(n-1) + dT_bar(n))
            end do
            dT(n_nodes) = dT_bar(n_nodes-1)

        end function dT_dt_inlet


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! This function computes the nodal temperature rate of change for the return header
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function dT_dt_return(m_dot_l, m_dot_r, m_dot_out, T_l, T_r, T, T_bar, Vols, L_cv, mc, n_nodes, n_loop, inds_header, fluid, heat_loss)result(dT)
            integer, intent(in) :: n_nodes, n_loop
            double precision, dimension(n_nodes-1), intent(in) :: m_dot_out, T_bar, Vols, L_cv
            double precision, dimension(n_nodes), intent(in) :: T
            double precision, dimension(n_loop/2), intent(in) :: T_l, T_r, m_dot_l, m_dot_r
            double precision, intent(in) :: fluid, mc, heat_loss
            integer :: n, jj, cc
            double precision, dimension(n_nodes-1) :: dT_bar
            double precision, dimension(n_nodes) :: dT
            integer, dimension(n_loop/2) :: inds_header
            double precision :: rho, h1, h2, c, T_in, h_r, h_l

            ! Compute control volume temperature rate of changes
            jj = 2
            do n = 1, n_nodes-1
                rho = den_NN(fluid, T_bar(n), 0.d0)
                c = 1000.d0*spec_NN(fluid, T_bar(n), 0.d0)
                h1 = H_Dowtherm_A(t(n))
                h2 = H_Dowtherm_A(t(n+1))
                if(n==1)then
                    dT_bar(n) = 1/mc/(Vols(n)*rho*c)*(m_dot_out(n)*(h1 - h2))
                else
                    if(n == inds_header(jj))then
                        h_r = H_Dowtherm_A(T_r(jj))
                        h_l = H_Dowtherm_A(T_l(jj))
                        dT_bar(n) = 1/mc/(Vols(n)*rho*c) * ( m_dot_l(jj)*h_l + m_dot_r(jj)*h_r + m_dot_out(n-1)*h1 - m_dot_out(n)*h2 - heat_loss*L_cv(n))
                        jj = jj + 1
                    else
                        dT_bar(n) = 1/Vols(n)/rho/c * ( m_dot_out(n)*(h1 - h2) - heat_loss*L_cv(n))
                    end if
                end if
            end do

            ! Compute nodal temperature rate of changes
            dT(1) = 0.d0
            do n = 2, n_nodes-1
                dT(n) = 1.0/2.0*(dT_bar(n-1) + dT_bar(n))
            end do
            dT(n_nodes) = dT_bar(n_nodes-1)
        end function dT_dt_return




        !  DENSITY AND SPECIFIC HEAT FUNCTIONS TAKEN FROM SAM
        ! ||||||||||||||||||||||||||||||||||||||||||||||||||||
        ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        double precision function den_NN(fnumd,T,P)
            !This function accepts as inputs temperature [K] and pressure [Pa]
            !This function outputs in units of [kg/m^3]
            double precision::xlo,xhi, Dens_fluid, Td,HTFPropsav
            double precision::T,P,fnumd
            !!double precision,dimension(size(fprop(1,:)))::dxx,dyy !Create dummy arrays
            integer::fnum,lb,ub,dum,t_warn
            !den_NN=1.
            fnum=nint(fnumd)
            Td=T-273.15             !Convert from K to C
                
            select case(fnum)
            case(1)   !    1.) Air
            den_NN = P/(287.*T)
            case(2)   !    2.) Stainless_AISI316
                den_NN=8349.38 - 0.341708*T - 0.0000865128*T*T  !EES
            case(3)   !    3.) Water (liquid)
                den_NN = 1000 
            case(4)   !    4.) Steam
                continue
            case(5)   !    5.) CO2
                continue
            case(6)   !    6.) Salt (68% KCl, 32% MgCl2)
            den_NN = 1E-10*T*T*T - 3E-07*T*T - 0.4739*T + 2384.2
            case(7)   !    7.) Salt (8% NaF, 92% NaBF4)
            den_NN = 8E-09*T*T*T - 2E-05*T*T - 0.6867*T + 2438.5
            case(8)   !    8.) Salt (25% KF, 75% KBF4)
            den_NN = 2E-08*T*T*T - 6E-05*T*T - 0.7701*T + 2466.1
            case(9)   !    9.) Salt (31% RbF, 69% RbBF4)
            den_NN = -1E-08*T*T*T + 4E-05*T*T - 1.0836*T + 3242.6
            case(10)   !    10.) Salt (46.5% LiF, 11.5%NaF, 42%KF)
            den_NN =  -2E-09*T*T*T + 1E-05*T*T - 0.7427*T + 2734.7
            case(11)   !    11.) Salt (49% LiF, 29% NaF, 29% ZrF4)
            den_NN = -2E-11*T*T*T + 1E-07*T*T - 0.5172*T + 3674.3
            case(12)   !    12.) Salt (58% KF, 42% ZrF4)
            den_NN =  -6E-10*T*T*T + 4E-06*T*T - 0.8931*T + 3661.3
            case(13)   !    13.) Salt (58% LiCl, 42% RbCl)
            den_NN = -8E-10*T*T*T + 1E-06*T*T - 0.689*T + 2929.5
            case(14)   !    14.) Salt (58% NaCl, 42% MgCl2)
            den_NN = -5E-09*T*T*T + 2E-05*T*T - 0.5298*T + 2444.1
            case(15)   !    15.) Salt (59.5% LiCl, 40.5% KCl)
            den_NN = 1E-09*T*T*T - 5E-06*T*T - 0.864*T + 2112.6
            case(16)   !    16.) Salt (59.5% NaF, 40.5% ZrF4)
            den_NN =  -5E-09*T*T*T + 2E-05*T*T - 0.9144*T + 3837.
            case(17)   !    17.) Salt (60% NaNO3, 40% KNO3)
            den_NN = dmax1(-1E-07*T*T*T + 0.0002*T*T - 0.7875*T + 2299.4,1000.d0)
            case(18)
            !den_NN of Nitrate Salt, [kg/m3]
            den_NN = dmax1(2090 - 0.636 * (T-273.15),1000.d0)
            case(19)
            !den_NN of Caloria HT 43 [kg/m3]
            den_NN = dmax1(885 - 0.6617 * Td - 0.0001265 * Td*Td,100.d0)
            case(20)
            !den_NN of HITEC XL Nitrate Salt, [kg/m3]
            den_NN = dmax1(2240 - 0.8266 * Td,800.d0)
            case(21)
            !den_NN of Therminol Oil [kg/m3]
            den_NN = dmax1(1074.0 - 0.6367 * Td - 0.0007762 * Td*Td,400.d0)
            case(22)
            !den_NN of HITEC Salt, [kg/m3]
            den_NN = dmax1(2080 - 0.733 * Td,1000.d0)
            case(23)
            !den_NN of Dowtherm Q [kg/m3]
            den_NN = dmax1(-0.757332 * Td + 980.787,100.d0)                               ! Russ 10-2-03
            case(24)
            !den_NN of Dowtherm RP [kg/m3]
            den_NN = dmax1(-0.000186495 * Td*Td - 0.668337 * Td + 1042.11,200.d0)  !Russ 10-2-03
            case(25)
            !den_NN of HITEC XL Nitrate Salt, [kg/m^3]
            den_NN = dmax1(2240 - 0.8266 * Td,800.d0)
            case(26) !Argon
            den_NN = dmax1(P/(208.13*T),1.e-10)
            case(27) !Hydrogen
            den_NN = dmax1(P/(4124.*T),1.e-10)
            case(28)    !T-91 Steel: "Thermo hydraulic optimisation of the EURISOL DS target" - Paul Scherrer Institut
            den_NN = -0.3289*Td + 7742.5
            case(29)    !Therminol 66: Reference: Therminol Reference Disk by Solutia: http://www.therminol.com/pages/tools/toolscd.asp
            den_NN = -0.7146*Td + 1024.8
            case(30)    !Therminol 59: Reference: Therminol Reference Disk by Solutia: http://www.therminol.com/pages/tools/toolscd.asp
            den_NN = -0.0003*Td*Td - 0.6963*Td + 988.44
            case(31:35) 
            continue !no informaion
            !!case(36:) !Any integer greater than 35
            !!!Call the user-defined property table
            !!lb=fl_bounds(fnum-35)
            !!ub=fl_bounds(fnum-35+1)-1
            !!if(ub.lt.lb) ub=size(fprop(1,:))
            !!dxx(:)=fprop(1,lb:ub)
            !!dyy(:)=fprop(3,lb:ub)
            !!call interp(Td,size(dxx),dxx,dyy,Gjsav,den_NN)
            !!if((Gjsav.eq.ub).or.(Gjsav.eq.lb)) dum=t_warn(Td,dxx(lb),dxx(ub),"User-specified fluid")
            ! case(36) !36-User defined SF HTF
            ! call NR_LINEAR_INTERPOLATION_00(Td,size(T_SF_HTF),T_SF_HTF,den_NN_SF_HTF,HTFPropsav,den_NN) !den_NN, Td is in [C], den_NN in [kg/m3]
            ! case(37) !37-User defined TES HTF
            ! call NR_LINEAR_INTERPOLATION_00(Td,size(T_TES_HTF),T_TES_HTF,den_NN_TES_HTF,HTFPropsav,den_NN) !den_NN, Td is in [C], den_NN in [kg/m3]
            case(40)
            !den_NN of Dowtherm A [kg/m3]
            den_NN = dmax1(1063.61 - 0.605235*Td - 0.000860877*Td*Td,400.d0)!den_NN, Td is in [C], den_NN in [kg/m3]
            end select

        end function

        double precision function spec_NN(fnumd,T,P)
            !This function accepts as inputs temperature [K] and pressure [Pa]
            !This function outputs in units of [kJ/kg-K]
            double precision::xlo,xhi, Td,HTFPropsav
            double precision,intent(in)::T,P,fnumd
            !!double precision,dimension(size(fprop(1,:)))::dxx,dyy !Create dummy arrays
            integer::fnum,lb,ub,dum,t_warn

            spec_NN=1.
            fnum=nint(fnumd)
            Td = T - 273.15
            select case(fnum)
            case(1)   !    1.) Air
                spec_NN = 1.03749 - 0.000305497*T + 7.49335E-07*T*T - 3.39363E-10*T*T*T
            !spec_NN = 1.03749 - 0.000305497*T + 7.49335E-07*T*T - 3.39363E-10*T*T*T
            case(2)   !    2.) Stainless_AISI316
                spec_NN = 0.368455 + 0.000399548*T - 1.70558E-07*T*T !EES
            case(3)   !    3.) Water (liquid)
                spec_NN = 4.181d0  !mjw 8.1.11 
            case(4)   !    4.) Steam
                continue
            case(5)   !    5.) CO2
                continue
            case(6)   !    6.) Salt (68% KCl, 32% MgCl2)
                spec_NN = 1.156
            case(7)   !    7.) Salt (8% NaF, 92% NaBF4)
                spec_NN = 1.507
            case(8)   !    8.) Salt (25% KF, 75% KBF4)
                spec_NN = 1.306
            case(9)   !    9.) Salt (31% RbF, 69% RbBF4)
                spec_NN = 9.127
            case(10)   !    10.) Salt (46.5% LiF, 11.5%NaF, 42%KF)
                spec_NN = 2.010
            case(11)   !    11.) Salt (49% LiF, 29% NaF, 29% ZrF4)
                spec_NN = 1.239
            case(12)   !    12.) Salt (58% KF, 42% ZrF4)
                spec_NN = 1.051
            case(13)   !    13.) Salt (58% LiCl, 42% RbCl)
                spec_NN = 8.918
            case(14)   !    14.) Salt (58% NaCl, 42% MgCl2)
                spec_NN = 1.080
            case(15)   !    15.) Salt (59.5% LiCl, 40.5% KCl)
                spec_NN = 1.202
            case(16)   !    16.) Salt (59.5% NaF, 40.5% ZrF4)
                spec_NN = 1.172
            case(17)   !    17.) Salt (60% NaNO3, 40% KNO3)
                spec_NN = -1E-10*T*T*T + 2E-07*T*T + 5E-06*T + 1.4387
            case(18) !Heat Capacity of Nitrate Salt, [J/kg/K]
                spec_NN = (1443. + 0.172 * (T-273.15))/1000.d0
            case(19)
            !Specific Heat of Caloria HT 43 [J/kgC]
                spec_NN = (3.88 * (T-273.15) + 1606.0)/1000.
            case(20)
            !Heat Capacity of HITEC XL Nitrate Salt, [J/kg/K]
                spec_NN = dmax1(1536 - 0.2624 * Td - 0.0001139 * Td * Td,1000.d0)/1000.
            case(21)
            !Specific Heat of Therminol Oil, kJ/kg/K
                spec_NN = (1.509 + 0.002496 * Td + 0.0000007888 * Td*Td)
            case(22)
            !Heat Capacity of HITEC Salt, [J/kg/K]
                spec_NN = (1560 - 0.0 * Td)/1000.
            case(23)
            !Specific Heat of Dowtherm Q, J/kg/K
                spec_NN = (-0.00053943 * Td*Td + 3.2028 * Td + 1589.2)/1000.               ! Russ 10-2-03
            case(24)
            !Specific Heat of Dowtherm RP, J/kg/K
                spec_NN = (-0.0000031915 * Td**2 + 2.977 * Td + 1560.8)/1000.       !Russ 10-2-03
            case(25)
            !Heat Capacity of HITEC XL Nitrate Salt, [J/kg/K]
                spec_NN = dmax1(1536 - 0.2624 * Td - 0.0001139 * Td * Td,1000.d0)/1000.
            case(26)    ! Argon
                spec_NN = 0.5203 !Cp only, Cv is different
            case(27)    ! Hydrogen
                spec_NN = dmin1(dmax1(-45.4022 + 0.690156*T - 0.00327354*T*T + 0.00000817326*T*T*T - 1.13234E-08*T*T*T*T + 8.24995E-12*T*T*T*T*T - 2.46804E-15*T*T*T*T*T*T,11.3d0),14.7d0)
            case(28)    !T-91 Steel: "Thermo hydraulic optimisation of the EURISOL DS target" - Paul Scherrer Institut
                spec_NN = 0.0004*Td*Td + 0.2473*Td + 450.08
            case(29)    !Therminol 66: Reference: Therminol Reference Disk by Solutia: http://www.therminol.com/pages/tools/toolscd.asp
                spec_NN = 0.0036*Td + 1.4801   
            case(30)    !Therminol 59: Reference: Therminol Reference Disk by Solutia: http://www.therminol.com/pages/tools/toolscd.asp
                spec_NN = 0.0033*Td + 1.6132
            case(31:35)	
            continue
            !!case(36:) !Any integer greater than 35
            !!!Call the user-defined property table
            !!lb=fl_bounds(fnum-35)
            !!ub=fl_bounds(fnum-35+1)-1
            !!if(ub.lt.lb) ub=size(fprop(1,:))
            !!dxx(:)=fprop(1,lb:ub)
            !!dyy(:)=fprop(2,lb:ub)
            !!call interp(Td,size(dxx),dxx,dyy,Gjsav,spec_NN)
            !!        if((Gjsav.eq.ub).or.(Gjsav.eq.lb)) dum=t_warn(Td,dxx(lb),dxx(ub),"User-specified fluid")
            ! case(36) !36-User defined SF HTF
            !     call NR_LINEAR_INTERPOLATION_00(Td,size(T_SF_HTF),T_SF_HTF,spec_NN_SF_HTF,HTFPropsav,spec_NN) !Specific heat, Td is in [C], cp in [kJ/kg/K]
            ! case(37) !37-User defined TES HTF
            !     call NR_LINEAR_INTERPOLATION_00(Td,size(T_TES_HTF),T_TES_HTF,spec_NN_TES_HTF,HTFPropsav,spec_NN) !Specific heat, Td is in [C], cp in [kJ/kg/K]
            case(40)
            !Specific Heat of Dowtherm A, kJ/kg/K
                spec_NN = 1.47524 + 0.00368606*Td - 0.00000516458*Td**2 + 8.99399E-09*Td**3 !Specific heat, Td is in [C], cp in [kJ/kg/K]
            end select

        end function

        !
        !*************** Dowtherm A **************************
        !
        !Enthalpy of Dowtherm A [J/kg]
        Double Precision Function H_Dowtherm_A(T) !T [K]
         ! chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.appliedthermalfluids.com/wp-content/uploads/2018/02/Dowtherm-A-heat-transfer-fluid-TDS.pdf
        implicit none
        Double Precision T, Td
        Td = T - 273.15 ! [C]
        !H_Dowtherm_A = (-19.8113 + 1.50647*T + 0.00144152*T**2) * 1000      !
        !H_Dowtherm_A = (-38.0792 + 1.50904*T + 0.00142671*T**2) * 1000
        H_Dowtherm_A = (-12.7078 + 1.481714*Td + 0.0014292857*Td**2) * 1000 ! [J/kg]      
        End Function

end module Header_functions


module NN_data
    implicit none

    Double precision, dimension(:,:,:), allocatable :: w1, w2, w3, w4, w5, w6, w7, w8
    Double precision, dimension(:,:), allocatable :: b1, b2, b3, b4, b5, b6, b7, B8
    Double precision, dimension(16, 4) :: minMax
    

    contains
        subroutine unload_NN()
            deallocate (B1)
            deallocate (B2)
            deallocate (B3)
            deallocate (B4)
            deallocate (B5)
            deallocate (B6)
            deallocate (B7)
            deallocate (B8)
            deallocate (W1)
            deallocate (W2)
            deallocate (W3)
            deallocate (W4)
            deallocate (W5)
            deallocate (W6)
            deallocate (W7)
            deallocate (W8)
        end subroutine unload_NN


        subroutine load_NN(base)
            character(len=20) :: base
            character(len=40), dimension(4) :: keys
            integer :: n, k, row, col, io, cc, curr_NN, nn
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!! Load in minMax data and NN's !!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            keys = ["pristine_NN", "vacuumLost_NN", "brokenGlass_NN", "H2_NN"]
            do n = 1, 4
                keys(n) = TRIM(keys(n))//TRIM(base)//".txt"
            end do
            ! keys = ["pristine_NNv2.txt", "vacuumLost_NNv2.txt", "brokenGlass_NNv2.txt", "H2_NNv2.txt"]
            !keys = ["pristine_NN_v3.txt", "vacuumLost_NN_v3.txt", "brokenGlass_NN_v3.txt", "H2_NN_v3.txt"]
            do n = 1,4
                nn = 1
                row = -1
                curr_NN = 0
            OPEN(1, File =keys(n))
                do k = 1, 1000
                    ! Read in min-max data for the state
                    IF (k==1) THEN 
                        READ(1, *, IOSTAT=io), minMax(:,n)

                    ! Read in the dimensions for the current NN
                    ELSE IF (nn==row+2) THEN
                        READ(1, *, IOSTAT=io), row, col
                        nn = 1
                        curr_NN = curr_NN + 1
                        IF (io < 0) THEN ! AT END OF FILE
                            EXIT
                        END IF

                    !!! First Layer
                    ELSE IF (curr_NN==1) THEN
                        ! allocate memory
                        If (nn == 1 .AND. n == 1) THEN
                            allocate (w1(col, row, 4)) ! Transpose data
                            allocate (b1(row, 4))
                        END IF

                        ! Load in biases
                        IF (nn==row+1) THEN
                            READ(1, *,IOSTAT=io), b1(:,n)
                            nn = nn + 1
                        ! Load in weights
                        ELSE
                            READ(1, *, IOSTAT=io), w1(:, nn, n)
                            nn = nn + 1
                        END IF

                    !!! Second Layer
                    ELSE IF (curr_NN==2) THEN
                        ! allocate memory
                        If (nn == 1 .AND. n == 1) THEN
                            allocate (w2(col, row, 4)) ! Transpose data
                            allocate (b2(row, 4))
                        END IF

                        ! Load in biases
                        IF (nn==row+1) THEN
                            READ(1, *,IOSTAT=io), b2(:,n)
                            nn = nn + 1
                        ! Load in weights
                        ELSE
                            READ(1, *, IOSTAT=io), w2(:, nn, n)
                            nn = nn + 1
                        END IF

                    !!! Third Layer
                    ELSE IF (curr_NN==3) THEN
                        ! allocate memory
                        If (nn == 1 .AND. n == 1) THEN
                            allocate (w3(col, row, 4)) ! Transpose data
                            allocate (b3(row, 4))
                        END IF

                        ! Load in biases
                        IF (nn==row+1) THEN
                            READ(1, *,IOSTAT=io), b3(:,n)
                            nn = nn + 1
                        ! Load in weights
                        ELSE
                            READ(1, *, IOSTAT=io), w3(:, nn, n)
                            nn = nn + 1
                        END IF

                    !!! Fourth Layer
                    ELSE IF (curr_NN==4) THEN
                        ! allocate memory
                        If (nn == 1 .AND. n == 1) THEN
                            allocate (w4(col, row, 4)) ! Transpose data
                            allocate (b4(row, 4))
                        END IF

                        ! Load in biases
                        IF (nn==row+1) THEN
                            READ(1, *,IOSTAT=io), b4(:,n)
                            nn = nn + 1
                        ! Load in weights
                        ELSE
                            READ(1, *, IOSTAT=io), w4(:, nn, n)
                            nn = nn + 1
                        END IF

                    !!! Fifth Layer
                    ELSE IF (curr_NN==5) THEN
                        ! allocate memory
                        If (nn == 1 .AND. n == 1) THEN
                            allocate (w5(col, row, 4)) ! Transpose data
                            allocate (b5(row, 4))
                        END IF

                        ! Load in biases
                        IF (nn==row+1) THEN
                            READ(1, *,IOSTAT=io), b5(:,n)
                            nn = nn + 1
                        ! Load in weights
                        ELSE
                            READ(1, *, IOSTAT=io), w5(:, nn, n)
                            nn = nn + 1
                        END IF

                    !!! Sixth Layer
                    ELSE IF (curr_NN==6) THEN
                        ! allocate memory
                        If (nn == 1 .AND. n == 1) THEN
                            allocate (w6(col, row, 4)) ! Transpose data
                            allocate (b6(row, 4))
                        END IF

                        ! Load in biases
                        IF (nn==row+1) THEN
                            READ(1, *,IOSTAT=io), b6(:,n)
                            nn = nn + 1
                        ! Load in weights
                        ELSE
                            READ(1, *, IOSTAT=io), w6(:, nn, n)
                            nn = nn + 1
                        END IF

                    !!! Seventh Layer
                    ELSE IF (curr_NN==7) THEN
                        ! allocate memory
                        If (nn == 1 .AND. n == 1) THEN
                            allocate (w7(col, row, 4)) ! Transpose data
                            allocate (b7(row, 4))
                        END IF

                        ! Load in biases
                        IF (nn==row+1) THEN
                            READ(1, *,IOSTAT=io), b7(:,n)
                            nn = nn + 1
                        ! Load in weights
                        ELSE
                            READ(1, *, IOSTAT=io), w7(:, nn, n)
                            nn = nn + 1
                        END IF

                    !!! Eighth Layer
                    ELSE 
                        ! allocate memory
                        If (nn == 1 .AND. n == 1) THEN
                            allocate (w8(col, row, 4)) ! Tranpspose data
                            allocate (b8(row, 4))
                        END IF

                        ! Load in biases
                        IF (nn==row+1) THEN
                            READ(1, *,IOSTAT=io), b8(:,n)
                            nn = nn + 1
                        ! Load in weights
                        ELSE
                            READ(1, *, IOSTAT=io), w8(:, nn, n)
                            nn = nn + 1
                        END IF


                    END IF

                end do
            CLOSE(1)
            end do
        end subroutine load_NN
end module NN_data




module NN_functions_static
    USE NN_data
    implicit none

    contains
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! THIS FUNCTION LOADS THE NEURAL NETWORK INTO THE SCEOPE OF NN_FUNCTIONS MODULE
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function load_NN_toMod(base)result(x)
            implicit none
            Double precision, dimension(16, 4) :: x
            character(len=20) :: base
            Call load_NN(base)
            x = minMax
        end function load_NN_toMod

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! THIS FUNCTION DEALLOCATES THE MEMEORY ASSOCIATED WITH STORING THE NEURAL NETWORK
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine unload_NN_toMod()
            Call unload_NN()
        end subroutine unload_NN_toMod

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! THIS FUNCTION PERFORMS THE COMUTATION OF FORWARD PROPOGATION USING A RELU ACTIVATION FUNCTION
        ! 
        !   val = relu(curr_feat*w + b)
        !
        ! Inputs: 
        !        curr_feat: Feature matrix (n_curr x n_in)
        !        n_curr: number of rows in feature matrix
        !        n_in: number of columns in feature matrix
        !        n_out: number of columns in output matrix
        !        w: weight matrix (n_in x n_out)
        !        b: bias vector (n_out)
        ! Output:
        !        vals: result of forward propogation with relu (n_curr x n_out)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function relu(curr_feat, n_curr, n_in, n_out, w, b)result(val)
            implicit none
            integer, intent(in) :: n_curr, n_in, n_out
            double precision, dimension(n_curr, n_in), intent(in) :: curr_feat
            double precision, dimension(n_in, n_out), intent(in) :: w
            double precision, dimension(n_out), intent(in) :: b
            integer :: n, i, j
            double precision, dimension(n_curr, n_out) :: val

            val = matmul(curr_feat, w)
            do n = 1, n_out
                val(:,n) = val(:,n) + b(n)
            end do

            do i = 1,n_curr 
                do j = 1,n_out
                    val(i,j) = Max(val(i,j), 0.d0)
                end do
            end do
        end function relu 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! THIS FUNCTION PERFORMS THE COMUTATION OF FORWARD PROPOGATION USING A Linear ACTIVATION FUNCTION
        ! 
        !   val = (curr_feat*w + b)
        !
        ! Inputs: 
        !        curr_feat: Feature matrix (n_curr x n_in)
        !        n_curr: number of rows in feature matrix
        !        n_in: number of columns in feature matrix
        !        n_out: number of columns in output matrix
        !        w: weight matrix (n_in x n_out)
        !        b: bias vector (n_out)
        ! Output:
        !        vals: result of forward propogation (n_curr x n_out)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function linear(curr_feat, n_curr, n_in, n_out, w, b)result(val)
            implicit none
            integer, intent(in) :: n_curr, n_in, n_out
            double precision, dimension(n_curr, n_in), intent(in) :: curr_feat
            double precision, dimension(n_in, n_out), intent(in) :: w
            double precision, dimension(n_out), intent(in) :: b
            integer :: n, i, j
            double precision, dimension(n_curr, n_out) :: val

            val = matmul(curr_feat, w)
            do n = 1, n_out
                val(:,n) = val(:,n) + b(n)
            end do
        end function linear

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! THIS FUNCTION FORWARD PROPOGATES AN ENTIRE STATE THROUGH THE NEURAL NETWORK
        ! 
        !   heat = NeuralNetwork(features)
        !
        ! Inputs: 
        !        feat: Feature matrix (c_curr, 7))
        !        n_curr: number of collectors in the current state
        !        state: state of collector (1 : pristine, 2 : lost vacuum, 3 : broken glass, 4 : Hydrogen))
        ! Output:
        !        heat: vector output of the neural network (n_curr) [W/m]
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function forward_Prop(feat_in, n_curr, state)result(heat)
            implicit none
            integer, intent(in) :: n_curr, state
            !double precision, dimension(n_curr, 7), intent(in) :: feat
            double precision, dimension(7, n_curr), intent(in) :: feat_in
            integer :: n_in, n_out
            double precision, dimension(n_curr) :: heat
            integer, dimension(2) :: dim
            double precision, dimension(n_curr, 7) :: curr_mat, new_mat
            double precision, dimension(n_curr, 7) :: feat
            
            feat = transpose(feat_in)
            ! Layer 1
            dim = Shape(w1(:,:,state))
            n_in = dim(1)
            n_out = dim(2)
            new_mat(:,1:n_out) = relu(feat, n_curr, n_in, n_out, w1(:,:,state), b1(:,state))
            curr_mat(:,1:n_out) = new_mat(:,1:n_out)
            

            ! Layer 2
            dim = Shape(w2(:,:,state))
            n_in = dim(1)
            n_out = dim(2)
            new_mat(:,1:n_out) = relu(curr_mat(:,1:n_in), n_curr, n_in, n_out, w2(:,:,state), b2(:,state))
            curr_mat(:,1:n_out) = new_mat(:,1:n_out)

            ! Layer 3
            dim = Shape(w3(:,:,state))
            n_in = dim(1)
            n_out = dim(2)
            new_mat(:,1:n_out) = relu(curr_mat(:,1:n_in), n_curr, n_in, n_out, w3(:,:,state), b3(:,state))
            curr_mat(:,1:n_out) = new_mat(:,1:n_out)

            ! Layer 4
            dim = Shape(w4(:,:,state))
            n_in = dim(1)
            n_out = dim(2)
            new_mat(:,1:n_out) = relu(curr_mat(:,1:n_in), n_curr, n_in, n_out, w4(:,:,state), b4(:,state))
            curr_mat(:,1:n_out) = new_mat(:,1:n_out)

            ! Layer 5
            dim = Shape(w5(:,:,state))
            n_in = dim(1)
            n_out = dim(2)
            new_mat(:,1:n_out) = relu(curr_mat(:,1:n_in), n_curr, n_in, n_out, w5(:,:,state), b5(:,state))
            curr_mat(:,1:n_out) = new_mat(:,1:n_out)

            ! Layer 6
            dim = Shape(w6(:,:,state))
            n_in = dim(1)
            n_out = dim(2)
            new_mat(:,1:n_out) = relu(curr_mat(:,1:n_in), n_curr, n_in, n_out, w6(:,:,state), b6(:,state))
            curr_mat(:,1:n_out) = new_mat(:,1:n_out)

            ! Layer 7
            dim = Shape(w7(:,:,state))
            n_in = dim(1)
            n_out = dim(2)
            new_mat(:,1:n_out) = relu(curr_mat(:,1:n_in), n_curr, n_in, n_out, w7(:,:,state), b7(:,state))
            curr_mat(:,1:n_out) = new_mat(:,1:n_out)

            ! Layer 8
            dim = Shape(w8(:,:,state))
            n_in = dim(1)
            n_out = dim(2)
            new_mat(:,1:n_out) = linear(curr_mat(:,1:n_in), n_curr, n_in, n_out, w8(:,:,state), b8(:,state)) ! NOTE: USING LINEAR FOR LAST LAYER
            curr_mat(:,1) = new_mat(:,1)

            heat = curr_mat(:,1)
            ! Undo normalization and get into units of W/m
            heat = heat - 0.1d0
            heat = heat * (minMax(15, state) - minMax(16, state)) + minMax(16, state)
        end function forward_Prop


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! THIS FUNCTION COMPUTES THE TEMPERATURE DERIVATIVE OF THE HTF NODES
        ! 
        !   
        !
        ! Inputs: 
        !        t: nodal temperatures
        !        t_bar: control volume average temperatures
        !        features: normalized NN features (n_nodes-1, 7)
        !        m_dot: mass flow rate [kg/s]
        !        L_segment: Length of the control volume [m]
        !        Vol: Volume of the control volume [m^3]
        !        n_nodes: Number of nodes in the loop [-] (should be one more than the number of control volumes)
        !        nCV_state: Vector with number of control volumes in each collector state (4) 
        !        inds_pristine: indices of control volumes in pristine state
        !        inds_lVacuum: indices of control volumes with lost vacuum
        !        inds_bGlass: indices of control volumes with broken glass
        !        inds_H2: indices of control volumes with hydrogen in annulus
        ! Output:
        !        deriv: Vector of temperature derivatives of all nodes [K/s]
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function dt_dtime_NN(t, t_bar, features, m_dot, mc_sf, L_segment, Vol, n_nodes, nCV_state, inds_pristine, inds_lVacuum, inds_bGlass, inds_H2)result(dt_dtheta)
            implicit none
            integer, intent(in) :: n_nodes
            Double precision, intent(in) :: L_segment, Vol, m_dot, mc_sf
            Double precision, dimension(n_nodes), intent(in) :: t
            integer, dimension(4), intent(in) :: nCV_state
            Double precision, dimension(n_nodes-1, 7), intent(in) :: features
            Double precision, dimension(n_nodes-1), intent(in) :: t_bar
            Integer, dimension(n_nodes-1), intent(in) :: inds_pristine
            Integer, dimension(n_nodes-1), intent(in) :: inds_lVacuum
            Integer, dimension(n_nodes-1), intent(in) :: inds_bGlass
            Integer, dimension(n_nodes-1), intent(in) :: inds_H2
            Double precision, dimension(7, n_nodes-1) :: feat_temp
            integer :: n, state, zz
            double precision:: start, finish, time_diff
            DOUBLE PRECISION :: fnumd = 40
            Double precision, dimension(n_nodes) :: dt_dtheta, h
            Double precision, dimension(n_nodes-1) :: c_bar, dt_dtheta_bar, q_in, rho
            Integer, dimension(n_nodes-1) :: curr_inds
            ! Forward propogate NN
            do state = 1,4
                n = nCV_state(state)
                if (state == 1)THEN
                    curr_inds(1:n) = inds_pristine(1:n)
                ELSE IF (state == 2)THEN
                    curr_inds(1:n) = inds_lVacuum(1:n)
                ELSE IF (state == 3)THEN
                    curr_inds(1:n) = inds_bGlass(1:n)
                ELSE
                    curr_inds(1:n) = inds_H2(1:n)
                END IF
                IF (n>0) THEN
                    feat_temp(:, 1:n) = Transpose(features(curr_inds(1:n),:))
                    q_in(curr_inds(1:n)) = forward_Prop(feat_temp(:, 1:n), n, state)*L_segment ! [W]
                END IF
            end do

            ! Solve for HTF Properties
            do n = 1,n_nodes-1
                rho(n) = den_NN(fnumd, t_bar(n), 0.d0)
                c_bar(n) = 1000.d0*spec_NN(fnumd, t_bar(n), 0.d0)
                h(n) = H_Dowtherm_A(t(n))
            end do
            h(n) = H_Dowtherm_A(t(n))
            n= 1
            ! Compute derivative
            !dt_dtheta_bar = 1.d0/mc_sf/(rho*Vol*c_bar)*(q_in + m_dot*(c_1*t(1:n_nodes-1) - c_2*t(2:n_nodes)))
            dt_dtheta_bar = 1.d0/mc_sf/(rho*Vol*c_bar)*(q_in + m_dot*c_bar*(t(1:n_nodes-1) - t(2:n_nodes)))
            dt_dtheta(1) = 0.d0
            dt_dtheta(2:n_nodes-1) = 1.d0/2.d0 * (dt_dtheta_bar(1:n_nodes-2) + dt_dtheta_bar(2:n_nodes-1))
            dt_dtheta(n_nodes) = dt_dtheta_bar(n_nodes-1)
        end function dt_dtime_NN



        !  DENSITY AND SPECIFIC HEAT FUNCTIONS TAKEN FROM SAM
        ! ||||||||||||||||||||||||||||||||||||||||||||||||||||
        ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        double precision function den_NN(fnumd,T,P)
            !This function accepts as inputs temperature [K] and pressure [Pa]
            !This function outputs in units of [kg/m^3]
            double precision::xlo,xhi, Dens_fluid, Td,HTFPropsav
            double precision::T,P,fnumd
            !!double precision,dimension(size(fprop(1,:)))::dxx,dyy !Create dummy arrays
            integer::fnum,lb,ub,dum,t_warn
            !den_NN=1.
            fnum=nint(fnumd)
            Td=T-273.15             !Convert from K to C
                
            select case(fnum)
            case(1)   !    1.) Air
            den_NN = P/(287.*T)
            case(2)   !    2.) Stainless_AISI316
                den_NN=8349.38 - 0.341708*T - 0.0000865128*T*T  !EES
            case(3)   !    3.) Water (liquid)
                den_NN = 1000 
            case(4)   !    4.) Steam
                continue
            case(5)   !    5.) CO2
                continue
            case(6)   !    6.) Salt (68% KCl, 32% MgCl2)
            den_NN = 1E-10*T*T*T - 3E-07*T*T - 0.4739*T + 2384.2
            case(7)   !    7.) Salt (8% NaF, 92% NaBF4)
            den_NN = 8E-09*T*T*T - 2E-05*T*T - 0.6867*T + 2438.5
            case(8)   !    8.) Salt (25% KF, 75% KBF4)
            den_NN = 2E-08*T*T*T - 6E-05*T*T - 0.7701*T + 2466.1
            case(9)   !    9.) Salt (31% RbF, 69% RbBF4)
            den_NN = -1E-08*T*T*T + 4E-05*T*T - 1.0836*T + 3242.6
            case(10)   !    10.) Salt (46.5% LiF, 11.5%NaF, 42%KF)
            den_NN =  -2E-09*T*T*T + 1E-05*T*T - 0.7427*T + 2734.7
            case(11)   !    11.) Salt (49% LiF, 29% NaF, 29% ZrF4)
            den_NN = -2E-11*T*T*T + 1E-07*T*T - 0.5172*T + 3674.3
            case(12)   !    12.) Salt (58% KF, 42% ZrF4)
            den_NN =  -6E-10*T*T*T + 4E-06*T*T - 0.8931*T + 3661.3
            case(13)   !    13.) Salt (58% LiCl, 42% RbCl)
            den_NN = -8E-10*T*T*T + 1E-06*T*T - 0.689*T + 2929.5
            case(14)   !    14.) Salt (58% NaCl, 42% MgCl2)
            den_NN = -5E-09*T*T*T + 2E-05*T*T - 0.5298*T + 2444.1
            case(15)   !    15.) Salt (59.5% LiCl, 40.5% KCl)
            den_NN = 1E-09*T*T*T - 5E-06*T*T - 0.864*T + 2112.6
            case(16)   !    16.) Salt (59.5% NaF, 40.5% ZrF4)
            den_NN =  -5E-09*T*T*T + 2E-05*T*T - 0.9144*T + 3837.
            case(17)   !    17.) Salt (60% NaNO3, 40% KNO3)
            den_NN = dmax1(-1E-07*T*T*T + 0.0002*T*T - 0.7875*T + 2299.4,1000.d0)
            case(18)
            !den_NN of Nitrate Salt, [kg/m3]
            den_NN = dmax1(2090 - 0.636 * (T-273.15),1000.d0)
            case(19)
            !den_NN of Caloria HT 43 [kg/m3]
            den_NN = dmax1(885 - 0.6617 * Td - 0.0001265 * Td*Td,100.d0)
            case(20)
            !den_NN of HITEC XL Nitrate Salt, [kg/m3]
            den_NN = dmax1(2240 - 0.8266 * Td,800.d0)
            case(21)
            !den_NN of Therminol Oil [kg/m3]
            den_NN = dmax1(1074.0 - 0.6367 * Td - 0.0007762 * Td*Td,400.d0)
            case(22)
            !den_NN of HITEC Salt, [kg/m3]
            den_NN = dmax1(2080 - 0.733 * Td,1000.d0)
            case(23)
            !den_NN of Dowtherm Q [kg/m3]
            den_NN = dmax1(-0.757332 * Td + 980.787,100.d0)                               ! Russ 10-2-03
            case(24)
            !den_NN of Dowtherm RP [kg/m3]
            den_NN = dmax1(-0.000186495 * Td*Td - 0.668337 * Td + 1042.11,200.d0)  !Russ 10-2-03
            case(25)
            !den_NN of HITEC XL Nitrate Salt, [kg/m^3]
            den_NN = dmax1(2240 - 0.8266 * Td,800.d0)
            case(26) !Argon
            den_NN = dmax1(P/(208.13*T),1.e-10)
            case(27) !Hydrogen
            den_NN = dmax1(P/(4124.*T),1.e-10)
            case(28)    !T-91 Steel: "Thermo hydraulic optimisation of the EURISOL DS target" - Paul Scherrer Institut
            den_NN = -0.3289*Td + 7742.5
            case(29)    !Therminol 66: Reference: Therminol Reference Disk by Solutia: http://www.therminol.com/pages/tools/toolscd.asp
            den_NN = -0.7146*Td + 1024.8
            case(30)    !Therminol 59: Reference: Therminol Reference Disk by Solutia: http://www.therminol.com/pages/tools/toolscd.asp
            den_NN = -0.0003*Td*Td - 0.6963*Td + 988.44
            case(31:35) 
            continue !no informaion
            !!case(36:) !Any integer greater than 35
            !!!Call the user-defined property table
            !!lb=fl_bounds(fnum-35)
            !!ub=fl_bounds(fnum-35+1)-1
            !!if(ub.lt.lb) ub=size(fprop(1,:))
            !!dxx(:)=fprop(1,lb:ub)
            !!dyy(:)=fprop(3,lb:ub)
            !!call interp(Td,size(dxx),dxx,dyy,Gjsav,den_NN)
            !!if((Gjsav.eq.ub).or.(Gjsav.eq.lb)) dum=t_warn(Td,dxx(lb),dxx(ub),"User-specified fluid")
            ! case(36) !36-User defined SF HTF
            ! call NR_LINEAR_INTERPOLATION_00(Td,size(T_SF_HTF),T_SF_HTF,den_NN_SF_HTF,HTFPropsav,den_NN) !den_NN, Td is in [C], den_NN in [kg/m3]
            ! case(37) !37-User defined TES HTF
            ! call NR_LINEAR_INTERPOLATION_00(Td,size(T_TES_HTF),T_TES_HTF,den_NN_TES_HTF,HTFPropsav,den_NN) !den_NN, Td is in [C], den_NN in [kg/m3]
            case(40)
            !den_NN of Dowtherm A [kg/m3]
            den_NN = dmax1(1063.61 - 0.605235*Td - 0.000860877*Td*Td,400.d0)!den_NN, Td is in [C], den_NN in [kg/m3]
            end select

        end function

        double precision function spec_NN(fnumd,T,P)
            !This function accepts as inputs temperature [K] and pressure [Pa]
            !This function outputs in units of [kJ/kg-K]
            double precision::xlo,xhi, Td,HTFPropsav
            double precision,intent(in)::T,P,fnumd
            !!double precision,dimension(size(fprop(1,:)))::dxx,dyy !Create dummy arrays
            integer::fnum,lb,ub,dum,t_warn

            spec_NN=1.
            fnum=nint(fnumd)
            Td = T - 273.15
            select case(fnum)
            case(1)   !    1.) Air
                spec_NN = 1.03749 - 0.000305497*T + 7.49335E-07*T*T - 3.39363E-10*T*T*T
            !spec_NN = 1.03749 - 0.000305497*T + 7.49335E-07*T*T - 3.39363E-10*T*T*T
            case(2)   !    2.) Stainless_AISI316
                spec_NN = 0.368455 + 0.000399548*T - 1.70558E-07*T*T !EES
            case(3)   !    3.) Water (liquid)
                spec_NN = 4.181d0  !mjw 8.1.11 
            case(4)   !    4.) Steam
                continue
            case(5)   !    5.) CO2
                continue
            case(6)   !    6.) Salt (68% KCl, 32% MgCl2)
                spec_NN = 1.156
            case(7)   !    7.) Salt (8% NaF, 92% NaBF4)
                spec_NN = 1.507
            case(8)   !    8.) Salt (25% KF, 75% KBF4)
                spec_NN = 1.306
            case(9)   !    9.) Salt (31% RbF, 69% RbBF4)
                spec_NN = 9.127
            case(10)   !    10.) Salt (46.5% LiF, 11.5%NaF, 42%KF)
                spec_NN = 2.010
            case(11)   !    11.) Salt (49% LiF, 29% NaF, 29% ZrF4)
                spec_NN = 1.239
            case(12)   !    12.) Salt (58% KF, 42% ZrF4)
                spec_NN = 1.051
            case(13)   !    13.) Salt (58% LiCl, 42% RbCl)
                spec_NN = 8.918
            case(14)   !    14.) Salt (58% NaCl, 42% MgCl2)
                spec_NN = 1.080
            case(15)   !    15.) Salt (59.5% LiCl, 40.5% KCl)
                spec_NN = 1.202
            case(16)   !    16.) Salt (59.5% NaF, 40.5% ZrF4)
                spec_NN = 1.172
            case(17)   !    17.) Salt (60% NaNO3, 40% KNO3)
                spec_NN = -1E-10*T*T*T + 2E-07*T*T + 5E-06*T + 1.4387
            case(18) !Heat Capacity of Nitrate Salt, [J/kg/K]
                spec_NN = (1443. + 0.172 * (T-273.15))/1000.d0
            case(19)
            !Specific Heat of Caloria HT 43 [J/kgC]
                spec_NN = (3.88 * (T-273.15) + 1606.0)/1000.
            case(20)
            !Heat Capacity of HITEC XL Nitrate Salt, [J/kg/K]
                spec_NN = dmax1(1536 - 0.2624 * Td - 0.0001139 * Td * Td,1000.d0)/1000.
            case(21)
            !Specific Heat of Therminol Oil, kJ/kg/K
                spec_NN = (1.509 + 0.002496 * Td + 0.0000007888 * Td*Td)
            case(22)
            !Heat Capacity of HITEC Salt, [J/kg/K]
                spec_NN = (1560 - 0.0 * Td)/1000.
            case(23)
            !Specific Heat of Dowtherm Q, J/kg/K
                spec_NN = (-0.00053943 * Td*Td + 3.2028 * Td + 1589.2)/1000.               ! Russ 10-2-03
            case(24)
            !Specific Heat of Dowtherm RP, J/kg/K
                spec_NN = (-0.0000031915 * Td**2 + 2.977 * Td + 1560.8)/1000.       !Russ 10-2-03
            case(25)
            !Heat Capacity of HITEC XL Nitrate Salt, [J/kg/K]
                spec_NN = dmax1(1536 - 0.2624 * Td - 0.0001139 * Td * Td,1000.d0)/1000.
            case(26)    ! Argon
                spec_NN = 0.5203 !Cp only, Cv is different
            case(27)    ! Hydrogen
                spec_NN = dmin1(dmax1(-45.4022 + 0.690156*T - 0.00327354*T*T + 0.00000817326*T*T*T - 1.13234E-08*T*T*T*T + 8.24995E-12*T*T*T*T*T - 2.46804E-15*T*T*T*T*T*T,11.3d0),14.7d0)
            case(28)    !T-91 Steel: "Thermo hydraulic optimisation of the EURISOL DS target" - Paul Scherrer Institut
                spec_NN = 0.0004*Td*Td + 0.2473*Td + 450.08
            case(29)    !Therminol 66: Reference: Therminol Reference Disk by Solutia: http://www.therminol.com/pages/tools/toolscd.asp
                spec_NN = 0.0036*Td + 1.4801   
            case(30)    !Therminol 59: Reference: Therminol Reference Disk by Solutia: http://www.therminol.com/pages/tools/toolscd.asp
                spec_NN = 0.0033*Td + 1.6132
            case(31:35)	
            continue
            !!case(36:) !Any integer greater than 35
            !!!Call the user-defined property table
            !!lb=fl_bounds(fnum-35)
            !!ub=fl_bounds(fnum-35+1)-1
            !!if(ub.lt.lb) ub=size(fprop(1,:))
            !!dxx(:)=fprop(1,lb:ub)
            !!dyy(:)=fprop(2,lb:ub)
            !!call interp(Td,size(dxx),dxx,dyy,Gjsav,spec_NN)
            !!        if((Gjsav.eq.ub).or.(Gjsav.eq.lb)) dum=t_warn(Td,dxx(lb),dxx(ub),"User-specified fluid")
            ! case(36) !36-User defined SF HTF
            !     call NR_LINEAR_INTERPOLATION_00(Td,size(T_SF_HTF),T_SF_HTF,spec_NN_SF_HTF,HTFPropsav,spec_NN) !Specific heat, Td is in [C], cp in [kJ/kg/K]
            ! case(37) !37-User defined TES HTF
            !     call NR_LINEAR_INTERPOLATION_00(Td,size(T_TES_HTF),T_TES_HTF,spec_NN_TES_HTF,HTFPropsav,spec_NN) !Specific heat, Td is in [C], cp in [kJ/kg/K]
            case(40)
            !Specific Heat of Dowtherm A, kJ/kg/K
                spec_NN = 1.47524 + 0.00368606*Td - 0.00000516458*Td**2 + 8.99399E-09*Td**3 !Specific heat, Td is in [C], cp in [kJ/kg/K]
            end select

        end function

        !
        !*************** Dowtherm A **************************
        !
        !Enthalpy of Dowtherm A [J/kg]
        Double Precision Function H_Dowtherm_A(T) !T [K]
         ! chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.appliedthermalfluids.com/wp-content/uploads/2018/02/Dowtherm-A-heat-transfer-fluid-TDS.pdf
        implicit none
        Double Precision T, Td
        Td = T - 273.15 ! [C]
        !H_Dowtherm_A = (-19.8113 + 1.50647*T + 0.00144152*T**2) * 1000      !
        !H_Dowtherm_A = (-38.0792 + 1.50904*T + 0.00142671*T**2) * 1000
        H_Dowtherm_A = (-12.7078 + 1.481714*Td + 0.0014292857*Td**2) * 1000 ! [J/kg]      
        End Function



end Module NN_functions_static


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module SF_data

    USE NN_functions_static
    implicit none

    Double precision, dimension (:,:,:), allocatable :: t_sf
    Integer, dimension(:,:), allocatable :: inds_header_in, inds_right, inds_left
    Double precision, dimension (:,:), allocatable :: Vol_inlet, t_header_inlet, L_cv_inlet, D_inlet
    Double precision, dimension (:,:), allocatable :: Vol_return, t_header_return, L_cv_return, D_return
    Double precision, dimension(:), allocatable :: m_dots_in, k1_ih, k2_ih, k3_ih, k4_ih, t_bar_inlet, t_bar_hat_inlet, t_hat_inlet
    Double precision, dimension(:), allocatable :: m_dots_return, k1_rh, k2_rh, k3_rh, k4_rh, t_bar_return, t_bar_hat_return, t_hat_return
    Double precision, dimension(16, 4) :: minMax2
    Integer, dimension(:, :, :), allocatable :: nCV_state, inds_pristine, inds_lVacuum, inds_bGlass, inds_H2, H2_pressure
    Double precision, dimension(:,:), allocatable :: m_dot_var, defocus_mode, Time_df
    Double precision, dimension(:), allocatable :: t_bar, t_bar_hat, t_hat, k1, k2, k3, k4, dni_array
    Double precision, dimension(:), allocatable :: m_left, m_right, t_hold_l, t_hold_r, t_hold, L_cv_hold, vol_hold
    Integer, dimension(:), allocatable :: inds_hold
    Double precision :: L_segment, Vol
    Double precision, dimension(:,:), allocatable :: features
    Double precision, dimension(:,:,:), allocatable :: r_number
    Integer, dimension(:), allocatable :: num_cv_header
    Double precision, dimension(:), allocatable :: mass_HTF_hold, t4Ave_hold
    Double precision, dimension(:,:), allocatable :: t_bar_sf
    
    contains 
        subroutine allocate_memory_sf(max_loops, n_sectors, n_nodes_per_loop, n_SCA, base)
            integer :: max_loops, n_sectors, n_nodes_per_loop, n_SCA
            integer :: cc, k, n_node_header, n_cv_header
            character(len=20) :: base

            ! Determine Inlet and Return Header Matrix Dimensions
            cc = 0
            do k = 1, max_loops/2
                if(k==1)then
                    cc = cc + 1
                else if(mod(k,2) == 0)then
                    cc = cc + 1
                else
                    cc = cc + 2
                end if
            end do
            n_node_header = cc
            n_cv_header = cc -  1

            ! Allocate Memory for Inlet Headers
            allocate (num_cv_header(n_sectors))
            allocate (Vol_inlet(n_cv_header, n_sectors))
            allocate (D_inlet(n_cv_header, n_sectors))
            allocate (t_header_inlet(n_node_header, n_sectors))
            allocate (m_dots_in(n_cv_header))
            allocate (k1_ih(n_node_header))
            allocate (k2_ih(n_node_header))
            allocate (k3_ih(n_node_header))
            allocate (k4_ih(n_node_header))
            allocate (t_bar_inlet(n_cv_header))
            allocate (t_bar_hat_inlet(n_cv_header))
            allocate (t_hat_inlet(n_node_header))
            allocate (inds_header_in(max_loops/2, n_sectors))
            allocate (L_cv_inlet(n_cv_header, n_sectors))
            

            ! Allocate Memory for Return Headers
            allocate (Vol_return(n_cv_header+1, n_sectors))
            allocate (D_return(n_cv_header+1, n_sectors))
            allocate (t_header_return(n_node_header+1, n_sectors))
            allocate (m_dots_return(n_cv_header+1))
            allocate (k1_rh(n_node_header+1))
            allocate (k2_rh(n_node_header+1))
            allocate (k3_rh(n_node_header+1))
            allocate (k4_rh(n_node_header+1))
            allocate (t_bar_return(n_cv_header+1))
            allocate (t_bar_hat_return(n_cv_header+1))
            allocate (t_hat_return(n_node_header+1))
            allocate (inds_right(max_loops/2, n_sectors))
            allocate (inds_left(max_loops/2, n_sectors))


            ! Allocate Memory for Return Headers
            allocate (L_cv_return(n_cv_header+1, n_sectors))
            
            ! Load in Neural Network Data
            minMax2 = load_NN_toMod(base)

            ! Allocate Memory for receiver states in each loop
            allocate (nCV_state(4, max_loops, n_sectors))
            allocate (inds_pristine(n_nodes_per_loop-1, max_loops, n_sectors))
            allocate (inds_lVacuum(n_nodes_per_loop-1, max_loops, n_sectors))
            allocate (inds_bGlass(n_nodes_per_loop-1, max_loops, n_sectors))
            allocate (inds_H2(n_nodes_per_loop-1, max_loops, n_sectors))
            allocate (H2_pressure(n_nodes_per_loop-1, max_loops, n_sectors))
            allocate (m_dot_var(max_loops, n_sectors))

            ! Allocate Memory for tracking strategy
            allocate (defocus_mode(max_loops, n_sectors))
            allocate (Time_df(max_loops, n_sectors))

            ! Allocate Memory for temperatures in solar field loop
            allocate (t_bar(n_nodes_per_loop-1))
            allocate (t_bar_hat(n_nodes_per_loop-1))
            allocate (t_hat(n_nodes_per_loop))
            allocate (k1(n_nodes_per_loop))
            allocate (k2(n_nodes_per_loop))
            allocate (k3(n_nodes_per_loop))
            allocate (k4(n_nodes_per_loop))
            allocate (dni_array(n_nodes_per_loop -1))
            allocate (t_sf(n_nodes_per_loop, max_loops, n_sectors))

            ! Allocate memory for features matrix 
            allocate (features(n_nodes_per_loop-1, 7))

            ! Allocate memory for SF_Avail mask
            allocate (r_number(n_SCA, max_loops, n_sectors))

            ! Allocate temporary arrays to avoid function call warnings
            allocate (m_left(max_loops/2))
            allocate (m_right(max_loops/2))
            allocate (t_hold_l(max_loops/2))
            allocate (t_hold_r(max_loops/2))
            allocate (t_hold(n_node_header+1))
            allocate (inds_hold(max_loops/2))
            allocate (L_cv_hold(n_cv_header+1))
            allocate (vol_hold(n_cv_header+1))

            ! Allocate memory for HTF mass holder, average sector temperature, and t4ave
            allocate (mass_htf_hold(n_sectors))
            allocate (t_bar_sf(n_nodes_per_loop-1, n_sectors))
            allocate (t4Ave_hold(n_sectors))


        end subroutine allocate_memory_sf


        subroutine deallocate_memory_sf()
            deallocate (t_sf)
            deallocate (num_cv_header)
            deallocate (Vol_inlet, D_inlet, t_header_inlet, inds_header_in, L_cv_inlet)
            deallocate (Vol_return, t_header_return, L_cv_return, D_return)
            deallocate (m_dots_in, k1_ih, k2_ih, k3_ih, k4_ih, t_bar_inlet, t_bar_hat_inlet, t_hat_inlet)
            deallocate (nCV_state, inds_pristine, inds_lVacuum, inds_bGlass, inds_H2, H2_pressure)
            deallocate (m_dot_var, defocus_mode, Time_df)
            deallocate (features)
            deallocate (r_number)
            deallocate (m_left, m_right, t_hold_l, t_hold_r, t_hold, inds_hold, L_cv_hold, vol_hold)
        end subroutine





end module SF_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Solar_Position   
    implicit none

    contains

    function Row_shadow(phi, row_distance, w_ap)result(eta_row)
        double precision :: phi, row_distance, w_ap, eta_row
        eta_row = Abs(cos(phi))*row_distance / w_ap

        if (eta_row > 1.0)then
            eta_row = 1.d0
        else if (eta_row < 0.d0)then
            eta_row = 0.d0
        end if

    end function
    
    function solar_tracking(TimeZone, LongD, LatD, julian_day, time)result(results)
        double precision :: time, TimeZone, LongD, LatD, pi, dec, dec_rad, julian_day, n, EOT, solartime, lat_rad, solaraz_rad
        double precision :: B, B_rad, HourAngle, HourAngle_rad, SolarAlt_rad, SolarAlt, SolzarAz_rad, SolarAz, ColAz_rad, ColAz
        double precision :: SolarZenith_rad, SolarZenith, L_st, phi, CosTh, Theta_rad, Theta
        double precision, dimension(2) :: results

        pi = 3.141592653
        n = julian_day
        L_st = 105.d0

        !3.3-New B per Duffie & Beckman 1.4.2
        B = (n - 1) * 360. / 365. 
        B_rad = B * pi/180.	

        !3.4-Equation of Time in minutes
        EOT = 229.2 * (0.000075 + 0.001868 * Cos(B_rad) - 0.032077 * Sin(B_rad) - 0.014615 * Cos(B_rad * 2) - 0.04089 * Sin(B_rad * 2.))

        !3.5- Declination   (per Duffie & Beckman 1.6.1a)"
        dec = 23.45 * Sin(360. * (284. + n) / 365. * Pi / 180.)
        dec_rad = dec * pi / 180.
        
        ! Solar Time in Hours
        SolarTime = time+ (4.d0*(L_st - LongD) + EOT)/60.d0
        


        !3.14- Calculation of Hour Angle in radians
        HourAngle = (SolarTime - 12.) * 15. 
        HourAngle_rad = HourAngle * pi / 180.

        !3.14- Solar Altitude   (Radians)
        Lat_rad = LatD * pi/180
        SolarAlt_rad = ASIN( Sin(dec_rad) * Sin(Lat_rad) + Cos(Lat_rad) * Cos(Dec_rad) * Cos(HourAngle_rad))
        SolarAlt = SolarAlt_rad * 180./pi
        
        !3.15- Solar azimuth
        SolarAz_rad = sign(1.d0, HourAngle_rad)*abs(acos(dmin1(1.d0,(cos(pi/2.-SolarAlt_rad)*sin(Lat_rad)-sin(Dec_rad))/(sin(pi/2.-SolarAlt_rad)*cos(Lat_rad)))))
        !SolarAz_rad = ACOS((Sin(Dec_rad) * Cos(Lat_rad) - Cos(Dec_rad) * Cos(HourAngle_rad) * Sin(Lat_rad)) / Cos(SolarAlt_rad))	
        !SolarAz_rad = sign(1.d0, HourAngle_rad)
        SolarAz = SolarAz_rad * 180./pi

        !3.16- Solar Zenith
        SolarZenith_rad = Pi / 2. - SolarAlt_rad
        SolarZenith = SolarZenith_rad * 180./pi
        
        !3.17- Calculation of Solar Incidence Angle for Trough (FROM YOUR OWN DERIVATION: ASSUMES HORIZONTAL COLLECTOR ROTATING ABOUT THE N/S AXIS)
        if (HourAngle_rad>0)then
            ColAz_rad = 90*pi/180
        else
            ColAz_rad = -90*pi/180
        end if

        ! Compute Tracking Angle
        phi = atan(tan(SolarZenith_rad)*sin(SolarAz_rad))
        
        if(PHI<0)then
            costh = 1.d0
        endif
        ! Compute cos(theta)
        CosTh = -sin(SolarZenith_rad)*sin(SolarAz_rad)*(-sin(phi)) + cos(SolarZenith_rad)*cos(phi)
        

        theta = ACOS(CosTh)

        results(1) = phi
        results(2) = theta

    end function solar_tracking

end module Solar_Position
