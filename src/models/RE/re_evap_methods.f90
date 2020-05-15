module re_evap_methods
 use typy
 
  
  public :: McGuinessBorder_fcn
  public :: makkink_fcn
  public :: makkinkHasem_fcn
  public :: DoorenbosPruit_fcn
  public :: Abtew_fcn
  public :: PriestleyTaylor_fcn
  public :: Turc_fcn
  public :: Hargreaves_fcn
  public :: HargreaveSamani_fcn
  public :: JensenHaise_fcn
  public :: Thornthwaite_fcn
  public :: Hamon_fcn
  public :: Linacre_fcn
  public :: BlaneyCriddle_fcn
  public :: Oudin_fcn
  public :: Romanenko_fcn
  public :: RomanenkOudin_fcn
  public :: Kharrufa_fcn
  
  !> if units [m/s] conv_const = 1e-3/86400
  !! if units [mm/day] conv_const = 1.0
  !! if different units find on your own
  !<
  
  real(kind=rkind), private :: conv_const = 1.0
  contains


    !> 1. McGuiness and Border Method
    !> val: McGuiness and Border 
    function McGuinessBorder_fcn(t,Rs) result (val) 
      use typy
      real (kind=rkind), intent(in) ::t,Rs
      !x: # of the day in the year 
      !t: average temperature
      !Rs: solar radiation
      real (kind=rkind) :: val
      !McGuiness and Border [mm/day]         
      val = ((0.0082*t)-0.19)*(Rs/1500)*2.54
      val = val*conv_const
    end function McGuinessBorder_fcn
    
   !#################################################################
   !> 2. Makkink Method
   !> val: Makkink Method
    function makkink_fcn(t,Rs) result(val)
      use typy
      real (kind=rkind), intent(in) :: t,Rs
      !x: # of the day in the year 
      !sl: slope of the saturation vapour pressure curve [kPa/°C]
      !Rs: solar radiation [MJ/m2 day]
      !Ps: psichromatic constant [kPa/°C]
      real (kind=rkind) :: val
      real(kind=rkind) :: sl,Ps,Lh
       !slope of the saturation vapour pressure curve
       sl = 33.8639_rkind*(0.05904_rkind*(0.00738_rkind*(t))+0.8072**7_rkind-0.0000342_rkind)
       !Latent heat [in calories per gram]
       Lh = (5.95_rkind-0.51_rkind*(t))
       !psicromatic constant [kPa/°C]
       Ps = (1.013**-3_rkind*101.325_rkind)/(0.622_rkind*Lh)
      !Makkink Method [mm/day]         
      val = 0.61_rkind*(sl/(sl-Ps))*(Rs/58.5_rkind)-0.012_rkind
    end function makkink_fcn
    
    !##################################################################
    !> 3. Makkink-Hasem Method
    !> val: Makkink-Hasem 
    function makkinkHasem_fcn (x,y,z,a,e,Rs,t,t_min,t_max) result (val)
    use typy
    use core_tools
    real (kind=rkind), intent(in) :: y,z,e,a,Rs,t,t_min,t_max
    !sl: slope of the saturation vapour pressure curve [kPa/°C]
    !Rs: solar radiation [MJ/m2 day]
    !Ps: psichromatic constant [mbar/°C]
    integer (kind=ikind), intent(in):: x
     real (kind=rkind) :: val
    real(kind=rkind) :: sl,Ps,Lh,Rn,omega,R_nl,R_so,R_ns,R_a,dr,delta
      !slope of the saturation vapour pressure curve
      sl = 33.8639*(0.05904*(0.00738*(t))+0.8072**7-0.0000342)
      !Latent heat [in calories per gram]
      Lh = (5.95-0.51*(t))
      !psicromatic constant [kPa/°C]
      Ps = (1.013**-3*101.325)/(0.622*Lh)
        !inverse distance Earth-sun [rad]
      dr = 1.0_rkind + 0.033_rkind*cos((2.0_rkind*pi()*x)/365.0_rkind)
      !solar declination [rad]
      delta = 0.409_rkind*sin(((2.0_rkind*pi()*x)/365.0_rkind) -1.39_rkind)
      !sunset hour angle [rad]
      omega = acos(-tan(y)*tan(delta))
      !Extraterrestrial radiation [MJ/m2 day]        
      R_a = ((24.0_rkind*60.0_rkind)/pi())*dr*0.0820_rkind*(omega*sin(y)*sin(delta)&
              + cos(y)*cos(delta)*sin(omega))
      !Clear-sky solar radiation [MJ/m2 day] 
      R_so = (0.75_rkind + z*2e-5)*R_a
      !net shortwave radiation MJ/m2 day]
      R_ns = (1.0_rkind-a)*Rs
      !Net longwave radiation [MJ/m2 day] 
      R_nl = 4.903e-9*((t_min**4.0_rkind + t_max**4.0_rkind)/2.0_rkind)*(0.34_rkind - & 
                0.14_rkind*sqrt(e))*(1.35_rkind*(Rs/R_so) -0.35_rkind)
      !Net radiation [mm/day]         
      Rn = R_ns - R_nl
       
    !Makkink-Hasem Method [mm/day] 
    val = 1.26_rkind*(sl/(sl+Ps))*(Rn/Lh)
    end function makkinkHasem_fcn 
    
    !#!!!!############################################################
    !> 4. Doorenbos and Pruit Method
    !> val: Doorenbos_Pruit
    function DoorenbosPruit_fcn (t,Rs,RH,Ud)  result (val)
     use typy
     real (kind=rkind), intent(in) :: t,Rs,RH,Ud
     !sl: slope of the saturation vapour pressure curve [kPa/°C]
     !Rs: solar radiation [MJ/m2 day]
     !Ps: psichromatic constant [kPa/°C]
     !b: adjustment factor that varies with mean relative humidity and daytime wind speed
     
     real(kind=rkind) :: val
     real(kind=rkind) :: b,sl,Ps,Lh,a
       b = 0.3
       !slope of the saturation vapour pressure curve [kPa/°C]
       sl = 33.8639_rkind*(0.05904_rkind*(0.00738_rkind*(t))+0.8072e7_rkind-0.0000342_rkind)
       !heat latent [kPa/°C]
       Lh = (5.95_rkind-0.51_rkind*(t))
       !psicromatic constant [kPa/°C]
       Ps = (1.013e-3_rkind*101.25_rkind)/(0.622_rkind*Lh)
       !adjustment factor that varies with mean relative humidity and daytime wind speed [m/s]
       a = 1.066_rkind-0.13e-2*(RH)+0.045_rkind*(Ud)-0.20e-3*(RH)*(Ud)-0.315e-4*(RH)**2-0.11e-2*(Ud)**2
     !Doorenbos and Pruit method [mm/day] 
     val = a*((sl/(sl+Ps))*Rs)+b
    end function DoorenbosPruit_fcn
    
    !#!!!!###########################################################
    !> 5. Abtew Method
    !> val: Abtew
    function Abtew_fcn (Rs,t) result (val)
      use typy
      real (kind=rkind), intent (in) :: Rs,t
      real (kind=rkind) :: val
      real(kind=rkind) :: Lh
      !Latent heat [in calories per gram]
       Lh = 5.95_rkind-0.51_rkind*(t)
      !Abtew [mm/day]
      val = 0.53_rkind*(Rs/Lh)
    end function Abtew_fcn
    
    !#################################################################
    !> 6. Priestley and Taylor 
    !> val: Priestley and Taylor 
    function PriestleyTaylor_fcn (t,Tmonth1,Tmonth2,Rs) result (val)
     use typy
     real (kind=rkind), intent(in) :: Rs,t,Tmonth1,Tmonth2
     real (kind=rkind) :: val
     real(kind=rkind) :: a,sl,PS,Lh,G
     a = 1.26
     !slope of the saturation vapour pressure curve
     sl = (33.8639_rkind*(0.05904_rkind*(0.00738_rkind*(t))+0.8072**7_rkind-0.0000342_rkind))
     !psicromatic constant [kPa/°C]
     Ps = ((1.013e-3)*10313.25_rkind)/(0.622_rkind*Lh)
     !Latent heat [in calories per gram]
     Lh = (5.95_rkind-0.51_rkind*(t))
     !soil heat flux density [MJm2/day]
     
     !Priestley and Taylor [mm/day]
     G = 0.07_rkind*(Tmonth1-Tmonth2)
      val = (a*(sl/(sl+Ps))*((Rs-G)/Lh))
    end function PriestleyTaylor_fcn
    
    
     !#################################################################
    !> 7. Turc Method
    !> val: Turc
    function Turc_fcn (t,Rs) result (val)
      use typy
      !> t: input temperature [°C], Rs: solar radiation [m/day]
      real (kind=rkind), intent(in) :: t,Rs
      real (kind=rkind) :: val
    
    !Turc [mm/day]
     val = 0.013_rkind*(23.88_rkind*(Rs)+ 50_rkind)*(t)*(t+15_rkind)**(-1) 
    end function Turc_fcn
    
    
    !#################################################################
    !> 7.2 Turc1 Method
    !> val: Turc1
    !function Turc1_fcn (t,Rs) result (val)
     !use typy
     !real (kind=rkind), intent(in) :: t,Rs
    !real (kind=rkind) :: val
    !Turc1
    !  val = 0.013_rkind*(t/(t+15_rkind))*(Rs+50_rkind)
    !end function Turc1_fcn
    
    !#################################################################
    !> 7.3 Turc2 Method
    !> val: Turc1
    !function Turc2_fcn (t,Rs,RH) result (val)
     !use typy
     !real (kind=rkind), intent(in) :: t,Rs,RH
     !real (kind=rkind) :: val
    !Turc1
     ! val = 0.013_rkind*((t/(t+15_rkind))*(Rs+50_rkind)*(1+((50_rkind-RH)/70)))
    !end function Turc2_fcn
    
    !#################################################################
    !> 8. Hargreaves Method
    !> val: Hargreaves
    function Hargreaves_fcn (t,Rs) result (val)
     use typy
     real (kind=rkind), intent(in) :: t,Rs
     !Temperature [°C]
     !Solar Radiation [m/day]
     real (kind=rkind) :: val
    !Hargreaves [mm/day]
      val = 0.0135_rkind*(t+17.8_rkind)*Rs
    end function Hargreaves_fcn
    
    !#################################################################
    !> 9. HargreaveSamani Method
    !> val: HargreaveSamani
    function HargreaveSamani_fcn (t,tmax,tmin,Rs) result (val)
     use typy
     real (kind=rkind), intent(in) :: t,tmax,tmin,Rs
     !t = average temperature [°C]
     !tmax = maximum temperature [°C]
     !tmin = minimum temperature [°C]
     !Rs = Solar radiation [m/day]
     real (kind=rkind) :: val
    !HargreaveSamani [mm/day]
      val = 0.0023_rkind*(tmax-tmin)*0.5_rkind*(t+17.8_rkind)*Rs
    end function HargreaveSamani_fcn
    
    !#################################################################
    !> 10. JensenHaise Method
    !> val: JensenHaise
    function JensenHaise_fcn (t,Rs)  result (val)
     use typy
     real (kind=rkind), intent(in) :: t,Rs
     !Temperature [°C]
     !Solar Radiation [m/day]
     real (kind=rkind) :: val
     real(kind=rkind) :: e,Lh,Ch,Ct,tx
     !tx = constant
     tx = -3
     !e = Saturation vapour pressure [kPa]
     e = 102.2
     !Lh = Latent heat [in calories per gram]
     Lh = (5.95_rkind-0.51_rkind*(t))
     !Ch = Coefficient 
     Ch = 50.0_rkind*(t/e)
     !Coefficient temperature 
     Ct = 1.0_rkind/(27_rkind+7.3_rkind*(Ch))
     
    !Jensen_Haise [mm/day]
      val = (Ct*(t-tx)*Rs)/Lh
    end function JensenHaise_fcn
    
    !###############################################################
    ! 11. Thornthwaite method
    !> val: Thornthwaite method
    function Thornthwaite_fcn (t) result (val)
    use typy
    use core_tools
    real (kind=rkind), intent(in) :: t
    !t = average temperature [°C]
    real (kind=rkind) :: val
     real(kind=rkind) :: a_th,I
    !Monthly heat index [-]
    I = (t/5.0_rkind)**1.514
    !Constant [-]
    a_th = 0.49239_rkind + 1792e-5*(I) - 771e-2*(I)**2 + 675e-9*(I)**3
    
    !Thornthwaite
    val = 16_rkind*((10_rkind*t)/I)
    end function Thornthwaite_fcn
    
    
    !###############################################################
    ! 12. Hamon method
    !> val: Hamon method
    function Hamon_fcn (Ld,KPEC,RHOSAT) result (val)
    use typy
    real (kind=rkind), intent(in) :: Ld,KPEC,RHOSAT
    real (kind=rkind) :: val
    !Daytime length [h]
    !Ld = 12.0_rkind
    !Calibration coefficient [-]
    !KPEC = 1.20_rkind
    !Saturated vapor density [g/m3]
    !RHOSAT = 17.30_rkind
    
    !Hamon [mm/day]
    val = 0.1651_rkind*(Ld)*(RHOSAT)*(KPEC)
    end function Hamon_fcn
    
    
    !###############################################################
    ! 13. Linacre method
    !> val: Linacre method
    function Linacre_fcn (z,t,tmin,Ls) result (val)
    use typy
    real (kind=rkind), intent(in) :: z,t,tmin,Ls
    !z = elevation (m)
    !t = average temperature [°C]
    !tmin = minimum temperature [°C]
    !Ls = = latitude of the station [rad]
    real (kind=rkind) :: val
    
    !Linacre [mm/day]
    val = (((500_rkind*(t+0.006_rkind+z))/(100.00_rkind-Ls))+(15.00_rkind*(t-tmin)))/(80-t)
    end function Linacre_fcn
    
    
    !###############################################################
    ! 14. BlaneyCriddle method
    !> val: BlaneyCriddle method
    function BlaneyCriddle_fcn (q,tmin) result (val)
    use typy
    real (kind=rkind), intent(in) :: q,tmin
    real (kind=rkind) :: val
    !q: daily percentage of annual daytime hours
    !tmin = minimum temperarure [°C]
    
    !BlaneyCriddle [mm/day]
    val = (q)*((0.46_rkind)*(tmin) + 8.00_rkind)
    end function BlaneyCriddle_fcn
    
    !###############################################################
    ! 15. Oudin method
    !> val: Oudin method
    function Oudin_fcn (Rs,t) result (val)
    use typy
    real (kind=rkind), intent(in) :: Rs,t
    !Rs = solar radiation [m/day]
    !t = average temperature [°C]
    real (kind=rkind) :: val
    
    !Oudin [mm/day]
    val = 0.408_rkind*Rs*(t+5.00_rkind)
    end function Oudin_fcn
    
      !###############################################################
    ! 16. Romanenko method
    !> val: Romanenko method
    function Romanenko_fcn (RH,t) result (val)
    use typy
    real (kind=rkind), intent(in) :: RH,t
    !RH = realtive humidity [-]
    !t = average temperature [°C]
    real (kind=rkind) :: val
    
    !Romanenko [mm/day]
    val = 0.0018_rkind*(25_rkind + t)**2*(100_rkind - RH)
    end function Romanenko_fcn
    
      !###############################################################
    ! 17. RomanenkOudin method
    !> val: Romanenko_Oudin method
    function RomanenkOudin_fcn (RH,t) result (val)
    use typy
    real (kind=rkind), intent(in) :: RH,t
    real (kind=rkind) :: val
    
    !RomanenkOudin
    val = 4.5_rkind*(1+(t/25_rkind))**2*(1.0_rkind-RH)
    end function RomanenkOudin_fcn
    
    !###############################################################
    ! 18. Kharrufa method
    !> val: Kharrufa method
    function Kharrufa_fcn (p,t) result (val)
    use typy
    real (kind=rkind), intent(in) :: p,t
    !p = % of total daytime hours [-]
    !t = average temperature [°C]
    real (kind=rkind) :: val
    
    !Kharrufa [mm/day]
    val = 0.34_rkind*(p)*(t)**1.3
    end function Kharrufa_fcn


end module re_evap_methods
