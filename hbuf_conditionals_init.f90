subroutine hbuf_conditionals_init(count,trcount)
  use vars, only: ncondavg, condavgname, condavglongname, &
       dowtg_blossey_etal_JAMES2009, use_scam_reference_sounding
  use rad, only: do_output_clearsky_heating_profiles
  use params, only: dodamping, docoriolis, donudging_uv, dosurface
  use grid, only: doShipTrackConditionals
  implicit none

  ! Initialize the list of UW statistics variables written in statistics.f90
  integer count,trcount, n

  call add_to_namelist(count,trcount,'CLDCUMU1', &
         'Upward Cumulative Shaded Cloud Fraction (threshold 20 g/m2)',' ',0)
  call add_to_namelist(count,trcount,'CLDCUMD1', &
         'Downward Cumulative Shaded Cloud Fraction (threshold 20 g/m2)',' ',0)

  call add_to_namelist(count,trcount,'CLDCUMU2', &
         'Upward Cumulative Shaded Cloud Fraction (threshold 0.2 g/m2)',' ',0)
  call add_to_namelist(count,trcount,'CLDCUMD2', &
         'Downward Cumulative Shaded Cloud Fraction (threshold 0.2 g/m2)',' ',0)

  call add_to_namelist(count,trcount,'TBIAS', &
         'Absolute temperature bias (model-OBS)','K',0)
  call add_to_namelist(count,trcount,'QBIAS', &
         'Water vapor mass mixing ratio bias (model-OBS)','g/kg',0)
  call add_to_namelist(count,trcount,'RELHBIAS', &
         'Relative humidity bias (model-OBS)','K',0)

  if(donudging_uv) then
    call add_to_namelist(count,trcount,'UBIAS', &
         'Zonal wind bias (model-OBS)','m/s',0)
    call add_to_namelist(count,trcount,'VBIAS', &
         'Meridional wind bias (model-OBS)','m/s',0)
  end if

  call add_to_namelist(count,trcount,'UXGRID', &
         'Mean Cross-grid flow in zonal direction (U - x translation velocity)','m/s',0)
  call add_to_namelist(count,trcount,'VXGRID', &
         'Mean Cross-grid flow in meridional direction (V - y translation velocity)','m/s',0)

  if(docoriolis) then
    call add_to_namelist(count,trcount,'UGEOSTR', &
         'Geostrophic wind in zonal direction','m/s',0)
    call add_to_namelist(count,trcount,'VGEOSTR', &
         'Geostrophic wind in meridional direction','m/s',0)
  end if

  if(use_scam_reference_sounding) then
    call add_to_namelist(count,trcount,'TABSREF', &
         'Referenxe absolute temperature sounding','K',0)
    call add_to_namelist(count,trcount,'QVREF', &
         'Reference water vapor mass mixing ratio','g/kg',0)
  end if

  if(dodamping) then
    call add_to_namelist(count,trcount,'TKEDAMP', &
         'Daming of TKE by sponge region at domain top','m2/s3',0)
  end if

  if(do_output_clearsky_heating_profiles) then
    call add_to_namelist(count,trcount,'RADQRCLW', &
         'Clearsky longwave heating rate','K/d',0)
    call add_to_namelist(count,trcount,'RADQRCSW', &
         'Clearsky shortwave heating rate','K/d',0)
  end if

  if(dowtg_blossey_etal_JAMES2009) then
    call add_to_namelist(count,trcount,'WWTG', &
         'Large-scale W induced by weak temperature gradient approx','m/s',0)
    call add_to_namelist(count,trcount,'TVPR_WTG', &
         'Virtual temperature anomaly wrt reference sounding, used to drive WWTG','K',0)
  end if

  if(dowtg_blossey_etal_JAMES2009) then
    call add_to_namelist(count,trcount,'WOBSREF', &
         'Reference Large-scale W Before Modifications by WTG/Scaling','m/s',0)
  end if

  if(dosurface) then
    call add_to_namelist(count,trcount,'SURFVARS', &
         'NOT A PROFILE (suface vars): wspd, ustar, taux, tauy, SHF/rho/Cp, LHF/rho/L, Ct*wspd, Cq*wspd, qsat_surf, tskin, SST','',0)
  end if

  !bloss: setup to add an arbitrary number of conditional statistics
  do n = 1,ncondavg

     !bloss: add all of the conditional statistics here, so that they don't
     !  have to be added to the lst file
     call add_to_namelist(count,trcount,TRIM(condavgname(n)), &
          TRIM(condavglongname(n))//' Fraction',' ',0)
     call add_to_namelist(count,trcount,'W'//TRIM(condavgname(n)), &
          'Mean W in '//TRIM(condavglongname(n)),'m/s',n)
     call add_to_namelist(count,trcount,'U'//TRIM(condavgname(n)), &
          'Mean U in '//TRIM(condavglongname(n)),'m/s',n)
     call add_to_namelist(count,trcount,'V'//TRIM(condavgname(n)), &
          'Mean V in '//TRIM(condavglongname(n)),'m/s',n)
     
     if(doShipTrackConditionals) then
     call add_to_namelist(count,trcount,'ZB'//TRIM(condavgname(n)), &
          'Inversion base height in '//TRIM(condavglongname(n)),'m',n)
     call add_to_namelist(count,trcount,'INVT'//TRIM(condavgname(n)), &
          'Inversion Tend. in '//TRIM(condavglongname(n)),'m/s',n)
     call add_to_namelist(count,trcount,'WEM'//TRIM(condavgname(n)), &
          'Entrain. Rate in '//TRIM(condavglongname(n)),'m/s',n)
     call add_to_namelist(count,trcount,'WLO'//TRIM(condavgname(n)), &
          'Local vert. vel. in '//TRIM(condavglongname(n)),'m/s',n)
     call add_to_namelist(count,trcount,'WLA'//TRIM(condavgname(n)), &
          'Large-scale vert. vel. in '//TRIM(condavglongname(n)),'m/s',n)
     call add_to_namelist(count,trcount,'RDIV'//TRIM(condavgname(n)), &
          'Rad. Div. in '//TRIM(condavglongname(n)),'W/m2',n)
     call add_to_namelist(count,trcount,'LDIV'//TRIM(condavgname(n)), &
          'LW Rad. Div. in '//TRIM(condavglongname(n)),'W/m2',n)
     call add_to_namelist(count,trcount,'SDIV'//TRIM(condavgname(n)), &
          'SW Rad. Div. in '//TRIM(condavglongname(n)),'W/m2',n)
     call add_to_namelist(count,trcount,'SH_F'//TRIM(condavgname(n)), &
          'Sfc. SHF in '//TRIM(condavglongname(n)),'W/m2',n)
     call add_to_namelist(count,trcount,'LH_F'//TRIM(condavgname(n)), &
          'Sfc. LHF in '//TRIM(condavglongname(n)),'W/m2',n)
     call add_to_namelist(count,trcount,'PR_S'//TRIM(condavgname(n)), &
          'Sfc. Precip. rate in '//TRIM(condavglongname(n)),'m/s',n)
     call add_to_namelist(count,trcount,'PR_Z'//TRIM(condavgname(n)), &
          'Inv. Base Precip. rate in '//TRIM(condavglongname(n)),'m/s',n)
     call add_to_namelist(count,trcount,'TL_T'//TRIM(condavgname(n)), &
          'TL BL-avg. tendency in '//TRIM(condavglongname(n)),'K/s',n)
     call add_to_namelist(count,trcount,'QT_T'//TRIM(condavgname(n)), &
          'QT BL-avg. tendency in '//TRIM(condavglongname(n)),'1/s',n)
     call add_to_namelist(count,trcount,'CL_T'//TRIM(condavgname(n)), &
          'Cloud Thick. tendency in '//TRIM(condavglongname(n)),'m/s',n)
     call add_to_namelist(count,trcount,'CB_T'//TRIM(condavgname(n)), &
          'Cloud Base tendency in '//TRIM(condavglongname(n)),'m/s',n)
     call add_to_namelist(count,trcount,'LWPT'//TRIM(condavgname(n)), &
          'LWP tendency in '//TRIM(condavglongname(n)),'kg/kg/m2/s',n)
     call add_to_namelist(count,trcount,'SW_U'//TRIM(condavgname(n)), &
          'Upwelling Shortwave in '//TRIM(condavglongname(n)),'W/m2',n)
     call add_to_namelist(count,trcount,'SW_D'//TRIM(condavgname(n)), &
          'Downwelling Shortwave in '//TRIM(condavglongname(n)),'W/m2',n)
     call add_to_namelist(count,trcount,'LW_U'//TRIM(condavgname(n)), &
          'Upwelling longwave in '//TRIM(condavglongname(n)),'W/m2',n)
     call add_to_namelist(count,trcount,'LW_D'//TRIM(condavgname(n)), &
          'Downwelling longwave in '//TRIM(condavglongname(n)),'W/m2',n)
     call add_to_namelist(count,trcount,'TS2'//TRIM(condavgname(n)), &
          'TL variance in '//TRIM(condavglongname(n)),'K2',n)
     call add_to_namelist(count,trcount,'QS2'//TRIM(condavgname(n)), &
          'QT variance in '//TRIM(condavglongname(n)),'',n)
     call add_to_namelist(count,trcount,'WS2'//TRIM(condavgname(n)), &
          'W variance in '//TRIM(condavglongname(n)),'m2/s2',n)
     call add_to_namelist(count,trcount,'US2'//TRIM(condavgname(n)), &
          'U variance in '//TRIM(condavglongname(n)),'m2/s2',n)
     call add_to_namelist(count,trcount,'VS2'//TRIM(condavgname(n)), &
          'V variance in '//TRIM(condavglongname(n)),'m2/s2',n)
     call add_to_namelist(count,trcount,'WTL'//TRIM(condavgname(n)), &
          'Vertical sensible heat flux in '//TRIM(condavglongname(n)),'W/m2',n)
     call add_to_namelist(count,trcount,'WQT'//TRIM(condavgname(n)), &
          'Vertical latent heat flux in '//TRIM(condavglongname(n)),'W/m2',n)
     call add_to_namelist(count,trcount,'WTV'//TRIM(condavgname(n)), &
          'Vertical buoyancy flux in '//TRIM(condavglongname(n)),'W/m2',n)
     call add_to_namelist(count,trcount,'PFS'//TRIM(condavgname(n)), &
          'Precipitation rate in '//TRIM(condavglongname(n)),'m/s',n)
     endif
     
     call add_to_namelist(count,trcount,'MSE'//TRIM(condavgname(n)), &
          'Mean moist static energy in '//TRIM(condavglongname(n)),'K',n)
     call add_to_namelist(count,trcount,'DSE'//TRIM(condavgname(n)), &
          'Mean dry static energy in '//TRIM(condavglongname(n)),'K',n)
     call add_to_namelist(count,trcount,'TL'//TRIM(condavgname(n)), &
          'Mean liquid-ice static energy in '//TRIM(condavglongname(n)),'K',n)
     call add_to_namelist(count,trcount,'TA'//TRIM(condavgname(n)), &
          'Mean TABS in '//TRIM(condavglongname(n)),'K',n)
     call add_to_namelist(count,trcount,'TV'//TRIM(condavgname(n)), &
          'Mean THETAV in '//TRIM(condavglongname(n)),'K',n)
     call add_to_namelist(count,trcount,'TV'//TRIM(condavgname(n))//'A', &
          'Mean THETAV anomaly in '//TRIM(condavglongname(n)),'K',n)
     call add_to_namelist(count,trcount,'QT'//TRIM(condavgname(n)), &
          'Mean QT in '//TRIM(condavglongname(n)),'g/kg',n)
     call add_to_namelist(count,trcount,'QN'//TRIM(condavgname(n)), &
          'Mean QN in '//TRIM(condavglongname(n)),'g/kg',n)
     !bloss: these conditional averages are now computed inside the microphysics
     !         routines.
     !bloss        call add_to_namelist(count,trcount,'QC'//TRIM(condavgname(n)), &
     !bloss             'Mean QC in '//TRIM(condavglongname(n)),'g/kg',n)
     !bloss        call add_to_namelist(count,trcount,'QI'//TRIM(condavgname(n)), &
     !bloss             'Mean QI in '//TRIM(condavglongname(n)),'g/kg',n)
     call add_to_namelist(count,trcount,'QP'//TRIM(condavgname(n)), &
          'Mean QP in '//TRIM(condavglongname(n)),'g/kg',n)
     call add_to_namelist(count,trcount,'S'//TRIM(condavgname(n)), &
          'Mean scalar in '//TRIM(condavglongname(n)),'K',n)
     call add_to_namelist(count,trcount,'W'//TRIM(condavgname(n))//'A', &
          'W in '//TRIM(condavglongname(n))//' averaged over the whole domain','m/s',0)
     call add_to_namelist(count,trcount,'TLW'//TRIM(condavgname(n)), &
          'TLW in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'Km/s',0)
     call add_to_namelist(count,trcount,'TVW'//TRIM(condavgname(n)), &
          'TVW in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'Km/s',0)
     call add_to_namelist(count,trcount,'SW'//TRIM(condavgname(n)), &
          'SW in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'Km/s',0)
     call add_to_namelist(count,trcount,'QTW'//TRIM(condavgname(n)), &
          'QTW in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'g/kg m/s',0)
     call add_to_namelist(count,trcount,'QCW'//TRIM(condavgname(n)), &
          'QCW in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'g/kg m/s',0)
     call add_to_namelist(count,trcount,'QIW'//TRIM(condavgname(n)), &
          'QIW in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'g/kg m/s',0)

     !bloss: frozen moist static energy statistics
     call add_to_namelist(count,trcount,'HF'//TRIM(condavgname(n)), &
          'Mean Frozen MSE in '//TRIM(condavglongname(n)),'K',n)
     call add_to_namelist(count,trcount,'HF'//TRIM(condavgname(n))//'A', &
          'Mean Frozen MSE anomaly in '//TRIM(condavglongname(n)),'K',n)

     !bloss: velocity anomalies
     call add_to_namelist(count,trcount,'U'//TRIM(condavgname(n))//'A', &
          'Mean U anomaly in '//TRIM(condavglongname(n)),'m/s',n)
     call add_to_namelist(count,trcount,'V'//TRIM(condavgname(n))//'A', &
          'Mean V anomaly in '//TRIM(condavglongname(n)),'m/s',n)

     !bloss: pressure gradients
     call add_to_namelist(count,trcount,'UPGF'//TRIM(condavgname(n)), &
          'Zonal pressure gradient in '//TRIM(condavglongname(n)),'m/s2',n)
     call add_to_namelist(count,trcount,'VPGF'//TRIM(condavgname(n)), &
          'Meridional pressure gradient in '//TRIM(condavglongname(n)),'m/s2',n)
     call add_to_namelist(count,trcount,'WPGF'//TRIM(condavgname(n)), &
          'Vertical pressure gradient in '//TRIM(condavglongname(n)),'m/s2',n)

     !bloss: momentum statistics
     call add_to_namelist(count,trcount,'UW'//TRIM(condavgname(n)), &
          'UW in '//TRIM(condavglongname(n)),'m2/s2',n)
     call add_to_namelist(count,trcount,'VW'//TRIM(condavgname(n)), &
          'VW in '//TRIM(condavglongname(n)),'m2/s2',n)
     call add_to_namelist(count,trcount,'UWSB'//TRIM(condavgname(n)), &
          'Subgrid UW in '//TRIM(condavglongname(n)),'m2/s2',n)
     call add_to_namelist(count,trcount,'VWSB'//TRIM(condavgname(n)), &
          'Subgrid VW in '//TRIM(condavglongname(n)),'m2/s2',n)

     !bloss: UW-added mass flux weighted statistics
     call add_to_namelist(count,trcount,'MF'//TRIM(condavgname(n)), &
          'Mass flux in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'kg/m2/s',0)
     call add_to_namelist(count,trcount,'MFH'//TRIM(condavgname(n)), &
          'RHO*W*HF in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'K kg/m2/s',0)
     call add_to_namelist(count,trcount,'MFH'//TRIM(condavgname(n))//'A', &
          'RHO*W*HF anomaly in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'K kg/m2/s',0)
     call add_to_namelist(count,trcount,'MFTL'//TRIM(condavgname(n)), &
          'RHO*W*TL in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'K kg/m2/s',0)
     call add_to_namelist(count,trcount,'MFTL'//TRIM(condavgname(n))//'A', &
          'RHO*W*TL anomaly in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'K kg/m2/s',0)
     call add_to_namelist(count,trcount,'MFTV'//TRIM(condavgname(n)), &
          'RHO*W*TV in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'K kg/m2/s',0)
     call add_to_namelist(count,trcount,'MFTV'//TRIM(condavgname(n))//'A', &
          'RHO*W*TV anomaly in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'K kg/m2/s',0)
     call add_to_namelist(count,trcount,'MFQT'//TRIM(condavgname(n)), &
          'RHO*W*QT in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'g/m2/s',0)
     call add_to_namelist(count,trcount,'MFQT'//TRIM(condavgname(n))//'A', &
          'RHO*W*QT anomaly in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'g/m2/s',0)
     call add_to_namelist(count,trcount,'RUW'//TRIM(condavgname(n)), &
          'RHOUW in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'kg/m/s2',0)
     call add_to_namelist(count,trcount,'RVW'//TRIM(condavgname(n)), &
          'RHOVW in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'kg/m/s2',0)
     call add_to_namelist(count,trcount,'RWW'//TRIM(condavgname(n)), &
          'RHOWW in '//TRIM(condavglongname(n))//' averaged over the whole domain', &
          'kg/m/s2',0)
  end do ! n = 1,ncondavg

end

subroutine add_to_namelist(count,trcount,varname,varlongname,varunits,varavg)
  use hbuffer, only: namelist,deflist,unitlist,status,average_type
  implicit none

  ! add variable to namelist
  integer count, trcount, ntr, n, varstatus, varavg
  character(*) varname
  character(*) varlongname
  character(*) varunits

  count = count + 1
  trcount = trcount + 1
  namelist(count) = trim(varname)
  deflist(count) = trim(varlongname)
  unitlist(count) = trim(varunits)
  status(count) = 1
  average_type(count) = varavg

end subroutine add_to_namelist
