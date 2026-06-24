pro Example_ep
 Nf=100 ;number of frequencies
 L=1d10 ;source depth, cm
 S=1d20 ;source area 

 libname='GSmulti_64.dll'
 
 ;--------------------------------------------------------------------------
 ;a) Compute the emission from electrons and positrons with equal parameters
 ;   Particle distributions: isotropic power-law
 
 Lparms=lonarr(4)
 Lparms[0]=Nf
 Lparms[1]=2 ;number of components
 Lparms[3]=1 ;FF is off
 
 Rparms=dblarr(3)
 Rparms[0]=S    ;area, cm^2
 Rparms[1]=1d9  ;starting frequency to calculate spectrum, Hz
 Rparms[2]=0.02 ;logarithmic step in frequency
 
 Sparms=dblarr(5)
 Sparms[0]=L   ;voxel depth, cm
 Sparms[1]=3d9 ;n_th - thermal electron density, cm^{-3}
 Sparms[2]=3d7 ;T_0, K
 Sparms[3]=180 ;B - magnetic field, G
 Sparms[4]=75  ;theta - viewing angle, /degrees
 
 Bparms=dblarr(14, 2)
 Bparms[0, 0]=1    ;q/e - electrons
 Bparms[1, 0]=1    ;m/e - electrons
 Bparms[2, 0]=1d6  ;n_b
 Bparms[3, 0]=3    ;energy distribution (PLW)
 Bparms[4, 0]=0.1  ;E_min
 Bparms[5, 0]=10.0 ;E_max
 Bparms[7, 0]=4    ;\delta_1
 Bparms[13, 0]=1   ;add n_b to n_th to compute the dispersion
 
 Bparms[0, 1]=-1   ;q/e - positrons
 Bparms[1, 1]=1    ;m/e - positrons
 Bparms[2, 1]=1d6  ;n_b
 Bparms[3, 1]=3    ;energy distribution (PLW)
 Bparms[4, 1]=0.1  ;E_min
 Bparms[5, 1]=10.0 ;E_max
 Bparms[7, 1]=4    ;\delta_1 
 
 Flux=dblarr(3, Nf)
 
 res=call_external(libname, 'GET_MW', Lparms, Rparms, Sparms, Bparms, Flux, /unload)
 
 f=reform(Flux[0, *])
 IXep=reform(Flux[1, *])
 IOep=reform(Flux[2, *])

 ;--------------------------------------------------------------------------
 ;b) Compute the emission from electrons only
 
 Lparms[1]=1 ;number of components
 
 Bparms=dblarr(14, 1)
 Bparms[0, 0]=1    ;q/e - electrons
 Bparms[1, 0]=1    ;m/e - electrons
 Bparms[2, 0]=1d6  ;n_b
 Bparms[3, 0]=3    ;energy distribution (PLW)
 Bparms[4, 0]=0.1  ;E_min
 Bparms[5, 0]=10.0 ;E_max
 Bparms[7, 0]=4    ;\delta_1
 Bparms[13, 0]=1   ;add n_b to n_th to compute the dispersion
 
 ;other parameters are the same as above
 
 res=call_external(libname, 'GET_MW', Lparms, Rparms, Sparms, Bparms, Flux, /unload)
 
 IXe=reform(Flux[1, *])
 IOe=reform(Flux[2, *])
 
 ;--------------------------------------------------------------------------
 ;c) Compute the emission from positrons only
 
 Lparms[1]=1 ;number of components
 
 Bparms=dblarr(14, 1)
 Bparms[0, 0]=-1   ;q/e - positrons
 Bparms[1, 0]=1    ;m/e - positrons
 Bparms[2, 0]=1d6  ;n_b
 Bparms[3, 0]=3    ;energy distribution (PLW)
 Bparms[4, 0]=0.1  ;E_min
 Bparms[5, 0]=10.0 ;E_max
 Bparms[7, 0]=4    ;\delta_1
 
 ;other parameters are the same as above
 
 res=call_external(libname, 'GET_MW', Lparms, Rparms, Sparms, Bparms, Flux, /unload)
 
 IXp=reform(Flux[1, *])
 IOp=reform(Flux[2, *])
 
 ;--------------------------------------------
 ;Visualize the results
 
 window, 1, title='Intensity, X+O'
 wset, 1
 plot, f, IXep+IOep, /xlog, /ylog, xtitle='Frequency, GHz', ytitle='Intensity, sfu', thick=2, linestyle=2
       ;dashed line: electrons + positrons
 oplot, f, IXe+IOe, thick=2 ;solid line: electrons
 oplot, f, IXp+IOp, thick=2, linestyle=1 ;dotted line: positrons
    
 window, 2, title='Polarization, (X-O)/(X+O)'
 wset, 2
 plot, f, (IXep-IOep)/(IXep+IOep), /xlog, yrange=[-0.2, 0.2], xtitle='Frequency, GHz', ytitle='Polarization degree', $
       thick=2, linestyle=2 ;dashed line: electrons + positrons
 oplot, f, (IXe-IOe)/(IXe+IOe), thick=2 ;solid line: electrons
 oplot, f, (IXp-IOp)/(IXp+IOp), thick=2, linestyle=1 ;dotted line: positrons
 plots, [1, 100], [0, 0], linestyle=5, /data
end