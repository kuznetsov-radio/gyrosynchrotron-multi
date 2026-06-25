pro gsProtons
 Nf=160    ;number of frequencies
 L=725d5*15    ;total source depth, cm
 S=(725d5*20)^2    ;area, cm^2
 
 libname='GSmulti_64.dll'
 
 Lparms=lonarr(4)
 Lparms[0]=Nf
 Lparms[1]=1 ;number of components
 Lparms[2]=32 ;N_nodes
 Lparms[3]=1 ;FF off
 
 Rparms=dblarr(3)
 Rparms[0]=S ;area, cm^2
 Rparms[1]=3d8  ;starting frequency to calculate spectrum, Hz
 Rparms[2]=0.02 ;logarithmic step in frequency
 
 Sparms=dblarr(5)
 Sparms[0]=L   ;voxel depth, cm
 Sparms[1]=1d8 ;n_0 - thermal electron density, cm^{-3}
 Sparms[2]=2d7 ;T_0, K
 Sparms[3]=1500 ;B - magnetic field, G
 Sparms[4]=80  ;theta - viewing angle, /degrees
 
 Bparms=dblarr(14, 1)
 Bparms[0, 0]=-1    ;q/e - protons
 Bparms[1, 0]=1836  ;m/e - protons
 Bparms[2, 0]=3d7  ;n_b
 Bparms[3, 0]=3    ;energy distribution (PLW)
 Bparms[4, 0]=2  ;E_min
 Bparms[5, 0]=7d4 ;E_max
 Bparms[6, 0]=1.0  ;E_break
 Bparms[7, 0]=1.5    ;\delta_1
 Bparms[8, 0]=6    ;\delta_2
 Bparms[9, 0]=0    ;pitch-angle distribution (ISO)
 Bparms[10, 0]=60  ;\alpha_c
 Bparms[11, 0]=0.1 ;\Delta\mu
 
 Flux=dblarr(3, Nf)
 
 res=call_external(libname, 'GET_MW', Lparms, Rparms, Sparms, Bparms, Flux, /unload)
 
 f=reform(Flux[0, *])
 IX=reform(Flux[1, *])
 IO=reform(Flux[2, *])
 
 ;----------------------------------------------------
 
 fp=8.98d3*sqrt(Sparms[1])
 fB=2.8d6*Sparms[3]
 fc=fB/2+sqrt(fp^2+fB^2/4)
 
 forward_function expnumber
 
 set_plot, 'ps'
 device, file='Spectrum_p+.eps', xsize=8.9, ysize=7, font_size=10, /encapsulated, /color
 !X.MARGIN=[7.5, 1.5]
 !Y.MARGIN=[3.2, 0.5]
 tvlct, [0,     255, 0,     0,    196,     0,    196, 196, 128], $
        [0,     0,   128,   0,    0,       196,  196, 196, 128], $
        [0,     0,   0,     255,  196,     196,  0,   196, 128] 

 plot, f, IX, /xlog, /ylog, xtitle='!6Frequency, GHz!3', ytitle='!6Intensity, sfu!3', /nodata, $
       yrange=[1d-4, 1d2], xrange=[0.3, 300], xstyle=1, ytickformat='expnumber'
 plots, fc/1d9*[1, 1], [1d-4, 1d2], /data, linestyle=2
 oplot, f, IX, color=1, thick=3
 oplot, f, IO, color=3, thick=3
 
 xyouts, 1830, 6250, '!6p!U+!N!3', /device
 xyouts, 5450, 6430, '!6X-mode!3', /device, color=1
 xyouts, 1700, 4600, '!6O-mode!3', /device, color=3
 xyouts, 200, 200, '!6(a)!3', /device
 
 device, /close
end

function expnumber, axis, index, value
 u=floor(alog10(value))
 return, '10!U'+string(u, format='(2I0)')+'!N'
end