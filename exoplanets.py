# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as scicon
import scipy.integrate as integrate
import pandas as pd
import seaborn as sb
import scipy.special as sp

import star_planet_data as data
 

#https://orbi.uliege.be/bitstream/2268/25280/1/Absil06_thesis.pdf page 151 of pdf
#confirmed
def phase_fct_refl():
    return 0.5

#https://orbi.uliege.be/bitstream/2268/25280/1/Absil06_thesis.pdf page 151 of pdf
#confirmed
def phase_fct_thermal():
    return 1
        
#confirmed
def var_phase_error(wl,sigma_OPD):
    phi_squared=(2*np.pi*sigma_OPD/wl)**2
    
    return phi_squared

#ToDo Check formula
def var_frac_intensity(wl,sigma_TipTilt):
    I_squared=np.sqrt(8)*(np.pi*0.25*sigma_TipTilt)**4/(64*wl**4)
    
    return I_squared


def inst_throughput(wl):
    wl_um=wl*10**6
    factor=-2/295*wl_um+297/590
    
    return factor


def blackbody(T,wl):
    B=2*scicon.h*scicon.c**2/wl**5/(np.exp(scicon.h*scicon.c/(scicon.k*wl*T))-1)
    
    return B


def local_zodiacal_emission_per_wl(star,wl):
    I=4*10**(-8)*(blackbody(265,wl)+NIR_dust_albedo*blackbody(data.T_Sun,wl)*(data.R_Sun/(1.5*data.AU))**2)*np.sqrt(np.pi/(np.arccos(np.cos(np.pi)*np.cos(star.lat)))/(np.sin(star.lat)**2+0.36*(wl/(11*10**(-6)))**(-0.8)*np.cos(star.lat)**2))

    return I


def local_zodiacal_emission(star,wl,wl_band,telescope):
    #integrate over wavelength-band
    zodi_W_perArea=integrate.quad(lambda x: local_zodiacal_emission_per_wl(star,x),wl-wl_band/2,wl+wl_band/2)[0]
    #multiply by area of both apertures
    zodi_W_apertures=zodi_W_perArea*2*(telescope.D/2)**2*np.pi
    #multiply by solid angle
    #matlab 0.61
    zodi_solid_angle=(0.514*wl/telescope.D)**2*np.pi*zodi_W_apertures
    #convert from W to ph/s
    zodi_photons_per_second=zodi_solid_angle/(scicon.h*scicon.c/wl)
    
    return zodi_photons_per_second
    

def exozodiacal_emission(star,wl):
    return 0


def thermal_emission_instrument_per_wl(wl,telescope):
    I=grey_body_emissivity*blackbody(telescope.temp,wl)
    
    return I


def thermal_emission_instrument(wl,wl_band,telescope):
    thermal_W_perArea=integrate.quad(lambda x: thermal_emission_instrument_per_wl(x,telescope),wl-wl_band/2,wl+wl_band/2)[0]
    #ToDo Implement telescope
    thermal_W_total=thermal_W_perArea*instrument_thermal_multiplicator
    thermal_photons_per_second=thermal_W_total/(scicon.h*scicon.c/wl)
    
    return thermal_photons_per_second


def output_planet_flux_per_wl(star,planet,wl,telescope):
    #matlab doesn't have /4 at end
    reflection=planet.albedo*np.pi*star.radius**2/star.dist**2*blackbody(star.temp,wl)*phase_fct_refl()*planet.radius**2/(1*planet.dist**2)
    thermal_emission=np.pi*planet.radius**2/star.dist**2*phase_fct_thermal()*blackbody(planet.temp,wl)
    Transmission=np.sin(np.pi*np.tan(planet.dist/star.dist)*telescope.b/wl)**2
    total_output=(reflection+thermal_emission)*Transmission
    
    return total_output


def output_planet_flux(star,planet,wl,wl_band,telescope):
    #integrate over wavelength-band
    OPF_W_perArea=integrate.quad(lambda x: output_planet_flux_per_wl(star,planet,x,telescope),wl-wl_band/2,wl+wl_band/2)[0]
    #multiply by projected solid angle
    #maybe 2*pi? matlab doesn't have it at all
    OPF_solid_angle=OPF_W_perArea#*np.pi
    #multiply by area of both apertures
    OPF_W_apertures=OPF_solid_angle*2*(telescope.D/2)**2*np.pi
    #convert from W to ph/s
    OPF_photons_per_second=OPF_W_apertures/(scicon.h*scicon.c/wl)
    #multiply by intrumental throughput
    OPF_throughput=OPF_photons_per_second*inst_throughput(wl)
    
    return OPF_throughput


def output_background_flux(star,planet,wl,wl_band,telescope):
    inst_thermal=thermal_emission_instrument(wl,wl_band,telescope)
    local_zodi=local_zodiacal_emission(star,wl,wl_band,telescope)
    exozodi=exozodiacal_emission(star,wl)
    OBF_throughput=(local_zodi+exozodi+inst_thermal)*inst_throughput(wl)
    
    return OBF_throughput


def output_stellar_flux_per_wl(star,wl,telescope,sigma_OPD,sigma_TipTilt):
    F=np.pi*star.radius**2/star.dist**2*blackbody(star.temp,wl)
    N_time_average=1/4*(var_phase_error(wl,sigma_OPD)+np.pi**2/4*(star.solid_angle/((wl)/telescope.b))**2+var_frac_intensity(wl,sigma_TipTilt))
    F_remaining=F*N_time_average
    
    return F_remaining


def output_stellar_flux(star,wl,wl_band,telescope,sigma_OPD,sigma_TipTilt):
    #without shortest wl
    F_W_perArea=integrate.quad(lambda x: output_stellar_flux_per_wl(star,x,telescope,sigma_OPD,sigma_TipTilt),wl-wl_band/2,wl+wl_band/2)[0]
    #maybe2*pi? matlab doesn't have it at all
    F_W_solid_angle=F_W_perArea#*np.pi
    F_W_apertures=F_W_solid_angle*2*(telescope.D/2)**2*np.pi
    F_photons_per_second=F_W_apertures/(scicon.h*scicon.c/wl)
    OSF_throughput=F_photons_per_second*inst_throughput(wl)
    
    return OSF_throughput

'''
def output_stellar_flux(star,wl,wl_band,telescope,sigma_OPD,sigma_TipTilt):
    #old version taking into account shortest wavelength
    N_time_average=1/4*(var_phase_error(wl,sigma_OPD)+np.pi**2/4*(star.solid_angle/((wl-wl_band/2)/telescope.b))**2+var_frac_intensity(wl,sigma_TipTilt))
    F_W_perArea=integrate.quad(lambda x: output_stellar_flux_per_wl(star,x),wl-wl_band/2,wl+wl_band/2)[0]
    F_W_solid_angle=F_W_perArea*np.pi
    F_W_apertures=F_W_solid_angle*2*(telescope.D/2)**2*np.pi
    F_photons_per_second=F_W_apertures/(scicon.h*scicon.c/wl)
    OSF_throughput=F_photons_per_second*N_time_average*inst_throughput(wl)
    
    return OSF_throughput
'''

def optimal_wl(star,planet,telescope,sigma_OPD,sigma_TipTilt):
    wl_array=np.linspace(0.5*10**-6,30*10**-6,100)

    pl_signal=np.empty_like(wl_array)
    for i in range(wl_array.size):
        pl_signal[i]=output_planet_flux(star,planet,wl_array[i],wl_array[i]/R,telescope)
    
    shot_n=np.empty_like(wl_array)
    for i in range(wl_array.size):
        shot_n[i]=np.sqrt(output_stellar_flux(star,wl_array[i],wl_array[i]/R,telescope,sigma_OPD,sigma_TipTilt)+pl_signal[i]+output_background_flux(star,planet,wl_array[i],wl_array[i]/R,telescope))
        
    instrumental_n=np.empty_like(wl_array)
    for i in range(wl_array.size):
        instrumental_n[i]=output_stellar_flux(star,wl_array[i],wl_array[i]/R,telescope,sigma_OPD,sigma_TipTilt)*np.sqrt((var_phase_error(wl_array[i],sigma_OPD)**2+var_frac_intensity(wl_array[i],sigma_TipTilt)**2)/8)
    
    SNR_array=np.empty_like(wl_array)
    for i in range(wl_array.size):
        SNR_array[i]=pl_signal[i]/np.sqrt(shot_n[i]**2+instrumental_n[i]**2)
        
    optimal_wavelength=wl_array[np.argmax(SNR_array)]
    
    return optimal_wavelength


def adjust_baseline(star,planet,telescope,wl):
    angular_sep=np.tan(planet.dist/star.dist)
    if(telescope.b>=wl/(2*angular_sep)):
        adjusted_bl=wl/(2*angular_sep)
    else:
        adjusted_bl=telescope.b
    
    return adjusted_bl

        
def integration_time(star,planet,wl,telescope,sigma_OPD,sigma_TipTilt,show_details,convert_to_sec,adjust_bl):
    
    if(wl=="optimal"):
        #max function to make sure you do not go outside the OWA of the telescope
        wl=np.max(np.array([optimal_wl(star,planet,telescope,sigma_OPD,sigma_TipTilt),telescope.D*np.tan(planet.dist/star.dist)/0.514]))
    
    wl_band=wl/R
    original_bl=telescope.b
    if(adjust_bl==True):
        adjusted_bl=adjust_baseline(star,planet,telescope,wl)
    else:
        adjusted_bl=original_bl
    telescope.adjust_baseline(adjusted_bl)
    
    #fluxes in ph/s
    OPF=output_planet_flux(star,planet,wl,wl_band,telescope)
    OBF=output_background_flux(star,planet,wl,wl_band,telescope)
    OSF=output_stellar_flux(star,wl,wl_band,telescope,sigma_OPD,sigma_TipTilt)
    N_s=np.sqrt(OPF+OBF+OSF)
    N_inst=OSF*np.sqrt((var_phase_error(wl,sigma_OPD)**2+var_frac_intensity(wl,sigma_TipTilt)**2)/8)
    SNR_1s=OPF/np.sqrt(N_s**2+N_inst**2)
    t=(SNR_req/SNR_1s)**2/3600
    if(convert_to_sec==True):
        t=t*3600
    if(show_details==True):
        print(telescope.name)
        print(planet.name,", ",star.name)
        print("adjusted telescope baseline:",telescope.b)
        print("OPF:",OPF)
        print("OBF:",OBF)
        print("OSF:",OSF)
        print("Shot noise:",N_s)
        print("Instrumental noise:",N_inst)
        print("SNR_1s:",SNR_1s)
        if(convert_to_sec==True):
            print("integration time:",t,"seconds at wavelength",np.round(wl*10**6,2),"\u03BCm")
        else:
            print("integration time:",t,"hours at wavelength",np.round(wl*10**6,2),"\u03BCm")
        print("-----------------------------------------------")
    telescope.adjust_baseline(original_bl)
    
    return t, adjusted_bl


# constants
SNR_req=5
R=1.2
grey_body_emissivity=0.25
NIR_dust_albedo=0.22
sigma_OPD_zero=0
sigma_TipTilt_zero=0
sigma_OPD_max=10*10**-9
sigma_TipTilt_max=300*np.pi/(180*3600)*10**-3
sigma_OPD_test=5*10**-9
sigma_TipTilt_test=50*np.pi/(180*3600)*10**-3

#ToDo
instrument_thermal_multiplicator=10**-10.5


#main function for star-planet system
def main(star,planet,wl,sigma_OPD,sigma_TipTilt,show_details,convert_to_sec,adjust_bl):
    angular_sep_planet=np.tan(planet.dist/star.dist)
    telescope_list=[data.CubeSat6U,data.CubeSat12U,data.PROBA]
    time_list=[]
    adjusted_bl_list=[]
    for i in range(len(telescope_list)):
        if(wl=="optimal"):
            used_wl=telescope_list[i].wl
        else:
            used_wl=wl
        if(angular_sep_planet>0.514*used_wl/telescope_list[i].D):
        #if(angular_sep_planet>telescope_list[i].OWA):
            if(show_details==True):
                print(planet.name,"outside OWA of",telescope_list[i].name)
            time_list.append("inf")
            adjusted_bl_list.append("n/a")
    
        else:
            t,bl=integration_time(star,planet,wl,telescope_list[i],sigma_OPD,sigma_TipTilt,show_details,convert_to_sec,adjust_bl)
            time_list.append(t)
            adjusted_bl_list.append(bl)
    
    return time_list,adjusted_bl_list



#final results; all directly imaged planets for optimized baseline and desired wavelength, sigma_OPD and sigma_TipTilt
def get_results(wl,sigma_OPD,sigma_TipTilt,save_fig,convert_to_sec,adjust_bl):
    directly_imaged_table=pd.DataFrame(columns=["CubeSat6U","CubeSat12U","PROBA"])
    for i in range(len(data.directly_imaged_systems_planets)):
        time_list=main(data.directly_imaged_systems_stars[i],data.directly_imaged_systems_planets[i],wl,sigma_OPD,sigma_TipTilt,False,convert_to_sec,adjust_bl)[0]
        new_row=pd.DataFrame([time_list],index=[data.directly_imaged_systems_planets[i].name],columns=["CubeSat6U","CubeSat12U","PROBA"])
        directly_imaged_table=directly_imaged_table.append(new_row)
    
    if(save_fig==True):
        directly_imaged_table.to_csv("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_files\\directly_imaged_table.csv")
    
    return directly_imaged_table


#call functions

#directly_imaged_final_table_1=get_results("optimal",sigma_OPD_max,sigma_TipTilt_max,False,True,False)
'''
directly_imaged_final_table_2=get_results("optimal",10*sigma_OPD_max,sigma_TipTilt_max,False,True,False)
directly_imaged_final_table_3=get_results("optimal",sigma_OPD_max,10*sigma_TipTilt_max,False,True,False)
all_sigmas = pd.concat([directly_imaged_final_table_1,directly_imaged_final_table_2,directly_imaged_final_table_3], axis=1)
all_sigmas.to_csv("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_files\\all_sigmas2.csv")
'''
#a,b=main(data.TWOMASS_J01225093_2439505,data.TWOMASS_J01225093_2439505_b,"optimal",sigma_OPD_max,sigma_TipTilt_max,True,True,False)
#a,b=main(data.kap_And,data.kap_A_b,"optimal",sigma_OPD_max,sigma_TipTilt_max,True,True,False)
#a,b=main(data.bet_Pic,data.bet_Pic_b,"optimal",sigma_OPD_max,sigma_TipTilt_max,True,True,False)
#a,b=main(data.TYC_8998_760_1,data.TYC_8998_760_1_b,"optimal",sigma_OPD_max,sigma_TipTilt_max,True,True,False)
#a,b=main(data.HD_95086,data.HD_95086_b,"optimal",sigma_OPD_max,sigma_TipTilt_max,True,True,False)

#a,b=main(data.tau_Cet,data.tau_Cet_f,1.5*10**-6,sigma_OPD_zero,sigma_TipTilt_zero,True,False,False)

############################-----------------testing-------------------##############################

'''
#Reproduce figure 3
x=np.linspace(0.1*10**-6,100*10**-6,10000)

def output_planet_flux_per_wl_test(star,planet,wl):
    reflection=planet.albedo*np.pi*star.radius**2/star.dist**2*blackbody(star.temp,wl)*phase_fct_refl()*planet.radius**2/(4*planet.dist**2)
    thermal_emission=np.pi*planet.radius**2/star.dist**2*phase_fct_thermal()*blackbody(planet.temp,wl)
    total_output=reflection+thermal_emission
    
    return total_output

def output_stellar_flux_per_wl_test(star,wl):
    F=np.pi*star.radius**2/star.dist**2*blackbody(star.temp,wl)
        
    return F

exoplanet=np.empty_like(x)
for i in range (x.size):
    exoplanet[i]=output_planet_flux_per_wl_test(data.proxima_centauri,data.proxima_synthetic_planet,x[i])*x[i]**2/scicon.c*10**26

star_=np.empty_like(x)
for i in range(x.size):
    star_[i]=output_stellar_flux_per_wl_test(data.proxima_centauri,x[i])*x[i]**2/scicon.c*10**26

thermal_emission=np.empty_like(x)
Omega=np.empty_like(x)
for i in range(x.size):
    thermal_emission[i]=thermal_emission_instrument_per_wl(x[i],data.test_telescope)*x[i]**2/scicon.c*10**26*instrument_thermal_multiplicator
    
local_zodi=np.empty_like(x)
Omega=np.empty_like(x)
for i in range(x.size):
    local_zodi[i]=local_zodiacal_emission_per_wl(data.proxima_centauri,x[i])*x[i]**2/scicon.c*10**26
    Omega[i]=(0.512*x[i]/data.test_telescope.D)**2
    local_zodi[i]=Omega[i]*local_zodi[i]*np.pi
   
plt.loglog(10**6*x,exoplanet,label='exoplanet')
plt.loglog(10**6*x,star_,label='star')
plt.loglog(10**6*x,thermal_emission,label='Thermal emission')
plt.loglog(10**6*x,local_zodi,label='Local zodi')
plt.xlabel("Wavelength (\u03BCm)")
plt.ylabel("Input signals (Jy)")
plt.grid(True)
plt.xlim(10**-1,10**2)
plt.ylim(10**-10,10**2)
plt.legend(loc='lower right')
plt.show()


#Reproduce Figure 4
y=np.linspace(0.5*10**-6,30*10**-6,1000)


planet_signal=np.empty_like(y)
for i in range(y.size):
    planet_signal[i]=output_planet_flux(data.proxima_centauri,data.proxima_synthetic_planet,y[i],y[i]/R,data.PROBA)

  
shot_noise=np.empty_like(y)
for i in range(y.size):
    shot_noise[i]=np.sqrt(output_stellar_flux(data.proxima_centauri,y[i],y[i]/R,data.PROBA,sigma_OPD_max,sigma_TipTilt_max)+planet_signal[i]+output_background_flux(data.proxima_centauri,data.proxima_synthetic_planet,y[i],y[i]/R,data.PROBA))
    
    
instrumental_noise=np.empty_like(y)
for i in range(y.size):
    instrumental_noise[i]=output_stellar_flux(data.proxima_centauri,y[i],y[i]/R,data.PROBA,sigma_OPD_max,sigma_TipTilt_max)*np.sqrt((var_phase_error(y[i],sigma_OPD_max)**2+var_frac_intensity(y[i],sigma_TipTilt_max)**2)/8)

SNR_test=np.empty_like(y)
for i in range(y.size):
    SNR_test[i]=planet_signal[i]/np.sqrt(shot_noise[i]**2+instrumental_noise[i]**2)
    
optimal_wavelength=y[np.argmax(SNR_test)]
print(optimal_wavelength)


plt.plot(10**6*y,planet_signal,label="Planet")
plt.plot(10**6*y,shot_noise,label="Shot noise")
plt.plot(10**6*y,instrumental_noise,label="Instrumental noise")
plt.legend(loc='best')
plt.yscale("log")
plt.xlabel("Wavelength (\u03BCm)")
plt.ylabel("Detected signals (ph/s/bin)")
plt.grid(True)
plt.ylim(10**-6,10**5)
plt.xlim(0,30)
plt.show()

plt.plot(10**6*y,SNR_test,label="SNR")
plt.legend(loc='best')
plt.yscale("log")
plt.xlabel("Wavelength (\u03BCm)")
plt.ylabel("signal/noise")
plt.grid(True)
plt.ylim(10**-6,10**-1)
plt.xlim(0,30)
plt.show()



#heatmap of local zodiacal function
def local_zodiacal_emission_per_wl_test(long,lat,wl):
    I=4*10**(-8)*(blackbody(265,wl)+NIR_dust_albedo*blackbody(data.T_Sun,wl)*(data.R_Sun/(1.5*data.AU))**2)*np.sqrt(np.pi/(np.arccos(np.cos(long)*np.cos(lat)))/(np.sin(lat)**2+0.36*(wl/(11*10**(-6)))**(-0.8)*np.cos(lat)**2))

    return I

wavelength=5*10**-6
long=np.linspace(np.pi/2,3*np.pi/2,1000)
lat=np.linspace(-np.pi/2,np.pi/2,1000)
l,beta=np.meshgrid(long,lat)
z=local_zodiacal_emission_per_wl_test(l,beta,wavelength)
z*=wavelength**2/scicon.c*10**26*10**-6
df=pd.DataFrame(z,index=-lat*180/np.pi,columns=long*180/np.pi)
ax=sb.heatmap(data=df,cbar_kws={'label': '$MJy$ $sr^{-1}$'})
ax.set_title('\u03BB=5\u03BCm')
ax.set_xticks([0,166,333,499,666,832,999])
ax.set_xticklabels([90,120,150,180,210,240,270],rotation='horizontal')
ax.set_yticks([0,166,333,499,666,832,999])
ax.set_yticklabels([90,60,30,0,-30,-60,-90],rotation='horizontal')
ax.set_ylabel('\u03B2 [\u00B0]')
ax.set_xlabel('$\lambda_{rel}$ [\u00B0]')
plt.show()


#to check transmission maxima
baseline=np.linspace(0,3,10000)
T=np.sin(np.pi*np.tan(data.TYC_8998_760_1_b.dist/data.TYC_8998_760_1.dist)*baseline/(2.59*10**-6))**2
plt.plot(baseline,T)
plt.legend(loc='best')
plt.title("Transmission function for Cubesat")
plt.ylabel("Transmission function []")
plt.xlabel("Baselength [m]")
plt.vlines(np.array([0.5,1]),0,1)
plt.show()


#plot transmission function
u=np.linspace(-data.PROBA.OWA,data.PROBA.OWA,1000)
transmission=(2*sp.j1(np.pi*u*data.PROBA.D/data.PROBA.wl)/(np.pi*u*data.PROBA.D/data.PROBA.wl))**2*np.sin(np.pi*u*data.PROBA.b/data.PROBA.wl)**2
plt.plot(u,transmission)
plt.vlines(data.PROBA.IWA,0,1)
plt.show()


#reproduce colin's results for integration time
test_table=pd.DataFrame(columns=["CubeSat6U","CubeSat12U","PROBA"])
for i in range(len(data.close_systems_planets)):
    time_list=main(data.close_systems_stars[i],data.close_systems_planets[i],1.5*10**-6,sigma_OPD_zero,sigma_TipTilt_zero,False,False,False)[0]
    new_row=pd.DataFrame([time_list],index=[data.close_systems_planets[i].name],columns=["CubeSat6U","CubeSat12U","PROBA"])
    test_table=test_table.append(new_row)
print(test_table)
#test_table.to_csv("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_files\\integration_times_table2.csv")


#reproduce figure 7
OPD_array=np.linspace(0,75*10**(-9),50)
TipTilt_array=np.linspace(0,500*np.pi/(180*3600)*10**-3,50)
t_int_CubeSat12U=np.empty_like(OPD_array)
t_int_PROBA=np.empty_like(OPD_array)
t_int2_CubeSat12U=np.empty_like(TipTilt_array)
t_int2_PROBA=np.empty_like(TipTilt_array)
for i in range(OPD_array.size):
    t_int_CubeSat12U[i]=integration_time(data.proxima_centauri,data.proxima_synthetic_planet,"optimal",data.CubeSat12U,OPD_array[i],sigma_TipTilt_zero,False,False,False)[0]
    t_int_PROBA[i]=integration_time(data.proxima_centauri,data.proxima_synthetic_planet,"optimal",data.PROBA,OPD_array[i],sigma_TipTilt_zero,False,False,False)[0]
    t_int2_CubeSat12U[i]=integration_time(data.proxima_centauri,data.proxima_synthetic_planet,1.5*10**-6,data.CubeSat12U,sigma_OPD_zero,TipTilt_array[i],False,False,False)[0]
    t_int2_PROBA[i]=integration_time(data.proxima_centauri,data.proxima_synthetic_planet,2.5*10**-6,data.PROBA,sigma_OPD_zero,TipTilt_array[i],False,False,False)[0]


plt.plot(10**9*OPD_array,t_int_CubeSat12U,label="CubeSat12U")
plt.plot(10**9*OPD_array,t_int_PROBA,label="PROBA")
plt.legend(loc='best')
plt.yscale("log")
plt.xlim(0,80*10**(-9))
plt.ylim(0,10**6)
plt.hlines(np.array([24]),0,10**6)
plt.show()

plt.plot(TipTilt_array*np.pi/(180*3600)*10**-3,t_int2_CubeSat12U,label="CubeSat12U")
plt.plot(TipTilt_array*np.pi/(180*3600)*10**-3,t_int2_PROBA,label="PROBA")
plt.legend(loc='best')
plt.yscale("log")
plt.xlim(0,500*np.pi/(180*3600)*10**-3)
plt.ylim(0,10**6)
plt.hlines(np.array([24]),0,10**6)
plt.show()
'''