# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import exoplanets as exo
import star_planet_data as data

path="C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\Colin_Model_directly_imaged\\"
wl_values=np.linspace(10**-6,6*10**-6,21)
row_names=["2MASS J01225093 2439505 b","kap A b","bet Pic b","TYC 8998-760-1 b","HD 95086 b"]
save_files=False


#####################################  CubeSat6U plots  #########################################

CubeSat_6U_t_file=pd.read_csv(path+"CubeSat 6U\\tipopd\\T_iCubeSat6Uopdtilt.csv",sep=",",names=wl_values)
CubeSat_6U_t_file.index=row_names

TWOMASS_t_6U_colin=CubeSat_6U_t_file.iloc[0].values
kap_t_6U_colin=CubeSat_6U_t_file.iloc[1].values
bet_Pic_t_6U_colin=CubeSat_6U_t_file.iloc[2].values
TYC_t_6U_colin=CubeSat_6U_t_file.iloc[3].values
HD_t_6U_colin=CubeSat_6U_t_file.iloc[4].values


plt.plot(wl_values*10**6,TWOMASS_t_6U_colin*3600)
plt.ylim(0,100)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat6U 2Mass (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_6U\\2MASS_J01225093_2439505_b_Colin_6U.png")
plt.show()

TWOMASS_t_6U=np.empty_like(wl_values)
for i in range(wl_values.size):
    TWOMASS_t_6U[i]=exo.integration_time(data.TWOMASS_J01225093_2439505,data.TWOMASS_J01225093_2439505_b,wl_values[i],data.CubeSat6U,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,TWOMASS_t_6U)
plt.ylim(0,100)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat6U 2Mass (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_6U\\2MASS_J01225093_2439505_b_Philipp_6U.png")
plt.show()
print("-----------------------------------------------")


plt.plot(wl_values*10**6,kap_t_6U_colin*3600)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat6U kap A b (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_6U\\kap_And_b_Colin_6U.png")
plt.show()

kap_t_6U=np.empty_like(wl_values)
for i in range(wl_values.size):
    kap_t_6U[i]=exo.integration_time(data.kap_And,data.kap_A_b,wl_values[i],data.CubeSat6U,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,kap_t_6U)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat6U Kap A b (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_6U\\kap_And_b_Philipp_6U.png")
plt.show()
print("-----------------------------------------------")


plt.plot(wl_values*10**6,bet_Pic_t_6U_colin*3600)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat6U bet Pic b (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_6U\\bet_Pic_b_Colin_6U.png")
plt.show()

bet_Pic_t_6U=np.empty_like(wl_values)
for i in range(wl_values.size):
    bet_Pic_t_6U[i]=exo.integration_time(data.bet_Pic,data.bet_Pic_b,wl_values[i],data.CubeSat6U,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,bet_Pic_t_6U)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat6U bet Pic (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_6U\\bet_Pic_b_Philipp_6U.png")
plt.show()
print("-----------------------------------------------")


plt.plot(wl_values*10**6,TYC_t_6U_colin*3600)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat6U TYC (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_6U\\TYC_8998_760_1_b_Colin_6U.png")
plt.show()

TYC_t_6U=np.empty_like(wl_values)
for i in range(wl_values.size):
    TYC_t_6U[i]=exo.integration_time(data.TYC_8998_760_1,data.TYC_8998_760_1_b,wl_values[i],data.CubeSat6U,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,TYC_t_6U)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat6U TYC (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_6U\\TYC_8998_760_1_b_Philipp_6U.png")
plt.show()
print("-----------------------------------------------")


plt.plot(wl_values*10**6,HD_t_6U_colin*3600)
plt.ylim(0,2000)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat6U HD (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_6U\\HD_95086_b_Colin_6U.png")
plt.show()

HD_t_6U=np.empty_like(wl_values)
for i in range(wl_values.size):
    HD_t_6U[i]=exo.integration_time(data.HD_95086,data.HD_95086_b,wl_values[i],data.CubeSat6U,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,HD_t_6U)
plt.ylim(0,2000)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat6U HD (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_6U\\HD_95086_b_Philipp_6U.png")
plt.show()
print("-----------------------------------------------")


########################################## CubeSat12U plots  #################################

CubeSat_12U_t_file=pd.read_csv(path+"CubeSat 12U\\opdtip\\T_iCubeSat12Uopdtilt.csv",sep=",",names=wl_values)
CubeSat_12U_t_file.index=row_names

TWOMASS_t_12U_colin=CubeSat_12U_t_file.iloc[0].values
kap_t_12U_colin=CubeSat_12U_t_file.iloc[1].values
bet_Pic_t_12U_colin=CubeSat_12U_t_file.iloc[2].values
TYC_t_12U_colin=CubeSat_12U_t_file.iloc[3].values
HD_t_12U_colin=CubeSat_12U_t_file.iloc[4].values


plt.plot(wl_values*10**6,TWOMASS_t_12U_colin*3600)
plt.ylim(0,100)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat12U 2Mass (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_12U\\2MASS_J01225093_2439505_b_Colin_12U.png")
plt.show()

TWOMASS_t_12U=np.empty_like(wl_values)
for i in range(wl_values.size):
    TWOMASS_t_12U[i]=exo.integration_time(data.TWOMASS_J01225093_2439505,data.TWOMASS_J01225093_2439505_b,wl_values[i],data.CubeSat12U,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,TWOMASS_t_12U)
plt.ylim(0,100)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat12U 2Mass (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_12U\\2MASS_J01225093_2439505_b_Philipp_12U.png")
plt.show()
print("-----------------------------------------------")


plt.plot(wl_values*10**6,kap_t_12U_colin*3600)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat12U kap A b (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_12U\\kap_A_b_Colin_12U.png")
plt.show()

kap_t_12U=np.empty_like(wl_values)
for i in range(wl_values.size):
    kap_t_12U[i]=exo.integration_time(data.kap_And,data.kap_A_b,wl_values[i],data.CubeSat12U,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,TWOMASS_t_12U)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat12U Kap A b (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_12U\\kap_A_b_Philipp_12U.png")
plt.show()
print("-----------------------------------------------")


plt.plot(wl_values*10**6,bet_Pic_t_12U_colin*3600)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat12U bet Pic b (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_12U\\bet_Pic_b_Colin_12U.png")
plt.show()

bet_Pic_t_12U=np.empty_like(wl_values)
for i in range(wl_values.size):
    bet_Pic_t_12U[i]=exo.integration_time(data.bet_Pic,data.bet_Pic_b,wl_values[i],data.CubeSat12U,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,bet_Pic_t_12U)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat12U bet Pic (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_12U\\bet_Pic_b_Philipp_12U.png")
plt.show()
print("-----------------------------------------------")


plt.plot(wl_values*10**6,TYC_t_12U_colin*3600)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat12U TYC (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_12U\\TYC_8998_760_1_b_Colin_12U.png")
plt.show()

TYC_t_12U=np.empty_like(wl_values)
for i in range(wl_values.size):
    TYC_t_12U[i]=exo.integration_time(data.TYC_8998_760_1,data.TYC_8998_760_1_b,wl_values[i],data.CubeSat12U,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,TYC_t_12U)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat12U TYC (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_12U\\TYC_8998_760_1_b_Philipp_12U.png")
plt.show()
print("-----------------------------------------------")


plt.plot(wl_values*10**6,HD_t_12U_colin*3600)
plt.ylim(0,20000)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat12U HD (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_12U\\HD_95086_b_Colin_12U.png")
plt.show()

HD_t_12U=np.empty_like(wl_values)
for i in range(wl_values.size):
    HD_t_12U[i]=exo.integration_time(data.HD_95086,data.HD_95086_b,wl_values[i],data.CubeSat12U,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,HD_t_12U)
plt.ylim(0,20000)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("CubeSat12U HD (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\CubeSat_12U\\HD_95086_b_Philipp_12U.png")
plt.show()
print("-----------------------------------------------")


########################################## PROBA plots  #################################

PROBA_t_file=pd.read_csv(path+"PROBA\\opdtip\\T_iPROBAopdtilt.csv",sep=",",names=wl_values)
PROBA_t_file.index=row_names

TWOMASS_t_PROBA_colin=PROBA_t_file.iloc[0].values
kap_t_PROBA_colin=PROBA_t_file.iloc[1].values
bet_Pic_t_PROBA_colin=PROBA_t_file.iloc[2].values
TYC_t_PROBA_colin=PROBA_t_file.iloc[3].values
HD_t_PROBA_colin=PROBA_t_file.iloc[4].values


plt.plot(wl_values*10**6,TWOMASS_t_PROBA_colin*3600)
plt.ylim(0,100)
plt.xlim(1,6)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("PROBA 2Mass (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\PROBA\\TWOMASS_J01225093_2439505_b_Colin_PROBA.png")
plt.show()

TWOMASS_t_PROBA=np.empty_like(wl_values)
for i in range(wl_values.size):
    TWOMASS_t_PROBA[i]=exo.integration_time(data.TWOMASS_J01225093_2439505,data.TWOMASS_J01225093_2439505_b,wl_values[i],data.PROBA,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,TWOMASS_t_PROBA)
plt.ylim(0,100)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("PROBA 2Mass (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\PROBA\\TWOMASS_J01225093_2439505_b_Philipp_PROBA.png")
plt.show()
print("-----------------------------------------------")


plt.plot(wl_values*10**6,kap_t_PROBA_colin*3600)
plt.ylim(0,2)
plt.xlim(1,6)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("PROBA kap A b (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\PROBA\\kap_A_b_Colin_PROBA.png")
plt.show()

kap_t_PROBA=np.empty_like(wl_values)
for i in range(wl_values.size):
    kap_t_PROBA[i]=exo.integration_time(data.kap_And,data.kap_A_b,wl_values[i],data.PROBA,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,kap_t_PROBA)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("PROBA kap A b (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\PROBA\\kap_A_b_Philipp_PROBA.png")
plt.show()
print("-----------------------------------------------")


plt.plot(wl_values*10**6,bet_Pic_t_PROBA_colin*3600)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("PROBA bet Pic b (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\PROBA\\bet_Pic_b_Colin_PROBA.png")
plt.show()

bet_Pic_t_PROBA=np.empty_like(wl_values)
for i in range(wl_values.size):
    bet_Pic_t_PROBA[i]=exo.integration_time(data.bet_Pic,data.bet_Pic_b,wl_values[i],data.PROBA,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,bet_Pic_t_PROBA)
plt.ylim(0,2)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("PROBA bet Pic (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\PROBA\\bet_Pic_b_Philipp_PROBA.png")
plt.show()
print("-----------------------------------------------")


plt.plot(wl_values*10**6,TYC_t_PROBA_colin*3600)
plt.ylim(0,100)
plt.xlim(1,6)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("PROBA TYC (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\PROBA\\TYC_8998_760_1_b_Colin_PROBA.png")
plt.show()

TYC_t_PROBA=np.empty_like(wl_values)
for i in range(wl_values.size):
    TYC_t_PROBA[i]=exo.integration_time(data.TYC_8998_760_1,data.TYC_8998_760_1_b,wl_values[i],data.PROBA,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,TYC_t_PROBA)
plt.ylim(0,100)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("PROBA TYC (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\PROBA\\TYC_8998_760_1_b_Philipp_PROBA.png")
plt.show()
print("-----------------------------------------------")


plt.plot(wl_values*10**6,HD_t_PROBA_colin*3600)
plt.ylim(0,2000)
plt.xlim(1,6)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("PROBA HD (Colin)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\PROBA\\HD_95086_b_Colin_PROBA.png")
plt.show()

HD_t_PROBA=np.empty_like(wl_values)
for i in range(wl_values.size):
    HD_t_PROBA[i]=exo.integration_time(data.HD_95086,data.HD_95086_b,wl_values[i],data.PROBA,exo.sigma_OPD_max,exo.sigma_TipTilt_max,False,True,False)[0]
plt.plot(wl_values*10**6,HD_t_PROBA)
plt.ylim(0,2000)
plt.ylabel("Integration time (s)")
plt.xlabel("Wavelength (\u03BCm)")
plt.title("PROBA HD (Philipp)")
if(save_files==True):
    plt.savefig("C:\\Users\\Phili\\OneDrive\\Desktop\\ETH\\Semesterarbeit_Exoplaneten\\output_Colin_vs_Philipp\\PROBA\\HD_95086_b_Philipp_PROBA.png")
plt.show()
print("-----------------------------------------------")
