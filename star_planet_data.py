# -*- coding: utf-8 -*-
import numpy as np

#constants
parsec=3.0857*10**16
R_Sun=6.96*10**8
T_Sun=5778
AU=1.496*10**11
m_jup=1.898*10**27
r_jup=7.149*10**7


#define temp in K, dist in parsec, radius in R_Sun
class star:
    def __init__(self,temp,dist,radius,longitude,latitude,name):
        self.temp = temp #in K
        self.dist = dist*parsec # in m
        self.radius = radius*R_Sun #in m
        self.solid_angle=np.tan(self.radius/self.dist) # in radian measure, this is actually the radius not the diameter
        self.long=(longitude[0]+longitude[1]/60+longitude[2]/3600)*np.pi/12 # in radian measure
        self.lat=(latitude[0]+latitude[1]/60+latitude[2]/3600)*np.pi/180 # in radian measure
        self.name=name
        
     
#define temp in K, distance to star in AU, radius in r_jupiter, albedo in [0,1] (dimensionless)
class planet:
    def __init__(self,temp,dist,radius,albedo,name):
        self.temp=temp # in K
        self.dist=dist*AU # in m
        self.radius=radius*r_jup # in m
        self.albedo=albedo
        self.name=name


#define all in SI-Units
class telescope:
    def __init__(self,baseline,diameter,temp,optimal_wl,name):
        self.b=baseline
        self.D=diameter
        self.temp=temp
        self.wl=optimal_wl
        self.IWA=self.wl/(2*self.b)
        self.OWA=0.514*self.wl/self.D
        self.name=name
        
    def adjust_baseline(self,adjusted_bl):
        self.b=adjusted_bl


def estimate_planet_temp(T_star,Albedo,R_star,planet_dist):
    T_planet=T_star*(1-Albedo)**(1/4)*(R_star*R_Sun/(2*planet_dist*AU))**(1/2)
    
    return T_planet


#define telescopes
CubeSat6U=telescope(0.5,0.08,150,5*10**-6,"CubeSat6U")
CubeSat12U=telescope(1,0.08,150,5*10**-6,"CubeSat12U")
PROBA=telescope(5,0.25,150,5*10**-6,"PROBA")
FKSI=telescope(12.5,0.5,60,5*10**-6,"FKSI")
test_telescope=telescope(14.4,1,150,1.5*10**-6,"Test telescope")




#define stars and planets

######################################   close systems (from colin's paper)   ######################

proxima_centauri=star(3054,1.29,0.14,np.array([14,29,43]),-1*np.array([62,40,46]),"Proxima Centauri")
proxima_synthetic_planet=planet(137,0.077,2.36/11.21,0.5,"Proxima synthetic planet")
proxima_centauri_b=planet(216,0.05,0.098,0.3,"Proxima centauri b")
proxima_centauri_c=planet(39,1.5,0.30,0.3,"Proxima centauri c")

Barnard=star(3278,1.83,0.194,np.array([17,57,49]),np.array([4,41,36]),"Barnard's star")
Barnard_b=planet(105,0.404,1.5/11.21,0.3,"Barnard's star b")

Wolf_359=star(2800,2.39,0.16,np.array([10,56,29]),np.array([7,0,52]),"Wolf 359")
Wolf_359_b=planet(36,1.845,5.7/11.21,0.3,"Wolf 359 b")
Wolf_359_c=planet(368,0.018,1.6/11.21,0.3,"Wolf 359 c")

Lalande_21185=star(3828,2.547,0.4,np.array([11,3,20]),np.array([35,58,12]),"Lalande 21185")
Lalande_21185_b=planet(370,0.0789,1.4/11.21,0.3,"Lalande 21185 b")

eps_Eridani=star(5084,3.212,0.735,np.array([3,32,56]),-1*np.array([9,27,30]),"eps Eridani")
eps_Eridani_b=planet(116,3.39,10/11.21,0.3,"eps Eridani b")

GJ887=star(3688,3.3,0.5,np.array([23,5,52]),np.array([-35,51,11]),"GJ887")
GJ887_b=planet(468,0.07,0.14,0.3,"GJ887 b")
GJ887_c=planet(352,0.12,0.238,0.3,"GJ887 c")

Ross_128=star(3192,3.37,0.21,np.array([11,47,44]),np.array([0,48,16]),"Ross 128")
Ross_128_b=planet(256,0.05,1.11/11.21,0.3,"Ross 128")

GJ_15A=star(3567,3.587,0.386,np.array([0,18,23]),np.array([44,1,23]),"GJ 15A")
GJ_15A_b=planet(364,0.072,1.4/11.21,0.3,"GJ 15A b")

YZ_Cet=star(3056,3.71,0.169,np.array([1,12,31]),-1*np.array([16,59,56]),"YZ Cet")
YZ_Cet_b=planet(432,0.02,0.9/11.21,0.3,"YZ Cet b")
YZ_Cet_c=planet(376,0.02,1/11.21,0.3,"YZ Cet c")
YZ_Cet_d=planet(327,0.03,1/11.21,0.3,"YZ Cet d")

eps_Ind_A=star(4630,3.6,0.7,np.array([22,3,22]),-1*np.array([56,47,10]),"eps Ind A")
eps_Ind_A_b=planet(51,11.6,16.2/11.21,0.3,"eps Ind A b")

tau_Cet=star(5344,3.7,0.8,np.array([1,44,4]),-1*np.array([15,56,15]),"tau Cet")
tau_Cet_e=planet(286,0.5,1.6/11.21,0.3,"tau Cet e")
tau_Cet_f=planet(182,1.3,1.6/11.21,0.3,"tau Cet f")
tau_Cet_g=planet(576,0.1,1.2/11.21,0.3,"tau Cet g")
tau_Cet_h=planet(426,0.2,1.2/11.21,0.3,"tau Cet h")

GJ1061=star(2953,3.7,0.2,np.array([3,36,0]),-1*np.array([44,30,46]),"GJ1061")
GJ1061_b=planet(355,0.02,1.1/11.21,0.3,"GJ1061 b")
GJ1061_c=planet(275,0.04,1.2/11.21,0.3,"GJ1061 c")
GJ1061_d=planet(226,0.05,1.2/11.21,0.3,"GJ1061 d")

GJ273=star(3382,3.8,0.3,np.array([7,27,25]),np.array([5,13,33]),"GJ273")
GJ273_b=planet(268,0.1,1.4/11.21,0.3,"GJ273 b")
GJ273_c=planet(426,0.04,1.1/11.21,0.3,"GJ273 c")
GJ273_d=planet(96,0.7,3.6/11.21,0.3,"GJ273 d")
GJ273_e=planet(88,0.8,3.4/11.21,0.3,"GJ273 e")

Teegarden=star(2637,3.8,0.107,np.array([2,53,1]),np.array([16,52,53]),"Teegarden's star")
Teegarden_b=planet(240,0.03,1/11.21,0.3,"Teegarden's star b")
Teegarden_c=planet(181,0.04,1/11.21,0.3,"Teegarden's star c")



#lists of all stars/planets
close_systems_stars=[proxima_centauri,proxima_centauri,Barnard,Wolf_359,Wolf_359,Lalande_21185,eps_Eridani,GJ887,GJ887,Ross_128,GJ_15A,YZ_Cet,YZ_Cet,YZ_Cet,eps_Ind_A,tau_Cet,tau_Cet,tau_Cet,tau_Cet,GJ1061,GJ1061,GJ1061,GJ273,GJ273,GJ273,GJ273,Teegarden,Teegarden]
close_systems_planets=[proxima_centauri_b,proxima_centauri_c,Barnard_b,Wolf_359_b,Wolf_359_c,Lalande_21185_b,eps_Eridani_b,GJ887_b,GJ887_c,Ross_128_b,GJ_15A_b,YZ_Cet_b,YZ_Cet_c,YZ_Cet_d,eps_Ind_A_b,tau_Cet_e,tau_Cet_f,tau_Cet_g,tau_Cet_h,GJ1061_b,GJ1061_c,GJ1061_d,GJ273_b,GJ273_c,GJ273_d,GJ273_e,Teegarden_b,Teegarden_c]



######################################## directly imaged planets ###################################

HIP_65426=star(8840,108.875,1.77,np.array([13,24,36]),np.array([51,30,34]),"HIP 65426")
HIP_65426_b=planet(1500,92,1.5,0.3,"HIP 65426 b")

#Radius unknown; assumption density like jupiter (also often done later)
USco_CTIO_108=star(2846,144,0.345,np.array([16,5,54]),-1*np.array([18,18,45]),"USco CTIO 108")
USco_CTIO_108_b=planet(2350,670,14**(1/3),0.3,"USco CTIO 108 b")

TWOMASS_J01225093_2439505=star(3309,33.83,0.366,np.array([1,22,51]),-1*np.array([24,39,53]),"2MASS J01225093 2439505")
TWOMASS_J01225093_2439505_b=planet(1600,52,1,0.3,"2MASS J01225093 2439505 b")

#temp from https://arxiv.org/pdf/1611.00364.pdf
TWOMASSJ22362452_4751425=star(4033,69.57,0.64,np.array([22,36,25]),np.array([47,51,42]),"2MASS J22362452+4751425")
TWOMASSJ22362452_4751425_b=planet(1100,230,12.5**(1/3),0.3,"2MASS J22362452+4751425 b")

kap_And=star(10839,50,2.31,np.array([23,40,25]),np.array([44,20,2]),"kap A")
kap_A_b=planet(1900,55,13.62**(1/3),0.3,"kap A b")

GSC_06214_00210=star(4432,108.5,0.905,np.array([16,12,55]),-1*np.array([20,43,10]),"GSC 06214 00210")
GSC_06214_00210_b=planet(2200,320,1.8,0.3,"GSC 06214 00210 b")

ONERXS_J1609291_210524=star(3993,139.135,1.313,np.array([16,9,30]),-1*np.array([21,4,59]),"1RXS J160929.1-210524")
ONERXS_J1609291_210524_b=planet(1800,330,1.664,0.3,"1RXS J160929.1-210524 b")

#temp from https://arxiv.org/pdf/astro-ph/0609464.pdf
HN_Peg=star(6186,18.119,1.011,np.array([21,44,32]),np.array([14,46,17]),"HN Peg")
HN_Peg_b=planet(1130,773,1.051,0.3,"HN Peg b")

HD_203030=star(5472,39.25,0.856,np.array([21,18,58]),np.array([26,13,50]),"HD 203030")
HD_203030_b=planet(1206,487,24.09**(1/3),0.3,"HD 203030 b")

GQ_Lup=star(4360,151.2,1.94,np.array([15,49,12]),-1*np.array([35,39,5]),"GQ Lup")
GQ_Lup_b=planet(2650,100,3,0.3,"GQ Lup b")

#temp from https://arxiv.org/pdf/1412.5173.pdf
HD_100546=star(10331,109.7,1.853,np.array([11,33,25]),-1*np.array([70,11,41]),"HD 100546")
HD_100546_b=planet(932,53,6.9,0.3,"HD 100546 b")

#temp using stefan boltzmann
TWOMASS_J04414489_2301513=star(2936,120.4,0.232,np.array([4,41,45]),np.array([23,1,51]),"2MASS J04414489+2301513")
TWOMASS_J04414489_2301513_b=planet(estimate_planet_temp(2936,0.3,0.232,15),15,7.5**(1/3),0.3,"2MASS J04414489+2301513 b")

Fomalhaut=star(8590,7.704,1.842,np.array([22,57,39]),-1*np.array([29,37,20]),"Fomalhaut")
Fomalhaut_b=planet(estimate_planet_temp(8590,0.3,1.842,117),117,0.208**(1/3),0.3,"Fomalhaut b")

ROXs_12=star(4059,136.65,1.066,np.array([16,26,28]),-1*np.array([25,26,48]),"ROXs 12")
ROXs_12_b=planet(estimate_planet_temp(4059,0.3,1.066,210),210,16**(1/3),0.3,"ROXs 12 b")

PDS_70=star(4152,113.1,1.186,np.array([14,8,10]),-1*np.array([41,23,53]),"PDS 70")
PDS_70_b=planet(1204,20,2.72,0.3,"PDS 70 b")
PDS_70_c=planet(995,34,2.04,0.3,"PDS 70 c")

Oph_11=star(2878,136.2,0.233,np.array([16,22,25]),-1*np.array([24,5,14]),"Oph 11")
Oph_11_b=planet(2175,243,14**(1/3),0.3,"Oph 11 b")

HR_8799=star(7339,41.24,1.493,np.array([23,7,29]),np.array([21,8,3]),"HR 8799")
HR_8799_b=planet(1200,68,1.2,0.3,"HR 8799 b")
HR_8799_c=planet(1200,38,1.2,0.3,"HR 8799 c")
HR_8799_d=planet(1300,24,1.2,0.3,"HR 8799 d")
HR_8799_e=planet(1150,16.4,1.17,0.3,"HR 8799 e")

HIP_78530=star(10690,136.76,2.109,np.array([16,1,55]),-1*np.array([21,58,50]),"HIP 78530")
HIP_78530_b=planet(2700,740,23**(1/3),0.3,"HIP 78530 b")

HD_95086=star(7883,86.23,1.487,np.array([10,57,3]),-1*np.array([68,40,2]),"HD 95086")
HD_95086_b=planet(916,55.7,4.5**(1/3),0.3,"HD 95086 b")

HD_106906=star(6798,103.03,2.124,np.array([12,17,53]),-1*np.array([55,58,32]),"HD 106906")
HD_106906_b=planet(1800,650,11**(1/3),0.3,"HD 106906 b")

GJ_504=star(6291,17.53,1.352,np.array([13,16,46]),np.array([9,25,30]),"GJ 504")
GJ_504_b=planet(510,43.5,4**(1/3),0.3,"GJ 504 b")

TWOMASS_J12073346_3932539=star(2825,64.31,0.2189,np.array([12,7,33]),-1*np.array([39,32,54]),"2MASS J12073346-3932539")
TWOMASS_J12073346_3932539_b=planet(1150,46,4**(1/3),0.3,"2MASS J12073346-3932539 b")

#no known planet distance from sun --> do not include in analysis (not in list)
#TWOMASS_J21402931_1625183_A=star(2300,33.197,0.1339,np.array([21,40,29]),np.array([16,25,17]),"2MASS J21402931+1625183 A")
#TWOMASS_J21402931_1625183_A_b=planet(2075,"??",0.92,0.3,"2MASS J21402931+1625183 A b")

#planet distance calculated from star distance and anglular separation (https://arxiv.org/pdf/astro-ph/0504658.pdf)
AB_Pic=star(5285,50.0475,0.999,np.array([6,19,13]),-1*np.array([58,3,15]),"AB Pic")
AB_Pic_b=planet(275,260,13.5**(1/3),0.3,"AB Pic b")

bet_Pic=star(8039,19.74,1.54,np.array([5,47,17]),-1*np.array([51,3,58]),"bet Pic")
bet_Pic_b=planet(1650,8.9,1.5,0.3,"bet Pic b")

#star radius and distance from https://en.wikipedia.org/wiki/CFBDSIR_J145829%2B101343
CFBDSIR_J145829_101343=star(581,32,0.15,np.array([164,58,29]),np.array([10,13,43]),"CFBDSIR J145829+101343")
CFBDSIR_J145829_101343_b=planet(370,2.6,10.5**(1/3),0.3,"CFBDSIR J145829+101343 b")

#star radius from https://en.wikipedia.org/wiki/CHXR_73; planet temp from stefan boltzmann
CHXR_73=star(3099,190,0.83,np.array([11,6,29]),-1*np.array([77,37,33]),"CHXR 73")
CHXR_73_b=planet(estimate_planet_temp(3099,0.3,0.83,210),210,12.569**(1/3),0.3,"CHXR 73 b")

CT_Cha=star(4403,190.72,1.68,np.array([11,4,9]),-1*np.array([76,27,19]),"CT Cha")
CT_Cha_b=planet(2600,440,2.2,0.3,"CT Cha b")

DH_Tau=star(4371,134.846,1.367,np.array([4,29,42]),np.array([26,32,58]),"DH Tau")
DH_Tau_b=planet(2200,330,2.7,0.3,"DH Tau b")

GU_Psc=star(3420,47.55,0.45,np.array([1,12,35]),np.array([17,3,54]),"GU Psc")
GU_Psc_b=planet(1050,2000,1.265,0.3,"GU Psc b")

#VHS J125601.92-125723.9 b appears to be in a binary system

#star radius not found
#WD_0806_661=star(9552,19.2444,"??",np.array([8,6,55]),-1*np.array([66,18,21]),"WD 0806-661")
#WD_0806_661_b=planet(322,2500,7.5**(1/3),0.3,"WD 0806-661 b")

#star distance from https://en.wikipedia.org/wiki/WISE_1217%2B1626
WISEP_J121756_91_162640_2_A=star(575,8.8,0.091,np.array([12,17,57]),np.array([16,26,40]),"WISEP J121756.91+162640.2 A")
WISEP_J121756_91_162640_2_A_b=planet(450,8,0.934,0.3,"WISEP J121756.91+162640.2 A b")

#No information on planet temperature or size; system is still in an accreting stage (see https://arxiv.org/pdf/1511.07456.pdf)
#LkCa_15=star(4589,158.152,1.61,np.array([4,39,18]),np.array([22,21,3]),"LkCa 15")
#LkCa_15_b=planet()
#LkCa_15_c=planet()

#planet temp from https://arxiv.org/pdf/1505.01747.pdf
TWOMASS_J02192210_3925225=star(3064,40.09,0.281,np.array([2,19,22]),-1*np.array([39,25,23]),"2MASS J02192210-3925225")
TWOMASS_J02192210_3925225_b=planet(1683,156,1.44,0.3,"2MASS J02192210-3925225 b")

#planet temp from https://arxiv.org/pdf/1912.04284.pdf
TYC_8998_760_1=star(4573,94.6177,0.954,np.array([13,25,12]),-1*np.array([64,56,21]),"TYC 8998-760-1")
TYC_8998_760_1_b=planet(1727,162,14**(1/3),0.3,"TYC 8998-760-1 b")


#lists of all stars/planets
directly_imaged_systems_stars=[HIP_65426,USco_CTIO_108,TWOMASS_J01225093_2439505,TWOMASSJ22362452_4751425,kap_And,GSC_06214_00210,ONERXS_J1609291_210524,HN_Peg,HD_203030,GQ_Lup,HD_100546,TWOMASS_J04414489_2301513,Fomalhaut,ROXs_12,PDS_70,PDS_70,Oph_11,HR_8799,HR_8799,HR_8799,HR_8799,HIP_78530,HD_95086,HD_106906,GJ_504,TWOMASS_J12073346_3932539,AB_Pic,bet_Pic,CFBDSIR_J145829_101343,CHXR_73,CT_Cha,DH_Tau,GU_Psc,WISEP_J121756_91_162640_2_A,TWOMASS_J02192210_3925225,TYC_8998_760_1]
directly_imaged_systems_planets=[HIP_65426_b,USco_CTIO_108_b,TWOMASS_J01225093_2439505_b,TWOMASSJ22362452_4751425_b,kap_A_b,GSC_06214_00210_b,ONERXS_J1609291_210524_b,HN_Peg_b,HD_203030_b,GQ_Lup_b,HD_100546_b,TWOMASS_J04414489_2301513_b,Fomalhaut_b,ROXs_12_b,PDS_70_b,PDS_70_c,Oph_11_b,HR_8799_b,HR_8799_c,HR_8799_d,HR_8799_e,HIP_78530_b,HD_95086_b,HD_106906_b,GJ_504_b,TWOMASS_J12073346_3932539_b,AB_Pic_b,bet_Pic_b,CFBDSIR_J145829_101343_b,CHXR_73_b,CT_Cha_b,DH_Tau_b,GU_Psc_b,WISEP_J121756_91_162640_2_A_b,TWOMASS_J02192210_3925225_b,TYC_8998_760_1_b]
