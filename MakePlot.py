import os
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

import math

HF_lim = []
W_lim = []
HFW_lim = []

lumi = [97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0]

ZTT_lim=[[2.1719, 1.75, 0.7891, 0.5098, 0.3008, 0.2295, 0.1895, 0.1641, 0.1445], [0.8823, 0.7109, 0.3298, 0.219, 0.1316, 0.1022, 0.0844, 0.0731, 0.0649], [1.3155, 1.06, 0.4877, 0.319, 0.1898, 0.1469, 0.1213, 0.1058, 0.0929], [3.7605, 3.0467, 1.3437, 0.856, 0.4979, 0.3756, 0.3091, 0.2654, 0.2345], [6.2706, 5.0579, 2.207, 1.3795, 0.7964, 0.5928, 0.4863, 0.4154, 0.3673]]
ZTT_lim_const=[[2.0156, 1.6406, 0.7773, 0.5156, 0.3193, 0.25, 0.2109, 0.1855, 0.167], [0.7874, 0.6537, 0.3279, 0.2256, 0.1447, 0.1152, 0.0989, 0.087, 0.0796], [1.1904, 0.9852, 0.4824, 0.3253, 0.2075, 0.1637, 0.1391, 0.1224, 0.111], [3.5282, 2.8562, 1.3016, 0.8487, 0.512, 0.3961, 0.3322, 0.2895, 0.2598], [5.8847, 4.7418, 2.0967, 1.3415, 0.7925, 0.6054, 0.502, 0.4374, 0.3895]]

ZTT_lim_1_sqrtL=[]
ZTT_lim_1_L=[]
for lumino,lumival in enumerate(lumi):
    ZTT_lim_1_L.append(ZTT_lim[0][0]*(lumi[0]/lumi[lumino]))
    ZTT_lim_1_sqrtL.append(ZTT_lim[0][0]*math.sqrt(lumi[0]/lumi[lumino]))

HF_lim=[[0.3398, 0.292, 0.1738, 0.125, 0.0864, 0.0708, 0.061, 0.0542, 0.0493], [0.1646, 0.1414, 0.0856, 0.0615, 0.0425, 0.0348, 0.03, 0.0267, 0.0243], [0.229, 0.1967, 0.1173, 0.0848, 0.0587, 0.048, 0.0412, 0.037, 0.0333], [0.5239, 0.4487, 0.2655, 0.1909, 0.1316, 0.1075, 0.0926, 0.0825, 0.0751], [0.7899, 0.6757, 0.3974, 0.2848, 0.1954, 0.1599, 0.1378, 0.1225, 0.1115]]

HF_lim_1_sqrtL=[]
HF_lim_1_L=[]
for lumino,lumival in enumerate(lumi):
    HF_lim_1_L.append(HF_lim[0][0]*(lumi[0]/lumi[lumino]))
    HF_lim_1_sqrtL.append(HF_lim[0][0]*math.sqrt(lumi[0]/lumi[lumino]))

W_lim=[[0.5918, 0.5039, 0.2764, 0.1982, 0.1323, 0.1069, 0.0923, 0.0825, 0.0747], [0.282, 0.2401, 0.1339, 0.0976, 0.0662, 0.0535, 0.0461, 0.0413, 0.0379], [0.3958, 0.3349, 0.1873, 0.1345, 0.0905, 0.0735, 0.0634, 0.0567, 0.0517], [0.9122, 0.7744, 0.4194, 0.2981, 0.1977, 0.1592, 0.1374, 0.1221, 0.1113], [1.3623, 1.1511, 0.6158, 0.4339, 0.2848, 0.2282, 0.1956, 0.1743, 0.1583]]
W_lim_const=[[0.7754, 0.6699, 0.3848, 0.2803, 0.1899, 0.1548, 0.1338, 0.1196, 0.1089], [0.3877, 0.335, 0.1924, 0.1412, 0.0965, 0.0786, 0.0679, 0.0607, 0.0553], [0.5301, 0.458, 0.2645, 0.1934, 0.1315, 0.1072, 0.0926, 0.0828, 0.0754], [1.1621, 1.0008, 0.5693, 0.4147, 0.2802, 0.2276, 0.1967, 0.1759, 0.1601], [1.6715, 1.442, 0.8129, 0.5878, 0.3978, 0.3225, 0.2777, 0.2483, 0.226]]


W_lim_1_sqrtL=[]
W_lim_1_L=[]
for lumino,lumival in enumerate(lumi):
    W_lim_1_L.append(W_lim[0][0]*(lumi[0]/lumi[lumino]))
    W_lim_1_sqrtL.append(W_lim[0][0]*math.sqrt(lumi[0]/lumi[lumino]))



plt.plot(lumi, HF_lim[0], 'b--')
plt.plot(lumi, HF_lim_1_sqrtL, 'b*')
plt.plot(lumi, HF_lim_1_L, 'b*')
plt.fill_between(lumi, HF_lim[1], HF_lim[4], alpha=1.0, facecolor='yellow',label='_nolegend_')
plt.fill_between(lumi, HF_lim[2], HF_lim[3], alpha=1.0, facecolor='lime',label='_nolegend_')


plt.plot(lumi, W_lim[0], 'c--')
plt.plot(lumi, W_lim_1_sqrtL, 'c*')
plt.plot(lumi, W_lim_1_L, 'c*')
plt.fill_between(lumi, W_lim[1], W_lim[4], alpha=1.0, facecolor='yellow',label='_nolegend_')
plt.fill_between(lumi, W_lim[2], W_lim[3], alpha=1.0, facecolor='lime',label='_nolegend_')

plt.plot(lumi, ZTT_lim[0], 'r--')
plt.plot(lumi, ZTT_lim_1_sqrtL, 'r*')
plt.plot(lumi, ZTT_lim_1_L, 'r*')
plt.fill_between(lumi, ZTT_lim[1], ZTT_lim[4], alpha=1.0, facecolor='yellow',label='_nolegend_')
plt.fill_between(lumi, ZTT_lim[2], ZTT_lim[3], alpha=1.0, facecolor='lime',label='_nolegend_')


plt.xlabel('Integrated Luminosity $fb^{-1}$')
plt.ylabel('Expected Limit ($10^{-7}$)')
plt.legend(['HF',r'$HF\;\frac{1}{\sqrt{L}}$',r'$HF\;\frac{1}{L}$', 'W',r'$W\;\frac{1}{\sqrt{L}}$',r'$W\;\frac{1}{L}$','ZTT',r'$ZTT\;\frac{1}{\sqrt{L}}$',r'$ZTT\;\frac{1}{L}$'],loc="lower left", ncol = 3)
plt.xscale("log")

plt.yscale("log")
plt.savefig("test.png", format="png")
plt.show()

