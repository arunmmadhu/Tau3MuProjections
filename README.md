# Projections

Run the following before using any code in this repository:

~~~
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
~~~

Run the following in the makeYield_ZTT directory in ROOT (for Limit/Significance scan):
~~~
.x makeYield_fromBDTFit_Combine_2018.C
~~~

Run the following to get projected signal/background yields at future luminosities:

~~~
makeYield_ZTT/makeYield_fromBDTFit_Combine.C
makeYield_W/makeYield_fromBDTFit_Combine_W.C
~~~

These give the limits at different luminosities:

~~~
makeYield_ZTT/LumiLimit.py
makeYield_HF/LumiLimit_HF.py
makeYield_W/LumiLimit_W.py
~~~

To make the final plot:

~~~
python3 MakePlot.py
~~~
