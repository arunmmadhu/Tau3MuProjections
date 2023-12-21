# Projections

Run the following to get projected signal/background yields at future luminosities:

~~~
makeYield_ZTT/makeYield_fromBDTFit_Combine.C
makeYield_HF/makeYield_fromBDTFit_Combine.C
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
