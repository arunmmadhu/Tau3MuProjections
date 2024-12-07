# Projections

To create datacards in the `makeYield/ZTT`, `makeYield/ZTT` and `makeYield/W` directories, run the `./makeCards.py` file in these directories to create the datacards in the `makeYield/*/datacards_modified` directory.

For the ZTT limit projections, a limit scan needs to be run at each luminosity wrt BDT cuts. The steps are:



0) cmsenv is done with CMSSW_14_1_0_pre4:
```sh


cmsrel CMSSW_14_1_0_pre4;
cd CMSSW_14_1_0_pre4/src;
cmsenv;
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit;
scramv1 b clean; scramv1 b;

```


1) Define the categories to be run here (do not run the 'combined' category at this stage, but all the other categories can be run here):

[categories = ['taue','taumu','tauhA','tauhB','all']](https://github.com/T3MuAnalysisTools/Projections/blob/cb5efdac12f291e91bc00f29275c390868a90af3/makeYield/ZTT/makeCards.py#L486)

2) Set WhetherFitBDTandMakeCards to True
3) Run makeCards.py (`./makeCards.py`) in the `makeYield/ZTT` directory to create the datacards for the limit scan wrt bdt cuts.
4) Once the datacards are made, Set WhetherFitBDTandMakeCards to False and uncomment the following:
[executeDataCards_onCondor](https://github.com/T3MuAnalysisTools/Projections/blob/cb5efdac12f291e91bc00f29275c390868a90af3/makeYield/ZTT/makeCards.py#L537) 
5) Run makeCards.py again to submit at all the datacards to condor.
6) Comment the executeDataCards_onCondor and uncomment the following:
[ReadAndCopyMinimumBDTCard](https://github.com/T3MuAnalysisTools/Projections/blob/cb5efdac12f291e91bc00f29275c390868a90af3/makeYield/ZTT/makeCards.py#L538)
7) Run makeCards.py again to read all the limits, finds the minimum BDT cut and copies the datacards to the `makeYield/ZTT/datacards_modified` folder.
8) To get the datacards for the 'combined' category, set WhetherFitBDTandMakeCards to True, uncomment executeDataCards_onCondor and ReadAndCopyMinimumBDTCard and run makeCards.py with the category name set to 'combined'. Run all the other categories first before setting the `categories` to `['combined']`.

Once the datacards are made, run the following to submit lumi scan to condor and to plot the result:

~~~
./readLimits.py
~~~

Help for running it is:

~~~
./readLimits.py -h
~~~

For the ZTT category, `./readLimits.py -c ZTT -m A -r run` will submit the datacards (in the `makeYield/ZTT/datacards_modified` folder) to condor. Then `./readLimits.py -c ZTT -m A -r plot` plots the result.

Vladimir: Arun, here it is not clear how I produce the projection curve only for, lets say, tau_e. 
