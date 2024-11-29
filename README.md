# Projections

To create datacards in the `makeYield/ZTT`, `makeYield/ZTT` and `makeYield/W` directories, run the `python makeCards.py` file in these directories to create the datacards in the `makeYield/*/datacards_modified` directory.

For the ZTT limit projections, a limit scan needs to be run at each luminosity wrt BDT cuts. To create the datacards for it, define the categories to run here:

~~~
[categories = ['taue','taumu','tauhA','tauhB','all']](https://github.com/T3MuAnalysisTools/Projections/blob/cb5efdac12f291e91bc00f29275c390868a90af3/makeYield/ZTT/makeCards.py#L486)
~~~
Set WhetherFitBDTandMakeCards to True and run makeCards.py to create the datacards for the limit scan wrt bdt cuts. Once the datacards are made, Set it to false and uncomment the following:

~~~
[executeDataCards_onCondor(lumi,categories,False,bdt_points)](https://github.com/T3MuAnalysisTools/Projections/blob/cb5efdac12f291e91bc00f29275c390868a90af3/makeYield/ZTT/makeCards.py#L537)
[ReadAndCopyMinimumBDTCard(lumi,categories,False,bdt_points)](https://github.com/T3MuAnalysisTools/Projections/blob/cb5efdac12f291e91bc00f29275c390868a90af3/makeYield/ZTT/makeCards.py#L538)
~~~

This submits the lumi limit runs to condor and then reads the limits, finds the minimum BDT cut and copies the datacards to the `makeYield/*/datacards_modified` folder.

This isn't needed for the 'combined' category. For that category, run all the other categories first before setting the `categories` to `['combined']` and run with WhetherFitBDTandMakeCards set to True.

Once the datacards are made, run the following to submit lumi scan to condor and to plot the result:

~~~
./readLimits.py
~~~

Help for running it is:

~~~
./readLimits.py -h
~~~


