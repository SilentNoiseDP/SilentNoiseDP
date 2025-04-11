# Usage (Frequency Estimation)

**(1) Install**

Install C/C++ codes as follows.
```
$ cd Frequency estimation/cpp
$ make
$ cd ../
```

**(2) Download datasets**

Download the [Foursquare dataset](https://sites.google.com/site/yangdingqi/home/foursquare-dataset) and the [USCensus dataset](https://archive.ics.uci.edu/dataset/116/us+census+data+1990).

Run "ReadFoursquare.exe" and "ReadUSCensus.exe" and place the results in "data/Foursquare_TIST15_MH/" and "data/USCensus1990/", respectively.

**(2) Get numerical bounds in shuffle protocols**

Place [ml-shuffling-amplification](https://github.com/apple/ml-shuffling-amplification) in "python/", run "NumericalBound.py", and place "numerical-bound_n359094_d10-8.csv" and "numerical-bound_n2458285_d10-8" in "data/Foursquare_TIST15_MH/" and "data/USCensus1990/", respectively.

**(4) Run protocols**

Run codes in "cpp/".
