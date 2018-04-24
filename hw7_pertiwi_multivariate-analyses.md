Homework 7 (Multivariate Analyses)
================
Cininta Pertiwi

------------------------------------------------------------------------

### <span style="color:blue">Notes on approach</span>

------------------------------------------------------------------------

**--Data wrangling--**

-   The problem have asked to *"develop a predictive model of the mean monthly Chl-a concentration in the California Bay Delta using other mean monthly water quality variables"*. I interprested the 'mean monthly' values to be the mean value for a single month meaning that for each month of each year, there is one value which is the mean of all the values in that month. In this dataset, for example, for each variable there is one value for January 2005, one value for February 2005, and so on. This mean monthly value, therefore, was calculated by taking the mean values of all stations since the measurement for each station is collected monthly.
-   Since the predictive model is to predict Chl-a concentration using other mean monthly water quality variables, all other variables (columns) that are not water quality variables were removed. This leaves a dataset that only contains the response variable (mean monthly Chl-a) and the explanatory/predictor variables (mean monthly water quality variables). This dataset is then used to prepare datasets for problems 1 and 2.

**--For problem 1--**

-   A stepwise regression was used as the first step in selecting a model. Even though stepwise regression is somewhat a 'black box' approach, it is used because I did not want to use expert judgment (since I do not have sufficient expertise in water quality assessment) so stepwise regression was used to reduce the number of variables.
-   To explore parsimony, however, I did additional multiple correlations (pairwise simple correlation) to further reduce the number of variables in the model.
-   To compare the models, I mostly used the combination of the following values: adjusted R-squared, test statistic, AIC, BIC, and p-value.

**--For problem 2--**

-   The grouping of values into dry and wet season was done according [to this article](https://www.water.ca.gov/LegacyFiles/floodmgmt/hafoo/csc/docs/CA_Precipitation_2pager.pdf) from the CA Department of Water Resources. The **wet season** is categorized as months where 90% of annual precipitation occurs which is from October to April. Meanwhile the **dry season** is categorized as the remaining months which is from May to September.

------------------------------------------------------------------------

### <span style="color:blue">Data Wrangling</span>

------------------------------------------------------------------------

``` r
# load packages

library(tidyverse)
library(corrplot)
library(lubridate)
library(broom)
library(reshape2)
```

``` r
# load and check dataset

wq <- read_csv("BayDeltaWQ.csv",
               col_names = TRUE,
               na = c("NA", "n/p", "n/a"),
               guess_max = 30000)
wq
```

    ## # A tibble: 30,520 x 103
    ##       X1 SampleDate StationCode Depth `1% Light Depth` `Conductance (EC)`
    ##    <int> <date>     <chr>       <int>            <dbl>              <dbl>
    ##  1     1 1975-01-07 D11             3             1.77                235
    ##  2     2 1975-01-07 D15             3             1.77                230
    ##  3     3 1975-01-07 D16             3             1.54                185
    ##  4     4 1975-01-07 D19             3             1.64                209
    ##  5     5 1975-01-07 D22             3             1.77                191
    ##  6     6 1975-01-07 D24             3             1.94                170
    ##  7     7 1975-01-07 D26             3             1.64                193
    ##  8     8 1975-01-07 D4              3             1.48                344
    ##  9     9 1975-01-08 D10             3             1.57               2500
    ## 10    10 1975-01-08 D10            37            NA                  4510
    ## # ... with 30,510 more rows, and 97 more variables: SiteDepth <dbl>,
    ## #   `Field Notes` <dbl>, Fluorescence <dbl>, Latitude <dbl>, Longitude
    ## #   <dbl>, Oxygen <dbl>, pH <dbl>, `Secchi Depth` <dbl>, Temperature
    ## #   <dbl>, `Tide Stage` <int>, `Tide Time` <chr>, Turbidity <dbl>,
    ## #   `Weather Observations` <chr>, `Wind Direction` <int>, `Wind Velocity`
    ## #   <int>, VarType <chr>, Alachlor <chr>, Aldrin <chr>, `Ammonia
    ## #   (Dissolved)` <dbl>, `Ammonia (Total)` <dbl>, `Arsenic (Dissolved)`
    ## #   <int>, `Arsenic (Total)` <int>, `Atra+Simazine (Atrazine & Simazine)`
    ## #   <dbl>, Atrazine <chr>, BHC <dbl>, `BHC-alpha` <chr>, `BHC-beta` <chr>,
    ## #   `BHC-delta` <chr>, `BHC-gamma (Lindane)` <chr>, `Biochemical Oxygen
    ## #   Demand (BOD)` <dbl>, `Cadmium (Dissolved)` <int>, `Cadmium (Total)`
    ## #   <int>, `Calcium (Dissolved)` <int>, Captafol <dbl>, Captan <chr>,
    ## #   Chlordane <chr>, `Chloride (Dissolved)` <dbl>, `Chloride (Total)`
    ## #   <dbl>, `Chlorophyll a` <dbl>, Chlorothalonil <chr>, Chlorpropham
    ## #   <chr>, Chlorpyrifos <chr>, `Chromium (Dissolved)` <int>, `Chromium
    ## #   (Total)` <int>, `Copper (Dissolved)` <int>, `Copper (Total)` <int>,
    ## #   `Dacthal (DCPA)` <dbl>, Dichloran <chr>, Dicofol <chr>, Dieldrin
    ## #   <dbl>, Diuron <dbl>, `Endosulfan (mixed isomers)` <dbl>,
    ## #   `Endosulfan-I` <chr>, `Endosulfan-II` <chr>, Endrin <dbl>, `Endrin
    ## #   aldehyde` <chr>, Heptachlor <chr>, `Iron (Dissolved)` <int>, `Iron
    ## #   (Total)` <int>, `Kjeldahl Nitrogen (Total)` <dbl>, `Lead (Dissolved)`
    ## #   <int>, `Lead (Total)` <int>, `Manganese (Dissolved)` <int>, `Manganese
    ## #   (Total)` <int>, `Mercury (Total)` <dbl>, Methoxychlor <dbl>, `Nitrate
    ## #   (Dissolved)` <int>, `Nitrite (Dissolved)` <dbl>, `Nitrite + Nitrate
    ## #   (Dissolved)` <dbl>, `Organic Carbon (Dissolved)` <dbl>, `Organic
    ## #   Carbon (Total)` <dbl>, `Organic Nitrogen (Dissolved)` <dbl>, `Organic
    ## #   Nitrogen (Total)` <dbl>, `Ortho-phosphate (Dissolved)` <dbl>,
    ## #   `p,p'-DDD` <chr>, `p,p'-DDE` <dbl>, `p,p'-DDT` <chr>, `PCB-1016`
    ## #   <chr>, `PCB-1221` <chr>, `PCB-1232` <chr>, `PCB-1242` <chr>,
    ## #   `PCB-1248` <chr>, `PCB-1254` <chr>, `PCB-1260` <chr>,
    ## #   `Pentachloronitrobenzene (PCNB)` <dbl>, `Pheophytin a` <dbl>,
    ## #   `Phosphorus (Total)` <dbl>, `Silica (SiO2) (Dissolved)` <dbl>,
    ## #   Simazine <chr>, `Solids (Total Dissolved)` <dbl>, `Solids (Total
    ## #   Suspended)` <dbl>, `Solids (Volatile Suspended)` <dbl>, Thiobencarb
    ## #   <chr>, Toxaphene <chr>, `Unknown hydrocarbon` <dbl>, `Zinc
    ## #   (Dissolved)` <int>, `Zinc (Total)` <int>

``` r
# pull out month and year from SampleDate--
# to allow for grouping by year and month

wq_bymonth <- wq %>%
  mutate(Year = lubridate::year(SampleDate)) %>%    # add Year column
  mutate(Month = lubridate::month(SampleDate)) %>%  # add Month column
  group_by(Year, Month) %>%                         # # group rows by year and month
  select(Month, everything()) %>%                   # move Month column as 1st column
  select(Year, everything())                        # move Year column as 1st column
wq_bymonth
```

    ## # A tibble: 30,520 x 105
    ## # Groups: Year, Month [455]
    ##     Year Month    X1 SampleDate StationCode Depth `1% Light Depth`
    ##    <dbl> <dbl> <int> <date>     <chr>       <int>            <dbl>
    ##  1  1975  1.00     1 1975-01-07 D11             3             1.77
    ##  2  1975  1.00     2 1975-01-07 D15             3             1.77
    ##  3  1975  1.00     3 1975-01-07 D16             3             1.54
    ##  4  1975  1.00     4 1975-01-07 D19             3             1.64
    ##  5  1975  1.00     5 1975-01-07 D22             3             1.77
    ##  6  1975  1.00     6 1975-01-07 D24             3             1.94
    ##  7  1975  1.00     7 1975-01-07 D26             3             1.64
    ##  8  1975  1.00     8 1975-01-07 D4              3             1.48
    ##  9  1975  1.00     9 1975-01-08 D10             3             1.57
    ## 10  1975  1.00    10 1975-01-08 D10            37            NA   
    ## # ... with 30,510 more rows, and 98 more variables: `Conductance (EC)`
    ## #   <dbl>, SiteDepth <dbl>, `Field Notes` <dbl>, Fluorescence <dbl>,
    ## #   Latitude <dbl>, Longitude <dbl>, Oxygen <dbl>, pH <dbl>, `Secchi
    ## #   Depth` <dbl>, Temperature <dbl>, `Tide Stage` <int>, `Tide Time`
    ## #   <chr>, Turbidity <dbl>, `Weather Observations` <chr>, `Wind Direction`
    ## #   <int>, `Wind Velocity` <int>, VarType <chr>, Alachlor <chr>, Aldrin
    ## #   <chr>, `Ammonia (Dissolved)` <dbl>, `Ammonia (Total)` <dbl>, `Arsenic
    ## #   (Dissolved)` <int>, `Arsenic (Total)` <int>, `Atra+Simazine (Atrazine
    ## #   & Simazine)` <dbl>, Atrazine <chr>, BHC <dbl>, `BHC-alpha` <chr>,
    ## #   `BHC-beta` <chr>, `BHC-delta` <chr>, `BHC-gamma (Lindane)` <chr>,
    ## #   `Biochemical Oxygen Demand (BOD)` <dbl>, `Cadmium (Dissolved)` <int>,
    ## #   `Cadmium (Total)` <int>, `Calcium (Dissolved)` <int>, Captafol <dbl>,
    ## #   Captan <chr>, Chlordane <chr>, `Chloride (Dissolved)` <dbl>, `Chloride
    ## #   (Total)` <dbl>, `Chlorophyll a` <dbl>, Chlorothalonil <chr>,
    ## #   Chlorpropham <chr>, Chlorpyrifos <chr>, `Chromium (Dissolved)` <int>,
    ## #   `Chromium (Total)` <int>, `Copper (Dissolved)` <int>, `Copper (Total)`
    ## #   <int>, `Dacthal (DCPA)` <dbl>, Dichloran <chr>, Dicofol <chr>,
    ## #   Dieldrin <dbl>, Diuron <dbl>, `Endosulfan (mixed isomers)` <dbl>,
    ## #   `Endosulfan-I` <chr>, `Endosulfan-II` <chr>, Endrin <dbl>, `Endrin
    ## #   aldehyde` <chr>, Heptachlor <chr>, `Iron (Dissolved)` <int>, `Iron
    ## #   (Total)` <int>, `Kjeldahl Nitrogen (Total)` <dbl>, `Lead (Dissolved)`
    ## #   <int>, `Lead (Total)` <int>, `Manganese (Dissolved)` <int>, `Manganese
    ## #   (Total)` <int>, `Mercury (Total)` <dbl>, Methoxychlor <dbl>, `Nitrate
    ## #   (Dissolved)` <int>, `Nitrite (Dissolved)` <dbl>, `Nitrite + Nitrate
    ## #   (Dissolved)` <dbl>, `Organic Carbon (Dissolved)` <dbl>, `Organic
    ## #   Carbon (Total)` <dbl>, `Organic Nitrogen (Dissolved)` <dbl>, `Organic
    ## #   Nitrogen (Total)` <dbl>, `Ortho-phosphate (Dissolved)` <dbl>,
    ## #   `p,p'-DDD` <chr>, `p,p'-DDE` <dbl>, `p,p'-DDT` <chr>, `PCB-1016`
    ## #   <chr>, `PCB-1221` <chr>, `PCB-1232` <chr>, `PCB-1242` <chr>,
    ## #   `PCB-1248` <chr>, `PCB-1254` <chr>, `PCB-1260` <chr>,
    ## #   `Pentachloronitrobenzene (PCNB)` <dbl>, `Pheophytin a` <dbl>,
    ## #   `Phosphorus (Total)` <dbl>, `Silica (SiO2) (Dissolved)` <dbl>,
    ## #   Simazine <chr>, `Solids (Total Dissolved)` <dbl>, `Solids (Total
    ## #   Suspended)` <dbl>, `Solids (Volatile Suspended)` <dbl>, Thiobencarb
    ## #   <chr>, Toxaphene <chr>, `Unknown hydrocarbon` <dbl>, `Zinc
    ## #   (Dissolved)` <int>, `Zinc (Total)` <int>

``` r
# filter rows to select only those sampled from--
# October 2004 to September 2012 (water years 2005-2012)

wq_bymonth <- wq_bymonth %>%
  #group_by(Year, Month) %>%               # group rows by year and month
  filter(Year > 2003) %>%                 # select observations after 2003
  filter(!(Year==2004 && Month<10)) %>%   # take out January-September 2004
  filter(!(Year==2012 && Month>9))        # take out October-December 2012
wq_bymonth
```

    ## # A tibble: 4,723 x 105
    ## # Groups: Year, Month [96]
    ##     Year Month    X1 SampleDate StationCode Depth `1% Light Depth`
    ##    <dbl> <dbl> <int> <date>     <chr>       <int>            <dbl>
    ##  1  2004  10.0 12573 2004-10-18 D16             3               NA
    ##  2  2004  10.0 12574 2004-10-18 D19             3               NA
    ##  3  2004  10.0 12575 2004-10-18 D26             3               NA
    ##  4  2004  10.0 12576 2004-10-18 D28A            3               NA
    ##  5  2004  10.0 12577 2004-10-19 MD10A           3               NA
    ##  6  2004  10.0 12578 2004-10-19 P8              3               NA
    ##  7  2004  10.0 12579 2004-10-20 C10             3               NA
    ##  8  2004  10.0 12580 2004-10-20 C3              3               NA
    ##  9  2004  10.0 12581 2004-10-20 D7              3               NA
    ## 10  2004  10.0 12582 2004-10-20 NZ032           3               NA
    ## # ... with 4,713 more rows, and 98 more variables: `Conductance (EC)`
    ## #   <dbl>, SiteDepth <dbl>, `Field Notes` <dbl>, Fluorescence <dbl>,
    ## #   Latitude <dbl>, Longitude <dbl>, Oxygen <dbl>, pH <dbl>, `Secchi
    ## #   Depth` <dbl>, Temperature <dbl>, `Tide Stage` <int>, `Tide Time`
    ## #   <chr>, Turbidity <dbl>, `Weather Observations` <chr>, `Wind Direction`
    ## #   <int>, `Wind Velocity` <int>, VarType <chr>, Alachlor <chr>, Aldrin
    ## #   <chr>, `Ammonia (Dissolved)` <dbl>, `Ammonia (Total)` <dbl>, `Arsenic
    ## #   (Dissolved)` <int>, `Arsenic (Total)` <int>, `Atra+Simazine (Atrazine
    ## #   & Simazine)` <dbl>, Atrazine <chr>, BHC <dbl>, `BHC-alpha` <chr>,
    ## #   `BHC-beta` <chr>, `BHC-delta` <chr>, `BHC-gamma (Lindane)` <chr>,
    ## #   `Biochemical Oxygen Demand (BOD)` <dbl>, `Cadmium (Dissolved)` <int>,
    ## #   `Cadmium (Total)` <int>, `Calcium (Dissolved)` <int>, Captafol <dbl>,
    ## #   Captan <chr>, Chlordane <chr>, `Chloride (Dissolved)` <dbl>, `Chloride
    ## #   (Total)` <dbl>, `Chlorophyll a` <dbl>, Chlorothalonil <chr>,
    ## #   Chlorpropham <chr>, Chlorpyrifos <chr>, `Chromium (Dissolved)` <int>,
    ## #   `Chromium (Total)` <int>, `Copper (Dissolved)` <int>, `Copper (Total)`
    ## #   <int>, `Dacthal (DCPA)` <dbl>, Dichloran <chr>, Dicofol <chr>,
    ## #   Dieldrin <dbl>, Diuron <dbl>, `Endosulfan (mixed isomers)` <dbl>,
    ## #   `Endosulfan-I` <chr>, `Endosulfan-II` <chr>, Endrin <dbl>, `Endrin
    ## #   aldehyde` <chr>, Heptachlor <chr>, `Iron (Dissolved)` <int>, `Iron
    ## #   (Total)` <int>, `Kjeldahl Nitrogen (Total)` <dbl>, `Lead (Dissolved)`
    ## #   <int>, `Lead (Total)` <int>, `Manganese (Dissolved)` <int>, `Manganese
    ## #   (Total)` <int>, `Mercury (Total)` <dbl>, Methoxychlor <dbl>, `Nitrate
    ## #   (Dissolved)` <int>, `Nitrite (Dissolved)` <dbl>, `Nitrite + Nitrate
    ## #   (Dissolved)` <dbl>, `Organic Carbon (Dissolved)` <dbl>, `Organic
    ## #   Carbon (Total)` <dbl>, `Organic Nitrogen (Dissolved)` <dbl>, `Organic
    ## #   Nitrogen (Total)` <dbl>, `Ortho-phosphate (Dissolved)` <dbl>,
    ## #   `p,p'-DDD` <chr>, `p,p'-DDE` <dbl>, `p,p'-DDT` <chr>, `PCB-1016`
    ## #   <chr>, `PCB-1221` <chr>, `PCB-1232` <chr>, `PCB-1242` <chr>,
    ## #   `PCB-1248` <chr>, `PCB-1254` <chr>, `PCB-1260` <chr>,
    ## #   `Pentachloronitrobenzene (PCNB)` <dbl>, `Pheophytin a` <dbl>,
    ## #   `Phosphorus (Total)` <dbl>, `Silica (SiO2) (Dissolved)` <dbl>,
    ## #   Simazine <chr>, `Solids (Total Dissolved)` <dbl>, `Solids (Total
    ## #   Suspended)` <dbl>, `Solids (Volatile Suspended)` <dbl>, Thiobencarb
    ## #   <chr>, Toxaphene <chr>, `Unknown hydrocarbon` <dbl>, `Zinc
    ## #   (Dissolved)` <int>, `Zinc (Total)` <int>

``` r
# summarize the column values to the mean monthly--
# so that each month in each year has a single value--
# which is the mean of all values in that month.
# this will result in 96 rows where each row is the--
# mean value for that month in that particular year--
# (12 months x 8 years = 96 months)

wq_meanmonthly <- wq_bymonth %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
       # take the mean values of each month of--
       # each year only if the column type is a numeric
  select(-X1) %>%
       # take out index column (not useful)
  select(`Chlorophyll a`, everything())
       # move Chl-a column as first column
wq_meanmonthly
```

    ## # A tibble: 96 x 69
    ## # Groups: Year [9]
    ##    `Chlorophyll a`  Year Month Depth `1% Light Depth` `Conductance (EC)`
    ##              <dbl> <dbl> <dbl> <dbl>            <dbl>              <dbl>
    ##  1           2.43   2004 10.0   3.00              NaN              12951
    ##  2           1.63   2004 11.0   3.00              NaN              12791
    ##  3           0.462  2004 12.0   3.00              NaN              11360
    ##  4           0.923  2005  1.00  3.00              NaN               6480
    ##  5           1.51   2005  2.00  3.00              NaN               6788
    ##  6           4.55   2005  3.00  3.00              NaN               4981
    ##  7           4.50   2005  4.00  3.00              NaN               5753
    ##  8           4.37   2005  5.00  3.00              NaN               3545
    ##  9           3.53   2005  6.00  3.00              NaN               4873
    ## 10           5.54   2005  7.00  3.00              NaN               6708
    ## # ... with 86 more rows, and 63 more variables: SiteDepth <dbl>, `Field
    ## #   Notes` <dbl>, Fluorescence <dbl>, Latitude <dbl>, Longitude <dbl>,
    ## #   Oxygen <dbl>, pH <dbl>, `Secchi Depth` <dbl>, Temperature <dbl>, `Tide
    ## #   Stage` <dbl>, Turbidity <dbl>, `Wind Direction` <dbl>, `Wind Velocity`
    ## #   <dbl>, `Ammonia (Dissolved)` <dbl>, `Ammonia (Total)` <dbl>, `Arsenic
    ## #   (Dissolved)` <dbl>, `Arsenic (Total)` <dbl>, `Atra+Simazine (Atrazine
    ## #   & Simazine)` <dbl>, BHC <dbl>, `Biochemical Oxygen Demand (BOD)`
    ## #   <dbl>, `Cadmium (Dissolved)` <dbl>, `Cadmium (Total)` <dbl>, `Calcium
    ## #   (Dissolved)` <dbl>, Captafol <dbl>, `Chloride (Dissolved)` <dbl>,
    ## #   `Chloride (Total)` <dbl>, `Chromium (Dissolved)` <dbl>, `Chromium
    ## #   (Total)` <dbl>, `Copper (Dissolved)` <dbl>, `Copper (Total)` <dbl>,
    ## #   `Dacthal (DCPA)` <dbl>, Dieldrin <dbl>, Diuron <dbl>, `Endosulfan
    ## #   (mixed isomers)` <dbl>, Endrin <dbl>, `Iron (Dissolved)` <dbl>, `Iron
    ## #   (Total)` <dbl>, `Kjeldahl Nitrogen (Total)` <dbl>, `Lead (Dissolved)`
    ## #   <dbl>, `Lead (Total)` <dbl>, `Manganese (Dissolved)` <dbl>, `Manganese
    ## #   (Total)` <dbl>, `Mercury (Total)` <dbl>, Methoxychlor <dbl>, `Nitrate
    ## #   (Dissolved)` <dbl>, `Nitrite (Dissolved)` <dbl>, `Nitrite + Nitrate
    ## #   (Dissolved)` <dbl>, `Organic Carbon (Dissolved)` <dbl>, `Organic
    ## #   Carbon (Total)` <dbl>, `Organic Nitrogen (Dissolved)` <dbl>, `Organic
    ## #   Nitrogen (Total)` <dbl>, `Ortho-phosphate (Dissolved)` <dbl>,
    ## #   `p,p'-DDE` <dbl>, `Pentachloronitrobenzene (PCNB)` <dbl>, `Pheophytin
    ## #   a` <dbl>, `Phosphorus (Total)` <dbl>, `Silica (SiO2) (Dissolved)`
    ## #   <dbl>, `Solids (Total Dissolved)` <dbl>, `Solids (Total Suspended)`
    ## #   <dbl>, `Solids (Volatile Suspended)` <dbl>, `Unknown hydrocarbon`
    ## #   <dbl>, `Zinc (Dissolved)` <dbl>, `Zinc (Total)` <dbl>

This `wq_meanmonthly` dataset is what will be used in preparing datasets for problems 1 and 2.

------------------------------------------------------------------------

### <span style="color:blue">1| Predictive Model</span>

------------------------------------------------------------------------

#### **--Data preparation**

``` r
# remove NaN values

wq_meanmonthly_mod <- wq_meanmonthly %>%
  select_if(~sum(!is.na(.)) > 0) %>%
       # remove columns with values that are all NaN
  select_if(~sum(!is.na(.)) == nrow(wq_meanmonthly))
       # remove columns which have NaN values
wq_meanmonthly_mod
```

    ## # A tibble: 96 x 23
    ## # Groups: Year [9]
    ##     Year `Chlorophyll a` Month Depth `Conductance (EC)` SiteDepth
    ##    <dbl>           <dbl> <dbl> <dbl>              <dbl>     <dbl>
    ##  1  2004           2.43  10.0   3.00              12951      33.4
    ##  2  2004           1.63  11.0   3.00              12791      34.6
    ##  3  2004           0.462 12.0   3.00              11360      34.3
    ##  4  2005           0.923  1.00  3.00               6480      36.5
    ##  5  2005           1.51   2.00  3.00               6788      35.0
    ##  6  2005           4.55   3.00  3.00               4981      35.1
    ##  7  2005           4.50   4.00  3.00               5753      36.5
    ##  8  2005           4.37   5.00  3.00               3545      33.0
    ##  9  2005           3.53   6.00  3.00               4873      32.7
    ## 10  2005           5.54   7.00  3.00               6708      32.9
    ## # ... with 86 more rows, and 17 more variables: Fluorescence <dbl>, Oxygen
    ## #   <dbl>, `Secchi Depth` <dbl>, Temperature <dbl>, Turbidity <dbl>,
    ## #   `Ammonia (Dissolved)` <dbl>, `Chloride (Dissolved)` <dbl>, `Kjeldahl
    ## #   Nitrogen (Total)` <dbl>, `Nitrite + Nitrate (Dissolved)` <dbl>,
    ## #   `Organic Nitrogen (Dissolved)` <dbl>, `Ortho-phosphate (Dissolved)`
    ## #   <dbl>, `Pheophytin a` <dbl>, `Phosphorus (Total)` <dbl>, `Silica
    ## #   (SiO2) (Dissolved)` <dbl>, `Solids (Total Dissolved)` <dbl>, `Solids
    ## #   (Total Suspended)` <dbl>, `Solids (Volatile Suspended)` <dbl>

``` r
# remove variables that are not water quality variables

wq_meanmonthly_mod <- wq_meanmonthly_mod %>%
  ungroup() %>%
  select(-Year, -Month, -Depth)
wq_meanmonthly_mod
```

    ## # A tibble: 96 x 20
    ##    `Chlorophyll a` `Conductance (EC)` SiteDepth Fluorescence Oxygen
    ##              <dbl>              <dbl>     <dbl>        <dbl>  <dbl>
    ##  1           2.43               12951      33.4         1.35   8.01
    ##  2           1.63               12791      34.6         1.33   8.32
    ##  3           0.462              11360      34.3         1.28   8.93
    ##  4           0.923               6480      36.5         1.77   9.43
    ##  5           1.51                6788      35.0         1.37   8.87
    ##  6           4.55                4981      35.1         2.27   8.84
    ##  7           4.50                5753      36.5         2.52   9.17
    ##  8           4.37                3545      33.0         2.12   8.44
    ##  9           3.53                4873      32.7         1.84   7.98
    ## 10           5.54                6708      32.9         1.82   7.71
    ## # ... with 86 more rows, and 15 more variables: `Secchi Depth` <dbl>,
    ## #   Temperature <dbl>, Turbidity <dbl>, `Ammonia (Dissolved)` <dbl>,
    ## #   `Chloride (Dissolved)` <dbl>, `Kjeldahl Nitrogen (Total)` <dbl>,
    ## #   `Nitrite + Nitrate (Dissolved)` <dbl>, `Organic Nitrogen (Dissolved)`
    ## #   <dbl>, `Ortho-phosphate (Dissolved)` <dbl>, `Pheophytin a` <dbl>,
    ## #   `Phosphorus (Total)` <dbl>, `Silica (SiO2) (Dissolved)` <dbl>, `Solids
    ## #   (Total Dissolved)` <dbl>, `Solids (Total Suspended)` <dbl>, `Solids
    ## #   (Volatile Suspended)` <dbl>

Now I will use the `wq_meanmonthly_mod` dataset for the model selection.

#### **--Model selection**

``` r
# use stepwise to reduce the number of variables

wq_step <- step(lm(`Chlorophyll a` ~ ., data = wq_meanmonthly_mod), trace = 0)
     # stepwise regression; trace = 0 returns only the final model of the stepwise
wq_step
```

    ## 
    ## Call:
    ## lm(formula = `Chlorophyll a` ~ `Conductance (EC)` + Oxygen + 
    ##     Temperature + `Ammonia (Dissolved)` + `Kjeldahl Nitrogen (Total)` + 
    ##     `Organic Nitrogen (Dissolved)` + `Pheophytin a` + `Solids (Total Dissolved)`, 
    ##     data = wq_meanmonthly_mod)
    ## 
    ## Coefficients:
    ##                    (Intercept)              `Conductance (EC)`  
    ##                     -2.818e+01                      -3.190e-04  
    ##                         Oxygen                     Temperature  
    ##                      2.288e+00                       5.800e-01  
    ##          `Ammonia (Dissolved)`     `Kjeldahl Nitrogen (Total)`  
    ##                     -8.995e+00                       6.765e+00  
    ## `Organic Nitrogen (Dissolved)`                  `Pheophytin a`  
    ##                     -4.271e+00                       8.332e-01  
    ##     `Solids (Total Dissolved)`  
    ##                      5.821e-04

``` r
# the stepwise regression resulted in a model with--
# 8 predictor variables instead of the initial 19 variables

# make a new dataset containing only the 8 selected variables

wq_model <- select(wq_meanmonthly_mod,
                   `Chlorophyll a`,
                   `Conductance (EC)`,
                   Oxygen,
                   Temperature,
                   `Ammonia (Dissolved)`,
                   `Kjeldahl Nitrogen (Total)`,
                   `Organic Nitrogen (Dissolved)`,
                   `Pheophytin a`,
                   `Solids (Total Dissolved)`)
wq_model
```

    ## # A tibble: 96 x 9
    ##    `Chlorophyll a` `Conductance (EC~ Oxygen Temperature `Ammonia (Dissolv~
    ##              <dbl>             <dbl>  <dbl>       <dbl>              <dbl>
    ##  1           2.43              12951   8.01       17.4              0.154 
    ##  2           1.63              12791   8.32       13.8              0.178 
    ##  3           0.462             11360   8.93       10.9              0.236 
    ##  4           0.923              6480   9.43        8.96             0.156 
    ##  5           1.51               6788   8.87       12.4              0.175 
    ##  6           4.55               4981   8.84       15.6              0.118 
    ##  7           4.50               5753   9.17       15.6              0.0927
    ##  8           4.37               3545   8.44       18.5              0.0736
    ##  9           3.53               4873   7.98       20.1              0.0562
    ## 10           5.54               6708   7.71       22.7              0.0836
    ## # ... with 86 more rows, and 4 more variables: `Kjeldahl Nitrogen (Total)`
    ## #   <dbl>, `Organic Nitrogen (Dissolved)` <dbl>, `Pheophytin a` <dbl>,
    ## #   `Solids (Total Dissolved)` <dbl>

``` r
# let's check if the model with 8 predictor variables is actually--
# better than with 19 predictor variables

lm_19v <- lm(`Chlorophyll a` ~ ., data = wq_meanmonthly_mod)
lm_8v <- lm(`Chlorophyll a` ~ ., data = wq_model)

lms <- list(w19v = lm_19v, w8v = lm_8v)
lms.stats <- mapply(glance, lms)
colnames(lms.stats) <- names(lms)
lms.stats
```

    ##               w19v        w8v         
    ## r.squared     0.6170207   0.5929971   
    ## adj.r.squared 0.5212759   0.5555716   
    ## sigma         1.594462    1.536287    
    ## statistic     6.444429    15.84471    
    ## p.value       1.66057e-09 3.569059e-14
    ## df            20          9           
    ## logLik        -169.7921   -172.7124   
    ## AIC           381.5842    365.4247    
    ## BIC           435.4355    391.0682    
    ## deviance      193.2155    205.3356    
    ## df.residual   76          87

Although maybe not by much, the model from the stepwise regression (with 8 predictor variables) has a higher adj-r-square and test statistic with lower p-value, AIC, and BIC. From these values, the model with 8 predictor variables do seem to be better.

#### **--Parsimony?**

Having 8 predictor variables still sounds like a lot of variables. Let's see if variable reduction can still be done by checking potential correlations between the predictor variables.

``` r
# --- (1) let's visualize the relationships between Chl-a and the 8 variables

plot(wq_model, pch=16, col="blue", cex = 0.5, main="Model: Chl-a ~ 8 variables")
```

![](hw7_pertiwi_multivariate-analyses_files/figure-markdown_github/unnamed-chunk-11-1.png)

From the plot, it looks like the following variables are potentially highly correlated:

-   Oxyen and Temperature
-   EC and Solids
-   Ammonia and Kjeldahl Nitrogen
-   Kjeldahl Nitrogen and Organic Nitrogen

``` r
# --- (2) let's check with correlation plot

mycor <- cor(wq_model)

cex.before <- par("cex")
par(cex = 0.7)

corrplot(mycor, method = "number", tl.cex = 1/par("cex"),
         cl.cex = 1/par("cex"))
```

![](hw7_pertiwi_multivariate-analyses_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
par(cex = cex.before)
```

The correlation plot verifies that the following variables are highly correlated:

-   Oxyen x Temperature
-   EC x Solids
-   Ammonia x Kjeldahl Nitrogen
-   Kjeldahl Nitrogen x Organic Nitrogen
-   Temperature x Ammonia

Should some variables be removed? But which ones?

``` r
# --- (3) check regression summary to determine which variable to remove

chla_all <- lm(`Chlorophyll a` ~ ., data = wq_model)
summary(chla_all)
```

    ## 
    ## Call:
    ## lm(formula = `Chlorophyll a` ~ ., data = wq_model)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5837 -0.7604 -0.1729  0.5220  8.1688 
    ## 
    ## Coefficients:
    ##                                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                    -2.818e+01  7.968e+00  -3.536 0.000653 ***
    ## `Conductance (EC)`             -3.190e-04  1.838e-04  -1.736 0.086158 .  
    ## Oxygen                          2.288e+00  6.500e-01   3.519 0.000690 ***
    ## Temperature                     5.800e-01  1.417e-01   4.092 9.52e-05 ***
    ## `Ammonia (Dissolved)`          -8.995e+00  5.672e+00  -1.586 0.116359    
    ## `Kjeldahl Nitrogen (Total)`     6.765e+00  2.741e+00   2.469 0.015522 *  
    ## `Organic Nitrogen (Dissolved)` -4.271e+00  2.764e+00  -1.545 0.125911    
    ## `Pheophytin a`                  8.332e-01  2.971e-01   2.804 0.006219 ** 
    ## `Solids (Total Dissolved)`      5.821e-04  3.324e-04   1.751 0.083415 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.536 on 87 degrees of freedom
    ## Multiple R-squared:  0.593,  Adjusted R-squared:  0.5556 
    ## F-statistic: 15.84 on 8 and 87 DF,  p-value: 3.569e-14

#### **--Model comparison**

Basing on the p-values and the variable correlations above:

-   The variables Oxygen, Organic Nitrogen, and Solids are removed while EC, Ammonia and Kjeldahl Nitrogen can also probably be removed.
-   Variables Temperature and Pheophytin a are both left in the models with 3 and 4 variables due to their low p-value and high correlation with Chl-a but they are not highly correlated with each other.
-   Since Temperature and Ammonia as well as Ammonia with Kjeldahl Nitrogen have relatively high correlations (-0.62 and 0.67 respectively), they are not included together in 3 and 4-variable models.

Let's see if additional reduction improves the model by comparing the 8-variable model and some other models derived from combinations of the variables.

``` r
# --- null
chla_null <- lm(`Chlorophyll a` ~ 1, data = wq_model)

# --- all 8 variables
chla_all <- lm(`Chlorophyll a` ~ ., data = wq_model)
```

``` r
# --- 4 variables: EC, Temp, Kj N, Pheo a

chla_ec.tmp.kjn.phe <- lm(`Chlorophyll a` ~
                              `Conductance (EC)` +
                              Temperature +
                              `Kjeldahl Nitrogen (Total)` +
                              `Pheophytin a`,
                            data = wq_model)
```

``` r
# --- 3 variables: EC, Temp, Pheo a
#                  EC, Ammonia, Pheo a
#                  Temp, Kj N, Pheo a

chla_ec.tmp.phe <- lm(`Chlorophyll a` ~
                         `Conductance (EC)` +
                         Temperature +
                         `Pheophytin a`,
                       data = wq_model)
chla_ec.amm.phe <- lm(`Chlorophyll a` ~
                         `Conductance (EC)` +
                         `Ammonia (Dissolved)` +
                         `Pheophytin a`,
                       data = wq_model)
chla_tmp.kjn.phe <- lm(`Chlorophyll a` ~
                         Temperature +
                         `Kjeldahl Nitrogen (Total)` +
                         `Pheophytin a`,
                       data = wq_model)
```

``` r
# --- 2 variables: Temp and Pheo a
#                  Ammonia and Pheo a
#                  Kj N and Pheo a
#                  Temp and Kj N

chla_tmp.phe <- lm(`Chlorophyll a` ~
                      Temperature +
                      `Pheophytin a`,
                    data = wq_model)
chla_amm.phe <- lm(`Chlorophyll a` ~
                      `Ammonia (Dissolved)` +
                      `Pheophytin a`,
                    data = wq_model)
chla_kjn.phe <- lm(`Chlorophyll a` ~
                      `Kjeldahl Nitrogen (Total)` +
                      `Pheophytin a`,
                    data = wq_model)
chla_tmp.kjn <- lm(`Chlorophyll a` ~
                      Temperature +
                      `Kjeldahl Nitrogen (Total)`,
                    data = wq_model)
```

``` r
lms_compare <- list(null=chla_null,
                    all=chla_all,
                    ec.tmp.kjn.phe=chla_ec.tmp.kjn.phe,
                    ec.tmp.phe=chla_ec.tmp.phe,
                    ec.amm.phe=chla_ec.amm.phe,
                    tmp.kjn.phe=chla_tmp.kjn.phe,
                    tmp.phe=chla_tmp.phe,
                    amm.phe=chla_amm.phe,
                    kjn.phe=chla_kjn.phe,
                    tmp.kjn=chla_tmp.kjn)

lms_compare.stats <- mapply(glance, lms_compare)
colnames(lms_compare.stats) <- names(lms_compare)
lms_compare.stats
```

    ##               null      all          ec.tmp.kjn.phe ec.tmp.phe  
    ## r.squared     0         0.5929971    0.4806624      0.4790654   
    ## adj.r.squared 0         0.5555716    0.4578344      0.4620784   
    ## sigma         2.304473  1.536287     1.696827       1.690173    
    ## statistic     NA        15.84471     21.0558        28.20188    
    ## p.value       NA        3.569059e-14 2.583643e-12   5.064252e-13
    ## df            1         9            5              4           
    ## logLik        -215.8612 -172.7124    -184.4116      -184.559    
    ## AIC           435.7225  365.4247     380.8232       379.118     
    ## BIC           440.8512  391.0682     396.2093       391.9397    
    ## deviance      504.5064  205.3356     262.0091       262.8148    
    ## df.residual   95        87           91             92          
    ##               ec.amm.phe  tmp.kjn.phe  tmp.phe      amm.phe   
    ## r.squared     0.4728065   0.4574591    0.4570602    0.4728041 
    ## adj.r.squared 0.4556154   0.4397676    0.4453841    0.4614666 
    ## sigma         1.700296    1.724867     1.716199     1.691133  
    ## statistic     27.50299    25.85749     39.14486     41.70251  
    ## p.value       8.71719e-13 3.212171e-12 4.633964e-13 1.1795e-13
    ## df            4           4            3            3         
    ## logLik        -185.1322   -186.5096    -186.5449    -185.1325 
    ## AIC           380.2645    383.0193     381.0898     378.2649  
    ## BIC           393.0862    395.841      391.3472     388.5223  
    ## deviance      265.9725    273.7153     273.9166     265.9737  
    ## df.residual   92          92           93           93        
    ##               kjn.phe      tmp.kjn     
    ## r.squared     0.3863068    0.2691238   
    ## adj.r.squared 0.3731091    0.253406    
    ## sigma         1.824599     1.991195    
    ## statistic     29.27076     17.12227    
    ## p.value       1.379546e-10 4.663771e-07
    ## df            3            3           
    ## logLik        -192.4248    -200.8127   
    ## AIC           392.8495     409.6254    
    ## BIC           403.1069     419.8828    
    ## deviance      309.6121     368.7317    
    ## df.residual   93           93

#### **--Selected multiple regression model**

By comparing the models above, ***the model with all 8 variables actually had the highest adj-r-squared value compared to the models with less predictor variables***. However, in the case that fewer variables are available for measurement, I would select the following based on comparison of values of the adj-r-square, AIC/BIC, p-value, and test statistic.

-   **the 3-variable model `Chl-a ~ Ammonia + Pheophytin a`**

This model is able to explain 46% of the variation in Chl-a compared to 55% when using all 8 variables.

Let's compare the residuals and predicted values for the model with 8 predictor variables and the model `Chl-a ~ Ammonia + Pheophtin a`.

``` r
# plot residuals

par(mfrow = c(2, 2))
plot(chla_all, pch = 16, which = 1)
plot(chla_all, pch = 16, which = 2)
plot(chla_amm.phe, pch = 16, which = 1)
plot(chla_amm.phe, pch = 16, which = 2)
```

![](hw7_pertiwi_multivariate-analyses_files/figure-markdown_github/unnamed-chunk-19-1.png)

``` r
par(mfrow = c(1, 1))
```

Top: all 8 variables; Bottom: Chl-a ~ Ammonia + Pheophytin a.

``` r
# plot predicted vs actual

par(mfrow = c(1, 2))
plot(predict(chla_all),wq_model$`Chlorophyll a`,
     xlab="predicted",ylab="actual", main ="Chl-a ~ .")
abline(a=0,b=1)
plot(predict(chla_amm.phe),wq_model$`Chlorophyll a`,
     xlab="predicted",ylab="actual", main ="Chl-a ~ Ammonia + Pheophytin a")
abline(a=0,b=1)
```

![](hw7_pertiwi_multivariate-analyses_files/figure-markdown_github/unnamed-chunk-20-1.png)

``` r
par(mfrow = c(1, 1))
```

#### **--Best predictor**

The selected multiple regression model has 2 variables. Is ther one best predictor for Chl-a concentration? I look at 3 variables: Temperature, Ammonia, and Pheophytin a. These 3 variables have the highest correlation value to Chl-a based on the correlation plot.

``` r
# lm for each of the 3 variables
chla_tmp <- lm(`Chlorophyll a` ~ Temperature, data = wq_model)
chla_amm <- lm(`Chlorophyll a` ~ `Ammonia (Dissolved)`, data = wq_model)
chla_phe <- lm(`Chlorophyll a` ~ `Pheophytin a`, data = wq_model)
```

``` r
# compare lms of the 3 variables
lms_best <- list(null=chla_null,
                    all=chla_all,
                    tmp=chla_tmp,
                    amm=chla_amm,
                    phe=chla_phe)

lms_best.stats <- mapply(glance, lms_best)
colnames(lms_best.stats) <- names(lms_best)
lms_best.stats
```

    ##               null      all          tmp         amm          phe         
    ## r.squared     0         0.5929971    0.2595211   0.1718766    0.3558155   
    ## adj.r.squared 0         0.5555716    0.2516436   0.1630668    0.3489624   
    ## sigma         2.304473  1.536287     1.993544    2.108225     1.859407    
    ## statistic     NA        15.84471     32.94487    19.50966     51.92092    
    ## p.value       NA        3.569059e-14 1.15385e-07 2.675906e-05 1.426426e-10
    ## df            1         9            2           2            2           
    ## logLik        -215.8612 -172.7124    -201.4393   -206.8088    -194.7523   
    ## AIC           435.7225  365.4247     408.8785    419.6176     395.5046    
    ## BIC           440.8512  391.0682     416.5716    427.3106     403.1976    
    ## deviance      504.5064  205.3356     373.5763    417.7935     324.9952    
    ## df.residual   95        87           94          94           94

From the comparison, I would say that **the most important variable explaining Chl-a is Pheophytin a** which is able to explain 35% of the variability in Chl-a compared to 25% with Temperature and only 16% with Ammonia.

``` r
# plot residuals

par(mfrow=c(1,2))
plot(chla_phe, pch=16, which=1)
plot(chla_phe, pch=16, which=2)
```

![](hw7_pertiwi_multivariate-analyses_files/figure-markdown_github/unnamed-chunk-23-1.png)

``` r
par(mfrow=c(1,1))
```

The are somewhat normal and distributed along the 0 line although improvements may need to be done, such as investigation potential outliers (57, 45, and 3).

``` r
# plot predicted vs actual

plot(predict(chla_phe),wq_model$`Chlorophyll a`,
     xlab="predicted",ylab="actual", main = "Chl-a ~ Pheophytin a")
abline(a=0,b=1)
```

![](hw7_pertiwi_multivariate-analyses_files/figure-markdown_github/unnamed-chunk-24-1.png)

``` r
# plot Chl-a vs Pheophytin a actual values with prediction line from model

ggplot(wq_model, aes(x = `Pheophytin a`, y = `Chlorophyll a`)) +
  geom_point() +
  geom_line(aes(y = predict(chla_phe)), shape = 1)
```

![](hw7_pertiwi_multivariate-analyses_files/figure-markdown_github/unnamed-chunk-25-1.png)

The predicted values from the model are somewhat predicting the trend of the actual values.

------------------------------------------------------------------------

### <span style="color:blue">2| Parallel Regression</span>

------------------------------------------------------------------------

#### **--Data preparation**

``` r
# add seasons

wq_meanmonthly_season <- wq_meanmonthly

wq_meanmonthly_season$Season <- ifelse(wq_meanmonthly_season$Month > 9 |
                                   wq_meanmonthly_season$Month < 5, "wet season", "dry season")

wq_season <- wq_meanmonthly_season %>%
  ungroup() %>%
  select(Season, everything()) %>%
  select(Month, everything()) %>%
  select(Year, everything())

wq_season$Season <- as.factor(wq_season$Season)
```

``` r
# create new dataset for parallel regression model

wq_prl_model <- select(wq_season,
                       `Chlorophyll a`,
                       Season,
                       `Pheophytin a`)
```

#### **--Model comparison**

``` r
# perform lms

# Chl-a ~ season
chla_season <- lm(`Chlorophyll a` ~ Season, data = wq_prl_model)

# Chl-a ~ Pheophytin a
chla_phea <- lm(`Chlorophyll a` ~ `Pheophytin a`, data = wq_prl_model)

# Chl-a ~ season + pheophytin a
chla_season.phea <- lm(`Chlorophyll a` ~ Season + `Pheophytin a`, data = wq_prl_model)
```

``` r
# compare models
lms_prl <- list(season=chla_season, phea=chla_phea, season.phea=chla_season.phea)

lms_prl.stats <- mapply(glance, lms_prl)
colnames(lms_prl.stats) <- names(lms_prl)
lms_prl.stats
```

    ##               season      phea         season.phea 
    ## r.squared     0.2360985   0.3558155    0.4413055   
    ## adj.r.squared 0.2279719   0.3489624    0.4292905   
    ## sigma         2.024828    1.859407     1.740921    
    ## statistic     29.05251    51.92092     36.72974    
    ## p.value       5.21002e-07 1.426426e-10 1.752311e-12
    ## df            2           2            3           
    ## logLik        -202.9341   -194.7523    -187.9179   
    ## AIC           411.8681    395.5046     383.8359    
    ## BIC           419.5612    403.1976     394.0933    
    ## deviance      385.3932    324.9952     281.865     
    ## df.residual   94          94           93

The addition of the season category does seem to improve the model. The parallel model is able to explain 43% of the variation while the model with only Pheophytin a as the variable is only able to explain 35% of the variation.

#### **--Residuals**

``` r
# plot residuals

par(mfrow=c(2,3))

plot(chla_season, pch=16, which=1)
plot(chla_phea, pch=16, which=1)
plot(chla_season.phea, pch=16, which=1)

plot(chla_season, pch=16, which=2)
plot(chla_phea, pch=16, which=2)
plot(chla_season.phea, pch=16, which=2)
```

![](hw7_pertiwi_multivariate-analyses_files/figure-markdown_github/unnamed-chunk-30-1.png)

``` r
par(mfrow=c(1,1))
```

Left: Chl-a ~ season; Middle: Chl-a ~ Pheophytin a; Right: Parallel model.

#### **--Plot**

``` r
# extract intercept and slope values to generate regression lines

season_coef <- data.frame(t(coef(chla_season)))
colnames(season_coef) <- c("intercept", "slope")
season_coef
```

    ##   intercept     slope
    ## 1  5.020391 -2.259392

``` r
phea_coef <- data.frame(t(coef(chla_phea)))
colnames(phea_coef) <- c("intercept", "slope")
phea_coef
```

    ##   intercept    slope
    ## 1 0.4458041 1.908262

``` r
season.phea_coef <- data.frame(t(coef(chla_season.phea)))
colnames(season.phea_coef) <- c("intercept", "slope.season", "slope.phe")
season.phea_coef
```

    ##   intercept slope.season slope.phe
    ## 1  1.904096    -1.455273  1.551185

``` r
# plot Chl-a vs Pheophytin a with regression lines

plot(wq_prl_model$`Pheophytin a`, wq_prl_model$`Chlorophyll a`, 
     xlab = "Pheophytin a",
     ylab = "Chlorophyll a",
     col = wq_prl_model$Season,
     pch = 16,
     main = "By Season (dry season: black, dry season: red)")
abline(reg = chla_phea, lty = 2)
abline(a = season.phea_coef$intercept, b = season.phea_coef$slope.phe, col = "blue")
```

![](hw7_pertiwi_multivariate-analyses_files/figure-markdown_github/unnamed-chunk-32-1.png) The blue line is the regression line from the parallel model while the dashed black line is the regression line for the `Chl-a ~ Pheophytin a` model. Although not entirely clear from the plot, it seems as though the parallel model is accounting for for the dry season points compared to the univariate model. This may be because of the influence of the season variable in the parallel model. The plot also shows how the dry season may have a higher mean value than the wet season as it's values seem to be higher than values from the wet season.

------------------------------------------------------------------------

### <span style="color:blue">3| GitHub Link</span>

------------------------------------------------------------------------

<https://github.com/cnpw/es207_hw7>
