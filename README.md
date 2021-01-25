# TestDimorph

<!-- badges: start -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/TestDimorph)](https://cran.r-project.org/package=TestDimorph)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/TestDimorph)](https://cran.r-project.org/package=TestDimorph)
[![Codecov test coverage](https://codecov.io/gh/bassam-abulnoor/TestDimorph/branch/master/graph/badge.svg)](https://codecov.io/gh/bassam-abulnoor/TestDimorph?branch=master)
<!-- badges: end -->

Analysis Of The Interpopulation Difference In Degree of Sexual Dimorphism Using Summary Statistics

## Installation

You can install the released version of TestDimorph from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("TestDimorph")
```
or the development version from GitHub:

``` r
devtools::install_github("bassam-abulnoor/TestDimorph")
```

## Data entry

The input can be in the form of wide format **raw data** or **summary statistics** where each row represents measurements of a single individual or population. In order to be recognized by the functions, columns need to have specific names (case sensitive):


1. **Sex** : either **M** for male or **F** for female (*Raw data*)
2. **Pop** : studied populations (*Summary/Raw data*)
3. **Trait** : studied traits or parameters (*Summary*)
4. **m** : male sample size (*Summary*)
5. **f** : female sample size (*Summary*)
6. **M.mu**: male mean (*Summary*)
7. **F.mu** : female mean (*Summary*)
8. **M.sdev** : male standard deviation (*Summary*)
9. **F.sdev** : female standard deviation (*Summary*)

**N.B**

* All non numerical values need to be factors.
* The list and the data frame formats for the `multivariate` and `RawGen` functions input, can be be found in the attached data sets under the names `baboon.parms_list` and `baboon.parms_df` respectively.

## Examples of data entry

#### Summary statistics for univariate analysis

##### Comparison of maximal femur length and mediolateral diameter in mid-shaft


|    Pop    |  m  | M.mu  | M.sdev |  f  | F.mu  | F.sdev |
| :-------: | :-: | :---: | :----: | :-: | :---: | :----: |
|  Turkish  | 150 | 49.39 |  3.01  | 150 | 42.91 |  2.90  |
| Bulgarian | 82  | 48.33 |  2.53  | 58  | 42.89 |  2.84  |



#### Summary statistics for multivariate analysis and raw data generation

##### Comparison of maximal femur length and mediolateral diameter in mid-shaft


|    Pop    | Parms |  m  |  M.mu  | M.sdev |  f  |  F.mu  | F.sdev |
| :-------: | :---: | :-: | :----: | :----: | :-: | :----: | :----: |
| Bulgarian | MXFL  | 82  | 461.80 | 19.90  | 58  | 411.70 | 23.20  |
|   Greek   | MXFL  | 36  | 440.40 | 19.60  | 34  | 409.80 | 21.40  |
| Bulgarian |  MLD  | 82  | 27.67  |  2.21  | 58  | 24.89  |  1.78  |
|   Greek   |  MLD  | 36  | 27.74  |  1.79  | 34  | 26.69  |  2.42  |



#### Raw data for summary extraction

##### The Howells' craniometric data

| Sex |  Pop  | GOL | NOL | BNL | BBH | XCB | XFB | ZYB | AUB |
| :-: | :---: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
|  M  | NORSE | 189 | 185 | 100 | 135 | 143 | 120 | 133 | 119 |
|  F  | NORSE | 179 | 177 | 97  | 127 | 140 | 114 | 122 | 118 |
|  M  | EGYPT | 192 | 190 | 104 | 132 | 140 | 116 | 134 | 125 |
|  F  | EGYPT | 174 | 172 | 93  | 123 | 140 | 118 | 121 | 114 |


#### Summary statistics for van_vark function

##### Howells summary data

| Trait |  Pop  |   M.mu   |  m |   F.mu   |  f |
|:-----:|:-----:|:--------:|:--:|:--------:|:--:|
|  GOL  | EGYPT | 185.6207 | 58 | 175.5849 | 53 |
|  GOL  | NORSE | 188.4727 | 55 | 179.9818 | 55 |
|  NOL  | EGYPT |  183.569 | 58 | 175.0377 | 53 |
|  NOL  | NORSE | 186.1818 | 55 | 178.6909 | 55 |
