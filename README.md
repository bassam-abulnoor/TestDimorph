# TestDimorph

<!-- badges: start -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/TestDimorph)](https://cran.r-project.org/package=TestDimorph)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/TestDimorph)](https://cran.r-project.org/package=TestDimorph)
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
## Acknowledgement

The authors would like to express their endless gratitude to **Prof. Konigsberg** for his contribution to this package with authentic code used in the `multivariate`, `extract_sum`,`RawGen` functions, data sets used in the examples and for his overall valuable work in the field of physical anthropology without which this package wouldn't have been possible.

## Data entry

The input can be in the form of wide format **raw data** or **summary statistics** where each row represents measurements of a single individual or population. In order to be recognized by the functions, columns need to have a specific name (case sensitive):


1. **Sex** : Either **M** for male or **F** for female (*Raw data*)
2. **Pop** : Populations of the measured individuals (*Summary/Raw data*)
3. **Parms** : Measured numerical parameters for multivariate analysis and raw data generation (*Summary*)
4. **m** : Male sample size (*Summary*)
5. **f** :Female sample size (*Summary*)
6. **M.mu**: Male mean (*Summary*)
7. **F.mu** :Female mean (*Summary*)
8. **M.sdev** : Male standard deviation (*Summary*)
9. **F.sdev** :Female standard deviation (*Summary*)

**N.B**

* All non numerical values need to be entered as factors.
* The list and the data frame formats for the `multivariate` and `RawGen` functions can be be found in the attached data sets under the names `baboon.parms_list` and `baboon.parms_df` respectively.
* In the input of most of the functions, columns are referred to by number to avoid confusion.

## Examples

#### Summary statistics for univariate analysis

##### Comparison of femur head diameters in two populations


|    Pop    |  m  | M.mu  | M.sdev |  f  | F.mu  | F.sdev |
| :-------: | :-: | :---: | :----: | :-: | :---: | :----: |
|  Turkish  | 150 | 49.39 |  3.01  | 150 | 42.91 |  2.90  |
| Bulgarian | 82  | 48.33 |  2.53  | 58  | 42.89 |  2.84  |



#### Summary statistics for multivariate analysis and raw data generation

##### Comparison of femur head diameters in two populations


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



