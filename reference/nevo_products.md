# Nevo (2000) Cereal Product Data

Product-level data from the Nevo (2000) study of the US ready-to-eat
cereal market. Contains 2256 product-market observations across 94
markets (47 cities x 2 quarters).

## Usage

``` r
nevo_products
```

## Format

A data frame with 2256 rows and the following columns:

- market_ids:

  Market identifier (city-quarter combination)

- city_ids:

  City identifier

- quarter:

  Quarter (1 or 2)

- product_ids:

  Product identifier within market

- firm_ids:

  Firm identifier

- brand_ids:

  Brand identifier

- shares:

  Market share

- prices:

  Product price (dollars per serving)

- sugar:

  Sugar content (grams per serving)

- mushy:

  Mushiness indicator (0/1)

- demand_instruments0:

  Excluded demand-side instrument 0

- demand_instruments1:

  Excluded demand-side instrument 1

- demand_instruments2:

  Excluded demand-side instrument 2

- demand_instruments3:

  Excluded demand-side instrument 3

- demand_instruments4:

  Excluded demand-side instrument 4

- demand_instruments5:

  Excluded demand-side instrument 5

- demand_instruments6:

  Excluded demand-side instrument 6

- demand_instruments7:

  Excluded demand-side instrument 7

- demand_instruments8:

  Excluded demand-side instrument 8

- demand_instruments9:

  Excluded demand-side instrument 9

- demand_instruments10:

  Excluded demand-side instrument 10

- demand_instruments11:

  Excluded demand-side instrument 11

- demand_instruments12:

  Excluded demand-side instrument 12

- demand_instruments13:

  Excluded demand-side instrument 13

- demand_instruments14:

  Excluded demand-side instrument 14

- demand_instruments15:

  Excluded demand-side instrument 15

- demand_instruments16:

  Excluded demand-side instrument 16

- demand_instruments17:

  Excluded demand-side instrument 17

- demand_instruments18:

  Excluded demand-side instrument 18

- demand_instruments19:

  Excluded demand-side instrument 19

## Source

<https://pyblp.readthedocs.io/en/stable/>

## References

Nevo, A. (2000). A Practitioner's Guide to Estimation of
Random-Coefficients Logit Models of Demand. *Journal of Economics &
Management Strategy*, 9(4), 513-548.

## Examples

``` r
nevo_products <- load_nevo_products()
head(nevo_products)
#>   market_ids city_ids quarter product_ids firm_ids brand_ids      shares
#> 1      C01Q1        1       1       F1B04        1         4 0.012417212
#> 2      C01Q1        1       1       F1B06        1         6 0.007809387
#> 3      C01Q1        1       1       F1B07        1         7 0.012994511
#> 4      C01Q1        1       1       F1B09        1         9 0.005769961
#> 5      C01Q1        1       1       F1B11        1        11 0.017934141
#> 6      C01Q1        1       1       F1B13        1        13 0.026601892
#>       prices sugar mushy demand_instruments0 demand_instruments1
#> 1 0.07208794     2     1          -0.2159728          0.04057341
#> 2 0.11417849    18     1          -0.2452393          0.05474226
#> 3 0.13239066     4     1          -0.1764587          0.04659597
#> 4 0.13034408     3     0          -0.1214013          0.04876037
#> 5 0.15482331    12     0          -0.1326114          0.03962835
#> 6 0.13704921    14     0          -0.1534998          0.04298842
#>   demand_instruments2 demand_instruments3 demand_instruments4
#> 1           -3.247948        -0.523937690         -0.23246005
#> 2          -19.832461        -0.180519690          0.01468859
#> 3           -2.878531        -0.284219000         -0.21553691
#> 4           -2.059918        -0.328412260         -0.22206995
#> 5           -6.137598        -0.138625100         -0.18936521
#> 6           -8.417332         0.007829087         -0.13850121
#>   demand_instruments5 demand_instruments6 demand_instruments7
#> 1        0.0068326605           3.1397395         -0.57478633
#> 2        0.0007988026           0.2876539          0.03293960
#> 3       -0.0318693280           2.8862741         -0.74976495
#> 4       -0.0314740400           4.4531096          0.25567529
#> 5       -0.0437471020          -3.5546508          0.13882114
#> 6       -0.0210582270          -2.7594799          0.05020052
#>   demand_instruments8 demand_instruments9 demand_instruments10
#> 1           0.2062201           0.1774656            2.1163580
#> 2           0.1051208          -0.2875618           -7.3740909
#> 3          -0.4789565           0.2147389            2.1878721
#> 4          -0.4729673           0.3560980            2.7045762
#> 5          -0.6886784           0.2602726            1.2612419
#> 6          -0.2734440           0.1273060            0.3375543
#>   demand_instruments11 demand_instruments12 demand_instruments13
#> 1          -0.15470824        -0.0057964065           0.01453801
#> 2          -0.57641176         0.0129908540           0.07614324
#> 3          -0.20734643         0.0035092777           0.09178117
#> 4           0.04074801        -0.0037242656           0.09473168
#> 5           0.03483558        -0.0005676374           0.10245147
#> 6           0.02351037         0.0002637777           0.08627983
#>   demand_instruments14 demand_instruments15 demand_instruments16
#> 1           0.12624398           0.06734464           0.06842261
#> 2           0.02973565           0.08786672           0.11050060
#> 3           0.16377308           0.11188073           0.10822551
#> 4           0.13527378           0.08809001           0.10176745
#> 5           0.13063951           0.08481820           0.10107461
#> 6           0.07233581           0.02225051           0.10564387
#>   demand_instruments17 demand_instruments18 demand_instruments19
#> 1           0.03480046           0.12634612           0.03548368
#> 2           0.08778380           0.04987192           0.07257905
#> 3           0.08643905           0.12234707           0.10184248
#> 4           0.10177748           0.11074119           0.10433204
#> 5           0.12516923           0.13346381           0.12111110
#> 6           0.11603699           0.09965064           0.10572660
```
