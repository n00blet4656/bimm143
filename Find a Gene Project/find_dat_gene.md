Find a Gene
================

``` r
library(bio3d)
```

``` r
muscle <- read.fasta("muscle alignment.fasta")
```

``` r
seqmatrix <- seqidentity(muscle)
seqmatrix
```

    ##                Human Cavia Urocitellus Ictidomys M.flaviventris Marmota
    ## Human          1.000 0.283       0.271     0.271          0.282   0.271
    ## Cavia          0.283 1.000       0.815     0.815          0.820   0.820
    ## Urocitellus    0.271 0.815       1.000     0.974          0.959   0.969
    ## Ictidomys      0.271 0.815       0.974     1.000          0.959   0.969
    ## M.flaviventris 0.282 0.820       0.959     0.959          1.000   0.984
    ## Marmota        0.271 0.820       0.969     0.969          0.984   1.000

``` r
heatmap(seqmatrix, cexRow = 0.7, cexCol = 0.7)
```

![](find_dat_gene_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
blast.pdb("5F6Z_A")
```

    ##  Searching ... please wait (updates every 5 seconds) RID = 8G2K5MR4014 
    ##  .
    ##  Reporting 66 hits

    ## $hit.tbl
    ##    queryid subjectids identity alignmentlength mismatches gapopens q.start
    ## 1   5F6Z_A     5F6Z_A  100.000             170          0        0       1
    ## 2   5F6Z_A     5EZ2_A   99.408             169          1        0       2
    ## 3   5F6Z_A     2HZQ_A   40.764             157         91        2       6
    ## 4   5F6Z_A     1GKA_B   26.220             164        119        2       1
    ## 5   5F6Z_A     1Z24_A   29.114             158        101        6       5
    ## 6   5F6Z_A     1H91_A   27.059             170        115        5       1
    ## 7   5F6Z_A     1S2P_A   27.059             170        115        5       1
    ## 8   5F6Z_A     1I4U_A   27.219             169        114        5       2
    ## 9   5F6Z_A     4ALO_A   26.627             169        115        5       2
    ## 10  5F6Z_A     2ACO_A   30.909             165        101        8       2
    ## 11  5F6Z_A     3MBT_A   30.380             158         99        6      13
    ## 12  5F6Z_A     1QWD_A   31.579             152         93        7      13
    ## 13  5F6Z_A     1IIU_A   31.757             148         84        6       7
    ## 14  5F6Z_A     1JYJ_A   30.714             140         84        6      13
    ## 15  5F6Z_A     2WR6_A   30.714             140         84        6      13
    ## 16  5F6Z_A     2WQA_E   30.714             140         84        6      13
    ## 17  5F6Z_A     1QAB_E   30.714             140         84        6      13
    ## 18  5F6Z_A     3BSZ_E   30.714             140         84        6      13
    ## 19  5F6Z_A     2WQ9_A   30.714             140         84        6      13
    ## 20  5F6Z_A     1JYD_A   30.714             140         84        6      13
    ## 21  5F6Z_A     1BRP_A   30.714             140         84        6      13
    ## 22  5F6Z_A     3FMZ_A   30.714             140         84        6      13
    ## 23  5F6Z_A     4O9S_A   30.137             146         89        6       7
    ## 24  5F6Z_A     1RLB_E   26.241             141         92        5      13
    ## 25  5F6Z_A     1AQB_A   28.467             137         86        5       7
    ## 26  5F6Z_A     1KT5_A   26.241             141         92        5      13
    ## 27  5F6Z_A     1HBQ_A   26.241             141         92        5      13
    ## 28  5F6Z_A     1KT3_A   26.241             141         92        5      13
    ## 29  5F6Z_A     1BBP_A   27.338             139         89        4       6
    ## 30  5F6Z_A     1T0V_A   27.660             141         86        7       6
    ## 31  5F6Z_A     3EBW_A   22.930             157        104        7      17
    ## 32  5F6Z_A     1N0S_A   27.660             141         86        7       6
    ## 33  5F6Z_A     1KXO_A   27.193             114         73        3       6
    ## 34  5F6Z_A     3EBK_A   42.857              35         20        0      92
    ## 35  5F6Z_A     4N7D_A   42.857              35         20        0      92
    ## 36  5F6Z_A     5WWL_M   45.161              31         17        0     140
    ## 37  5F6Z_A     2IVD_A   34.375              64         35        2     111
    ## 38  5F6Z_A     1EPA_A   55.556              18          8        0      12
    ## 39  5F6Z_A     2WWP_A   31.884              69         39        2      96
    ## 40  5F6Z_A     3O19_A   31.884              69         39        2      96
    ## 41  5F6Z_A     4OS0_A   23.529             170         85        6      12
    ## 42  5F6Z_A     1OCP_A   38.889              36         20        1       7
    ## 43  5F6Z_A     4ORR_A   31.884              69         39        2      96
    ## 44  5F6Z_A     4ORW_A   31.884              69         39        2      96
    ## 45  5F6Z_A     2B0C_A   28.049              82         36        3      16
    ## 46  5F6Z_A     3VE7_A   36.111              36         23        0     126
    ## 47  5F6Z_A     4OS8_A   33.898              59         38        1      99
    ## 48  5F6Z_A     4OS3_A   33.898              59         38        1      99
    ## 49  5F6Z_A     2XST_A   19.632             163        108        4      13
    ## 50  5F6Z_A     5N8Y_M   55.000              20          9        0      70
    ## 51  5F6Z_A     5C5E_A   55.000              20          9        0      70
    ## 52  5F6Z_A     4G86_A   55.000              20          9        0      70
    ## 53  5F6Z_A     2A2G_A   21.477             149        101        5      16
    ## 54  5F6Z_A     1R8J_A   55.000              20          9        0      70
    ## 55  5F6Z_A     3ZQ3_A   21.192             151         99        6      16
    ## 56  5F6Z_A     3KZA_A   24.528             106         68        3      11
    ## 57  5F6Z_A     1M2E_A   55.000              20          9        0      70
    ## 58  5F6Z_A     3I3L_A   25.000             104         58        4       1
    ## 59  5F6Z_A     2L9C_A   19.672              61         49        0      96
    ## 60  5F6Z_A     3JYV_G   30.769              39         27        0     112
    ## 61  5F6Z_A     1S1H_G   30.769              39         27        0     112
    ## 62  5F6Z_A     4UER_G   30.769              39         27        0     112
    ## 63  5F6Z_A     4CUY_F   30.769              39         27        0     112
    ## 64  5F6Z_A     5IT9_F   30.769              39         27        0     112
    ## 65  5F6Z_A     3IZB_F   30.769              39         27        0     112
    ## 66  5F6Z_A     4UJF_G   30.769              39         27        0     112
    ##    q.end s.start s.end    evalue bitscore positives  mlog.evalue pdb.id
    ## 1    170       1   170 8.88e-128    356.0    100.00 292.54709035 5F6Z_A
    ## 2    170       2   170 1.00e-126    353.0     99.41 290.12572172 5EZ2_A
    ## 3    160       4   160  2.27e-38    129.0     58.60  86.67845370 2HZQ_A
    ## 4    162       5   168  3.11e-12     62.4     39.02  26.49639839 1GKA_B
    ## 5    153       6   161  4.07e-12     62.0     48.10  26.22737812 1Z24_A
    ## 6    163       4   171  6.28e-12     61.6     44.12  25.79365114 1H91_A
    ## 7    163       5   172  6.79e-12     61.6     44.12  25.71557017 1S2P_A
    ## 8    163       6   172  9.07e-12     61.2     43.79  25.42604885 1I4U_A
    ## 9    163       6   172  1.48e-10     57.8     43.20  22.63380884 4ALO_A
    ## 10   160      13   170  1.20e-09     55.1     44.24  20.54094428 2ACO_A
    ## 11   166      12   162  1.71e-09     54.7     44.30  20.18677247 3MBT_A
    ## 12   160      30   174  1.85e-09     54.7     45.39  20.10808020 1QWD_A
    ## 13   139       4   149  4.87e-09     53.5     43.24  19.14017190 1IIU_A
    ## 14   139      12   151  1.07e-08     52.8     45.00  18.35302210 1JYJ_A
    ## 15   139      12   151  3.80e-08     51.2     45.00  17.08567968 2WR6_A
    ## 16   139      12   151  3.86e-08     51.2     45.00  17.07001356 2WQA_E
    ## 17   139       8   147  4.05e-08     51.2     45.00  17.02196386 1QAB_E
    ## 18   139      11   150  4.05e-08     51.2     45.00  17.02196386 3BSZ_E
    ## 19   139      11   150  4.11e-08     51.2     45.00  17.00725772 2WQ9_A
    ## 20   139      12   151  4.28e-08     51.2     45.00  16.96672773 1JYD_A
    ## 21   139      11   150  4.77e-08     50.8     45.00  16.85833444 1BRP_A
    ## 22   139      40   179  5.27e-08     51.2     45.00  16.75865038 3FMZ_A
    ## 23   139      37   182  6.28e-08     51.2     43.84  16.58331076 4O9S_A
    ## 24   141      11   151  1.04e-06     47.4     45.39  13.77628984 1RLB_E
    ## 25   131       5   141  1.05e-06     47.4     43.07  13.76672039 1AQB_A
    ## 26   141      11   151  1.18e-06     47.0     45.39  13.64999612 1KT5_A
    ## 27   141      11   151  1.40e-06     47.0     45.39  13.47903832 1HBQ_A
    ## 28   141      11   151  1.43e-06     47.0     45.39  13.45783611 1KT3_A
    ## 29   136       6   140  7.94e-06     44.7     41.01  11.74359728 1BBP_A
    ## 30   136       6   140  5.96e-04     39.7     46.10   7.42526989 1T0V_A
    ## 31   165      13   160  6.58e-04     39.3     42.68   7.32630563 3EBW_A
    ## 32   136       6   140  6.70e-04     39.3     46.10   7.30823285 1N0S_A
    ## 33   113       6   115  1.80e-02     35.0     40.35   4.01738352 1KXO_A
    ## 34   126      91   125  1.20e-01     32.7     54.29   2.12026354 3EBK_A
    ## 35   126      91   125  1.20e-01     32.7     54.29   2.12026354 4N7D_A
    ## 36   170      85   115  4.80e-01     31.2     64.52   0.73396918 5WWL_M
    ## 37   170     199   259  9.50e-01     30.8     46.88   0.05129329 2IVD_A
    ## 38    29       1    18  1.10e+00     30.0     83.33  -0.09531018 1EPA_A
    ## 39   156      89   157  1.10e+00     30.0     46.38  -0.09531018 2WWP_A
    ## 40   156      82   150  1.10e+00     29.6     46.38  -0.09531018 3O19_A
    ## 41   156      29   178  1.30e+00     29.6     37.65  -0.26236426 4OS0_A
    ## 42    40      29    64  1.40e+00     28.1     63.89  -0.33647224 1OCP_A
    ## 43   156     110   178  2.10e+00     29.3     46.38  -0.74193734 4ORR_A
    ## 44   156     110   178  2.10e+00     29.3     46.38  -0.74193734 4ORW_A
    ## 45    97      20    78  3.20e+00     28.9     36.59  -1.16315081 2B0C_A
    ## 46   161       5    40  3.80e+00     28.5     58.33  -1.33500107 3VE7_A
    ## 47   156     120   178  4.20e+00     28.1     49.15  -1.43508453 4OS8_A
    ## 48   156     120   178  4.20e+00     28.1     49.15  -1.43508453 4OS3_A
    ## 49   163       8   159  5.50e+00     27.7     39.26  -1.70474809 2XST_A
    ## 50    89      81   100  5.60e+00     28.1     70.00  -1.72276660 5N8Y_M
    ## 51    89      81   100  5.70e+00     28.1     70.00  -1.74046617 5C5E_A
    ## 52    89      81   100  5.70e+00     28.1     70.00  -1.74046617 4G86_A
    ## 53   156      28   168  5.80e+00     27.7     36.91  -1.75785792 2A2G_A
    ## 54    89      86   105  5.80e+00     28.1     70.00  -1.75785792 1R8J_A
    ## 55   156      22   162  6.80e+00     27.7     39.07  -1.91692261 3ZQ3_A
    ## 56   109       4   104  6.90e+00     27.3     33.96  -1.93152141 3KZA_A
    ## 57    89      81   100  7.30e+00     27.3     70.00  -1.98787435 1M2E_A
    ## 58    92     434   529  7.30e+00     28.1     41.35  -1.98787435 3I3L_A
    ## 59   156     103   163  7.40e+00     27.3     40.98  -2.00148000 2L9C_A
    ## 60   150     134   172  8.00e+00     27.3     58.97  -2.07944154 3JYV_G
    ## 61   150      98   136  9.20e+00     26.9     58.97  -2.21920348 1S1H_G
    ## 62   150     154   192  9.50e+00     27.3     58.97  -2.25129180 4UER_G
    ## 63   150     154   192  9.60e+00     27.3     58.97  -2.26176310 4CUY_F
    ## 64   150     154   192  9.70e+00     27.3     58.97  -2.27212589 5IT9_F
    ## 65   150     173   211  9.80e+00     27.3     58.97  -2.28238239 3IZB_F
    ## 66   150     172   210  9.90e+00     27.3     58.97  -2.29253476 4UJF_G
    ##       acc
    ## 1  5F6Z_A
    ## 2  5EZ2_A
    ## 3  2HZQ_A
    ## 4  1GKA_B
    ## 5  1Z24_A
    ## 6  1H91_A
    ## 7  1S2P_A
    ## 8  1I4U_A
    ## 9  4ALO_A
    ## 10 2ACO_A
    ## 11 3MBT_A
    ## 12 1QWD_A
    ## 13 1IIU_A
    ## 14 1JYJ_A
    ## 15 2WR6_A
    ## 16 2WQA_E
    ## 17 1QAB_E
    ## 18 3BSZ_E
    ## 19 2WQ9_A
    ## 20 1JYD_A
    ## 21 1BRP_A
    ## 22 3FMZ_A
    ## 23 4O9S_A
    ## 24 1RLB_E
    ## 25 1AQB_A
    ## 26 1KT5_A
    ## 27 1HBQ_A
    ## 28 1KT3_A
    ## 29 1BBP_A
    ## 30 1T0V_A
    ## 31 3EBW_A
    ## 32 1N0S_A
    ## 33 1KXO_A
    ## 34 3EBK_A
    ## 35 4N7D_A
    ## 36 5WWL_M
    ## 37 2IVD_A
    ## 38 1EPA_A
    ## 39 2WWP_A
    ## 40 3O19_A
    ## 41 4OS0_A
    ## 42 1OCP_A
    ## 43 4ORR_A
    ## 44 4ORW_A
    ## 45 2B0C_A
    ## 46 3VE7_A
    ## 47 4OS8_A
    ## 48 4OS3_A
    ## 49 2XST_A
    ## 50 5N8Y_M
    ## 51 5C5E_A
    ## 52 4G86_A
    ## 53 2A2G_A
    ## 54 1R8J_A
    ## 55 3ZQ3_A
    ## 56 3KZA_A
    ## 57 1M2E_A
    ## 58 3I3L_A
    ## 59 2L9C_A
    ## 60 3JYV_G
    ## 61 1S1H_G
    ## 62 4UER_G
    ## 63 4CUY_F
    ## 64 5IT9_F
    ## 65 3IZB_F
    ## 66 4UJF_G
    ## 
    ## $raw
    ##    queryid subjectids identity alignmentlength mismatches gapopens q.start
    ## 1   5F6Z_A     5F6Z_A  100.000             170          0        0       1
    ## 2   5F6Z_A     5EZ2_A   99.408             169          1        0       2
    ## 3   5F6Z_A     2HZQ_A   40.764             157         91        2       6
    ## 4   5F6Z_A     1GKA_B   26.220             164        119        2       1
    ## 5   5F6Z_A     1Z24_A   29.114             158        101        6       5
    ## 6   5F6Z_A     1H91_A   27.059             170        115        5       1
    ## 7   5F6Z_A     1S2P_A   27.059             170        115        5       1
    ## 8   5F6Z_A     1I4U_A   27.219             169        114        5       2
    ## 9   5F6Z_A     4ALO_A   26.627             169        115        5       2
    ## 10  5F6Z_A     2ACO_A   30.909             165        101        8       2
    ## 11  5F6Z_A     3MBT_A   30.380             158         99        6      13
    ## 12  5F6Z_A     1QWD_A   31.579             152         93        7      13
    ## 13  5F6Z_A     1IIU_A   31.757             148         84        6       7
    ## 14  5F6Z_A     1JYJ_A   30.714             140         84        6      13
    ## 15  5F6Z_A     2WR6_A   30.714             140         84        6      13
    ## 16  5F6Z_A     2WQA_E   30.714             140         84        6      13
    ## 17  5F6Z_A     1QAB_E   30.714             140         84        6      13
    ## 18  5F6Z_A     3BSZ_E   30.714             140         84        6      13
    ## 19  5F6Z_A     2WQ9_A   30.714             140         84        6      13
    ## 20  5F6Z_A     1JYD_A   30.714             140         84        6      13
    ## 21  5F6Z_A     1BRP_A   30.714             140         84        6      13
    ## 22  5F6Z_A     3FMZ_A   30.714             140         84        6      13
    ## 23  5F6Z_A     4O9S_A   30.137             146         89        6       7
    ## 24  5F6Z_A     1RLB_E   26.241             141         92        5      13
    ## 25  5F6Z_A     1AQB_A   28.467             137         86        5       7
    ## 26  5F6Z_A     1KT5_A   26.241             141         92        5      13
    ## 27  5F6Z_A     1HBQ_A   26.241             141         92        5      13
    ## 28  5F6Z_A     1KT3_A   26.241             141         92        5      13
    ## 29  5F6Z_A     1BBP_A   27.338             139         89        4       6
    ## 30  5F6Z_A     1T0V_A   27.660             141         86        7       6
    ## 31  5F6Z_A     3EBW_A   22.930             157        104        7      17
    ## 32  5F6Z_A     1N0S_A   27.660             141         86        7       6
    ## 33  5F6Z_A     1KXO_A   27.193             114         73        3       6
    ## 34  5F6Z_A     3EBK_A   42.857              35         20        0      92
    ## 35  5F6Z_A     4N7D_A   42.857              35         20        0      92
    ## 36  5F6Z_A     5WWL_M   45.161              31         17        0     140
    ## 37  5F6Z_A     2IVD_A   34.375              64         35        2     111
    ## 38  5F6Z_A     1EPA_A   55.556              18          8        0      12
    ## 39  5F6Z_A     2WWP_A   31.884              69         39        2      96
    ## 40  5F6Z_A     3O19_A   31.884              69         39        2      96
    ## 41  5F6Z_A     4OS0_A   23.529             170         85        6      12
    ## 42  5F6Z_A     1OCP_A   38.889              36         20        1       7
    ## 43  5F6Z_A     4ORR_A   31.884              69         39        2      96
    ## 44  5F6Z_A     4ORW_A   31.884              69         39        2      96
    ## 45  5F6Z_A     2B0C_A   28.049              82         36        3      16
    ## 46  5F6Z_A     3VE7_A   36.111              36         23        0     126
    ## 47  5F6Z_A     4OS8_A   33.898              59         38        1      99
    ## 48  5F6Z_A     4OS3_A   33.898              59         38        1      99
    ## 49  5F6Z_A     2XST_A   19.632             163        108        4      13
    ## 50  5F6Z_A     5N8Y_M   55.000              20          9        0      70
    ## 51  5F6Z_A     5C5E_A   55.000              20          9        0      70
    ## 52  5F6Z_A     4G86_A   55.000              20          9        0      70
    ## 53  5F6Z_A     2A2G_A   21.477             149        101        5      16
    ## 54  5F6Z_A     1R8J_A   55.000              20          9        0      70
    ## 55  5F6Z_A     3ZQ3_A   21.192             151         99        6      16
    ## 56  5F6Z_A     3KZA_A   24.528             106         68        3      11
    ## 57  5F6Z_A     1M2E_A   55.000              20          9        0      70
    ## 58  5F6Z_A     3I3L_A   25.000             104         58        4       1
    ## 59  5F6Z_A     2L9C_A   19.672              61         49        0      96
    ## 60  5F6Z_A     3JYV_G   30.769              39         27        0     112
    ## 61  5F6Z_A     1S1H_G   30.769              39         27        0     112
    ## 62  5F6Z_A     4UER_G   30.769              39         27        0     112
    ## 63  5F6Z_A     4CUY_F   30.769              39         27        0     112
    ## 64  5F6Z_A     5IT9_F   30.769              39         27        0     112
    ## 65  5F6Z_A     3IZB_F   30.769              39         27        0     112
    ## 66  5F6Z_A     4UJF_G   30.769              39         27        0     112
    ##    q.end s.start s.end    evalue bitscore positives
    ## 1    170       1   170 8.88e-128    356.0    100.00
    ## 2    170       2   170 1.00e-126    353.0     99.41
    ## 3    160       4   160  2.27e-38    129.0     58.60
    ## 4    162       5   168  3.11e-12     62.4     39.02
    ## 5    153       6   161  4.07e-12     62.0     48.10
    ## 6    163       4   171  6.28e-12     61.6     44.12
    ## 7    163       5   172  6.79e-12     61.6     44.12
    ## 8    163       6   172  9.07e-12     61.2     43.79
    ## 9    163       6   172  1.48e-10     57.8     43.20
    ## 10   160      13   170  1.20e-09     55.1     44.24
    ## 11   166      12   162  1.71e-09     54.7     44.30
    ## 12   160      30   174  1.85e-09     54.7     45.39
    ## 13   139       4   149  4.87e-09     53.5     43.24
    ## 14   139      12   151  1.07e-08     52.8     45.00
    ## 15   139      12   151  3.80e-08     51.2     45.00
    ## 16   139      12   151  3.86e-08     51.2     45.00
    ## 17   139       8   147  4.05e-08     51.2     45.00
    ## 18   139      11   150  4.05e-08     51.2     45.00
    ## 19   139      11   150  4.11e-08     51.2     45.00
    ## 20   139      12   151  4.28e-08     51.2     45.00
    ## 21   139      11   150  4.77e-08     50.8     45.00
    ## 22   139      40   179  5.27e-08     51.2     45.00
    ## 23   139      37   182  6.28e-08     51.2     43.84
    ## 24   141      11   151  1.04e-06     47.4     45.39
    ## 25   131       5   141  1.05e-06     47.4     43.07
    ## 26   141      11   151  1.18e-06     47.0     45.39
    ## 27   141      11   151  1.40e-06     47.0     45.39
    ## 28   141      11   151  1.43e-06     47.0     45.39
    ## 29   136       6   140  7.94e-06     44.7     41.01
    ## 30   136       6   140  5.96e-04     39.7     46.10
    ## 31   165      13   160  6.58e-04     39.3     42.68
    ## 32   136       6   140  6.70e-04     39.3     46.10
    ## 33   113       6   115  1.80e-02     35.0     40.35
    ## 34   126      91   125  1.20e-01     32.7     54.29
    ## 35   126      91   125  1.20e-01     32.7     54.29
    ## 36   170      85   115  4.80e-01     31.2     64.52
    ## 37   170     199   259  9.50e-01     30.8     46.88
    ## 38    29       1    18  1.10e+00     30.0     83.33
    ## 39   156      89   157  1.10e+00     30.0     46.38
    ## 40   156      82   150  1.10e+00     29.6     46.38
    ## 41   156      29   178  1.30e+00     29.6     37.65
    ## 42    40      29    64  1.40e+00     28.1     63.89
    ## 43   156     110   178  2.10e+00     29.3     46.38
    ## 44   156     110   178  2.10e+00     29.3     46.38
    ## 45    97      20    78  3.20e+00     28.9     36.59
    ## 46   161       5    40  3.80e+00     28.5     58.33
    ## 47   156     120   178  4.20e+00     28.1     49.15
    ## 48   156     120   178  4.20e+00     28.1     49.15
    ## 49   163       8   159  5.50e+00     27.7     39.26
    ## 50    89      81   100  5.60e+00     28.1     70.00
    ## 51    89      81   100  5.70e+00     28.1     70.00
    ## 52    89      81   100  5.70e+00     28.1     70.00
    ## 53   156      28   168  5.80e+00     27.7     36.91
    ## 54    89      86   105  5.80e+00     28.1     70.00
    ## 55   156      22   162  6.80e+00     27.7     39.07
    ## 56   109       4   104  6.90e+00     27.3     33.96
    ## 57    89      81   100  7.30e+00     27.3     70.00
    ## 58    92     434   529  7.30e+00     28.1     41.35
    ## 59   156     103   163  7.40e+00     27.3     40.98
    ## 60   150     134   172  8.00e+00     27.3     58.97
    ## 61   150      98   136  9.20e+00     26.9     58.97
    ## 62   150     154   192  9.50e+00     27.3     58.97
    ## 63   150     154   192  9.60e+00     27.3     58.97
    ## 64   150     154   192  9.70e+00     27.3     58.97
    ## 65   150     173   211  9.80e+00     27.3     58.97
    ## 66   150     172   210  9.90e+00     27.3     58.97
    ## 
    ## $url
    ##                                                                                                                                                        8G2K5MR4014 
    ## "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=Alignment&ALIGNMENT_VIEW=Tabular&RESULTS_FILE=on&FORMAT_TYPE=CSV&ALIGNMENTS=20000&RID=8G2K5MR4014" 
    ## 
    ## attr(,"class")
    ## [1] "blast"
