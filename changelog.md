# Changelog

## 06 june 2018
Benchmarking results :

order_1
| instances |           POP_build |           pb_build_SDPMosekstruct |           pb_construction |           pb_mosek_solve |
| --------- | -------------------:| ---------------------------------:| -------------------------:| ------------------------:|
|       WB2 |           45.130 ms |                          5.915 ms |                 33.937 ms |                32.317 ms |
|     LMBM3 |          147.365 ms |                         11.696 ms |                 80.194 ms |                48.073 ms |


order_2
| instances |           POP_build |           pb_build_SDPMosekstruct |           pb_construction |           pb_mosek_solve |
| --------- | -------------------:| ---------------------------------:| -------------------------:| ------------------------:|
|       WB2 |           42.331 ms |                         48.805 ms |                341.357 ms |                88.134 ms |
|     LMBM3 |          126.692 ms |                        221.792 ms |                   1.813 s |               194.405 ms |
