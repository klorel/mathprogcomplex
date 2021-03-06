# Benchmark SOShierarchy

Results of the hierarchy benchmark script `dev/benchmark_hierarchy`.

Timing will be done on the REX server. Julia system image is precompiled, julia version is 0.6.3.
Cpuinfo is as follows :

```
processor       : 0
vendor_id       : GenuineIntel
cpu family      : 6
model           : 42
model name      : Intel Xeon E312xx (Sandy Bridge)
stepping        : 1
microcode       : 0x1
cpu MHz         : 2397.222
cache size      : 4096 KB
fpu             : yes
fpu_exception   : yes
cpuid level     : 13
wp              : yes
flags           : fpu de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx fxsr sse sse2 syscall nx rdtscp lm constant_tsc arch_perfmon rep_good nopl pni pclmulqdq ssse3 cx16 pcid sse4_1 sse4_2 x2apic popcnt tsc_deadline_timer aes xsave avx hypervi
sor lahf_lm xsaveopt
bogomips        : 4794.44
clflush size    : 64
cache_alignment : 64
address sizes   : 40 bits physical, 48 bits virtual
```

## 06 june 2018
Benchmarking results :

order_1

| instances |           1.POP_build |           2.pb_construction |           3.pb_build_SDPMosekstruct |           4.pb_mosek_solve |
| --------- | ---------------------:| ---------------------------:| -----------------------------------:| --------------------------:|
|       WB2 |             56.285 ms |                   28.250 ms |                            7.372 ms |                  29.601 ms |
|     LMBM3 |            133.133 ms |                   94.086 ms |                           16.390 ms |                  45.613 ms |
|       WB5 |            264.598 ms |                  147.907 ms |                           23.777 ms |                  73.157 ms |
|   case6ww |            485.896 ms |                  325.370 ms |                           35.288 ms |                  64.492 ms |
|     case9 |            224.405 ms |                  419.656 ms |                           45.620 ms |                 111.820 ms |

order_2

| instances |           1.POP_build |           2.pb_construction |           3.pb_build_SDPMosekstruct |           4.pb_mosek_solve |
| --------- | ---------------------:| ---------------------------:| -----------------------------------:| --------------------------:|
|       WB2 |             52.534 ms |                  407.634 ms |                           55.889 ms |                 131.921 ms |
|     LMBM3 |            139.650 ms |                     2.571 s |                          256.090 ms |                 468.952 ms |
|       WB5 |            227.937 ms |                    13.783 s |                          954.469 ms |                    3.274 s |
|   case6ww |            485.768 ms |                    36.592 s |                             2.053 s |                    7.373 s |
|     case9 |            267.784 ms |                    89.101 s |                             4.347 s |                  141.749 s |

### Specific benchmarks: MomentRelaxation()

#### Before any improvements

MomentRelaxation(), case30, d=1
```
  memory estimate:  374.77 MiB
  allocs estimate:  8117070
  --------------
  minimum time:     1.077 s (17.42% GC)
  median time:      1.187 s (22.11% GC)
  mean time:        1.213 s (24.76% GC)
  maximum time:     1.346 s (31.57% GC)
  --------------
  samples:          5
  evals/sample:     1
```

MomentRelaxation(), case89pegase, d=1
```
BenchmarkTools.Trial:
  memory estimate:  3.61 GiB
  allocs estimate:  80377019
  --------------
  minimum time:     14.887 s (22.51% GC)
  median time:      14.887 s (22.51% GC)
  mean time:        14.887 s (22.51% GC)
  maximum time:     14.887 s (22.51% GC)
  --------------
  samples:          1
  evals/sample:     1
```

#### Removing Sorted Structures

```
MomentRelaxation(), case30, d=1
  memory estimate:  83.77 MiB
  allocs estimate:  1707556
  --------------
  minimum time:     254.627 ms (15.71% GC)
  median time:      292.712 ms (24.97% GC)
  mean time:        309.915 ms (29.00% GC)
  maximum time:     416.560 ms (48.47% GC)
  --------------
  samples:          17
  evals/sample:     1
```

```
  memory estimate:  509.83 MiB
  allocs estimate:  10403116
  --------------
  minimum time:     2.032 s (37.80% GC)
  median time:      2.078 s (39.65% GC)
  mean time:        2.233 s (43.12% GC)
  maximum time:     2.590 s (50.08% GC)
  --------------
  samples:          3
  evals/sample:     1
```

### Specific benchmarks: build_SDPInstance()

#### Before improvements

build_SDPInstance(), case30, d=1
```
  memory estimate:  123.22 MiB
  allocs estimate:  2846439
  --------------
  minimum time:     397.124 ms (10.93% GC)
  median time:      416.186 ms (11.33% GC)
  mean time:        417.275 ms (11.44% GC)
  maximum time:     444.241 ms (13.60% GC)
  --------------
  samples:          13
  evals/sample:     1
```

build_SDPInstance(), case89pegase, d=1
```
  memory estimate:  1.06 GiB
  allocs estimate:  25103390
  --------------
  minimum time:     3.824 s (15.38% GC)
  median time:      3.865 s (15.81% GC)
  mean time:        3.865 s (15.81% GC)
  maximum time:     3.905 s (16.23% GC)
  --------------
  samples:          2
  evals/sample:     1
```

#### Removing Sorted Structures:

build_SDPInstance(), case30, d=1
```
  memory estimate:  2.13 MiB
  allocs estimate:  63165
  --------------
  minimum time:     11.314 ms (0.00% GC)
  median time:      13.707 ms (0.00% GC)
  mean time:        15.624 ms (4.61% GC)
  maximum time:     40.133 ms (0.00% GC)
  --------------
  samples:          320
  evals/sample:     1
```

build_SDPInstance(), case89pegase, d=1
```
  memory estimate:  15.04 MiB
  allocs estimate:  448344
  --------------
  minimum time:     116.857 ms (0.00% GC)
  median time:      135.347 ms (11.83% GC)
  mean time:        133.295 ms (8.86% GC)
  maximum time:     152.965 ms (9.19% GC)
  --------------
  samples:          38
  evals/sample:     1
```

### End of day timings
order_1

| instances |           1.POP_build |           2.pb_construction |           3.pb_build_SDPMosekstruct |           4.pb_mosek_solve |
| --------- | ---------------------:| ---------------------------:| -----------------------------------:| --------------------------:|
|       WB2 |             51.372 ms |                   16.967 ms |                            8.793 ms |                  31.642 ms |
|     LMBM3 |            157.561 ms |                   41.722 ms |                           16.289 ms |                  44.674 ms |
|       WB5 |            281.963 ms |                   60.528 ms |                           41.304 ms |                 111.494 ms |
|   case6ww |            635.079 ms |                  115.655 ms |                          225.262 ms |                 189.157 ms |
|     case9 |            327.359 ms |                   74.992 ms |                           88.304 ms |                 289.448 ms |

order_2

| instances |           1.POP_build |           2.pb_construction |           3.pb_build_SDPMosekstruct |           4.pb_mosek_solve |
| --------- | ---------------------:| ---------------------------:| -----------------------------------:| --------------------------:|
|       WB2 |             51.632 ms |                  108.792 ms |                           60.252 ms |                 150.116 ms |
|     LMBM3 |            141.495 ms |                  737.637 ms |                          280.620 ms |                 455.741 ms |
|       WB5 |            246.962 ms |                     2.804 s |                             1.001 s |                    3.781 s |
|   case6ww |            543.881 ms |                     6.506 s |                             2.245 s |                    9.069 s |
|     case9 |            253.338 ms |                    15.039 s |                             5.176 s |                  148.176 s |

## 7 june 2018

changes : addindex! function.

order_1

| instances |           1.POP_build |           2.pb_construction |           3.pb_build_SDPMosekstruct |           4.pb_mosek_solve |
| --------- | ---------------------:| ---------------------------:| -----------------------------------:| --------------------------:|
|       WB2 |             46.554 ms |                   14.502 ms |                            6.057 ms |                  28.743 ms |
|     LMBM3 |            132.754 ms |                   38.195 ms |                           12.925 ms |                  56.226 ms |
|       WB5 |            213.880 ms |                   55.178 ms |                           21.882 ms |                  69.439 ms |
|   case6ww |            466.185 ms |                  103.533 ms |                           28.754 ms |                  68.627 ms |
|     case9 |            238.251 ms |                   65.787 ms |                           47.890 ms |                 118.387 ms |

order_2

| instances |           1.POP_build |           2.pb_construction |           3.pb_build_SDPMosekstruct |           4.pb_mosek_solve |
| --------- | ---------------------:| ---------------------------:| -----------------------------------:| --------------------------:|
|       WB2 |             47.870 ms |                   85.025 ms |                           46.952 ms |                 123.663 ms |
|     LMBM3 |            127.208 ms |                  580.757 ms |                          212.954 ms |                 404.964 ms |
|       WB5 |            207.997 ms |                     1.924 s |                          814.236 ms |                    3.243 s |
|   case6ww |            513.525 ms |                     5.076 s |                             1.898 s |                    6.827 s |
|     case9 |            208.660 ms |                    11.547 s |                             4.141 s |                  146.763 s |
