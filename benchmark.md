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
| instances |           POP_build |           pb_build_SDPMosekstruct |           pb_construction |           pb_mosek_solve |
| --------- | -------------------:| ---------------------------------:| -------------------------:| ------------------------:|
|       WB2 |           45.130 ms |                          5.915 ms |                 33.937 ms |                32.317 ms |
|     LMBM3 |          147.365 ms |                         11.696 ms |                 80.194 ms |                48.073 ms |


order_2
| instances |           POP_build |           pb_build_SDPMosekstruct |           pb_construction |           pb_mosek_solve |
| --------- | -------------------:| ---------------------------------:| -------------------------:| ------------------------:|
|       WB2 |           42.331 ms |                         48.805 ms |                341.357 ms |                88.134 ms |
|     LMBM3 |          126.692 ms |                        221.792 ms |                   1.813 s |               194.405 ms |

### Specific benchmarks:

#### Initial benchmark results
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

```
MomentRelaxation(), case89pegase, d=1
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

#### Removing Sorted Structures:

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
MomentRelaxation(), case89pegase, d=1
  memory estimate:  545.63 MiB
  allocs estimate:  11182870
  --------------
  minimum time:     2.110 s (36.22% GC)
  median time:      2.175 s (38.66% GC)
  mean time:        2.309 s (42.06% GC)
  maximum time:     2.640 s (49.54% GC)
  --------------
  samples:          3
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