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
