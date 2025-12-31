# Benchmarks

Timings based on `make run` vs. `make run_legacy`.

System     | Original   | Optimized | Threads | Speedup
---------- | ---------- | --------- | ------- | -------
Apple M4   | 110        |  28       | 1       | 3.9
Apple M4   | 110        |  17       | 2       | 6.4
Apple M4   | 110        |  12       | 4       | 9.2
AMD Ryzen  | 134        |  36       | 1       | 3.7
AMD Ryzen  | 134        |  20       | 2       | 6.7
AMD Ryzen  | 134        |  12       | 4       | 12
AMD Ryzen  | 134        |  7.1      | 8       | 19
AMD Ryzen  | 134        |  3.9      | 16      | 35
AMD Ryzen  | 134        |  3.0      | 24      | 45
AMD Ryzen  | 134        |  3.2      | 32      | 42
