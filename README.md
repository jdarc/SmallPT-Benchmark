# Path Tracing Benchmark
Based on [SmallPT](https://www.kevinbeason.com/smallpt/) global illumination renderer developed by Kevin Beason.

To run the benchmarks:
```
./benchmark.sh
```
It will execute each program at 5000 samples per pixel.

## Results
Razer Blade 15 Advanced Edition - Super Max-Q - Turbo Boost

CPU | MEM | OS
--- | --- | ---
i7-10875H | 16GB | Ubuntu 20.10 |

### C++
Total | User | Kernel | CPU | MEM
--- | --- | --- | --- | --- 
10:07.28 | 9619.36 | 5.58 | 1584% | 24,520K

### Rust 1.49.0
Total | User | Kernel | CPU | MEM
--- | --- | --- | --- | --- 
10:31.51 | 9997.22 | 11.05 | 1584% | 38,152K

### .NET Core 5.0.101
Total | User | Kernel | CPU | MEM
--- | --- | --- | --- | ---
15:44.71 | 14026.58 | 0.79 | 1484% | 110,040K

### Java 15.0.1
Total | User | Kernel | CPU | MEM
--- | --- | --- | --- | ---
12:35.02 | 11178.18 | 13.48 | 1482% | 1,675,952K

### Kotlin version 1.4.21-release-351 (JRE 15.0.1+9)
Total | User | Kernel | CPU | MEM
--- | --- | --- | --- | ---
14:06.16 | 12715.47 | 71.27 | 1511% | 342,232K

### JavaScript - Node v15.4.0
Total | User | Kernel | CPU | MEM
--- | --- | --- | --- | ---
53:39.25 | 50181.13 | 19.82 | 1559% | 3,199,180K

### WebGL2
Total | User | Kernel | CPU | MEM
--- | --- | --- | --- | ---
00:02.05 ||||

## License
MIT License

Copyright © 2021 Jean d'Arc

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
