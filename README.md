# SmallPT - Language Benchmark
Based on the [SmallPT](https://www.kevinbeason.com/smallpt/) global illumination renderer developed by Kevin Beason.

## C++
```
g++ -O3 -fopenmp smallpt.cpp -o smallpt && time ./smallpt 5000
```

## .Net Core
```
dotnet run -c release 5000
```

## Java 11
```
java SmallPT.java 5000
```

## Rust
```
cargo run --release 5000
```

## JavaScript
```
node smallpt.js 5000
```












--------------------------------
C++
--------------------------------
Rendering (5000 spp) 100.00%

Total: 10:16.36
User: 9471.74
Kernel: 7.84

CPU: 1537%
MEM: 47136K

--------------------------------
Rust 1.48.0
--------------------------------
Rendering (5000 spp) 73.40%
thread '<unknown>' has overflowed its stack
fatal runtime error: stack overflow
Command terminated by signal 6
Total: 7:44.99
User: 7372.63
Kernel: 9.17

CPU: 1587%
MEM: 40704K

--------------------------------
.NET Core 5.0.101
--------------------------------
Rendering (5000 spp) 100.00%

Total: 15:12.10
User: 14142.30
Kernel: 3.64

CPU: 1550%
MEM: 109776K

--------------------------------
Java 15.0.1
--------------------------------
Rendering (5000 spp) 100.00%

Total: 11:33.88
User: 10929.85
Kernel: 12.79

CPU: 1577%
MEM: 1684836K

--------------------------------
JavaScript v15.4.0
--------------------------------
Rendering (5000 spp) 100.00%

Total: 53:39.25
User: 50181.13
Kernel: 19.82

CPU: 1559%
MEM: 3199180K
