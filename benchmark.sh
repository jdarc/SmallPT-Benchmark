#! /bin/sh
spp=5000
format="Total: %E\nUser: %U\nKernel: %S\n\nCPU: %P\nMEM: %MK"

echo "--------------------------------"
echo "C++"
echo "--------------------------------"
g++ -O3 -fopenmp smallpt.cpp -o smallpt
/usr/bin/time -f "$format" ./smallpt $spp
echo ""

echo "--------------------------------"
echo "Rust" $( cargo --version 2>&1 | head -n 1 | awk -F ' ' '{print $2}')
echo "--------------------------------"
cargo build --release --quiet
/usr/bin/time -f "$format" cargo run --release --quiet $spp
echo ""

echo "--------------------------------"
echo ".NET Core" $(dotnet --version)
echo "--------------------------------"
dotnet build -c release > /dev/null
/usr/bin/time -f "$format" dotnet run -c release $spp
echo ""

echo "--------------------------------"
echo "Java" $(java -version 2>&1 | head -n 1)
echo "--------------------------------"
javac SmallPT.java
/usr/bin/time -f "$format" java SmallPT $spp
echo ""

echo "--------------------------------"
echo $(kotlin -version)
echo "--------------------------------"
kotlinc -no-reflect -jvm-target 11 SmallPT.kt -d SmallPT.jar 2> /dev/null
/usr/bin/time -f "$format" kotlin -cp SmallPT.jar SmallPT $spp
echo ""

echo "--------------------------------"
echo "JavaScript" $(node --version)
echo "--------------------------------"
/usr/bin/time -f "$format" node smallpt.js $spp

rm -rf target bin obj
rm *.class *.jar
rm smallpt
