rem PysCeS: Metatool compile script

g++ meta4.3_double_gcc4.3.2.cpp -Wno-deprecated -Wno-write-strings -O0 -o meta43_double.exe
g++ meta4.3_int_gcc4.3.2.cpp -Wno-deprecated -Wno-write-strings -fpermissive -O0 -o meta43_int.exe

exit
