#!/USR/BIN/ENV sh
rm main.exe
ifort -fpp -standard-semantics -O3 -Wl,-rpath,../../../lib -I../../../inc main.F90 ../../../lib/libparamonte* -o main.exe
./main.exe