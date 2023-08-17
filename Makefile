
1d_rk1_fou_serial:
	g++ -O3 -std=c++17 1d_rk1_fou_serial.cc -lfmt -o 1d_rk1_fou_serial

clean:
	rm 1d_rk1_fou_serial
	rm *.h5 tpH5.txt log.txt *.dat
	
