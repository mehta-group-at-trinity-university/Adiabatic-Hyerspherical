ThreeBody.sech.x:	ThreeBody.sech.f Bsplines.f ThreeBody.sech.o  
	g95   -O4 ThreeBody.sech.o -L/sw/lib/ -framework vecLib  -I/opt/local/lib/ -L/opt/local/lib/ -I/usr/local/lib/ -larpack -lgfortran -lg2c -o ThreeBody.sech.x

ThreeBody.sech.o:	ThreeBody.sech.f
	g95   -ftrace=full -ffixed-line-length-132 -O4 -c ThreeBody.sech.f
