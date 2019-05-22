ThreeBody.sech.x:	ThreeBody.sech.f Bsplines.f ThreeBody.sech.o  
	gfortran   -O4 ThreeBody.sech.o  -framework accelerate  -I/usr/local/lib/ -larpack -o ThreeBody.sech.x

ThreeBody.sech.o:	ThreeBody.sech.f
	gfortran   -ffixed-line-length-132 -O4 -c ThreeBody.sech.f
