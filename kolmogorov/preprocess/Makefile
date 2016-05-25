F90 = gfortran
INC = /usr/include/

OBJ = readbin.o cpl2f90.o
flags =  -malign-double -fall-intrinsics -ffree-line-length-none -I$(INC) -cpp  -g -fcheck=all -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow 
#flags =  -malign-double -fall-intrinsics -ffree-line-length-none  -I$(INC) -cpp -Ofast

cpl2f90: $(OBJ)
	$(F90) $(flags) -o $@ $(OBJ)
%.o : %.f90
	$(F90) $(flags) -c $<
clean: 
	rm *.mod *.o
