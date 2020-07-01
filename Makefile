packmod.o: packmod.f90
	gfortran -Wunused -O3 -c $^

studentpack.o: studentpack.f90 packmod.o
	gfortran -Wunused -O3 -c $<

studentpack: studentpack.o packmod.o
	gfortran -Wunused -O3 -L$(ALGENCAN_LIB_PATH) $^ -lalgencan -o $@

clean:
	rm -f packmod.o packmod.mod studentpack studentpack.o
