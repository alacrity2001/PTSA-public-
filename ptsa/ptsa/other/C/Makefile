CC=gcc
CFLAGS=

OBJ = main.o VALMOD.o FiboKLAMD.o

# $@ says to put the output of the compilation in the file named on the left side of the:
#the $< is the first item in the dependencies list

%.o: %.c
	$(CC) $(CFLAGS) -c -fPIC -o $@ $<


MATRIXPROFILE: $(OBJ)
	$(CC) $(CFLAGS) -fPIC -shared -o $@.so $^ -lfftw3 -lm
	rm $(OBJ)


clean:
	rm MATRIXPROFILE.*
