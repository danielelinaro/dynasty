CC = gcc
CXX = g++
CFLAGS = -Wall -g -I/home/daniele/local/sundials/include
LDFLAGS = -L/home/daniele/local/sundials/lib
LIBS = -lsundials_cvode -lsundials_nvecserial -lm
OBJS = dynasty.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

dynasty : $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)
	
sharedlib :
	$(CC) -g -c -fPIC -DPYTHON $(shell python3-config --cflags) \
		$(CFLAGS) dynasty.c -o dynasty.o
	$(CXX) -g -o dynasty.so -fPIC -shared \
		-static-libgcc -static-libstdc++ \
		dynasty.o \
		-Wl,-zmuldefs \
		-Wl,--whole-archiv \
		-Wl,--no-whole-archiv \
		$(shell python3-config --ldflags) \
		-Wl,--library-path=/home/daniele/local/sundials/lib \
		-lsundials_cvode \
		-lsundials_nvecserial \
		-lm \
		-ldl \
		-lpthread

clean :
	rm -f dynasty *.o *.so
