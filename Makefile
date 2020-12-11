CC = gcc
CXX = g++
CFLAGS = -I/Users/daniele/local/include -I/usr/include -Wall -g -O3 -DDEBUG
LDFLAGS = -L/Users/daniele/local/lib
LIBS = -lsundials_cvode -lsundials_nvecserial -lm -ldl -lpthread
OBJS = dynasty.o

%.o : %.c
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY : linuxlib

sharedlib : $(OBJS)
	$(CXX) -g -o dynasty.so -fPIC -shared \
		-static-libgcc -static-libstdc++ \
		$(OBJS) \
		-Wl,-zmuldefs \
		-Wl,--whole-archiv \
		-Wl,--no-whole-archiv \
		$(shell python3-config --ldflags) \
		-Wl,--library-path=/home/daniele/local/sundials/lib \
		$(LIBS)

.PHONY : maclib

maclib : $(OBJS)
	$(CXX) -g -o dynasty.dylib -fPIC -shared \
		$(OBJS) \
		$(LDFLAGS) \
		$(LIBS)

clean :
	rm -f dynasty *.o *.dylib *.so
