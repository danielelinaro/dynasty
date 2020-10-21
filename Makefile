CC = gcc
CXX = g++
CFLAGS = -g -Wall -fPIC -I/home/daniele/local/sundials/include $(shell python3-config --cflags)
OBJS = dynasty.o

%.o : %.c
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY : sharedlib

sharedlib : $(OBJS)
	$(CXX) -g -o dynasty.so -fPIC -shared \
		-static-libgcc -static-libstdc++ \
		$(OBJS) \
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
	rm -f *.o *.so
