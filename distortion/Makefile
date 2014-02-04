SRCS=distortion.cpp
OBJS=$(SRCS:.cpp=.o) m5op_arm.o

include ../compile.mk

all: libar.a

libar.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $?

clean:
	rm -f *.o *.a
