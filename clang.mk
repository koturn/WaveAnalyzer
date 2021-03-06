### This Makefile was written for GNU Make. ###
CC      = clang
CFLAGS  = -pipe -O3 -Wall -Wextra -Wno-unknown-pragmas $(if $(STD), $(addprefix -std=, $(STD)),)
LDFLAGS = -pipe -O3 -s
TARGET  = waveAnalyzer
SRCDIR  = src
OBJ1    = $(SRCDIR)/main.o
OBJ2    = $(SRCDIR)/fourier.o
OBJ3    = $(SRCDIR)/gplot.o
OBJ4    = $(SRCDIR)/wav.o
HEADER  = $(SRCDIR)/compatibility.h


ifeq ($(OS),Windows_NT)
    TARGET := $(addsuffix .exe, $(TARGET))
endif

%.exe :
	$(CC) $(LDFLAGS) $(filter %.c %.o, $^) $(LDLIBS) -o $@


.PHONY : all
all : $(TARGET)

$(TARGET) : $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(HEADER)

$(OBJ1) : $(OBJ1:%.o=%.c) $(OBJ2:%.o=%.h) $(OBJ3:%.o=%.h) $(OBJ4:%.o=%.h) $(HEADER)

$(OBJ2) : $(OBJ2:%.o=%.c) $(OBJ2:%.o=%.h) $(HEADER)

$(OBJ3) : $(OBJ3:%.o=%.c) $(OBJ3:%.o=%.h) $(HEADER)

$(OBJ4) : $(OBJ4:%.o=%.c) $(OBJ4:%.o=%.h) $(HEADER)


.PHONY : clean
clean :
	$(RM) $(TARGET) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4)
.PHONY : objclean
objclean :
	$(RM) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4)
