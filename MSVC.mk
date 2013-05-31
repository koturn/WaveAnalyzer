### This Makefile was written for nmake. ###
CC      = cl
CFLAGS  = /O2 /Wall /c
LDFLAGS = /O2
RM      = del /F
TARGET  = waveAnalyzer.exe
SRCDIR  = src
OBJ1    = $(SRCDIR)\$(TARGET:.exe=.obj)
OBJ2    = $(SRCDIR)\fourier.obj
OBJ3    = $(SRCDIR)\gplot.obj
OBJ4    = $(SRCDIR)\wav.obj
SRC1    = $(OBJ1:.obj=.c)
SRC2    = $(OBJ2:.obj=.c)
SRC3    = $(OBJ3:.obj=.c)
SRC4    = $(OBJ4:.obj=.c)
HEADER  = $(SRCDIR)\compatibility.h


.SUFFIXES : .exe .obj
.exe.obj :
	$(CC) $(LDFLAGS) $**

.SUFFIXES : .in .out
.c.obj :
	$(CC) $(CFLAGS) $** /Fo$@



all : $(TARGET)

$(TARGET) : $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4)
	$(CC) $(LDFLAGS) $**

$(OBJ1) : $(SRC1) $(OBJ2:.obj=.h) $(OBJ3:.obj=.h) $(OBJ4:.obj=.h) $(HEADER)

$(OBJ2) : $(SRC2) $(OBJ2:.obj=.h) $(HEADER)

$(OBJ3) : $(SRC3) $(OBJ3:.obj=.h) $(HEADER)

$(OBJ4) : $(SRC4) $(OBJ4:.obj=.h) $(HEADER)


clean :
	$(RM) $(TARGET) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4)
objclean :
	$(RM) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4)
