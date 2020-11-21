CC = g++
OPT = -O3 -DNDEBUG
PTH = ../include
PTHEIGEN = /usr/local/include/eigen3
PTHKDTR = /usr/local/include
PTHCGAL = /usr/local/include
PTHBIN = ../bin
PTHLIB = /usr/local/lib
LIB = -L$(PTHLIB)
INC = -I$(PTH) -I$(PTHEIGEN) -I$(PTHKDTR)
CFLAGS = -Wall -std=c++17 $(OPT)
LFLAGS = -lgmp

all : BeadPack BeadPack_Reactive BeadPack_Fp BeadPack_Return BeadPack_Strips BeadPack_ReactionMap

BeadPack : BeadPack.o
	$(CC) $(CFLAGS) $(LIB) $(LFLAGS) -o BeadPack BeadPack.o
	rm BeadPack.o
	mv BeadPack $(PTHBIN)/BeadPack

BeadPack_Reactive : BeadPack_Reactive.o
	$(CC) $(CFLAGS) $(LIB) $(LFLAGS) -o BeadPack_Reactive BeadPack_Reactive.o
	rm BeadPack_Reactive.o
	mv BeadPack_Reactive $(PTHBIN)/BeadPack_Reactive

BeadPack_Fp : BeadPack_Fp.o
	$(CC) $(CFLAGS) $(LIB) $(LFLAGS) -o BeadPack_Fp BeadPack_Fp.o
	rm BeadPack_Fp.o
	mv BeadPack_Fp $(PTHBIN)/BeadPack_Fp

BeadPack_Return : BeadPack_Return.o
	$(CC) $(CFLAGS) $(LIB) $(LFLAGS) -o BeadPack_Return BeadPack_Return.o
	rm BeadPack_Return.o
	mv BeadPack_Return $(PTHBIN)/BeadPack_Return

BeadPack_Strips : BeadPack_Strips.o
	$(CC) $(CFLAGS) $(LIB) $(LFLAGS) -o BeadPack_Strips BeadPack_Strips.o
	rm BeadPack_Strips.o
	mv BeadPack_Strips $(PTHBIN)/BeadPack_Strips

BeadPack_ReactionMap : BeadPack_ReactionMap.o
	$(CC) $(CFLAGS) $(LIB) $(LFLAGS) -o BeadPack_ReactionMap BeadPack_ReactionMap.o
	rm BeadPack_ReactionMap.o
	mv BeadPack_ReactionMap $(PTHBIN)/BeadPack_ReactionMap

BeadPack.o : BeadPack.cpp
	$(CC) $(CFLAGS) $(INC) -c $<

BeadPack_Reactive.o : BeadPack_Reactive.cpp
	$(CC) $(CFLAGS) $(INC) -c $<

BeadPack_Fp.o : BeadPack_Fp.cpp
	$(CC) $(CFLAGS) $(INC) -c $<

BeadPack_Return.o : BeadPack_Return.cpp
	$(CC) $(CFLAGS) $(INC) -c $<

BeadPack_Strips.o : BeadPack_Strips.cpp
	$(CC) $(CFLAGS) $(INC) -c $<

BeadPack_ReactionMap.o : BeadPack_ReactionMap.cpp
	$(CC) $(CFLAGS) $(INC) -c $<

clean :
	rm -f BeadPack BeadPack_Reactive BeadPack_Fp BeadPack_Strips BeadPack_Return BeadPack_ReactionMap BeadPack.o BeadPack_Reactive.o BeadPack_Fp.o BeadPack_Strips.o BeadPack_Return.o BeadPack_ReactionMap.o
