CC = g++
OPT = -O3 -DNDEBUG
PTH = ../include
PTHINC = /usr/local/include
PTHBIN = ../bin
PTHLIB = /usr/local/lib

LIB = -L$(PTHLIB)
INC = -I$(PTH) -I$(PTHINC)
CFLAGS = -Wall -std=c++17 $(OPT)
LFLAGS = -lgmp -lmpfr

all : BeadPack_Conservative BeadPack_Reactive BeadPack_Fp BeadPack_Return BeadPack_Strips BeadPack_ReactionMap BeadPack_Statistics

BeadPack_Conservative : BeadPack_Conservative.o
	$(CC) $(CFLAGS) $(LIB) -o BeadPack_Conservative BeadPack_Conservative.o $(LFLAGS)
	rm BeadPack_Conservative.o
	mv BeadPack_Conservative $(PTHBIN)/BeadPack_Conservative

BeadPack_Reactive : BeadPack_Reactive.o
	$(CC) $(CFLAGS) $(LIB) -o BeadPack_Reactive BeadPack_Reactive.o $(LFLAGS)
	rm BeadPack_Reactive.o
	mv BeadPack_Reactive $(PTHBIN)/BeadPack_Reactive

BeadPack_Fp : BeadPack_Fp.o
	$(CC) $(CFLAGS) $(LIB) -o BeadPack_Fp BeadPack_Fp.o $(LFLAGS)
	rm BeadPack_Fp.o
	mv BeadPack_Fp $(PTHBIN)/BeadPack_Fp

BeadPack_Return : BeadPack_Return.o
	$(CC) $(CFLAGS) $(LIB) -o BeadPack_Return BeadPack_Return.o $(LFLAGS)
	rm BeadPack_Return.o
	mv BeadPack_Return $(PTHBIN)/BeadPack_Return

BeadPack_Strips : BeadPack_Strips.o
	$(CC) $(CFLAGS) $(LIB) -o BeadPack_Strips BeadPack_Strips.o $(LFLAGS)
	rm BeadPack_Strips.o
	mv BeadPack_Strips $(PTHBIN)/BeadPack_Strips

BeadPack_ReactionMap : BeadPack_ReactionMap.o
	$(CC) $(CFLAGS) $(LIB) -o BeadPack_ReactionMap BeadPack_ReactionMap.o $(LFLAGS)
	rm BeadPack_ReactionMap.o
	mv BeadPack_ReactionMap $(PTHBIN)/BeadPack_ReactionMap

BeadPack_Statistics : BeadPack_Statistics.o
	$(CC) $(CFLAGS) $(LIB) -o BeadPack_Statistics BeadPack_Statistics.o $(LFLAGS)
	rm BeadPack_Statistics.o
	mv BeadPack_Statistics $(PTHBIN)/BeadPack_Statistics

BeadPack_Conservative.o : BeadPack_Conservative.cpp
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

BeadPack_Statistics.o : BeadPack_Statistics.cpp
	$(CC) $(CFLAGS) $(INC) -c $<

clean :
	rm -f BeadPack_Conservative BeadPack_Reactive BeadPack_Fp BeadPack_Strips BeadPack_Return BeadPack_ReactionMap BeadPack_Statistics BeadPack_Conservative.o BeadPack_Reactive.o BeadPack_Fp.o BeadPack_Strips.o BeadPack_Return.o BeadPack_ReactionMap.o BeadPack_Statistics.o
