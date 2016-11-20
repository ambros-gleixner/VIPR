COMP = g++
CFLAGS = -O1 -g -Wall -std=c++11 -Wshadow -W -Wmissing-declarations -Wpointer-arith -Wundef -Wcast-qual -Wcast-align -Wwrite-strings -DNDEBUG

all: viprchk vipr2html viprttn

clean:
	rm *.o

viprchk: viprchk.o 
	$(COMP) $(CFLAGS) -o viprchk viprchk.o -lgmpxx -lgmp

viprchk.o:  viprchk.cpp
	$(COMP) $(CFLAGS) -c viprchk.cpp


vipr2html: vipr2html.o 
	$(COMP) $(CFLAGS) -o vipr2html vipr2html.o

vipr2html.o:  vipr2html.cpp
	$(COMP) $(CFLAGS) -c vipr2html.cpp


viprttn: viprttn.o 
	$(COMP) $(CFLAGS) -o viprttn viprttn.o

viprttn.o:  viprttn.cpp
	$(COMP) $(CFLAGS) -c viprttn.cpp
