CC=g++
CFLAGS=-g -O3 -I /home/cbil/Downloads/eigen-eigen-bdd17ee3b1b3 -I /home/cbil/Henry/bamtools/include -L /home/cbil/Henry/bamtools/lib -lbamtools -lboost_system -lboost_filesystem -lpthread -lboost_thread
all: IntAPT
IntAPT: IntAPT.o filter_bam.o readinstance.o utility.o Info.o options.o merge_instances.o bedio.o rangeset.o
		$(CC) $(CFLAGS) IntAPT.o filter_bam.o readinstance.o utility.o Info.o options.o merge_instances.o bedio.o rangeset.o -o IntAPT

IntAPT.o: IntAPT.cpp process_junc.cpp mergerange.cpp
						$(CC) $(CFLAGS) -c IntAPT.cpp process_junc.cpp mergerange.cpp

filter_bam.o: filter_bam.cpp filter_bam.h
						$(CC) $(CFLAGS) -c filter_bam.cpp

merge_instances.o: merge_instances.cpp merge_instances.h
						$(CC) $(CFLAGS) -c merge_instances.cpp
						
readinstance.o: readinstance.cpp readinstance.h
						$(CC) $(CFLAGS) -c readinstance.cpp

utility.o: utility.cpp utility.h
						$(CC) $(CFLAGS) -c utility.cpp				

options.o: options.cpp options.h
						$(CC) $(CFLAGS) -c options.cpp	

Info.o: Info.cpp Info.h
						$(CC) $(CFLAGS) -c Info.cpp
						
bedio.o: bedio.cpp bedio.h
						$(CC) $(CFLAGS) -c bedio.cpp
						
rangeset.o: rangeset.cpp rangeset.h
						$(CC) $(CFLAGS) -c rangeset.cpp
						
											
clean:
	rm -rf *.o IntAPT
