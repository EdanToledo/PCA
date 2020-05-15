EXECUTABLE_NAME=RunPCA
PROGRAM_Driver=RunPCA
PREFILE1=PCA





CC=g++
COMPILE=-c -o
LINK=-o
LIBRARY=-I eigen-3.3.7/

FLAGS=-std=c++11 -g

$(PROGRAM_Driver): $(PROGRAM_Driver).o
	$(CC) $(FLAGS) $(LIBRARY) $(LINK) $(EXECUTABLE_NAME) $(PROGRAM_Driver).o $(PREFILE1).o
$(PROGRAM_Driver).o: $(PREFILE1).o $(PREFILE1).h
	$(CC) $(FLAGS) $(LIBRARY) $(COMPILE) $(PROGRAM_Driver).o $(PROGRAM_Driver).cpp

$(PREFILE1).o: $(PREFILE1).h
	$(CC) $(FLAGS) $(LIBRARY) $(COMPILE) $(PREFILE1).o $(PREFILE1).cpp


clean:
	rm -f *.o
	rm -f $(EXECUTABLE_NAME)
