#------------------------------------------------------------
# Target file to be compiled by default
#------------------------------------------------------------
MAIN = main
#------------------------------------------------------------
#CC is the compiler to be used
#------------------------------------------------------------
CC = gcc
#------------------------------------------------------------
#CFLAGS are the options passed to the compiler
#------------------------------------------------------------
CFLAGS = -Wall -lpthread -lrt -lm
#------------------------------------------------------------
# OBJS are the object files to be linked
#------------------------------------------------------------
OBJ1 = graphics_lib
OBJ2 = task_lib
OBJS = $(MAIN).o $(OBJ1).o $(OBJ2).o
#------------------------------------------------------------
# LIBS are the external libraries to be used
#------------------------------------------------------------
LIBS = `allegro-config --libs`
#------------------------------------------------------------
# Dependencies
#------------------------------------------------------------
$(MAIN): $(OBJS)
	$(CC) -o $(MAIN) $(OBJS) $(LIBS) $(CFLAGS)
$(MAIN).o: $(MAIN).c
	$(CC) -c $(MAIN).c  
$(OBJ1).o: $(OBJ1).c 
	$(CC) -c $(OBJ1).c 
$(OBJ2).o: $(OBJ2).c 
	$(CC) -c $(OBJ2).c
#------------------------------------------------------------
# Command that can be specified inline: make clean. Removes all the files
# created before with "make" command.
#------------------------------------------------------------
clean:
	rm -rf *.o $(MAIN) 
#------------------------------------------------------------
# Command that can be specified inline: make run. 
#------------------------------------------------------------
run:
	sudo ./$(MAIN) 
