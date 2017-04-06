# BayeScan makefile

bayescan_2.1: start.o beta.o dirichlet.o RJupdates.o MHupdates.o likelihood.o read_write.o anyoption.o 
	g++ -fopenmp -lpthread -lgomp -o bayescan_2.1 start.o beta.o dirichlet.o RJupdates.o MHupdates.o likelihood.o read_write.o anyoption.o 

start.o: start.cpp errors.cpp anyoption.h global_defs.h
	g++ -fopenmp -c start.cpp errors.cpp 

beta.o: beta.cpp global_defs.h
	g++ -fopenmp -c beta.cpp 
      
dirichlet.o: dirichlet.cpp global_defs.h
	g++ -fopenmp -c dirichlet.cpp 

RJupdates.o: RJupdates.cpp global_defs.h
	g++ -fopenmp -c RJupdates.cpp 

MHupdates.o: MHupdates.cpp global_defs.h
	g++ -fopenmp -c MHupdates.cpp 

likelihood.o: likelihood.cpp global_defs.h
	g++ -fopenmp -c likelihood.cpp 

read_write.o: read_write.cpp errors.cpp global_defs.h
	g++ -fopenmp -c read_write.cpp errors.cpp 

anyoption.o: anyoption.cpp anyoption.h 
	g++ -fopenmp -c anyoption.cpp 

clean: 
	rm *.o bayescan_2.1
