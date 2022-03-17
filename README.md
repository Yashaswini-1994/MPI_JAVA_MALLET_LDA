# DataProvenance.LDA-GA
Version of the LDA-GA algorithm for data provenance in C++

## Requirements
GCC Compiler
lnuma Library
MPICH2 implementation of Message Passing Interface
JDK 11.0.14

## MPI Cluster set up
1.    Get the remote machines (CSSMPI)
2. Need to be able to ssh without a password to other machines in order to use MPICH.
        	Run setup_mpi_cluster.sh. this creates ssh-keygen under ./.shh 
        	./setup_mpi_cluster.sh
 
        	The first time a new host is added to the first "ring"; it needs to be
        	established by a "yes" response to "continue connecting".  Every time the
        	output hangs, type "yes".
           
        	Note that setup_mpi_cluster.sh must have created this ring but just in case login    each  remote machine to check if you no longer need to type anything to get there.
        	ssh cssmpi1h
        	ssh cssmpi2h
        	ssh cssmpi3h
        	ssh cssmpi4h
        	ssh cssmpi5h
        	ssh cssmpi6h
        	ssh cssmpi7h
        	ssh cssmpi8h
 
3. Make file .mpd.conf 
        	vi/emacs/pico .mpd.conf
        	in it write one line:
        	secretword=<secretword>
        	where <secretword> is a secure key you create but not your normal password save        	the file
        	set the correct permissions on this file (other permissions won't work)
        	chmod 600 .mpd.conf
 
4. Create the mpd.hosts file in your home directory. The file should
        	include a list of cssmpi machines:
        	cssmpi2h
        	cssmpi3h
        	cssmpi4h
        	Note that you should not include cssmpi1h where you are logging in.
 
5. Edit .bash_profile file as follows:
 
 
        	For Java setup
        	export PATH=/usr/apps/mpich121-`uname -p`/bin:$PATH
        	export JAVAPATH=/usr/java/latest
        	export CLASSPATH=$CLASSPATH:/usr/apps/mpiJava-`uname -p`/lib/classes:.
        	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/apps/mpiJava-`uname -p`/lib
        	export PATH=/usr/apps/mpiJava-`uname -p`/src/scripts:$JAVAPATH/bin:$PATH
 
        	For C++ setup
        	PATH=/usr/apps/mpich121-`uname -p`/bin:$PATH
        	export PATH
 
        	either relogin or type at the command line:
        	source .bash_profile
 
6. Test that your set-up works on the current host
        	mpd &
        	mpiexec -n 1 /bin/hostname
        	mpdallexit
 
7. If you get an error or warning this is a problem. 
        	You should get the hostname of your current host
        	CSSmpdboot -n 4 -v
        	Mpdallexit
 
8. Note that, if you want to use all 8 machines, you have to list 7
        	machine names in your ~/mpd.hosts file. Then, type
        	CSSmpdboot -n 8 -v
        	Mpdallexit
 
9. To stop MPI cluster
        	Mpdallexit

## Running Java Parallel Solution
bash compileScript.sh 

To execute 
bash runScript.sh 
Number of nodes can be edited in this script to run for different number of nodes

Different data sets can be added to txtData, rawData directives along with truthfile. 
Instead processed data can be added in input0.txt file. 

# MPI_Java_Mallet_LDA
# MPI_LDA_MALLET_JAVA1
# MPI_JAVA_MALLET_LDA
