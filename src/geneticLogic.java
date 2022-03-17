
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Hashtable;
import java.util.List;
import java.util.Scanner;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import mpi.*;
import java.io.Serializable;
import java.util.*;
import java.util.Enumeration;
import java.util.concurrent.TimeUnit;
import cc.mallet.types.*;
import cc.mallet.pipe.*;
import cc.mallet.pipe.iterator.*;
import cc.mallet.topics.*;
import java.util.regex.*;
import java.io.*;

public class geneticLogic {

	static // a hashtable to save the fitness for values that have already been saved
	Hashtable<String, Double> fitnessTable = new Hashtable<String, Double>();

// static ArrayList<Mapclass> fitnessmapclass = new ArrayList<>();
	// the initial population of size 9
	static Integer[][] initialPopulation = new Integer[9][2];

	// to get the fitness values
	static double[] fitnessValues = new double[9];

	/**
	 * the total number of documents that are being processed. Put them in a folder
	 * and add the folder path here.
	 */
	static int numberOfDocuments; 
	static int min;
	static int max;
	static int mpi_size;
	static int mpi_rank;
	static Integer[] canBreak = new Integer[1];
	static long avg = 0;
	static int currentbesttopics =0;
	static int curentbestitr =0;
	static int numberofthreads =1;

	// public static void main(String[] args) throws IOException,
	// InterruptedException {
	public static void geneticLogic(int noOfArticles, int numdoc, int ms, int mr, int nt)
			throws IOException, InterruptedException, MPIException {

		numberOfDocuments = numdoc;
		numberofthreads = nt;
		mpi_rank = mr;
		mpi_size = ms;

		boolean maxFitnessFound = false;
	
		
			//min = (int) noOfArticles / 2;
			//max = numberOfDocuments - min;
			min = 2;
			max = 15;
			
			
			// populating the initial population
			for (int i = 0; i < initialPopulation.length; i++) {

				// the first value is the number of topics. Assign a range which you think is
				// reasonable
				// the second value is the number of iterations
				initialPopulation[i][0] = (int) Math.floor(Math.random() * max + min);
				initialPopulation[i][1] = (int) Math.floor(Math.random() * 1500 + 500);
				//ip[y++] = initialPopulation[i][0];
				//ip[y++] = initialPopulation[i][1];

				// initialPopulation[i][0] = 2;
				// initialPopulation[i][1] = 500;
			}

		
		// stop when you reach 50 iterations and take the best chromosome found till now
		int iterationNo = 0;
		int idleitr =0;
		double bestfitfoundallitr = Double.MIN_VALUE;
		while (!maxFitnessFound && iterationNo < 100) {

			// Integer[] mfc = new Integer[1];
			System.out.println("Interation No : " + iterationNo);
			Integer[][] newPopulation = new Integer[initialPopulation.length][2];
			// to get the fitness values
			fitnessValues = new double[9];
			int stripe = 9 / mpi_size;
			int start = stripe * mpi_rank;
			int end = (stripe * (mpi_rank + 1)) - 1;

			long startavg = System.currentTimeMillis();
			
			// make it map of chromosome number and fitnessval
			ArrayList<Mapclass> gtfitnessmapclass = new ArrayList<>();
			for (int i = start; i <= end; i++) {
				//System.out.println("fff loop" + i);
				Gt g = new Gt(i);
				//System.out.println("fff loop" + i + "completed");
				Mapclass mcgt = new Mapclass();
				// just recv double value of fitness
				mcgt = g.func(numberofthreads);
				// add to map i,double value
				if (mcgt != null)
					gtfitnessmapclass.add(mcgt);
			}


			//System.out.println("all ranks completed");
			Mapclass fitarray[] = new Mapclass[gtfitnessmapclass.size()];

			int h = 0;

			for (Mapclass r : gtfitnessmapclass) {
				
				fitarray[h++] = r;
				fitnessTable.put(r.s, r.d);
				if(mpi_rank==0) {
					fitnessValues[r.chromosomeNumber] = r.d;
				}
			}

			for (int r = 0; r < mpi_size; r++) {
				if (r != mpi_rank) {

					MPI.COMM_WORLD.Send(fitarray, 0, stripe, MPI.OBJECT, r, 0);
				}
			}


		

			for (int r = 0; r < mpi_size; r++) {

				if (r != mpi_rank) {
					Mapclass fitarrayrecv[] = new Mapclass[stripe];
					MPI.COMM_WORLD.Recv(fitarrayrecv, 0, stripe, MPI.OBJECT, r, 0);
					if (mpi_rank == 0) {
						for (int d = 0; d < fitarrayrecv.length; d++) {
							if (fitarrayrecv[d] != null) {
								fitnessValues[fitarrayrecv[d].chromosomeNumber] = fitarrayrecv[d].d;
							}
						}
					}

					for (Mapclass mrr : fitarrayrecv) {
						fitnessTable.put(mrr.s, mrr.d);
					}

				}

			}

			if((9%mpi_size>0)){
				if(mpi_rank==0)
				{
					Mapclass fitarrayextra[] = new Mapclass[1];
					Gt g = new Gt(8);
				
				Mapclass mcgt = new Mapclass();
				
				mcgt = g.func(numberofthreads);
				
				if (mcgt != null)
					fitarrayextra[0] = mcgt;

				fitnessTable.put(fitarrayextra[0].s,fitarrayextra[0].d);
				fitnessValues[fitarrayextra[0].chromosomeNumber] = fitarrayextra[0].d;

				for (int r = 1; r < mpi_size; r++) {
					MPI.COMM_WORLD.Send(fitarrayextra, 0, 1, MPI.OBJECT, r, 0);
				}

				} else {

					Mapclass fitarrayrecv1[] = new Mapclass[1];
					MPI.COMM_WORLD.Recv(fitarrayrecv1, 0, 1, MPI.OBJECT, 0, 0);
					fitnessTable.put(fitarrayrecv1[0].s, fitarrayrecv1[0].d);

				}
			}


			long endavg   = System.currentTimeMillis();
			
			long totalavg = endavg - startavg;
		
			avg+=totalavg;

			

			if(mpi_rank==0) {
				

				double bestfitfound = Double.MIN_VALUE;

				
				
				int idx = 0;
				for(int i = 0;i<fitnessValues.length;i++) {
					
					if(fitnessValues[i]>bestfitfound) {
						bestfitfound = fitnessValues[i];
						idx =i;
					}

				}

				

				if(bestfitfound > bestfitfoundallitr) {
					bestfitfoundallitr = bestfitfound; 
					currentbesttopics = initialPopulation[idx][0];
					curentbestitr = initialPopulation[idx][0];

				}
				else
					idleitr++;

				if(idleitr==5) {
					
					TopicModelling tm = new TopicModelling();
					tm.LDA(currentbesttopics, curentbestitr, true, -1);
					System.out.println("the best distribution is " + currentbesttopics + " topics and "
										+ curentbestitr + "iterations and fitness is " + bestfitfound);

					Ranking.createDocTopMatrix(numberOfDocuments, currentbesttopics);


					canBreak[0] = 1;
					maxFitnessFound = true;
					
				}
				else 
					canBreak[0] = 0;
				

				if(canBreak[0]== 0) {

				// copy only the top 1/3rd of the chromosomes to the new population
				for (int i = 0; i < (initialPopulation.length / 4); i++) {
					double maxFitness = Integer.MIN_VALUE;
					int maxFitnessChromosome = -1;
					for (int j = 0; j < initialPopulation.length; j++) {
						if (fitnessValues[j] > maxFitness) {
							maxFitness = fitnessValues[j];

							// CHANGE QUALITY THRESHOLD HERE
							if (maxFitness > 0.5) {
								
								TopicModelling tm = new TopicModelling();
								tm.LDA(initialPopulation[j][0], initialPopulation[j][1], true, -1);
								System.out.println("the best distribution is " + initialPopulation[j][0] + " topics and "
												+ initialPopulation[j][1] + "iterations and fitness is " + maxFitness);

								Ranking.createDocTopMatrix(numberOfDocuments, initialPopulation[j][0]);
								
								maxFitnessFound = true;
								break;
							}
							maxFitnessChromosome = j;
						}
					}

					if (maxFitnessFound) {
						break;
					}

					// copy the chromosome with high fitness to the next generation
					newPopulation[i] = initialPopulation[maxFitnessChromosome];
					fitnessValues[maxFitnessChromosome] = Double.valueOf(Integer.MIN_VALUE);
					//mfc[0] = maxFitnessChromosome;
				}

				if (maxFitnessFound) {
					canBreak[0] = 1;
				} else {
					canBreak[0] = 0;
				}

			}
			}
			if (mpi_rank == 0) {

				//System.out.println("can break sent " + canBreak[0]);
				// send canbreak values
				for (int f = 1; f < mpi_size; f++) {
					MPI.COMM_WORLD.Send(canBreak, 0, 1, MPI.OBJECT, f, 0);
				}

				//System.out.println("canbreak sent by master node");

				if (canBreak[0] == 1) {
					long a = avg/(iterationNo+1);
					System.out.println("a="+ a);
					System.out.println("TOTAL TIME FOR AVERAGE: " + a + "ms");
					break;
				} else {

					
					// perform crossover - to fill the rest of the 2/3rd of the initial Population
					for (int i = 0; i < initialPopulation.length / 4; i++) {
						newPopulation[(i + 1) * 2][0] = newPopulation[i][0];
						newPopulation[(i + 1) * 2][1] = (int) Math.floor(Math.random() * 1500 + 500);
						newPopulation[(i + 1) * 2 + 1][0] = (int) Math.floor(Math.random() * max + min);
						newPopulation[(i + 1) * 2 + 1][1] = newPopulation[i][1];
					}

					// performing crossovers based on the best combinations of the previous
					// generation
					int bestTopicNo1 = initialPopulation[0][0];
					int bestTopicNo2 = initialPopulation[1][0];
					int bestIteration = (int) Math.ceil((initialPopulation[0][1] + initialPopulation[1][1]) / 2);

					newPopulation[6][0] = (int) Math.ceil((bestTopicNo1 + bestTopicNo2) / 2);
					newPopulation[6][1] = newPopulation[7][1] = newPopulation[8][1] = bestIteration;
					newPopulation[7][0] = (bestTopicNo1 > bestTopicNo2)
							? bestTopicNo1 + (int) Math.ceil(bestTopicNo1 - bestTopicNo2) / 2
							: bestTopicNo2 + (int) Math.ceil(bestTopicNo2 - bestTopicNo1) / 2;
					newPopulation[7][0] = newPopulation[7][0] > 0 ? newPopulation[7][0] : min;
					newPopulation[7][0] = newPopulation[7][0] < max + min ? newPopulation[7][0] : max + min;
					newPopulation[8][0] = (bestTopicNo1 > bestTopicNo2)
							? bestTopicNo2 - (int) Math.ceil(bestTopicNo1 - bestTopicNo2) / 2
							: bestTopicNo1 - (int) Math.ceil(bestTopicNo2 - bestTopicNo1) / 2;
					newPopulation[8][0] = newPopulation[8][0] > 0 ? newPopulation[8][0] : min;
					newPopulation[8][0] = newPopulation[8][0] < max + min ? newPopulation[8][0] : max + min;

					// substitute the initial population with the new population and continue
					initialPopulation = newPopulation;

					
					
					  
					  Integer[] ip = new Integer[18]; int y = 0; for (int w = 0; w <
					  initialPopulation.length; w++) { for (int e = 0; e < 2; e++) { ip[y++] =
					  initialPopulation[w][e]; } }
					  
					  for (int f = 1; f < mpi_size; f++) { MPI.COMM_WORLD.Send(ip, 0, 18,
					  MPI.OBJECT, f, 0); }
					 
					 

					iterationNo++;

				}

			} else {
				MPI.COMM_WORLD.Recv(canBreak, 0, 1, MPI.OBJECT, 0, 0);
				//System.out.println("canbreak  received by slave");
				if (canBreak[0] == 1) {
					//iterationNo = 200; 
					maxFitnessFound = true;
					break;
				} else {
					

					
					  Integer[] ipr = new Integer[18]; 
					  MPI.COMM_WORLD.Recv(ipr, 0, 18, MPI.OBJECT,0, 0); 
					  int y = 0; for (int w = 0; w < 9; w++) { 
					  	for (int e = 0; e < 2; e++) {
					  initialPopulation[w][e] = ipr[y++]; 
							} 
						}
					
					 
					
					iterationNo++;
				}

			}

			

			/**
			 * The genetic algorithm loop will not exit until the required fitness is
			 * reached. For some cases, we might expect a very high fitness that will never
			 * be reached. In such cases add a variable to check how many times the GA loop
			 * is repeated. Terminate the loop in predetermined number of iterations.
			 */

		}

		if (mpi_rank == 0) {
			
			if (!maxFitnessFound) {
				long a = avg/iterationNo;
				System.out.println("TOTAL TIME FOR AVERAGE: " + a + "ms");
				// create an instance of the topic modeling class
				TopicModelling tm = new TopicModelling();
				tm.LDA(initialPopulation[0][0], initialPopulation[0][1], true, -1);
				System.out.println("the best distribution is " + initialPopulation[0][0] + " topics and "
						+ initialPopulation[0][1] + "iterations ");
				
				Ranking.createDocTopMatrix(numberOfDocuments, initialPopulation[0][0]);
				

			}

		}

	}

	public static class Gt {
		int chromosomeNo;

		public Gt(int chromosomeNo) {
			this.chromosomeNo = chromosomeNo;
		}

		public Mapclass func(int numberofthreads) {

			String combo = initialPopulation[chromosomeNo][0] + "+" + initialPopulation[chromosomeNo][1];
			
			  boolean comboPresentBefore = false;
			if (fitnessTable.containsKey(combo)) {
				Mapclass mcr = new Mapclass();
				mcr.s = combo;
				mcr.d = fitnessTable.get(combo);
				mcr.chromosomeNumber = chromosomeNo;
				return mcr;
			}
			 
			int numberOfTopics = initialPopulation[chromosomeNo][0];
			int numberOfIterations = initialPopulation[chromosomeNo][1];
			int num_th = numberofthreads;
			int divitr = numberOfIterations / num_th;
			long starttm = System.currentTimeMillis();
			


			try {

				TopicThread tml[] = new TopicThread[num_th];
				for (int th = 0; th < num_th; th++) {
					String filename = "distribution" + chromosomeNo + th;
					
					if ((numberOfIterations % num_th != 0) && th == (num_th-1)) {
						tml[th] = new TopicThread(numberOfTopics, divitr+(numberOfIterations % num_th), false, chromosomeNo, 
								filename);

					} else {
						tml[th] = new TopicThread(numberOfTopics, divitr, false, chromosomeNo,  filename);
					}

					tml[th].start();
				}

				for (int th = 0; th < num_th; th++) {
					tml[th].join();
				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			long endtm = System.currentTimeMillis();
			long tttm = endtm - starttm;
			


			int threaditr =0;
			double bestfitness = Double.MIN_VALUE;
			int documentCount = 0;
			Multimap<Integer, Integer> clusterMap = ArrayListMultimap.create();
			Scanner fileRead = null;
			double[][] clusterMatrix = new double[numberOfDocuments - 1][numberOfTopics];
			double[][] clusterCentroids = new double[numberOfTopics][numberOfTopics];
			double[] maxDistanceInsideCluster = new double[numberOfDocuments - 1];
			int rowNumber = 0, columnNumber = 0;
			double[] minDistanceOutsideCluster = new double[numberOfDocuments - 1];
			double[] silhouetteCoefficient = new double[numberOfDocuments - 1];
			double total = 0.0;
			double calc_fit = 0.0;
			while(threaditr<num_th) {

			
			clusterMatrix = new double[numberOfDocuments - 1][numberOfTopics];
			// reading the values from distribution.txt and populating the cluster matrix
			rowNumber = 0;
			columnNumber = 0;
			
			try {
				fileRead = new Scanner(new File("distribution" + chromosomeNo + threaditr +".txt"));
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			fileRead.nextLine();
			
			// Map to save the documents that belong to each cluster
			// An arraylist multimap allows to save <Key, Value[]> combination
			// each topic will have a cluster
			clusterMap = ArrayListMultimap.create();

			documentCount = 0;
			// read the values for every document
			while (documentCount < (numberOfDocuments - 1)) {

				rowNumber = fileRead.nextInt();
				fileRead.next();

				for (int z = 0; z < numberOfTopics; z++) {
					columnNumber = fileRead.nextInt();
					if (z == 0) {
						clusterMap.put(columnNumber, rowNumber);
					}
					clusterMatrix[rowNumber][columnNumber] = fileRead.nextDouble();
				}
				documentCount = documentCount + 1;
			}
			fileRead.close();
			
			// getting the centroid of each cluster by calculating the average of their
			// cluster distribution
			clusterCentroids = new double[numberOfTopics][numberOfTopics];
			for (int k : clusterMap.keySet()) {
				List<Integer> values = (List<Integer>) clusterMap.get(k);

				for (int j = 0; j < values.size(); j++) {
					int docNo = values.get(j);
					for (int y = 0; y < numberOfTopics; y++) {
						clusterCentroids[k][y] = clusterCentroids[k][y] + clusterMatrix[docNo][y];
					}
				}
				for (int y = 0; y < numberOfTopics; y++) {
					clusterCentroids[k][y] = clusterCentroids[k][y] / values.size();
				}
			}

			// finding the distance of each documents in each cluster finding max distance
			// from other documents in the same cluster
			maxDistanceInsideCluster = new double[numberOfDocuments - 1];
			for (int k : clusterMap.keySet()) {
				List<Integer> values = (List<Integer>) clusterMap.get(k);

				// for each of the documents find the maxDistance from other cluster members
				for (int y = 0; y < values.size(); y++) {
					int docNo = values.get(y);
					maxDistanceInsideCluster[docNo] = 0;
					for (int z = 0; z < values.size(); z++) {
						int otherDocNo = values.get(z);
						if (otherDocNo == docNo) {
							continue;
						}

						// finding euclidean distance between the two points/docuemnts
						double distance = 0;
						for (int h = 0; h < numberOfTopics; h++) {
							distance = distance + Math.pow((clusterMatrix[otherDocNo][h] - clusterMatrix[docNo][h]), 2);
						}
						distance = Math.sqrt(distance);
						if (distance > maxDistanceInsideCluster[docNo]) {
							maxDistanceInsideCluster[docNo] = distance;
						}
					}
				}
			}

			// finding each documents minimum distance to the centroids of other clusters
			 minDistanceOutsideCluster = new double[numberOfDocuments - 1];
			for (int k : clusterMap.keySet()) {
				List<Integer> values = (List<Integer>) clusterMap.get(k);

				// find the documents min distance from the centroid of other clusters
				for (int y = 0; y < values.size(); y++) {
					int docNo = values.get(y);
					minDistanceOutsideCluster[docNo] = Integer.MAX_VALUE;
					for (int z = 0; z < numberOfTopics; z++) {

						// don't calculate the distance to the same cluster
						if (z == k) {
							continue;
						}
						double distance = 0;
						for (int h = 0; h < numberOfTopics; h++) {
							distance = distance + Math.pow((clusterCentroids[z][h] - clusterMatrix[docNo][h]), 2);
						}
						distance = Math.sqrt(distance);
						if (distance < minDistanceOutsideCluster[docNo]) {
							minDistanceOutsideCluster[docNo] = distance;
						}
					}
				}
			}

			// calculate the Silhouette coefficient for each document
			silhouetteCoefficient = new double[numberOfDocuments - 1];
			for (int m = 0; m < (numberOfDocuments - 1); m++) {
				silhouetteCoefficient[m] = (minDistanceOutsideCluster[m] - maxDistanceInsideCluster[m])
						/ Math.max(minDistanceOutsideCluster[m], maxDistanceInsideCluster[m]);
			}

			// find the average of the Silhouette coefficient of all the documents - fitness
			// criteria
			total = 0;
			for (int m = 0; m < (numberOfDocuments - 1); m++) {
				total = total + silhouetteCoefficient[m];
			}

			calc_fit = total / (numberOfDocuments - 1);
			if(calc_fit>bestfitness)
				bestfitness = calc_fit;

			threaditr++;

		}

			

			
			
			fitnessValues[chromosomeNo] = bestfitness;
			

			// save the value in the fitnessTable to prevent LDA running for the same
			// combination again
			Mapclass mc = new Mapclass();
			mc.s = combo;
			mc.d = fitnessValues[chromosomeNo];
			mc.chromosomeNumber = chromosomeNo;

			return mc;
			

		}

	}

	public static class Mapclass implements Serializable {

		public String s;
		public Double d;
		public int chromosomeNumber;

		public void MapClass(String s, Double d) {
			this.s = s;
			this.d = d;
		}
	}

		public static class TopicThread extends Thread {
		int numberOfTopics;
		int numberOfIterations;
		boolean topicFile;
		int chromosomeNo;
		String filename;
		

		TopicThread(int numberOfTopics, int numberOfIterations, boolean topicFile, int chromosomeNo,
				String filename) {
			this.numberOfTopics = numberOfTopics;
			this.numberOfIterations = numberOfIterations;
			this.topicFile = topicFile;
			this.chromosomeNo = chromosomeNo;
			this.filename = filename;
			
		}

		public void run() {

			
			
				

			// import documents from the texts to Mallet format
		ArrayList<Pipe> pipeList = new ArrayList<Pipe>();

		// to prepare the documents for topic modelling
		// convert all the letters to lowercase
		pipeList.add(new CharSequenceLowercase());

		// tokenizing the words
		pipeList.add(new CharSequence2TokenSequence(Pattern.compile("\\p{L}[\\p{L}\\p{P}]+\\p{L}")));

		/**
		 * the stopwords file that is to be used is to be set in as the file value here.
		 */
		pipeList.add(new TokenSequenceRemoveStopwords(new File("stopwords.txt"), "UTF-8", false, false, false));

		// convert the token sequence into feature sequence
		pipeList.add(new TokenSequence2FeatureSequence());

		InstanceList instances = new InstanceList(new SerialPipes(pipeList));

			try {

			Reader fileReader = new InputStreamReader(new FileInputStream(new File("input1.txt")), "UTF-8");
			instances.addThruPipe(new CsvIterator(fileReader, Pattern.compile("^(\\S*)[\\s,]*(\\S*)[\\s,]*(.*)$"), 3, 2, 1)); 
			} catch (Exception e) {
			e.printStackTrace();
			}

			ParallelTopicModel model = new ParallelTopicModel(numberOfTopics, 0.01, 0.01);
			//System.out.println("numberOfTopics inside run" + numberOfTopics);
			if (!instances.isEmpty())
				model.addInstances(instances);

			
			//model.setNumThreads(1);
			//System.out.println("numThreads in run" + model.numThreads);
			// Run the model for 50 iterations and stop (this is for testing only,
			// for real applications, use 1000 to 2000 iterations)
			model.setNumIterations(numberOfIterations);

			
			try{
				long startes = System.currentTimeMillis();
				//System.out.println("numThreads in run" + model.numThreads);
				model.estimate();
				long endes = System.currentTimeMillis();
				long ttes = endes - startes;
				//System.out.println("TOTAL TIME func  ttes: " + ttes + "ms");

				// Show the words and topics in the first instance
				long startpr = System.currentTimeMillis();
				model.printDocumentTopics(new File(filename + ".txt"));
				long endpr = System.currentTimeMillis();
				long ttpr = endpr - startpr;
				//System.out.println("TOTAL TIME func  ttpr: " + ttpr + "ms");
			
			}
			catch (Exception e) {
			e.printStackTrace();
			}

			
	}
}
}


