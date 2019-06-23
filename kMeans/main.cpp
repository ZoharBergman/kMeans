#include <iostream>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <vector>

#include "mainService.h"
#include "CudaProto.h"
#include "Cluster.h"
#include "Product.h"

using namespace std;
int main(int argc, char *argv[])
{
	// Variables
	int myid, numprocs, errorCode = 999, isContinue;
	MPI_Datatype inputDataMPIType;
	double t1, t2;

	vector<Cluster*> clusters;
	Product* products = nullptr;
	int numOfClusters = 2;
	sInputData inputData;
	int numOfProducts;
	double qualityMeasure;
	
	Point* sumOfProductsValuesPerCluster;
	int* sumOfNumOfProductsPerCluster;

	double* distancesProductsToClusters;
	sDistanceFromCluster* minimumDistancesProductsToClusters;
	double* radiusOfClusters;

	int numOfTransfers, sumOfNumOfTransfers;

	double* allClustersCentersValues, *allProductsValues;
	double* clustersDistancesFromEachOther;

	// Init MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	// Creating MPI data type of sInputData
	createMyInputDataMpiDataType(inputDataMPIType);

	if (myid == 0)
	{
		const char* filePath = "D:\\Zohar\\15.11.17\\Sales_Transactions_Dataset_Weekly_part.txt";

		// Reading the products from file
		inputData = readProductsFromFile(filePath, products, numOfClusters);		
	}

	// Broad casting the input data to all the processes
	MPI_Bcast(&inputData, 1, inputDataMPIType, 0, MPI_COMM_WORLD);

	// Init the quality Measure and the num of products
	qualityMeasure = inputData.qualityMeasure;
	numOfProducts = inputData.myNumOfProducts;

	// While we didn't find qualitative clusters and didn't reach for the max number of clusters yet
	while (qualityMeasure >= inputData.qualityMeasure && numOfClusters <= inputData.maxNumOfClusters)
	{
		clusters.clear();

		if (myid == 0)
		{
			t1 = MPI_Wtime();

			// Init the clusters
			initClusters(clusters, numOfClusters, products, numOfProducts, inputData.dimension);

			// Dividing the number of products to be handle between the processes
			int* numOfProductsPerProcess = divideNumberOfProductsPerProcess(numprocs, numOfProducts);

			// Sending to all the other processes the number of products that they should be handle with
			for (int i = 1; i < numprocs; i++)
			{
				MPI_Send(numOfProductsPerProcess + i, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			}

			// Setting the number of products that process 0 should be handle with
			inputData.myNumOfProducts = numOfProductsPerProcess[0];

			// Packing the clusters and the products and send them to the other processes
			packClustersAndProductsAndSend(numprocs, inputData.dimension, numOfProductsPerProcess, 
										   products, clusters, numOfClusters);
			
			delete[]numOfProductsPerProcess;
		}
		else // Other processes
		{
			delete[]products;
			
			// Recieving my number of products from process 0
			MPI_Recv(&inputData.myNumOfProducts, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

			// Recieving the clusters and my products from process 0 and unpacking it
			unpackClustersAndProducts(inputData.dimension, inputData.myNumOfProducts, clusters, products, numOfClusters, status);
		}

		isContinue = YES;

		// Trying to find clusters
		for (int iterationNum = 0; iterationNum < inputData.iterationLimit && isContinue == YES; iterationNum++)
		{
			distancesProductsToClusters = new double[inputData.myNumOfProducts * numOfClusters];
			allProductsValues = createAllProductsValuesArray(products, inputData.myNumOfProducts, inputData.dimension);
			allClustersCentersValues = createAllClustersCentersValuesArray(clusters, numOfClusters, inputData.dimension);			
			
			// Calculate the distances of all products from all of the clusters' centers
			calcDistancesFromClusters(inputData.myNumOfProducts, inputData.dimension, numOfClusters, 
									  allProductsValues, allClustersCentersValues, distancesProductsToClusters);
			
			delete[]allProductsValues;
			delete[]allClustersCentersValues;
			
			// For each product, find the cluster that is in the minimum distance from it
			minimumDistancesProductsToClusters = findMinimumDistancesProductsToClusters(
														inputData.myNumOfProducts, distancesProductsToClusters, numOfClusters);

			delete[]distancesProductsToClusters;

			// Transfer products that should be transfer to other clusters
			numOfTransfers = transferProductsToClusters(clusters, inputData.myNumOfProducts, products, minimumDistancesProductsToClusters);

			delete[]minimumDistancesProductsToClusters;

			if (myid != 0)
			{
				// Packing the local data of the clusters and send it to process 0
				packLocalClustersDataAndSend(inputData.dimension, numOfClusters, numOfTransfers, clusters);

				// Getting a message from process 0 - if we should continue searching for the clusters
				MPI_Recv(&isContinue, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
				
				if (isContinue == YES)
				{
					// Recieving the new clusters' centers and setting them
					unpackAndSetClustersCenters(clusters, inputData.dimension, numOfClusters, status);
				}
			}
			else // Process 0
			{
				// Allocating the sums: products values per cluster and the number of products per cluster
				sumOfProductsValuesPerCluster = new Point[numOfClusters];
				sumOfNumOfProductsPerCluster = new int[numOfClusters];

				// Recieving the local clusters data from the other processes and sum it up
				sumOfNumOfTransfers = recieveAndHandleLocalClustersData(sumOfProductsValuesPerCluster, sumOfNumOfProductsPerCluster,
					clusters, numOfClusters, inputData.dimension, numprocs, status);

				// Adding the number of transfers of process 0
				sumOfNumOfTransfers += numOfTransfers;

				// Checking if there were transfers of products between the clusters
				if (sumOfNumOfTransfers == 0)
				{
					isContinue = NO;

					// Sending to the other processes a flag that they should not continue to search for the clusters 
					for (int procIndex = 1; procIndex < numprocs; procIndex++)
					{
						MPI_Send(&isContinue, 1, MPI_INT, procIndex, 0, MPI_COMM_WORLD);
					}
				}
				else
				{
					// Calculating the new clusters' centers
					calculateClustersCenters(numOfClusters, inputData.dimension, clusters, 
													sumOfProductsValuesPerCluster, sumOfNumOfProductsPerCluster);

					// Packing and sending the clusters to the other processes
					packAndsendClusterCenters(numprocs, numOfClusters, inputData.dimension, clusters, isContinue);
				}

				delete[]sumOfProductsValuesPerCluster;
				delete[]sumOfNumOfProductsPerCluster;
			}

			if (isContinue == NO)
			{
				// For each cluster, finding the product that is the farest from the cluster's center and 
				// setting that distance as the local cluster radius
				radiusOfClusters = findRadiusOfClusters(numOfClusters, inputData.dimension, clusters, products, inputData.myNumOfProducts);
				
				if (myid != 0)
				{
					// Sending the radius of the local clusters to process 0 
					MPI_Send(radiusOfClusters, numOfClusters, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
				}
				else // Process 0
				{
					// Recieving the radius of the clusters from the other processes and finding the radius of each cluster
					recieveAndFindRadiusOfClusters(numOfClusters, numprocs, radiusOfClusters, status);

					allClustersCentersValues = createAllClustersCentersValuesArray(clusters, numOfClusters, inputData.dimension);
					clustersDistancesFromEachOther = new double[numOfClusters * numOfClusters];

					// Calculating the distance of each cluster from the other clusters
					calcDistancesFromClusters(numOfClusters, inputData.dimension, numOfClusters,
						allClustersCentersValues, allClustersCentersValues, clustersDistancesFromEachOther);

					delete[]allClustersCentersValues;

					// Calculating the clusters quality
					qualityMeasure = calculateClustersQuality(clustersDistancesFromEachOther, radiusOfClusters, numOfClusters);

					delete[]clustersDistancesFromEachOther;
				}
				
				MPI_Bcast(&qualityMeasure, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				numOfClusters++;

				delete[]radiusOfClusters;
			}
		}
	}

	if (myid == 0)
	{
		t2 = MPI_Wtime();
		writeClustersToFile("D:\\Zohar\\15.11.17\\Result.txt", clusters, qualityMeasure);
		printf("Time measured: %1.6f\n", t2 - t1); fflush(stdout);
	}

	delete[]products;

	MPI_Finalize();
	return 0;
}