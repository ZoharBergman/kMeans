#include "Cluster.h"
#include "Product.h"
#include "CudaProto.h"

#include <mpi.h>
#include <vector>

// Structs
struct sInputData
{
	int myNumOfProducts;
	int dimension;
	int maxNumOfClusters;
	int iterationLimit;
	double qualityMeasure;
};

struct sDistanceFromCluster
{
	int clusterId;
	double distance;
};

// Consts
const int YES = 1;
const int NO = 0;

sInputData readProductsFromFile(const char* filePath, Product*& products, const int numOfClusters);
void createMyInputDataMpiDataType(MPI_Datatype& inputDataMPIType);
Product* getProductById(Product* products, int numOfProducts, const char* id);
double* createAllProductsValuesArray(const Product* products, const int numOfProducts, const int dimension);
double* createAllClustersCentersValuesArray(const vector<Cluster*> clusters, const int numOfClusters, const int dimension);
void packClusters(vector<Cluster*>& clusters, int clusterBufferSize, char*& clusterBuffer, int& clusterPosition);
double calculateClustersQuality(const double* clustersDistancesFromEachOther, const double* radiusOfClusters, const int numOfClusters);
void toFsCentersOfClusters(const vector<Cluster*>& clusters, ofstream& file);
Product* getProductsByClusterId(Product* products, const int numOfProducts, const int numOfProductsInCluster, const int clusterId);
void initClusters(vector<Cluster*>& clusters, const int numOfClusters, Product* products, const int numOfProducts, const int dimension);
void writeClustersToFile(const char* filePath, const vector<Cluster*>& clusters, const double qualityMeasure);
int* divideNumberOfProductsPerProcess(const int numprocs, const int numOfProducts);
void packClustersAndProductsAndSend(const int numprocs, const int dimension, const int* numOfProductsPerProcess, 
									Product* products, vector<Cluster*>& clusters, const int numOfClusters);
void unpackClustersAndProducts(const int dimension, const int myNumOfProducts, vector<Cluster*>& clusters,
							   Product*& products, const int numOfClusters, MPI_Status& status);
sDistanceFromCluster* findMinimumDistancesProductsToClusters(const int myNumOfProducts, const double* distancesProductsToClusters, const int numOfClusters);
int transferProductsToClusters(vector<Cluster*>& clusters, const int myNumOfProducts, Product* products, sDistanceFromCluster* minimumDistancesProductsToClusters);
void packLocalClustersDataAndSend(const int dimension, const int numOfClusters, int numOfTransfers, vector<Cluster*>& clusters);
void unpackAndSetClustersCenters(vector<Cluster*>& clusters, const int dimension, const int numOfClusters, MPI_Status& status);
void calculateClustersCenters(const int numOfClusters, const int dimension, vector<Cluster*>& clusters, Point* sumOfProductsValuesPerCluster, int* sumOfNumOfProductsPerCluster);
void packAndsendClusterCenters(const int numprocs, const int numOfClusters, const int dimension, vector<Cluster*>& clusters, int isContinue);
int recieveAndHandleLocalClustersData(Point*& sumOfProductsValuesPerCluster, int*& sumOfNumOfProductsPerCluster, vector<Cluster*>& clusters, 
									  const int numOfClusters, const int dimension, const int numprocs, MPI_Status& status);
double* findRadiusOfClusters(const int numOfClusters, const int dimension, vector<Cluster*>& clusters, Product* products,const int myNumOfProducts);
void recieveAndFindRadiusOfClusters(const int numOfClusters, const int numprocs, double*& radiusOfClusters, MPI_Status& status);