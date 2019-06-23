#include "mainService.h"
#include <fstream>

sInputData readProductsFromFile(const char* filePath, Product*& products, const int numOfClusters)
{
	FILE *file;
	int productIndex = 0;

	// Open the products' file
	file = fopen(filePath, "r");

	if (file)
	{
		sInputData sInputData;
		double value;
		char id[Product::ID_SIZE];

		// Read the input data
		fscanf(file, "%d", &sInputData.myNumOfProducts);
		fscanf(file, "%d", &sInputData.dimension);
		fscanf(file, "%d", &sInputData.maxNumOfClusters);
		fscanf(file, "%d", &sInputData.iterationLimit);
		fscanf(file, "%lf", &sInputData.qualityMeasure);

		// Creating the products
		products = new Product[sInputData.myNumOfProducts];

		// Reading all the products from the file
		while (fscanf(file, "%6s", &id) != EOF)
		{
			products[productIndex].setId(id);
			products[productIndex].getPoint().setSize(sInputData.dimension);

			// Add the values to the product
			for (int i = 0; i < sInputData.dimension; i++)
			{
				fscanf(file, "%lf", &value);
				products[productIndex].addValue(value);
			}

			productIndex++;
		}

		fclose(file);

		return sInputData;
	}
}

void createMyInputDataMpiDataType(MPI_Datatype& inputDataMPIType)
{
	int numOfProducts;
	int dimension;
	int maxNumOfClusters;
	int iterationLimit;
	double qualityMeasure;

	struct sInputData inputData;
	MPI_Datatype type[5] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE};
	int blocklen[5] = {1, 1, 1, 1, 1};
	MPI_Aint disp[5];

	disp[0] = (char *)&inputData.myNumOfProducts - (char *)&inputData;
	disp[1] = (char *)&inputData.dimension - (char *)&inputData;
	disp[2] = (char *)&inputData.maxNumOfClusters - (char *)&inputData;
	disp[3] = (char *)&inputData.iterationLimit - (char *)&inputData;
	disp[4] = (char *)&inputData.qualityMeasure - (char *)&inputData;

	MPI_Type_create_struct(5, blocklen, disp, type, &inputDataMPIType);
	MPI_Type_commit(&inputDataMPIType);
}

Product* getProductById(Product* products, int numOfProducts, const char* id)
{
	for (int i = 0; i < numOfProducts; i++)
	{
		if (strcmp(products[i].getId(), id) == 0)
			return &products[i];
	}

	return nullptr;
}

double* createAllProductsValuesArray(const Product* products, const int numOfProducts, const int dimension)
{
	double* allProductsValues = new double[numOfProducts * dimension];
	int position = 0;

	for (int i = 0; i < numOfProducts; i++)
	{
		const double* values = products[i].getPoint().getValues();

		for (int j = 0; j < dimension; j++)
		{
			allProductsValues[position++] = values[j];
		}
	}

	return allProductsValues;
}

double* createAllClustersCentersValuesArray(const vector<Cluster*> clusters, const int numOfClusters, const int dimension)
{
	double* allClustersCentersValues = new double[numOfClusters * dimension];
	vector<Cluster*>::const_iterator endIndex = clusters.end();
	vector<Cluster*>::const_iterator itr = clusters.begin();
	int position = 0;

	for (; itr < endIndex; ++itr)
	{
		const double* values = (*itr)->getCenter().getValues();

		for (int j = 0; j < dimension; j++)
		{
			allClustersCentersValues[position++] = values[j];
		}
	}

	return allClustersCentersValues;
}

void packClusters(vector<Cluster*>& clusters, int clusterBufferSize, char*& clusterBuffer, int& clusterPosition)
{
	double* point;
	int clusterId;

	clusterPosition = 0;

	for (int i = 0; i < clusters.size(); i++)
	{
		point = clusters.at(i)->getCenter().getValues();
		MPI_Pack(point, clusters.at(i)->getCenter().getSize(), MPI_DOUBLE, clusterBuffer, clusterBufferSize, &clusterPosition, MPI_COMM_WORLD);

		clusterId = clusters.at(i)->getId();
		MPI_Pack(&clusterId, 1, MPI_INT, clusterBuffer, clusterBufferSize, &clusterPosition, MPI_COMM_WORLD);

		//cout << "cluster " << clusters[i]->getId() << ": " << clusters[i]->getCenter() << endl;
	}
}

double calculateClustersQuality(const double* clustersDistancesFromEachOther, const double* radiusOfClusters, const int numOfClusters)
{
	double quality = 0;

	for (int clusterIndex = 0; clusterIndex < numOfClusters; clusterIndex++)
	{
		for (int otherClusterIndex = 0; otherClusterIndex < numOfClusters; otherClusterIndex++)
		{
			if (otherClusterIndex != clusterIndex)
			{
				quality += radiusOfClusters[clusterIndex] / clustersDistancesFromEachOther[clusterIndex * numOfClusters + otherClusterIndex];
			}
		}
	}

	quality /= numOfClusters * (numOfClusters - 1);

	return quality;
}

void toFsCentersOfClusters(const vector<Cluster*>& clusters, ofstream& file)
{
	vector<Cluster*>::const_iterator end = clusters.end();

	for (vector<Cluster*>::const_iterator itr = clusters.begin(); itr < end; ++itr)
	{
		file << "C" << (*itr)->getId() + 1 << "\t" << (*itr)->getCenter() << endl;
	}
}

Product* getProductsByClusterId(Product* products, const int numOfProducts, const int numOfProductsInCluster, const int clusterId)
{
	if (numOfProductsInCluster > 0)
	{
		Product* productsPerCluster = new Product[numOfProductsInCluster];
		int productsPerClusterIndex = 0;

		for (int i = 0; i < numOfProducts; i++)
		{
			if (products[i].getClusterId() == clusterId)
			{
				productsPerCluster[productsPerClusterIndex++] = products[i];
			}
		}

		return productsPerCluster;
	}
	
	return nullptr;
}

void initClusters(vector<Cluster*>& clusters, const int numOfClusters, Product* products, const int numOfProducts, const int dimension)
{
	clusters.clear();

	// Creating the clusters
	for (int i = 0; i < numOfClusters; i++)
	{
		Cluster* cluster = new Cluster(i, dimension);
		clusters.push_back(cluster);
	}

	// Putting all the products in the first cluster
	clusters.at(0)->addProducts(products, numOfProducts);
	clusters.at(0)->setCenter(getProductById(products, numOfProducts, clusters.at(0)->getProducts()[0])->getPoint());

	for (int clusterIndex = 1; clusterIndex < numOfClusters; clusterIndex++)
	{
		// Removing a product from the first cluster
		clusters.at(0)->removeProduct(products[clusterIndex]);

		// Putting the removed product in the current cluster
		clusters.at(clusterIndex)->addProduct(products[clusterIndex]);

		// Setting centers of the cluster (point of the first product)
		clusters.at(clusterIndex)->setCenter(getProductById(products, numOfProducts, clusters.at(clusterIndex)->getProducts()[0])->getPoint());
	}
}

void writeClustersToFile(const char* filePath, const vector<Cluster*>& clusters, const double qualityMeasure)
{
	ofstream file(filePath);

	if (file)
	{
		file << "K = " << clusters.size() << " QM = " << qualityMeasure << endl << "Centers of the clusters: " << endl;
		toFsCentersOfClusters(clusters, file);
		file.close();
	}
}

int* divideNumberOfProductsPerProcess(const int numprocs, const int numOfProducts)
{
	int* numOfProductsPerProcess = new int[numprocs];

#pragma omp parallel for num_threads(numprocs)
	for (int i = 0; i < numprocs; i++)
	{
		numOfProductsPerProcess[i] = numOfProducts / numprocs;
	}

#pragma omp parallel for
	for (int i = 0; i < numOfProducts % numprocs; i++)
	{
		numOfProductsPerProcess[i]++;
	}

	return numOfProductsPerProcess;
}

void packClustersAndProductsAndSend(const int numprocs, const int dimension, const int* numOfProductsPerProcess,
									Product* products, vector<Cluster*>& clusters, const int numOfClusters)
{
	int productBufferSize, clusterBufferSize, clusterPosition;
	char* productBuffer, *clusterBuffer;
	int productPosition;
	Product* currProduct;
	double* point;
	char* id;
	int clusterId;

	// Pack clusters
	clusterBufferSize = (sizeof(double) * dimension + sizeof(int)) * numOfClusters;
	clusterBuffer = new char[clusterBufferSize];
	packClusters(clusters, clusterBufferSize, clusterBuffer, clusterPosition);

	int startIndex = numOfProductsPerProcess[0];

	// Sending the products and the clusters to the processes
	for (int procIndex = 1; procIndex < numprocs; procIndex++)
	{
		productBufferSize = (sizeof(double) * dimension + sizeof(char) * Product::ID_SIZE + sizeof(int)) *  numOfProductsPerProcess[procIndex];
		productBuffer = new char[productBufferSize];
		productPosition = 0;

		// Packing the data of the products that are relevant to the current process
		for (int productIndex = startIndex; productIndex < startIndex + numOfProductsPerProcess[procIndex]; productIndex++)
		{
			currProduct = &products[productIndex];

			// Packing the point values
			point = currProduct->getPoint().getValues();
			MPI_Pack(point, dimension, MPI_DOUBLE, productBuffer, productBufferSize, &productPosition, MPI_COMM_WORLD);

			// Packing the id
			id = currProduct->getId();
			MPI_Pack(id, Product::ID_SIZE, MPI_CHAR, productBuffer, productBufferSize, &productPosition, MPI_COMM_WORLD);

			// Packing the cluster id
			clusterId = currProduct->getClusterId();
			MPI_Pack(&clusterId, 1, MPI_INT, productBuffer, productBufferSize, &productPosition, MPI_COMM_WORLD);

			clusters.at(clusterId)->removeProduct(*currProduct);
		}

		startIndex += numOfProductsPerProcess[procIndex];

		// Sending the clusters and the products
		MPI_Send(clusterBuffer, clusterPosition, MPI_PACKED, procIndex, 0, MPI_COMM_WORLD);
		MPI_Send(productBuffer, productPosition, MPI_PACKED, procIndex, 0, MPI_COMM_WORLD);

		delete[]productBuffer;
	}

	delete[]clusterBuffer;
}

void unpackClustersAndProducts(const int dimension, const int myNumOfProducts, vector<Cluster*>& clusters,
							   Product*& products, const int numOfClusters, MPI_Status& status)
{
	// Variables initialize
	int productPosition = 0;
	int productBufferSize = (sizeof(double) * dimension + sizeof(char) * Product::ID_SIZE + sizeof(int)) * myNumOfProducts;
	char* productBuffer = new char[productBufferSize];

	int clusterPosition = 0;
	int clusterBufferSize = (sizeof(double) * dimension + sizeof(int)) * numOfClusters;
	char* clusterBuffer = new char[clusterBufferSize];

	double* point = new double[dimension];
	products = new Product[myNumOfProducts];
	char* id = new char[Product::ID_SIZE];
	int clusterId;

	// Recieving the packed clusters
	MPI_Recv(clusterBuffer, clusterBufferSize, MPI_PACKED, 0, 0, MPI_COMM_WORLD, &status);

	// Recieving the packed products
	MPI_Recv(productBuffer, productBufferSize, MPI_PACKED, 0, 0, MPI_COMM_WORLD, &status);

	// Creating the clusters
	for (int i = 0; i < numOfClusters; i++)
	{
		// Unpacking the clusters message
		MPI_Unpack(clusterBuffer, clusterBufferSize, &clusterPosition, point, dimension, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Unpack(clusterBuffer, clusterBufferSize, &clusterPosition, &clusterId, 1, MPI_INT, MPI_COMM_WORLD);

		Cluster* cluster = new Cluster(clusterId, dimension);
		Point center(dimension);
		center.addValues(point, dimension);
		cluster->setCenter(center);
		clusters.push_back(cluster);
	}

	// Creating the products
	for (int i = 0; i < myNumOfProducts; i++)
	{
		// Unpacking the product message
		MPI_Unpack(productBuffer, productBufferSize, &productPosition, point, dimension, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Unpack(productBuffer, productBufferSize, &productPosition, id, Product::ID_SIZE, MPI_CHAR, MPI_COMM_WORLD);
		MPI_Unpack(productBuffer, productBufferSize, &productPosition, &clusterId, 1, MPI_INT, MPI_COMM_WORLD);

		// Setting the product
		products[i].setId(id);
		products[i].getPoint().setSize(dimension);
		products[i].addValues(point, dimension);

		// Adding the product to the cluster that it belongs to
		clusters.at(clusterId)->addProduct(products[i]);
	}

	delete[]id;
	delete[]point;
	delete[]productBuffer;
	delete[]clusterBuffer;
}

sDistanceFromCluster* findMinimumDistancesProductsToClusters(const int myNumOfProducts, const double* distancesProductsToClusters, const int numOfClusters)
{
	sDistanceFromCluster* minimumDistancesProductsToClusters = new sDistanceFromCluster[myNumOfProducts];

	// Put in the array the distance of the first cluster from each product
#pragma omp parallel for	
	for (int i = 0; i < myNumOfProducts; i++)
	{
		minimumDistancesProductsToClusters[i].distance = distancesProductsToClusters[i * numOfClusters];
		minimumDistancesProductsToClusters[i].clusterId = 0;
	}

	// For each product, find the cluster that is in the minimum distance from it
	for (int productIndex = 0; productIndex < myNumOfProducts; productIndex++)
	{
		for (int clusterIndex = 1; clusterIndex < numOfClusters; clusterIndex++)
		{
			if (distancesProductsToClusters[productIndex * numOfClusters + clusterIndex] < minimumDistancesProductsToClusters[productIndex].distance)
			{
				minimumDistancesProductsToClusters[productIndex].distance = distancesProductsToClusters[productIndex * numOfClusters + clusterIndex];
				minimumDistancesProductsToClusters[productIndex].clusterId = clusterIndex;
			}
		}
	}

	return minimumDistancesProductsToClusters;
}

int transferProductsToClusters(vector<Cluster*>& clusters, const int myNumOfProducts, Product* products, sDistanceFromCluster* minimumDistancesProductsToClusters)
{
	int numOfTransfers = 0;
	int currentClusterId, newClusterId;

	// Going over the minimum distances of the products from the clusters, and checking if we should transfer products to other clusters
	for (int productIndex = 0; productIndex < myNumOfProducts; productIndex++)
	{
		currentClusterId = products[productIndex].getClusterId();
		newClusterId = minimumDistancesProductsToClusters[productIndex].clusterId;
		products[productIndex].setDistanceFromClusterCenter(minimumDistancesProductsToClusters[productIndex].distance);
		
		if (currentClusterId != newClusterId)
		{
			// Removing the product from the old cluster and putting it in the new cluster				
			clusters.at(currentClusterId)->removeProduct(products[productIndex]);
			clusters.at(newClusterId)->addProduct(products[productIndex]);
			numOfTransfers++;
		}
	}

	return numOfTransfers;
}

void packLocalClustersDataAndSend(const int dimension, const int numOfClusters, int numOfTransfers, vector<Cluster*>& clusters)
{
	int clusterId, numOfProductsPerCluster;
	double* sumOfProductsValues;

	// Building a package
	int clustersSumProductsSize = sizeof(int) + (sizeof(int) + sizeof(int) + sizeof(double) * dimension) * numOfClusters;
	char* clustersSumProductsBuffer = new char[clustersSumProductsSize];
	int clustersSumProductsPosition = 0;

	// Packing the number of transfers from clusters 
	MPI_Pack(&numOfTransfers, 1, MPI_INT, clustersSumProductsBuffer, clustersSumProductsSize, &clustersSumProductsPosition, MPI_COMM_WORLD);

	// Packing the clusters' IDs, number of products per cluster and the clusters' products' values' sum
	for (int i = 0; i < numOfClusters; i++)
	{
		clusterId = clusters.at(i)->getId();
		MPI_Pack(&clusterId, 1, MPI_INT, clustersSumProductsBuffer, clustersSumProductsSize, &clustersSumProductsPosition, MPI_COMM_WORLD);

		numOfProductsPerCluster = clusters.at(i)->getProducts().size();
		MPI_Pack(&numOfProductsPerCluster, 1, MPI_INT, clustersSumProductsBuffer, clustersSumProductsSize, &clustersSumProductsPosition, MPI_COMM_WORLD);

		sumOfProductsValues = clusters.at(i)->getSumOfProductsValues();
		MPI_Pack(sumOfProductsValues, dimension, MPI_DOUBLE, clustersSumProductsBuffer, clustersSumProductsSize, &clustersSumProductsPosition, MPI_COMM_WORLD);
	}

	// Sending the package to process 0
	MPI_Send(clustersSumProductsBuffer, clustersSumProductsPosition, MPI_PACKED, 0, 0, MPI_COMM_WORLD);

	delete[]clustersSumProductsBuffer;
}

void unpackAndSetClustersCenters(vector<Cluster*>& clusters, const int dimension, const int numOfClusters, MPI_Status& status)
{
	// Init variables
	int clusterPosition = 0;
	int clusterBufferSize = (sizeof(double) * dimension + sizeof(int)) * numOfClusters;
	char* clusterBuffer = new char[clusterBufferSize];
	
	double* point = new double[dimension];
	int clusterId;

	// Recieving the new clusters' centers from process 0
	MPI_Recv(clusterBuffer, clusterBufferSize, MPI_PACKED, 0, 0, MPI_COMM_WORLD, &status);

	// Setting the clusters' centers
	for (int clusterIndex = 0; clusterIndex < numOfClusters; clusterIndex++)
	{
		// Unpacking the clusters message
		MPI_Unpack(clusterBuffer, clusterBufferSize, &clusterPosition, point, dimension, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Unpack(clusterBuffer, clusterBufferSize, &clusterPosition, &clusterId, 1, MPI_INT, MPI_COMM_WORLD);

		clusters.at(clusterId)->getCenter().setValues(point);
	}

	delete[]clusterBuffer;
	delete[]point;
}

void calculateClustersCenters(const int numOfClusters, const int dimension, vector<Cluster*>& clusters, Point* sumOfProductsValuesPerCluster, int* sumOfNumOfProductsPerCluster)
{
#pragma omp parallel for
	for (int clusterIndex = 0; clusterIndex < numOfClusters; clusterIndex++)
	{
#pragma omp parallel for
		for (int valueIndex = 0; valueIndex < dimension; valueIndex++)
		{
			sumOfProductsValuesPerCluster[clusterIndex][valueIndex] /= sumOfNumOfProductsPerCluster[clusterIndex];
		}

		// Setting the new center of the cluster
		clusters.at(clusterIndex)->setCenter(sumOfProductsValuesPerCluster[clusterIndex]);
	}
}

void packAndsendClusterCenters(const int numprocs, const int numOfClusters, const int dimension, vector<Cluster*>& clusters, int isContinue)
{
	int clusterPosition;
	int clusterBufferSize = (sizeof(double) * dimension + sizeof(int)) * numOfClusters;
	char* clusterBuffer = new char[clusterBufferSize];

	// Packing the clusters
	packClusters(clusters, clusterBufferSize, clusterBuffer, clusterPosition);

	// Sending to the other processes a flag that they should continue to search for the clusters 
	// and the new centers of the clusters
	for (int procIndex = 1; procIndex < numprocs; procIndex++)
	{
		MPI_Send(&isContinue, 1, MPI_INT, procIndex, 0, MPI_COMM_WORLD);
		MPI_Send(clusterBuffer, clusterPosition, MPI_PACKED, procIndex, 0, MPI_COMM_WORLD);
	}

	delete[]clusterBuffer;
}

int recieveAndHandleLocalClustersData(Point*& sumOfProductsValuesPerCluster, int*& sumOfNumOfProductsPerCluster, vector<Cluster*>& clusters, 
									  const int numOfClusters, const int dimension, const int numprocs, MPI_Status& status)
{
	int sumOfNumOfTransfers = 0;
	double* sumOfProductsValues;
	int numOfTransfers, clusterId, numOfProductsPerCluster;

	// Init the sum of values of product per cluster  and 
	// the sum of num of products per cluster  with the data of process 0
#pragma omp parallel for
	for (int clusterIndex = 0; clusterIndex < numOfClusters; clusterIndex++)
	{
		// Init the sum of values of product per cluster
		sumOfProductsValuesPerCluster[clusterIndex].setSize(dimension);
		sumOfProductsValues = clusters.at(clusterIndex)->getSumOfProductsValues();
		sumOfProductsValuesPerCluster[clusterIndex].addValues(sumOfProductsValues, dimension);

		// Init the sum of num of products per cluster
		sumOfNumOfProductsPerCluster[clusterIndex] = clusters.at(clusterIndex)->getProducts().size();
	}

	int clustersSumProductsPosition;
	int clustersSumProductsSize = sizeof(int) + (sizeof(int) + sizeof(int) + sizeof(double) * dimension) * numOfClusters;
	char* clustersSumProductsBuffer = new char[clustersSumProductsSize];
	sumOfProductsValues = new double[dimension];

	// Recieving packages from all the processes
	for (int procIndex = 1; procIndex < numprocs; procIndex++)
	{
		clustersSumProductsPosition = 0;

		// Recieving the number of transfers from clusters, the clusters' IDs and the clusters' products' values' sum
		MPI_Recv(clustersSumProductsBuffer, clustersSumProductsSize, MPI_PACKED, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

		// Unpacking the number of transfers
		MPI_Unpack(clustersSumProductsBuffer, clustersSumProductsSize, &clustersSumProductsPosition, &numOfTransfers, 1, MPI_INT, MPI_COMM_WORLD);
		
		sumOfNumOfTransfers += numOfTransfers;

		for (int clusterIndex = 0; clusterIndex < numOfClusters; clusterIndex++)
		{
			// Unpacking the clusters' IDs, number of products per cluster and the clusters' products' values' sum
			MPI_Unpack(clustersSumProductsBuffer, clustersSumProductsSize, &clustersSumProductsPosition, &clusterId, 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(clustersSumProductsBuffer, clustersSumProductsSize, &clustersSumProductsPosition, &numOfProductsPerCluster, 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(clustersSumProductsBuffer, clustersSumProductsSize, &clustersSumProductsPosition, sumOfProductsValues, dimension, MPI_DOUBLE, MPI_COMM_WORLD);

			// Summing up the number of products in the current cluster
			sumOfNumOfProductsPerCluster[clusterIndex] += numOfProductsPerCluster;

			// Summing up the values of products in the current cluster						
#pragma omp parallel for
			for (int valueIndex = 0; valueIndex < dimension; valueIndex++)
			{
				sumOfProductsValuesPerCluster[clusterIndex][valueIndex] += sumOfProductsValues[valueIndex];
			}
		}
	}

	delete[]sumOfProductsValues;
	delete[]clustersSumProductsBuffer;

	return sumOfNumOfTransfers;
}

double* findRadiusOfClusters(const int numOfClusters, const int dimension, vector<Cluster*>& clusters, Product* products, const int myNumOfProducts)
{
	double* radiusOfClusters = new double[numOfClusters];
	int numOfProductsPerCluster;
	
	// In each cluster, searching for the product that is the farest from the cluster center
#pragma omp parallel for private(numOfProductsPerCluster)
	for (int clusterIndex = 0; clusterIndex < numOfClusters; clusterIndex++)
	{
		numOfProductsPerCluster = clusters.at(clusterIndex)->getProducts().size();

		// Checking if there are products in the current cluster
		if (numOfProductsPerCluster == 0)
		{
			radiusOfClusters[clusterIndex] = 0;
		}
		else
		{
			// Getting the products of the current cluster
			Product* productsOfCluster = 
				getProductsByClusterId(products, myNumOfProducts, numOfProductsPerCluster, clusters.at(clusterIndex)->getId());
			
			// Init the max distance to the first product
			radiusOfClusters[clusterIndex] = productsOfCluster[0].getDistanceFromClusterCenter();

			// Finding the max distance - the radius of the cluster
			for (int i = 1; i < numOfProductsPerCluster; i++)
			{
				if (productsOfCluster[i].getDistanceFromClusterCenter() > radiusOfClusters[clusterIndex])
				{
					radiusOfClusters[clusterIndex] = productsOfCluster[i].getDistanceFromClusterCenter();
				}
			}

			delete[]productsOfCluster;
	}

	return radiusOfClusters;
}

void recieveAndFindRadiusOfClusters(const int numOfClusters, const int numprocs, double*& radiusOfClusters, MPI_Status& status)
{
	double* radiusOfClusterOfProcesses = new double[numOfClusters];

	// Recieving from the other processes the radius of each of their local clusters
	for (int procIndex = 1; procIndex < numprocs; procIndex++)
	{
		MPI_Recv(radiusOfClusterOfProcesses, numOfClusters, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

		// Find the max radius for each of the clusters
#pragma omp parallel for
		for (int clusterIndex = 0; clusterIndex < numOfClusters; clusterIndex++)
		{
			if (radiusOfClusterOfProcesses[clusterIndex] > radiusOfClusters[clusterIndex])
			{
				radiusOfClusters[clusterIndex] = radiusOfClusterOfProcesses[clusterIndex];
			}
		}
	}

	delete[]radiusOfClusterOfProcesses;
}