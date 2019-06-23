#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>

#include <math.h>

#include <stdio.h>
#include "CudaProto.h"

#include <iostream>

using namespace std;

// Functions signatures
__global__ void initDistancesProductsToClusters(double* dev_distancesProductsToClusters, const int numOfClusters);
__global__ void calcSubtractionsPower(const double* dev_products, const double* dev_clusters, double* centersOfClustersPerProduct, const int numOfClusters, const int productSize);
__global__ void calcSumOfSubtractionsPower(double* dev_distancesProductsToClusters, const double* dev_centersOfClustersPerProduct, const int numOfClusters, const int productSize);
__global__ void doSquareRoot(double* dev_distancesProductsToClusters, const int numOfClusters);
void print(double* arr, int size);
void free(double*& dev_products, double*& dev_clusters);

cudaError_t calcDistancesFromClusters(const int numOfProducts, const int productSize, const int numOfClusters, 
													const double* products, double* clusters, double* distancesProductsToClusters)
{
	double* dev_products; // Array of all the values of all the products
	double* dev_clusters; // Array of all the values of all the clusters' centers
	double* dev_centersOfClustersPerProduct; // Array of all the clusters per product (for all the products)
	double* dev_distancesProductsToClusters; // Array of all distances from each product to each cluster's center

	cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		free(dev_products, dev_clusters);
	}

	// Allocate GPU buffers for products
	cudaStatus = cudaMalloc((void**)&dev_products, numOfProducts * productSize * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc products failed!");
		free(dev_products, dev_clusters);
	}

	// Allocate GPU buffers for clusters
	cudaStatus = cudaMalloc((void**)&dev_clusters, numOfClusters * productSize * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc clusters failed!");
		free(dev_products, dev_clusters);
	}

	// Allocate GPU buffers for centers of clusters per product
	cudaStatus = cudaMalloc((void**)&dev_centersOfClustersPerProduct, productSize * numOfProducts * numOfClusters * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc centers of clusters per product failed!");
		free(dev_products, dev_clusters);
	}

	// Allocate GPU buffers for distance from products to clusters centers
	cudaStatus = cudaMalloc((void**)&dev_distancesProductsToClusters, numOfProducts * numOfClusters * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc distance from products to clusters centers failed!");
		free(dev_products, dev_clusters);
	}

	// Copy input of the products from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_products, products, numOfProducts * productSize * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy products failed!");
		free(dev_products, dev_clusters);
	}

	// Copy input of the clusters from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_clusters, clusters, numOfClusters * productSize * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy clusters failed!");
		free(dev_products, dev_clusters);
	}

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching initExtendHistograma!\n", cudaStatus);
		free(dev_products, dev_clusters);
	}
	
	// Init the array of distances from each product to each cluster's center
	initDistancesProductsToClusters << <numOfProducts, numOfClusters >> > (dev_distancesProductsToClusters, numOfClusters);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "initDistancesProductsToClusters launch failed: %s\n", cudaGetErrorString(cudaStatus));
		free(dev_products, dev_clusters);
	}

	// For each product, calc the subtraction power of each value of the product and the value at the same index of each cluster's center
	// (product[i] - clusterCenter[i])^2 --> do this calculation for all the clusters' centers
	calcSubtractionsPower << < numOfProducts, productSize >> > (dev_products, dev_clusters, dev_centersOfClustersPerProduct, numOfClusters, productSize);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "calcSubtractionsPower launch failed: %s\n", cudaGetErrorString(cudaStatus));
		free(dev_products, dev_clusters);
	}

	// Calc the sum of all the subtracions power
	calcSumOfSubtractionsPower << < numOfProducts, numOfClusters >> > (dev_distancesProductsToClusters, dev_centersOfClustersPerProduct, numOfClusters, productSize);	

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "calcSumOfSubtractionsPower launch failed: %s\n", cudaGetErrorString(cudaStatus));
		free(dev_products, dev_clusters);
	}

	// Do square root on each sum of subtraction power (to finish the calculation of distance from each product to each cluster's center)
	doSquareRoot << <numOfProducts, numOfClusters >> > (dev_distancesProductsToClusters, numOfClusters);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "doSquareRoot launch failed: %s\n", cudaGetErrorString(cudaStatus));
		free(dev_products, dev_clusters);
	}

	// Copy the distances of each product from each cluster's center from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(distancesProductsToClusters, dev_distancesProductsToClusters, numOfProducts * numOfClusters * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy distances of each product from each cluster's center failed!");
		free(dev_products, dev_clusters);
	}

	return cudaStatus;
}


// Functions
__global__ void initDistancesProductsToClusters(double* dev_distancesProductsToClusters, const int numOfClusters)
{
	int productIndex = blockIdx.x;
	int clusterIndex = threadIdx.x;

	dev_distancesProductsToClusters[productIndex * numOfClusters + clusterIndex] = 0;
}

__global__ void calcSubtractionsPower(const double* dev_products, const double* dev_clusters, double* centersOfClustersPerProduct, const int numOfClusters, const int productSize)
{
	int productIndex = blockIdx.x;
	int valueIndex = threadIdx.x;

	for (int clusterIndex = 0; clusterIndex < numOfClusters; clusterIndex++)
	{
		centersOfClustersPerProduct[productIndex * productSize * numOfClusters + clusterIndex * productSize + valueIndex] =
			pow(dev_products[productIndex * productSize + valueIndex] - dev_clusters[clusterIndex * productSize + valueIndex], 2);
	}
}

__global__ void calcSumOfSubtractionsPower(double* dev_distancesProductsToClusters, const double* dev_centersOfClustersPerProduct,
	const int numOfClusters, const int productSize)
{
	int productIndex = blockIdx.x;
	int clusterIndex = threadIdx.x;;

	for (int valueIndex = 0; valueIndex < productSize; valueIndex++)
	{
		dev_distancesProductsToClusters[productIndex * numOfClusters + clusterIndex] += 
			dev_centersOfClustersPerProduct[productIndex * productSize * numOfClusters + productSize * clusterIndex + valueIndex];
	}
}

__global__ void doSquareRoot(double* dev_distancesProductsToClusters, const int numOfClusters)
{
	int productIndex = blockIdx.x;
	int clusterIndex = threadIdx.x;

	dev_distancesProductsToClusters[productIndex * numOfClusters + clusterIndex] = sqrt(dev_distancesProductsToClusters[productIndex * numOfClusters + clusterIndex]);
}

void free(double*& dev_products, double*& dev_clusters)
{
	cudaFree(dev_products);
	cudaFree(dev_clusters);
}

void print(double* arr, int size)
{
	for (int i = 0; i < size; i++)
	{
		printf("%d ", arr[i]);
	}
}