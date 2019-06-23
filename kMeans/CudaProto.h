#include "cuda_runtime.h"

cudaError_t calcDistancesFromClusters(const int numOfProducts, const int productSize, const int numOfClusters,
	const double* products, double* clusters, double* distancesProductsToClusters);
