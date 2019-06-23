#include "Cluster.h"

Cluster::Cluster(int id, int dimension) : id(id)
{
	center.setSize(dimension);
	sumOfProductsValues = new double[dimension];
	initSumOfProductsValues(dimension);
}

Cluster::~Cluster()
{
	delete[]sumOfProductsValues;
}

void Cluster::addProduct(Product& product)
{
	product.setClusterId(id);
	products.push_back(product.getId());	
	addToSumOfProductsValues(product);
}

void Cluster::addProducts(Product* products, const int numOfProducts)
{
	for (int i = 0; i < numOfProducts; i++)
	{
		addProduct(products[i]);
	}
}

void Cluster::removeProduct(Product& product)
{
	minusFromSumOfProductsValues(product);
	products.erase(getProductItrById(product.getId()));
}

vector<char*>::const_iterator Cluster::getProductItrById(const char* id) const
{
	vector<char*>::const_iterator end = products.end();

	for (vector<char*>::const_iterator itr = products.begin(); itr < end; ++itr)
	{
		if (strcmp(*itr, id) == 0)
			return itr;
	}

	return products.end();
}

void Cluster::initSumOfProductsValues(const int dimension)
{
#pragma omp parallel for
	for (int i = 0; i < dimension; i++)
	{
		sumOfProductsValues[i] = 0;
	}
}

void Cluster::addToSumOfProductsValues(Product& product)
{
	double* point = product.getPoint().getValues();
#pragma omp parallel for
	for (int i = 0; i < product.getPoint().getSize(); i++)
	{
		sumOfProductsValues[i] += point[i];
	}
}

void Cluster::minusFromSumOfProductsValues(Product& product)
{
	double* point = product.getPoint().getValues();
#pragma omp parallel for
	for (int i = 0; i < product.getPoint().getSize(); i++)
	{
		sumOfProductsValues[i] -= point[i];
	}
}

ostream& operator<<(ostream& os, const Cluster& cluster)
{
	os << "ID: " << cluster.getId() << ", Center: {" << cluster.getCenter() << "}, Products: ";

	vector<char*> products = cluster.getProducts();
	vector<char*>::const_iterator end = products.end();

	for (vector<char*>::const_iterator itr = products.begin(); itr < end; ++itr)
	{
		os << *itr << ", ";
	}

	os << "}";

	return os;
}