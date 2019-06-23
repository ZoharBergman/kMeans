#ifndef __CLUSTER_H
#define __CLUSTER_H

#include <vector>
#include "Point.h"
#include "Product.h"

using namespace std; 

class Cluster
{
private:
	vector<char*> products;
	Point center;
	int id;
	double* sumOfProductsValues;
	Cluster(const Cluster& other);

	vector<char*>::const_iterator getProductItrById(const char* id) const;
	void initSumOfProductsValues(const int dimension);
	void addToSumOfProductsValues(Product& product);
	void minusFromSumOfProductsValues(Product& product);

public:
	Cluster(int id, int dimension);
	~Cluster();

	const vector<char*>& getProducts() const { return products; }
	vector<char*>& getProducts() { return products; }
	const Point& getCenter() const { return center; }	
	Point& getCenter() { return center; }
	void setCenter(const Point& center) { this->center = center; }
	const int getId() const { return id;  }
	const double* getSumOfProductsValues() const { return sumOfProductsValues; }
	double* getSumOfProductsValues() { return sumOfProductsValues; }
	void addProduct(Product& product);
	void addProducts(Product* products, const int numOfProducts);
	void removeProduct(Product& product);

	friend ostream& operator<<(ostream& os, const Cluster& cluster);
};

#endif // __CLUSTER_H