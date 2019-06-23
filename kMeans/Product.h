#ifndef __PRODUCT_H
#define __PRODUCT_H

#include "Point.h"
#include <iostream>

using namespace std;

class Cluster;

class Product
{
public:
	const static int ID_SIZE = 6;
private:
	Point point;
	char id[ID_SIZE];
	int clusterId;
	double distanceFromClusterCenter;

	Product(const Product& other);

public:
	Product(const char* id = nullptr, const int dimension = 0);

	void addValue(double value);
	void addValues(const double* values, int size);

	const Point& getPoint() const { return point; }
	Point& getPoint() { return point; }
	
	int getClusterId() const { return clusterId; }
	void setClusterId(int clusterId){ this->clusterId = clusterId; }
	
	const char* getId() const { return id; }
	char* getId() { return id; }
	void setId(char* id) { strcpy(this->id, id); }

	const double getDistanceFromClusterCenter() const { return distanceFromClusterCenter; }
	void setDistanceFromClusterCenter(double distance) { this->distanceFromClusterCenter = distance; }
	
	double& operator[](int index);

	friend ostream& operator<<(ostream& os, const Product& product);
};

#endif // __PRODUCT_H