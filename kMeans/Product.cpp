#include "Product.h"

Product::Product(const char* id, const int dimension) : distanceFromClusterCenter(0)
{
	if (id != nullptr)
	{
		strcpy(this->id, id);
	}
	
	if (dimension > 0)
	{
		point.setSize(dimension);
	}
}

void Product::addValue(double value)
{
	point.addValue(value);
}

void Product::addValues(const double* values, int size)
{
	point.addValues(values, size);
}

double& Product::operator[](int index)
{
	return point[index];
}

ostream& operator<<(ostream& os, const Product& product)
{
	return os << "ID: " << product.getId() << ", Cluster ID: " << product.getClusterId() << ", Distance from cluster center: "
		<< product.getDistanceFromClusterCenter() << ", Point: {" << product.getPoint() << "}";
}