#include "Point.h"

Point::Point(int size) : physicalSize(size), logicalSize(0), values(nullptr)
{
	if (size > 0)
	{
		values = new double[size];
	}
}

void Point::addValue(double value)
{
	if (logicalSize < physicalSize)
		values[logicalSize++] = value;
	else
		throw ("Point is full");
}

void Point::addValues(const double* values, int size)
{
	for (int i = 0; i < size; i++)
	{
		addValue(values[i]);
	}
}

void Point::setSize(int size)
{
	physicalSize = size;
	logicalSize = 0;

	delete []values;

	values = new double[size];
}

double& Point::operator[](int index)
{
	if (index < physicalSize)
		return values[index];
	else
		throw ("Index is out of range");
}

void Point::setValues(const double* values)
{
	for (int i = 0; i < logicalSize; i++)
	{
		this->values[i] = values[i];
	}
}

const Point& Point::operator=(const Point& other)
{
	if (this != &other)
	{	
		delete []values;

		physicalSize = other.physicalSize;
		logicalSize = other.logicalSize;		

		values = new double[physicalSize];

		for (int i = 0; i < logicalSize; i++)
		{
			values[i] = other.values[i];
		}
	}

	return *this;
}

ostream& operator<<(ostream& os, const Point& point)
{
	const double* values = point.getValues();

	for (int i = 0; i < point.getLogicalSize(); i++)
	{
		os << values[i] << "\t";
	}

	return os;
}
