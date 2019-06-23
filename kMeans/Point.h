#ifndef __POINT_H
#define __POINT_H

#include <iostream>

using namespace std;

class Point
{
private:
	double* values;
	int physicalSize;
	int logicalSize;
	
public:
	Point(int size = 0);
	Point(const Point& other) { *this = other; }
	virtual ~Point(){ delete []values; }

	void addValue(double value);
	void addValues(const double* values, int size);
	int getSize() const { return logicalSize; }
	const int getLogicalSize() const { return logicalSize; }
	void setSize(int size);
	const double* getValues() const { return values; }
	double* getValues() { return values; }
	void setValues(const double* values);

	const Point& operator=(const Point& other);
	double& operator[](int index);
	friend ostream& operator<<(ostream& os, const Point& point);
};

#endif // __POINT_H