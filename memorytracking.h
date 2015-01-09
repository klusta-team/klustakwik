/*
Memory tracking utilities: used to warn the user that the memory they are requesting goes over certain limits.
*/
#ifndef _MEMORY_TRACKING_H
#define _MEMORY_TRACKING_H

#include "numerics.h"
#include<vector>

size_t total_physical_memory();

class MemoryUsage
{
public:
	long long num_bytes;
	long long bytes_per_element;
	long long num_elements;
	char *name_of_type;
	char *name_of_array;
	char *expr;
	double min_multiplier, max_multiplier;
	MemoryUsage(char *_name_of_array, char *_name_of_type, long long _bytes_per_element, long long _num_elements, char *_expr, double _min_multiplier, double _max_multiplier)
	{
		name_of_array = _name_of_array;
		name_of_type = _name_of_type;
		bytes_per_element = _bytes_per_element;
		num_elements = _num_elements;
		expr = _expr;
		min_multiplier = _min_multiplier;
		max_multiplier = _max_multiplier;
		num_bytes = bytes_per_element*num_elements;
	}
};

void check_memory_usage(vector<MemoryUsage> &usages, scalar limit_gb, integer nPoints, integer nDims, integer MaxPossibleClusters);

#endif
