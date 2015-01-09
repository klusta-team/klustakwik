#include "memorytracking.h"
#include "log.h"
#include<algorithm>

// Platform independent way to get available memory

#ifdef __linux__
// Unix way
#include <unistd.h>

size_t total_physical_memory()
{
	long pages = sysconf(_SC_PHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	return pages * page_size;
}

#endif

#ifdef __APPLE__
#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/sysctl.h>

size_t total_physical_memory()
{
	int mib [] = { CTL_HW, HW_MEMSIZE };
	uint64_t value = 0;
	size_t length = sizeof(value);

	sysctl(mib, 2, &value, &length, NULL, 0);
	// Physical memory is now in value

	return (size_t) value;
}

#endif

#ifdef _WIN32

// Windows way

#include <windows.h>

size_t total_physical_memory()
{
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	GlobalMemoryStatusEx(&status);
	return status.ullTotalPhys;
}

#endif

bool memory_usage_comparison(const MemoryUsage &lhs, const MemoryUsage &rhs)
{
	return lhs.num_bytes > rhs.num_bytes;
}

void check_memory_usage(vector<MemoryUsage> &usages, scalar limit_gb, integer nPoints, integer nDims, integer MaxPossibleClusters)
{
	sort(usages.begin(), usages.end(), memory_usage_comparison);
	double total_min = 0.0, total_max = 0.0;
	Output("\nExpected memory usage by array (largest first):\n\n");
	for (int i = 0; i < usages.size(); i++)
	{
		MemoryUsage &m = usages[i];
		double base_usage = (double)m.num_bytes / (1024.0*1024.0*1024.0);
		double min_usage = m.min_multiplier*base_usage;
		double max_usage = m.max_multiplier*base_usage;
		total_min += min_usage;
		total_max += max_usage;
		Output("Array %s uses %lld bytes per element (%s) and has %lld elements when full. Total usage will be between %.2f and %.2f GB. Memory usage scales as %s.\n",
			m.name_of_array, m.bytes_per_element, m.name_of_type, m.num_elements, min_usage, max_usage, m.expr);
	}
	Output("\nNote that nPoints=%d, nDims=%d, MaxPossibleClusters=%d", (int)nPoints, (int)nDims, (int)MaxPossibleClusters);
	Output("\nTotal memory usage will be between %.2f and %.2f GB.\n", total_min, total_max);
	Output("RAM limit is set at %.2f GB.\n", limit_gb);
	if (total_min > limit_gb)
	{
		Error("The RAM limit does not cover the minimum possible memory usage: KlustaKwik definitely cannot run this.\n");
		exit(EXIT_FAILURE);
	} else if(total_max > limit_gb)
	{
		Error("The RAM limit covers the minimum but not maximum possible memory usage, so it cannot be guaranteed not to crash. Call KlustaKwik with -RamLimitGB -1 to force it to run.\n");
		exit(EXIT_FAILURE);
	}
}
