#include "memorytracking.h"
#include "log.h"

void KKMemoryTracker::request(integer num_bytes)
{
	double reqsize = (double)num_bytes / (1024.0*1024.0*1024.0);
	num_bytes_allocated += num_bytes;
	double totalsize = (double)num_bytes_allocated / (1024.0*1024.0*1024.0);
	Output("Memory request: %.2f GB would take us to total of %.2f GB.\n", reqsize, totalsize);
	if ((limit_gb > 0) && ((integer)totalsize >= limit_gb))
	{
		Error("Memory request exceeds limit.\n");
		exit(EXIT_FAILURE);
	}
};

void KKMemoryTracker::free(integer num_bytes)
{
	num_bytes_allocated -= num_bytes;
};

KKMemoryTracker memory_tracker;
