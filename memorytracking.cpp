#include "memorytracking.h"
#include "log.h"

// Platform independent way to get available memory

#ifdef __linux__
// Unix way
#include <unistd.h>

size_t available_physical_memory()
{
	long pages = sysconf(_SC_AVPHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	return pages * page_size;
}

#endif

#ifdef __APPLE__
// Mac way only returns total, not available physical memory because of the
// way the Mac uses memory to cache some data meaning that almost all memory
// is used at all times
#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/sysctl.h>

size_t available_physical_memory()
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

size_t available_physical_memory()
{
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	GlobalMemoryStatusEx(&status);
	return status.ullAvailPhys;
}

#endif


void KKMemoryTracker::request(integer num_bytes)
{
	double reqsize = (double)num_bytes / (1024.0*1024.0*1024.0);
	num_bytes_allocated += num_bytes;
	double totalsize = (double)num_bytes_allocated / (1024.0*1024.0*1024.0);
	Output("Memory request: %.2f GB would take us to total of %.2f GB.\n", reqsize, totalsize);
	if (totalsize >= limit_gb)
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
