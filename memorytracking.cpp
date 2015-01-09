#include "memorytracking.h"
#include "log.h"

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
