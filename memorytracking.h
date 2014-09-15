/*
Memory tracking utilities: used to warn the user that the memory they are requesting goes over certain limits.
*/
#ifndef _MEMORY_TRACKING_H
#define _MEMORY_TRACKING_H

#include "numerics.h"

size_t available_physical_memory();

class KKMemoryRequest;

class KKMemoryTracker
{
public:
	integer num_bytes_allocated;
	scalar limit_gb;
	KKMemoryTracker() { num_bytes_allocated = 0; limit_gb = 0.0; };
	void request(integer num_bytes);
	void free(integer num_bytes);
};

extern KKMemoryTracker memory_tracker;

class KKMemoryRequest
{
public:
	KKMemoryRequest() { num_bytes = 0; };
	KKMemoryRequest(integer N) : num_bytes(N) {};
	~KKMemoryRequest() { memory_tracker.free(num_bytes); }
	void add(integer N) { num_bytes += N; memory_tracker.request(N); }
	integer num_bytes;
};

#endif