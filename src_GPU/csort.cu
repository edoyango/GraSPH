#include <thrust/device_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>	

extern "C" {
	
	void sortbykey_int_wrapper( int *keys, int *values, int N)
	{
		thrust::device_ptr<int> dev_ptr_keys(keys);
		thrust::device_ptr<int> dev_ptr_values(values);
		thrust::sort_by_key(thrust::device, dev_ptr_keys, dev_ptr_keys+N, dev_ptr_values);
	}
	
} 
