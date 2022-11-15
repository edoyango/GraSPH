#include <thrust/device_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/extrema.h>
#include <thrust/count.h>

struct is_zero
{
    __host__ __device__
    bool operator()(const int &x)
    {
        return (x == 0);
    }
};

extern "C" {
	
	void sortbykey_int_wrapper( int *keys, int *values, int N)
	{
		thrust::device_ptr<int> dev_ptr_keys(keys);
		thrust::device_ptr<int> dev_ptr_values(values);
		thrust::sort_by_key(thrust::device, dev_ptr_keys, dev_ptr_keys+N, dev_ptr_values);
	}

	void max_int_wrapper( int *data, int *maxvalue, int *maxindex, int N)
	{
		thrust::device_ptr<int> dev_ptr_data(data);
		thrust::device_ptr<int> max_ptr = thrust::max_element(thrust::device, dev_ptr_data, dev_ptr_data+N);
		*maxvalue = max_ptr[0];
		*maxindex = &max_ptr[0] - &dev_ptr_data[0] + 1; //Fortran indexing
	}

	void min_int_wrapper( int *data, int *minvalue, int *minindex, int N)
	{
		thrust::device_ptr<int> dev_ptr_data(data);
		thrust::device_ptr<int> min_ptr = thrust::min_element(thrust::device, dev_ptr_data, dev_ptr_data+N);
		*minvalue = min_ptr[0];
		*minindex = &min_ptr[0] - &dev_ptr_data[0] + 1; //Fortran indexing
	}

	void countifzero_int_wrapper( int *data, int *nzero, int N)
    {
        thrust::device_ptr<int> dev_ptr_data(data);
        *nzero = thrust::count_if(thrust::device, dev_ptr_data, dev_ptr_data+N, is_zero());
    }
	
} 