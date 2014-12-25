#include<vector>

using namespace std;

template<class index_type>
class Mask
{
public:
	vector<index_type> unmasked_indices;
	vector<index_type> default_indices;
	Mask() {}
	template<class data_type>
	Mask(data_type *x, index_type n)
	{
		from_flat_mask(x, n);
	}
	template<class data_type>
	Mask(vector<data_type> x, index_type n, index_type offset=0)
	{
		from_flat_mask(x, n, offset);
	}
	template<class data_type>
	void from_flat_mask(data_type *x, index_type n)
	{
		for(index_type i=0; i<n; i++)
		{
			if(x[i]!=(data_type)0)
				unmasked_indices.push_back(i);
			else
				default_indices.push_back(i);
		}
	}
	template<class data_type>
	void from_flat_mask(vector<data_type> x, index_type n, index_type offset=0)
	{
		from_flat_mask<data_type>(&(x[offset]), n);
	}
};

template<class data_type, class index_type>
class SparseArray
{
public:
	vector<data_type> unmasked_data_array;
	vector<data_type> *default_data_array;
	Mask<index_type> *mask;
	index_type num_all, num_masked, num_unmasked;

	SparseArray(Mask<index_type> &_mask, vector<data_type> &_default)
	{
		mask = &_mask;
		default_data_array = &_default;
		num_masked = mask->default_indices.size();
		num_unmasked = mask->unmasked_indices.size();
		num_all = num_masked+num_unmasked;
		unmasked_data_array.resize(num_unmasked);
	}

	// access unmasked part only
	inline data_type& unmasked_data(index_type i) { return unmasked_data_array[i]; }
	inline index_type& unmasked_index(index_type i) { return mask->unmasked_indices[i]; }
	// access masked part only
	inline data_type& masked_data(index_type i) { return (*default_data_array)[masked_index(i)]; }
	inline index_type& masked_index(index_type i) { return mask->default_indices[i]; }
	// access either part (less efficient)
	inline data_type& all_data(index_type i)
	{
		if(i<num_unmasked)
			return unmasked_data(i);
		else
			return masked_data(i-num_unmasked);
	}
	inline index_type& all_index(index_type i)
	{
		if(i<num_unmasked)
			return unmasked_index(i);
		else
			return masked_index(i-num_unmasked);
	}
};
