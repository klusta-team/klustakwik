#include "sparsearray.h"
#include<iostream>

using namespace std;

int main(int argc, char *argv[])
{
	vector<double> floatmask(3);
	floatmask[0] = 0.0;
	floatmask[1] = 1.0;
	floatmask[2] = 0.0;

	vector<double> da(3);
	da[0] = 0.0;
	da[1] = 1.0;
	da[2] = 2.0;

	Mask<int> m(floatmask, 3);
	SparseArray<double, int> sa(m, da);
	sa.unmasked_data(0) = 10.0;

	cout << "Iterate over all unmasked" << endl;
	for(int i=0; i<sa.num_unmasked; i++)
	{
		cout << sa.unmasked_index(i) << " " << sa.unmasked_data(i) << endl;
	}

	cout << "Iterate over all masked" << endl;
	for(int i=0; i<sa.num_masked; i++)
	{
		cout << sa.masked_index(i) << " " << sa.masked_data(i) << endl;
	}

	cout << "Iterate over all" << endl;
	for(int i=0; i<sa.num_all; i++)
	{
		cout << sa.all_index(i) << " " << sa.all_data(i) << endl;
	}

	return 0;
}
