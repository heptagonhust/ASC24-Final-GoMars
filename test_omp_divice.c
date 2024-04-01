#include <omp.h>
#include <stdio.h>
int main() {
	printf("omp_get_num_devices() returns: %d\n", omp_get_num_devices());
	return 0;
}
