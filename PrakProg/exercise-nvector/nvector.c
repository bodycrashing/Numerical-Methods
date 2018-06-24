#include "nvector.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/*--------- Initiliazing vectors -----------*/
nvector *nvector_alloc(int n)
{
	nvector *v = (nvector *) malloc(sizeof(nvector));
	(*v).size = n;
	(*v).data = (double *)malloc(n * sizeof(double));
	if (v == NULL)
		fprintf(stderr, "error in nvector_alloc\n");
	return v;
}

void nvector_free(nvector * v)
{
	free(v->data);
	free(v);
}


/*--------- Accessing and setting vector elements -----------*/
void nvector_set(nvector * v, int i, double value)
{
	assert(0 <= i && i < (*v).size);
	(*v).data[i] = value;
}

double nvector_get(nvector * v, int i)
{
	assert(0 <= i && i < (*v).size);
	return v->data[i];
}


/*--------- Vector print utility -----------*/
void nvector_print(char *s, nvector * v)
{
	printf("%s", s);
	for (int i = 0; i < v->size; i++)
		printf("%9.3g ", v->data[i]);
	printf("\n");
}



/*--------- Vector Operations -----------*/
void nvector_add(nvector * a, nvector * b)
{
	assert(a->size == b->size);
	for (int i = 0; i < a->size; i++) {
		double s = nvector_get(a, i) + nvector_get(b, i);
		nvector_set(a, i, s);
	}
}

int double_equal(double a, double b)
{
	double TAU = 1e-6, EPS = 1e-6;
	if (fabs(a - b) < TAU)
		return 1;
	if (fabs(a - b) / (fabs(a) + fabs(b)) < EPS / 2)
		return 1;
	return 0;
}

int nvector_equal(nvector * a, nvector * b)
{
	if (a->size != b->size) return 0;
	for (int i = 0; i < a->size; i++)
		if (!double_equal(a->data[i], b->data[i]))
			return 0;
	return 1;
}
