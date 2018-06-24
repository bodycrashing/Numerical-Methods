#include<stdio.h>
int i=2; /* file scope */
void f(){printf("i=%i\n",i);}
int main(){
	int i=1; /* function scope */
	{
		int i=0; /* block scope */
		printf("i=%i\n",i);
	}
	printf("i=%i\n",i);
	f();
	return 0; }
