#include<stdio.h>
void f1(int i){i=0;}
void f2(int* i){*i=0;}
void f3(int* i){i=NULL;}
int main(){
	int i=1; f1(i); printf("i=%i\n",i);

  i=1; f2(&i); printf("i=%i\n",i);

	i=1; f3(&i); printf("i=%i\n",i);
              
	return 0; }
