#include "stdio.h"
/*Dette er en funktion som kopiere inpu til output. */
int main(void) {
  int c;
//  c = getchar();
  while (c != EOF) {
    putchar(c);
    c = getchar();
  }


  return 0;
}
