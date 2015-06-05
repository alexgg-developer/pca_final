#include <stdio.h>

#define max(a,b) (a) - (((a) - (b)) & -((a) < (b)))

int main()
{
    int a = 5;
    int b = 10;
    int c = -1;
    int d = -5;
    float e = 5.0f;
    float f = 10.0f;
    float g = -1.0f;
    float h = -5.0f;

    printf("\nmax:::%i", max(a,b));
    printf("\nmax:::%i", max(b,c));
    printf("\nmax:::%i", max(c,d));

    //printf("\nmax:::%f", max(f,e));
    //printf("\nmax:::%f", max(g,f));
    //printf("\nmax:::%f", max(h,g));

    //float i = g & (unsigned )h;

}
