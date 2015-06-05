#include <stdio.h>

#define max(a,b) (a) - (((a) - (b)) & -((a) < (b)))

typedef  union {
      int integer;
      float real;
    } t_intreal;

int main()
{
    int a = 5;
    int b = 10;
    int c = -1;
    int d = -5;
    float e = 5.5f;
    float f = 10.5f;
    float g = -1.0f;
    float h = -5.0f;

    printf("\nmax:::%i", max(a,b));
    printf("\nmax:::%i", max(b,c));
    printf("\nmax:::%i", max(c,d));

    //printf("\nmax:::%f", max(f,e));
    //printf("\nmax:::%f", max(g,f));
    //printf("\nmax:::%f", max(h,g));

    //float i = g & (unsigned )h;
    t_intreal u_e; //5
    u_e.real = e;
    //u_e.integer = e;
    t_intreal u_f;
    u_f.real = f; //11
    //u_f.x.integer = f; //10
    //printf("\nmax:::%i\n", max(u_e.integer, u_f.integer));

    printf("\n u_f.int:::%i\n", u_f.integer);

}
