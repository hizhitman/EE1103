 
#include <stdio.h>
#include <math.h>
#include <time.h>
/*
double fact(double n) //FOR EXTRA PRECISION WE CAN UNCOMMENT THIS
{
    int factorial=1;
    for (int i=1; i<=n;i++){
        factorial=factorial*i;
    }
    return factorial;
}
*/
double f(double x, double y)
{
    return  -100000*y + 99999*pow(M_E,-x); //RETURNS F(X) VALUE
}
/*
double f1(double x, double y)//FOR EXTRA PRECISION WE CAN UNCOMMENT THIS
{
    return  3*y*pow(x,2);
}
double f2(double x, double y)//FOR EXTRA PRECISION WE CAN UNCOMMENT THIS
{
    return  y*6*x;
}
double f3(double x, double y)//FOR EXTRA PRECISION WE CAN UNCOMMENT THIS
{
    return  6*y;
}
*/
int main()
{
    double h;
    printf("BY EXPLICIT EULER'S METHOD\n");
    printf("Enter step size value h: "); //THAT VALUE OF H FOR WHICH THE ANSWER IS STABLE,LESS THAN 0.00002
    
    scanf("%lf", &h);
    printf("h=%lf\n",h);
    printf("\n");
    double x_low=0;//EQUIVALENT TO X_I
    double y_low=0,y_high;//EQUIVALENT TO Y_I AND Y_I+1
    FILE *fptr2;
    fptr2 = fopen("eee.txt", "a");//OPENING FILE TO INPUT VALUES AS ORDERED PAIRS
    printf("X\t\tY\n");
    printf("\n");
    fprintf(fptr2,"%lf \t%lf\n",x_low,y_low);
    for (int i=1;i<=2/h;i++){
        y_high=y_low+h*f(x_low,y_low);// for higher orders +f1(x_low,y_low)*h*h/fact(2)+f2(x_low,y_low)*h*h*h/fact(3)+f3(x_low,y_low)*h*h*h*h/fact(4)
        
        fprintf(fptr2,"%lf \t%0.12lf\n",x_low+h,y_high);//TO WRITE INTO THE FILE ORDERED PAIRS OF X AND Y AS OBTAINED FROM EXPLICIT CODE
        y_low=y_high;//RECURSIVE RELATION
        x_low=x_low+h; //RECURSIVE RELATION
        
    }
    fclose(fptr2);

    clock_t begin = clock();//START TIME
    clock_t end = clock();//END TIME
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;//TIME OF EXECUTION IN SECONDS
    printf("EXECUTION TIME: %lf",time_spent);
    return 0;
}






