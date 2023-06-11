//C PROGRAM BY AMIZHTHNI TO CALCULATE THE VARIOUS PARAMETERS ASSOCIATED WITH QUADRATIC REGRESSION AND NEWTON DIVIDED DIFFERENCE INTERPOLATION

//SPECIFYING HEADERS REQUIRED
#include <stdio.h>
#include <math.h>
#define N 3//DEFINING THE SIZE OF MATRIX

//E DENOTES SUMMATIONS IN NAMES OF FUNCTIONS

void BackwardSubstitution(double matrik[N][N+1]); //FOR PARABOLA MATRIX EQUATION SOLVING
int ForwardElimination(double matrik[N][N+1]);

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF DATAPOINTS X
double Exi(double x[6])
{
    double temp1=0;//COUNTER FOR THE SUM OF DATAPOINTS X
    for (int a=0; a<6; a++)
    {
        temp1+=x[a];
    }
    return temp1;//THE FUNCTION THUS RETURNS THE SUM OF DATAPOINTS X
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF DATAPOINTS Y
double Eyi(double y[6])
{
    double temp2=0;//COUNTER FOR THE SUM OF DATAPOINTS Y
    for (int b=0; b<6; b++)
    {
        temp2+=y[b];
    }
    return temp2;//THE FUNCTION THUS RETURNS THE SUM OF DATAPOINTS Y
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF SQUARES OF DATAPOINTS X
double Exi2(double x[6])
{
    double temp3=0;//COUNTER FOR THE SUM OF SQUARES OF DATAPOINTS X
    for (int c=0; c<6; c++)
    {
        temp3+=pow(x[c],2);
    }
    return temp3;//THE FUNCTION THUS RETURNS THE SUM OF SQUARES OF DATAPOINTS X
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF CUBES OF DATAPOINTS X
double Exi3(double x[6])
{
    double temp3=0;//COUNTER FOR THE SUM OF CUBES OF DATAPOINTS X
    for (int c=0; c<6; c++)
    {
        temp3+=pow(x[c],3);
    }
    return temp3;//THE FUNCTION THUS RETURNS THE SUM OF CUBES OF DATAPOINTS X
}
//USER DEFINED FUNCTION TO CALCULATE THE SUM OF 4THPOWERS OF DATAPOINTS X
double Exi4(double x[6])
{
    double temp3=0;//COUNTER FOR THE SUM OF 4THPOWERS OF DATAPOINTS X
    for (int c=0; c<6; c++)
    {
        temp3+=pow(x[c],4);
    }
    return temp3;//THE FUNCTION THUS RETURNS THE SUM OF 4THPOWERS OF DATAPOINTS X
}
//USER DEFINED FUNCTION TO CALCULATE THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y
double Exiyi(double x[6],double y[6])
{
    double temp4=0;//COUNTER FOR THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y
    for (int d=0; d<6; d++)
    {
        temp4+=x[d]*y[d];
    }
    return temp4;//THE FUNCTION THUS RETURNS THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y
}
double Exi2yi(double x[6],double y[6])
{
    double temp4=0;//COUNTER FOR THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y
    for (int d=0; d<6; d++)
    {
        temp4+=x[d]*x[d]*y[d];
    }
    return temp4;//THE FUNCTION THUS RETURNS THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y
}

void GaussianElimination(double matrik[N][N+1])
{
    //reduction to row reduced echelon form
    int flag = ForwardElimination(matrik);

    //if matrix is singular
    if (flag != 1)
    {
        printf("Singular Matrix.\n");

        /* if the RHS of equation corresponding to
        zero row is 0, * system has infinitely
        many solutions, else inconsistent*/
        if (matrik[flag][N])
            printf("Inconsistent System."); //from cramer's rule
        else
            printf("May have infinitely many solutions.");

        return;
    }

    /* get solution to system and print it using
    backward substitution */
    BackwardSubstitution(matrik);
}
//TO SWAP ROWS OF A MATRIX
void swap_row(double matrik[N][N+1], int i, int j)
{
    //printf("Swapped rows %d and %d\n", i, j);

    for (int k=0; k<=N; k++)
    {
        double temp = matrik[i][k];
        matrik[i][k] = matrik[j][k];
        matrik[j][k] = temp;
    }
}
// function to reduce matrix to ROW REDUCED ECHELON FORM
int ForwardElimination(double matrik[N][N+1])
{
    for (int k=0; k<N; k++)
    {
        // Initialize maximum value and index for pivot
        int i_max = k;
        int v_max = matrik[i_max][k];

        /* find greater amplitude for pivot if any */
        for (int i = k+1; i < N; i++)
            if (fabs(matrik[i][k]) > v_max)
                v_max = matrik[i][k], i_max = i;

        /* if a prinicipal diagonal element is zero,
        * it denotes that matrix is singular, and
        * will lead to a division-by-zero later. */
        if (!matrik[k][i_max])
            return k; // Matrix is singular

        /* Swap the greatest value row with current row */
        if (i_max != k)
            swap_row(matrik, k, i_max);


        for (int i=k+1; i<N; i++)
        {
            /* factor f to set current row kth element to 0,
            * and subsequently remaining kth column to 0 */
            double f = matrik[i][k]/matrik[k][k];

            /* subtract fth multiple of corresponding kth
            row element*/
            for (int j=k+1; j<=N; j++)
                matrik[i][j] -= matrik[k][j]*f;

            /* filling lower triangular matrix with zeros*/
            matrik[i][k] = 0;
        }

    }

    return 1;
}

// function to calculate the values of the unknowns
void BackwardSubstitution(double matrik[N][N+1])
{
    double yy[N]; // An array to store solution

    /* Start calculating from last equation up to the
    first */
    for (int i = N-1; i >= 0; i--)
    {
        /* start with the RHS of the equation */
        yy[i] = matrik[i][N];

        /* Initialize j to i+1 since matrix is upper
        triangular*/
        for (int j=i+1; j<N; j++)
        {
            /* subtract all the lhs values
            * except the coefficient of the variable
            * whose value is being calculated */
            yy[i] -= matrik[i][j]*yy[j];
        }

        /* divide the RHS by the coefficient of the
        unknown being calculated */
        yy[i] = yy[i]/matrik[i][i];
    }

    printf("\nSolution for the system of equations involving a0, a1 and a2:\n");
    for (int i=0; i<N; i++)
        printf("a%d\t=\t%lf\n", i,yy[i]);
    double a0=yy[0];
    double a1=yy[1];
    double a2=yy[2];
    printf("The dynamic viscosity of water μ (10−3 N·s/m2) at T=7.5 C is:");
    printf("%lf",a0+a1*7.5+a2*7.5*7.5);
}

//INTERPOLATION CODE
//WE USE  LONG INT HERE BECAUSE THE VALUES ATTAIN UTMOST PRECISION AND AT A POINT UPTO SIX DECIMALS WE GET JUST ZEROS. TO PREVENT THIS, WE USE LONG DOUBLE
long double f(int x,int  x0,int  x1,int x2,int x3,int x4,int x5)
{
    if (x==x0)      //THIS FUNCTION RETURNS THE VALUE OF F(X) FOR VALUE OF X=X0
    { return 1.787 ;
    }
    if (x==x1)      //THIS FUNCTION RETURNS THE VALUE OF F(X) FOR VALUE OF X=X1
    { return 1.519;
    }
    if (x==x2)      //THIS FUNCTION RETURNS THE VALUE OF F(X) FOR VALUE OF X=X2
    { return 1.307;
    }
    if (x==x3)      //THIS FUNCTION RETURNS THE VALUE OF F(X) FOR VALUE OF X=X3
    { return 1.002;
    }
    if (x==x4)      //THIS FUNCTION RETURNS THE VALUE OF F(X) FOR VALUE OF X=X4
    { return 0.7975;
    }
    if (x==x5)      //THIS FUNCTION RETURNS THE VALUE OF F(X) FOR VALUE OF X=X5
    { return 0.6529;
    }
    return 0;
}

long double f1(int x,int  y,int x0,int x1,int x2,int x3,int x4,int x5)//FUNCTION TO CALCULATE THE FIRST DIVIDED DIFFERENCE
{
    return (f(x,x0,x1,x2,x3,x4,x5)-f(y,x0,x1,x2,x3,x4,x5))/(x-y);
    
}
long double f2(int a,int b,int c,int x0,int x1,int x2,int x3,int x4,int x5)//FUNCTION TO CALCULATE THE SECOND DIVIDED DIFFERENCE
{
    return (f1(a,b,x0,x1,x2,x3,x4,x5)-f1(b,c,x0,x1,x2,x3,x4,x5))/(a-c);
    
}
long double f3(int a,int b,int c,int d,int x0,int x1,int x2,int x3,int x4,int x5)//FUNCTION TO CALCULATE THE THIRD DIVIDED DIFFERENCE
{
    return (f2(a,b,c,x0,x1,x2,x3,x4,x5)-f2(b,c,d,x0,x1,x2,x3,x4,x5))/(a-d);
    
}
long double f4(int a,int b,int c,int d,int e,int x0,int x1,int x2,int x3,int x4,int x5)//FUNCTION TO CALCULATE THE FOURTH DIVIDED DIFFERENCE
{
    return (f3(a,b,c,d,x0,x1,x2,x3,x4,x5)-f3(b,c,d,e,x0,x1,x2,x3,x4,x5))/(a-e);
  
}
long double f5(int a,int b,int c,int d,int e,int f,int x0,int x1,int x2,int x3,int x4,int x5)//FUNCTION TO CALCULATE THE FIFTH DIVIDED DIFFERENCE
{
    return (f4(a,b,c,d,e,x0,x1,x2,x3,x4,x5)-f4(b,c,d,e,f,x0,x1,x2,x3,x4,x5))/(40);

}
long double ffinal(double x,int x0,int x1,int x2,int x3,int x4,long double b0,long double b1, long double b2,long double b3,long double b4,long double b5)
//FUNCTION TO RETURN THE Y VALUE FOR GIVEN value of X
{
    return b0 + b1*(x-x0) +b2*(x-x0)*(x-x1) + b3*(x-x0)*(x-x1)*(x-x2)+b4*(x-x0)*(x-x1)*(x-x2)*(x-x3)+b5*(x-x0)*(x-x1)*(x-x2)*(x-x3)*(x-x4);
}
//DRIVER-MAIN PROGRAM
int main(){
    printf("GIVEN DATA POINTS ARE:\n");
    printf("\n");//newline for output formatting
    printf("x : 0     5     10    20    30     40\n");
    //PRINTING GIVEN DATASET FOR REFERENCE
    printf("y : 1.787 1.519 1.307 1.002 0.7975 0.6529\n");
    //PRINTING GIVEN DATASET FOR REFERENCE
    printf("\n");//newline for output formatting
    double x[6]={0,5,10, 20, 30, 40}, y[6]={1.787,1.519,1.307,1.002,0.7975,0.6529};
    //INITIALISING ARRAYS FOR EACH DATASET
    
    //A
    printf("(A) POLYNOMIAL INTERPOLATION\n");
    //WE HAVE 6 DATASETS, THEREFORE DEGREE OF POLYNOMIAL IS 5
    /*f5(x)=b0 + b1(x-x0) +b2(x-x0)(x-x1) + b3(x-x0)(x-x1)(x-x2)
    +b4(x-x0)(x-x1)(x-x2)(x-x3)+b5(x-x0)(x-x1)(x-x2)(x-x3)(x-x4)*/
    int  x0=x[0]; //variables for datatype x
    int  x1=x[1];
    int  x2=x[2];
    int  x3=x[3];
    int  x4=x[4];
    int  x5=x[5];
    long double b1,b2,b3,b4,b5;
   
    //b0=f(x0,x0,x1,x2,x3,x4,b0,b1,b2,b3,b4,b5);
    long double b0=f(x0,x0,x1,x2,x3,x4,x5);     //assigning coefficients for equation of f(x) polynomial interpolation
    b1=f1(x1,x0,x0,x1,x2,x3,x4,x5);             //assigning coefficients for equation of f(x) polynomial interpolation
    b2=f2(x2,x1,x0,x0,x1,x2,x3,x4,x5);          //assigning coefficients for equation of f(x) polynomial interpolation
    b3=f3(x3,x2,x1,x0,x0,x1,x2,x3,x4,x5);       //assigning coefficients for equation of f(x) polynomial interpolation
    b4=f4(x4,x3,x2,x1,x0,x0,x1,x2,x3,x4,x5);    //assigning coefficients for equation of f(x) polynomial interpolation
    b5=f5(x5,x4,x3,x2,x1,x0,x0,x1,x2,x3,x4,x5); //assigning coefficients for equation of f(x) polynomial interpolation
    
    printf("EQUATION OF POLYNOMIAL INTERPOLATION IS:\n");
    printf(" f(x)= %Le+ %Le(x-%d) +%Le(x-%d)(x-%d) + %Le(x-%d)(x-%d)(x-%d)+%Le(x-%d)(x-%d)(x-%d)(x-%d)+%Le(x-%d)(x-%d)(x-%d)(x-%d)(x-%d)\n\n",b0,b1,x0,b2,x0,x1,b3,x0,x1,x2,b4,x0,x1,x2,x3,b5,x0,x1,x2,x3,x4);//to print the  f(x) function
    printf("The dynamic viscosity of water μ (10−3 N·s/m2) at T=7.5 C is:"); //finding the interpolation value of t=7.5
    printf("%Lf",ffinal(7.5,x0,x1,x2,x3,x4,b0,b1,b2,b3,b4,b5));
    printf("\n");//newline
    printf("\n");//newline
    printf("(B) QUADRATIC REGRESSION - PARABOLA\n");
    double matrik[N][N+1] = {{6,Exi(x) ,Exi2(x), Eyi(y)},       //inititalising the matrix
                            {Exi(x),Exi2(x),Exi3(x), Exiyi(x,y)},
                            {Exi2(x),Exi3(x),Exi4(x),Exi2yi(x,y)}
                            };

    GaussianElimination(matrik); // to find the roots of matrix
    FILE *fptr; //entering values of the discrete values to a text file
    fptr = fopen("datasett.txt", "a");
    for (int u=0;u<6;u++){
        fprintf(fptr, "%lf \t %lf\n",x[u],y[u]);
    }
    fclose(fptr);
    //WE HAVE COMMENTED OUT THE ABOVE WRITING INTO FILE CODE AS WE HAVE ALREADY PERFORMED THE OPERATION
    //PLOTTING THE POINTS OF DATASETS AND THE LINE OF BEST FIT
    FILE *pipek = popen("gnuplot --persist", "w");
    fprintf(pipek, "set title \"ACTUAL DATA IN DISCRETE POINTS, QUADRATIC REGRESSION and POLYNOMIAL INTERPOLATION\"\n");//GIVING TITLE TO THE GNUPLOT
    fprintf(pipek, "set xlabel \"X VALUES \"\n");//GIVING X AXIS TITLE TO THE GNUPLOT
    fprintf(pipek, "set ylabel \"Y VALUES\"\n");//GIVING Y AXIS TITLE TO THE GNUPLOT
    fprintf(pipek, "set xrange [0:40]\n");//SETTING THE RANGE OF VALUES FOR X AXIS
    fprintf(pipek, "set yrange [0.6:1.9]\n");//SETTING THE RANGE OF VALUES FOR Y AXIS
    fprintf(pipek, "set grid\n");//TO MAKE THE GRID LINES VISIBLE FOR BETTER CLARITY
    fprintf(pipek, "plot 'datasett.txt' with points pointtype 7 lt rgb  \"red\" title 'DISCRETE DATA POINTS',1.787000 - 5.360000 * 10 ** (-02) * (x-0) +(1.120000 * 10 **  (-03)) * (x-0) * (x-5) -1.633333 * 10 ** (-5) * (x-0) * (x-5) * (x-10)+1.566667 * 10 ** (-7) * (x-0) * ( x-5) * ( x-10) * (x-20)-4.404762 * 10 ** (-10) * (x-0) * (x-5) * (x-10) * (x-20) * (x-30) lt rgb \"blue\" title 'POLYNOMIAL INTERPOLATION',1.767245-0.049493*x+0.000548*x**2 lt rgb \"black\" title 'PARABOLA - QUADRATIC  REGRESSION' \n");
    /*PLOTTING THE TWO REQUIRED GRAPHS: 1. THE quadratic REGRESSION  POTTED IN BLack
                                        2. THE DISCRETE DATAPOINTS OF THE GIVEN DATA IN RED
                                        3. INTERPOLATION PLOT LINES IN BLUE
     THE STATEMENT POINTTYPE 7 GIVES US THE UNIQUE FILLED CIRCLE WHICH IS PREFERRED OVER THE USUAL CROSS
     LT RGB COLOUR GIVES US THE DESIRED COLOUR FOR THE POINT
    */
    
}




