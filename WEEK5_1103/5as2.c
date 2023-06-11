//C PROGRAM BY AMIZHTHNI TO CALCULATE THE VARIOUS PARAMETERS ASSOCIATED WITH LINEAR REGRESSION

//SPECIFYING HEADERS REQUIRED
#include <stdio.h>
#include <math.h>
#define N 3

int ForwardElimination(double matrik[N][N+1]);
void BackwardSubstitution(double matrik[N][N+1],double x[10],double y[10]);

//E DENOTES SUMMATIONS IN NAMES OF FUNCTIONS

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF DATAPOINTS X

double Exi(double x[10])
{
    int temp1=0;//COUNTER FOR THE SUM OF DATAPOINTS X
    for (int a=0; a<10; a++)
    {
        temp1+=x[a];
    }
    return temp1;//THE FUNCTION THUS RETURNS THE SUM OF DATAPOINTS X
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF DATAPOINTS Y
double Eyi(double y[10])
{
    int temp2=0;//COUNTER FOR THE SUM OF DATAPOINTS Y
    for (int b=0; b<10; b++)
    {
        temp2+=y[b];
    }
    return temp2;//THE FUNCTION THUS RETURNS THE SUM OF DATAPOINTS Y
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF SQUARES OF DATAPOINTS X
double Exi2(double x[10])
{
    int temp3=0;//COUNTER FOR THE SUM OF SQUARES OF DATAPOINTS X
    for (int c=0; c<10; c++)
    {
        temp3+=pow(x[c],2);
    }
    return temp3;//THE FUNCTION THUS RETURNS THE SUM OF SQUARES OF DATAPOINTS X
}
//USER DEFINED FUNCTION TO CALCULATE THE SUM OF CUBES OF DATAPOINTS X
double Exi3(double x[10])
{
    int temp3=0;//COUNTER FOR THE SUM OF CUBES OF DATAPOINTS X
    for (int c=0; c<10; c++)
    {
        temp3+=pow(x[c],3);
    }
    return temp3;//THE FUNCTION THUS RETURNS THE SUM OF CUBES OF DATAPOINTS X
}
//USER DEFINED FUNCTION TO CALCULATE THE SUM OF 4TH POWERS OF DATAPOINTS X
double Exi4(double x[10])
{
    int temp3=0;//COUNTER FOR THE SUM OF 4TH POWERS OF DATAPOINTS X
    for (int c=0; c<10; c++)
    {
        temp3+=pow(x[c],4);
    }
    return temp3;//THE FUNCTION THUS RETURNS THE SUM OF 4TH POWERS OF DATAPOINTS X
}


//USER DEFINED FUNCTION TO CALCULATE THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y
double Exiyi(double x[10],double y[10])
{
    int temp4=0;//COUNTER FOR THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y
    for (int d=0; d<10; d++)
    {
        temp4+=x[d]*y[d];
    }
    return temp4;//THE FUNCTION THUS RETURNS THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X^2 AND Y
double Exi2yi(double x[10],double y[10])
{
    int temp4=0;//COUNTER FOR THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X ^2 AND Y
    for (int d=0; d<10; d++)
    {
        temp4+=x[d]*x[d]*y[d];
    }
    return temp4;//THE FUNCTION THUS RETURNS THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X^2 AND Y
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF SQUARES OF RESIDUALS FOR PARABOLA
double bSr(double x[10],double y[10],double a0,double a1,double a2)
{
    double temp5=0;//COUNTER TO CALCULATE THE SUM OF SQUARES OF RESIDUALS
    for (int e=0; e<10; e++)
    {
        temp5+=pow((y[e]-a0-a1*x[e]-a2*(pow(x[e],2))),2);//THIS IS A DIFFERENT EQUATION AS WE USE QUADRATIC REGRESSION AS OPPOSED TO NORMAL LINEAR REGRESSION
    }
    return temp5;//THE FUNCTION THUS RETURNS THE SUM OF SQUARES OF RESIDUALS
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF SQUARES OF RESIDUALS
double Sr(double x[10],double y[10],double intercept,double slope)
{
    double temp5=0;//COUNTER TO CALCULATE THE SUM OF SQUARES OF RESIDUALS
    for (int e=0; e<10; e++)
    {
        temp5+=pow(y[e]-intercept-slope*x[e],2);
    }
    return temp5;//THE FUNCTION THUS RETURNS THE SUM OF SQUARES OF RESIDUALS
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF SQUARES AROUND MEAN
double St(double x[10],double y[10])
{
    double temp6=0;//COUNTER TO CALCULATE THE SUM OF SQUARES AROUND MEAN
    for (int f=0; f<10; f++)
    {
        temp6+=pow(y[f]-Eyi(y)/10,2);//NOTE THAT THIS IS THE SAME FOR ALL TYPES OF REGRESSION
    }
    return temp6;//THE FUNCTION THUS RETURNS THE SUM OF SQUARES AROUND MEAN
}

//FOR SOLVING PARABOLA USING MATRIX GAUSS ELIMINATION
void GaussianElimination(double matrik[N][N+1],double x[10],double y[10])
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
    BackwardSubstitution(matrik,x,y);
}
//FUNCTION TO SWAP ROWS. THIS IS LATER USED IN BACK  WARD SUBSTITUTION
void swap_row(double matrik[N][N+1], int i, int j)
{
    //printf("Swapped rows %d and %d\n", i, j);

    for (int k=0; k<=N; k++)
    {
        double temp = matrik[i][k];//TEMPORARY  VARIABLE THAT HELPS  IN TRANSFERRING REQUIRED VALUE
        matrik[i][k] = matrik[j][k];//SWAP
        matrik[j][k] = temp;
    }
}
// function to reduce matrix to r.e.f.
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
void BackwardSubstitution(double matrik[N][N+1],double x[10],double y[10])
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

    printf("\nSolution for the system of equations involving a0,a1 and a2 for y=a0+a1x+a2x^2:\n");
    for (int i=0; i<N; i++)
        printf("a%d\t=\t%0.10lf\n", i,yy[i]);
    double a0=yy[0];
    double a1=yy[1];
    double a2=yy[2];
    double bsr=bSr(x,y,a0,a1,a2);
    printf("coefficient of correlation r:%lf\n",sqrt(1-bsr/601.6));
}


//USER DEFINED FUNCTION TO CALCULATE THE SUM OF DATAPOINTS X FOR POWER EQUATION

double PExi(double x[10])
{
    double temp1=0;//COUNTER FOR THE SUM OF DATAPOINTS X FOR POWER EQUATION
    for (int a=0; a<10; a++)
    {
        temp1+=log(x[a]);
    }
    return temp1;//THE FUNCTION THUS RETURNS THE SUM OF DATAPOINTS X FOR POWER EQUATION
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF DATAPOINTS Y FOR POWER EQUATION
double PEyi(double y[10])
{
    double temp2=0;//COUNTER FOR THE SUM OF DATAPOINTS Y FOR POWER EQUATION
    for (int b=0; b<10; b++)
    {
        temp2+=log(y[b]);
    }
    return temp2;//THE FUNCTION THUS RETURNS THE SUM OF DATAPOINTS Y FOR POWER EQUATION
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF SQUARES OF DATAPOINTS X FOR POWER EQUATION
double PExi2(double x[10])
{
    double temp3=0;//COUNTER FOR THE SUM OF SQUARES OF DATAPOINTS X FOR POWER EQUATION
    for (int c=0; c<10; c++)
    {
        temp3+=pow(log(x[c]),2);
    }
    return temp3;//THE FUNCTION THUS RETURNS THE SUM OF SQUARES OF DATAPOINTS X FOR POWER EQUATION
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y FOR POWER EQUATION
double PExiyi(double x[10],double y[10])
{
    double temp4=0;//COUNTER FOR THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y FOR POWER EQUATION
    for (int d=0; d<10; d++)
    {
        temp4+=log(x[d])*log(y[d]);
    }
    return temp4;//THE FUNCTION THUS RETURNS THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y FOR POWER EQUATION
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF SQUARES OF RESIDUALS FOR POWER EQUATION
double PSr(double x[10],double y[10],double Pintercept,double Pslope)
{
    double temp5=0;//COUNTER TO CALCULATE THE SUM OF SQUARES OF RESIDUALS FOR POWER EQUATION
    for (int e=0; e<10; e++)
    {
        temp5+=pow(y[e]-pow(M_E,Pintercept)*(pow(x[e],Pslope)),2);
    }
    return temp5;//THE FUNCTION THUS RETURNS THE SUM OF SQUARES OF RESIDUALS FOR POWER EQUATION
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF SQUARES AROUND MEAN  FOR POWER EQUATION
double PSt(double x[10],double y[10])
{
    double temp6=0;//COUNTER TO CALCULATE THE SUM OF SQUARES AROUND MEAN  FOR POWER EQUATION
    for (int f=0; f<10; f++)
    {
        temp6+=pow(y[f]-Eyi(y)/10,2);
    }
    return temp6;//THE FUNCTION THUS RETURNS THE SUM OF SQUARES AROUND MEAN FOR POWER EQUATION
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF DATAPOINTS X FOR SATURATION EQUATION
double SExi(double x[10])
{
    double temp1=0;//COUNTER FOR THE SUM OF DATAPOINTS X FOR SATURATION EQUATION
    for (int a=0; a<10; a++)
    {
        temp1+=pow(x[a],-1);
    }
    return temp1;//THE FUNCTION THUS RETURNS THE SUM OF DATAPOINTS X FOR SATURATION EQUATION
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF DATAPOINTS Y FOR SATURATION EQUATION
double SEyi(double y[10])
{
    double temp2=0;//COUNTER FOR THE SUM OF DATAPOINTS Y FOR SATURATION EQUATION
    for (int b=0; b<10; b++)
    {
        temp2+=pow(y[b],-1);
    }
    return temp2;//THE FUNCTION THUS RETURNS THE SUM OF DATAPOINTS Y FOR SATURATION EQUATION
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF SQUARES OF DATAPOINTS X FOR SATURATION EQUATION
double SExi2(double x[10])
{
    double temp3=0;//COUNTER FOR THE SUM OF SQUARES OF DATAPOINTS X FOR SATURATION EQUATION
    for (int c=0; c<10; c++)
    {
        temp3+=pow(pow(x[c],-1),2);
    }
    return temp3;//THE FUNCTION THUS RETURNS THE SUM OF SQUARES OF DATAPOINTS X FOR SATURATION EQUATION
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y FOR SATURATION EQUATION
double SExiyi(double x[10],double y[10])
{
    double temp4=0;//COUNTER FOR THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y FOR SATURATION EQUATION
    for (int d=0; d<10; d++)
    {
        temp4+=pow(x[d],-1)*pow(y[d],-1);
    }
    return temp4;//THE FUNCTION THUS RETURNS THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y FOR SATURATION EQUATION
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF SQUARES OF RESIDUALS FOR SATURATION EQUATION
double SSr(double x[10],double y[10],double Sintercept,double Sslope)
{
    double temp5=0;//COUNTER TO CALCULATE THE SUM OF SQUARES OF RESIDUALS FOR SATURATION EQUATION
    for (int e=0; e<10; e++)
    {
        temp5+=pow(y[e]-pow(Sintercept,-1)*x[e]/(Sslope/Sintercept+x[e]),2);
    }
    return temp5;//THE FUNCTION THUS RETURNS THE SUM OF SQUARES OF RESIDUALS FOR SATURATION EQUATION
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF SQUARES AROUND MEAN FOR SATURATION EQUATION
double SSt(double x[10],double y[10])
{
    double temp6=0;//COUNTER TO CALCULATE THE SUM OF SQUARES AROUND MEAN FOR SATURATION EQUATION
    for (int f=0; f<10; f++)
    {
        temp6+=pow(y[f]-Eyi(y)/10,2);
    }
    return temp6;//THE FUNCTION THUS RETURNS THE SUM OF SQUARES AROUND MEAN FOR SATURATION EQUATION
}
//DRIVER-MAIN PROGRAM
int main(){
    printf("GIVEN DATA POINTS ARE:\n");
    printf("\n");//newline for output formatting
    printf("x : 5  10 15 20 25 30 35 40 45 50\n");
    //PRINTING GIVEN DATASET FOR REFERENCE
    printf("y : 17 24 31 33 37 37 40 40 42 41\n");
    //PRINTING GIVEN DATASET FOR REFERENCE
    printf("\n");//newline for output formatting
    double x[10]={5,10,15,20,25,30,35,40,45,50}, y[10]={17,24,31,33,37,37,40,40,42,41};
    //INITIALISING ARRAYS FOR EACH DATASET
    
    //A
    printf("(A) LINEAR REGRESSION - STRAIGHT LINE\n");
    printf("\n");//newline for output formatting
    printf("\n");//newline for output formatting
    printf("SLOPE OF LINEAR REGRESSION LINE(a1):");
    double slope=(10*Exiyi(x,y)-(Exi(x)*Eyi(y)))/(10*Exi2(x)-Exi(x)*Exi(x));
    //CALCULATING THE SLOPE DIRECTLY USING FORMULA INSTEAD OF GAUSS ELIMINATION
    //a1=[n(Exiyi)-(Exi)(Eyi)]/[n(Exi^2)-(Exi)^2]
    printf("%lf\n",slope);
    printf("INTERCEPT OF LINEAR REGRESSION LINE(a0):");
    //CALCULATING THE INTERCEPT DIRECTLY USING FORMULA INSTEAD OF GAUSS ELIMINATION
    //a0=1/n(Eyi-a1(Exi))
    double intercept=(Eyi(y)/10-slope*Exi(x)/10);
    printf("%lf\n",intercept);
    printf("\n");//newline for output formatting
    printf("SQUARES OF RESIDUALS(Sr):");
    //WE OBTAIN THE PREVIOUS TWO RESULTS BY MINIMISING THIS PARTICULAR EXPRESSION
    //SR = (Y-a0-a1X)^2
    double sr=Sr(x,y,intercept,slope);
    printf("%lf\n",sr);
    printf("SQUARES AROUND MEAN(St):");
    //St=E{(yi-(Eyi/n))^2}
    double st=St(x,y);
    printf("%lf\n",st);
    printf("\n");//newline for output formatting
    printf("COEFFICIENT OF CORRELATION(R):");
    //r=(1-Sr/St)^0.5
    //where  r is the coefficient of determination
    //where r^2 is the coefficient if correlation
    double r=sqrt(1-(sr/st));
    printf("%lf\n",r);
    printf("STANDARD DEVIATION OF REGRESSION(Sxy):");
    //Sxy=(Sr/(n-2))^0.5
    double sxy=sqrt(sr/8);
    printf("%lf\n",sxy);
    printf("\n");//newline for output formatting
    printf("\n");//newline for output formatting
    
    //B
    printf("(B) LEAST SQUARES REGRESSION - A POWER EQUATION\n");
    printf("y=a*x^b, let this be the power equation \n");
    printf("log y=log a+ b*log x, let this be the power equation expressed in linear format \n");
    printf("\n");//newline for output formatting
    printf("\n");//newline for output formatting
    printf("SLOPE OF LINEAR REGRESSION LINE(b):");
    double Pslope=(10*PExiyi(x,y)-(PExi(x)*PEyi(y)))/(10*PExi2(x)-PExi(x)*PExi(x));
    //CALCULATING THE SLOPE DIRECTLY USING FORMULA INSTEAD OF GAUSS ELIMINATION
    //a1=[n(Exiyi)-(Exi)(Eyi)]/[n(Exi^2)-(Exi)^2]
    printf("%lf\n",Pslope);//b
    printf("E^INTERCEPT OF LINEAR REGRESSION LINE(a):");
    //CALCULATING THE INTERCEPT DIRECTLY USING FORMULA INSTEAD OF GAUSS ELIMINATION
    //a0=1/n(Eyi-a1(Exi))
    double Pintercept=(PEyi(y)/10-Pslope*PExi(x)/10);
    printf("%lf\n",pow(M_E,Pintercept));//a
    printf("\n");//newline for output formatting
    printf("SQUARES OF RESIDUALS(Sr):");
    //WE OBTAIN THE PREVIOUS TWO RESULTS BY MINIMISING THIS PARTICULAR EXPRESSION
    //SR = (Y-a0-a1X)^2
    double Psr=PSr(x,y,Pintercept,Pslope);
    printf("%lf\n",Psr);
    printf("SQUARES AROUND MEAN(St):");
    //St=E{(yi-(Eyi/n))^2}
    double Pst=PSt(x,y);
    printf("%lf\n",Pst);
    printf("\n");//newline for output formatting
    printf("COEFFICIENT OF CORRELATION(R):");
    //r2=1-Sr/St
    //where  r is the coefficient of correlation
    //where r^2 is the coefficient of determination
    double Pr=sqrt(1-(Psr/Pst));
    printf("%lf\n",Pr);
    printf("STANDARD DEVIATION OF REGRESSION(Sxy):");
    //Sxy=(Sr/(n-2))^0.5
    double Psxy=sqrt(Psr/8);
    printf("%lf\n",Psxy);
    printf("\n");//newline for output formatting
    printf("\n");//newline for output formatting
    //C
    printf("(C) LEAST SQUARES REGRESSION - SATURATION GROWTH RATE\n");
    printf(" y=ax/(b+x), let this be the saturation equation \n");
    printf(" 1/y = 1/a +(b/a)(1/x)let this be the saturation equation expressed in linear format \n");
    printf("\n");//newline for output formatting
    printf("\n");//newline for output formatting
    printf("SLOPE OF LINEAR REGRESSION LINE(b/a):");
    double Sslope=(10*SExiyi(x,y)-(SExi(x)*SEyi(y)))/(10*SExi2(x)-SExi(x)*SExi(x));
    //CALCULATING THE SLOPE DIRECTLY USING FORMULA INSTEAD OF GAUSS ELIMINATION
    //a1=[n(Exiyi)-(Exi)(Eyi)]/[n(Exi^2)-(Exi)^2]
    printf("%lf\n",Sslope);
    printf("INTERCEPT OF LINEAR REGRESSION LINE(1/a):");
    //CALCULATING THE INTERCEPT DIRECTLY USING FORMULA INSTEAD OF GAUSS ELIMINATION
    //a0=1/n(Eyi-a1(Exi))
    double Sintercept=(SEyi(y)/10-Sslope*SExi(x)/10);
    printf("%lf\n",Sintercept);
    printf("a:%lf\n",pow(Sintercept,-1));
    printf("b:%lf\n",Sslope/Sintercept);
    printf("\n");//newline for output formatting
    printf("SQUARES OF RESIDUALS(Sr):");
    //WE OBTAIN THE PREVIOUS TWO RESULTS BY MINIMISING THIS PARTICULAR EXPRESSION
    //SR = (Y-a0-a1X)^2
    double Ssr=SSr(x,y,Sintercept,Sslope);
    printf("%lf\n",Ssr);
    printf("SQUARES AROUND MEAN(St):");
    //St=E{(yi-(Eyi/n))^2}
    double Sst=SSt(x,y);
    printf("%lf\n",Sst);
    printf("\n");//newline for output formatting
    printf("COEFFICIENT OF CORRELATION(R):");
    //r2=1-Sr/St
    //where  r is the coefficient of correlation
    //where r^2 is the coefficient of determination
    double Sr=sqrt(1-(Ssr/Sst));
    printf("%lf\n",Sr);
    printf("STANDARD DEVIATION OF REGRESSION(Sxy):");
    //Sxy=(Sr/(n-2))^0.5
    double Ssxy=sqrt(Ssr/8);
    printf("%lf\n",Ssxy);
    printf("\n");//newline for output formatting
    printf("\n");//newline for output formatting
    
    //note: it is highly important to find variables in the  given order: a1 -> a0 -> Sr -> St -> r^2
    
    //A
    printf("(D) QUADRATIC REGRESSION - PARABOLA\n");
    double matrik[N][N+1] = {{10,Exi(x) ,Exi2(x), Eyi(y)},
                        {Exi(x),Exi2(x),Exi3(x), Exiyi(x,y)},
                        {Exi2(x),Exi3(x),Exi4(x),Exi2yi(x,y)}
                        };
    
    GaussianElimination(matrik,x,y);
    //WRITING DATA FOR PLOTTING
    
    /*FILE *fptr;
    fptr = fopen("dataseta.txt", "a");
    for (int u=0;u<10;u++){
        fprintf(fptr, "%lf \t %lf\n",x[u],y[u]);
    }
    fclose(fptr);*/
    //WE HAVE COMMENTED OUT THE ABOVE WRITING INTO FILE CODE AS WE HAVE ALREADY PERFORMED THE OPERATION OF WRITING INTO THE FILE
    //PLOTTING THE POINTS OF DATASETS AND THE LINE OF BEST FIT
    FILE *pipek = popen("gnuplot --persist", "w");
    fprintf(pipek, "set title \"ACTUAL DATA IN DISCRETE POINTS AND LINEAR AND QUADRATIC REGRESSION\"\n");//GIVING TITLE TO THE GNUPLOT
    fprintf(pipek, "set xlabel \"X VALUES \"\n");//GIVING X AXIS TITLE TO THE GNUPLOT
    fprintf(pipek, "set ylabel \"Y VALUES\"\n");//GIVING Y AXIS TITLE TO THE GNUPLOT
    fprintf(pipek, "set xrange [0:50]\n");//SETTING THE RANGE OF VALUES FOR X AXIS
    fprintf(pipek, "set yrange [0:50]\n");//SETTING THE RANGE OF VALUES FOR Y AXIS
    fprintf(pipek, "set grid\n");//TO MAKE THE GRID LINES VISIBLE FOR BETTER CLARITY
    fprintf(pipek, "plot 'dataseta.txt' with points pointtype 7 lt rgb  \"red\" title 'DISCRETE DATA POINTS',0.494545*x+20.6 lt rgb \"blue\" title 'STRAIGHT LINE FIT',9.952915*x**0.385077 lt rgb \"violet\" title 'POWER EQUATION FIT',50.092118*x/(x+9.891369)  lt rgb \"green\"  title 'SATURATION EQUATION',11.766667+1.3778787879*x-0.0160606061*x**2 lt rgb \"black\" title 'PARABOLA ' \n");
    /*PLOTTING THE  REQUIRED GRAPHS: 1. THE LINEAR REGRESSION LINE POTTED IN BLUE
                                        2. THE DISCRETE DATAPOINTS OF THE GIVEN DATA IN RED
                                        3. THE  POWER EQUATION IN VIOLET
                                        4. SATURATION EQUATION IN GREEN
                                        5. PARABOLA QUADRATIC REGRESSION IN BLACK
     THE STATEMENT POINTTYPE 7 GIVES US THE UNIQUE FILLED CIRCLE WHICH IS PREFERRED OVER THE USUAL CROSS
     LT RGB COLOUR GIVES US THE DESIRED COLOUR FOR THE POINT
    */
    
    return 0;
}





