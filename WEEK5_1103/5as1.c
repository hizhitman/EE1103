//C PROGRAM BY AMIZHTHNI TO CALCULATE THE VARIOUS PARAMETERS ASSOCIATED WITH LINEAR REGRESSION

//SPECIFYING HEADERS REQUIRED
#include <stdio.h>
#include <math.h>

//E DENOTES SUMMATIONS IN NAMES OF FUNCTIONS
//USER DEFINED FUNCTION TO CALCULATE THE SUM OF DATAPOINTS X

double Exi(double x[11])
{
    int temp1=0;//COUNTER FOR THE SUM OF DATAPOINTS X
    for (int a=0; a<11; a++)
    {
        temp1+=x[a];
    }
    return temp1;//THE FUNCTION THUS RETURNS THE SUM OF DATAPOINTS X
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF DATAPOINTS Y
double Eyi(double y[11])
{
    int temp2=0;//COUNTER FOR THE SUM OF DATAPOINTS Y
    for (int b=0; b<11; b++)
    {
        temp2+=y[b];
    }
    return temp2;//THE FUNCTION THUS RETURNS THE SUM OF DATAPOINTS Y
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF SQUARES OF DATAPOINTS X
double Exi2(double x[11])
{
    int temp3=0;//COUNTER FOR THE SUM OF SQUARES OF DATAPOINTS X
    for (int c=0; c<11; c++)
    {
        temp3+=pow(x[c],2);
    }
    return temp3;//THE FUNCTION THUS RETURNS THE SUM OF SQUARES OF DATAPOINTS X
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y
double Exiyi(double x[11],double y[11])
{
    int temp4=0;//COUNTER FOR THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y
    for (int d=0; d<11; d++)
    {
        temp4+=x[d]*y[d];
    }
    return temp4;//THE FUNCTION THUS RETURNS THE SUM OF PRODUCT OF RESPECTIVE DATAPOINTS IN X AND Y
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF SQUARES OF RESIDUALS
double Sr(double x[11],double y[11],double intercept,double slope)
{
    double temp5=0;//COUNTER TO CALCULATE THE SUM OF SQUARES OF RESIDUALS
    for (int e=0; e<11; e++)
    {
        temp5+=pow(y[e]-intercept-slope*x[e],2);
    }
    return temp5;//THE FUNCTION THUS RETURNS THE SUM OF SQUARES OF RESIDUALS
}

//USER DEFINED FUNCTION TO CALCULATE THE SUM OF SQUARES AROUND MEAN
double St(double x[11],double y[11])
{
    double temp6=0;//COUNTER TO CALCULATE THE SUM OF SQUARES AROUND MEAN
    for (int f=0; f<11; f++)
    {
        temp6+=pow(y[f]-Eyi(y)/11,2);
    }
    return temp6;//THE FUNCTION THUS RETURNS THE SUM OF SQUARES AROUND MEAN
}

//DRIVER-MAIN PROGRAM
int main(){
    printf("GIVEN DATA POINTS ARE:\n");
    printf("\n");//newline for output formatting
    printf("x : 6  7  11 15 17 21 23 29 29 37 39\n");
    //PRINTING GIVEN DATASET FOR REFERENCE
    printf("y : 29 21 29 14 21 15 7  7  13  0  3\n");
    //PRINTING GIVEN DATASET FOR REFERENCE
    printf("\n");//newline for output formatting
    double x[11]={6,7,11,15,17,21,23,29,29,37,39}, y[11]={29,21,29,14,21,15,7,7,13,0,3};
    //INITIALISING ARRAYS FOR EACH DATASET
    printf("SLOPE OF LINEAR REGRESSION LINE(a1):");
    double slope=(11*Exiyi(x,y)-(Exi(x)*Eyi(y)))/(11*Exi2(x)-Exi(x)*Exi(x));
    //CALCULATING THE SLOPE DIRECTLY USING FORMULA INSTEAD OF GAUSS ELIMINATION
    //a1=[n(Exiyi)-(Exi)(Eyi)]/[n(Exi^2)-(Exi)^2]
    printf("%lf\n",slope);
    printf("INTERCEPT OF LINEAR REGRESSION LINE(a0):");
    //CALCULATING THE INTERCEPT DIRECTLY USING FORMULA INSTEAD OF GAUSS ELIMINATION
    //a0=1/n(Eyi-a1(Exi))
    double intercept=(Eyi(y)/11-slope*Exi(x)/11);
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
    //r2=1-Sr/St
    //where  r is the coefficient of correlation
    //where r^2 is the coefficient if  determination
    double r=sqrt(1-(sr/st));
    printf("%lf\n",r);
    printf("STANDARD DEVIATION OF REGRESSION(Sxy):");
    //Sxy=(Sr/(n-2))^0.5
    double sxy=sqrt(sr/9);
    printf("%lf\n",sxy);
    //note: it is highly important to find variables in the  given order: a1 -> a0 -> Sr -> St -> r^2
    //WRITING DATA FOR PLOTTING
    /*
    FILE *fptr;
    fptr = fopen("datasets.txt", "a");
    for (int u=0;u<11;u++){
        fprintf(fptr, "%lf \t %lf\n",x[u],y[u]);
    }
    fclose(fptr);*/
    //WE HAVE COMMENTED OUT THE ABOVE WRITING INTO FILE CODE AS WE HAVE ALREADY PERFORMED THE OPERATION
    //PLOTTING THE POINTS OF DATASETS AND THE LINE OF BEST FIT
    FILE *pipek = popen("gnuplot --persist", "w");
    fprintf(pipek, "set title \"LINEAR REGRESSION FOR GIVEN DATASET\"\n");//GIVING TITLE TO THE GNUPLOT
    fprintf(pipek, "set xlabel \"X VALUES \"\n");//GIVING X AXIS TITLE TO THE GNUPLOT
    fprintf(pipek, "set ylabel \"Y VALUES\"\n");//GIVING Y AXIS TITLE TO THE GNUPLOT
    fprintf(pipek, "set xrange [0:40]\n");//SETTING THE RANGE OF VALUES FOR X AXIS
    fprintf(pipek, "set yrange [-5:35]\n");//SETTING THE RANGE OF VALUES FOR Y AXIS
    fprintf(pipek, "set grid\n");//TO MAKE THE GRID LINES VISIBLE FOR BETTER CLARITY
    fprintf(pipek, "plot 'datasets.txt' with points pointtype 7 lt rgb  \"red\",-0.780547*x+31.058898 lt rgb \"blue\" \n");
    /*PLOTTING THE TWO REQUIRED GRAPHS: 1. THE LINEAR REGRESSION LINE POTTED IN BLUE (y=-0.780547*x+31.058898)
                                        2. THE DISCRETE DATAPOINTS OF THE GIVEN DATA IN RED
     THE STATEMENT POINTTYPE 7 GIVES US THE UNIQUE FILLED CIRCLE WHICH IS PREFERRED OVER THE USUAL CROSS
     LT RGB COLOUR GIVES US THE DESIRED COLOUR FOR THE POINT
    */
    
    return 0;
}

