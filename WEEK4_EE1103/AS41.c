//GaussianElimination
#include <stdio.h>
#include <math.h>

#define N 3     // Number of unknowns

// function to reduce matrix to row reduced echelon form.
int ForwardElimination(double matrik[N][N+1]);

// function to calculate the values of the unknowns
void BackwardSubstitution(double matrik[N][N+1]);

// function to get matrix content
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

// function for elementary operation of swapping two rows
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

// function to print matrix content at any stage
void print(double matrik[N][N+1])
{
    for (int i=0; i<N; i++, printf("\n"))
        for (int j=0; j<=N; j++)
            printf("%lf ", matrik[i][j]);

    printf("\n");
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
void BackwardSubstitution(double matrik[N][N+1])
{
    double x[N]; // An array to store solution

    /* Start calculating from last equation up to the
    first */
    for (int i = N-1; i >= 0; i--)
    {
        /* start with the RHS of the equation */
        x[i] = matrik[i][N];

        /* Initialize j to i+1 since matrix is upper
        triangular*/
        for (int j=i+1; j<N; j++)
        {
            /* subtract all the lhs values
            * except the coefficient of the variable
            * whose value is being calculated */
            x[i] -= matrik[i][j]*x[j];
        }

        /* divide the RHS by the coefficient of the
        unknown being calculated */
        x[i] = x[i]/matrik[i][i];
    }

    printf("\nSolution for the system of equations involving V1 ,V2 and V3:\n");
    for (int i=0; i<N; i++)
        printf("v%d\t=\t%lf\n", i+1,x[i]);
}

// Driver program
int main()
{
    /* input matrix */
    double matrik[N][N+1] = {{0.5833, -0.3333,-0.25, -11}, //these values were obtained by solving the kirchhoff's equations manually
                        {-0.3333, 1.4762, -0.1429, 3},
                        {-0.25, -0.1429, 0.5929, 25}
                        };

    GaussianElimination(matrik);
    //print(matrik); (can be used to print as and when required)
    printf("where v1, v2 and v3 are the potentials of junctions 1, 2 and 3 respectively");
    return 0;
}



