//To find the inverse of a matrix using LU decomposition
#include<math.h>
#include<stdio.h>

int main ()
{
  //Variable declarations
  int i, j, n, m;
  float M[3][3], d[3], C[3][3];
  //declaring two matrices, one 3X3 matrix for taking Matrix M and another column vector :M=MATRIX M of question,d= is the matrix D
  float x, s[3][3], y[3];
  
  //declaring two matrices, one 3X3 matrix for taking Matrix M and another column vector.s takes the inverse numbers, y is the intermediate d matrix in inverse calculation
  //Function def type
  void LU();

  n = 2;
  /* the matrix to be inverted is initialised here*/
  M[0][0] = 1.0; //initialised from 0 till N-1 since we use pointers
  M[0][1] = 4.0;
  M[0][2] = -3.0;
  M[1][0] = -2.0;
  M[1][1] = 8.0;
  M[1][2] = 5.0;
  M[2][0] = 3.0;
  M[2][1] = 4.0;
  M[2][2] = 7.0;


/* Call a sub-function to calculate the LU decomposed matrix. Note that
we pass the two dimensional array [D] to the function and get it back */
  LU (M, n);

  printf (" \n");
  printf
    ("The matrix LU decomposed which stores a_ij and f_ij values together i.e. L and U together \n");
  for (m = 0; m <= 2; m++)
    {
      printf (" %f  %f   %f \n", M[m][0], M[m][1], M[m][2]);
      //we reduce another loop by typing the column values
    }


  printf (" \n");
//Printing the U MATRIX from the LU matrix. we selectively choose the required elements. however this holds good if and only if N=3.
  printf
    ("The matrix U which is the upper triangular matrix which satisfies UX=D AND LU=A \n");
  printf (" %f  %f   %f \n", M[0][0], M[0][1], M[0][2]);
  printf (" %f  %f   %f \n", 0.00, M[1][1], M[1][2]);
  printf (" %f  %f   %f \n", 0.00, 0.00, M[2][2]);

  printf (" \n");
//Printing the L MATRIX from the LU matrix. we selectively choose the required elements. however this holds good if and only if N=3.
  printf
    ("The matrix L which is the lower triangular matrix which satisfies LD=B AND LU=A \n");
  printf (" %f  %f   %f \n", 1.0, 0.0, 0.0);
  printf (" %f  %f   %f \n", M[1][0], 1.0, 0.0);
  printf (" %f  %f   %f \n", M[2][0], M[2][1], 1.00);
printf (" \n");
//TO FIND THE INVERSE

/* to find the inverse we solve [D][y]=[d] with only one element in
the [d] array put equal to one at a time, for example:
in the first iteration: the first row of the column vector is 1 and rest are zeros. Similarly for the i th iteration,
the i th row of the column vector iz one while others are zeros*/

  for (m = 0; m <= 2; m++)
    {
      d[0] = 0.0; //all the non m th elements of the column vector are zeros
      d[1] = 0.0; //all the non m th elements of the column vector are zeros
      d[2] = 0.0; //all the non m th elements of the column vector are zeros
      d[m] = 1.0; //the m th iteration takes value 1
      for (i = 0; i <= n; i++)//row iteration
    {
      x = 0.0;
      for (j = 0; j <= i - 1; j++)//column iteration
        {
          x = x + M[i][j] * y[j];

        }
      y[i] = (d[i] - x);//intermediate d value
    }

      for (i = n; i >= 0; i--)
    {
      x = 0.0;
      for (j = i + 1; j <= n; j++)
        {
          x = x + M[i][j] * s[j][m];
        }
      s[i][m] = (y[i] - x) / M[i][i]; //whichever row is 1, the corresponding column gets the values. while rest are zeros.
    }
    }

//Print the inverse matrix
  printf ("The Inverse Matrix\n");
  for (m = 0; m <= 2; m++)
    {
      printf (" %f %f %f \n", s[m][0], s[m][1], s[m][2]);
    }
}

//The function that calcualtes the LU decomposed matrix.
//we use the concept of pointers here
void LU (float (*M)[3][3], int n) //we use pointers here so that the changes here are reflected in Original M matrix too and no need to return any value
{ //this function gives us the lu matrix with both upper and lower combined i.e. all aijs and fijs.
  int i, j, k, m;
  float multi;
  for (k = 0; k <= n - 1; k++) //fixes the column
    {
      for (j = k + 1; j <= n; j++)//fixes the row
    {
       multi= (*M)[j][k] / (*M)[k][k]; //the factor by which we need to multiply
      for (i = k; i <= n; i++)
        {
          (*M)[j][i] = (*M)[j][i] - multi * (*M)[k][i]; //this operation subtracts every element of the j th row from multi times corresponding element of preceding row
        }
      (*M)[j][k] = multi; //this is essentially because after all subtractions the leading term must remain intact and must not vanish
    }
    }
}
