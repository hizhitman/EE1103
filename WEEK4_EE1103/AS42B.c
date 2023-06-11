//ENCRYPTION-DECRYPTION
#include<math.h>
#include<stdio.h>
 
double main ()
{

  int i, j, n, m; //variable declarations
  float M[3][3] = { {1, 4, -3}, {-2, 8, 5}, {3, 4, 7} }, d[3]; //matrix definition
  float x, s[3][3], y[3]; //extra matrix definitons, s takes the inverse numbers, y is the intermediate d matrix in inverse calculation
  void LU ();

  n = 2;

  LU (M, n);//function type def

  printf (" \n");


/*  TO FIND THE INVERSE */

  for (m = 0; m <= 2; m++)
    {
      d[0] = 0.0; //all the non m th elements of the column vector are zeros
      d[1] = 0.0;//all the non m th elements of the column vector are zeros
      d[2] = 0.0;//all the non m th elements of the column vector are zeros
      d[m] = 1.0;//the m th iteration takes value 1
      for (i = 0; i <= n; i++) //row ieration
    {
      x = 0.0;
      for (j = 0; j <= i - 1; j++) //column iteration
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
      s[i][m] = (y[i] - x) / M[i][i];//whichever row is 1, the corresponding column gets the values. while rest are zeros.
    }
    }

/* Print the inverse matrix */
  printf ("The Inverse Matrix of M is as follows\n");
  for (m = 0; m <= 2; m++)
    {
      printf (" %f %f %f \n", s[m][0], s[m][1], s[m][2]);
    }
    printf("\n");
    
//decoding the hidden encrypted message
  int number;
  printf ("Enter the number of letters in the supposed secret message\n");
  scanf ("%i", &number);
  //we assign this length to the variable number
  char ami[number];
  printf ("Enter the numbers which are to be encoded. Also hit enter after inserting each value\n");
  for (int t = 0; t <number / 3; t++)
  //t is the number of column vectors each containing three rws
    {
    float a[3][1];
    //initialise a matrix a with the encoding values i.e the column vector
    for (int ank = 0; ank < 3; ank++)
      {
        for (int anj = 0; anj < 1; anj++)
        scanf ("%f", &a[ank][anj]);
        //each value entered by us is added to the matrix to give a column vector
      }
      
      for (int z = 0; z <= 2; z++)
    {
      double sum = s[z][0] * a[0][0] + s[z][1] * a[1][0] + s[z][2] * a[2][0];
      //the elements of the column vector [3x1] are post multiplied with those of s [3x3] to get another colum vector [3x1]
      int c = 64 + round (sum);
      //we add 64 to the values of the column vector to get the ASCII Values of characters
      if (c==91){
          c=32;
          //this addresses the discrepancy when a left bracket is typed instead of a space. space has ascii 32
      }
      ami[3*t+z]=c;
//this array takes up values of all characters as elements. 3*t+z is a special combination that yields 0,1,2...number methodically.
//this facilitates number and element assignment
      
    }

    }
    printf("THE SECRET MESSAGE IS:\n"); //this prints values of the array a individually
    for(i=0;i<number;i++){
    printf("%c",ami[i]); }
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
