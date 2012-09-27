#include <cmath>
#include <iostream>
#include <lpsolve/lp_lib.h>

using std::cout;
using std::endl;

long fact(long n) {
  long i, f = 1;

  for(i=1; i<=n; ++i)
     f *= i;

  return f;
}

long binom(long n, long k) {
    if(k>n) return 0;
    return ( fact(n)/(fact(k)*fact(n-k)) );
}


long kraw(long j, long x, long n){
    double res =0;
    long s;
    for(s=0;s<=j;s++)
        {
            if( ((s%2) == 0))
            {
                res += (long)pow(3,j-s) * binom(x,s) * binom(n-x,j-s);
            }
        else
            {
                res += (-1) * (long)pow(3,j-s) * binom(x,s) * binom(n-x,j-s);
            }
    }
    return res;
}

double b_(long j, long n, long K, long *aj)
{
    double res = 0;
    long r;
    for(r=0; r <= n; r++)
        {
            res+= kraw(j,r,n)*aj[r];
        }
    res /= pow(2,n);
    res *= K;
    return res;
}


double s_(long j, long n, long K, long *aj)
{
    double res = 0;
    long r;
    for(r=0; r <= n; r++)
        {
            if((r % 2) == 0)
                {
                    res+= kraw(j,r,n)*aj[r];
                }
            else
                {
                    res+= (-1)* kraw(j,r,n)*aj[r];
                }
        }
    res /= pow(2,n);
    res *= K;
    return res;
}

int lp_shadow(long n, long K, long d, long coeff)
{
  lprec *lp;
  long Ncol = n+1, i,j, ret = 0;  
  REAL *row = NULL;
  lp = make_lp(0,Ncol);
  if(lp == NULL) ret=1;

  if(ret == 0)
      {
          row = (REAL *) malloc((Ncol + 1) * sizeof(*row));
          if(row == NULL) ret = 2;
      }

  //Make all variables integers
  for(i=1;i<=Ncol; i++)
      {
          if(!set_int(lp,i,1)) ret = 5;
      }

  if(ret == 0) set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */

  //Initialize colno, row.
  for(j=0;j<Ncol;j++)
      {
          row[j] = 0;
      }

  // Eqn 76

  if(ret == 0)
      {
          row[1] = 1;
          if(!add_constraintex(lp,Ncol+1,row,NULL, EQ, 1)) ret = 3;
          
          for(j=1;j<=n;j++)
              {
                  for(i=0;i<Ncol;i++) row[i+1] = 0;
                  row[j+1] = 1;
                  if(!add_constraintex(lp,Ncol,row, NULL, GE, 0)) ret = 3;
              }
      }


  // Eqn 78: \Leftrightarrow B_j - A_j = 0

  if(ret == 0)
      {
          for(j=0;j<d;j++)
              {
                  for(i=0;i<Ncol;i++) row[i+1] = 0;
                  for(i=0;i<Ncol;i++) row[i+1]= K/(double)pow(2,n)*kraw(j,i,n);
                  row[j+1] -= 1;
                  if(!add_constraintex(lp,Ncol,row, NULL, EQ, 0)) ret = 3;
              }
      }

  // Eqn 79: \Leftrightarrow B_j - A_j >= 0

  if(ret == 0)
      {
          for(j=d;j<=n;j++)
              {
                  for(i=0;i<Ncol;i++) row[i+1] = 0;
                  for(i=0;i<Ncol;i++) row[i+1]= K/(double)pow(2,n)*kraw(j,i,n);
                  row[j+1] -= 1;
                  if(!add_constraintex(lp,Ncol,row, NULL, GE, 0)) ret = 3;
              }
      }

  // Eqn 81: Shadow

  if(ret == 0)
      {
          for(j=0;j<=n;j++)
              {
                  for(i=0;i<Ncol;i++) row[i] = 0;
                  for(i=0;i<Ncol;i++)
                      {
                          row[i+1]= K/(double)pow(2,n)*kraw(j,i,n);
                          if((i%2) != 0) row[i+1] *= -1;
                      }
                  if(!add_constraintex(lp,Ncol,row, NULL, GE, 0)) ret = 3;
              }
      }


  if(ret == 0) {
    set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

    // Set objective function to coeff-th coefficient

    /*
    for(i=0;i<Ncol;i++) row[i] = 0;
    if(coeff < 0) row[(-coeff)] = -1;
    else  row[coeff] = 1;

    if(!set_obj_fnex(lp, j, row, NULL)) ret = 4;    

    */
    for(i=0;i<Ncol;i++) row[i+1] = 1;
    if(!set_obj_fnex(lp, j, row, NULL)) ret = 4;    
  }

  if(ret == 0) {
      set_maxim(lp);

      cout << endl << "n = " << n << ", K = " << K << ", d = " << d << endl << "---------------------------------------------" << endl;

      write_LP(lp, stdout);
      /* I only want to see important messages on screen while solving */
      //      set_verbose(lp, IMPORTANT);

      
      ret = solve(lp);
      if(ret != INFEASIBLE)
          {
              get_variables(lp, row);
              cout << endl;
              for(j = 0; j < Ncol; j++)
                  cout << "a[" << j << "] = " << row[j+1] << ", ";
          }
      else cout << endl << "infeasible. ";
      ret = 0;
  }


  cout  << endl << "---------------------------------------------" << endl;
  cout << std::flush;

  if(row != NULL) free(row);
  if(lp != NULL) delete_lp(lp);

  return ret;
}




int main()
{

    /*
    int i,j,k;
    for(i=1;i<6;i++)
        {
            for(j=1;j<6;j++)
                {
                    for(k=1;k<6;k++)
                        {
                            cout << endl << kraw(i,j,k);
                        }
                }
        }
    cout << endl;
    return 0;
    */
    long r;
    return lp_shadow(7,2,1,2);
    
}
