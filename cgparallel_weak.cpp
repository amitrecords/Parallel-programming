#include<iostream>
#include<math.h>
#include<iomanip>
#include<array>
#include<mpi.h>

using namespace std;

//----------------------------------Global Variables----------------------------------------------
#define pi 3.1459
#define MAX 10
float u[MAX][MAX] ={0};
float ut[MAX][MAX] ={0};
const int nsquare=(MAX-2)*(MAX-2);
double x[nsquare] ={0};
int rn,NP;

//-------------------------------FUNCTION DEFINITIONS FOLLOW---------------------------------------
//Function 1:
// Applying the boundary conditions
void boundary_conditions()
   {  int i;
      double dx=1.0/MAX;
      for(i=0;i<MAX;i++)
  	{
	  
	 u[i][0]   = sin(pi*i*dx);
	 u[i][MAX-1] = sin(pi*i*dx)* exp(-pi);
         u[0][i]  = 0;
	 u[MAX-1][i] = 0;
	 
         }
         u[MAX-1][0]=0;
   }
//Function2:----------------------------------------------------------------------------------------
// Congugate_gradient function definition
double* Conjugate_gradient(int start , int en, int n, int rowptac[] ,int colac[] ,double valac[] ,double b[])
   {
	  const int size=n;
	  double r[size],res[size], x[size],d[size],z[size],minv[size][size];
   	  double alpha,alphaden,alden,alnum,alphanum,beta,btnum,betanum;

//----------------------------------Initial declarations-------------------------------------------
      int iter,i,j;
      res[size] ={0};d[size] ={0};
      r[size] ={0};z[size] ={0};minv[size][size]={0};
     for(i=0;i<size;i++)
     { x[i]=0;res[i]=0;r[i]=0;d[i]=0;z[i]=0;
        for(j=0;j<size;j++)
        	{minv[i][j]=0;}

    	 }

    	    for(i=start;i<en;i++)
    	    {

    	    	r[i]= b[i]-res[i];
    	    //	cout<<"before r[]"<<r[i]<<"\n";
    	    }   // Got Ro=b-Ax!


    	   // Creating M-inverse
    	    for(i = start; i < en;i++){
    	    		minv[i][i] = -0.25;
    	    	}
    	   //Computation of z= Minv*r=d
    	    	    for(int i = start ; i < en ;i++){
    	    		    for(int j = 0; j < n; j++){
    	    		       z[i] = z[i] + minv[i][j] * r[j];
    	    			   d[i]=z[i];
    	    			   }                         // got z and d

    	    	    }

//------------------------------------------Iterative search begins!-----------------------------------------------
      for(iter=0;iter<1000;iter++)
  	{
	      // making alpha
           for(i=start;i<en;i++)
           {
         	res[i]=0;
      	        for(j=rowptac[i]-1;j<rowptac[i+1]-1;j++)
      	         {
      	           	res[i]=res[i]+valac[j]*d[colac[j]-1];

      	         }
      		}
        // calclulating the dot product of the denominator
        alphaden=0;
      	for(i=start;i<en;i++)
      	{
      		alphaden=alphaden+res[i]*d[i];
      	}
         MPI_Allreduce(&alphaden,&alden,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
         alphanum=0;

        for(i=start;i<en;i++)
      	{
      		alphanum=alphanum+z[i]*r[i];
      	}

        MPI_Allreduce(&alphanum,&alnum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        alpha=alnum/alden;  //got alpha

        for(i=start;i<en;i++)
      	{
      		x[i]=x[i]+alpha*d[i];
      		r[i]=r[i]-alpha*res[i];
      	}                    // got new values of x and r

        for(i = start ; i < en ;i++){
          for(j = start; j < en; j++){
            z[i] = z[i] + minv[i][j] * r[j];
   		     }// got new z
        }

        // making beta
        betanum=0;
        for(i=start;i<en;i++)
        {
    	  betanum=betanum + z[i]*r[i];
        }
        MPI_Allreduce(&betanum,&btnum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    	beta=btnum/alnum ;

    	for(i=start;i<en;i++)
    	{
    		d[i]=z[i]+beta*d[i];
    	}  // got new d value

  	} // Once this is over, we would have our final value of the unknowns which we name as 'x'


      //take the answer generated from conjugate gradient and store it back to u matrix;
      	 	int k=1;int m=1;
            for(i = 0; k <MAX-1; i++)
              { //cout << "b[" << i << "]\t" << b[0][i] << "\n";
                   u[k][m]=x[i];
                   m++;
                   if(m==MAX-1)
                   { k++;
                     m=1;
  	 	            }
              }

 }
//----------------------------------------MAIN FUNCTION STARTS HERE-----------------------------
int main(int argc, char *argv[])
{ 
      // Values needed for dgesv

    int n;
    int nrhs = 1;
    int k,kmax;
    kmax=MAX-2;
    int i,j;
    double a[nsquare][nsquare]={ 0 },AT[nsquare][nsquare]= { 0 };
    double b[nsquare]={0.0};
    double val[nsquare*nsquare]={ 0 },col_ind[nsquare*nsquare]= { 0 },row_ptr[nsquare*nsquare]= { 0 };
    n=nsquare;

    boundary_conditions();
//Displaying the matrix u	
    	/*cout<<" The matrix u:"<<endl;
    	for(i=0;i<MAX;i++)
    	{  for(j=0;j<MAX;j++)
    	      cout<<setprecision(2)<< u[i][j]<< "\t ";
    	 cout<<endl;
    	}        */
  // Assigning the coefficient matrix values 
//Transposing U
    	for(i=0;i<MAX;i++)
    	{	for(j=0;j<MAX;j++)
    			ut[j][i]=u[i][j];
    	}  
  
	  
// k=i+(j-1)*kmax;

    	for(k=0;k<nsquare;k++)
    	{  
    	  	a[k][k]=-4;
    	  	a[k][k-1]=1;
    	  	a[k][k+1]=1;
    	  	if((k+1) % 4==0 && k != 0)
    	  	 a[k][k+1]=0;
    	  	if((k) % 4==0 && k != 0)
    	  	 a[k][k-1]=0;
    	  	if(k>3)
    	  	 a[k][k-4]=1; 
    	  	if(k<nsquare-4)
    	  	 a[k][k+4]=1;
	 	  
    	}
    	a[nsquare-1][nsquare-1]=-4;
    	cout<< endl;
//--------------------------------------Display of coefficient matrix if required----------------
  /*//Displaying the coefficient matrix
    	   	cout<<"\n The coefficient matrix:"<<endl;
    	    	for(i=0;i<nsquare;i++)
    	    	{  for(j=0;j<nsquare;j++)
    	    	      cout<< a[i][j]<< " ";
    	   	  cout<<endl;
    	   	}*/
//------------------------------------------------------------------------------------------------

 //Compressed Row storage of sparse matrix
    	 k=0;
    	 int ctr=0;
    	int p=0,q=0;
    	int flag=0;
    	for(i=0;i<nsquare;i++)
    	{	for(j=0;j<nsquare;j++)
    		{	if(a[i][j] != 0)
    			 { ctr++;
    			   val[k]=a[i][j];
    			   col_ind[q]=j+1;
    			   if(flag==0)
    			   { row_ptr[p]=k+1;			//creation of row_ptr vector
    			     flag=1;
    			     p++;
    			   }
    			   k++;
    			   q++;
    			 }
    		}

    	 flag=0;
    			 	 // val and col_ind vector created here
    	}
    			double *valac = new double[k];
    	    	int *colac = new int[k];
    	    	int *rowptac = new int[p];
    	        for(i=0;i<k;i++)
    	        {  valac[i]=val[i];
    	           colac[i]=col_ind[i];
    	        }
    	    	for(i=0;i<p;i++)
    	    	        {  rowptac[i]=row_ptr[i];
    	    	        }

  // Assigning values to the matrix B (Ax=B)
     //starting from the bottom inner row
    	{  for(i=0;i<MAX-2;i++)
    	   b[i]= -ut[0][i];
    	   b[0]= -ut[1][0]-ut[0][1]   ;
           b[MAX-3]= -ut[1][MAX-1] - ut[0][MAX-2];
        }  
     //the rows further upwards
    	 

    	 for(j=2;j<=MAX-3;j++) 
    	{
    		 	 
    	  b[MAX-4+j]= -u[j][0];
    	   for(i=1;i<MAX-3;i++)
    	   { 
    	    b[MAX-2+i]= 0   ;
           } 
          b[j*(MAX-2)-1]= -u[j][MAX-1];
        }  
     //Ending at the topmost inner row
       	{  
        	   b[(MAX-2)*(MAX-3)]= -u[MAX-2][0];
         	   b[(MAX-2)*(MAX-2)-1]= -u[MAX-2][MAX-1];
         	   
        }  	 
       	MPI_Init(&argc,&argv);
       	  	MPI_Comm_size( MPI_COMM_WORLD, &NP );
       	    MPI_Comm_rank(MPI_COMM_WORLD, &rn);
       	    double star= MPI_Wtime();
       	    int MAX2=MAX*NP
       	 // ROW BLOCKING-*******
       	         int start=rn*((MAX2-2)/NP);
       	         int en=(int)(rn+1) * (MAX2-2)/NP);
       	            if(rn==(NP-1))
       	          en=en+((MAX2-2) % NP);
       	         //**********************

  // Conjugate gradient solver
	 Conjugate_gradient(start,en,(MAX2-2)*(MAX2-2),rowptac,colac,valac,b);

 	  double diff= MPI_Wtime()- star;
 	  MPI_Barrier(MPI_COMM_WORLD);
	  MPI_Finalize();
	  //Displaying the matrix u
	  	      	cout<<" The matrix u:"<<endl;
	  	       	for(i=0;i<MAX2;i++)
	  	        	{  for(j=0;j<MAX2;j++)
	  	       	      cout<<setprecision(2)<< u[i][j]<< "\t ";
	  	           	 cout<<endl;
	  	           	}
	  cout<< "The Time taken to do this in parallel with "<<NP<<"processors is"<< diff<<endl;
	  cin.get();

      return 0; 
}      
