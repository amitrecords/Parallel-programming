#include<iostream>
#include<math.h>
#include<iomanip>
#include<array>

using namespace std;


#define pi 3.1459
#define MAX 5
float u[MAX][MAX] ={0};
float ut[MAX][MAX] ={0};
const int nsquare=(MAX-2)*(MAX-2);
double x[nsquare] ={0};
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

// Congugate_gradient function definition
double* Conjugate_gradient(int n, double (*a)[nsquare] , double b[9])
   {  const int size=n;
	  double r[size],res[size], x[size],d[size],z[size],minv[size][size];
   	  double alpha,alphaden,alphanum,beta,betanum;

 //---------------------------------Initial declarations-------------------------------------------
      int iter,i,j;
      res[size] ={0};d[size] ={0};
      r[size] ={0};z[size] ={0};minv[size][size]={0};
      // making alpha
                       for(i=0;i<size;i++)
                       {
                     	      for(j=0;j<size;j++)
                  	         {
                  	           	res[i]=res[i]+a[i][j]*x[j];
                  		      }
                  		}
    	    for(i=0;i<size;i++)
    	    { r[i]= b[i]-res[i];
    	    //cout<<" res [i]"<< res[i];
    	    }   // Got Ro=b-Ax!
    	    res[size] ={0};

    	   // Creating M-inverse
    	    for(i = 0; i < size;i++){
    	    		minv[i][i] = 1/a[i][i];
    	    	}
    	    //Computation of z= Minv*r=d
    	    	    for(int i = 0 ; i < n ;i++){
    	    		    for(int j = 0; j < n; j++){
    	    			   z[i] = z[i] + minv[i][j] * r[j];
    	    			   d[i]=z[i];
    	    		     }// got z and d
    	    	     }
 //----------------------------------------------------------------------------------------
      for(iter=0;iter<10000;iter++)
  	{    if (r[0]==0)
  	     { cout<< "Solution found"; break;
  	     }

	    // making alpha
           for(i=0;i<(size);i++)
           {
         	res[i]=0;
      	        for(j=0;j<size;j++)
      	         {
      	           	res[i]=res[i]+a[i][j]*d[j];
      		      }
      		}
      // calclulating the dot product of the denominator
       alphaden=0;
      	for(i=0;i<size;i++)
      	{
      		alphaden=alphaden+res[i]*d[i];
      	}
       alphanum=0;
      	for(i=0;i<size;i++)
      	{
      		alphanum=alphanum+z[i]*r[i];
      	}
      alpha=alphanum/alphaden;  //got alpha

      for(i=0;i<size;i++)
      	{
      		x[i]=x[i]+alpha*d[i];
      		r[i]=r[i]-alpha*res[i];
      	}// got new values of x and r

      for(i = 0 ; i < n ;i++){
          for(j = 0; j < n; j++){
            z[i] = z[i] + minv[i][j] * r[j];
   		     }// got new z
      }
      // making beta
      betanum=0;
      for(i=0;i<size;i++)
       {
    	  betanum=betanum + z[i]*r[i];
       }
    	beta=betanum/alphanum ;
      for(i=0;i<size;i++)
    	{
    		d[i]=z[i]+beta*d[i];
    	}// got new d value

  	}
      cout<<"size of b:"<<sizeof(x)/sizeof(x[0])<<endl;
      cout<<" The matrix r:"<<endl;
          	for(i=0;i<size;i++)
          	{
          	      cout<< r[i]<< " ";
                   }
      //take the answer generated from conjugate gradient;
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

int main()
{ 
      // Values needed for dgesv

      int n;
      int nrhs = 1;
      int k,kmax;
      kmax=MAX-2;

      double a[nsquare][nsquare]={ 0 },AT[nsquare][nsquare]= { 0 };
      double b[nsquare]={0};
      double val[nsquare*nsquare]={ 0 },col_ind[nsquare*nsquare]= { 0 },row_ptr[nsquare*nsquare]= { 0 };
      
      
      int i,j;
      n=nsquare;
      
            
      boundary_conditions();	    
//Displaying the matrix u	
    	cout<<" The matrix u:"<<endl;
    	for(i=0;i<MAX;i++)
    	{  for(j=0;j<MAX;j++)
    	      cout<<setprecision(2)<< u[i][j]<< "\t ";
    	 cout<<endl;
    	}        
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

  //Displaying the coefficient matrix
    	   	cout<<"\n The coefficient matrix:"<<endl;
    	    	for(i=0;i<nsquare;i++)
    	    	{  for(j=0;j<nsquare;j++)
    	    	      cout<< a[i][j]<< " ";
    	   	  cout<<endl;
    	   	}

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
    	    	double *colac = new double[k];
    	    	double *rowptac = new double[p];
    	        for(i=0;i<k;i++)
    	        {  valac[i]=val[i];
    	           colac[i]=col_ind[i];
    	        }
    	    	for(i=0;i<p;i++)
    	    	        {  rowptac[i]=row_ptr[i];
    	    	        }

    	cout<<endl<<"the ctr value is :"<< ctr;
    	//Display the vecors
    	cout<<" \nthe vector val: \n";
    	for(i=0;i<k;i++)
    	{  cout<< valac[i]<<" ";}
    	cout<<" \nthe vector col_ind: \n";
    	for(i=0;i<k;i++)
    	{  cout<< colac[i]<<" ";}
    	cout<<" \nthe vector row_ptr: \n";
    	for(i=0;i<p;i++)
    	{  cout<< rowptac[i]<<" ";}
    	cout<< endl;


  // Assigning values to the matrix B (Ax=B)
     //starting from the bottom inner row
    	{  for(i=1;i<MAX-2;i++)
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
       	
//Displaying the matrix B
    	/*cout<<" The matrix B:"<<endl;
    	for(i=0;i<nsquare;i++)
    	{  
    	      cout<< b[i]<< " ";
    	 }     */

	      cout<<"size of u:"<<sizeof(u)/sizeof(u[0][0])<<endl;      
	      cout<<"size of A:"<<sizeof(a)/sizeof(a[0][0])<<endl;  
	      cout<<"size of b:"<<sizeof(b)/sizeof(b[0])<<endl;

  // Conjugate gradient solver
	 Conjugate_gradient(nsquare,a,b);



	//Displaying the matrix u
	      	cout<<" The matrix u:"<<endl;
	       	for(i=0;i<MAX;i++)
	        	{  for(j=0;j<MAX;j++)
	       	      cout<<setprecision(1)<< u[i][j]<< "\t ";
	           	 cout<<endl;
	           	}

	  cin.get();
      return 0; 
}      
