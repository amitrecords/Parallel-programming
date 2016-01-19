#include<iostream>
#include<math.h>
#include<iomanip>
#include<array>
extern "C" {
     void dgesv_(int *n, int *nrhs,  double *a,  int  *lda,  
           int *ipivot, double *b, int *ldb, int *info) ;
}

using namespace std;

#define pi 3.1459
#define MAX 10
float u[MAX][MAX] ={0};
float ut[MAX][MAX] ={0};
// Applying the boundary conditions
void boundary_conditions()
   {  int i;
      double dx=1.0/MAX;
      for(i=0;i<MAX;i++)
  	{
	  
	 u[i][0]   = sin(pi*i*dx);
	 u[i][MAX-1] = sin(pi*i*dx)* exp(-pi);
         u[0][i]  = 0.0;
	 u[MAX-1][i] = 0.0;
	 
         }
         u[MAX-1][0]=0.0;
   }

int main()
{ 
      // Values needed for dgesv
	cout<<"hi";
      int n;
      int nrhs = 1;
      int k,kmax;
      kmax=MAX-2;
      const int nsquare=(MAX-2)*(MAX-2);
      double a[nsquare][nsquare]={ 0.0 },AT[nsquare][nsquare]= { 0.0 };
      double b[nsquare]={0.0};
      int lda = nsquare;
      int ldb = nsquare;
      int ipiv[nsquare];
      int info;     
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
    	
    	for(k=1;k<nsquare;k++)
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
    	a[0][0]=-4;
    	a[0][4]=1;
    	a[0][1]=1;
    	a[nsquare-1][nsquare-1]=-4;
    	cout<< endl;
//Transposing A
    	for(i=0;i<nsquare;i++)
    	{	for(j=0;j<nsquare;j++)
    			AT[j][i]=a[i][j];
    	}
//Displaying the coefficient matrix	
   	cout<<"\n The coefficient matrix:"<<endl;
    	for(i=0;i<nsquare;i++)
    	{  for(j=0;j<nsquare;j++)
    	      cout<< a[i][j]<< " ";
   	  cout<<endl;
   	}      	
    	
  // Assigning values to the matrix B (Ax=B)
     //starting from the bottom inner row
    	{  for(i=1;i<MAX-2;i++)
    	   b[i]= -ut[0][i];
    	   b[0]= -ut[1][0]-ut[0][1]   ;
           b[MAX-3]= -ut[1][MAX-1] - ut[0][MAX-2];
        }  
     //the rows further upwards
    	 int ctr=0;
    	 
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
    	cout<<" The matrix B:"<<endl;
    	for(i=0;i<nsquare;i++)
    	{  
    	      cout<< b[i]<< " ";
    	}  
	      cout<<"size of u:"<<sizeof(u)/sizeof(u[0][0])<<endl;      
	      cout<<"size of A:"<<sizeof(a)/sizeof(a[0][0])<<endl;  
	      cout<<"size of b:"<<sizeof(b)/sizeof(b[0])<<endl;

// Solve the linear system
	         dgesv_(&n, &nrhs, &AT[0][0], &lda, ipiv, &b[0], &ldb, &info);

	         // Check for success
	            if(info == 0)
	            {
	               // Write the answer
	            	k=1;int m=1;
	               cout << "The answer is\n";
	               for(i = 0; k <MAX-1; i++)
	                 { //cout << "b[" << i << "]\t" << b[0][i] << "\n";
	                   u[k][m]=b[i];
	                   m++;
	                   if(m==MAX-1)
	                   { k++;
	                    m=1;
	                   }
	                 } 
	            }
	            else
	            {
	            // Write an error message
	               cerr << "dgesv returned error " << info << "\n";
	            }

	            //Displaying the matrix u	
	                	cout<<" The matrix u:"<<endl;
	                	for(i=0;i<MAX;i++)
	                	{  for(j=0;j<MAX;j++)
	                	      cout<<setprecision(2)<< u[i][j]<< "\t ";
	                	 cout<<endl;
	                	}        	            
	  cin.get();
      return 0; 
}      
