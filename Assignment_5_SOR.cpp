#include<stdio.h>
#include<conio.h>
#include<math.h>

#define imax 100
#define jmax 100
#define PI 3.14159265
#define omg 1.33

FILE *fT;
char *x="'",*y_="STREAM",*z="'";

double l,dX,dY,X[imax+1][jmax+1],Y[imax+1][jmax+1],phi[imax+1][jmax+1],phi_new[imax+1][jmax+1],phi_new_SOR[imax+1][jmax+1];
double epsilon,f,e;

main()
{
	l=2*PI;
	dX=l/imax;
	dY=l/jmax;
	
	for(int i=0;i<=imax;i++)			//grid formation
	{
		for(int j=0;j<=jmax;j++)
		{
			X[i][j]=i*dX;
			Y[i][j]=j*dY;
		}
	}
	
	for(int i=1;i<=imax-1;i++)			//initialization of phi
	{
		for(int j=1;j<=jmax-1;j++)
		{
			phi[i][j]=0;
		}
	}
	
	for(int j=0;j<=jmax;j++)			//west boundary
	{
		phi[0][j]=0;
	}
	for(int j=0;j<=jmax;j++)			//east boundary
	{
		phi[imax][j]=0;
	}
	for(int i=1;i<=imax-1;i++)			//north boundary
	{
		phi[i][jmax]=0;
	}
	for(int i=1;i<=imax-1;i++)			//south boundary
	{
		phi[i][0]=-0.1*dY*cos(X[i][0]);
	}
	
	for(int j=0;j<=jmax;j++)			//west boundary
	{
		phi_new[0][j]=0;
	}
	for(int j=0;j<=jmax;j++)			//east boundary
	{
		phi_new[imax][j]=0;
	}
	for(int i=1;i<=imax-1;i++)			//north boundary
	{
		phi_new[i][jmax]=0;
	}
	for(int i=1;i<=imax-1;i++)			//south boundary
	{
		phi_new[i][0]=-0.1*dY*cos(X[i][0]);
	}
	
	epsilon=1;
	while(epsilon>=0.01)
	{
		for(int i=1;i<=imax-1;i++)
		{
			for(int j=1;j<=jmax-1;j++)
			{
				phi_new[i][j]=0.5*(((0.75/(dX*dX))*(phi[i+1][j]+phi_new[i-1][j]))+((1/(dY*dY))*(phi[i][j+1]+phi_new[i][j-1])))/((0.75/(dX*dX))+(1/(dY*dY)));
			}
		}
		for(int j=0;j<=jmax;j++)			//west boundary
		{
			phi_new[0][j]=phi_new[1][j];
		}
		for(int j=0;j<=jmax;j++)			//east boundary
		{
			phi_new[imax][j]=phi_new[imax-1][j];
		}
		for(int i=1;i<=imax-1;i++)			//north boundary
		{
			phi_new[i][jmax]=phi_new[i][jmax-1];
		}
		for(int i=1;i<=imax-1;i++)			//south boundary
		{
			phi_new[i][0]=phi_new[i][1]-0.1*dY*cos(X[i][0]);
		}
		
		for(int i=0;i<=imax;i++)
		{
			for(int j=0;j<=jmax;j++)
			{
				phi_new_SOR[i][j]=((1-omg)*phi[i][j]+omg*phi_new[i][j]);
			}
		}
		
		f=0;								//epsilon calculation
		e=0;
		
		for(int i=1;i<=imax-1;i++)
		{
			for(int j=1;j<=jmax-1;j++)
			{
				e=e+abs(phi[i][j]);
				//printf("%lf\n",phi[i][j]);
				f=f+fabs(phi_new_SOR[i][j]-phi[i][j]);
			}
		}
		epsilon=f;
		printf("%lf\t%lf\t%lf\n",f,e,epsilon);
		
		for(int i=0;i<=imax;i++)
		{
			for(int j=0;j<=jmax;j++)
			{
				phi[i][j]=phi_new_SOR[i][j];
			}
		}
	}
	
	fT=fopen("PHI_DISTRIBUTION_SOR.dat","a");		//.dat for PHI by ADI
	//fprintf(fT,"TITLE=%s%s%s\n",x,y_,z);
	//fprintf(fT,"VARIABLES=X,Y,PHI\n");
	//fprintf(fT,"ZONE I=%d,J=%d,F=POINT\n",imax,jmax);

	for(int j=jmax;j>=0;j--)
    {
       	for(int i=imax;i>=0;i--)
       	{
  	   		fprintf(fT,"%lf\t %lf\t %lf\n",X[i][j],Y[i][j],phi[i][j]);
  		}
    }
    fclose(fT);
}
