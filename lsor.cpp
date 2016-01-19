#include<iostream>
#include<cmath>
using namespace std;
int main()
{
//constants
float m=0.5,l=2*3.142,ep=0.1,E=1e-5,w=1.25;
//cout<<m<<l<<ep;
int n=20; 

double dx,dy,bet;
dx=(double)2*3.142/(n);
dy=(double)2*3.142/(n);
bet=double(dx*dx)/((1-m*m)*(dy*dy));
//cout<<bet;
float phi[n][n],phi_n[n][n];
int i,j,k;
//making the old values and new values  zero
for(i=0;i<n;i++)
{
for(j=0;j<n;j++)
{
phi[i][j]=0;
phi_n[i][j]=0; 
//cout<<phi_n[i][j];
}
//cout<<endl;
}

//boundary condtions
for(i=0;i<n;i++)
	{
	phi[i][0]=phi[i][1]-dy*cos(dx*i);
	phi_n[i][0]=phi_n[i][1]-dy*cos(dx*i);
//	cout<<phi_n[i][0]<<endl;	
	}
// running across the time
for(k=0;k<10000;k++)
	{
		for(i=1;i<n-1;i++)
			{
				for(j=1;j<n-1;j++)
				{
			phi_n[i][j]=(double)(1/(2*(1+bet)))*(phi[i+1][j]+phi_n[i-1][j]+bet*(phi[i][j+1]+phi_n[i][j-1]));
//			cout<<phi_n[i][j]<<endl;
				}
		//		cout<<phi_n[i][j];
			}
		for(j=0;j<n;j++)
			{
			phi_n[0][j]=phi_n[1][j];
		//	cout<<phi_n[0][j]<<"  ";
			phi_n[n-1][j]=phi_n[n-2][j];
			phi_n[j][n-1]=phi_n[j][n-2];
		//	cout<<phi_n[j][n-1]<<"yo"<<endl;
			}
		
		for(i=0;i<n;i++)
			{
			phi_n[i][0]=phi_n[i][1]-(double)dy*ep*cos(i*dx);
		//	cout<<phi_n[i][0];
			}
			for(i=0;i<n;i++)
			{
			for(j=0;j<n;j++)
			{
			phi_n[i][j]=(1-w)*phi[i][j]+w*phi_n[i][j];
			}
			}
	}
//intializing grids
for(i=0;i<n;i++)
{
for(j=0;j<n;j++)
{
cout<<phi_n[i][j]<< " ";
}
cout<<endl;
}
return 0;
}
