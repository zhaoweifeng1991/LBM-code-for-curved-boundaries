//Couette flow, convex scheme

#include<iostream> 
#include<cmath> 
#include<cstdlib> 
#include<iomanip> 
#include<fstream> 
#include<sstream> 
#include<string> 
#include<stdio.h> 
 
using namespace std; 
const int Q=9; 
const int NX = 66;   // the mesh size should be 64 x 64, two meshes are added for the convenience of computation
const int NY = 66; 
const double U=0.1; 
const double pi=3.1415926;
const double beta=0.2;


 
int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};    //9 directions
int ne[Q]={0,3,4,1,2,7,8,5,6};                                                //back directions
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36}; 
double rho[NX+1][NY+1],u[NX+1][NY+1][2],u0[NX+1][NY+1][2],f[NX+1][NY+1][Q],ff[NX+1][NY+1][Q],F[NX+1][NY+1][Q],xlabel[NX+1][NY+1],ylabel[NX+1][NY+1]; 
int i,j,k,ip,jp,n,q_flag; 
double c,Re,dx,dy,Lx,Ly,D,dt,rho0,p0,tau_f,niu,error,E_r,y,yy1,yy2,kk,b,cc,x1,x2,x,q; 

double iq,jq,AA,BB,CC,DD,EE,rr,uu,vv,Center_x,Center_y;
double R1,R2;  //radius of the circles 


double ell=1.0;  // a parameter to control the shape of the boundary: ell*(x-x0)^2+(y-y0)^2=R^2, ell=1 for circle
double s_nu,s_q,SS,cs_2;

double abs( double i);
void comput_q (int i, int j, int ip, int jp, double R);




 
void init(); 
double feq(int k,double rho,double u[2]); 
void evolution(); 
void output(int m); 
void Error(); 
int flag[NX+1][NY+1];


 
int main() 
{ 
  using namespace std; 

  

  init();   //initiate


  for(n=0; ;n++) 
  { 
    evolution(); 
    if(n%100==0) 
    { 
      Error(); 
	}

	if(n>10 && n%100==0)
	{
      cout<<"The"<<n<<"th computational result:"<<endl<<"The u,v of point (NX/2,NY/2)is : " 
     
      <<setprecision(6)<<u[NX/2][NY/2][0]<<","<<u[NX/2][NY/2][1]<<endl; 
      cout<<"The max relative error of uv is:" 
        <<setiosflags(ios::scientific)<<error<<" "<<E_r<<endl;         
	} 


//	if(n%1000==0)  
//		output(n);

	if(error<1.0e-10) 
	{
        Error();
		cout<<"The max relative error of uv is:" 
        <<setiosflags(ios::scientific)<<error<<" "<<E_r<<endl; 

		output(n+1);
		break;
	}

  } 
 
  return 0; 
} 



void init() 
{ 
  
  Lx=1.6; 
  Ly=1.6; 
  dx=Lx/(NX-2.0); 
  dy=dx; 
  niu=0.02;
  SS=1.0;   //tau
  s_nu=-1.0/SS;
  s_q=-8.0*(2+s_nu)/(8+s_nu);             
  dt=(SS -0.5)/3.0 *dx*dx /niu;
  c=dx/dt;

  cout<<"U/c = "<<U/c<<"\n";


  R2=Lx/2.0;  // radius 
  R1=R2*beta;

  rho0=1.0;
  cs_2=c*c/3.0;
  
  for(i=0;i<=NX;i++)    //corrdinates of the lattice nodes
  for(j=0;j<=NY;j++) 
  {
	  xlabel[i][j]=(i-NX/2)*dx ;
	  ylabel[i][j]=(j-NY/2)*dy ;
  }

  Center_x=0.0;
  Center_y=0.0;

  

  for (i=0;i<=NX;i++)    //tag different lattice nodes
	  for (j=0;j<=NY;j++)
	  {
		  flag[i][j]=0;           

		  if(    ( (  ell*(xlabel[i][j]-Center_x)*(xlabel[i][j]-Center_x) + (ylabel[i][j]-Center_y)*(ylabel[i][j]-Center_y) ) < R2*R2 )
			  && ( (  ell*(xlabel[i][j]-Center_x)*(xlabel[i][j]-Center_x) + (ylabel[i][j]-Center_y)*(ylabel[i][j]-Center_y) ) > R1*R1 ))
		  {
		        flag[i][j]=1;	   //fluid lattice nodes between the two circles
		  }

		  else if(    ( (  ell*(xlabel[i][j]-Center_x)*(xlabel[i][j]-Center_x) + (ylabel[i][j]-Center_y)*(ylabel[i][j]-Center_y) ) >= R2*R2 ) )

			     flag[i][j] = 2;   //outside of the circle with radius R2
	  }

  
 
  for(i=0;i<=NX;i++)    //initialization of DF
    for(j=0;j<=NY;j++) 
	if(flag[i][j]==1)
    {
		
      u[i][j][0]= 0.0; 
      u[i][j][1]= 0.0; 
      rho[i][j]= rho0; 

      for(k=0;k<Q;k++) 
      { 
        f[i][j][k]=feq(k,rho[i][j],u[i][j]); 
      } 
    }
} 



 
double feq(int k,double rho,double u[2])   // equilibrium distribution
{ 
  double eu,uv,feq; 
  eu=(e[k][0]*u[0]+e[k][1]*u[1]); 
  uv=(u[0]*u[0]+u[1]*u[1]); 
  feq=w[k]*(rho + rho0* (3.0*eu/c+4.5*eu*eu/c/c-1.5*uv/c/c) ); 
  return feq; 
} 


 
void evolution() 
{ 


  
 for(i=0;i<=NX;i++)    
    for(j=0;j<=NY;j++) 
		if(flag[i][j]==1)
		{
			for(k=0;k<Q;k++) 
				F[i][j][k]=f[i][j][k]
				           +
				           s_nu*(  0.5*(f[i][j][k]+f[i][j][ne[k]]) - 0.5*(feq(k,rho[i][j], u[i][j])+feq(ne[k],rho[i][j], u[i][j]))   )
						   +
				           s_q*(  0.5*(f[i][j][k]-f[i][j][ne[k]]) - 0.5*(feq(k,rho[i][j], u[i][j])-feq(ne[k],rho[i][j], u[i][j]))   );    
			//collision
		}
	

 for(i=0;i<=NX;i++)    
    for(j=0;j<=NY;j++) 
		if(flag[i][j]==1)
		{
				for(k=0;k<Q;k++) 
				{
					ip=i-e[k][0]; 
					jp=j-e[k][1];
				

					if( flag[ip][jp]==0 )
					{
					   
						comput_q(ip,jp,i,j,R1);   //the distance between the node (ip,jp) and the boundary

					   iq=xlabel[i][j]-(dx-q)*double(e[k][0]);
					   jq=ylabel[i][j]-(dx-q)*double(e[k][1]);


						q = 1.0 - q/dx;     // the ratio $gamma$ in the paper


						uu =  U* jq / sqrt(iq*iq+jq*jq) ; //velocity at the boundary
						vv = -U* iq / sqrt(iq*iq+jq*jq) ;


                        AA= 2.0*q/(1.0+2.0*q);
						BB = 1.0-AA;
						CC = 1.0-AA+BB;

						ff[i][j][k] = AA*F[i][j][k] +BB*f[i][j][ne[k]] + CC*w[k]*rho0*3.0/c*(e[k][0]*uu+e[k][1]*vv);

					}
					
					
					else if(flag[ip][jp]==2)
						{
					   
						comput_q(i,j,ip,jp,R2);   //the distance between the node (i,j) and the boundary
						

						q = q/dx;     // the ratio $gamma$ in the paper

                        AA= 2.0*q/(1.0+2.0*q);
						BB = 1.0-AA;
						CC = 1.0-AA+BB;

						ff[i][j][k] = AA*F[i][j][k] +BB*f[i][j][ne[k]] ; // + CC*w[k]*rho0*3.0/c*(e[k][0]*uu+e[k][1]*vv) = 0

					}

					else

					ff[i][j][k]=F[ip][jp][k];					

				}
			
		}

    


for(i=0;i<=NX;i++)    
    for(j=0;j<=NY;j++) 
		if(flag[i][j]==1)
        { 
			
          rho[i][j]=0; 
          u[i][j][0]=0; 
          u[i][j][1]=0; 
		
          for(k=0;k<Q;k++) 
          { 
            f[i][j][k] = ff[i][j][k];
            rho[i][j]+=f[i][j][k]; 
            u[i][j][0]+=c*e[k][0]*f[i][j][k]; 
            u[i][j][1]+=c*e[k][1]*f[i][j][k]; 
          } 

		  u[i][j][0]/=rho0; 
          u[i][j][1]/=rho0; 

        
		} 

}



double abs( double i)
{
	if(i>=0.0) return i;
	else 
		return -i;
}



void comput_q (int i, int j, int ip, int jp, double R)  //Compute the distance between the node (i,j) and the boundary
{
     if (ip==i)
	 {   
		 yy1  = abs( Center_y+sqrt( (R*R-(xlabel[i][j]-Center_x)*(xlabel[i][j]-Center_x))/ell )-ylabel[i][j] );
		 yy2  = abs( Center_y-sqrt( (R*R-(xlabel[i][j]-Center_x)*(xlabel[i][j]-Center_x))/ell )-ylabel[i][j] );

		 if(yy1<=yy2) q=yy1;
		 else q=yy2;
	 }


	 else
	 {
         kk  =  (ylabel[ip][jp]-ylabel[i][j])/(xlabel[ip][jp]-xlabel[i][j]);
		 
		 b   =  (  2.0*Center_x - 2.0*ell*kk*(ylabel[i][j]-kk*xlabel[i][j]-Center_y)  ) / (ell*kk*kk+1.0);

		 cc  =  ( ell*(ylabel[i][j]-kk*xlabel[i][j]-Center_y)*(ylabel[i][j]-kk*xlabel[i][j]-Center_y)+Center_x*Center_x-R*R  ) / (ell*kk*kk+1.0);  //×¢Òâderta>=0

         x1  =  abs( (b+sqrt(b*b-4.0*cc))/2.0-xlabel[i][j] );
         x2  =  abs( (b-sqrt(b*b-4.0*cc))/2.0-xlabel[i][j] );

		 if(x1<=x2) q=x1;
		 else q=x2;
		 
	 }
}




void output(int m)    //output the data
        { 
          ostringstream name; 
          name<<"TaylorGreen"<<"beta_"<<beta<<"Mesh"<<NX<<"_"<<NY<<"tau_"<<SS<<".dat"; 
          ofstream out(name.str().c_str()); 
          out<<"Title= \"TaylorGreen\"\n"<<"VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"U0\",\"V0\",\"p\" ,\"flag\" \n"<<"ZONE T=\"BOX\",I=" 
            <<NX+1<<",J="<<NY+1<<",F=POINT"<<endl;
            for(j=0;j<=NY;j++) 
              for(i=0;i<=NX;i++) 
              { 
                out<<setprecision(15)<<xlabel[i][j]<<" "<<ylabel[i][j]<<" "<<u[i][j][0]<<" "<<  u[i][j][1]<<" "<<u0[i][j][0]<<" "<<  u0[i][j][1]<<" "<<rho[i][j]<<" "<<flag[i][j]<<endl; 
		
              } 


		  ostringstream name2; 
          name2<<"Computation_X_TaylorGreen"<<"beta_"<<beta<<"NX_"<<NX<<"NY_"<<NY<<"tau_"<<SS<<".dat"; 
          ofstream out2(name2.str().c_str()); 
              for(i=ceil(R1/dx)+NX/2;  i<=NX-3 ; i++) 
              { 
                out2<<setprecision(15)<<(xlabel[i][NY/2]-R1)/(R2-R1)<<" "<<  -u[i][NY/2][1]/U<<endl; 		
              } 

			  
        } 
 
      
void Error()   //compute the error
        { 
          double temp1,temp2,temp3,temp4; 
		  double R_temp, U_temp;
          temp1=0; 
          temp2=0;
		  temp3=temp4=0;
		  
          for(i=1;i<NX;i++) 
            for(j=1;j<NY;j++)
			{
			   if(flag[i][j]==1)
			   {

                  temp1+=( (u[i][j][0]-u0[i][j][0])*(u[i][j][0]-u0[i][j][0])+(u[i][j][1]-u0[i][j][1])*(u[i][j][1]-u0[i][j][1])); 
                  temp2+=(u0[i][j][0]*u0[i][j][0]+u0[i][j][1]*u0[i][j][1]); 	


				  u0[i][j][0] = u[i][j][0] ;
				  u0[i][j][1] = u[i][j][1] ;


				  R_temp = sqrt(xlabel[i][j]*xlabel[i][j] + ylabel[i][j]*ylabel[i][j]);
				  U_temp = U*beta/(1.0-beta*beta) * ( R2/R_temp - R_temp/R2  );
				  uu =  U_temp* ylabel[i][j] / R_temp ;
				  vv = -U_temp* xlabel[i][j] / R_temp ; 
				  temp3+=( (u[i][j][0]-uu)*(u[i][j][0]-uu)+(u[i][j][1]-vv)*(u[i][j][1]-vv)); 
                  temp4+=(uu*uu+vv*vv); 
				
			   }
			}


            temp1=sqrt(temp1); 
            temp2=sqrt(temp2); 
            error=temp1/(temp2+1e-30); 

			temp3=sqrt(temp3); 
            temp4=sqrt(temp4); 
            E_r=temp3/(temp4+1e-30); 

        } 
 