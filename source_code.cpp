#include <iostream>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>//数组计算
#include <algorithm>//数组最大最小函数(未使用）
using namespace std;

//定义初始变量
#define N 201
#define M 201
#define R 30
#define Pi 3.1415926
#define V0 10
#define C 1


int main()
{
	//1.定义全场上的x值，y值数组
	float x [201][201];
	float y [201][201];

	//计算物面x
	for(int i=0;i<201;i++)
	{
		x[i][0]=0.5*(1+cos(1.8*i*Pi/180));
	}
	//计算物面y
	for(int i=0;i<201;i++)
	{
		y[i][0]=0.6*(-0.1015*pow(x[i][0],4)+0.2843*pow(x[i][0],3)
			 -0.3576*pow(x[i][0],2)-0.1221*x[i][0]+0.2969*pow(x[i][0],0.5));
	}

	for(int i=100;i<201;i++)
	{
		y[i][0]=-y[i][0];
	}

	/*
	验证程序
	for(int i=0;i<101;i++)
	{
		for(int j=0;j<101;j++)
		{
			cout<<x0[i]<<" "<<y0[j]<<" ";
		}
		cout<<endl;
	}
	*/

	//2.计算远场
	for(int i=0;i<201;i++)
	{
		x[i][200]=30*cos(1.8*i*Pi/180);
		y[i][200]=30*sin(1.8*i*Pi/180);
	}

	//3.插值中间点
	for(int i=0;i<201;i++)
	{
		for(int j=0;j<201;j++)
		{
			x[i][j]=(x[i][200]-x[i][0])*j/200+x[i][0];
			y[i][j]=(y[i][200]-y[i][0])*j/200+y[i][0];
		}
	}
	//由上网格初始化完成

	
	/*//4.输出初始网格，检查
	ofstream file;
	file.open("01.dat");
	file << "variables = X,Y" << endl;
	file << "zone t = naca, i = " << 201 << ", j = " << 201 << endl;   
	for(int j=0;j<201;j++)    
	{       
		for(int i=0;i<201;i++)       
		{          
			file << x[i][j] << " " <<  y[i][j]<<endl;     
		}  
	}  
	file.close();
	*/

	//5.内场网格迭代(简单迭代收敛速度过慢）
	float x_max=1;
	float y_max=1;
	float a;
	float b;
	float c;
	int n=0;
	float x_res [201][201];
	float y_res [201][201];
	//float x0_max=1;
	//float y0_max=1;
	//网格长j均取1，边界条件外

	while(x_max>10e-6)
	{
			n++;
			
		//5.1 将上一次内场的数据保存在了残差数组中
			for(int i=0;i<201;i++)
			{
				for(int j=0;j<201;j++)
				{
					x_res[i][j]=x[i][j];
					y_res[i][j]=y[i][j];
				}
			}

		
		//5.3 内场网格迭代
			for(int i=1;i<200;i++)
			{
				for(int j=1;j<200;j++)
				{
					a=pow((x[i][j+1]-x[i][j-1])/2,2)+pow((y[i][j+1]-y[i][j-1])/2,2);
					b=(x[i+1][j]-x[i-1][j])*(x[i][j+1]-x[i][j-1])/4+(y[i+1][j]-y[i-1][j])*(y[i][j+1]-y[i][j-1])/4;
					c=pow((x[i+1][j]-x[i-1][j])/2,2)+pow((y[i+1][j]-y[i-1][j])/2,2);
					x[i][j]=0.5*(
									a*(x[i+1][j]+x[i-1][j])
									+c*(x[i][j+1]+x[i][j-1])
									-0.5*b*(x[i+1][j+1]+x[i-1][j-1]-x[i+1][j-1]-x[i-1][j+1])
								)/(a+c);
					y[i][j]=0.5*(
									a*(y[i+1][j]+y[i-1][j])
									+c*(y[i][j+1]+y[i][j-1])
									-0.5*b*(y[i+1][j+1]+y[i-1][j-1]-y[i+1][j-1]-y[i-1][j+1])
								)/(a+c);

				}
			}
		//5，4 边界网格迭代
			for(int j=1;j<200;j++)
			{
					
					a=pow((x[0][j+1]-x[0][j-1])/2,2)+pow((y[0][j+1]-y[0][j-1])/2,2);
					b=(x[1][j]-x[199][j])*(x[0][j+1]-x[0][j-1])/4+(y[1][j]-y[199][j])*(y[0][j+1]-y[0][j-1])/4;
					c=pow((x[1][j]-x[199][j])/2,2)+pow((y[1][j]-y[199][j])/2,2);
					x[0][j]=0.5*(
									a*(x[1][j]+x[199][j])
									+c*(x[0][j+1]+x[0][j-1])
									-0.5*b*(x[1][j+1]+x[199][j-1]-x[1][j-1]-x[199][j+1])
								)/(a+c);
					y[0][j]=0.5*(
									a*(y[1][j]+y[199][j])
									+c*(y[0][j+1]+y[0][j-1])
									-0.5*b*(y[1][j+1]+y[199][j-1]-y[1][j-1]-y[199][j+1])
								)/(a+c);

					//割缝
					x[200][j]=x[0][j];
					y[200][j]=y[0][j];
					
			}


		/*//5.5 原点迭代
			a=pow((-x[0][2]+4*x[0][1]-3*x[0][0])/2,2)+pow((-y[0][2]+4*y[0][1]-3*y[0][0])/2,2);
			b=(x[1][0]-x[199][0])*(-x[0][2]+4*x[0][1]-3*x[0][0])/4+(y[1][0]-y[199][0])*(-y[0][2]+4*y[0][1]-3*y[0][0])/4;
			c=pow((x[1][0]-x[199][0])/2,2)+pow((y[1][0]-y[199][0])/2,2);
			x[0][0]=0.5*(
								a*(x[1][0]+x[199][0])
								+c*(-x[0][2]+4*x[0][1]-3*x[0][0])
								-0.5*b*(x[1][1]+x[199][2]-3*x[199][1]+3*x[199][0]-(x[1][2]-3*x[1][1]+3*x[1][0])-x[199][1])
							)/(a+c);
			y[0][0]=0.5*(
								a*(y[1][0]+y[199][0])
								+c*(-y[0][2]+4*y[0][1]-3*y[0][0])
								-0.5*b*(y[1][1]+y[199][2]-3*y[199][1]+3*y[199][0]-(y[1][2]-3*y[1][1]+3*y[1][0])-y[199][1])
							)/(a+c);

			//割缝
			x[200][0]=x[0][0];
			y[200][0]=y[0][0];
		*/

		//5.7残差
			for(int i=0;i<201;i++)
			{
				for(int j=0;j<201;j++)
				{
					x_res[i][j]=x_res[i][j]-x[i][j];
					y_res[i][j]=y_res[i][j]-y[i][j];
				}
			}
			x_max=x_res[0][0];
			y_max=y_res[0][0];
			for(int i=0;i<201;i++)
			{
				for(int j=0;j<201;j++)
				{
					if(abs(x_res[i][j])>x_max)
					{
						x_max=abs(x_res[i][j]);
					}
					if(abs(y_res[i][j])>y_max)
					{
						y_max=abs(y_res[i][j]);
					}
				}
			}
			if(y_max>x_max)
			{
				x_max=y_max;
			}


			if(n%100==0)
			{
			cout<<n<<" "<<x_max<<" "<<endl; //（监测迭代过程）
			}
	}
	cout<<n<<endl; //监测迭代次数

	
	// 6. 输出网格，检查
	/*
	ofstream file;
	file.open("03.dat");
	file << "variables = X,Y" << endl;
	file << "zone t = naca, i = " << 201 << ", j = " << 201 << endl;   
	for(int j=0;j<201;j++)    
	{       
		for(int i=0;i<201;i++)       
		{          
			file << x[i][j] << " " <<  y[i][j]<<endl;     
		}  
	}  
	file.close();
	*/


	//以上网格生成完成，下面进行流场计算

	//7 解流场
	float phi [201][201];
	float a_phi;
	float b_phi;
	float c_phi;
	float phi_res [201][201];
	float phi_max=1;
	int m=0;

	while(phi_max>10e-6)
	{
		m++;

		//7.1 残差
		for(int i=0;i<201;i++)
		{
			for(int j=0;j<201;j++)
			{
					phi_res[i][j]=phi[i][j];
			}
		}

		//7.3 迭代
		for(int i=1;i<200;i++)
		{
			for(int j=1;j<200;j++)
			{
				
				a_phi=pow((x[i][j+1]-x[i][j-1])/2,2)+pow((y[i][j+1]-y[i][j-1])/2,2);
				b_phi=(x[i+1][j]-x[i-1][j])*(x[i][j+1]-x[i][j-1])/4+(y[i+1][j]-y[i-1][j])*(y[i][j+1]-y[i][j-1])/4;
				c_phi=pow((x[i+1][j]-x[i-1][j])/2,2)+pow((y[i+1][j]-y[i-1][j])/2,2);
				
				phi[i][j]=0.5*(
									a_phi*(phi[i+1][j]+phi[i-1][j])
									+c_phi*(phi[i][j+1]+phi[i][j-1])
									-0.5*b_phi*(phi[i+1][j+1]+phi[i-1][j-1]-phi[i+1][j-1]-phi[i-1][j+1])
								)/(a_phi+c_phi);

			}
		}

		//7.4 远场
		for(int i=0;i<201;i++)
		{
		phi[i][200]=V0*x[i][200];
		}

		//7.5 割缝
		for(int j=1;j<200;j++)
		{
			a_phi=pow((x[0][j+1]-x[0][j-1])/2,2)+pow((y[0][j+1]-y[0][j-1])/2,2);
			b_phi=(x[1][j]-x[199][j])*(x[0][j+1]-x[0][j-1])/4+(y[1][j]-y[199][j])*(y[0][j+1]-y[0][j-1])/4;
			c_phi=pow((x[1][j]-x[199][j])/2,2)+pow((y[1][j]-y[199][j])/2,2);

			phi[0][j]=0.5*(
									a_phi*(phi[1][j]+phi[199][j])
									+c_phi*(phi[0][j+1]+phi[0][j-1])
									-0.5*b_phi*(phi[1][j+1]+phi[199][j-1]-phi[1][j-1]-phi[199][j+1])
							)/(a_phi+c_phi);
			phi[200][j]=phi[0][j];
		}

		//7.6 物面
		//隐格子：x[i][j+1]-x[i][j-1]=-x[i][2]+4*x[i][1]-3*x[i][0]
		//x[i][-1]=x[i][2]-3*x[i][1]+3*x[i][0]
		//phi[i][0]=((b/c)*(phi[i-1][0]-phi[i+1][0])+4*phi[i][1]-phi[i][2])/3
		for(int i=1;i<200;i++)
		{
			a_phi=pow((-x[i][2]+4*x[i][1]-3*x[i][0])/2,2)+pow((-y[i][2]+4*y[i][1]-3*y[i][0])/2,2);
			b_phi=(x[i+1][0]-x[i-1][0])*(-x[i][2]+4*x[i][1]-3*x[i][0])/4+(y[i+1][0]-y[i-1][0])*(-y[i][2]+4*y[i][1]-3*y[i][0])/4;
			c_phi=pow((x[i+1][0]-x[i-1][0])/2,2)+pow((y[i+1][0]-y[i-1][0])/2,2);

			phi[i][0]=(
						(b_phi/c_phi)*(phi[i-1][0]-phi[i+1][0])
						+4*phi[i][1]
						-phi[i][2]
						)/3;
		}

		//7.6 原点

		a_phi=pow((-x[0][2]+4*x[0][1]-3*x[0][0])/2,2)+pow((-y[0][2]+4*y[0][1]-3*y[0][0])/2,2);
		b_phi=(x[1][0]-x[199][0])*(-x[0][2]+4*x[0][1]-3*x[0][0])/4+(y[1][0]-y[199][0])*(-y[0][2]+4*y[0][1]-3*y[0][0])/4;
		c_phi=pow((x[1][0]-x[199][0])/2,2)+pow((y[1][0]-y[199][0])/2,2);

		phi[0][0]=(
					(b_phi/c_phi)*(phi[199][0]-phi[1][0])
					+4*phi[0][1]
					-phi[0][2]
					)/3;

		phi[200][0]=phi[0][0];



		//7.7获得残差
		for(int i=0;i<201;i++)
		{
			for(int j=0;j<201;j++)
			{
					phi_res[i][j]=phi_res[i][j]-phi[i][j];
			}
		}

		phi_max=phi_res[0][0];
		for(int i=0;i<201;i++)
		{
			for(int j=0;j<201;j++)
			{
				if(abs(phi_res[i][j])>phi_max)
				{
					phi_max=abs(phi_res[i][j]);
				}
			}
		}
		if(m%100==0)
			{
			cout<<m<<" "<<phi_max<<" "<<endl; //（监测迭代过程）
			}
		
	}
	cout<<m<<endl;//监测迭代次数
	//速度势迭代完成


	
	//8.求解速度u,v

	//8.1 求解雅可比行列式

	float J [201][201];


	for(int i=1;i<200;i++)
	{
		for(int j=1;j<200;j++)
		{
			J[i][j]=(
					(x[i+1][j]-x[i-1][j])*(y[i][j+1]-y[i][j-1])
					-(x[i][j+1]-x[i][j-1])*(y[i+1][j]-y[i-1][j])
					)/4;
		}
	}
	
	
	for(int j=1;j<200;j++)
	{
		J[0][j]=((x[1][j]-x[199][j])*(y[0][j+1]-y[0][j-1]))/4
					-((x[0][j+1]-x[0][j-1])*(y[1][j]-y[199][j]))/4;
		J[200][j]=J[0][j];
	}

	for(int i=1;i<200;i++)
	{
		J[i][0]=((x[i+1][0]-x[i-1][0])*(-y[i][2]+4*y[i][1]-3*y[i][0]))/4
				-((-x[i][2]+4*x[i][1]-3*x[i][0])*(y[i+1][0]-y[i-1][0]))/4;
	}

	J[0][0]=((x[1][0]-x[199][0])*(-y[0][2]+4*y[0][1]-3*y[0][0]))/4
			-((-x[0][2]+4*x[0][1]-3*x[0][0])*(y[1][0]-y[199][0]))/4;
	J[200][0]=J[0][0];

	
	//8.2 求解uv


	float u [201][201];
	float v [201][201];

	for(int i=1;i<200;i++)
	{
		for(int j=1;j<200;j++)
		{
			u[i][j]=(
					(phi[i+1][j]-phi[i-1][j])*(y[i][j+1]-y[i][j-1])/4
					-(phi[i][j+1]-phi[i][j-1])*(y[i+1][j]-y[i-1][j])/4
					)/J[i][j];
			v[i][j]=(
					(phi[i][j+1]-phi[i][j-1])*(x[i+1][j]-x[i-1][j])/4
					-(phi[i+1][j]-phi[i-1][j])*(x[i][j+1]-x[i][j-1])/4
					)/J[i][j];
		}
	}

	for(int j=1;j<200;j++)
	{
		u[0][j]=(
					(phi[1][j]-phi[199][j])*(y[0][j+1]-y[0][j-1])/4
					-(phi[0][j+1]-phi[0][j-1])*(y[1][j]-y[199][j])/4
					)/J[0][j];
		v[0][j]=(
					(phi[0][j+1]-phi[0][j-1])*(x[1][j]-x[199][j])/4
					-(phi[1][j]-phi[199][j])*(x[0][j+1]-x[0][j-1])/4
					)/J[0][j];
		u[200][j]=u[0][j];
		v[200][j]=v[0][j];
	}

	for(int i=1;i<200;i++)
	{
		u[i][0]=(
				(phi[i+1][0]-phi[i-1][0])*(-y[i][2]+4*y[i][1]-3*y[i][0])/4
				-(-phi[i][2]+4*phi[i][1]-3*phi[i][0])*(y[i+1][0]-y[i-1][0])/4
				)/J[i][0];

		v[i][0]=(
				(-phi[i][2]+4*phi[i][1]-3*phi[i][0])*(x[i+1][0]-x[i-1][0])/4
				-(phi[i+1][0]-phi[i-1][0])*(-x[i][2]+4*x[i][1]-3*x[i][0])/4
				)/J[i][0];
	}

	u[0][0]=(
			(phi[1][0]-phi[199][0])*(-y[0][2]+4*y[0][1]-3*y[0][0])/4
			-(-phi[0][2]+4*phi[0][1]-3*phi[0][0])*(y[1][0]-y[199][0])/4
			)/J[0][0];
	v[0][0]=(
			(-phi[0][2]+4*phi[0][1]-3*phi[0][0])*(x[1][0]-x[199][0])/4
			-(phi[1][0]-phi[199][0])*(-x[0][2]+4*x[0][1]-3*x[0][0])/4
			)/J[0][0];
	u[200][0]=u[0][0];
	v[200][0]=v[0][0];

	//9. 求解Cp
	float Cp[201][201];
	for(int i=0;i<201;i++)
	{
		for(int j=0;j<201;j++)
		{
			Cp[i][j]=1-(pow(u[i][j],2)+pow(v[i][j],2))/pow(V0,2);
		}
	}

	// 10. 输出结果
	ofstream file;
	file.open("04.dat");
	file << "variables = x,y,phi,u,v,Cp" << endl;
	file << "zone t = naca, i = " << 201 << ", j = " << 201 << endl;   
	for(int j=0;j<201;j++)    
	{       
		for(int i=0;i<201;i++)       
		{          
			file << x[i][j]<<" "<<y[i][j]<<" "<< phi[i][j] << " " <<  u[i][j]<<" "<< v[i][j]<<" "<<Cp[i][j]<<endl;     
		}  
	}  
	file.close();

	file.open("05.dat");
	file<<"variables = x, Cp"<<endl;
      
	for(int i=0;i<101;i++)       
	{          
		file << x[i][0]<<" "<<Cp[i][0]<<endl;     
	}  

	file.close();



system ("pause");
return 0;
}
