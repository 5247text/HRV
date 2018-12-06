/*
 *--------------------------------------------------------------------------------------------------------- 
 *功能：运用AR模型求功率谱密度 
 *设置的参数：
 *			  RRcount,输入序列的长度； 
 *            fft,FFT的点数 ； 
 *			  degree,AR模型的阶数； 
 *计算结果：结果保存于f,PSD数组（f,PSD为全局变量），计算结束后读取f,PSD即可，读取长度为fft/2(FFT是左右对称的，读取一半数据即可) ；
 *使用说明：调用AR_PSD()函数即可，返回1表示计算完成即可读取结果，返回0表示等待数据；
 *			init_AR_PSD()函数对全局变量进行初始化。 
 *------------------------------------------- ---------------------------------------------------------------
 *
 *note: write by Vigro,@2017,snarc
 * 
 */
#include <iostream>
#include <cmath>
#include "stdlib.h"
#include <vector>
#include "marple.cpp"
#include "string.h"
#include "spline.h"
#define RRcount 60 //进行谱估计的RR间期个数 
#define fft 1024 //FFT的点数，设置的点数为256*freq 
#define degree 16 //AR模型的阶数 

const int N = 1024; 
const float PI = 3.1416;
int RR_count = 0;
float freq = 4; // 重采样的频率 
double RR[RRcount];
double xi[RRcount];
double parama[fft];
double parama_imag[fft] ={0};
double PSD[fft] = {0};
double f[fft] = {0};

/*
 *功能：准备数据,当数据个数超过RRcount时，自动更新RR间期的序列（队列的形式） 
 *参数：rr,单个RR值
 *返回值：1表示数据已经准备好；0表示数据还未准备好； 
 */ 
int getdata(float rr){
	int i;
	if(RR_count < RRcount-1){
		RR[RR_count] =rr;
		RR_count++;
		return 0;
	}else if(RR_count == RRcount-1){
		RR[RR_count] =rr;
		RR_count++;
		return 1;
	}else{
		for(i=0;i<RRcount-1;i++){
			RR[i] = RR[i+1];
		}
		RR[RRcount-1] = rr;
		RR_count++;
		return 1;
	}

}
/*
 *功能：计算RR间期的总时长
 *参数：无
 *返回值：输入的RR间期总时长 
 */
double sum(){
	int i;
	double sum = 0;
	for(i = 0; i<RRcount;i++){
		sum += RR[i];
	}
	return sum;
}
/*
 *功能：计算需要插值的点数 
 *参数：无
 *返回值：插值序列的点数 
 *
 */ 
int numinterp(){
	double num =0;
	int num_interp =0;
	double t_sum;
	t_sum = sum();
	num = (t_sum - RR[1])*freq;
	num_interp = floor(num);
	return num_interp; 
}
/*
 *功能：获取RR间期的x坐标，存于数组xi 
 * 参数：无
 *返回值：无 
 */
void getx(){
	int i;
	double sum = 0;
	for(i = 0;i<RRcount;i++){
		sum = sum+RR[i];
		xi[i] = sum;
	}
}
/*
 *功能： 获取插值后的x坐标，存于数组xo
 *参数： *xo,用于存储插值后的x坐标；num,插值的个数 
 *返回值：无 
 */

void getxo(double *xo,int num){
	int i = 0;
	double s = RR[0];
	for(i = 0;i < num;i++ ){
		xo[i] = s;
		s += 1/freq;
	}
} 

/*
 *功能：交换两个变量的值 
 *参数：&a,&b 
 *返回值：无 
 */  
inline void swap (double &a, double &b)
{
    double t;
    t = a;
    a = b;
    b = t;
}
/*
 *功能：位反转置换 
 *参数：xreal,输入序列的实部；
 *		ximag,输入序列的虚部；
 *		n,个数；
 *返回值：无 
 *
 */
void bitrp (double xreal [], double ximag [], int n)
{
    // 位反转置换 Bit-reversal Permutation
    int i, j, a, b, p;

    for (i = 1, p = 0; i < n; i *= 2)
        {
        p ++;
        }
    for (i = 0; i < n; i ++)
        {
        a = i;
        b = 0;
        for (j = 0; j < p; j ++)
            {
            b = (b << 1) + (a & 1);    // b = b * 2 + a % 2;
            a >>= 1;        // a = a / 2;
            }
        if ( b > i)
            {
            swap (xreal [i], xreal [b]);
            swap (ximag [i], ximag [b]);
            }
        }
}
/*
 *功能：FFT 
 *参数：xreal,输入序列的实部；
 *		ximag,输入序列的虚部；
 *		n,个数；
 *返回值：无 
 */
void FFT(double xreal [], double ximag [], int n)
{
    // 快速傅立叶变换，将复数 x 变换后仍保存在 x 中，xreal, ximag 分别是 x 的实部和虚部
    double wreal [N / 2], wimag [N / 2], treal, timag, ureal, uimag, arg;
    int m, k, j, t, index1, index2;
    
    bitrp (xreal, ximag, n);

    // 计算 1 的前 n / 2 个 n 次方根的共轭复数 W'j = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
    arg = - 2 * PI / n;
    treal = cos (arg);
    timag = sin (arg);
    wreal [0] = 1.0;
    wimag [0] = 0.0;
    for (j = 1; j < n / 2; j ++)
        {
        wreal [j] = wreal [j - 1] * treal - wimag [j - 1] * timag;
        wimag [j] = wreal [j - 1] * timag + wimag [j - 1] * treal;
        }

    for (m = 2; m <= n; m *= 2)
        {
        for (k = 0; k < n; k += m)
            {
            for (j = 0; j < m / 2; j ++)
                {
                index1 = k + j;
                index2 = index1 + m / 2;
                t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
                treal = wreal [t] * xreal [index2] - wimag [t] * ximag [index2];
                timag = wreal [t] * ximag [index2] + wimag [t] * xreal [index2];
                ureal = xreal [index1];
                uimag = ximag [index1];
                xreal [index1] = ureal + treal;
                ximag [index1] = uimag + timag;
                xreal [index2] = ureal - treal;
                ximag [index2] = uimag - timag;
                }
            }
        }
}

/*
 *功能：求解复数的绝对值 
 *参数：x_real,复数的实部；
 *		x_imag,复数的虚部； 
 *返回值：复数的绝对值； 
 */
double complex_abs(double x_real,double x_imag){
	double x;
	x = sqrt((pow(x_real,2)+pow(x_imag,2)));
	return x;		
}
/*
 *功能：调整AR模型参数ak序列，添加a0 = 1，将原有参数后移一位； 
 *参数：无； 
 *返回值：无 
 */
void exchange(){
	int i = 0;
	for(i =16;i>0;i--){
		parama[i] = parama[i-1];
		
	}
	printf("\n");
	parama[0] = 1;
}

/*
 *功能：求功率谱密度 
 *参数：无 
 *返回值：无 
 */
void psd(double FY){
	exchange();
	
	FFT(parama,parama_imag,fft);
	int i;
	double x_real,x_imag,fft_mag;
	for(i = 0;i<fft;i++){
		x_real = parama[i];
		
		x_imag = parama_imag[i];
		
		fft_mag = complex_abs(x_real,x_imag);
		if (i == 0)
			PSD[i] = FY/(fft_mag*fft_mag*4);
		else
		
		PSD[i] = FY/(fft_mag*fft_mag*2);
		
	}
}
/*
 *功能：功率谱的频率序列即功率谱的x轴 
 *参数：无 
 *返回值：无 
 */
void getFreq(){
	int i;
	for(i = 0; i<fft;i++){
		f[i] = freq/fft*i;
	}
}
/*
 *功能：去除直流分量（对插值后的RR间期序列减去均值）
 *参数：num_interp, RR间期的个数；*yo,插值后的RR间期序列；
 *返回值：无 
 */
void minusmean(int num_interp,double *yo){
	int i;
	double mean = sum() / RRcount;
	for(i=0;i<num_interp;i++){
		yo[i] = yo[i] - mean;
	}
}

/*
 *功能：运用AR模型求解功率谱，调用此函数即可 
 *参数：rr,单个RR序列； 
 *返回值：计算完成返回1，数据等待返回0； 
 */
int AR_PSD(double rr){
	int i,j;
	if(getdata(rr)){
		
		int num_interp;
		num_interp = numinterp() + 1;
		double *xo = new double[num_interp];
		double *yo = new double[num_interp];
		getx();
		getxo(xo,num_interp);
	 	std::vector<double> x(RRcount),y(RRcount);
	 	for(i = 0; i<RRcount;i++){
	 	x[i] = xi[i];
	 	y[i] = RR[i];
	 	}
	 	tk::spline s;
		// s.set_boundary(tk::spline::second_deriv,0.0,tk::spline::first_deriv,0.0,false);
	 	s.set_points(x,y);
	 	for(j = 0; j<num_interp;j++){
	 	yo[j] = s(xo[j]);
		}
	 	minusmean(num_interp,yo);
	 	int max;
	 	double FPE,FY;
		Marple(yo,num_interp,degree,&max,parama,&FPE,&FY);
		psd(FY);
		getFreq();
		free(xo);
		free(yo); 
		return 1;
	}else{
		return 0;
	} 
	
}
/*
 *功能：初始化函数，初始化全局变量； 
 *参数；无 
 *返回值：无 
 */
    
void init_AR_PSD(){
	RR_count = 0;
	freq = 4;
 	memset(RR,0,RRcount*sizeof(double));
 	memset(xi,0,RRcount*sizeof(double));
 	memset(parama,0,fft*sizeof(double));
 	memset(parama_imag,0,fft*sizeof(double));
 	memset(PSD,0,fft*sizeof(double));
	memset(f,0,fft*sizeof(double));	
} 

int main(){
//	init_AR_PSD();
	int i;

	int n = 0;
	double rr[60] = {0.880000000000000,0.900000000000000,0.870000000000000,0.855000000000000,0.935000000000000,0.870000000000000,0.850000000000000,0.870000000000000,0.825000000000000,0.940000000000000,0.895000000000000,0.615000000000000,0.650000000000000,0.905000000000000,0.995000000000000,1.05499990000000,1,1.00500000000000,0.990000000000000,0.945000000000000,0.980000000000000,1,0.945000000000000,1,1.03000000000000,0.955000000000000,0.970000000000000,0.965000000000000,0.930000000000000,0.945000000000000,0.935000000000000,0.800000000000000,0.785000000000000,0.815000000000000,0.930000000000000,0.800000000000000,0.700000000000000,0.675000000000000,0.690000000000000,0.735000000000000,0.945000000000000,1.03000000000000,1.09000000000000,1.09500000000000,0.980000000000000,0.950000000000000,1.06500010000000,1.14500000000000,1.08000010000000,0.910000000000000,0.990000000000000,1.01500000000000,1.09500000000000,1.01500000000000,1.02000000000000,1.13500000000000,0.905000000000000,0.800000000000000,0.745000000000000,0.995000000000000};
	while(n<60){
		if(AR_PSD(rr[n])){
			for(i = 0;i<RRcount;i++){
			printf("%f ",f[i]);
			printf("%lf\n",PSD[i]);
		}
		printf("\n");
		}
		
		
		
		n++;
	}


		
}
