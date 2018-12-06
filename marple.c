#include <math.h>
#include <stdio.h>

void Marple( double *X,int N,int MMAX,int *MAX,double *A,double *FPE,double *FY)
{
	double E,DEN,EO,Q1,Q2,Q3,G,W,Q4,Q5,F,B,H,S,V,U,Q6,SAVE1,SAVE2;
	double YN,XN,GR,SAVE3,SAVE4,FPL,ALPHA,C1,C2,C3,C4,C5,C6;
	int SSS,NM,N1,K11,M1,NMK,M1K,NK,NMK1,MK,M2,N1M1;
	int M,MM,I,N1K1,K1,K,KM11;
	double Q7,Y1,Y2,Y3,Y4,DEL;
	double ZZ[1000],C[1000],D[1000],R[1000];
	SSS=0;
	
	XN=N;
	EO=0;
	for( K=1;K<=N;K++)
		EO=EO+X[K-1]*X[K-1];
	
	EO=EO*2;
	Q1=1/EO;Q2=Q1*X[0];
	G=Q1*(X[0]*X[0]);W=Q1*(X[N-1]*X[N-1]);
	DEN=1-G-W;Q4=1/DEN;
	Q5=1-G;Q6=1-W;F=X[0];B=X[N-1];
	H=Q2*X[N-1];S=Q2*X[N-1];U=Q1*(X[N-1]*X[N-1]);V=Q2*X[0];
	E=EO*DEN;Q1=1/E;C[0]=Q1*X[0];D[0]=Q1*X[N-1];
	M=1;
	SAVE1=0;N1=N+1;NM=N-1;
	for(K=1;K<=NM;K++)
	{ K11=K+1;
	SAVE1=SAVE1+X[K11-1]*X[K-1];
	}
	R[0]=2*SAVE1;  A[0]=-Q1*R[0]; E=E*(1-A[0]*A[0]);
	FPL=E/(2*(XN-1));FPL=FPL*((XN+1)/(XN-1));
	
L3: if( M>=MMAX ) goto L9;
	
	
	M1=M+1;F=X[M1-1];B=X[NM-1];
	for(K=1;K<=M;K++)
	{
		NMK=NM+K;M1K=M1-K;
		F=F+X[M1K-1]*A[K-1];B=B+X[NMK-1]*A[K-1];
	}
	Q1=1/E;Q2=Q1*F;Q3=Q1*B;
	for(I=1;I<=M;I++)
	{
		K=M+1-I;K1=K+1;C[K1-1]=C[K-1]+Q2*A[K-1];
		D[K1-1]=D[K-1]+Q3*A[K-1];
	}
	C[0]=Q2;D[0]=Q3;Q7=S*S;Y1=F*F;
	Y3=B*B;Y2=V*V;Y4=U*U;
	G=G+Y1*Q1+Q4*(Y2*Q6+Q7*Q5+2*V*S*H);
	W=W+Y3*Q1+Q4*(Y4*Q5+Q7*Q6+2*S*H*U);
	S=0;H=0;V=0;U=0;
	for(K=0;K<=M ;K++)
	{
		K1=K+1;NK=N-K;NMK1=NM+K;
		H=H+X[NMK1-1]*C[K1-1];S=S+X[NK-1]*C[K1-1];
		U=U+X[NK-1]*D[K1-1];V=V+X[K1-1]*C[K1-1];
	}
	
	Q5=1-G;Q6=1-W;
	DEN=Q5*Q6-H*H;
	if( DEN<=0 ) goto L9;
	Q4=1/DEN;
	
	Q1=Q1*Q4;
	ALPHA=1/(1+(Y1*Q6+Y3*Q5+2*H*F*B)*Q1);
	E=ALPHA*E;C1=Q4*(F*Q6+B*H);C2=Q4*(B*Q5+H*F);
	C3=Q4*(V*Q6+H*S);C4=Q4*(S*Q5+V*H);
	C5=Q4*(S*Q6+H*U);C6=Q4*(U*Q5+S*H);
	for(K=1 ;K<= M ;K++)
	{
		K1=K+1;
		A[K-1]=ALPHA*(A[K-1]+C1*C[K1-1]+C2*D[K1-1]);
	}
	
	M2=M / 2+1;
	for(K=1;K<= M2;K++)
	{
		MK=M+2-K;SAVE1=C[K-1];
		SAVE2=D[K-1];SAVE3=C[MK-1];
		SAVE4=D[MK-1];
		C[K-1]=C[K-1]+C3*SAVE3+C4*SAVE4;
		D[K-1]=D[K-1]+C5*SAVE3+C6*SAVE4;
		if( K!=MK )
		{
			C[MK-1]=C[MK-1]+C3*SAVE1+C4*SAVE2;
			D[MK-1]=D[MK-1]+C5*SAVE1+C6*SAVE2;
		}
	}
	
	M=M+1;NM=N-M;M1=M-1;DEL=0;
	N1M1=N1-M;C1=X[N1M1-1];C2=X[M-1];
	for(I=1;I<=M1 ;I++)
	{
		K=M1+1-I;K11=K+1;N1K1=N1-K;
		R[K11-1]=R[K-1]-X[N1K1-1]*C1-X[K-1]*C2;
		DEL=DEL+R[K11-1]*A[K-1];
	}
	SAVE1=0;
	for(K=1;K<=NM;K++)
	{
		KM11=K+M;SAVE1=SAVE1+X[KM11-1]*X[K-1];
	}
	R[0]=2*SAVE1;DEL=DEL+R[0];
	Q2=-DEL/E;A[M-1]=Q2;M2=M / 2;
	for(K=1;K<=M2 ;K++)
	{
		MK=M-K;SAVE1=A[K-1];A[K-1]=A[K-1]+Q2*A[MK-1];
		if( MK!=K ) A[MK-1]=A[MK-1]+Q2*SAVE1;
	}
	Y1=Q2*Q2;E=E*(1-Y1);
	
	if( Y1>=1 )
		goto L9;
	YN=M;
	*FPE=E/(2*(XN-YN))*(XN+YN+1)/(XN-YN-1);
//	if( *FPE>=FPL ) goto L3;
	if( M<MMAX ) goto L3; 
	MM=M;GR=E;FPL=*FPE;
	for (I=1;I<=M ;I++)
		ZZ[I-1]=A[I-1];goto L3;
L9: for( I=1;I<=MM ;I++)
		A[I-1]=ZZ[I-1];
	*MAX=MM;*FY=GR/(2*(XN-MM));
	*FPE=FPL;
	
}
