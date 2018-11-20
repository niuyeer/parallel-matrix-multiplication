#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>
//# define MAX 100000000
MPI_Status status;
int **A,**B, **C; //C=A*A
int *a,*b,*c; //the buffer of each processpr
int n=8;
int p1=3,p2=4,nk=2; //
int np; //the size of each proceesor
int p,rank; //the num of proceesor; thea rank
int *tempa, *tempb;
void ProduceAC(); //
void PrintAC();//print
void ScatterA();// scaater A matrix to each processor
void MainProcess(); //cannon
void collectC(); //Collect matrix C
void Mutiply(); //
void Printa();
void Printc();
void matrix_multiplication();


int det1( int nn, int **aa)
{
  // int bb [n][n] = {{0}};
 int **bb;
  // bb = (int**)malloc(n*sizeof(int*));
     bb = (int**)malloc(n*sizeof(int*));
    int i=0;
   for(i = 0; i < n; i++)
    {
        bb[i] = (int*)malloc(n*sizeof(int));
    }
    
    int j = 0;
    int sum=0;
    int x = 0, cc = 0, pp = 0;
    if (nn == 1)
        return aa[0][0];
    for (i = 0; i < nn; i++)
    {
        for (cc = 0; cc < nn - 1; cc++)
        {
            for (j = 0; j < nn - 1; j++)
            {
                if (cc < i) {
                    pp = 0;
                }
                else {
                    pp = 1;
                }
                bb[cc][j] = aa[cc + pp][j + 1];
                //printf("bb is %d\n",bb[cc][j]);
            }
        }
        if (i % 2 == 0) {
            x = 1;
        }
        else {
            x = (-1);
        }
        int de=det1(nn-1, bb);
       // printf("bb is %d",bb[cc][j]);
       //  printf("det is %d\n",de);
      //  printf("sum is %d,%d,%d,%d\n",sum,aa[i][0],de,x);
        sum += aa[i][0] * de* x;
    }
   // printf("sum is %d",sum);
    return sum;
}


int main(int argc, char *argv[])
{
    int i;
    double starttime,endtime;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
    {
        printf("Please input the size n= ");
        fflush(stdout);
        scanf("%d", &n);
        printf("%d\n",n);
        
        printf("Please input the number of 1 in 10 (0~9)= ");
        fflush(stdout);
        scanf("%d", &p1);
        printf("%d\n",p1);
        
        printf("Please input the number of 1 in 10 (0~9)= ");
        fflush(stdout);
        scanf("%d", &p2);
        printf("%d\n",p2);
        
        printf("problem input parameter k= ");
        fflush(stdout);
        scanf("%d", &nk);
        printf("%d\n",nk);
        
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // n = atoi(argv[1]);
    np = n/(int)sqrt(p);
    a = (int*)malloc(np*np*sizeof(int));
    b = (int*)malloc(np*np*sizeof(int));
    c = (int*)malloc(np*np*sizeof(int));
    memset(c, 0, np*np*sizeof(int));
    tempa = (int*)malloc(np*np*sizeof(int));
    tempb = (int*)malloc(np*np*sizeof(int));
    if(rank == 0)
    {
        
        A = (int**)malloc(n*sizeof(int*));
        B = (int**)malloc(n*sizeof(int*));
        C = (int**)malloc(n*sizeof(int*));
        
        for(i = 0; i < n; i++)
        {
            A[i] = (int*)malloc(n*sizeof(int));
            B[i] = (int*)malloc(n*sizeof(int));
            C[i] = (int*)malloc(n*sizeof(int));
        }
        ProduceAC();
        ScatterA();
    }
    else
    {
        MPI_Recv(a, np*np, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(b, np*np, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
    }
    //printf("My rank= %d\n",rank);
   
    starttime=MPI_Wtime();
    MainProcess(); //cannon
    
    if(rank == 0)
    {
        printf("My rank= %d\n",rank);
        collectC();
        
        PrintAC();
        matrix_multiplication();
        printf("My rank= %d\n",rank);
        printf("My det is %d\n", det1( n, C ) );
        
        endtime=MPI_Wtime();
        printf("time used: %lf\n",endtime - starttime);
        for(i = 0; i < n; i++)
        {
            free(A[i]);
            free(B[i]);
            free(C[i]);
        }
        free(A);
        free(B);
        free(C);
    }
    else
    {
        MPI_Send(c, np*np, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }
    free(a);
    free(b);
    free(c);
    free(tempa);
    free(tempb);
    MPI_Finalize();
    return 0;
}
void ProduceAC()
{
    srand((unsigned)time(NULL));
    int i;
    for(i=0; i<n; i++)
    {
        int j;
        for( j=0; j<n; j++)
        {
            // A[i][j]=1;
            int a=(int)rand()%10;
            int b=a;
            
            // printf("AA %d\t",b);
            if (p1>b) {
                A[i][j]=1;
                B[i][j]=1;
            }
            else
            {
                if(p1+p2>b)
                {
                    A[i][j]=-1;
                    B[i][j]=-1;
                }
                else
                {
                    A[i][j]=0;
                    B[i][j]=0;
                }
            }
            C[i][j]=0;
        }
    }
    
}
void PrintAC()
{
    printf("A matrix:\n");
    int i;int j;
    for( i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%d\t",A[i][j]);
        }
        printf("\n");
    }
    printf("B matrix:\n");
    for( i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%d\t",B[i][j]);
        }
        printf("\n");
    }
    
    printf("C matrix:\n");
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%d\t",C[i][j]);
        }
        printf("\n");
    }
    
}
void ScatterA()
{
    int imin,imax,jmin,jmax;
    int sp;
    int i,j,k,m;
    for(k=0; k<p; k++)
    {
       
        sp = (int)sqrt(p);
        imin = (k / sp) * np;
        imax = imin + np - 1;
        jmin = (k % sp) * np;
        jmax = jmin + np -1;
       
        m = 0;
        for(i=imin; i<=imax; i++)
        {
            for(j=jmin; j<=jmax; j++)
            {
                tempa[m] = A[i][j];
                 tempb[m] = B[i][j];
              //  tempb[m] = B[j][i];
                m++;
            }
        }
        
       
        if(k==0)
        {
            memcpy(a, tempa, np*np*sizeof(int));
            memcpy(b, tempb, np*np*sizeof(int));
        }
        else
        {
            MPI_Send(tempa, np*np, MPI_INT, k, 1, MPI_COMM_WORLD);
            MPI_Send(tempb, np*np, MPI_INT, k, 2, MPI_COMM_WORLD);
        }
    }
}
void MainProcess() //cannon
{
    MPI_Comm comm;
    int crank;
    int dims[2],periods[2], coords[2];
    int source, dest, up, down, right, left;
    int i;
    dims[0] = dims[1] = (int)sqrt(p);
    periods[0] = periods[1] = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm);
    MPI_Comm_rank(comm, &crank);
    MPI_Cart_coords(comm, crank, 2, coords);
    MPI_Cart_shift(comm, 1, -1, &right, &left);
    MPI_Cart_shift(comm, 0, -1, &down, &up);
    MPI_Cart_shift(comm, 1, -coords[0], &source, &dest);
    MPI_Sendrecv_replace(a, np*np, MPI_INT, dest, 1, source, 1, comm, &status);
    MPI_Cart_shift(comm, 0, -coords[1], &source, &dest);
    MPI_Sendrecv_replace(b, np*np, MPI_INT, dest, 1, source, 1, comm, &status);
   // printf("rank1 My rank= %d\n",rank);
    Mutiply();
    for(i = 1; i < dims[0]; i++)
    {
        // Mutiply();
        MPI_Sendrecv_replace(a, np*np, MPI_INT, left, 1, right, 1, comm, &status);
        MPI_Sendrecv_replace(b, np*np, MPI_INT, up, 1, down, 1, comm, &status);
     //   printf("rank2 My rank= %d\n",rank);
        Mutiply();
        //printf("rank = %d\n", rank);
        //   PrintAC();
    }
    MPI_Comm_free(&comm);
}
void collectC()
{
    int i,j,k,s,m;
    int imin,imax,jmin,jmax;
    int sp= (int)sqrt(p);
    for (i=0;i<np;i++)
    {
        for(j=0;j<np;j++)
            C[i][j]=c[i*np+j];
    }
    for (k=1;k<p;k++)
    {
       
        MPI_Recv(c, np*np, MPI_INT, k, 1, MPI_COMM_WORLD, &status);
        //printf("rank = %d\n", rank);PrintAC();
        imin = (k / sp) * np;
        imax = imin + np - 1;
        jmin = (k % sp) * np;
        jmax = jmin + np -1;
        
        for(i=imin,m=0; i<=imax; i++,m++)
        {
            for(j=jmin,s=0; j<=jmax; j++,s++)
                C[i][j]=c[m*np+s];
        }
    }
}
void Mutiply()
{
    int i,j,k;
    for(i=0; i<np; i++)
    {
        for(j=0; j<np; j++)
        {
            for(k=0; k<np; k++)
            {
                c[i*np+j] += a[i*np+k]*b[k*np+j];
              //  printf("c =%d ,a,b=(%d，%d) =(%d ,%d)\n",c[i*np+j],i*np+k,k*np+j,a[i*np+k],b[k*np+j]);
            }
        }
    }
}

void matrix_multiplication()
{
    int i, j, k;
    // int **arr3;
    int ret;
    printf("mult matrix:\n");
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ret = 0;
            for (k = 0; k < n; k++)
            {
                ret += A[i][k] * B[k][j];
            }
            //arr3[i][j] = ret;
            printf("%d\t",ret);
        }
        printf("\n");
    }
}

