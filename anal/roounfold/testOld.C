void testOld()
{
    constexpr int M = 2;
    constexpr int N = 2;
    double A[M][N];
    A[0][0]=0.6;A[0][1]=0;
    A[1][0]=0;A[1][1]=0.4;
    std::cout << "Response matrix" << std::endl;
    for(int i=0;i<M;i++){
    for(int j=0;j<N;j++)std::cout << A[i][j] << " ";
    std::cout << std::endl;
    }
    double m[N];
    m[0] = 19.2;
    m[1] = 15.6;
    // Bayes
    double dk[N],dk1[N];
    dk[0] = 1; dk[1] =2;
    int Niter=20;
    for(int iter=0;iter<Niter;iter++){
        //cout << iter << " iter " ;
        double denom[N];
        for(int i = 0; i < N; i++){
            double sum = 0;
            for(int j = 0; j < M; j++) {
                sum += A[i][j]*dk[j];
            }
            denom[i] = sum;
        }
        for(int j=0;j<M;j++){
            double sum=0.;
            for(int i=0;i<N;i++){
                //double sum2=0.;
                //for(int l=0;l<N;l++) sum2=sum2+A[i][l]*dk[l];
                sum=sum+A[i][j]*m[i]/denom[i];
            }
            double eps=0.;
            for(int i=0;i<N;i++)eps=eps+A[i][j];
            dk1[j]=dk[j]*sum/eps;
        }
        //cout << endl;
        for(int i=0;i<N;i++) {
            std::cout << dk[i] << " ";
            dk[i]=dk1[i];
        }
        std::cout << std::endl;
    }

}
