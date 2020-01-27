#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

double Wilsonloop(SuperArrayDim4<Link<MatrixSU2>> & arr, double & error, double & beta, TestClass &randomm){
    double loop[10] = {0,0,0,0,0,0,0,0,0,0};
    double Loop = 0;
    int t = 0;
    for(int q = 0; q < 10; q ++){
    for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                for(int k = 0; k < n; k++)
                {
                    for(int m = 0; m < n; m++)
                    {
                        loop[q] += (arr[i][j][k][m].get(0)*arr[i + 1][j][k][m].get(1)*arr[i][j + 1][k][m].get(0).Inverse()*arr[i][j][k][m].get(1).Inverse()).trace()+
                        (arr[i][j][k][m].get(0)*arr[i + 1][j][k][m].get(2)*arr[i][j][k + 1][m].get(0).Inverse()*arr[i][j][k][m].get(2).Inverse()).trace() +
                        (arr[i][j][k][m].get(0)*arr[i + 1][j][k][m].get(3)*arr[i][j][k][m + 1].get(0).Inverse()*arr[i][j][k][m].get(3).Inverse()).trace() +
                        (arr[i][j][k][m].get(1)*arr[i][j + 1][k][m].get(2)*arr[i][j][k + 1][m].get(1).Inverse()*arr[i][j][k][m].get(2).Inverse()).trace() +
                        (arr[i][j][k][m].get(1)*arr[i][j + 1][k][m].get(3)*arr[i][j][k][m + 1].get(1).Inverse()*arr[i][j][k][m].get(3).Inverse()).trace() +
                        (arr[i][j][k][m].get(2)*arr[i][j][k + 1][m].get(3)*arr[i][j][k][m + 1].get(2).Inverse()*arr[i][j][k][m].get(3).Inverse()).trace();
                    }
                }
            }
        }
        loop[q] = loop[q] /6/ n / n/ n / n;
        iteraction(2, arr, beta, randomm);
    }

    for(int q = 0; q < 10; q ++){
        Loop +=loop[q];
    }
    Loop = Loop/10;

        for(t = 0; t < 10; t ++){
            error += (Loop - loop[t])*(Loop - loop[t]);
        }
        error = sqrt(error / 10/9);
        error = 0.5*error;
        Loop = 1 - 0.5*Loop;

        return Loop;
}

template <typename MatrixSU>
double WilsonLoop(SuperArrayDim4<Link<MatrixSU2>> & m, int I, int J) {
    double sum = 0;
    double b = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                for (int l = 0; l < n; l++) {
                    MatrixSU U(0);

                    /*if((i + I < n - 1) &&(j + J < n - 1)){
                        b ++;
                    for (int x = 0; x < I; x++) {
                        U *= m[i+x][j][k][l].get(0);
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i+I][j+y][k][l].get(1);
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i+I-1-x][j+J][k][l].get(0).Inverse();
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j+J-1-y][k][l].get(1).Inverse();
                    }
                    sum += U.trace();
                    U = MatrixSU(0); }

                    if((i + J < n - 1) && (j + I < n - 1)){
                        b ++;
                    for (int y = 0; y < J; y++) {
                        U *= m[i+y][j][k][l].get(0);
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i+J][j+x][k][l].get(1);
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i+J-1-y][j+I][k][l].get(0).Inverse();
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j+I-1-x][k][l].get(1).Inverse();
                    }
                    sum += U.trace();
                    U = MatrixSU(0);}*/




                    /*if((i + I < n - 1) && (k + J < n - 1)){
                        b++;
                    for (int x = 0; x < I; x++) {
                        U *= m[i+x][j][k][l].get(0);
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i+I][j][k+y][l].get(2);
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i+I-1-x][j][k+J][l].get(0).Inverse();
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j][k+J-1-y][l].get(2).Inverse();
                    }
                    sum += U.trace();
                    U = MatrixSU(0);}

                    if((i + J < n - 1) && (k + I < n - 1)){
                        b++;
                    for (int y = 0; y < J; y++) {
                        U *= m[i+y][j][k][l].get(0);
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i+J][j][k+x][l].get(2);
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i+J-1-y][j][k+I][l].get(0).Inverse();
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j][k+I-1-x][l].get(2).Inverse();
                    }
                    sum += U.trace();
                    U = MatrixSU(0);}*/




                   /* if((i + I < n - 1) && (l + J < n - 1)){
                        b++;
                    for (int x = 0; x < I; x++) {
                        U *= m[i+x][j][k][l].get(0);
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i+I][j][k][l+y].get(3);
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i+I-1-x][j][k][l+J].get(0).Inverse();
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j][k][l+J-1-y].get(3).Inverse();
                    }
                    sum += U.trace();
                    U = MatrixSU(0);}

                    if((i + J < n - 1) && (l + I < n - 1)){
                        b++;
                    for (int y = 0; y < J; y++) {
                        U *= m[i+y][j][k][l].get(0);
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i+J][j][k][l+x].get(3);
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i+J-1-y][j][k][l+I].get(0).Inverse();
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j][k][l+I-1-x].get(3).Inverse();
                    }
                    sum += U.trace();
                    U = MatrixSU(0);}*/




                   /* if((j + I < n -1) && (k + J < n - 1)){
                        b++;
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j+x][k][l].get(1);
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j+I][k+y][l].get(2);
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j+I-1-x][k+J][l].get(1).Inverse();
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j][k+J-1-y][l].get(2).Inverse();
                    }
                    sum += U.trace();
                    U = MatrixSU(0);}

                    if((j + J < n - 1) && (k + I < n - 1)){
                        b ++;
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j+y][k][l].get(1);
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j+J][k+x][l].get(2);
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j+J-1-y][k+I][l].get(1).Inverse();
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j][k+I-1-x][l].get(2).Inverse();
                    }
                    sum += U.trace();
                    U = MatrixSU(0);}*/



                    /*if((j + I < n - 1) && (l + J < n - 1)){
                        b++;
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j+x][k][l].get(1);
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j+I][k][l+y].get(3);
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j+I-1-x][k][l+J].get(1).Inverse();
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j][k][l+J-1-y].get(3).Inverse();
                    }
                    sum += U.trace();
                    U = MatrixSU(0);}

                    if((j + J < n - 1) && (l + I < n - 1)){
                        b++;
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j+y][k][l].get(1);
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j+J][k][l+x].get(3);
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j+J-1-y][k][l+I].get(1).Inverse();
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j][k][l+I-1-x].get(3).Inverse();
                    }
                    sum += U.trace();
                    U = MatrixSU(0);}*/




                    if(((k + I) < (n - 1)) && ((l + J) < (n - 1))){
                        b++;
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j][k+x][l].get(2);
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j][k+I][l+y].get(3);
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j][k+I-1-x][l+J].get(2).Inverse();
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j][k][l+J-1-y].get(3).Inverse();
                    }
                    sum += U.trace();
                    U = MatrixSU(0);}

                    if(((k + J) < (n - 1))&& ((l + I) < (n - 1))){
                        b++;
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j][k+y][l].get(2);
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j][k+J][l+x].get(3);
                    }
                    for (int y = 0; y < J; y++) {
                        U *= m[i][j][k+J-1-y][l+I].get(2).Inverse();
                    }
                    for (int x = 0; x < I; x++) {
                        U *= m[i][j][k][l+I-1-x].get(3).Inverse();
                    }
                    sum += U.trace();
                    U = MatrixSU(0);}
                }
            }
        }
    }
    sum /= b;
    //sum /= 2*n*n*n*n;
    return sum;
}

double Average(SuperArrayDim4<Link<MatrixSU2>> & arr){
    return 1 - WilsonLoop<MatrixSU2>(arr, 1, 1)/2;
}

void ForWilsonloop (int i, int j, SuperArrayDim4<Link<MatrixSU2>> & arr, double & beta, TestClass &randomm, double & res, double &Error){
    int t = 50;
    double a[t + 1];
    double err = 0;
    for(int q = 0; q < t; q ++){
        iteraction(2, arr, beta, randomm);
        a[q] = WilsonLoop<MatrixSU2>(arr, i, j);
    }
    a[t + 1] = 0;
    for(int q = 0; q < t; q ++){
        a[t] += a[q];
    }
    a[t] = a[t]/t;
    for(int q = 0; q < t; q ++){
        err +=(a[t] - a[q])*(a[t] - a[q]);
    }
    err = err/t/(t - 1);
    res = a[t];
    Error = err;
}

void ReturnChi(int i, int j, SuperArrayDim4<Link<MatrixSU2>> & arr, double & beta, TestClass &randomm, double & res, double &Error){
    int t = ErrorNum;
    double a[t + 1], b[t + 1], c[t + 1], d[t + 1];
    for(int q = 0; q < t; q ++){
    //iteraction(5, arr, beta, randomm);
   // SuperArrayDim4<Link<MatrixSU2>>  arr;
    GetArray(arr, q);
    a[q] = WilsonLoop<MatrixSU2>(arr, i, j);
    b[q] = WilsonLoop<MatrixSU2>(arr, i - 1, j - 1);
    c[q] = WilsonLoop<MatrixSU2>(arr, i, j - 1);
    d[q] = WilsonLoop<MatrixSU2>(arr, i - 1, j);
    }
    a[t] = 0;
    b[t] = 0;
    c[t] = 0;
    d[t] = 0;
    for(int q = 0; q < t; q ++){
        a[t] += a[q];
        b[t] += b[q];
        c[t] += c[q];
        d[t] += d[q];
    }
    a[t] = a[t]/t;
    b[t] = b[t] /t;
    c[t] = c[t] /t;
    d[t] = d[t] /t;
    double erra = 0;
    double errb = 0;
    double errc = 0;
    double errd = 0;
    for(int q = 0; q < t; q ++){
        erra += (a[t] - a[q])*(a[t] - a[q]);
        errb += (b[t] - b[q])*(b[t] - b[q]);
        errc += (c[t] - c[q])*(c[t] - c[q]);
        errd += (d[t] - d[q])*(d[t] - d[q]);
    }
    erra = sqrt(erra/t/(t - 1));
    errb = sqrt(errb/t/(t - 1));
    errc = sqrt(errc/t/(t - 1));
    errd = sqrt(errd/t/(t - 1));
    res = -log(a[t]*b[t]/c[t]/d[t]);
    // res = -log(WilsonLoop<MatrixSU2>(arr, i, j)*WilsonLoop<MatrixSU2>(arr, i - 1, j - 1)/WilsonLoop<MatrixSU2>(arr, i - 1, j)/WilsonLoop<MatrixSU2>(arr, i, j- 1));
     cout << res << " "<< beta << endl;
    //Error = sqrt(1/a[t]/a[t]/log(a[t]*b[t]/c[t]/d[t])/log(a[t]*b[t]/c[t]/d[t]) * erra*erra + 1/b[t]/b[t]/log(a[t]*b[t]/c[t]/d[t])/log(a[t]*b[t]/c[t]/d[t]) *errb*errb
    //             + 1/c[t]/c[t]/log(a[t]*b[t]/c[t]/d[t])/log(a[t]*b[t]/c[t]/d[t])*errc*errc + 1/d[t]/d[t]/log(a[t]*b[t]/c[t]/d[t])/log(a[t]*b[t]/c[t]/d[t])*errd*errd);
    Error = sqrt(1/a[t]/a[t] * erra*erra + 1/b[t]/b[t] *errb*errb + 1/c[t]/c[t]*errc*errc + 1/d[t]/d[t]*errd*errd);
    //return ;
}

void Chi(SuperArrayDim4<Link<MatrixSU2>> & arr, TestClass & randomm, ofstream & fout1, ofstream & fout2, ofstream & fout3, ofstream & fout4){
    int j;
    int i;
    int iter = 10;
    double beta = 0.01;
    double res ;
    double Error;
    for(int q = 0; q < 20; q ++){
            //MakeArrayGreatAgain(arr);
            //iteraction(iter, arr, beta, randomm);
            res = 0;
            Error = 0;
            /*ReturnChi(1,1, arr, beta, randomm, res, Error);
            fout1 << "{{" << beta << "," << log(res) << "}, ErrorBar[" << Error << "]},";
            ReturnChi(2,2, arr, beta, randomm, res, Error);
            fout2 << "{{" << beta << "," << log(res) << "}, ErrorBar[" << Error << "]},";*/
            ReturnChi(3,3, arr, beta, randomm, res, Error);
            fout3 << "{{" << beta << "," << log(res) << "}, ErrorBar[" << Error << "]},";
            ReturnChi(4,4, arr, beta, randomm, res, Error);
            fout4 << "{{" << beta << "," << log(res) << "}, ErrorBar[" << Error << "]},";

            beta +=0.2;
   }

}
    void Sigma (SuperArrayDim4<Link<MatrixSU2>> & arr, TestClass & randomm, ofstream & fout, ofstream & foutt, ofstream & out){
        double res = 0;
        double Error = 0;
        double beta = 2.8;
        //MakeArrayGreatAgain(arr);
       // iteraction(15, arr, beta, randomm);
        for(int i = 1; i < 23; i ++){
            for(int j = 1; j < 23; j ++){
                    if((i*j > 14) && (i * j < 36)&& (i < 10) && (j < 10)&& (i > 3)&& (j > 3)){
                        cout << "yoy"<< endl;
                        ReturnChi(i,j, arr, beta, randomm, res, Error);
                        fout << "{{" << i*j << "," << res << "}, ErrorBar[" << Error << "]},";
                        foutt << "{" << i*j << "," << res << "},";
                        out << res <<",";
                        cout << res << endl;}
            }
        }
    }

    void Average (double a[24], ofstream & fout){
        double average = 0;
        double error = 0;
        for(int i = 0; i < 24; i ++){
            average +=a[i];
        }
        average = average / 24;
        for(int i = 0; i < 24; i ++){
            error +=(average - a[i])*(average - a[i]);
        }
        error = sqrt(error/24/23);
        fout << average << ", " << error << endl;
    }

#endif // FUNCTIONS_H_INCLUDED
