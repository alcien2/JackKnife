#ifndef ITERACTIONS_H_INCLUDED
#define ITERACTIONS_H_INCLUDED

void MakeArrayGreatAgain (SuperArrayDim4<Link<MatrixSU2>> & m){
    double arr[4] = {1, 0, 0, 0};
    for (int i = 0; i < n; i ++){
        for (int j = 0; j < n; j ++){
                for (int k = 0; k < n; k ++){
                    for (int l = 0; l < n; l++){
                            for(int t = 0; t < 4; t ++){
                                m[i][j][k][l].get(t).set_matrix(arr);
                            }
                    }
                }
        }
    }
   /* double b = 0;
    for(int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                    for (int k = 0; k < 2; k ++){
                           for (int q = 0; q < 2; q++){
                                    for (int l = 0; l < 4; l++){
                                            for(int t = 0; t < 4; t++){
                                                cout << m[i][j][k][q].get(l).get(t) << endl;
                                                b += (m[i][j][k][q].get(l).get(t))*(m[i][j][k][q].get(l).get(t));
                                            }
                                            cout<<"sum " << b << endl;
                                            b = 0;
                                    }
                            }
            }
        }
    }*/
}

void kcounter (SuperArrayDim4<Link<MatrixSU2>> & arr, int i, int j, int k, int m, int l, double & coef, MatrixSU2 & Matt){
    MatrixSU2 M[7];
    Matrix matrix ;
    double f[4] = {cos(phi), 0, 0, sin(phi)};
    MatrixSU2 phase;
    phase.set_matrix(f);
    switch(i){
        case n - 1:
        if (l == 0){
            M[0] = arr[i + 1][j][k][m].get(1)* phase *arr[i][j + 1][k][m].get(0).Inverse()*arr[i][j][k][m].get(1).Inverse();
        //for(int v = 0; v < 4; v ++){ffout << "M[0]" << M[0].get(v) <<endl;};
            M[1] = arr[i + 1][j][k][m].get(2) * arr[i][j][k + 1][m].get(0).Inverse()*arr[i][j][k][m].get(2).Inverse();
            M[2] = arr[i + 1][j][k][m].get(3) * arr[i][j][k][m + 1].get(0).Inverse()*arr[i][j][k][m].get(3).Inverse();
            M[3] = phase.Inverse() * arr[i + 1][j - 1][k][m].get(1).Inverse() * arr[i][j - 1][k][m].get(0).Inverse() * arr[i][j - 1][k][m].get(1);
            M[4] = arr[i + 1][j][k - 1][m].get(2).Inverse() *  arr[i][j][k - 1][m].get(0).Inverse() * arr[i][j][k - 1][m].get(2);
            M[5] = arr[i + 1][j][k][m - 1].get(3).Inverse() *  arr[i][j][k][m - 1].get(0).Inverse() * arr[i][j][k][m - 1].get(3);
        }
        if(l == 1){
            M[0] = arr[i][j + 1][k][m].get(0)* phase.Inverse() * arr[i + 1][j][k][m].get(1).Inverse() * arr[i][j][k][m].get(0).Inverse();
            M[1] = arr[i][j + 1][k][m].get(2)*arr[i][j][k + 1][m].get(1).Inverse()*arr[i][j][k][m].get(2).Inverse();
            M[2] = arr[i][j + 1][k][m].get(3)*arr[i][j][k][m + 1].get(1).Inverse()*arr[i][j][k][m].get(3).Inverse();
            M[3] = arr[i - 1][j + 1][k][m].get(0).Inverse() * arr[i - 1][j][k][m].get(1).Inverse() * arr[i - 1][j][k][m].get(0);
            M[4] = arr[i][j + 1][k - 1][m].get(2).Inverse() * arr[i][j][k - 1][m].get(1).Inverse() * arr[i][j][k - 1][m].get(2);
            M[5] = arr[i][j + 1][k][m - 1].get(3).Inverse() * arr[i][j][k][m - 1].get(1).Inverse() * arr[i][j][k][m - 1].get(3);
        }
        if(l == 2){
            M[0] = arr[i][j][k + 1][m].get(0)*arr[i + 1][j][k][m].get(2).Inverse()*arr[i][j][k][m].get(0).Inverse();
            M[1] = arr[i][j][k + 1][m].get(1)*arr[i][j + 1][k][m].get(2).Inverse()*arr[i][j][k][m].get(1).Inverse();
            M[2] = arr[i][j][k + 1][m].get(3)*arr[i][j][k][m + 1].get(2).Inverse()*arr[i][j][k][m].get(3).Inverse();
            M[3] = arr[i - 1][j][k + 1][m].get(0).Inverse() * arr[i - 1][j][k][m].get(2).Inverse() * arr[i - 1][j][k][m].get(0);
            M[4] = arr[i][j - 1][k + 1][m].get(1).Inverse() * arr[i][j - 1][k][m].get(2).Inverse() * arr[i][j - 1][k][m].get(1);
            M[5] = arr[i][j][k + 1][m - 1].get(3).Inverse() * arr[i][j][k][m - 1].get(2).Inverse() * arr[i][j][k][m - 1].get(3);
        }
        if(l == 3){
            M[0] = arr[i][j][k][m + 1].get(0)*arr[i + 1][j][k][m].get(3).Inverse()*arr[i][j][k][m].get(0).Inverse();
            M[1] = arr[i][j][k][m + 1].get(1)*arr[i][j + 1][k][m].get(3).Inverse()*arr[i][j][k][m].get(1).Inverse();
            M[2] = arr[i][j][k][m + 1].get(2)*arr[i][j][k + 1][m].get(3).Inverse()*arr[i][j][k][m].get(2).Inverse();
            M[3] = arr[i - 1][j][k][m + 1].get(0).Inverse() * arr[i - 1][j][k][m].get(3).Inverse() * arr[i - 1][j][k][m].get(0);
            M[4] = arr[i][j - 1][k][m + 1].get(1).Inverse() * arr[i][j - 1][k][m].get(3).Inverse() * arr[i][j - 1][k][m].get(1);
            M[5] = arr[i][j][k - 1][m + 1].get(2).Inverse() * arr[i][j][k - 1][m].get(3).Inverse() * arr[i][j][k - 1][m].get(2);
        }
        M[6] = M[0] + M[1] + M[2] + M[3] + M[4] + M[5];

        M[6].groupElement(matrix);
        break;

        case 0:
        if (l == 0){
            M[0] = arr[i + 1][j][k][m].get(1) *arr[i][j + 1][k][m].get(0).Inverse()*arr[i][j][k][m].get(1).Inverse();
            M[1] = arr[i + 1][j][k][m].get(2)* arr[i][j][k + 1][m].get(0).Inverse()*arr[i][j][k][m].get(2).Inverse();
            M[2] = arr[i + 1][j][k][m].get(3)* arr[i][j][k][m + 1].get(0).Inverse()*arr[i][j][k][m].get(3).Inverse();
            M[3] = arr[i + 1][j - 1][k][m].get(1).Inverse() * arr[i][j - 1][k][m].get(0).Inverse() * arr[i][j - 1][k][m].get(1);
            M[4] = arr[i + 1][j][k - 1][m].get(2).Inverse() * arr[i][j][k - 1][m].get(0).Inverse() * arr[i][j][k - 1][m].get(2);
            M[5] = arr[i + 1][j][k][m - 1].get(3).Inverse() * arr[i][j][k][m - 1].get(0).Inverse() * arr[i][j][k][m - 1].get(3);
        }
        if(l == 1){
            M[0] = arr[i][j + 1][k][m].get(0)*arr[i + 1][j][k][m].get(1).Inverse() *  arr[i][j][k][m].get(0).Inverse();
            M[1] = arr[i][j + 1][k][m].get(2)*arr[i][j][k + 1][m].get(1).Inverse()*arr[i][j][k][m].get(2).Inverse();
            M[2] = arr[i][j + 1][k][m].get(3)*arr[i][j][k][m + 1].get(1).Inverse()*arr[i][j][k][m].get(3).Inverse();
            M[3] = phase * arr[i - 1][j + 1][k][m].get(0).Inverse() * arr[i - 1][j][k][m].get(1).Inverse() * arr[i - 1][j][k][m].get(0);
            M[4] = arr[i][j + 1][k - 1][m].get(2).Inverse() * arr[i][j][k - 1][m].get(1).Inverse() * arr[i][j][k - 1][m].get(2);
            M[5] = arr[i][j + 1][k][m - 1].get(3).Inverse() * arr[i][j][k][m - 1].get(1).Inverse() * arr[i][j][k][m - 1].get(3);
        }
        if(l == 2){
            M[0] = arr[i][j][k + 1][m].get(0)*arr[i + 1][j][k][m].get(2).Inverse()* arr[i][j][k][m].get(0).Inverse();
            M[1] = arr[i][j][k + 1][m].get(1)*arr[i][j + 1][k][m].get(2).Inverse()*arr[i][j][k][m].get(1).Inverse();
            M[2] = arr[i][j][k + 1][m].get(3)*arr[i][j][k][m + 1].get(2).Inverse()*arr[i][j][k][m].get(3).Inverse();
            M[3] = arr[i - 1][j][k + 1][m].get(0).Inverse() * arr[i - 1][j][k][m].get(2).Inverse() * arr[i - 1][j][k][m].get(0);
            M[4] = arr[i][j - 1][k + 1][m].get(1).Inverse() * arr[i][j - 1][k][m].get(2).Inverse() * arr[i][j - 1][k][m].get(1);
            M[5] = arr[i][j][k + 1][m - 1].get(3).Inverse() * arr[i][j][k][m - 1].get(2).Inverse() * arr[i][j][k][m - 1].get(3);
        }
        if(l == 3){
            M[0] = arr[i][j][k][m + 1].get(0)*arr[i + 1][j][k][m].get(3).Inverse() * arr[i][j][k][m].get(0).Inverse();
            M[1] = arr[i][j][k][m + 1].get(1)*arr[i][j + 1][k][m].get(3).Inverse()*arr[i][j][k][m].get(1).Inverse();
            M[2] = arr[i][j][k][m + 1].get(2)*arr[i][j][k + 1][m].get(3).Inverse()*arr[i][j][k][m].get(2).Inverse();
            M[3] = arr[i - 1][j][k][m + 1].get(0).Inverse() * arr[i - 1][j][k][m].get(3).Inverse() * arr[i - 1][j][k][m].get(0);
            M[4] = arr[i][j - 1][k][m + 1].get(1).Inverse() * arr[i][j - 1][k][m].get(3).Inverse() * arr[i][j - 1][k][m].get(1);
            M[5] = arr[i][j][k - 1][m + 1].get(2).Inverse() * arr[i][j][k - 1][m].get(3).Inverse() * arr[i][j][k - 1][m].get(2);
        }
        M[6] = M[0] + M[1] + M[2] + M[3] + M[4] + M[5];

        M[6].groupElement(matrix);
        break;

        default:
        if (l == 0){
        M[0] = arr[i + 1][j][k][m].get(1) *arr[i][j + 1][k][m].get(0).Inverse()*arr[i][j][k][m].get(1).Inverse();
        M[1] = arr[i + 1][j][k][m].get(2)* arr[i][j][k + 1][m].get(0).Inverse()*arr[i][j][k][m].get(2).Inverse();
        M[2] = arr[i + 1][j][k][m].get(3) * arr[i][j][k][m + 1].get(0).Inverse()*arr[i][j][k][m].get(3).Inverse();
        M[3] = arr[i + 1][j - 1][k][m].get(1).Inverse()  * arr[i][j - 1][k][m].get(0).Inverse() * arr[i][j - 1][k][m].get(1);
        M[4] = arr[i + 1][j][k - 1][m].get(2).Inverse()  * arr[i][j][k - 1][m].get(0).Inverse() * arr[i][j][k - 1][m].get(2);
        M[5] = arr[i + 1][j][k][m - 1].get(3).Inverse()  * arr[i][j][k][m - 1].get(0).Inverse() * arr[i][j][k][m - 1].get(3);
    }
    if(l == 1){
        M[0] = arr[i][j + 1][k][m].get(0)*arr[i + 1][j][k][m].get(1).Inverse()*arr[i][j][k][m].get(0).Inverse();
        M[1] = arr[i][j + 1][k][m].get(2)*arr[i][j][k + 1][m].get(1).Inverse()*arr[i][j][k][m].get(2).Inverse();
        M[2] = arr[i][j + 1][k][m].get(3)*arr[i][j][k][m + 1].get(1).Inverse()*arr[i][j][k][m].get(3).Inverse();
        M[3] = arr[i - 1][j + 1][k][m].get(0).Inverse() * arr[i - 1][j][k][m].get(1).Inverse()  * arr[i - 1][j][k][m].get(0);
        M[4] = arr[i][j + 1][k - 1][m].get(2).Inverse() * arr[i][j][k - 1][m].get(1).Inverse() * arr[i][j][k - 1][m].get(2);
        M[5] = arr[i][j + 1][k][m - 1].get(3).Inverse() * arr[i][j][k][m - 1].get(1).Inverse() * arr[i][j][k][m - 1].get(3);
    }
    if(l == 2){
        M[0] = arr[i][j][k + 1][m].get(0)*arr[i + 1][j][k][m].get(2).Inverse()* arr[i][j][k][m].get(0).Inverse();
        M[1] = arr[i][j][k + 1][m].get(1)*arr[i][j + 1][k][m].get(2).Inverse()*arr[i][j][k][m].get(1).Inverse();
        M[2] = arr[i][j][k + 1][m].get(3)*arr[i][j][k][m + 1].get(2).Inverse()*arr[i][j][k][m].get(3).Inverse();
        M[3] = arr[i - 1][j][k + 1][m].get(0).Inverse() * arr[i - 1][j][k][m].get(2).Inverse() * arr[i - 1][j][k][m].get(0);
        M[4] = arr[i][j - 1][k + 1][m].get(1).Inverse() * arr[i][j - 1][k][m].get(2).Inverse() * arr[i][j - 1][k][m].get(1);
        M[5] = arr[i][j][k + 1][m - 1].get(3).Inverse() * arr[i][j][k][m - 1].get(2).Inverse() * arr[i][j][k][m - 1].get(3);
    }
    if(l == 3){
        M[0] = arr[i][j][k][m + 1].get(0)*arr[i + 1][j][k][m].get(3).Inverse() * arr[i][j][k][m].get(0).Inverse();
        M[1] = arr[i][j][k][m + 1].get(1)*arr[i][j + 1][k][m].get(3).Inverse()*arr[i][j][k][m].get(1).Inverse();
        M[2] = arr[i][j][k][m + 1].get(2)*arr[i][j][k + 1][m].get(3).Inverse()*arr[i][j][k][m].get(2).Inverse();
        M[3] = arr[i - 1][j][k][m + 1].get(0).Inverse() * arr[i - 1][j][k][m].get(3).Inverse() * arr[i - 1][j][k][m].get(0);
        M[4] = arr[i][j - 1][k][m + 1].get(1).Inverse() * arr[i][j - 1][k][m].get(3).Inverse() * arr[i][j - 1][k][m].get(1);
        M[5] = arr[i][j][k - 1][m + 1].get(2).Inverse() * arr[i][j][k - 1][m].get(3).Inverse() * arr[i][j][k - 1][m].get(2);
    }
    M[6] = M[0] + M[1] + M[2] + M[3] + M[4] + M[5];
    M[6].groupElement(matrix);
    }

    complex<double>coeff = sqrt(matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]);
    coef = coeff.real();
    for(int q = 0; q < 4; q ++){
        Matt.get(q) = M[6].get(q) / coef;
    }

}

void probability1(SuperArrayDim4<Link<MatrixSU2>>& arr, int i, int j, int k, int m, int l, double & beta, TestClass & randomm, double & coef){
    double probability;
	double a = exp(-(coef * beta));
	double b = exp(coef * beta);
    double z = randomm.returnRandom(a, b);
    probability = sqrt(1 - (log(z) * log(z)) / (beta * beta * coef * coef ));
    double prob = randomm.returnRandom(0, 1);
    while (prob > probability){
        z = randomm.returnRandom(a, b);
        probability = sqrt(1 - ( log(z) * log(z)) / (beta * beta * coef * coef ) );
        prob = randomm.returnRandom(0, 1);
    }
        arr[i][j][k][m].get(l).get(0) = log(z) / coef / beta ;

}

void probability2(SuperArrayDim4<Link<MatrixSU2>>& arr, int i, int j, int k, int m, int l, TestClass & randomm)
{
    double g[3];
    double sum;
    while(true){
        sum = 0;
        for(int q = 0; q < 3; q ++){
                g[q] = randomm.returnRandom(-1,1);
                sum += g[q]*g[q];
        }
        if(sum > 0 || sum <= 1){
            break;
            }
    }
    double sum1 = sqrt(sum);
    for(int q = 0; q < 3; q ++){
        g[q] = g[q] / sum1;
    }
    for(int q = 0; q < 3; q ++){
        arr[i][j][k][m].get(l).get(q + 1) = g[q]*sqrt(1 - arr[i][j][k][m].get(l).get(0)*arr[i][j][k][m].get(l).get(0));
    }
}

void iteraction(int iter, SuperArrayDim4<Link<MatrixSU2>> &arr, double & beta, TestClass & randomm)
{
    for (int v = 0; v < iter; v++)
    {
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                for(int k = 0; k < n; k++)
                {
                    for(int m = 0; m < n; m++)
                    {
                        for(int l = 0; l < 4; ++l)
                            {
                                MatrixSU2 M;
                                double coef = 0;
                                kcounter(arr, i, j, k, m, l, coef, M);
                                probability1(arr, i, j, k, m, l, beta, randomm, coef);
                                probability2(arr, i, j, k, m, l, randomm);
                                arr[i][j][k][m].get(l) = arr[i][j][k][m].get(l) * M.Inverse();
                                //double z = randomm.returnRandom(0,1);
                                //cout << z << endl;

                        }
                    }
                }
            }
        }
    }
}

#endif // ITERACTIONS_H_INCLUDED
