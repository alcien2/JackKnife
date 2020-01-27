#ifndef CLASS_SU2_H_INCLUDED
#define CLASS_SU2_H_INCLUDED

#define PI 3.14159265

using namespace std ;
const int n = 14;
const int ErrorNum = 500;
//const double beta = 1;
   using Matrix = complex<double>[2][2];
   using TMatrix = double[n][n][n][n][4][4];

    Matrix sigma0 = {
        {{1.0, 0.0}, {0.0, 0.0}},
        {{0.0, 0.0}, {1.0, 0.0}}
    };
    Matrix sigma1 = {
        {{0.0, 0.0}, {1.0, 0.0}},
        {{1.0, 0.0}, {0.0, 0.0}}
    };
    Matrix sigma2 = {
        {{0.0, 0.0}, {0.0, -1.0}},
        {{0.0, 1.0}, {0.0, 0.0}}
    };
    Matrix sigma3 = {
        {{1.0, 0.0}, {0.0, 0.0}},
        {{0.0, 0.0}, {-1.0, 0.0}}
        };


class MatrixSU2 {
public:
    MatrixSU2() {};
    MatrixSU2(int i);

    double trace();
    double det();
    MatrixSU2 Inverse();
    void set_matrix(double b[4]);
    double& get(int i);
    void groupElement(Matrix U);

    double& operator[] (int k);
    MatrixSU2 operator*= (MatrixSU2 B);
    MatrixSU2 operator/= (MatrixSU2 B);
    friend MatrixSU2 operator+(MatrixSU2 A, MatrixSU2 B);
    friend MatrixSU2 operator*(MatrixSU2 A, MatrixSU2 B);
    friend MatrixSU2 operator/(MatrixSU2 A, MatrixSU2 B);

private:
    double a[4];
    };

MatrixSU2::MatrixSU2(int i){
        double b[4] = {1,0,0,0};
        if (i == 0) {
            (*this).set_matrix(b);
        }
        if (i == 1){
            double g[4];
            double sum;
           // while(true){
                sum = 0;
                for(int j = 0; j < 4; j ++){
                    g[j] = (double)rand()/RAND_MAX;
                    sum += g[j]*g[j];
                }
              /*  if(sum == 0 || sum > 1){
                    continue;
                }*/
           // }
            double sum1 = sqrt(sum);
            for(int j = 0; j < 4; j ++){
                a[j] = g[j]/sum1;
                if(rand()%2){
                    a[j] = -a[j];
                }
            }

        }
}


double& MatrixSU2::get(int i){
        return a[i];
    }

void MatrixSU2::set_matrix(double b[4]){
        for(int i = 0; i < 4; i++){
            a[i] = b[i];
        }
    }

double MatrixSU2:: trace(){
    return 2*a[0];
}

double MatrixSU2::det(){

}

double& MatrixSU2::operator[] (int k) {
    if(k < 0 || k > 3){
        cout << "Выход за пределы"  << endl;
        exit(1);
    }
    return a[k];
}

MatrixSU2 MatrixSU2:: Inverse(){
        MatrixSU2 inversed;
        inversed[0] = a[0];
        inversed[1] = -a[1];
        inversed[2] = -a[2];
        inversed[3] = -a[3];
    return inversed;
}

void MatrixSU2::groupElement(Matrix U){
    complex<double> b = {0.0, 1.0};
    for ( int r = 0; r < 2; ++r){
        for(int v = 0; v < 2; ++v){
            U[r][v] = {0.0, 0.0};
            U[r][v] += a[0] * sigma0[r][v] + a[1] * b * sigma1[r][v]
            + a[2] * b * sigma2[r][v] + a[3]* b * sigma3[r][v];
        }
    }
}

/*MatrixSU2::MatrixSU2() {
    return;
}*/

MatrixSU2 operator*(MatrixSU2 A, MatrixSU2 B) {
    MatrixSU2 C;

    C[0] = A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    C[1] = A[0]*B[1] + B[0]*A[1] - A[2]*B[3] + A[3]*B[2];
    C[2] = A[0]*B[2] + B[0]*A[2] - A[3]*B[1] + A[1]*B[3];
    C[3] = A[0]*B[3] + B[0]*A[3] - A[1]*B[2] + A[2]*B[1];

    return C;
}

MatrixSU2 operator+(MatrixSU2 A, MatrixSU2 B){
    MatrixSU2 C;
    C[0] = A[0] + B[0];
    C[1] = A[1] + B[1];
    C[2] = A[2] + B[2];
    C[3] = A[3] + B[3];
    return C;
}

MatrixSU2 operator/(MatrixSU2 A, MatrixSU2 B) {
    return A*B.Inverse();
}

MatrixSU2 MatrixSU2::operator*=(MatrixSU2 B) {
    MatrixSU2 temp = (*this)*B;
    *this = temp;

    return temp;
}

MatrixSU2 MatrixSU2::operator/=(MatrixSU2 B) {
    MatrixSU2 temp = (*this)/B;
    *this = temp;

    return *this;
}


class TestClass {
public:
  TestClass() {
        mt19937 mersenne(static_cast<unsigned int>(time(NULL)));
    //std::random_device device;
    random_generator_.seed(mersenne());  }
    double returnRandom(){
        uniform_real_distribution<> range(0, RAND_MAX);
        return range(random_generator_);
    }
  double returnRandom(double min, double max) {
    uniform_real_distribution<> range(min, max);
    return range(random_generator_);
  }

private:
  mt19937 random_generator_;

};

string GetFileName (int num ){
    stringstream filename_constructor;
    filename_constructor << "/home/itep/petrova/Iteraction/files_beta_" << beta  << "_H_" << phi  <<"/file_" << num << ".txt";
    return filename_constructor.str();
};

string GetFileNameResult (){
    stringstream filename_constructor;
    filename_constructor << "/home/itep/petrova/JackKnife/Results_H" << phi << "/Result_alpha" << alpha << ".txt";
    return filename_constructor.str();
};

#endif // CLASS_SU2_H_INCLUDED
