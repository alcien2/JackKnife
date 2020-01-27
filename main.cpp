#include <random>
//#include <ctime>
#include <chrono>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
//#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <complex>
#include <new>
#include <iomanip>
#include <vector>
#include <fstream>
double phi;
double beta;
double alpha;
#include "Class_SU2.h"
#include "Arrays.h"
#include "Iteractions.h"
#include "ToGet.h"
#include "functions.h"
#include "ToSave.h"
#include "JackKnife.h"



int main(int argc, char** argv)
{
	phi = atof(argv[1]);
	beta = atof(argv[2]);
	alpha = atof(argv[3]);
    TestClass randomm;
    //srand(static_cast<unsigned int>(time(0)));
    srand (time(NULL));
    ofstream fout;
    string outfile = GetFileNameResult();
    fout.open(outfile);
    /*ofstream ffout, fout, fout1;
    ffout.open("/home/itep/petrova/JackKnife/H_=_0res2.txt");
    fout.open("/home/itep/petrova/JackKnife/H_=_0res2.txt");
    fout1.open("/home/itep/petrova/JackKnife/H_=_0res2.txt");*/
    //cout << randomm.returnRandom(0,1) << endl;
    //cout << randomm.returnRandom(0,1) << endl;
    SuperArrayDim4<Link<MatrixSU2>> *ar;
    try{
        ar = new SuperArrayDim4<Link<MatrixSU2>>(n);
    }catch(bad_alloc xa){
        cout << "yps" << endl;
        return 1;
    }
    //fout << phi << " " << beta << " " << alpha << endl;
    //double beta = 2.8;
    //MakeArrayGreatAgain(*ar);
    JackKnife(*ar, randomm, fout);
    //iteraction(15, *ar, beta, randomm);

    //Saver(*ar, beta, randomm, ffout); //Это для создания файлов
   // Sigma(*ar, randomm, fout, ffout, fout1); // Посчитать сигму

   /*
    double beta = 0.01;
    for(int q = 0; q < 16; q ++){
    MakeArrayGreatAgain(*ar);
    iteraction(iter, *ar, beta, randomm);
    int i,j;
    double res = 0;
    double error;
    ForWilsonloop(4,4, *ar, beta, randomm, res, error);
    cout << res << " " << error << endl;
    for(i = 1; i < 6; i ++){
        j = i;

    }
    double error = 0;
    //double loop = Average(*ar);
    double loop = Wilsonloop(*ar, error, beta, randomm);
    //double loop = Wilsonloop(ar);
    fout << "{{" << beta << ","<< loop << "}, ErrorBar[" << error << "]},";
    cout << loop <<" " << error << endl;
    beta += 0.25;
    //iter ++;
    }*/
    /*fout1.close();
    fout2.close();
    fout3.close();
    fout4.close();
    ffout.close();*/
    fout.close();
    return 0;
}
