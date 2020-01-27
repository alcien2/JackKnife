#ifndef JACKKNIFE_H_INCLUDED
#define JACKKNIFE_H_INCLUDED

void ReturnChi2(int i, int j, SuperArrayDim4<Link<MatrixSU2>> & arr, TestClass &randomm, double & res){
	int t = ErrorNum;
	int count = 0;
	double a[t + 1], b[t + 1], c[t + 1], d[t + 1];
	for(int q = 0; q < t; q ++){
		if((q >=  ErrorNum - alpha*50) && (q < ErrorNum -(alpha-1)*50)){ continue;}
            GetArray(arr, q);
		count ++;
            a[q] = WilsonLoop<MatrixSU2>(arr, i, j); 
	//cout << a[q] << " a[q]" << i << " i" << j << " j" << q << " q" <<  endl; 
            b[q] = WilsonLoop<MatrixSU2>(arr, i - 1, j - 1);
            c[q] = WilsonLoop<MatrixSU2>(arr, i, j - 1);
            d[q] = WilsonLoop<MatrixSU2>(arr, i - 1, j);
	}
	cout << count << " count configurations" << endl;
	a[t] = 0;
	b[t] = 0;
	c[t] = 0;
        d[t] = 0;
	for(int q = 0; q < t; q ++){
		if((q >=  ErrorNum - alpha*50) && (q < ErrorNum -(alpha-1)*50)){ continue;}
		a[t] += a[q];
		b[t] += b[q];
		c[t] += c[q];
		d[t] += d[q];
	}
	a[t] = a[t] /count;
	b[t] = b[t] /count;
	c[t] = c[t] /count;
	d[t] = d[t] /count;
	/*double erra = 0;
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
	errd = sqrt(errd/t/(t - 1));*/
	cout << a[t] << b[t] << c[t] << d[t] << "a,b,c,d" << endl;
	res = -log(a[t]*b[t]/c[t]/d[t]);
}

void JackKnife(SuperArrayDim4<Link<MatrixSU2>> & arr, TestClass & randomm, ofstream & fout){
	double Result;
	double points[100];
	//int alpha = 1;
    int count = 0;
	//for(alpha = 1; alpha < 11; alpha ++){
            for(int i = 1; i < 23; i++){
                for(int j = 1; j < 23; j ++){
                    if((i*j > 14) && (i * j < 36)&&(j <10)&&(i < 10)&&(i >= 4)&& (j >= 4)){
                        ReturnChi2(i, j, arr, randomm, points[count]);
                        count++;
                    }
                }
            }
	points[count + 1] = 0;
    for(int i = 0; i < count; i++){
        cout << points[i] << endl;    
	points[count + 1] += points[i];
    }
    Result = points[count + 1]/count;
    fout << Result << endl;
	}
#endif
