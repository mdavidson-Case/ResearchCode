#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <random>

using namespace std;
double random_numbers_test(){
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution <double> dist(0.0, 100.0);
    for(int i = 0; i < 100; i++){
        cout << dist(mt) << ", ";
    }
    return 0;
}

double Ising(double T, int N, double J1){
    //int N = 12;
    //double J1 = 0.50;
    double J2 = 1;
    double E = 0;
    double M = 0;
    double E1 = 0;
    double M1 = 0;
    double M2 = 0;
    double M4 = 0;
    double X = 0;
    double B = 0;
    double eq_sweep = pow(10,5);
    double mc_sweep = pow(10,4);
    double Beta = 1/T;
    double Beta2 = Beta*Beta;
    double n1 = 1/(mc_sweep*N*N);
    double n2 = 1/(mc_sweep*mc_sweep*N*N);
    vector <double> energy_vector;
    vector <double> mag_vector;
    vector <double> mag2;
    vector <double> mag4;
    
    //instantiate RNG
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution <double> dist(0.0, 100.0);
    
    //fill lattice
    double spin[N][N];
    for(int i = 0; i <= N-1; i++){
        for (int j = 0; j <= N-1; j++){
            spin[i][j] = -1;
        }
    }

    //begin equilibration runs
    for(int b=0; b<=eq_sweep; b ++){
        for(int i = 0; i <= N-1; i++){
            for (int j = 0; j <= N-1; j++){
                double cost = 0;
                double si = spin[i][j];
                double nb = (spin[(i+1)%(N)][j] + spin[i][(j+1)%(N)]+ spin[(i-1)%(N)][j] + spin[i][(j-1)%(N)]);
                double r = 0;
                
                //cost += 2*si*nb;
                r+= dist(mt);

                if(j < N/2 - 1){
                    cost = 2*(J2+(2/N)*j*(J1-J2))*si*nb;

                    if(cost < 0){
                        spin[i][j] = -spin[i][j];
                    }
                    else{
                        if(r <= 100*exp(-cost/T)){
                            spin[i][j] = -spin[i][j];
                        }
                        else{
                            spin[i][j] = spin[i][j];
                        }
                    }
                }
                else{
                    cost = 2*((2/N)*(J2-J1)*j-J2+2*J1)*si*nb;

                     if(cost < 0){
                        spin[i][j] = -spin[i][j];
                    }
                    else{
                        if(r <= 100*exp(-cost/T)){
                            spin[i][j] = -spin[i][j];
                        }
                        else{
                            spin[i][j] = spin[i][j];
                        }
                    }

                }

                /*if(b >= eq_sweep -2){
                    cout << cost << ", ";
                }*/
            }
        }
    }
    //end of equilibration runs

    //begin data collection
    for(int b = 0; b <= mc_sweep; b++){
        double count = 0 ;
        for(int i = 0; i <= N-1; i++){
            for (int j = 0; j <= N-1; j++){

                double cost = 0;
                double si = spin[i][j];
                double nb = (spin[(i+1)%(N)][j] + spin[i][(j+1)%(N)]+ spin[(i-1)%(N)][j] + spin[i][(j-1)%(N)]);
                double r = 0;

                //cost += 2*si*nb;
                r+= dist(mt);
                
                if(j < N/2 - 1){
                    cost = 2*(J2+(2/N)*j*(J1-J2))*si*nb;

                    if(cost < 0){
                        spin[i][j] = -spin[i][j];
                    }
                    else{
                        if(r <= 100*exp(-cost/T)){
                            spin[i][j] = -spin[i][j];
                        }
                        else{
                            spin[i][j] = spin[i][j];
                        }
                    }
                }
                else{
                    cost = 2*((2/N)*(J2-J1)*j-J2+2*J1)*si*nb;

                     if(cost < 0){
                        spin[i][j] = -spin[i][j];
                    }
                    else{
                        if(r <= 100*exp(-cost/T)){
                            spin[i][j] = -spin[i][j];
                        }
                        else{
                            spin[i][j] = spin[i][j];
                        }
                    }

                }

                count  += spin[i][j];

                energy_vector.push_back(-spin[i][j]*nb/2);
                mag_vector.push_back(spin[i][j]);
            }
        }
        mag2.push_back(count*count);
        mag4.push_back(count*count*count*count);

        //cout << "count: " << count << endl;
        /*double Msum = 0;
        for(int i =0; i <+ mag_vector.size(); i++){
            mag_vector[i]
        }*/
    }

    for(int i= 0; i <= energy_vector.size(); i++){
        E1 += energy_vector[i];
    }
    for(int i= 0; i <= energy_vector.size(); i++){
        M1 += mag_vector[i];
    }
    for(int i = 0; i<=mag2.size(); i++){
        M2 += mag2[i];
        M4 += mag4[i];
    }
    
    E1 = E1*n1;
    M = abs(M1*n1);
    
    X = (n1*M2 - n2*M1*M1)*Beta;
    B = 1 - n1*M4/(3*n2*M2*M2); 
    //for some reason J1 = 0.5 keeps returning values in reverse.
    //This might have to do with how the lattice is communicating with itself. 
    //could check to see if a value of J1 = 0 would return the same thing as a J1 = 1.
    //this would indicate the lattice is communicating primarily through the left and right boundary. 
   
   // cout << M1 << ", " << M2 << endl;
    return E1;
}



int main() {

    cout << Ising(1,9,0.5) << endl;
   cout << Ising(1.5,9,0.5) << endl;
    cout << Ising(2,9,0.5) << endl;
    cout << Ising(10,9,0.5) << endl;

    for(double J1 = 1; J1 >= 0.59; J1-= 0.10){
        cout << "j2 = 1, J1 = " << J1 << endl;
        for(int N = 8; N <= 40; N+=8){
            cout << " N: " <<  N << endl;
            for(double i = 0.01; i <= 4.0; i+=0.2){
                cout << "(" << i << ", " <<  Ising(i, N, J1) <<  "), " << endl;
            }
        }
    }
    return 0;
}