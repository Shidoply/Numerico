#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
using namespace std;

const double G=9.81;
const double A=M_PI*(7.0*100.0+37.0)/1000.0;
const double L=1.0;
ofstream fout("salida.csv");

double f1(double theta, double omega, double t){
	return omega;
}
double f2(double theta, double omega, double t){
	return -G/L*sin(theta);
}
double energia(double theta, double omega){
	return 0.5*omega*L*omega*L + G*L*(1-cos(theta));
}

double euler(double k, int n, int pasosPorPeriodo, double periodoIntegrado, double theta, double omega){
	int signoOmega = 1, per=0;

	for(int i=0;i<n;i++){
		if(i%pasosPorPeriodo == 0) cout << "\tPeriodo " << i/pasosPorPeriodo << " - Amplitud:" << theta << " ERel%:" << (theta-A)/A *100 << endl;
				
		theta = theta + k*f1(theta,omega,k*i);
		omega = omega + k*f2(theta,omega,k*i);
		fout << i*k << "\t" << theta << endl;
		
		if(omega<0 && signoOmega == 1){
			//Paso de positivo a negativo
			cout << "\tMaxima amplitud para T=" << k*i << " Esperado: " << periodoIntegrado*per << " ERel%:" << 100*(k*i-periodoIntegrado*per)/(k*i) << endl;
			per++;
		}
		signoOmega = omega>0?1:0;

	}
	return energia(theta,omega);
}
double rk2(double k, int n, int pasosPorPeriodo, double periodoIntegrado, double theta, double omega){
	int signoOmega = 1, per=0;
	for(int i=0;i<n;i++){
		if(i%pasosPorPeriodo == 0) cout << "\tPeriodo " << i/pasosPorPeriodo << " - Amplitud:" << theta << " ERel%:" << (theta-A)/A *100 << endl;
	
		double q1_theta = k*f1(theta,omega,k*i);
		double q1_omega = k*f2(theta,omega,k*i);
		double q2_theta = k*f1(theta+q1_theta,omega+q1_omega,k*(i+1.0));
		double q2_omega = k*f2(theta+q1_theta,omega+q1_omega,k*(i+1.0));
		theta = theta + (q1_theta+q2_theta)*0.5;
		omega = omega + (q1_omega+q2_omega)*0.5;
		
		fout << i*k << "\t" << theta << endl;
		
		if(omega<0 && signoOmega == 1){
			//Paso de positivo a negativo
			cout << "\tMaxima amplitud para T=" << k*i << " Esperado: " << periodoIntegrado*per << " ERel%:" << 100*(k*i-periodoIntegrado*per)/(k*i) << endl;
			per++;
		}
		signoOmega = omega>0?1:0;
	}
	return energia(theta,omega);
}
double rk4(double k, int n, int pasosPorPeriodo, double periodoIntegrado, double theta, double omega){
	int signoOmega = 1, per=0;
	for(int i=0;i<n;i++){
		if(i%pasosPorPeriodo == 0) cout << "\tPeriodo " << i/pasosPorPeriodo << " - Amplitud:" << theta << " ERel%:" << (theta-A)/A *100 << endl;

		double q1_theta = k*f1(theta,omega,k*i);
		double q1_omega = k*f2(theta,omega,k*i);
		double q2_theta = k*f1(theta+q1_theta/2.0,omega+q1_omega/2.0,k*(i+0.5));
		double q2_omega = k*f2(theta+q1_theta/2.0,omega+q1_omega/2.0,k*(i+0.5));
		double q3_theta = k*f1(theta+q2_theta/2.0,omega+q2_omega/2.0,k*(i+0.5));
		double q3_omega = k*f2(theta+q2_theta/2.0,omega+q2_omega/2.0,k*(i+0.5));
		double q4_theta = k*f1(theta+q3_theta,omega+q3_omega,k*(i+1.0));
		double q4_omega = k*f2(theta+q3_theta,omega+q3_omega,k*(i+1.0));
		theta = theta + (q1_theta+2.0*q2_theta+2.0*q3_theta+q4_theta)/6.0;
		omega = omega + (q1_omega+2.0*q2_omega+2.0*q3_omega+q4_omega)/6.0;
		
		fout << i*k << "\t" << theta << endl;
		
		if(omega<0 && signoOmega == 1){
			//Paso de positivo a negativo
			cout << "\tMaxima amplitud para T=" << k*i << " Esperado: " << periodoIntegrado*per << " ERel%:" << 100*(k*i-periodoIntegrado*per)/(k*i) << endl;
			per++;
		}
		signoOmega = omega>0?1:0;		
	}
	return energia(theta,omega);
}

double integrando(double x, double amplitud){
	return 4*sqrt(L/G)*pow(1-pow(sin(amplitud/2.0),2.0)*pow(sin(x),2.0),-0.5);
}

double calcularPeriodo(double amplitud){
	//Trapecio compuesto
	double acc = 0;
	const int n = 350;
	const double b = M_PI/2.0;
	const double HINT = b/n;
	const double fddMax= 17.0; //Cota de la derivada segunda entre 0 y PI/2

	acc += integrando(0,amplitud) * HINT * 0.5;
	acc += integrando(b,amplitud) * HINT * 0.5;

	for(int i=1;i<n;i++)
			acc += integrando(i*HINT,amplitud) * HINT;

	cout << "Error estimado: " << b*HINT*HINT*fddMax/12 << endl;
	return acc;
}

int main(){
	double theta_0 = A;
	double omega_0 = 0;
	double periodoIntegrado = calcularPeriodo(A);
	int pasosPorPeriodo=75;
	double paso = periodoIntegrado/pasosPorPeriodo;
	int cantPasos = pasosPorPeriodo*10;

	cout.precision(5);
	cout << "Periodo aprox: " << 2*M_PI*sqrt(L/G) << endl;
	cout << "Periodo integrado: " << periodoIntegrado << endl;
	cout << "Periodo para A=0.01: " << calcularPeriodo(0.01) << endl;

	cout.precision(8);
	cout << endl << "Amplitud inicial: " << A << endl;
	cout << "Energia inicial del pendulo: " << energia(theta_0,omega_0) << endl;
	cout << "Euler" << endl;
	cout << "\tEnergia final Euler:" << euler(paso, cantPasos, pasosPorPeriodo, periodoIntegrado, theta_0, omega_0) << endl;
	cout << "RK2" << endl;
	cout << "\tEnergia final RK2:" << rk2(paso, cantPasos, pasosPorPeriodo, periodoIntegrado, theta_0, omega_0) << endl;
	cout << "RK4" << endl;
	cout << "\tEnergia final RK4:" << rk4(paso, cantPasos, pasosPorPeriodo, periodoIntegrado, theta_0, omega_0) << endl;

	return 0;
}
