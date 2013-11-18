#include <iostream>
#include <cstdio>
#include <cmath>
using namespace std;

const double G=9.81;
const double A=M_PI*(7.0*100.0+37.0)/1000.0;
const double L=1.0;
const double HINT = 0.00001;

double f1(double theta, double omega, double t){
	return omega;
}
double f2(double theta, double omega, double t){
	return -G/L*sin(theta);
}
double energia(double theta, double omega){
	return 0.5*omega*L*omega*L + G*L*(1-cos(theta));
}

double euler(double k, int n, double theta, double omega){
	for(int i=0;i<n;i++){
		theta = theta + k*f1(theta,omega,k*i);		
		omega = omega + k*f2(theta,omega,k*i);
	}
	return energia(theta,omega);
}
double rk2(double k, int n, double theta, double omega){
	for(int i=0;i<n;i++){
		double q1_theta = k*f1(theta,omega,k*i);
		double q1_omega = k*f2(theta,omega,k*i);
		double q2_theta = k*f1(theta+q1_theta,omega+q1_omega,k*(i+1.0));
		double q2_omega = k*f2(theta+q1_theta,omega+q1_omega,k*(i+1.0));
		theta = theta + (q1_theta+q2_theta)*0.5;		
		omega = omega + (q1_omega+q2_omega)*0.5;
	}
	return energia(theta,omega);
}
double rk4(double k, int n, double theta, double omega){
	for(int i=0;i<n;i++){
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
	}
	return energia(theta,omega);
}

double calcularPeriodo(double amplitud){
	double acc = 0;
	//Rectángulo
	for(double ix=0;ix<=M_PI/2.0;ix+=HINT){
		acc += pow(1-pow(sin(amplitud/2.0),2.0)*pow(sin(ix),2.0),-0.5)*HINT;
	}
	//y=(1-(sin(737/2000*pi)^2)*(sin(x)^2))^(-0.5)
	//y''max = -13
	cout << 13/24.0*(M_PI/2.0)*HINT*HINT << endl;
	cout << acc << endl; 
	return 4*sqrt(L/G)*acc;
}

int main(){
	double theta_0 = A;
	double omega_0 = 0;
	double periodoIntegrado = calcularPeriodo(A);
	
	cout.precision(5);
	cout << "Periodo aprox: " << 2*M_PI*sqrt(L/G) << endl;
	cout << "Periodo integrado: " << periodoIntegrado << endl;
	cout << "Periodo para A=0.01: " << calcularPeriodo(0.01) << endl;
	cout << "Energia inicial del pendulo: " << energia(theta_0,omega_0) << endl;
	cout << "Energia final de Euler:" << euler(0.3, 10.0*periodoIntegrado, theta_0, omega_0) << endl;
	cout << "Energia final de RK2:" << rk2(0.3, 10.0*periodoIntegrado, theta_0, omega_0) << endl;
	cout << "Energia final de RK4:" << rk4(0.3, 10.0*periodoIntegrado, theta_0, omega_0) << endl;
	
	return 0;
}