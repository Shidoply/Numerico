#include <iostream>
#include <fstream>
using namespace std;

#define PADRON1 94950
#define PADRON2 93410

int cargarDatos(const char *archivo, double *&x, double *&y){
	int n;
	ifstream in(archivo);
	if(!in.is_open()) return -1;
	
	//Primera línea: cantidad de puntos
	in >> n;
	x = new double[n];
	y = new double[n];
	
	for(int i=0;i<n;i++){
		in >> x[i] >> y[i]; //Cada línea tiene un par (x,y)
		y[i] *= (PADRON1+PADRON2)/(2.0*100000.0);
	}
	
	return n;
}

void prepararMatriz(double *a, double *h, int n){
	for(int y=0;y<n;y++)
	for(int x=0;x<n;x++){
		int valor = 0;
		
		if(x == y){
			valor = 2.0*(h[y]+h[y+1]); 	//Diagonal
		}else if(x == y+1){
			valor = h[y+1]; 			//Arriba de la diagonal
		}else if(x == y-1){
			valor = h[y];				//Abajo de la diagonal
		}
		a[x+y*n] = valor;
	}
}

//Resuelve Ax = b por SOR, Siendo n la dimensión
//w la constante de relajación (?) y rtol el error mínimo deseado
int sor(double *a, double *x, double *b, int n, double w, double rtol){
	int it=0;
	while(it<100){ ///Por hacer: reemplazar por cálculo de error en cada iteración
		for(int j=0;j<n;j++){
			double suma = 0.0;
			
			for(int k=0;k<n;k++)
				if(j!=k)
					suma += a[j*n+k] * x[k];
			
			//Supongo que la diagonal no es cero
			x[j] = w*(b[j]-suma)/a[j*n+j] + (1.0-w)*x[j];
		}
		it++;
	}
	return it;
}

//Dados los C_k, h_k los pares de puntos y los h, 
//imprime todos los polinomios con sus intervalos
void poly(double *c, double *x, double *y, double *h, int n){
	for(int k=0;k<n;k++){
		double a = y[k];
		double b = (y[k+1]-y[k])/h[k] - (h[k]/3.0) * (2.0*c[k]+c[k+1]); ///Arreglo temporal
		double d = (c[k+1]-c[k])/(3.0*h[k]); 
		
		cout << "y = "
			 <<	a << "+" 
			 << b << "*(x-"<<x[k]<<")+" 
			 << c[k] << "*(x-"<<x[k]<<")^2+" 
			 << d << "*(x-"<<x[k]<<")^3"
			 << " {" << x[k] << ";" << x[k+1] << "}" << endl;
			 
		//Al graphmatica no le gusta la notación 2e5, hay que reemplazarlo por 2*10^5
	}
}

int main(int argc, char **args){
	double *x=NULL, *y=NULL;
	int n = cargarDatos("d1.csv", x, y)-1;
	if(n<0){
		cerr << "No se pudo abrir el archivo." << endl;
		return 1;
	}
	
	//N es la cantidad de polinomios, no de puntos totales!
	cout << "N = " << n << endl;
	
	double *h = new double[n];
	for(int k=0;k<n;k++)
		h[k] = x[k+1] - x[k];
	
	//Creo el sistema a resolver
	double *aSist = new double[(n-1)*(n-1)];
	double *xSist = new double[n+1];
	double *bSist = new double[n-1];
	
	//Semilla inicial para el SOR
	for(int i=0;i<=n;i++)
		xSist[i] = 0;
	
	prepararMatriz(aSist, h, n-1);
	
	//Preparo el B del sistema
	for(int k=0;k<n-1;k++)
		bSist[k] = (3.0/h[k+1]) * (y[k+2]-y[k+1]) - (3.0/h[k+1]) * (y[k+1]-y[k]);
	
	//Resuelvo el sistema, teniendo en cuenta que xSist son los C_k, por lo que 'salteo' el primer elemento
	sor(aSist,xSist+1,bSist,n-1,1.0,0.001);
	poly(xSist, x, y, h, n);
	
	delete []aSist;
	delete []xSist;
	delete []bSist;
	delete []x;
	delete []y;
	return 0;
}
