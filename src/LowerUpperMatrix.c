/*
 * LowerUpperMatrix.c
 *
 *  Created on: 13/02/2012
 *      Author: John Haiber Osorio, Lina Maria Perez Perez
 */

/* El siguiente codigo sirve para dar solucion a la ley de enfriamiento de newton definida de la siguiente forma
 *
 *     dT/dt = -K (T(t) - Ta)
 *
 *     Ta		= Temperatura ambiente
 *     T0 		= Temperatura Inicial
 *     T  		= Temperatura Final
 *     K  		= Coeficiente de conductividad
 *     DeltaX 	= Resolucion
 *     t		= Tiempo final
 *
 *
 *
 */

// Inclusion de librerias
#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/*-----------------------------------------------------------------------------------------------
 * Funcion para la solucion de un sistema de ecuaciones lineales de la forma A=LU
 *
 * donde:
 *
 * a	 = matriz que tiene la descomposicion LU
 * b 	 = Vector que tiene almacenado los terminos independientes de cada ecuacion, b ya ha sido permutado de
 * 			acuerdo al vector de permutacion pivot
 * pivot = Vector que contiene las permutaciones realizadas en la funcion que calcula la descomposicion LU
 * n	 = numero de ecuaciones a solucionar
 * x 	 = Vector respuesta
 *
 *
 * La solucion al sistema puede ser encontrado teniendo en cuenta lo siguiente
 *
 * primero se soluciona:
 * 							L*Y = B			(4) Se hace por forward substitution
 *
 * y despues solucionando
 * 							U*x = Y			(5) Se hace por back substitution
 *
 *
*----------------------------------------------------------------------------------------------
*/

int solveSystemEquationsLU(float *aModificada,float *bModificado, int *pivot,int n, float *x){

	int i,j;
	float *y;
	float sumatoria;

	y = (float*)malloc(sizeof(float)*n);

	for (i = 0; i < n; i++) {
		sumatoria = 0.0;										//Forward sustitution ecuacion (4)
		for(j = 0; j <= i - 1 ;j++)
			sumatoria = sumatoria + (aModificada[i*n+j]*y[j]);
		y[i] = bModificado[i] - sumatoria;
	}


	for (i = n-1 ; i >= 0 ; i--){
		sumatoria = 0.0;
		for(j = i+1 ; j < n ;j++)
			sumatoria = sumatoria + (aModificada[i*n+j]*x[j]);//back sustitution ecuacion (5)
		x[i] = (1/aModificada[i*n+i])*(y[i] - sumatoria);
	}

	return 0;
}


/*---------------------------------------------------------------------------------------------------------
 * Funcion que realiza la permutacion de un vector teniendo en cuenta el vector pivot.
 *
 * b 		= vector a permutar
 * pivot	= vector de permutacion
 * n		= tamaño del vector
 * --------------------------------------------------------------------------------------------------------
 */

int permutarVector(float *b, int *pivot, int n){
	int indice;
	float temp;
	for(indice = 0;indice < n; indice++){
		if(indice != pivot[indice]){
			temp = b[indice];
			b[indice] = b[pivot[indice]];
			b[pivot[indice]] = temp;
		}
	}
	return 0;

}

/*-----------------------------------------------------------------------------------------------------
 * Funcion que permite imprimir un vector
 * ----------------------------------------------------------------------------------------------------
 */

int imprimirVector(int *vector,int n,char *mensaje){
	int i;
	printf("\n %s: \n",mensaje);
		for(i = 0; i < n; i++){
			printf("%d  ",vector[i]);
		}
}

/*-----------------------------------------------------------------------------------------------------
 * Funcion que permite imprimir una matriz
 * ----------------------------------------------------------------------------------------------------
 */

int imprimirMatrix(float *matrix, int n, int m, char *mensaje){
	int i,j;
	printf("\n %s: \n",mensaje);
	for(i = 0; i < n; i++){
		for(j = 0;j < m ; j++)
			printf("%f  ",matrix[i*n+j]);
		printf("\n");
	}
	return 0;
}

/*--------------------------------------------------------------------------------------
 *Funcion que permite realizar la descomposicion LU de una matriz n*n
 *
 *El siguiente codigo realiza la solucion de uns sitema de ecuaciones lineales en la notacion
 *
 * 									Ax=b
 * utilizando el metodo de descomposicion LU con pivote parcial y haciendo uso de un vector de permutacion
 *
 * donde se tiene que
 * 									PA=LU.
 *
 *
 * a 	 = matriz de tamaño n*n a la cual se le va a realizar la descomposicion
 * n 	 = tamaño de la matriz (directamente relacionado con el numero de ecuaciones)
 * pivot = matriz que contendra las permutaciones realizadas por el algoritmo para
 * 			posteriormente ser utilizadas para dar la solucion al sistema de ecuaciones
 *
 * ------------------------------------------------------------------------------------
 */

int LUDecomposition(float *a, int n, int *pivot){

	float sumatoria = 0.0, big = 0.0, temp = 0.0, temp1 = 0.0;
	int i=0,j=0,k=0,rowMax = -1;
	float *permutationMatrix,*lMatrix,*uMatrix;
	permutationMatrix = (float *) malloc(sizeof(float)*n*n);
	lMatrix = (float *) malloc(sizeof(float)*n*n);
	uMatrix = (float *) malloc(sizeof(float)*n*n);


	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			if(i == j){
				permutationMatrix[i*n+j] = 1.0;
				lMatrix[i*n+j] = 1.0;
			}
			else{
				permutationMatrix[i*n+j] = 0.0;
				lMatrix[i*n+j] = 0.0;
			}
		}
	}


	for(j = 0; j < n; j++){

		pivot[j] = j;
		big = 0.0;
		rowMax = j;
		for(i = j;i < n; i++){
			k=j;
			if(fabs(a[i*n+k]) > big){
				big = fabs(a[i*n+k]);
				rowMax = i ;
			}
		}

///////////////////// Pivote medio//////////////////////////////////////////////////////////////
		if(rowMax != j){
			for(i=0;i<n;i++){
				temp = a[rowMax*n+i];
				a[rowMax*n+i] = a[j*n+i];
				a[j*n+i] = temp;

				temp1 = permutationMatrix[rowMax*n+i];
				permutationMatrix[rowMax*n+i] = permutationMatrix[j*n+i];
				permutationMatrix[j*n+i] = temp1;
			}
		}
		pivot[j] = rowMax;

//////////////////////// Calculo Betha ///////////////////////////////////////////////////////

		for(i = j; i < n; i++){
			sumatoria = 0.0;
			for(k = 0; k <= j - 1 ;k++)
				sumatoria = sumatoria + (a[j*n+k]*a[k*n+i]);
			a[j*n+i] = a[j*n+i] - sumatoria;
		}

////////////////////////Calculo Alpha ////////////////////////////////////////////////////////
		for ( i = j+1 ; i < n ; i++ ){
			sumatoria = 0.0;
			for(k = 0; k <= j - 1 ;k++)
				sumatoria = sumatoria + (a[i*n+k]*a[k*n+j]);
			a[i*n+j] = (1/a[j*n+j])*(a[i*n+j] - sumatoria);
		}

	}

	return 0;

}


/*------------------------------------------------------------------------------------------------------------------------
 * Funcion principal del programa
 * -----------------------------------------------------------------------------------------------------------------------
 */


int main(){



	float temperaturaInicial, temperaturaFinal, temperaturaAmbiente, tiempoFinal, tiempoInicial, coeficienteConductividad, deltaX;
	float *a, *b, *x;
	int *pivot;
	int numeroEcuaciones=0,i,j;


	temperaturaInicial = 80;
	temperaturaFinal = 28.12011;
	temperaturaAmbiente = 20;
	deltaX = 0.25;
	tiempoInicial = 0.0;
	tiempoFinal = 2;
	coeficienteConductividad = 1;

	//Reserva de memoria para la matriz y los vectores del programa

	numeroEcuaciones= (int)((tiempoFinal - tiempoInicial)/deltaX) - 1;
	a = (float *)malloc(sizeof(float)*numeroEcuaciones*numeroEcuaciones);
	b = (float *)malloc(sizeof(float)*numeroEcuaciones);
	x = (float *)malloc(sizeof(float)*numeroEcuaciones);
	pivot = (int *)malloc(sizeof(int)*numeroEcuaciones);

	printf("El numero de Ecuaciones es de : %d",numeroEcuaciones);

	for (i = 0; i < numeroEcuaciones*numeroEcuaciones; i++){
		a[i] = 0.0;
	}

	// LLenado de los vectores y matrices para la solucion de la ecuacion de enfriamiento de newton

	for (i = 0; i < numeroEcuaciones; i++) {
		b[i] = 2*coeficienteConductividad*deltaX*temperaturaAmbiente;
		pivot[i] = 0;
		x[i] = 0.0;
		for(j = 0; j < numeroEcuaciones; j++){
			if (i == j){
				a[i*numeroEcuaciones+j] = 2*coeficienteConductividad*deltaX;
				if (i < numeroEcuaciones-1){
				    a[(i+1)*numeroEcuaciones+j] = -1;
				    a[(i*numeroEcuaciones)+(j+1)] = 1;
				}
			}
		}
	}

	b[0] = 2*coeficienteConductividad*deltaX*temperaturaAmbiente + temperaturaInicial;
	b[numeroEcuaciones-1] = 2*coeficienteConductividad*deltaX*temperaturaAmbiente - temperaturaFinal;

	imprimirMatrix(a,numeroEcuaciones,numeroEcuaciones,"La matriz armada es");

	LUDecomposition(a,numeroEcuaciones,pivot);

	imprimirMatrix(a,numeroEcuaciones,numeroEcuaciones,"La matriz A conteniendo la descomposición LU es");
	imprimirVector(pivot,numeroEcuaciones,"El vector pivote es");
	imprimirMatrix(b,1,numeroEcuaciones,"El vector b es");
	permutarVector(b,pivot,numeroEcuaciones);
	imprimirMatrix(b,1,numeroEcuaciones,"El vector b' es");
	solveSystemEquationsLU(a,b,pivot,numeroEcuaciones,x);
	imprimirMatrix(x,1,numeroEcuaciones,"El vector x es");

	return 0;
}
