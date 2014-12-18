#include <iostream>
#include <iomanip>
#include <cstring> 
#include <Math.h>
#include <cmath>
#include <set>
#include <algorithm>
#include <fstream>
#include "Bmp.h" 
const double PI = 3.141592;
using namespace std;

Image::Bmp bmp;
int dataSize = 512 * 512 * 8 / 8;
unsigned char *ori = new unsigned char[dataSize];
unsigned char *comp = new unsigned char[dataSize];
int matInt[8][8];
int DCTMatrix[8][8];
int quantI[8][8];
int zigzag[8][8];
int izigzag[8][8];
int IquantI[8][8];
int IDCTMatrix[8][8];

int probData5[64 * 64][10];
int histoEntropy[64 * 64][10];
int entropylimit = 5;
int counterData = 0;

bool quantizebool = true;
int quality,re;
char filename[20];

void calcEntropy(int mode){
	for (int a = 0; a < mode; a++){
		for (int i = 0; i < 64 * 64; i++){
			if (probData5[i][a] >= 0){
				histoEntropy[probData5[i][a]][a] += 1;
			}
			else{
				histoEntropy[probData5[i][a] + 64 * 64-1][a] += 1;
			}
		}
	}
	
	double entropy=0.0;
	double sumEntropy = 0.0;
	double check = 0.0;
	for (int a = 0; a < mode; a++){
		for (int i = 0; i < 64 * 64; i++){
			if (histoEntropy[i][a] != 0){
				check += histoEntropy[i][a] / (64.0 * 64.0);
				entropy += (histoEntropy[i][a] / (64.0 * 64.0)) * log10(histoEntropy[i][a] / (64.0 * 64.0)) / log(2.0);
				//cout << histoEntropy[i][a] / (64.0 * 64.0) << endl;
			}
		}
		entropy = -1 * entropy;
		cout << "Total Prob : " << check;
		check = 0.0;
		cout << " Entropy : " << entropy << endl;
		sumEntropy += entropy;
		entropy = 0.0;
	}
	cout << "Total Entropy : " << sumEntropy << endl;
}

int getValue(int x,int y,int w){
	return w * y + x;
}


void DCT(int(&matInt)[8][8]){
	//DCT T Element
	for (int a = 0; a < 8; a++) {
		for (int b = 0; b < 8; b++) {
			matInt[a][b] -= 128;	
		}
	}
	float cv, cu;
	int i, j, u, v;
	for (u = 0; u < 8; ++u) {
		for (v = 0; v < 8; ++v) {
			if (u == 0)cu = 1 / sqrt(2); else cu = 1.0;
			if (v == 0)cv = 1 / sqrt(2); else cv = 1.0;
			DCTMatrix[u][v] = 0;
			for (i = 0; i < 8; i++) {
				for (j = 0; j < 8; j++) {
					DCTMatrix[u][v] += (int)(matInt[i][j] * cos(((2 * i + 1)*u*PI) / 16) * cos(((2 * j + 1)*v*PI) / 16));
				}
			}
			DCTMatrix[u][v] = (int)(DCTMatrix[u][v] * 0.25 * cu *cv);
			if (DCTMatrix[u][v] > 2047)DCTMatrix[u][v] = 2047;
			if (DCTMatrix[u][v] < -2048)DCTMatrix[u][v] = -2048;
		}
		
	}
}

int sgn(float val) {
	if(val < 0) return -1;
	if(val == 0) return 0;
	if(val > 0) return 1;
}

void Quantize(int(&zigzag)[8][8], float step, int retain){
	int u, v,n;
	n = 0;
	for (u = 0; u < 8; ++u) {
		for (v = 0; v < 8; ++v) {
			if ((n) < (retain)){
				if (quantizebool){
					quantI[u][v] = 0;
					quantI[u][v] = round((zigzag[u][v] + sgn(zigzag[u][v] * (step / 2))) / step);
				}else{
					quantI[u][v] = zigzag[u][v];
				}
			}
			else{
				quantI[u][v] = 0;
			}
			n++;
			//cout << quantI[u][v] << " ";
		}
		//cout << "\n";
	}
	//cout << "\n";
}
void IQuantize(int(&quantI)[8][8], int step, int retain){
	int u, v,n;
	n = 0;
	for (u = 0; u < 8; ++u) {
		for (v = 0; v < 8; ++v) {
			if ((n) < (retain)){
				IquantI[u][v] = 0;
				if (quantizebool){
					IquantI[u][v] = quantI[u][v] * step;
				}
				else{
					IquantI[u][v] = quantI[u][v];
				}
				
			}else{
				IquantI[u][v] = 0;
			}
			n++;
			//cout << IquantI[u][v] << " ";
		}
		//cout << "\n";
	}
}

void IDCT(int(&izigzag)[8][8]){
	//DCT T Element
	float cv, cu;
	int u, v, m, n;
	for (m = 0; m < 8; ++m) {
		for (n = 0; n < 8; ++n) {

			IDCTMatrix[m][n] = 0;
			for (u = 0; u < 8; u++) {
				for (v = 0; v < 8; v++) {
					if (u == 0)cu = 1 / sqrt(2); else cu = 1.0;
					if (v == 0)cv = 1 / sqrt(2); else cv = 1.0;
					IDCTMatrix[m][n] += (int)(cu *cv* izigzag[u][v] * cos(((2 * m + 1)*u*PI) / 16) * cos(((2 * n + 1)*v*PI) / 16));
				}
			}
			IDCTMatrix[m][n] = (int)(round(IDCTMatrix[m][n] * 0.25) + 128);
			if (IDCTMatrix[m][n] > 255)IDCTMatrix[m][n] = 255;
			if (IDCTMatrix[m][n] < 0)IDCTMatrix[m][n] = 0;
		}

	}

}

void zigzagOrder(int(&DCTMatrix)[8][8])
{
	int i = 0;
	int j = 0;
	int n = 0;
	//for upper triangle of matrix
	zigzag[n / 8][n % 8] = DCTMatrix[i][j]; n++;
	do
	{
		j++;
		zigzag[n / 8][n % 8] = DCTMatrix[i][j]; n++;

		while (j != 0)
		{
			i++;
			j--;

			zigzag[n / 8][n % 8] = DCTMatrix[i][j]; n++;
		}
		i++;
		if (i>7)
		{
			i--;
			break;
		}

		zigzag[n / 8][n % 8] = DCTMatrix[i][j]; n++;

		while (i != 0)
		{
			i--;
			j++;
			zigzag[n / 8][n % 8] = DCTMatrix[i][j]; n++;
		}
	} while (true);

	//for lower triangle of matrix
	do
	{
		j++;
		zigzag[n / 8][n % 8] = DCTMatrix[i][j]; n++;

		while (j != 7)
		{
			j++;
			i--;

			zigzag[n / 8][n % 8] = DCTMatrix[i][j]; n++;
		}
		i++;
		if (i>7)
		{
			i--;
			break;
		}

		zigzag[n / 8][n % 8] = DCTMatrix[i][j]; n++;

		while (i != 7)
		{
			i++;
			j--;
			zigzag[n / 8][n % 8] = DCTMatrix[i][j]; n++;
		}
	} while (true);

	//for (int k = 0; k < 8; k++){
	//	for (int l = 0; l < 8; l++){
	//		//Extract int Value from Image
	//		cout << zigzag[k][l]<<" ";

	//		//cout << matInt[l][k] << " ";
	//	}
	//	cout << "\n";
	//}//end 8x8 block
}

void izigzagOrder(int(&IquantI)[8][8])
{
	int i = 0;
	int j = 0;
	int n = 0;
	int temp[64];
	for (int k = 0; k < 8; k++){
		for (int l = 0; l < 8; l++){
			//Extract int Value from Image
			temp[n] = IquantI[k][l]; n++;
			//cout << izigzag[k][l] << " ";
			//cout << matInt[l][k] << " ";
		}
		//cout << "\n";
	}
	n = 0;
	//for upper triangle of matrix
	izigzag[i][j] = temp[n]; n++;
	do
	{
		j++;
		izigzag[i][j] = temp[n]; n++;

		while (j != 0)
		{
			i++;
			j--;

			izigzag[i][j] = temp[n]; n++;
		}
		i++;
		if (i>7)
		{
			i--;
			break;
		}

		izigzag[i][j] = temp[n]; n++;

		while (i != 0)
		{
			i--;
			j++;
			izigzag[i][j] = temp[n]; n++;
		}
	} while (true);

	//for lower triangle of matrix
	do
	{
		j++;
		izigzag[i][j] = temp[n]; n++;

		while (j != 7)
		{
			j++;
			i--;

			izigzag[i][j] = temp[n]; n++;
		}
		i++;
		if (i>7)
		{
			i--;
			break;
		}

		izigzag[i][j] = temp[n]; n++;

		while (i != 7)
		{
			i++;
			j--;
			izigzag[i][j] = temp[n]; n++;
		}
	} while (true);

	n = 0;
	//for (int k = 0; k < 8; k++){
	//	for (int l = 0; l < 8; l++){
	//		//Extract int Value from Image
	//		cout << izigzag[k][l] << " ";
	//		//cout << matInt[l][k] << " ";
	//	}
	//	cout << "\n";
	//}//end 8x8 block
}

double mse(unsigned char asli[512 * 512], unsigned char hasil[512 * 512], int size){
	int selisih = 0.0; // different
	int total = 0.0; // total=mse 

	for (int i = 0; i < size*size; i++){
		selisih = (int)asli[i] - (int)hasil[i];
		total += selisih*selisih;
	}

	return (double)(total / (size*size));
}

int main(){
	ifstream inFile;
	re = 0;
	cout << "Enter Retain Quantize : ";  cin >> re;
	char quan = 'y';
	cout << "Do you want to quantize ? (y/n) : ";  cin >> quan;
	if (quan == 'n'){
		quantizebool = false;
		quality = 1;
	}else{
		cout << "Enter Step Size : ";  cin >> quality;
	}
	entropylimit = re;
	inFile.open("LENA.raw", ios::binary);         // binary mode
	if (!inFile.good()){

		return false;
	}

	
	inFile.read((char*)ori, dataSize);
	bmp.save("LENA.bmp", 512, 512, 1, ori); // save new image 
	/*int aMatrix[8][8] = {
		{52,55,61,66,70,61,64,73},
		{63,59,55,90,109,85,69,72},
		{62,59,68,113,144,104,66,73},
		{63,58,71,122,154,106,70,69},
		{67,61,68,104,126,88,68,70},
		{79,65,60,70,77,68,58,75},
		{85,71,64,59,55,61,65,83},
		{87,79,69,68,65,76,78,94}
	};*/

	for (int x = 0; x < 512 * 512; x++){
		comp[x] = 0;
	}
	for (int i = 0; i < 512;){
		for (int j = 0; j < 512;){
			//Build 8x8 block
			for (int k = 0; k < 8; k++){
				for (int l = 0; l < 8; l++){
					//Extract int Value from Image
					matInt[l][k] = (int)ori[getValue(l + j, k + i,512)];
					//cout << matInt[l][k] << " ";
				}
				//cout << "\n";
			}//end 8x8 block
			//return 0;
			DCT(matInt);
			zigzagOrder(DCTMatrix);
			Quantize(zigzag, quality, re);

			int localcounter = 0;
			for (int k = 0; k < 8; k++){
				for (int l = 0; l < 8; l++){
					if (localcounter<entropylimit){
						probData5[counterData][localcounter] = quantI[k][l];
						localcounter++;
					}
					//cout << matInt[l][k] << " ";
				}
				//cout << "\n";
			}//end 8x8 block
			counterData++;

			IQuantize(quantI, quality, re);
			izigzagOrder(IquantI);
			IDCT(izigzag);
		
			for (int k = 0; k < 8; k++){
				for (int l = 0; l < 8; l++){
					//Extract int Value from Image
					comp[getValue(l + j, k + i, 512)] = IDCTMatrix[l][k];
					//cout << matInt[l][k] << " ";
				}
				//cout << "\n";
			}//end 8x8 block
			
			
			
			
			j = (j + 8);
		}
		i = (i + 8);
	}
	calcEntropy(entropylimit);
	/*DCT(aMatrix);
	Quantize(DCTMatrix, 16);
	zigzagOrder(quantI);
	cout << "\n";
	izigzagOrder(zigzag);*/
	//IQuantize(izigzag, 16);
	//IDCT(IquantI);
	
	double m = mse(ori, comp, 512);
	cout << "MSE : "<< m<<endl;
	bmp.save("LENAcom.bmp", 512, 512, 1, comp); // save new image 
	return 0;
}
