#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int main() {
	int pause;
	char *str;
	char *fileName;
	char *preFile = "ExFSNPs.txt";
	char *path = "./";


	str = (char*)malloc(80 * sizeof(char));
	fileName = (char*)malloc(80 * sizeof(char));

// Load SNPs
	strcpy(fileName, path);
	strcat(fileName, preFile);

	FILE *myData;

	myData = fopen(fileName, "r");
	if (myData == NULL)
	{
		printf("No file found\n");
		exit(1);
	}

	char *snp[310];
	char *snp1[310];
	char *snp2[310];
	char ch;
	char *ril[102][311];
	char tempsnp;

	int cnt = 0;
	int snpTot = -1;
	int rilnum = 0;
	int linenum = 0;
	char *riln;
	char *snpName;

	for (int i = 0; i < 102; i++)
	{
		for (int j = 0; j < 311; j++) {
			ril[i][j] = (char*)malloc(sizeof(char));
			*ril[i][j] = '0';
		}
	}

	while ((ch = getc(myData)) != EOF)
	{
		while ((ch != '\n'))
		{
			if (ch == '>')
			{
				snpTot++;
				snp[snpTot] = (char*)malloc(20 * sizeof(char));
				snp1[snpTot] = (char*)malloc(sizeof(char));
				snp2[snpTot] = (char*)malloc(sizeof(char));

				while (ch != '\t')
				{
					*(snp[snpTot] + cnt) = ch;
					ch = getc(myData);
					cnt++;
				}
				ch = getc(myData);
				cnt = 0;

				while (ch != ':')
				{
					ch = getc(myData);
					cnt++;
				}
				ch = getc(myData);

				*(snp1[snpTot]) = ch;
				cnt = 0;
				ch = getc(myData);

				while (ch != ':')
				{
					ch = getc(myData);
					cnt++;
				}
				ch = getc(myData);
				*(snp2[snpTot]) = ch;
				cnt = 0;
				linenum = 0;

				while (ch != '\n')
				{
					ch = getc(myData);
				}
			}

			else {
				cnt = 0;
				while (cnt < 2)
				{
					ch = getc(myData);

					if (ch == '\t')
					{
						cnt++;
					}
				}
				ch = getc(myData);
				tempsnp = ch;

				cnt = 0;
				while (cnt < 2)
				{
					ch = getc(myData);
					if (ch == '\t')
					{
						cnt++;
					}
				}

				ch = getc(myData);

				if (linenum > 2)
				{
					while (ch != '-')
					{
						ch = getc(myData);
					}

					cnt = 0;
					riln = (char*)malloc(2 * sizeof(char));
					ch = getc(myData);

					while (ch != '\t')
					{
						*(riln + cnt) = ch;
						ch = getc(myData);
						cnt++;
					}
					if (cnt == 2) {
						rilnum = *(riln)-'0';
						rilnum *= 10;
						rilnum += *(riln + 1) - '0';
						free(riln);
					}
					else if (cnt == 1) {
						rilnum = *(riln)-'0';
						free(riln);
					}
					cnt = 0;
				}

				else if (linenum < 3)
				{
					if (ch == 'E') {
						rilnum = 98;
					}
					else if (ch == 'F') {
						rilnum = 99;
					}
				}

				if (rilnum > 0 && rilnum < 100)
				{
					ril[rilnum][snpTot] = (char*)malloc(sizeof(char));
					*ril[rilnum][snpTot] = tempsnp;
				}

				while (ch != '\n')
				{
					ch = getc(myData);
				}
			}

			linenum++;
			rilnum = 0;
		}
	}

	printf("\nReached the end of the SNPs with %d SNPs\n\nYour SNPs are now loaded.\n\n\n", snpTot);

// Load Phenotypes
	FILE *myPhen;
	preFile = "ExFphen.txt";

	strcpy(fileName, path);
	strcat(fileName, preFile);

	myPhen = fopen(fileName, "r");
	if (myPhen == NULL)
	{
		printf("File not found\n");
		exit(1);
	}

	char buf[100];
	char *buff;
	float MNDI[100];
	float MNDS[100];
	float MNDX[100];
	float FI[100];
	float YLDSDS[100];

	buff = (char*)malloc(10 * sizeof(char));
	linenum = 0;

	while (fgets(buf, 100, myPhen) != NULL)
	{
		if (linenum > 0)
		{
			riln = strtok(buf, "\t");
			sscanf(riln, "%d", &rilnum);

			buff = strtok(NULL, "\t");
			sscanf(buff, "%f", &MNDI[rilnum]);

			buff = strtok(NULL, "\t");
			sscanf(buff, "%f", &MNDS[rilnum]);

			buff = strtok(NULL, "\t");
			sscanf(buff, "%f", &MNDX[rilnum]);

			buff = strtok(NULL, "\t");
			sscanf(buff, "%f", &FI[rilnum]);

			buff = strtok(NULL, "\t");
			sscanf(buff, "%f", &YLDSDS[rilnum]);
		}
		linenum++;
	}





	//START NEURNET

	printf("\n\nWelcome to NeurNet!\n\n");
	printf("\nYour SNPs and Phenotypes have been loaded.\n\n");

	cnt = 0;

	int datTot = 10000;
	int n = 308;
	int r = 101;
	int l = 2;
	int q = 2;
	int P = 97;
	int M = 6;

	double z[n][r];
	double t[n][r];
	double Y[n][r];
	double D[n][r];

	double Wij[n][r];
	double Wjk[n][r];
	double Wkm[n][r];

	double tWij[n];
	double tWjk[n];
	double tWkm[n];

	double dWij[n][r];
	double dWjk[n][r];
	double dWkm[n][r];

	double alpha[n][r];
	double beta[n][r];
	double gamma[n][r];

	double x[n][r];
	double e[n][r];
	double ep[r];
	double N = 0.000001;

	printf("\n\nVariables have their allocated memory\n");


	for (int u = 1; u <= 308; u++)
	{
		for (int v = 1; v < 102; v++)
		{
			Wij[u][v] = 0.001;
			Wjk[u][v] = 1.0;
			Wkm[u][v] = 1.0;
			x[u][v] = 0;
		}
	}

	for (int u0 = 0; u0 < 1; u0++)
	{
		for (int v = 0; v < 102; v++)
		{
			Wij[u0][v] = 1.0;
			Wjk[u0][v] = 1.0;
			Wkm[u0][v] = 1.0;
		}
	}

	printf("Weights have been assigned.\n");


	for (int u = 1; u <= 308; u++)
	{
		for (int v = 1; v < 101; v++)
		{
			if (*ril[v][u] != '0')
			{
				//ESSEX INITIAL
				if (*ril[v][u] == *ril[98][u])
				{
					x[u][v] += 0.9;
				}

				//FORREST INITIAL
				else if (*ril[v][u] == *ril[99][u])
				{
					x[u][v] += 0.1;
				}
			}
			else
			{
				x[u][v] = 0.5;
			}
		}
	}

	printf("RILs have been assigned\n");

	for (int v = 1; v < 100; v++)
	{
		D[1][v] = (FI[v] / 100);
//		D[2][v] = MNDS[v];
//		D[3][v] = MNDX[v];
//		D[4][v] = MNDI[v];
//		D[5][v] = YLDSDS[lll];
	}

	printf("Phenotypes have been assigned\n");
	//	free(FI);




	printf("\n\nBegin Neurnet\n\n");

	for (int b = 1; b < datTot; b++)
	{
		for (int p = 1; p <= P; p++)
		{
			ep[p] = 0;

			for (int i = 1; i < n; i++)
			{
				tWij[i] = 0;

				for (int j = 1; j < l; j++)
				{
					tWjk[j] = 0;
					alpha[i][j] = 0;
					alpha[i][j] += Wij[0][j];

					for (int u = 1; u < n; u++)
					{
						if (x[u][p] != 0.0)
						{
							alpha[i][j] += (Wij[u][j] * x[u][p]);
						}
					}

					if (((exp(alpha[i][j]) + exp(-alpha[i][j])) != 0))
					{
						z[j][p] = ((1 ) / (1 + exp(-alpha[i][j])));
					}

					for (int k = 1; k < q; k++)
					{
						tWkm[k] = 0;
						beta[j][k] = 0;
						beta[j][k] += Wjk[0][k];

						for (int u = 1; u < l; u++)
						{
							beta[j][k] += (Wjk[u][k] * z[u][p]);
						}

						if ((exp(beta[j][k]) + exp(-beta[j][k])) != 0)
						{
							t[k][p] = ((1 ) / (1 + exp(-beta[j][k])));
						}

						for (int m = 1; m < 2; m++)
						{
							gamma[k][m] = 0;
							gamma[k][m] += Wkm[0][m];

							for (int u = 1; u < q; u++)
							{
								gamma[k][m] += (Wkm[u][m] * t[u][p]);
							}

							if ((exp(gamma[k][m]) + exp(-gamma[k][m])) != 0)
							{
								Y[m][p] = ((1 ) / ((1 + exp(-gamma[k][m]))));
							}

							if (D[m][p] != 0 && Y[m][p] != 0)
							{
								e[m][p] = pow(((Y[m][p] - D[m][p])), 2);
								ep[p] += e[m][p];

								Wkm[0][m] = Wkm[0][m] + (-N) * Y[m][p] * (1 - Y[m][p]) * (Y[m][p] - D[m][p]);
								dWkm[k][m] = (-N) * t[k][p] * Y[m][p] * (1 - Y[m][p]) * (Y[m][p] - D[m][p]);
								Wkm[k][m] = Wkm[k][m] + dWkm[k][m];
								tWkm[k] += Wkm[k][m];
							}
						}

						if (tWkm[k] != 0)
						{
							Wjk[0][k] = Wjk[0][k] + (-N) * t[k][p] * (1 - t[k][p]) * tWkm[k];
							dWjk[j][k] = (-N) * z[j][p] * t[k][p] * (1 - t[k][p]) * tWkm[k];
							Wjk[j][k] = Wjk[j][k] + dWjk[j][k];
							tWjk[j] += Wjk[j][k];
						}
					}

					if (tWjk[j] != 0)
					{
						Wij[0][j] = Wij[0][j] + (-N) * z[j][p] * (1 - z[j][p]) * tWjk[j];

						if (x[i][j] != 0)
						{
							dWij[i][j] = (-N) * x[i][j] * z[j][p] * (1 - z[j][p]) * tWjk[j];
							Wij[i][j] = Wij[i][j] + dWij[i][j];
							tWij[i] += Wij[i][j];
						}
					}
				}
			}

			if (p == 1)
			{
				cnt++;
				if (cnt % 20 == 0)
				{
					for (int i = 0; i < 62; i++)
					{
						for (int j = 1; j < l; j++)
						{
							for (int iii = 0; iii < 5; iii++) {
								if (iii < 3) printf("[%d] %f\t", (i + (iii * 62)), (Wij[i + (iii * 62)][j] * 1000));
								if (iii >= 3 && i < 61) { printf("[%d] %f\t", (i + (3 * 62) + ((iii - 3) * 61)), (Wij[(i + (3 * 62) + ((iii - 3) * 61))][j] * 1000)); }
							}
						}
					}
					printf("\n\t\t\t\t\t\t   Error over p is %lf\n\n", ep[p]);
				}


				if (ep[p] < 0.004)
				{
					preFile = "weights.txt";
					strcpy(fileName, path);
					strcat(fileName, preFile);

					FILE *weights = fopen(fileName, "w");
					if (weights == NULL)
					{
						printf("Error");
						exit(1);
					}


					for (int i = 0; i < 300; i++)
					{
						for (int j = 1; j < 2; j++)
						{
							fprintf(weights, "Wij[%d][%d] = %f\n", i, j, Wij[i][j]);
						}
					}
				}
			}
		}
	}


	printf("END OF PROGRAM");
}