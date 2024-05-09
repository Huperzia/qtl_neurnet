import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.io.*;
import java.util.Random;

public class NeuralNet {

	public static void main(String[] args) {
		int y = 0;
		String xy;
		String delimiter = "[\t]";
		String[][] lines = new String[100][1000];
		String[][] snp = new String[100][1000];
		String[][] allele = new String[100][1000];
		String[][] cov = new String[100][1000];
		List<String> lbl = new ArrayList<String>();
		
	try
	{BufferedReader inn = new BufferedReader(new FileReader("./ExFSNPs.txt"));

		while ((xy = inn.readLine()) != null)
		{	if (xy.contains(">"))
				{	y = y + 1;
					lbl.add(xy);	}
			
			else
				{	String[] to = xy.split(delimiter);
					int too = Integer.parseInt(to[0]);
					lines[too][y] = to[4];
					cov[too][y] = to[1];
					snp[too][y] = to[2];
					allele[too][y] = to[3];
				}		
		}
		
			String[] Snp = lbl.toArray(new String[lbl.size()]);

		inn.close();
	}	
	catch (Exception eee)
	{System.out.println(eee);}
	
	Scanner scan = new Scanner(System.in);
//	String er = scan.nextLine();
	
	String xy1;
	List<String> lbl1 = new ArrayList<String>();
	double[] MNDI = new double[100];
	double[] MNDX = new double[100];
	double[] SCN = new double[100];
	double[] YLDSDS = new double[100];
	double[] MNDS = new double[100];
	int L = 0;
	
	try
	{	BufferedReader in = new BufferedReader(new FileReader("./ExFphen.txt"));
	
		while ((xy1 = in.readLine()) != null)
		{ 	if (L == 0)
				{ 	L = L + 1;	
					lbl1.add(xy1);	
				}
			else
			{	try
				{	String[] to = xy1.split(delimiter);
					int too = Integer.valueOf(to[0]);
					
					if(to[1] != null)
					{	MNDI[too] = Double.valueOf(to[1]);	}
					else
					{	MNDI[too] = 0;	}
							
					if(to[2] != null)
					{	MNDS[too] = Double.valueOf(to[2]);	}
					else
					{	MNDS[too] = 0;	}
							
					if(to[3] != null)
					{	MNDX[too] = Double.valueOf(to[3]);	}
					else
					{	MNDX[too] = 0;	}
							
					if(to[4] != null)
					{	SCN[too] = Double.valueOf(to[4]);	}
					else
					{	SCN[too] = 0;	}
					
					try
					{	if(to[5] != null)
						{	YLDSDS[too] = Double.valueOf(to[5]);	}
						else
						{	YLDSDS[too] = 0;	}
					}
					catch(Exception ex)
					{	System.out.println(ex);	}
					
				}
				catch(Exception ex)
				{ }
			}		
		}
		
		//String[] Lbl = lbl1.toArray(new String[lbl1.size()]);
		in.close();
	}	
	catch (Exception e)
	{	System.out.println(e);	}

	
		int n = 308;
		int l = 3;
		int q = 2;
		int P = 100;
		int M = 6;
		
		double[][] z = new double [1000][115];
		double[][] t = new double [1000][115];
		double[][] Y = new double [1000][115];
		double[][] D = new double [1000][115];

		double[][] Wij = new double [500][115];
		double[][] Wjk = new double [500][115];
		double[][] Wkm = new double [500][115];
		
		double[] tWij = new double [500];
		double[] tWjk = new double [500];
		double[] tWkm = new double [500];
		
		double[][] dWij = new double [500][115];
		double[][] dWjk = new double [500][115];
		double[][] dWkm = new double [500][115];
	
		double alpha[][] = new double[500][115];
		double beta[][] = new double[500][115];
		double gamma[][] = new double[500][115];
		
		double x[][] = new double[500][115];
		double e[][] = new double[500][115];
		double ep[] = new double[115];
		double N = 0.000001;
		
		
		System.out.println("Press enter to begin neural network");
		String err = scan.nextLine();
		Random weight = new Random();
		
		for (int u = 1; u < 500; u++)
		{	for (int v = 1; v < 115; v++)	
			{	Wij[u][v] = ((weight.nextInt(100) + 1));
				Wij[u][v] = Wij[u][v] * 0.0001;
				//Wij[u][v] = 0.001;
				Wjk[u][v] = 0.005;
				Wkm[u][v] = 1.0;
				x[u][v] = 0;
			}
		}
		
		for (int u = 1; u <= 300; u++)
		{	for (int v = 3; v < 100; v++)	
			{	if(lines[v][u] != null )
				{	String linee = lines[v][u].replaceAll("[^0-9]+", "");
					
					if(linee != null)
					{	int vvv = Integer.parseInt(linee);
						linee = null;
									
						if(allele[v][u].equals("PARENT_A"))
						{	x[u][vvv] = 10;	}

						else if(allele[v][u].equals("PARENT_B"))
						{	x[u][vvv] = -10; }	
						
						else
						{ 	x[u][vvv] = 0;	}
					}
				}
			}
		}
				
		for (int v = 1; v < 100; v++)	
			{	D[1][v] = SCN[v];
				D[2][v] = MNDS[v];
				D[3][v] = MNDX[v];
				D[4][v] = MNDI[v];
				D[5][v] = YLDSDS[v];	}
		
		for (int v = 100; v < 115; v++)	
		{	D[1][v] = 0;
			D[2][v] = 0;
			D[3][v] = 0;
			D[4][v] = 0;
			D[5][v] = 0;	}
		
		for (int u0 = 0; u0 < 1; u0++)
		{	for (int v = 0; v < 115; v++)	
			{	Wij[u0][v] = 1.0;
				Wjk[u0][v] = 1.0;
				Wkm[u0][v] = 1.0;	
			}
		}
		
		for (int b = 1; b < 1000; b++)
		{	for(int p = 1; p <= P; p++)
			{	ep[p] = 0;
				
				for	(int i = 1; i < n; i++)
				{	tWij[i] = 0;
					
					for (int j = 1; j < l; j++)
					{	tWjk[j] = 0;
						alpha[i][j] = 0;
						alpha[i][j] += Wij[0][j];	
						
						for (int u = 0; u < n; u++)
						{	if(x[u][p] != 0)
							{	alpha[i][j] += (Wij[u][j] * x[u][p]);	}
						}
						
						if ((Math.exp(alpha[i][j]) + Math.exp(-alpha[i][j]) != 0))
						{z[j][p] = (Math.exp(alpha[i][j]) - Math.exp(-alpha[i][j])) / ((Math.exp(alpha[i][j]) + Math.exp(-alpha[i][j])));}
						
						for (int k = 1; k < q; k++)
						{	tWkm[k] = 0;
							beta[j][k] = 0;
							beta[j][k] += Wjk[0][k];	
							
							for (int u = 0; u < l; u++)
							{	beta[j][k] += (Wjk[u][k] * z[u][p]);	}
							
							if (Math.exp(beta[j][k]) + Math.exp(-beta[j][k]) != 0)
							{t[k][p] = (Math.exp(beta[j][k]) - Math.exp(-beta[j][k])) / ((Math.exp(beta[j][k]) + Math.exp(-beta[j][k])));}	
							
							for (int m = 1; m < 2; m++)
							{	gamma[k][m] = 0;
								gamma[k][m] += Wkm[0][m];	
										
								for (int u = 0; u < q; u++)
								{	gamma[k][m] += (Wkm[u][m] * t[u][p]);	}
								
								if (Math.exp(gamma[k][m]) + Math.exp(-gamma[k][m]) != 0)
								{Y[m][p] = (Math.exp(gamma[k][m]) - Math.exp(-gamma[k][m])) / ((Math.exp(gamma[k][m]) + Math.exp(-gamma[k][m])));}
								
								if(Y[m][p] == Double.NaN)
								{int sdhf = scan.nextInt();}
								
								if(D[m][p] != 0 && Y[m][p] !=0)
								{	e[m][p] = Math.pow((Y[m][p] - D[m][p]), 2);
									ep[p] += e[m][p];
									
									Wkm[0][m] = Wkm[0][m] + (-N) * Y[m][p] * (1 - Y[m][p]) * (Y[m][p] - D[m][p]);
									dWkm[k][m] =  (-N) * t[k][p] * Y[m][p] * (1 - Y[m][p]) * (Y[m][p] - D[m][p]);
									Wkm[k][m] = Wkm[k][m] + dWkm[k][m];
									tWkm[k] += Wkm[k][m];
								}
							}
							
							if (tWkm[k] != 0)
							{	Wjk[0][k] =  Wjk[0][k] + (-N) * t[k][p] * (1 - t[k][p]) * tWkm[k];
								dWjk[j][k] =  (-N) * z[j][p] * t[k][p] * (1 - t[k][p]) * tWkm[k];
								Wjk[j][k] = Wjk[j][k] + dWjk[j][k];
								tWjk[j] += Wjk[j][k];
							}
						}
						
						if (tWjk[i] != 0)
						{	Wij[0][i] = Wij[0][i] + ((-N) * z[j][p] * (1 - z[j][p]) * tWjk[j]);
							
							if(x[i][j] != 0)
							{	dWij[i][j] =  (-N) * x[i][j] * z[j][p] * (1 - z[j][p]) * tWjk[j];
								Wij[i][j] = Wij[i][j] + dWij[i][j];
								tWij[i] += Wij[i][j];	
								//System.out.println("Wij: " + i + " " + j + " " + Wij[i][j]);
								if(dWij[i][j] == Double.NaN)
								{int sdhf = scan.nextInt();}
							}
						}
					}
				}
				
				if (p == 1)
				{System.out.println("Error over p: " + b + " " + ep[p]);}
			}
			//System.out.println(b);
		}
		
		for	(int ii = 1; ii < n; ii++)
		{	for (int jj = 1; jj < l; jj++)
			{	try (Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("Output.txt",true), "UTF-8")))
				{	writer.write("weight\t" + ii + "\t" + jj + "\t" + Wij[ii][jj] + "\n");
					writer.close();	}
				catch (Exception exc)
				{	System.out.println(exc);	}
			}
		}
		
		for (int u0 = 0; u0 < 1; u0++)
		{	for (int v = 0; v < 115; v++)	
			{	System.out.println("Wij " + u0 + " " + v + ": " + Wij[u0][v]);
				
			}
		System.out.println("Wjk " + Wjk[0][1]);
		System.out.println("Wjk: " + Wkm[0][1]);
		}
		
		scan.close();
	}
}
