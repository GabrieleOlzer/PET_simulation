/*
Programma per simulazione Montecarlo di esperimento di PET monodimensionale; fissata la posizione della sorgente radioattiva simula l'emissione di coppie di fotoni back to back (lungo direzioni casuali) e determina quanti di questi vengono assorbiti in entrambi i rivelatori (producendo un conteggio). La simulazione è ripetuta per diversi angoli di rotazione del rivelatore e per diverse posizioni della sorgente; per ogni posizione della sorgente viene calcolato l'angolo in corrispondenza del quale si ha un massimo dei conteggi, il programma restituisce un grafico in cui è plottato l'angolo dei conteggi massimi in funzione della posizione della sorgente (nel grafico è anche plottato il modello analitico  dell'andamento del massimo in funzione della posizione se i rivelatori fossero considerabili come puntiformi)
*/


#include<iostream>
#include<vector>
#include<cstdlib>
#include<fstream>
#include<string>
#include<time.h>
#include<cmath>

#include"TGraph.h"
#include"TCanvas.h"
#include"TString.h"
#include"TApplication.h"
#include"TH1F.h"
#include"TF1.h"

#define _USE_MATH_DEFINES
#define r 2.54 //raggio rivelatore cm
#define l 5.08 //lunghezza rivelatore cm
#define k 0.45 //inverso della lunghezza di radiazione (cm)-1
#define grad 0.0174533 //radianti per grado
#define rad 57.29578 //gradi per radiante

using namespace std;

//modello analitico, restiruisce l'angolo associato al massimo dei conteggi in funzione della poszione della sorgente lungo l'asse y (rivelatori approssimati puntiformi)
double alpha_max (double* x, double* par)
{
	double max = rad*(atan(x[0]/par[0]) + asin((x[0]/par[0])*cos(atan(x[0]/par[0]))));
	
	return max;
}

//genera numeri casuali secondo una distribuzione uniforme definita tra min e max
double rand_range (double min, double max)
{
	return min + (max - min)*(rand()/static_cast<double>(RAND_MAX));
}

//calcola la distaznza dalla sorgente (in cm) a cui si ha l'ingresso di un fotone nello spettrometro (a meno di dx/2); s viene usato per parametrizzare la retta lungo la quale si ha la propagazione del fotone
double bisezione_in (double s_min, double s_max, double x0, double y0, double theta, double phi, double dx, double d, double alpha)
{
	double s_avg = 0.5*(s_min + s_max);
	
	if (s_max - s_min < dx) return s_avg;
	
	double x_s = x0 + sin(theta)*cos(phi)*s_avg;
	double y_s = y0 + sin(theta)*sin(phi)*s_avg;
	double z_s = cos(theta)*s_avg;
	
	//cordinate rispetto ad un sistema ruotato attorno a z di un angolo alpha
	double x_p = cos(alpha)*x_s + sin(alpha)*y_s;
	double y_p = cos(alpha)*y_s - sin(alpha)*x_s;
	
	if (pow(y_p, 2) + pow(z_s, 2) <= pow(r, 2) && x_p > d && x_p < l + d)
	{
		s_avg = bisezione_in(s_min, s_avg, x0, y0, theta, phi, dx, d, alpha);
	}
	else s_avg = bisezione_in(s_avg, s_max, x0, y0, theta, phi, dx, d, alpha);
	
	return s_avg;	
}

//distanza dalla sorgente a cui si ha l'uscita dallo spettrometro
double bisezione_out (double s_min, double s_max, double x0, double y0, double theta, double phi, double dx, double d, double alpha)
{
	double s_avg = 0.5*(s_min + s_max);
	
	if (s_max - s_min < dx) return s_avg;
	
	double x_s = x0 + sin(theta)*cos(phi)*s_avg;
	double y_s = y0 + sin(theta)*sin(phi)*s_avg;
	double z_s = cos(theta)*s_avg;
	
	double x_p = cos(alpha)*x_s + sin(alpha)*y_s;
	double y_p = cos(alpha)*y_s - sin(alpha)*x_s;
	
	if (pow(y_p, 2) + pow(z_s, 2) <= pow(r, 2) && x_p > d && x_p < l + d)
	{
		s_avg = bisezione_out(s_avg, s_max, x0, y0, theta, phi, dx, d, alpha);
	}
	else s_avg = bisezione_out(s_min, s_avg, x0, y0, theta, phi, dx, d, alpha);
	
	return s_avg;	
}

//distanza dalla sorgente a cui si ha l'ingresso nel rivelatore di gate
double bisezione_in_g (double s_min, double s_max, double x0, double y0, double theta, double phi, double dx, double d)
{
	double s_avg = 0.5*(s_min + s_max);
	
	if (s_max - s_min < dx) return s_avg;
	
	double x_s = x0 - sin(theta)*cos(phi)*s_avg;
	double y_s = y0 - sin(theta)*sin(phi)*s_avg;
	double z_s = 0. - cos(theta)*s_avg;
	
	if (pow(y_s, 2) + pow(z_s, 2) <= pow(r, 2) && x_s < -d && x_s > -(l + d))
	{
		s_avg = bisezione_in_g(s_min, s_avg, x0, y0, theta, phi, dx, d);
	}
	else s_avg = bisezione_in_g(s_avg, s_max, x0, y0, theta, phi, dx, d);
	
	return s_avg;	
}

//distanza a cui si ha l'uscita dal rivelatore di gate
double bisezione_out_g (double s_min, double s_max, double x0, double y0, double theta, double phi, double dx, double d)
{
	double s_avg = 0.5*(s_min + s_max);
	
	if (s_max - s_min < dx) return s_avg;
	
	double x_s = x0 - sin(theta)*cos(phi)*s_avg;
	double y_s = y0 - sin(theta)*sin(phi)*s_avg;
	double z_s = 0. - cos(theta)*s_avg;
	
	if (pow(y_s, 2) + pow(z_s, 2) <= pow(r, 2) && x_s < -d && x_s > -(l + d))
	{
		s_avg = bisezione_out_g(s_avg, s_max, x0, y0, theta, phi, dx, d);
	}
	else s_avg = bisezione_out_g(s_min, s_avg, x0, y0, theta, phi, dx, d);
	
	return s_avg;	
}
	


int main (int argc, char** argv)
{
	
	//distanza dei rivelatori dall'origine del sistema di riferimento
	cout << "Inserire la distanza d: ";
	double d;
	cin >> d;
	
	//quaddrato della massima distanza dal centro alla quale un fotone si può trovare all'interno di un rivelatore
	double R_squared = pow(d + l, 2) + pow(r, 2); 
	
	//vettori delle posizioni della sorgente lungo l'asse y e degli angoli di massimo conteggio associati
	vector<double> Y;
	vector<double> alpha_max_Y;
	
	//posizione iniziale della sorgente
	double x0 = 0;
	double y0 = 0.;
		
	//ciclo sulle varie posizione della sorgente (movimenti a passo di 0.5 cm)
	while (y0 <= 12.)
	
	{
		//vettori contenenti gli angoli ed i conteggi associati a ciascuno di essi
		vector<double> alphas;
		vector<int> N_alpha; 
		
		srand(time(NULL));
		
		//il range di angoli considerato va da -20 a 90 gradi
		double alpha =  -20.* grad;
	
		while (alpha <= 90.*grad)
		{
		
			int N = 0; //inizializzazione numero conteggi
			
			
			for (int i = 0; i < atoi(argv[1]); i++)
			{
			
				double theta; //angolo formato con asse z
				double phi; //angolo formato dalla proiezione su xy con l'asse x
				
				theta = rand_range(0., M_PI);
				phi = rand_range(-M_PI, M_PI);
				
				/*CALCOLO DELLA DISTANZA PERCORSA ALL'INTERNO DELLO SPETTROMETRO*/
				
				double s = 0.; //parametrizza la retta
				
				double x = x0 + sin(theta)*cos(phi)*s;
				double y = y0 + sin(theta)*sin(phi);
				double z = cos(theta)*s;
				
				double x_p; //coordinate ruotate x'
				double y_p;
				
				//valutazione preliminare della posizione di entrata nello spetttrometro
				while (true)
				{
					x += sin(theta)*cos(phi);
					y += sin(theta)*sin(phi);
					z += cos(theta);
					
					x_p = cos(alpha)*x + sin(alpha)*y;
					y_p = cos(alpha)*y - sin(alpha)*x;
					
					if ( pow(y_p, 2) + pow(z, 2) <= pow(r, 2) && x_p > d && x_p < l + d)
					{
						s += 1.;
						
						break;	
					}
					
					//valuta se il fotone è uscito dalla regione di spazio in cui posso trovare un rivelatore
					if (pow(x, 2) + pow(y, 2) > R_squared || abs(z) > r)
					{
						s = 100.;
						break;
					}
					
					s += 1.;
					
				}
				
				double s_in = s;
				
				//interrompe la propagazione del fotone e passa al fotone successivo
				if (s_in == 100.) continue;
				
				//valutazione preliminare posizione uscita dallo spettrometro
				while (true)
				{
					x += sin(theta)*cos(phi);
					y += sin(theta)*sin(phi);
					z += cos(theta);
					
					x_p = cos(alpha)*x + sin(alpha)*y;
					y_p = cos(alpha)*y - sin(alpha)*x;
					
					if ( pow(y_p, 2) + pow(z, 2) > pow(r, 2) || x_p < d || x_p > l + d)
					{
						s += 1.;
						
						break;	
					}
					
					s += 1.;
				}
				
				/*CALCOLO DELLA DISTNZA PERCORSA ALL'INTERNO DEL RIVELATORE DI GATE, effettuata allo stesso modo*/
				
				double s_out = s;
				
				s = 0.;
				
				x = x0 - sin(theta)*cos(phi)*s;
				y = y0 - sin(theta)*sin(phi)*s;
				z = -cos(theta)*s;
				
				while (true)
				{
					x -= sin(theta)*cos(phi);
					y -= sin(theta)*sin(phi);
					z -= cos(theta);
					
					if ( pow(y, 2) + pow(z, 2) <= pow(r, 2) && x < -d && x > -(l + d))
					{
						s += 1.;
						
						break;	
					}
					
					if (pow(x, 2) + pow(y, 2) > R_squared || abs(z) > r)
					{
						s = 100.;
						break;
					}
					
					s += 1.;
					
				}
				
				double s_in_gate = s;
				
				if (s_in_gate == 100.) continue;
				
				while (true)
				{
					x -= sin(theta)*cos(phi);
					y -= sin(theta)*sin(phi);
					z -= cos(theta);
					
					if ( pow(y, 2) + pow(z, 2) > pow(r, 2) || x > -d || x < -(l + d))
					{
						s += 1.;
						
						break;	
					}
					
					s += 1.;
				}
				
				double s_out_gate = s;
				
				//calcolo delle posizoni precise di entrata e uscita dallo spettrometro 
				s_in = bisezione_in(s_in - 1., s_in, x0, y0, theta, phi, 0.01, d, alpha);
				s_out = bisezione_out(s_out - 1., s_out, x0, y0, theta, phi, 0.01, d, alpha);
				
				//distanza percorsa all'interno dello spettrometro
				double dist_riv = s_out - s_in;
				
				//proababilità che il fotone sia assorbito
				double P_assorbimento_riv = 1 - exp(-k * dist_riv);
				
				//calcolo delle posizoni precise di entrata e uscita dal rivelatore di gate
				s_in_gate = bisezione_in_g(s_in_gate - 1., s_in_gate, x0, y0, theta, phi, 0.01, d);
				s_out_gate = bisezione_out_g(s_out_gate - 1., s_out_gate, x0, y0, theta, phi, 0.01, d);
				
				double dist_gate = s_out_gate - s_in_gate;
				
				double P_assorbimento_gate = 1 - exp(-k * dist_gate);
				
				//calcolo della probabilità di assorbimento di enttrambi i fotoni (conteggio)
				double P = P_assorbimento_gate * P_assorbimento_riv;
				
				//si ha probabilità P di aggiunta di un conteggio
				if (rand_range(0., 1.) <= P) N++;
			
			}
			
			alphas.push_back(alpha*rad);
			N_alpha.push_back(N);
			
			alpha = alpha + 1.*grad;
		
		}
	
		/*
		Per ogni posizione si effettua un fit del grafico conteggi-angolo con un modello gaussiano, si usa la media della 			gaussiana di fit come stima per l'angolo che massimizza i conteggi
		*/
		TGraph graph;
		
		double alpha_max = alphas[0];
		double N_alpha_max = N_alpha[0];
		
		for (int i = 0; i < alphas.size(); i++)
		{
			graph.SetPoint(graph.GetN(), alphas[i], N_alpha[i]);
			
			if (N_alpha[i] > N_alpha_max)
			{
				alpha_max  = alphas[i];
				N_alpha_max = N_alpha[i];
			}
		}
		
		TF1 * fit_model = new TF1 ("fit_model", "gaus(0)", -20., 90.) ;
		fit_model->SetParameter(1, alpha_max);
		
		graph.Fit ("fit_model", "Q", "", -20., 90.);
		
		Y.push_back(y0);
		alpha_max_Y.push_back(fit_model->GetParameter(1));
		
	  	y0 += 0.5;
  	
  	}
  	
  	//Creazione e stampa del grafico alpha_N_max in fuznione della posizione y della sorgente
  	TGraph grafico;
  	grafico.SetMarkerStyle(20);
	grafico.SetMarkerColor(kBlue);
	grafico.SetMarkerSize(0.5);
	grafico.SetTitle ("AlphaMax - y0");
	grafico.GetHistogram ()->GetYaxis ()->SetTitle ("alpha max [°]");
	grafico.GetHistogram ()->GetXaxis ()->SetTitle ("y0 [cm]");
  	
  	for (int i = 0; i < Y.size(); i++)
  	{
  		grafico.SetPoint(grafico.GetN(), Y[i], alpha_max_Y[i]);
  	}
  	
  	TF1 modello_analitico("modello_analitico", alpha_max, 0., 12., 1);
  	modello_analitico.SetParameter(0, d);
  	
  	TApplication theApp("theApp", &argc, argv);
  	TCanvas c1;
  	grafico.Draw("AP");
  	modello_analitico.Draw("same");
  	c1.Print("grafico", "png");
  	theApp.Run();
		
	return 0;
}
