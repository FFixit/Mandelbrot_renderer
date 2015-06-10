#include "Main.h"
#include <iostream>
#include <fstream>
#include "Mandelbrot.h"
#include <limits>
#include <xmmintrin.h>
#include <complex>

using namespace std;
using namespace mandelbrot;




int main()
{
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	_mm_setcsr(_mm_getcsr() | 0x8040);

	/*
	mpfr_float xmin("-0.6367543465823899786437398318");
	mpfr_float xmax("-0.6367543465823899786437398320");
	mpfr_float ymin("0.68503129708367730147608783913");
	mpfr_float ymax("0.68503129708367730147608783915");



	mpfr_float xmin("-0.9223327810370947027656057193752719767735");
	mpfr_float xmax("-0.9223327810370947027656057193752719747735");
	mpfr_float ymin("0.3102598350874576432708737495917724826010");
	mpfr_float ymax("0.3102598350874576432708737495917724846010");




	mpfr_float xmin("0.01343887053201212902835491900");
	mpfr_float xmax("0.01343887053201212902837491900");
	mpfr_float ymin("0.65561421876946506225131002766");
	mpfr_float ymax("0.65561421876946506225133002766");



	mpfr_float xmin("-0.0056306862431773201255495207714275917884892103802234");
	mpfr_float xmax("0.032253864520761570210093726941663773934633200252854");
	mpfr_float ymin("0.63046253084936263298371639446880296399787487533285");
	mpfr_float ymax("0.66834708161330152331935964218189432972099728596593");

	mpfr_float xmin("0.01343887053201212902736491900");
	mpfr_float xmax("0.01343887053201212902936491900");
	mpfr_float ymin("0.65561421876946506225032002766");
	mpfr_float ymax("0.65561421876946506225232002766");



	mpfr_float xmin = -2.0;
	mpfr_float xmax = 1.0;
	mpfr_float ymin = -1.5;
	mpfr_float ymax = 1.5;
	*/



	system("mkdir \"animation\"");
	createAnimation();
	//panning();

	/*
	Mandelbrot mb(500, 500, 200000, xmin, xmax, ymin, ymax);
	mb.createPerturbation();
	mb.saveBitmap(string("mandelbrot"));
	*/

	system("pause");
}

#define N_FRAMES 5400
//mpfr_float P_X("0.01343887053201212902836491900");
//mpfr_float P_Y("0.65561421876946506225132002766");


//mpfr_float P_X("-0.1530676364065576286904807509243619153470203124158419432007797242524298416437269823619332285088962814515788961084586756693716720988397345635379187862366819034541545825442133825796657905172967598241623630174261227911955514660715572434320637912705787076679687");
//mpfr_float P_Y("1.02985410229893199881242903853907937978390247709547828169063255987514732525692441319928174495315016245522160632413811648857159519629570354825037579665366858956778548725156897337727664245912188390909350285813038025836887536749418909977401147654788781861858655");

#define RES_X 1280
#define RES_Y 720
#define MIN_ITER 200
#define MAX_ITER 600000
void createAnimation()
{
	mpfr_float::default_precision(270);

	mpfr_float P_X;
	P_X.precision(1000);
	P_X = mpfr_float_1000("-0.1530676364065576286904807509243619153470203124158419432007797242524298416437269823619332285088962814515788961084586756693716720988397345635379187862366819034541545825442133825796657905172967598241623630174261227911955514660715572434320637912705787076679687");


	mpfr_float P_Y;
	P_Y.precision(1000);
	P_Y = mpfr_float_1000("1.02985410229893199881242903853907937978390247709547828169063255987514732525692441319928174495315016245522160632413811648857159519629570354825037579665366858956778548725156897337727664245912188390909350285813038025836887536749418909977401147654788781861858655");


	mpfr_float X_ENDRANGE;
	X_ENDRANGE.precision(1000);
	X_ENDRANGE = mpfr_float_1000(8.051435920)*pow(mpfr_float_1000(10), -233);



	mpfr_float zoom_base = pow(((P_X + X_ENDRANGE) - (P_X - X_ENDRANGE)) / (2.0 - (-2.0)), mpfr_float(1.0) / N_FRAMES);

	for (size_t i = 5359; i < N_FRAMES; i++)
	{
		mpfr_float factor = pow(zoom_base, i);

		mpfr_float xmin = (P_X)-((P_X)-(P_X - 2.0))*factor;
		mpfr_float xmax = (P_X)+((P_X + 2.0) - (P_X))*factor;

		mpfr_float ymin = (P_Y)-((P_Y)-((P_Y - 2.0)*(9.0 / 16.0)))*factor;
		mpfr_float ymax = (P_Y)+((P_Y + 2.0)*(9.0 / 16.0) - (P_Y))*factor;
		
		cout << "rendering frame " << i + 1 << "...\n";

		size_t iter = MIN_ITER + ((MAX_ITER - MIN_ITER) / N_FRAMES)*i;

		Mandelbrot mb(RES_X, RES_Y, iter, xmin, xmax, ymin, ymax);
		mb.createPerturbation();
		mb.saveBitmap(string("animation\\mandelbrot_frame_") + to_string(i + 1));
	}
}


void panning()
{
	mpfr_float::default_precision(270);

	mpfr_float P_X;
	P_X.precision(1000);
	P_X = mpfr_float_1000("-0.1530676364065576286904807509243619153470203124158419432007797242524298416437269823619332285088962814515788961084586756693716720988397345635379187862366819034541545825442133825796657905172967598241623630174261227911955514660715572434320637912705787076679687");


	mpfr_float P_Y;
	P_Y.precision(1000);
	P_Y = mpfr_float_1000("1.02985410229893199881242903853907937978390247709547828169063255987514732525692441319928174495315016245522160632413811648857159519629570354825037579665366858956778548725156897337727664245912188390909350285813038025836887536749418909977401147654788781861858655");


	mpfr_float d0i, d0r;
	d0r = -P_X;
	d0i = -P_Y;

	for (size_t i = 0; i < N_FRAMES; i++)
	{
		double x = (1.0 / (N_FRAMES - 1))*i;
		double fx = x;// *(cbrt(x - 1.0) + 1);
		cout << fx << endl;

		mpfr_float dni, dnr;
		dnr = d0r*fx;
		dni = d0i*fx;

		mpfr_float pnr, pni;
		pnr = P_X + dnr;
		pni = P_Y + dni;

		mpfr_float xmin = (pnr)-((pnr)-(pnr - 2.0));
		mpfr_float xmax = (pnr)+((pnr + 2.0) - (pnr));

		mpfr_float ymin = (pni)-((pni)-((pni - 2.0)*(9.0 / 16.0)));
		mpfr_float ymax = (pni)+((pni + 2.0)*(9.0 / 16.0) - (pni));


		Mandelbrot mb(RES_X, RES_Y, 200, xmin, xmax, ymin, ymax);
		mb.createPerturbation();
		mb.saveBitmap(string("mandelbrot_moving_frame__") + to_string(i + 1));
	}


}