#include "Mandelbrot.h"
#include <iostream>
#include <algorithm>
#include <stack>
#include <omp.h>

using namespace std;

namespace mandelbrot
{

	void Mandelbrot::generatePalette()
	{
		const double START = 0.0;
		const double END = 360.0;

		hue_palette = new double[n_colors];

		for (size_t i = 0; i < n_colors; i++)
		{
			hue_palette[i] = ((END - START) / n_colors) * i + START;
		}
	}


	bool Mandelbrot::createPerturbation()
	{
		logss << getTimeString() << ": Setting up.\n";

		calculateReference(resx / 2, resy / 2);
		size_t approx_iter = calculateApproximation();
		calculatePoints(approx_iter);

		for (size_t i = 0; i < 15; i++) // tries to correct glitches up to 15 times, breaks if there are no more blobs with area greater than 5
		{
			logss << getTimeString() << ": Started to correct glitches (" << i + 1 << ").\n";

			if (!correctGlitches())
			{
				logss << getTimeString() << ": No more glitched regions with area over 5.\n";
				break;
			}

			logss << getTimeString() << ": Finished correcting glitches.\n\n";
		}

		pointsToColor();
		logss << getTimeString() << ": Render complete.\n";

		return true;
	}


	void Mandelbrot::calculatePoints(size_t approx_iter)
	{
		logss << getTimeString() << ": Starting to calculate the picture with " << omp_get_num_threads() << " threads.\n";

#pragma omp parallel for schedule(dynamic)
		for (__int64 y = 0; y < (__int64)resy; y++)
			for (__int64 x = 0; x < (__int64)resx; x++)
			{
				calculatePoint(x, y, approx_iter);
			}

		logss << getTimeString() << ": Finished calculating the picture.\n\n";
	}


	void Mandelbrot::calculateGlitchedPoints(size_t approx_iter)
	{
		logss << getTimeString() << ": Starting to calculate the glitched pixels with " << omp_get_num_threads() << " threads.\n";

#pragma omp parallel for schedule(dynamic)
		for (__int64 y = 0; y < (__int64)resy; y++)
			for (__int64 x = 0; x < (__int64)resx; x++)
			{
				if (points[y*resx + x].glitch)
					calculatePoint(x, y, approx_iter);
			}

		logss << getTimeString() << ": Finished calculating the glitched pixels.\n";
	}


	void Mandelbrot::calculatePoint(__int64& x, __int64& y, size_t& approx_iter)
	{
		double D0r, D0i;	// Δ0
		double Dnr, Dni;	// Δn
		double Xnr, Xni;	// ref-point
		double Yr, Yi;		// current point
		double d0r2, d0i2;	// Δ0^2
		double d0r3, d0i3;	// Δ0^3
		double d0r4, d0i4;	// Δ0^4
#ifdef SMOOTH_COLORING
		double r_test{};
#endif

		mpfr_float dx_scaled = (x_scale_factor*x + xmin) - ref_x_arr_mf[0];
		D0r = (double)dx_scaled;

		mpfr_float dy_scaled = (y_scale_factor*y + ymin) - ref_y_arr_mf[0];
		D0i = (double)dy_scaled;

		size_t i = approx_iter - 1;

		d0r2 = D0r*D0r - D0i*D0i;
		d0i2 = 2.0*D0r*D0i;
		d0r3 = d0r2*D0r - d0i2*D0i;
		d0i3 = d0i2*D0r + d0r2*D0i;
		d0r4 = d0r2*d0r2 - d0i2*d0i2;
		d0i4 = 2.0*d0r2*d0i2;

		// An * Δ0
		double tr = ar_arr[i] * D0r - ai_arr[i] * D0i;
		double ti = ai_arr[i] * D0r + ar_arr[i] * D0i;
		// + Bn * Δ0²
		tr += br_arr[i] * d0r2 - bi_arr[i] * d0i2;
		ti += bi_arr[i] * d0r2 + br_arr[i] * d0i2;
		// + Cn * Δ0³
		tr += cr_arr[i] * d0r3 - ci_arr[i] * d0i3;
		ti += ci_arr[i] * d0r3 + cr_arr[i] * d0i3;
		// + Dn * Δ0^4
		tr += dr_arr[i] * d0r4 - di_arr[i] * d0i4;
		ti += dr_arr[i] * d0i4 + di_arr[i] * d0r4;

		Dnr = tr;
		Dni = ti;

		bool glitch = false;
		for (; i < max_iterations; i++) // iterations
		{
			Xnr = ref_x_arr[i];
			Xni = ref_y_arr[i];

			Yr = Xnr + Dnr;
			Yi = Xni + Dni;

			double test = (Yr*Yr + Yi*Yi);
			double* test_ptr = &ref_test[i];

			if (test > 4.0)
			{
#ifdef SMOOTH_COLORING
				r_test = test;
#endif
				break;
			}

			if (test < *test_ptr)
			{
#ifdef SMOOTH_COLORING
				r_test = test;
#endif
				glitch = true;
				break;
			}

			double d0r_temp = 2.0*(Xnr*Dnr - Xni*Dni) + (Dnr*Dnr - Dni*Dni) + D0r;
			double d0i_temp = 2.0*(Dni*(Xnr + Dnr) + Xni*Dnr) + D0i;

			Dnr = d0r_temp;
			Dni = d0i_temp;
		}

#ifdef SMOOTH_COLORING
		*(points + y*resx + x) = Point{ r_test, (uint32_t)i, glitch, false };
#else
		// some hacking to replace the implicit pointer multiplication with a shift operation since sizeof(Point) == 2³
		*reinterpret_cast<Point*>(((y*resx + x) << 3) + reinterpret_cast<size_t>(points)) = Point{ (uint32_t)i, glitch, false };
#endif
	}


	void Mandelbrot::calculateReference(size_t x, size_t y)
	{
		logss << getTimeString() << ": Starting to calculate reference point at X:" << x << " Y:" << y << ".\n";

		mpfr_float x_scaled = xmin + x*x_scale_factor;
		mpfr_float y_scaled = ymin + y*y_scale_factor;

		mpfr_float ref_zx = 0.0;
		mpfr_float ref_zy = 0.0;

		for (size_t it = 0; it < max_iterations; it++)
		{
			mpfr_float xtemp = ref_zx*ref_zx - ref_zy*ref_zy + x_scaled;
			ref_zy = 2.0 * ref_zx*ref_zy + y_scaled;
			ref_zx = xtemp;

			ref_x_arr[it] = (double)ref_zx;
			ref_y_arr[it] = (double)ref_zy;
			ref_x_arr_mf[it] = ref_zx;
			ref_y_arr_mf[it] = ref_zy;
			ref_test[it] = (ref_x_arr[it] * ref_x_arr[it] + ref_y_arr[it] * ref_y_arr[it])*0.0000001;
		}

		logss << getTimeString() << ": Finished calulating reference point.\n";
	}


	size_t Mandelbrot::calculateApproximation()
	{
		logss << getTimeString() << ": Starting to calculate approximation.\n";

		double ar{ 1 }, ai{ 0 };
		double br{ 0 }, bi{ 0 };
		double cr{ 0 }, ci{ 0 };
		double dr{ 0 }, di{ 0 };

		mpfr_float* comp_px = new mpfr_float[4];
		mpfr_float* comp_py = new mpfr_float[4];
		comp_px[0] = xmax;
		comp_py[0] = ymax;
		comp_px[1] = xmax;
		comp_py[1] = ymin;
		comp_px[2] = xmin;
		comp_py[2] = ymax;
		comp_px[3] = xmin;
		comp_py[3] = ymin;
		double dif = 0.001;
		double* comp_d0r = new double[4];
		double* comp_d0i = new double[4];
		double* comp_dnr_eq1 = new double[4];
		double* comp_dni_eq1 = new double[4];
		double* comp_dnr_eq2 = new double[4];
		double* comp_dni_eq2 = new double[4];
		
		for (size_t i = 0; i < 4; i++)
		{
			mpfr_float comp_d0r_t = comp_px[i] - ref_x_arr_mf[0];
			mpfr_float comp_d0i_t = comp_py[i] - ref_y_arr_mf[0];
			comp_d0r[i] = (double)comp_d0r_t;
			comp_d0i[i] = (double)comp_d0i_t;
			comp_dnr_eq1[i] = comp_d0r[i];
			comp_dni_eq1[i] = comp_d0i[i];
			comp_dnr_eq2[i] = comp_d0r[i];
			comp_dni_eq2[i] = comp_d0i[i];
		}

		size_t iteration;
		for (iteration = 0; iteration < max_iterations; iteration++)
		{
			ar_arr[iteration] = ar;
			ai_arr[iteration] = ai;
			br_arr[iteration] = br;
			bi_arr[iteration] = bi;
			cr_arr[iteration] = cr;
			ci_arr[iteration] = ci;
			dr_arr[iteration] = dr;
			di_arr[iteration] = di;

			// term 1
			mpfr_float artemp = 2.0*(ref_x_arr_mf[iteration] * ar - ref_y_arr_mf[iteration] * ai) + 1.0;
			mpfr_float aitemp = 2.0*(ref_y_arr_mf[iteration] * ar + ref_x_arr_mf[iteration] * ai);
			// term 2
			mpfr_float brtemp = 2.0*(ref_x_arr_mf[iteration] * br - ref_y_arr_mf[iteration] * bi) + ar*ar - ai*ai;
			mpfr_float bitemp = 2.0*(ref_y_arr_mf[iteration] * br + ref_x_arr_mf[iteration] * bi) + 2.0*ar*ai;
			// term 3
			mpfr_float crtemp = 2.0*(ref_x_arr_mf[iteration] * cr - ref_y_arr_mf[iteration] * ci) + 2.0*(ar*br - ai*bi);
			mpfr_float citemp = 2.0*(ref_y_arr_mf[iteration] * cr + ref_x_arr_mf[iteration] * ci) + 2.0*(ai*br + ar*bi);
			// term 4
			mpfr_float drtemp = 2.0*(ref_x_arr_mf[iteration] * dr - ref_y_arr_mf[iteration] * di) + 2.0*(ar*cr - ai*ci) - br*br - bi*bi;
			mpfr_float ditemp = 2.0*(ref_y_arr_mf[iteration] * dr + ref_x_arr_mf[iteration] * di) + 2.0*(ai*cr + ar*ci) + 2.0*br*bi;

			size_t j = 0;
			for (; j < 4; j++)
			{
				// eq2
				double d0r2 = comp_d0r[j] * comp_d0r[j] - comp_d0i[j] * comp_d0i[j];
				double d0i2 = 2.0*comp_d0r[j] * comp_d0i[j];
				double d0r3 = d0r2*comp_d0r[j] - d0i2*comp_d0i[j];
				double d0i3 = d0i2*comp_d0r[j] + d0r2*comp_d0i[j];
				double d0r4 = d0r2*d0r2 - d0i2*d0i2;
				double d0i4 = 2.0*d0r2*d0i2;

				comp_dnr_eq2[j] = ar * comp_d0r[j] - ai * comp_d0i[j];
				comp_dni_eq2[j] = ai * comp_d0r[j] + ar * comp_d0i[j];
				// + Bn * Δ0²
				comp_dnr_eq2[j] += br * d0r2 - bi * d0i2;
				comp_dni_eq2[j] += bi * d0r2 + br * d0i2;
				// + Cn * Δ0³
				comp_dnr_eq2[j] += cr * d0r3 - ci * d0i3;
				comp_dni_eq2[j] += ci * d0r3 + cr * d0i3;
				// + Dn * Δ0^4
				comp_dnr_eq2[j] += dr * d0r4 - di * d0i4;
				comp_dni_eq2[j] += dr * d0i4 + di * d0r4;


				double diff_r = ((comp_dnr_eq2[j] - comp_dnr_eq1[j]) / comp_dnr_eq1[j]);
				if ((diff_r > dif) || (diff_r < -dif))
					break;

				double diff_i = ((comp_dni_eq2[j] - comp_dni_eq1[j]) / comp_dni_eq1[j]);
				if ((diff_i > dif) || (diff_i < -dif))
					break;

				// eq1
				double dnr_temp = 2.0*(ref_x_arr[iteration] * comp_dnr_eq1[j] - ref_y_arr[iteration] * comp_dni_eq1[j]) + (comp_dnr_eq1[j] * comp_dnr_eq1[j] - comp_dni_eq1[j] * comp_dni_eq1[j]) + comp_d0r[j];
				double dni_temp = 2.0*(ref_x_arr[iteration] * comp_dni_eq1[j] + ref_y_arr[iteration] * comp_dnr_eq1[j]) + 2.0*comp_dnr_eq1[j] * comp_dni_eq1[j] + comp_d0i[j];
				comp_dnr_eq1[j] = dnr_temp;
				comp_dni_eq1[j] = dni_temp;
			}
			if (j < 4)
				break;

			ar = (double)artemp;
			ai = (double)aitemp;
			br = (double)brtemp;
			bi = (double)bitemp;
			cr = (double)crtemp;
			ci = (double)citemp;
			dr = (double)drtemp;
			di = (double)ditemp;
		}

		delete[] comp_px;
		delete[] comp_py;
		delete[] comp_d0r;
		delete[] comp_d0i;
		delete[] comp_dnr_eq1;
		delete[] comp_dni_eq1;
		delete[] comp_dnr_eq2;
		delete[] comp_dni_eq2;

		logss << getTimeString() << ": Finished calcluating approximation. Skipping " << iteration - 1 << " iterations.\n";
		return iteration;
	}


	// this tries to correct the biggest blob by calculating a new reference point in the estimated center of the blob
	// returns false if there are no blobs with an area bigger than 5
	bool Mandelbrot::correctGlitches()
	{
		size_t curr_area{}, ref_x{}, ref_y{}, cnt{};
		for (size_t y = 1; y < resy - 1; y++)
			for (size_t x = 1; x < resx - 1; x++)
			{
				size_t area = getArea(x, y, points[y*resx + x].iterations);
				if (area > curr_area)
				{
					curr_area = area;
					ref_x = x;
					ref_y = y;
				}

				if (points[y*resx + x].glitch)
					cnt++;
			}

		for (size_t i = 0; i < resx*resy; i++)
			points[i].processed = false;

		logss << getTimeString() << ": Detected " << cnt << " glitched pixels.\n";

		if (curr_area < 6)
			return false;

		__int64 maxdist{};
		size_t rx{}, ry{};
		for (__int64 y = 1; y < (__int64)resy - 1; y++)
			for (__int64 x = 1; x < (__int64)resx - 1; x++)
			{
				if (points[y*resx + x].iterations != points[ref_y*resx + ref_x].iterations)
					continue;

				__int64 t = 0, to, c, ct;
				ct = c = 0;
				for (to = 0; x - to >= 0 && (points[y*(__int64)resx + (x - to)].iterations == points[ref_y*resx + ref_x].iterations); to++)
				{
					t++;
					ct++;
				}
				for (to = 0; x + to < (__int64)resx && points[y*resx + (x + to)].iterations == points[ref_y*resx + ref_x].iterations; to++)
				{
					t++;
					ct--;
				}
				c += (ct < 0 ? -ct : ct);
				ct = 0;

				for (to = 0; y - to >= 0 && (points[(y - to)*resx + x].iterations == points[ref_y*resx + ref_x].iterations); to++)
				{
					t++;
					ct++;
				}
				for (to = 0; y + to < (__int64)resy && points[(y + to)*resx + x].iterations == points[ref_y*resx + ref_x].iterations; to++)
				{
					t++;
					ct--;
				}
				c += (ct < 0 ? -ct : ct);
				ct = 0;

				for (to = 0; (y - to >= 0) && (x - to >= 0) && points[(y - to)*resx + (x - to)].iterations == points[ref_y*resx + ref_x].iterations; to++)
				{
					t++;
					ct++;
				}
				for (to = 0; (y + to < (__int64)resy) && (x + to >= (__int64)resx) && points[(y + to)*resx + (x + to)].iterations == points[ref_y*resx + ref_x].iterations; to++)
				{
					t++;
					ct--;
				}
				c += (ct < 0 ? -ct : ct);
				ct = 0;

				for (to = 0; (y + to < (__int64)resy) && (x - to >= 0) && points[(y + to)*resx + (x - to)].iterations == points[ref_y*resx + ref_x].iterations; to++)
				{
					t++;
					ct++;
				}
				for (to = 0; (y - to >= 0) && (x + to < (__int64)resx) && points[(y - to)*resx + (x + to)].iterations == points[ref_y*resx + ref_x].iterations; to++)
				{
					t++;
					ct--;
				}
				c += (ct < 0 ? -ct : ct);
				ct = 0;

				if (maxdist < t)
				{
					maxdist = t;
					rx = x;
					ry = y;
				}
			}

		ref_x = rx;
		ref_y = ry;
		calculateReference(ref_x, ref_y);

		size_t approx_iter = calculateApproximation();
		calculateGlitchedPoints(approx_iter);
		return true;
	}


	// stack based flood-fill algorithm to determine the area if a glitched region around the point P(x, y)
	// this also finds the minimum/maximum x/y-values of the area
	size_t Mandelbrot::getArea(size_t x, size_t y, uint32_t iteration)
	{
		size_t area{};
		stack<PointQ> queue;
		queue.push(PointQ{ x, y });

		while (queue.size() > 0)
		{
			x = queue.top().x;
			y = queue.top().y;

			size_t left = 0;
			while ((__int64)x - (__int64)left >= 0 && !points[y*resx + (x - left)].processed && points[y*resx + (x - left)].iterations == iteration && points[y*resx + (x - left)].glitch)
			{
				area++;
				points[y*resx + (x - left)].processed = true;

				if (y + 1 < resy)
				{
					queue.push(PointQ{ x - left, y + 1 });
				}
				if ((__int64)y - (__int64)1 >= 0)
				{
					queue.push(PointQ{ x - left, y - 1 });
				}

				left++;
			}

			size_t right = 1;
			while (x + right < resx && !points[y*resx + (x + right)].processed && points[y*resx + (x + right)].iterations == iteration && points[y*resx + (x + right)].glitch)
			{
				area++;
				points[y*resx + (x + right)].processed = true;

				if (y + 1 < resy)
				{
					queue.push(PointQ{ x + right, y + 1 });
				}
				if ((__int64)y - (__int64)1 >= 0)
				{
					queue.push(PointQ{ x + right, y - 1 });
				}

				right++;
			}

			queue.pop();
		}

		return area;
	}


	void Mandelbrot::pointsToColor()
	{
		for (size_t y = 0; y < resy; y++)
			for (size_t x = 0; x < resx; x++)
			{
#ifndef MANDEL_DEBUG
#ifdef SMOOTH_COLORING
				static const double LOG2 = log(2);

				if (points[y*resx + x].iterations == max_iterations)
					*(colorVals + y*resx + x) = Color{ 0, 0, 0 };
				else
				{
					double si = points[y*resx + x].iterations + 1 - (log(0.5 * log(points[y*resx + x].test)) / LOG2);
					double h1 = hue_palette[static_cast<int>(floor(si)) % n_colors];
					double h2 = hue_palette[(static_cast<int>(floor(si)) + 1) % n_colors];
					double t = fmod(si, 1);

					*(colorVals + y*resx + x) = HSVtoRGB(hueLinearInterpolate(h1, h2, t));
				}
#else
				if (points[y*resx + x].iterations == max_iterations)
					*(colorVals + y*resx + x) = Color{ 0, 0, 0 };
				else
					*(colorVals + y*resx + x) = HSVtoRGB(hue_palette[points[y*resx + x].iterations % n_colors]);
#endif
#else
				if (points[y*resx + x].glitch)
					*(colorVals + y*resx + x) = Color{ 255, 255, 255 };
				else
					*(colorVals + y*resx + x) = Color{ 0, 0, 0 };
#endif					
			}
	}


	bool Mandelbrot::saveBitmap(string& filename) const
	{
		if (colorVals == nullptr)
			return false;

		ofstream bmfs(filename + string(".bmp"), ios::out | ios::binary | ios::trunc);

		if (!bmfs)
			return false;

		writeBitmapHeader(bmfs, resx, resy);
		bmfs.write(reinterpret_cast<char*>(colorVals), 3 * resx*resy);
		bmfs.close();

#ifndef MANDEL_DEBUG
		ofstream logfs(filename + string(".log"), ios::out | ios::trunc);

		if (!logfs)
			return false;

		logfs.write(logss.str().c_str(), logss.str().size());
		logfs.close();
#endif
		return true;
	}


	string getTimeString()
	{
		char buffer[80];
		time_t rawtime;
		struct tm timeinfo;

		time(&rawtime);
		localtime_s(&timeinfo, &rawtime);

		strftime(buffer, 80, "%d.%m.%Y %I:%M:%S", &timeinfo);
		return string("[") + string(buffer) + string("]");
	}

}