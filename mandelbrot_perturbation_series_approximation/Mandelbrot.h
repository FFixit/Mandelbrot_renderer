#pragma once

#include <cmath>
#include <fstream>
#include <stdint.h>
#include <boost/multiprecision/mpfr.hpp>
#include <algorithm>
#include <iostream>

#ifdef _DEBUG
#define MANDEL_DEBUG
#endif
//#define MANDEL_DEBUG
#define SMOOTH_COLORING

using namespace boost::multiprecision;


namespace mandelbrot
{

#ifdef MANDEL_DEBUG
#define logss std::cout
#endif

	typedef unsigned char u_char;

	std::string getTimeString();

	struct Color
	{
		u_char blue;
		u_char green;
		u_char red;
	};

	struct Point
	{
#ifdef SMOOTH_COLORING
		double test;
#endif
		uint32_t iterations;
		bool glitch;
		bool processed;
	};

	struct PointQ
	{
		size_t x;
		size_t y;
	};

	// headers for the writeBitmapHeader() function
#pragma pack(push, 1)
	struct BitmapFileHeader
	{
		int8_t b = 'B';
		int8_t m = 'M';
		uint32_t size; // size of the file in bytes
		uint16_t reserved1 = 0; // generally depends on application
		uint16_t reserved2 = 0; // ^  
		uint32_t offset; // the offset at which the pixel-array starts
	};

	struct BitmapInfoHeader
	{
		uint32_t size = 40; // size of this header in bytes
		int32_t width;
		int32_t height;
		uint16_t planes = 1; // number of color planes
		uint16_t bpp; // bits per pixel
		uint32_t compression = 0; // compression method (0 for none)
		uint32_t dataSize = 0; // size in bytes of the pixel-array
		int32_t xPixelsPerMeter = 0;
		int32_t yPixelsPerMeter = 0;
		uint32_t clrsUsed = 0; // number of colors in the palette
		uint32_t importantClrs = 0; // important colors
	};
#pragma pack(pop)


	class Mandelbrot
	{
	public:
		Mandelbrot(const Mandelbrot&) = delete;
		Mandelbrot& operator=(const Mandelbrot&) = delete;

		Mandelbrot(size_t _resx, size_t _resy, size_t _max_iterations, mpfr_float _xmin, mpfr_float _xmax, mpfr_float _ymin, mpfr_float _ymax)
			: max_iterations{ _max_iterations }
			, xmin{ _xmin }
			, xmax{ _xmax }
			, ymin{ _ymin }
			, ymax{ _ymax }
			, resx{ _resx }
			, resy{ _resy }
			, x_scale_factor{ (_xmax - _xmin) / _resx }
			, y_scale_factor{ (_ymax - _ymin) / _resy }
			, n_colors{ 256 }
		{
			generatePalette();

			colorVals = new Color[resx*resy];
			points = new Point[resx*resy];

			ref_x_arr = new double[max_iterations];
			ref_y_arr = new double[max_iterations];
			ref_x_arr_mf = new mpfr_float[max_iterations];
			ref_y_arr_mf = new mpfr_float[max_iterations];
			ref_test = new double[max_iterations];

			ar_arr = new double[max_iterations];
			ai_arr = new double[max_iterations];
			br_arr = new double[max_iterations];
			bi_arr = new double[max_iterations];
			cr_arr = new double[max_iterations];
			ci_arr = new double[max_iterations];
			dr_arr = new double[max_iterations];
			di_arr = new double[max_iterations];

			size_t maxprec = std::max(std::max(xmin.precision(), xmax.precision()), std::max(ymin.precision(), ymax.precision()));
			//mpfr_float::default_precision(maxprec);
			logss.precision(maxprec);

			logss << "XMIN: " << xmin << "\n";
			logss << "XMAX: " << xmax << "\n";
			logss << "YMIN: " << ymin << "\n";
			logss << "YMAX: " << ymax << "\n";
			logss << "MAXITER: " << max_iterations << "\n\n";
		}

		~Mandelbrot()
		{
			delete[] points;
			delete[] colorVals;
			delete[] hue_palette;

			delete[] ref_x_arr;
			delete[] ref_y_arr;
			delete[] ref_x_arr_mf;
			delete[] ref_y_arr_mf;
			delete[] ref_test;

			delete[] ar_arr;
			delete[] ai_arr;
			delete[] br_arr;
			delete[] bi_arr;
			delete[] cr_arr;
			delete[] ci_arr;
			delete[] dr_arr;
			delete[] di_arr;
		}

		bool createPerturbation();
		bool saveBitmap(std::string& filename) const;

	private:
		void generatePalette(); // this generates the color table; will be replaced by custom color table support
		void pointsToColor();

		void calculatePoints(size_t approx_iter);
		void calculateGlitchedPoints(size_t approx_iter);
		void calculatePoint(__int64& x, __int64& y, size_t& approx_iter);
		bool correctGlitches();
		size_t getArea(size_t x, size_t y, uint32_t iteration);

		void calculateReference(size_t x, size_t y);
		size_t calculateApproximation();

#ifndef MANDEL_DEBUG
		std::stringstream logss;
#endif

		Point* points;
		double* hue_palette;
		Color* colorVals;

		const size_t resx;
		const size_t resy;

		const size_t n_colors;

		const size_t max_iterations;
		const double bailout_radius{ 1 << 8 };

		const mpfr_float ymin;
		const mpfr_float ymax;
		const mpfr_float y_scale_factor;

		const mpfr_float xmin;
		const mpfr_float xmax;
		const mpfr_float x_scale_factor;

		// series approximation
		double* ar_arr;
		double* ai_arr;
		double* br_arr;
		double* bi_arr;
		double* cr_arr;
		double* ci_arr;
		double* dr_arr;
		double* di_arr;
		double* ref_x_arr;
		double* ref_y_arr;
		mpfr_float* ref_x_arr_mf;
		mpfr_float* ref_y_arr_mf;
		double* ref_test; // = (Xnr*Xnr+Xni*Xni)*10^-7
	};

	double hueLinearInterpolate(double v1, double v2, double t);
	void writeBitmapHeader(std::ofstream& ofs, size_t resx, size_t resy);


	inline Color HSVtoRGB(double h, double s = 1.0, double v = 1.0) // http://www.cs.rit.edu/~ncs/color/t_convert.html
	{
		int i, tmp;
		double f, p, q, t, vs;
		if (s == 0) {
			tmp = static_cast<int>(v * 255.0);
			return Color{ u_char(tmp), u_char(tmp), u_char(tmp) };
		}
		static const double oos = 1.0 / 60.0;
		h *= oos;			// sector 0 to 5
		i = static_cast<int>(h);
		f = h - i;			// factorial part of h
		vs = v * s;
		p = v - vs;
		q = v - vs*f;
		t = v - vs + vs*f;
		switch (i) {
		case 0:
			return Color{ u_char(p * 255.0), u_char(t * 255.0), u_char(v * 255.0) };
		case 1:
			return Color{ u_char(p * 255.0), u_char(v * 255.0), u_char(q * 255.0) };
		case 2:
			return Color{ u_char(t * 255.0), u_char(v * 255.0), u_char(p * 255.0) };
		case 3:
			return Color{ u_char(v * 255.0), u_char(q * 255.0), u_char(p * 255.0) };
		case 4:
			return Color{ u_char(v * 255.0), u_char(p * 255.0), u_char(t * 255.0) };
		default:
			return Color{ u_char(q * 255.0), u_char(p * 255.0), u_char(v * 255.0) };
		}
	}


	inline void writeBitmapHeader(std::ofstream& ofs, size_t resx, size_t resy)
	{
		BitmapFileHeader bfh;
		bfh.size = uint32_t(resy*resx * 3 + 54);
		bfh.offset = uint32_t(0x36);

		ofs.write(reinterpret_cast<char*>(&bfh), 14);

		BitmapInfoHeader bih;
		bih.width = int32_t(resx);
		bih.height = int32_t(resy);
		bih.bpp = uint16_t(24);

		ofs.write(reinterpret_cast<char*>(&bih), 40);
	}

#ifdef SMOOTH_COLORING
	inline double hueLinearInterpolate(double v1, double v2, double t)
	{
		double v2_p360 = v2 + 360.0;

		if (abs(v2 - v1) < abs(v2_p360 - v1))
		{
			return v1 + t*(v2 - v1);
		}

		return v1 + t*(v2_p360 - v1);
	}
#endif

}