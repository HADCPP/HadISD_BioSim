#pragma once

#include "station.h"
#include "python_function.h"
#include "Utilities.h"

#include<vector>
#include <algorithm>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <valarray>
#include <dlib\optimization.h>

typedef dlib::matrix<double, 1, 1> input_vector;
typedef dlib::matrix<double, 2, 1> parameter_vector;
typedef dlib::matrix<double, 3, 1> parameter_vector1;
typedef dlib::matrix<double, 5, 1> parameter_vector2;


namespace UTILS
{
	void  month_starts(boost::gregorian::date start, boost::gregorian::date end, std::vector<int>& month_locs);
	inline std::vector<std::pair<int, int>> month_starts_in_pairs(boost::gregorian::date start, boost::gregorian::date end);
	std::valarray<bool> create_fulltimes(CStation& station, std::vector<std::string>& var_list, boost::gregorian::date start,
		boost::gregorian::date end, std::vector<std::string>& opt_var_list, bool do_input_station_id = true,
bool do_qc_flags = true, bool do_flagged_obs = true);
	/*
	Return reporting accuracy & reporting frequency for variable

	: param obj st_var :	CStation variable object
	: param datetime start : start of data series
	: param datatime end :  end of data series

	: returns :
	 reporting_stats - Nx2 array, one pair for each month
	'''
	*/
	std::valarray<std::pair<float, float>> monthly_reporting_statistics(CMetVar& st_var, boost::gregorian::date start, boost::gregorian::date end);
	
	/*
	 Apply these flags to all variables

    :param object station: the station object to be processed
    :param list all_variables: the variables where the flags are to be applied
    :param list flag_col: which column in the qc_flags array to work on
    :param file logfile: logfile to store outputs
    :param str test_name: test name for printing/loggin
        :returns:
	*/
	void apply_flags_all_variables(CStation& station, const std::vector<std::string>& full_variable_list, int flag_col, std::ofstream & logfile, std::string test);
	void apply_windspeed_flags_to_winddir(CStation& station);
	void append_history(CStation& station, std::string text);
	CMaskedArray<float>  apply_filter_flags(CMetVar& st_var);
	void print_flagged_obs_number(std::ofstream& logfile, std::string test, std::string variable, int nflags);
	inline float __cdecl MyApplyRoundFunc(float n){return std::floor(n);}
	inline float __cdecl MyApplyCeilFunc(float n){ return std::ceil(n); }
	inline float __cdecl MyApplyCastFloat(int n){ return static_cast<float>(n); }
	template<typename T>
	inline T __cdecl MyApplyAbsFunc(T n){return std::abs(n);}
	inline float __cdecl Log10Func(float n){return std::log10f(n);}
	/*
	Following reporting_accuracy.pro method.
	Uses histogram of remainders to look for special values

	:param array indata: masked array
	:output: frequency - reporting frequency of data
	*/
	inline float reporting_accuracy(CMaskedArray<float>& indata, bool winddir = false);
	inline float reporting_accuracy(varrayfloat& good_values, bool winddir=false);
	
	/*
	Uses histogram of remainders to look for special values

	: param array good_values :valarray (masked array.compressed())
	: param bool winddir : true if processing wind directions

	: output : resolution - reporting accuracy(resolution) of data
	*/
	inline float reporting_frequency(CMaskedArray<float>& good_indata);

	void create_bins(varrayfloat& indata, float binwidth, varrayfloat &bins, varrayfloat  &bincenters);
	void create_bins(std::valarray<CMaskedArray<float>>& indata, float binwidth, varrayfloat& bins, varrayfloat& bincenters);
	template<typename T>
	inline T idl_median(std::valarray<T> indata)
	{
		std::nth_element(std::begin(indata), std::begin(indata) + indata.size() / 2, std::end(indata));
		return indata[indata.size() / 2];
	}
	template<typename T, typename S>
	inline T Cast(const S &data){return static_cast<T>(data);}
	template<typename T>
	inline T Cast(const std::string &data)
	{
		string s = typeid(T).name();
		if (s == "int")
			return static_cast<T>(std::atoi(data.c_str()));

		else
			return static_cast<T>(std::atof(data.c_str()));

	}

	/*
		Plot histogram on log-y scale and fit 1/x decay curve to set threshold
		 :param array indata: input data to bin up
		:param int binmin: minimum bin value
		:param int binwidth: bin width
		:param str line_label: label for plotted histogram
		:param float old_threshold: (spike) plot the old threshold from IQR as well

		:returns:   threshold value
	   */
	inline double model(const input_vector&  input, const parameter_vector& params)
	{
		const double p0 = params(0);
		const double p1 = params(1);

		const double i0 = input(0);
		const double temp = p0*i0 + p1;

		return temp;
	}
	
	// ----------------------------------------------------------------------------------------

	// This function is the "residual" for a least squares problem.   It takes an input/output
	// pair and compares it to the output of our model and returns the amount of error.  The idea
	// is to find the set of parameters which makes the residual small on all the data pairs.
	// Source : dlib documentation
	inline double residual(const std::pair<input_vector, double>& data, const parameter_vector& params)
	{
		return model(data.first, params) - data.second;
	}

	// ----------------------------------------------------------------------------------------
	// This function is the derivative of the residual() function with respect to the parameters.
	inline parameter_vector residual_derivative(const std::pair<input_vector, double>& data, const parameter_vector& params)
	{
		parameter_vector der(2);

		const double p0 = params(0);
		const double p1 = params(1);


		const double x = data.first(0);

		const double temp = p0*x + p1;

		der(0) = x; // i0 * 2 * temp;
		der(1) = 1; // 2 * temp;

		return der;
	}
	
	/*
	decay function for line fitting
	p[0] = intercept
	p[1] = slope
	*/
	inline varrayfloat linear(const varrayfloat& X, parameter_vector& p)
	{
		return float(p(0)) * X + float(p(1));
	}
	
	template<typename T>
	inline float get_critical_values(std::vector<T>& indata, T binmin = T(0), T binwidth = T(1), float old_threshold = 0.)
	{
		float threshold;
		if (indata.size() > 0)  //sort indata ?
		{
			std::vector<T> dummy = np_abs(indata);
			dummy = np_ceil<T>(dummy);
			varrayfloat full_edges = arange<float>(float(3 * (*std::max_element(dummy.begin(), dummy.end()))), T(binmin), T(binwidth));
			dummy.clear();
			varrayfloat full_hist = histogram<int>(np_abs(indata), full_edges);
			if (full_hist.size() > 1)
			{
				// use only the central section (as long as it's not just 2 bins)
				int i = 0;
				int limit = 0;
				while (limit < 2)
				{
					if (npwhere<float>(full_hist, "=", float(0)).size()>0)
					{
						limit = npwhere<float>(full_hist, "=", float(0))[i];  //where instead of argwhere??
						i++;
					}
					else
					{
						limit = full_hist.size();
						break;
					}
				}
				varrayfloat edges = full_edges[std::slice(0, limit, 1)];
				varrayfloat hist = full_hist[std::slice(0, limit, 1)];
				hist=hist.apply(Log10Func);
				//Initialiser les parametres

				// DATA (X,Y)


				std::vector<std::pair<input_vector, double>> data;
				input_vector input;
				for (size_t i = 0; i < edges.size(); ++i)
				{
					input(0) = edges[i];
					data.push_back(make_pair(input, hist[i]));
				}
				// Now let's use the solve_least_squares_lm() routine to figure out what the
				// parameters are based on just the data.
				//Use the Levenberg-Marquardt method to determine the parameters which
				// minimize the sum of all squared residuals
				parameter_vector  fit;
				
				fit(0) = 1;
				fit(1) = hist.max();
				
				dlib::solve_least_squares_lm(dlib::objective_delta_stop_strategy(1e-7).be_verbose(),
					residual,
					residual_derivative,
					data,
					fit);
				varrayfloat fit_curve = linear(full_edges, fit);
				
				if (fit(0) < 0)
				{
					//in case the fit has a positive slope
					//where does fit fall below log10(-0.1)
					if (npwhere<float>(fit_curve, "<", float(-1)).size()>0)
					{
						size_t fit_below_point1 = npwhere<float>(fit_curve, "<", float(-1))[0];
						varrayfloat dummy = full_hist[slice(fit_below_point1, full_hist.size() - fit_below_point1, 1)];
						size_t first_zero_bin = npwhere<float>(dummy, "=", 0)[0] + 1;
						threshold = binwidth*(fit_below_point1 + first_zero_bin);
					}
					else
					{
						//too shallow a decay - use default maximum
						threshold = full_hist.size();
					}
					//find first empty bin after that
				}
				else threshold = full_hist.size();

			}
			else
				threshold = *std::max_element(indata.begin(), indata.end()) + binwidth;
		}
		else
			threshold = *std::max_element(indata.begin(), indata.end()) + binwidth;

		return threshold;
	}

	template<typename T>
	inline float get_critical_values(std::valarray<T>& indata, T binmin = T(0), T binwidth=T(1), float old_threshold = 0.)
	{

		float threshold;
		if (indata.size() > 0)  //sort indata ?
		{
			std::valarray<T> dummy = indata.apply(MyApplyAbsFunc);
			dummy = np_ceil<T>(dummy);
			varrayfloat full_edges = arange<float>(float(3 * (*std::max_element( std::begin(dummy), std::end(dummy) ) ) ), T(binmin), T(binwidth));
			dummy.free();
			varrayfloat full_hist = histogram<T>(indata, full_edges);
			if (full_hist.size() > 1)
			{
				// use only the central section (as long as it's not just 2 bins)
				int i = 0;
				int limit = 0;
				while (limit < 2)
				{
					if (npwhere<float>(full_hist, "=", float(0)).size()>0)
					{
						limit = npwhere<float>(full_hist, "=", float(0))[i];  
						i++;
					}
					else
					{
						limit = full_hist.size();
						break;
					}
				}
				varrayfloat edges = full_edges[std::slice(0, limit, 1)];
				varrayfloat hist = full_hist[std::slice(0, limit, 1)];
				hist = hist.apply(Log10Func);
				//Initialiser les parametres

				// DATA (X,Y)


				std::vector<std::pair<input_vector, double>> data;
				input_vector input;
				for (size_t i = 0; i < edges.size(); ++i)
				{
					input(0) = edges[i];
					data.push_back(make_pair(input, hist[i]));
				}
				// Now let's use the solve_least_squares_lm() routine to figure out what the
				// parameters are based on just the data.
				//Use the Levenberg-Marquardt method to determine the parameters which
				// minimize the sum of all squared residuals
				parameter_vector  fit;

				fit(0) = 1;
				fit(1) = hist.max();

				dlib::solve_least_squares_lm(dlib::objective_delta_stop_strategy(1e-7).be_verbose(),
					residual,
					residual_derivative,
					data,
					fit);
				varrayfloat fit_curve = linear(full_edges, fit);

				if (fit(0) < 0)
				{
					//in case the fit has a positive slope
					//where does fit fall below log10(-0.1)
					if (npwhere<float>(fit_curve, "<", float(-1)).size() > 0)
					{
						size_t fit_below_point1 = npwhere<float>(fit_curve, "<", float(-1))[0];
						varrayfloat dummy = full_hist[std::slice(fit_below_point1, full_hist.size() - fit_below_point1, 1)];
						size_t first_zero_bin = npwhere<float>(dummy, "=", 0)[0] + 1;
						threshold = binwidth*(fit_below_point1 + first_zero_bin);
					}
					else
					{
						//too shallow a decay - use default maximum
						threshold = full_hist.size();
					}
					//find first empty bin after that
				}
				else threshold = full_hist.size();

			}
			else
				threshold = *std::max_element(std::begin(indata), std::end(indata)) + binwidth;
		}
		else
			threshold = *std::max_element(std::begin(indata), std::end(indata)) + binwidth;

		return threshold;
	}

	//Sum up a single month across all years(e.g.all Januaries)
	std::valarray<int> concatenate_months(std::valarray<std::pair<int, int>>& month_ranges, CMaskedArray<float>& data, std::vector<CMaskedArray<float>>& this_month,
		std::vector<int>& year_ids, float missing_value, bool hours = true);

	inline float percentiles(varrayfloat& data, float percent, bool idl=false);
	void winsorize(varrayfloat& data, float percent, bool idl = false);
	float IQR(varrayfloat data, double percentile = 0.25);
	// array de taille 12*nombre d'années où chaque ligne correspond à un mois
	void reshapeMonth(std::vector<std::valarray<std::pair<int, int>>>& month_ranges_years, std::map<int, int>&  month_ranges);
	// array de taille (nombre d'années*12 où chaque ligne correspond à une année
	void reshapeYear(std::vector<std::vector<std::pair<int, int>>>& month_ranges_years, std::map<int, int>&  month_ranges);

	/*  Fit data with a gaussian distribution*/
	inline parameter_vector1 residual_derivativeGaussian(const std::pair<input_vector, double>& data, const parameter_vector1& params)
	{
		parameter_vector1 der(3);

		const double p0 = params(0);
		const double p1 = params(1);
		const double p2 = params(2);

		const double X = data.first(0);

		const double temp = std::exp(-(X - p0)*(X - p0) / (2 * p1*p1));

		der(0) = (p2*(X - p0) / (p1*p1))*temp;

		der(1) = p2*(((X - p0)*(X - p0)) / (p1*p1*p1))*temp;
		der(2) = temp;

		return der;
	}
	inline double modelGaussian(const input_vector&  input, const parameter_vector1& params)
	{
		/*	p[0] = mean
		p[1] = sigma
		p[2] = normalisation
		*/
		const double p0 = params(0);
		const double p1 = params(1);
		const double p2 = params(2);

		const double X = input(0);
		const double temp = p2*(std::exp(-(X - p0)*(X - p0) / (2 * p1*p1)));

		return temp;
	}
	inline double residualGaussian(const std::pair<input_vector, double>& data, const parameter_vector1& params)
	{
		return modelGaussian(data.first, params) - data.second;
	}
	inline varrayfloat Gaussian(varrayfloat& input, varrayfloat& params)
	{
		/*	p[0] = mean
		p[1] = sigma
		p[2] = normalisation
		*/
		const double p0 = params[0];
		const double p1 = params[1];
		const double p2 = params[2];

		varrayfloat temp(input.size());
		for (size_t i = 0; i < input.size(); i++)
		{
			 double X = input[i];
			temp[i] = p2*(std::exp(-(X - p0)*(X - p0) / (2 * p1*p1)));
		}

		return temp;
	}
	
	/*Fit a straight line to the data provided
		Inputs:
		  x - x-data
		  y - y-data
		  norm - norm
		Outputs:
		  fit - array of [mu,sigma,norm]
	*/
	inline varrayfloat fit_gaussian(varrayfloat& x, varrayfloat& y, double norm, double mu, double sig)
	{
		parameter_vector1  fit;
		fit(0) = mu;
		fit(1) = sig;
		fit(2) = norm;
		// DATA (X,Y)

		std::vector<std::pair<input_vector, double>> data;
		input_vector input;
		for (size_t i = 0; i < x.size(); ++i)
		{
			input(0) = x[i];
			data.push_back(std::make_pair(input, y[i]));
		}

		dlib::solve_least_squares(dlib::objective_delta_stop_strategy(1e-9).be_verbose(),
			residualGaussian,
			residual_derivativeGaussian,
			data,
			fit);
		varrayfloat fits(3);
		fits[0] = fit(0);
		fits[1] = fit(1);
		fits[2] = fit(2);
		return fits;
	}

	inline float invert_gaussian(float Y, varrayfloat& p)
	{
		return p[0] + (p[1] * std::sqrt(-2 * std::log(Y / p[2])));
	}

	inline double modelGaussianHermite(const input_vector&  input, const parameter_vector2& params)
	{
		//Model is a Gauss-Hermite function
		
		const double A = params(0);
		const double x0 = params(1);
		const double s = params(2);
		const double h3 = params(3);
		const double h4 = params(4);
		/* The Gauss-Hermite function is a superposition of functions of the form
		 F = (x-xc)/s                                            
			E =  A.Exp[-1/2.F^2] * {1 + h3[c1.F+c3.F^3] + h4[c5+c2.F^2+c4.F^4]} 
			From http://www.astro.rug.nl/software/kapteyn-beta/plot_directive/EXAMPLES/kmpfit_gausshermite.py
		*/
		const double c0 = std::sqrt(6.0) / 4.0;
		const double c1 = -std::sqrt(3.0);
		const double c2 = -std::sqrt(6.0);
		const double c3 = 2.0*(std::sqrt(3.0)) / 3.0;
		const double c4 = std::sqrt(6.0) / 3.0;

		const double X = input(0);

		const double F = (X - x0) / s;

		const double E = A*std::exp(-0.5*F*F)*(1.0 + h3*F*(c3*F*F + c1) + h4*(c0 + F*F*(c2 + c4*F*F)));

		return E;
	}

	inline double residualGaussienHermite(const std::pair<input_vector, double>& data, const parameter_vector2& params)
	{
		return modelGaussianHermite(data.first, params) - data.second;
	}
	inline parameter_vector2 residual_derivativeGaussienHermite(const std::pair<input_vector, double>& data, const parameter_vector2& params)
	{
		/* Derivée de la fonction de gauss hermite par rapport aux paramètres*/
		parameter_vector2 der(5);

		const double A = params(0);
		const double x0 = params(1);
		const double s = params(2);
		const double h3 = params(3);
		const double h4 = params(4);

		const double c0 = std::sqrt(6.0) / 4.0;
		const double c1 = -std::sqrt(3.0);
		const double c2 = -std::sqrt(6.0);
		const double c3 = 2.0*(std::sqrt(3.0)) / 3.0;
		const double c4 = std::sqrt(6.0) / 3.0;

		const double X = data.first(0);

		const double F = (X - x0) / s;

		const double E = A*std::exp(-0.5*F*F)*(1.0 + h3*F*(c3*F*F + c1) + h4*(c0 + F*F*(c2 + c4*F*F)));

		

		/*const double temp = std::exp(-(X - p0)*(X - p0) / (2 * p1*p1));

		der(0) = (p2*(X - p0) / (p1*p1))*temp;

		der(1) = p2*(((X - p0)*(X - p0)) / (p1*p1*p1))*temp;
		der(2) = temp;*/

		return der;
	}

	inline varrayfloat fit_gaussian_Hermite(const varrayfloat& initial_values, varrayfloat& x, varrayfloat& y, varrayfloat& err)
	{
		parameter_vector2  fit;
		fit(0) = initial_values[0];
		fit(1) = initial_values[1];
		fit(2) = initial_values[2];
		fit(3) = initial_values[3];
		fit(4) = initial_values[4];
		// DATA (X,Y)

		std::vector<std::pair<input_vector, double>> data;
		input_vector input;
		for (size_t i = 0; i < x.size(); ++i)
		{
			input(0) = x[i];
			data.push_back(std::make_pair(input, y[i]));
		}

		/*dlib::solve_least_squares_lm(dlib::objective_delta_stop_strategy(1e-7).be_verbose(),
			residualGaussienHermite,
			residual_derivativeGaussian,
			data,
			fit);*/
		varrayfloat fits(3);
		fits[0] = fit(0);
		fits[1] = fit(1);
		fits[2] = fit(2);
		return fits;
	}
	/*  Fit data with  distribution of Gauss_Hermite*/
	
	template<typename T>
	T Somme(std::vector<T> vec)
	{
		T sum = 0;
		for (size_t i = 0; i < vec.size(); i++)
		{
			sum += vec[i];
		}
		return sum;
	}

	float mean_absolute_deviation(const varrayfloat& data, bool median = false);

}