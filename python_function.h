#pragma once

#include "Utilities.h"

#include <dlib/matrix.h>
#include<algorithm>
#include <vector>
#include <valarray>



template<typename T> class CMaskedArray
{
public:
	/****** Constructeurs*****************************************************/
	CMaskedArray()
	{
		m_fill_value = 0;
		m_data = 0;
		m_mask = false;
	}
	CMaskedArray(std::valarray<T> data)
	{
		m_fill_value = static_cast<T>(1e20);
		m_data=data;
		m_mask.resize(data.size(), false);
	}
	CMaskedArray(size_t n, T *data, T missing_value)
	{
		m_fill_value = missing_value;
		valarray<T> dum(data, n);
		m_data= dum;
		m_mask.resize(n, false);
		masked(missing_value);
	}
	CMaskedArray(T missing_value,std::valarray<T> data)
	{
		m_fill_value = missing_value;
		m_data = data;
		m_mask.resize(data.size(), false);
		masked(missing_value);
	}
	CMaskedArray(T missing_value, size_t n,bool missing)
	{
		m_fill_value = missing_value;
		m_data.resize(n);
		m_mask.resize(n);
		masked(missing_value);
	}
	CMaskedArray(T data, size_t size)
	{
		m_fill_value = T(1e20);
		m_data.resize(size);
		m_data = data;
		m_mask.resize(size);
		m_mask = false;
	}
	CMaskedArray(size_t  size)
	{
		m_fill_value = T(1e20);
		m_data.resize(size);
		m_mask.resize(size);
		m_mask= false;
	}

	/*********************  END constructeurs ***************/
	CMaskedArray(CMaskedArray<T> const& mask_array_copy)
	{
		m_fill_value = mask_array_copy.m_fill_value;
		m_data = mask_array_copy.m_data;
		m_mask = mask_array_copy.m_mask;
		m_masked_indices = mask_array_copy.m_masked_indices;
	}
	
	//T fill_value() const{ return m_fill_value; }
	void addData(const varraysize& indices, const std::valarray<T>& M_array)
	{
		m_data[indices] += M_array;
		
	}
	void add(const varraysize& indices, const std::valarray<T>& M_array)
	{
		m_masked_indices.clear();
		
		m_data[indices] += M_array;
		masked(m_fill_value);
	}
	void fill(T value, bool mask = false)
	{
		if (mask)
		{
			
			for (size_t i = 0; i < m_data.size(); ++i)
			{
				if (m_data[i] == m_fill_value)
				{
					m_data[i] = value;
				}
			}
		}
		else m_data = value;

		m_fill_value = value;
	}
	void fill(const varraysize& indices,const CMaskedArray<T>& M_array)
	{
		m_masked_indices.clear();
		
		m_data[indices] = M_array;
		masked(m_fill_value);
	}

	void fill(const varraysize& indices,const std::valarray<T>& M_array)
	{
		m_masked_indices.clear();
		int j = 0;
		m_data[indices] = M_array;
		masked(m_fill_value);
	}
	void fillData(const varraysize& indices, const std::valarray<T>& M_array)
	{
		
		int j = 0;
		m_data[indices] = M_array;
		
	}
	void fill(const std::valarray<bool>& indices, const CMaskedArray<T>& M_array)
	{
		m_masked_indices.clear();
		int j = 0;
		for (size_t i = 0; i < indices.size(); ++i)
		{
			if (indices[i] == true)
			{
				m_data[i] = M_array.m_data[j++];
			}
		}
		masked(m_fill_value);
	}
	void fill(const std::valarray<bool>& indices, const std::valarray<T>& M_array)
	{
		m_masked_indices.clear();
		int j = 0;
		for (size_t i = 0; i < indices.size();i++)
		{
			if (indices[i] == true)
			{
				m_data[i] = M_array[j++];
			}	
		}
		masked(m_fill_value);
	}
	T ma_sum()
	{
		T sum = 0;
		for (int i = 0; i < m_data.size(); ++i)
		{
			if (!m_mask[i]) sum += m_data[i];
		}
		return sum;
	}
	T ma_mean()
	{
		T sum = 0;
		for (int i = 0; i < m_data.size(); ++i)
		{
			if (!m_mask[i]) sum += m_data[i];
		}
		return sum/m_data.size();
	}
	//std::valarray<T> data(){ return m_data; }
	//std::valarray<bool>	mask(){ return m_mask; }
	//void Setmask(std::valarray<bool> bool_array){ m_mask = bool_array; }
	bool nomask()
	{ 
		if (m_mask.max() == true)
			return false;
		else return true;
	}
	
	std::valarray<T> compressed(){ return m_data[!m_mask]; }

	void masked(size_t index)
	{
		m_mask[index] = true;
		m_masked_indices.push_back(index);
	}
	void masked(T missing_value,bool ma=false)
	{
		
		if (!ma)
		{
			for (size_t i = 0; i < m_data.size(); ++i)
			{
				if (m_data[i] == missing_value)
				{
					m_mask[i] = true;
					m_masked_indices.push_back(i);
				}
				else m_mask[i] = false;
			}
		}
		else /* If ma== true changer la valeur manquante par défaut(m_fill_value) par missing_value*/
		{
			for (size_t i = 0; i < m_data.size(); ++i)
			{
				if (m_data[i] == m_fill_value)
				{
					m_data[i] = missing_value;
					m_mask[i] = true;
					m_masked_indices.push_back(i);
				}
				else m_mask[i] = false;
			}
		}
		m_fill_value = missing_value;
	}
	void masked()
	{
		m_masked_indices.clear();
		for (size_t i = 0; i < m_data.size(); ++i)
		{
			if (m_data[i] == m_fill_value)
			{
				m_mask[i] = true;
				m_masked_indices.push_back(i);
			}
			else m_mask[i] = false;
		}
	}
	size_t size() const{return m_data.size();}
	
	void resize(size_t val,T value= T(0))
	{
		m_data.resize(val,value);
		m_mask.resize(val,false);
		m_masked_indices.clear();
	}
	//void fill_value(T fill){m_fill_value = fill;}

	//void setData(varraysize indices, std::valarray<T> new_data){m_data[indices] = new_data;}
	//void setData(size_t indice, T new_data){m_data[indice] = new_data;}
	void free()
	{
		m_mask.free();
		m_data.free();
		m_masked_indices.clear();
	}
	/***************	Surchage d'opérateurs	******/

	T operator[](int indice){return m_data[indice];}
		
	CMaskedArray<T> CMaskedArray<T>::operator()(size_t lb, size_t ub)
	{
		CMaskedArray<T> ma;
		ma.m_mask = m_mask[std::slice(lb, ub - lb, 1)];
		ma.m_data = m_data[std::slice(lb, ub - lb, 1)];
		ma.m_fill_value = m_fill_value;
		ma.masked(m_fill_value);
		return ma;
	}
	CMaskedArray<T> operator[](varraysize indices)
	{
		CMaskedArray<T> ma;
		ma.m_mask = m_mask[indices];
		ma.m_data = m_data[indices]; 
		ma.m_fill_value = m_fill_value;
		//ma.masked(m_fill_value);
		return ma; 
	}
	T operator[](size_t index){ return m_data[index]; }
	CMaskedArray<T> operator[](std::valarray<bool> indices)
	{
		CMaskedArray<T> ma;
		ma.m_mask = m_mask[indices];
		ma.m_data = m_data[indices];
		ma.m_fill_value = m_fill_value;
		//ma.masked(m_fill_value);
		return ma;
	}
	CMaskedArray<T>& CMaskedArray<T>::operator=(CMaskedArray<T> const& mask_array_copy)
	{
		if (this != &mask_array_copy)
		{
			m_fill_value = mask_array_copy.m_fill_value;
			m_data = mask_array_copy.m_data;
			m_mask = mask_array_copy.m_mask;
			m_masked_indices = mask_array_copy.m_masked_indices;
		}
		return *this;
	}
	CMaskedArray<T> CMaskedArray<T>::operator-(CMaskedArray<T> const& m1)
	{
		CMaskedArray<T> ma;
		ma.m_mask = m_mask*m1.m_mask;
		ma.m_data = m_data - m1.m_data;
		return ma;
	}	
	CMaskedArray<T> CMaskedArray<T>::operator+=(CMaskedArray<T> const& m1)
	{
		CMaskedArray<T> ma;
		ma.m_mask = m_mask*m1.m_mask;
		ma.m_data = m_data + m1.m_data;
		return ma;
	}
	CMaskedArray<T>& CMaskedArray<T>::operator/=(float val)
	{
		m_data /= val;
		m_data[m_mask] = m_fill_value;
		return *this;
	}
	/***************	Fin surchage d'opérateurs	******/
public:

	T m_fill_value;
	std::valarray<T> m_data;
	std::valarray<bool> m_mask; 
	std::vector<size_t> m_masked_indices;
	
}; 




namespace PYTHON_FUNCTION
{
	/*Return evenly spaced numbers over a specified interval.
	//
	//Returns `num` evenly spaced samples, calculated over the
	//interval[`start`, `stop`].*/
	template<typename T>
	inline void linspace(std::vector<T> *t, int start, int stop, int number)
	{

		int delta = int((stop - start) / (number - 1));
		for (int i = 0; i < number; i++)
			(*t).push_back(start + i*delta);
	}
	
	inline void arange(std::vector<int>& t, int stop, int start = 0)
	{

		for (int i = start; i < stop; i++)
			t.push_back(i);
	}
	template<typename T>
	inline std::valarray<T> arange(T stop, T start = 0, T step = 1)
	{
		int taille = int(ceil((stop - start) / step));
		std::valarray<T> val(taille);

		for (int j = 0; j < taille; j++, start = start + step)
			val[j] = start;
		return val;
		/*int taille = int((stop - start) / step);
		std::valarray<T> val(taille);
		int j = 0;
		for (int i = start, j = 0; i < stop; i += step, j++)
			val[j] = static_cast<T>(i);
		return val;*/
	}
	template<typename T>
	inline varraysize Arange(T stop, T start=0, T step=1)
	{
		int taille = int((stop - start) / step);
		varraysize val(taille);

		for (int j = 0; j < taille; j++, start = start + step)
			val[j] = int(ceil(start));
		return val;
	}
	/*Test whether each element of a 1 - D array is also present in a second array
	Returns a boolean array the same length as array 1 that is true where
	an element of array 1 is in array 2 and false otherwise
	@Docstring numpy.in1D python*/
	template<typename T,typename S>
	inline std::valarray<bool> in1D(std::valarray<T> v1, std::valarray<S> v2)
	{
		std::valarray<bool> vec(v1.size());
	
		for (int i = 0; i < v1.size(); i++)
		{
			if (std::binary_search(std::begin(v2), std::end(v2), static_cast<S>(v1[i])))
				vec[i]=true;
			else vec[i]=false;
		}
		return vec;
	}


	/*template<typename T>
	inline vind_matrice np_where(matrice& v1, T value, const std::string condition)
	{

		 vind_matrice indices;
		 ind_matrice ligne(v1.size1(), 1);
		 ind_matrice col(v1.size2(), 1);
		 size_t pos = 0;

		if (condition == "=")
			for (size_t i = 0; i < v1.size1(); i++)
				for (size_t j = 0; j < v1.size2(); j++)
			{
				if (v1(i,j) == value)
				{
					ligne(pos, 1) = i;
					col(pos, 1) = j;
					pos++;
				}
			}

		if (condition == "!")
			for (size_t i = 0; i < v1.size(); i++)
				for (size_t j = 0; j < v1.size2(); j++)
			{
				if (v1(i, j) != value)
				{
					ligne(pos, 1) = i;
					col(pos, 1) = j;
					pos++;
				}
			}

		if (condition == "<=")
			for (size_t i = 0; i < v1.size(); i++)
				for (size_t j = 0; j < v1.size2(); j++)
			{
				if (v1(i, j) <= value) 
				{
					ligne(pos, 1) = i;
					col(pos, 1) = j;
					pos++;
				}
			}

		if (condition == "<")
			for (size_t i = 0; i < v1.size(); i++)
				for (size_t j = 0; j < v1.size2(); j++)
			{
				if (v1(i, j) < value) 
				{
					ligne(pos, 1) = i;
					col(pos, 1) = j;
					pos++;
				}
			}

		if (condition == ">")
			for (size_t i = 0; i < v1.size(); i++)
				for (size_t j = 0; j < v1.size2(); j++)
			{
				if (v1(i, j) > value) 
				{
					ligne(pos, 1) = i;
					col(pos, 1) = j;
					pos++;
				}
			}

		if (condition == ">=")
			for (size_t i = 0; i < v1.size(); i++)
				for (size_t j = 0; j < v1.size2(); j++)
			{
				if (v1(i, j) >= value) 
				{
					ligne(pos, 1) = i;
					col(pos, 1) = j;
					pos++;
				}
			}

		ligne = subm(ligne, range(0, 0), range(0, pos));
		col = subm(col, range(0, 0), range(0, pos));

		indices.push_back(ligne);
		indices.push_back(col);
		return indices;


	}*/
	template<typename T>
	inline std::valarray<T> np_where(std::valarray<T> v1, T value)
	{
		std::vector<size_t> index;
		for (int i = 0; i < v1.size(); i++)
		{
			if (v1[i] != value) index.push_back(i);
		}
		std::valarray<T> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));
		return vec;
	}
	
	template<typename T>
	inline std::valarray<std::pair<std::size_t,std::size_t>> np_ma_where(std::valarray<CMaskedArray<T>>& v1, const std::string condition, T value)
	{

		std::vector<std::pair<std::size_t, std::size_t>> index;
		if (condition == "=")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				for (size_t j = 0; j < v1[i].size(); j++)
				{
					if (v1[i][j] == value && v1[i].m_mask[j] == false) index.push_back(make_pair(i, j));
				}
			}
		}
		if (condition == "!")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				for (size_t j = 0; j < v1[i].size(); j++)
				{
					if (v1[i][j] != value && v1[i].m_mask[j]== false) index.push_back(make_pair(i, j));
				}
			}
		}

		if (condition == "<=")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				for (size_t j = 0; j < v1[i].size(); j++)
				{
					if (v1[i][j] <= value && v1[i].m_mask[j]== false) index.push_back(make_pair(i, j));
				}
			}
		}

		if (condition == "<")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				for (size_t j = 0; j < v1[i].size(); j++)
				{
					if (v1[i][j] < value && v1[i].m_mask[j]== false) index.push_back(make_pair(i, j));
				}
			}
		}

		if (condition == ">")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				for (size_t j = 0; j < v1[i].size(); j++)
				{
					if (v1[i][j] > value  && v1[i].m_mask[j]== false) index.push_back(make_pair(i, j));
				}
			}
		}

		if (condition == ">=")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				for (size_t j = 0; j < v1[i].size(); j++)
				{
					if (v1[i][j] >= value  && v1[i].m_mask[j]== false) index.push_back(make_pair(i, j));
				}
			}
		}

		std::valarray<std::pair<std::size_t,std::size_t>> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));
		return vec;


	}


	template<typename T>
	inline std::valarray<std::size_t> npwhere(const std::valarray<T>& v1, const std::string condition,T value )
	{
		
		std::vector<size_t> index;
		if (condition == "=")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] == value) index.push_back(i);
			}
		}
			
		if (condition == "!")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] != value) index.push_back(i);
			}
		}
		
		if (condition == "<=")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] <= value) index.push_back(i);
			}
		}
			
		if (condition == "<")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] < value) index.push_back(i);
			}
		}
			
		if (condition == ">")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] > value) index.push_back(i);
			}
		}
			
		if (condition == ">=")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] >= value) index.push_back(i);
			}
		}

		std::valarray<std::size_t> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));
		return vec;
		
		
	}
	template<typename T>
	inline std::valarray<std::size_t> npwhere(const std::vector<T>& v1, const std::string condition, T value)
	{

		std::vector<size_t> index;
		if (condition == "=")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] == value) index.push_back(i);
			}
		}

		if (condition == "!")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] != value) index.push_back(i);
			}
		}

		if (condition == "<=")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] <= value) index.push_back(i);
			}
		}

		if (condition == "<")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] < value) index.push_back(i);
			}
		}

		if (condition == ">")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] > value) index.push_back(i);
			}
		}

		if (condition == ">=")
		{
			for (size_t i = 0; i < v1.size(); i++)
			{
			}
		}

		std::valarray<std::size_t> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));
		return vec;


	}
	template<typename T>
	inline std::valarray<std::size_t> npwhere(const std::valarray<T>& v1, const std::string condition, const std::valarray<T>& v2)
	{

		std::vector<size_t> index;
		if (condition == "=")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] == v2[i]) index.push_back(i);
			}

		if (condition == "!")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] != v2[i]) index.push_back(i);
			}

		if (condition == "<=")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] <= v2[i]) index.push_back(i);
			}

		if (condition == "<")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] < v2[i]) index.push_back(i);
			}

		if (condition == ">")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] > v2[i]) index.push_back(i);
			}

		if (condition == ">=")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] >= v2[i]) index.push_back(i);
			}
		std::valarray<std::size_t> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));
		return vec;


	}
	template<typename T>
	inline std::valarray<std::size_t> npwhere(const std::vector<T>& v1, const std::string condition, const std::valarray<T>& v2)
	{

		std::vector<size_t> index;
	
		if (condition == "=")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] == v2[i]) index.push_back(i);
			}

		if (condition == "!")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] != v2[i]) index.push_back(i);
			}

		if (condition == "<=")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] <= v2[i]) index.push_back(i);
			}

		if (condition == "<")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] < v2[i]) index.push_back(i);
			}

		if (condition == ">")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] > v2[i]) index.push_back(i);
			}

		if (condition == ">=")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] >= v2[i]) index.push_back(i);
			}
		std::valarray<std::size_t> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));
		return vec;


	}
	template<typename T>
	inline std::valarray<std::size_t> npwhere(const std::valarray<T>& v1, const std::valarray<T>& v2, const std::string condition)
	{
		
		std::vector<size_t> index;
		if (condition == "=")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] == v2[i]) index.push_back(i);
			}

		if (condition == "!")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] != v2[i]) index.push_back(i);
			}

		if (condition == "<=")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] <= v2[i]) index.push_back(i);
			}

		if (condition == "<")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] < v2[i]) index.push_back(i);
			}

		if (condition == ">")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] > v2[i]) index.push_back(i);
			}

		if (condition == ">=")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1[i] >= v2[i]) index.push_back(i);
			}
		std::valarray<std::size_t> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));
		return vec;


	}
	template<typename T>
	inline std::valarray<std::size_t> npwhereAnd(std::vector<T> v1, T value, T value2, const std::string condition1, const std::string condition2)
	{
		std::vector<size_t> index;
		if (condition1 == "<=")
		{
			if (condition2 == ">")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] <= value && v1[i]> value2) index.push_back(i);
				}
		}
		if (condition1 == "<")
		{
			if (condition2 == ">")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value && v1[i]> value2) index.push_back(i);
				}
		}
		if (condition1 == "<")
		{
			if (condition2 == ">=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value && v1[i] >= value2) index.push_back(i);
				}
		}
		if (condition1 == "<=")
		{
			if (condition2 == ">=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value && v1[i] >= value2) index.push_back(i);
				}
		}

		std::valarray<std::size_t> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));
		return vec;


	}
	template<typename T>
	inline std::valarray<std::size_t> npwhereAnd(std::valarray<T> v1, T value, T value2, const std::string condition1, const std::string condition2)
	{

		std::vector<size_t> index;
		if (condition1 == "<=")
		{
			if (condition2 == ">")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] <= value && v1[i]> value2) index.push_back(i);
				}
		}
		if (condition1 == "<")
		{
			if (condition2 == ">")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value && v1[i]> value2) index.push_back(i);
				}
		}
		if (condition1 == "<")
		{
			if (condition2 == ">=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value && v1[i] >= value2) index.push_back(i);
				}
		}
		if (condition1 == "<=")
		{
			if (condition2 == ">=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value && v1[i] >= value2) index.push_back(i);
				}
		}
		std::valarray<std::size_t> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));
		return vec;


	}
	template<typename T,typename S>
	inline std::valarray<std::size_t> npwhereAnd(std::valarray<T> v1,const std::string condition1, T value, std::valarray<S> v2, const std::string condition2,S value2)
	{
		std::vector<size_t> index;
		if (condition1 == "=")
		{
			if (condition2 == "=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] == value && v2[i]== value2) index.push_back(i);
				}
		}
		if (condition1 == "=")
		{
			if (condition2 == "!")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] == value && v2[i] != value2) index.push_back(i);
				}
		}
		if (condition1 == "<=")
		{
			if (condition2 == ">")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] <= value && v2[i]> value2) index.push_back(i);
				}
		}
		if (condition1 == "<")
		{
			if (condition2 == ">")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value && v2[i]> value2) index.push_back(i);
				}
		}
		if (condition1 == "<")
		{
			if (condition2 == ">=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value && v2[i] >= value2) index.push_back(i);
				}
		}
		if (condition1 == "<")
		{
			if (condition2 == "!")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value && v2[i] != value2) index.push_back(i);
				}
		}
		if (condition1 == "<=")
		{
			if (condition2 == ">=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value && v2[i] >= value2) index.push_back(i);
				}
		}
		std::valarray<std::size_t> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));
		return vec;


	}
	template<typename T>
	inline std::valarray<std::size_t> npwhereOr(std::valarray<T> v1, T value, T value2, const std::string condition1, const std::string condition2)
	{
		std::vector<size_t> index;
		std::valarray<bool> indices(false,v1.size());
		if (condition1 == "=")
		{
			if (condition2 == "=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] == value || v1[i]== value2) index.push_back(i);
				}

		}
		if (condition1 == "<=")
		{
			if (condition2 == ">")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] <= value || v1[i]> value2) index.push_back(i);
				}
		}
		if (condition1 == "<")
		{
			if (condition2 == ">")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value || v1[i]> value2) index.push_back(i);
				}
		}
		if (condition1 == "<")
		{
			if (condition2 == ">=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value || v1[i] >= value2) index.push_back(i);
				}
		}

		std::valarray<std::size_t> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));

		return vec;


	}
	template<typename T, typename S>
	inline std::valarray<std::size_t> npwhereOr(std::valarray<T> v1, const std::string condition1, T value, std::valarray<S> v2, const std::string condition2, S value2)
	{

		std::vector<size_t> index;
		if (condition1 == "=")
		{
			if (condition2 == "=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] == value || v2[i] == value2) index.push_back(i);
				}
		}
		if (condition1 == "=")
		{
			if (condition2 == "!")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] == value || v2[i] != value2) index.push_back(i);
				}
		}
		if (condition1 == "<=")
		{
			if (condition2 == ">")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] <= value || v2[i]> value2) index.push_back(i);
				}
		}
		if (condition1 == "<")
		{
			if (condition2 == ">")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value || v2[i]> value2) index.push_back(i);
				}
		}
		if (condition1 == "<")
		{
			if (condition2 == ">=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value || v2[i] >= value2) index.push_back(i);
				}
		}
		if (condition1 == "<=")
		{
			if (condition2 == ">=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1[i] < value || v2[i] >= value2) index.push_back(i);
				}
		}
		std::valarray<std::size_t> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));
		return vec;


	}
	template<typename T>
	inline std::valarray<std::size_t> npwhere(CMaskedArray<T>& v1, const std::string condition, T value)
	{
		std::vector<size_t> index;
		if (condition == "=")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1.m_mask[i] == false && v1.m_data[i] == value)
					index.push_back(i);
			}

		if (condition == "!")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1.m_mask[i] == false && v1.m_data[i] != value) index.push_back(i);
			}

		if (condition == "<=")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1.m_mask[i] == false && v1.m_data[i] <= value) index.push_back(i);
			}

		if (condition == "<")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1.m_mask[i] == false && v1.m_data[i] < value) index.push_back(i);
			}

		if (condition == ">")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1.m_mask[i] == false && v1.m_data[i] > value) index.push_back(i);
			}

		if (condition == ">=")
			for (size_t i = 0; i < v1.size(); i++)
			{
				if (v1.m_mask[i] == false && v1.m_data[i] >= value) index.push_back(i);
			}
		std::valarray<std::size_t> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));
		return vec;

		
	}
	template<typename T>
	inline std::vector<int> np_where_vec(std::valarray<T> v1, T value)
	{
		std::vector<int> vec;
		for (int i = 0; i < v1.size(); i++)
		{
			if (v1[i] != value) vec.push_back(i);
		}
		return vec;
	}
	template<typename T,typename S>
	CMaskedArray<S> ma_masked_where(std::valarray<T> data,  T value, std::valarray<S> data2,S missing_value)
	{
		CMaskedArray<S> dummy = CMaskedArray<S>(data2);
		dummy.masked(missing_value);
		for (size_t i = 0; i < data.size(); i++)
		{
			if (data[i] == value)
			{
				dummy.masked(i);
			}
		}
		return dummy;
	}
	template<typename T>
	CMaskedArray<T> ma_masked_where(std::valarray<T> data, const std::string condition, T value, std::valarray<T> data2)
	{
		CMaskedArray<T> dummy = CMaskedArray<S>(data2);
		if (condition == "=")
		{
			for (size_t i = 0; i < data.size(); i++)
			{
				if (data[i] == value) 
				{
					dummy.masked(i);
					
				}
			}
		}
		dummy.m_fill_value = value;
		return dummy;
	}
	template<typename T>
	varraysize ma_masked_where(CMaskedArray<T> data, const std::string condition,  T value)
	{
		std::vector<size_t> index;
		if (condition == "=")
		{
			for (size_t i = 0; i < data.size(); i++)
			{
				if (data.m_mask[i] == false && data.m_data[i] == value)
					index.push_back(i);
			}
		}

		if (condition == "!")
		{
			for (size_t i = 0; i < data.size(); i++)
			{
				if (data.m_mask[i] == false && data.m_data[i] != value) index.push_back(i);
			}
		}

		if (condition == "<=")
		{
			for (size_t i = 0; i < data.size(); i++)
			{
				if (data.m_mask[i] == false && data.m_data[i] <= value) index.push_back(i);
			}
		}

		if (condition == "<")
		{
			for (size_t i = 0; i < data.size(); i++)
			{
				if (data.m_mask[i] == false && data.m_data[i] < value) index.push_back(i);
			}
		}

		if (condition == ">")
		{
			for (size_t i = 0; i < data.size(); i++)
			{
				if (data.m_mask[i] == false && data.m_data[i] > value) index.push_back(i);
			}
		}

		if (condition == ">=")
		{
			for (size_t i = 0; i < data.size(); i++)
			{
				if (data.m_mask[i] == false && data.m_data[i] >= value) index.push_back(i);
			}
		}
		std::valarray<std::size_t> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));
		return vec;

	}
	template<typename T, typename S>
	inline varraysize ma_masked_whereAnd(CMaskedArray<T> v1, const std::string condition1, T value, CMaskedArray<T> v2, const std::string condition2, S value2)
	{

		std::vector<size_t> index;
		if (condition1 == "=")
		{
			if (condition2 == "=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1.m_mask[i] == false && v2.m_mask[i] == false)
						if (v1[i] == value && v2[i] == value2) index.push_back(i);
				}
		}
		if (condition1 == "=")
		{
			if (condition2 == "!")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1.m_mask[i]==false && v2.m_mask[i] == false)
						if (v1[i] == value && v2[i] != value2)
							index.push_back(i);
				}
		}
		if (condition1 == "<=")
		{
			if (condition2 == ">")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1.m_mask[i] == false && v2.m_mask[i] == false)
						if (v1[i] <= value && v2[i]> value2) index.push_back(i);
				}
		}
		if (condition1 == "<")
		{
			if (condition2 == ">")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1.m_mask[i] == false && v2.m_mask[i] == false)
						if (v1[i] < value && v2[i]> value2) index.push_back(i);
				}
		}
		if (condition1 == "<")
		{
			if (condition2 == ">=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1.m_mask[i] == false && v2.m_mask[i] == false)
						if (v1[i] < value && v2[i] >= value2) index.push_back(i);
				}
		}
		if (condition1 == "<=")
		{
			if (condition2 == ">=")
				for (size_t i = 0; i < v1.size(); i++)
				{
					if (v1.m_mask[i] == false && v2.m_mask[i] == false)
						if (v1[i] < value && v2[i] >= value2) index.push_back(i);
				}
		}
		std::valarray<std::size_t> vec(index.size());
		std::copy(index.begin(), index.end(), std::begin(vec));
		return vec;


	}
	template<typename T>
	inline std::vector<int> np_where(std::slice_array<T> v1, std::slice_array<T> v2)
	{
		std::vector<int> vec;
		std::valarray<T> vl1 = v1, vl2=v2;
		for (int i = 0; i < v1.size(); i++)
		{
			if (vl1[i]== vl2[i]) vec.push_back(i);
		}
		return vec;
	}
	/*
	    Return a valarray where value is deleted
		equivalent to compressed,masqued_equals here
	*/
	template<typename T>
	inline std::valarray <T> masked_values(std::valarray<T>& v1, T value)
	{
		std::valarray<bool> masque(true,v1.size());
		for (int i = 0; i < v1.size(); i++)
		{
			if (v1[i] == value) masque[i]=false ;
		}
		return v1[masque];
	}
	template<typename T>
	inline CMaskedArray<T> ma_masked_values(std::valarray<T>& v1, T value)
	{
		CMaskedArray<T> r_array(value, v1);
		r_array.masked(value);

		return r_array;
	}
	template<typename T>
	inline CMaskedArray<T> ma_masked_values(const CMaskedArray<T>& v1, T value)
	{
		CMaskedArray<T> r_array=v1;
		r_array.fill(value,true);
		r_array.masked(value);
		return r_array;
	}
	template<typename T>
	inline varrayfloat histogram(std::valarray<T>& data, varrayfloat& bin,bool density=false)
	{
		varrayfloat hist(bin.size() - 1);

		bool lastbin = true;
		for (int i = 0; i < data.size(); i++)
		{
			lastbin = true;
			for (int j = 0; j < bin.size()-1; j++)
			{
				if (data[i] >= bin[j] && data[i] < bin[j + 1])
					hist[j]++;
				if (data[i] == bin[bin.size() - 1] && lastbin)
				{
					hist[bin.size() - 2]++;
					lastbin = false;
				}
			}
		}
		if (density)
		{
			for (int j = 0; j < bin.size() - 1; j++)
			{
				 hist[j] /= ((bin[j+1]-bin[j])*data.size());
			}
		}
		return hist;

	}
	/*template<typename T>
	inline varrayfloat histogram(std::valarray<CMaskedArray<T>>& data, varrayfloat& bin, bool density = false)
	{
		varrayfloat hist(bin.size() - 1);

		bool lastbin = true;
		for (int i = 0; i < data.size(); i++)
			for (int j = 0; j < bin.size(); j++)
			{
				if (data[i] >= bin[j] && data[i] < bin[j + 1])
					hist[j]++;
				if (data[i] == bin[bin.size() - 1] && lastbin)
				{
					hist[bin.size() - 2]++;
					lastbin = false;
				}
			}
		if (density)
		{
			if (data.size() <= bin.size()) hist /= data.size();
			else hist /= bin.size();
		}
		return hist;

	}*/
	template<typename T>
	inline varrayfloat histogram(std::vector<T> data, varrayfloat bin)
	{
		varrayfloat hist(bin.size() - 1);

		bool lastbin = true;
		for (int i = 0; i < data.size(); i++)
			for (int j = 0; j < bin.size(); j++)
			{
				if (T(data[i]) >= bin[j] && (data[i]) < bin[j + 1])
					hist[j]++;
				if (T(data[i]) == bin[bin.size() - 1] && lastbin)
				{
					hist[bin.size() - 2]++;
					lastbin = false;
				}
			}
		return hist;
	}
	inline void Concatenate(varrayfloat &val1, const varrayfloat& val2)
	{
		varrayfloat val(val1.size() + val2.size());
		for (size_t i = 0; i < val1.size(); i++)
		{
			val[i] = val1[i];
		}
		int j = 0;
		for (size_t i = val1.size(); i < val.size(); i++)
		{
			val[i] = val2[j++];
		}
		val1.resize(val.size());
		val1 = val;
	}
	inline void Concatenate(CMaskedArray<float>& val1, const CMaskedArray<float>& val2)
	{
		CMaskedArray<float>  val(val1.size() + val2.size());
		for (size_t i = 0; i < val1.size(); i++)
		{
			val.m_data[i] = val1.m_data[i];
			val.m_mask[i] = val1.m_mask[i];
		}
		int j = 0;
		for (size_t i = val1.size(); i < val.size(); i++)
		{
			val.m_data[i] = val2.m_data[j];
			val.m_mask[i] = val2.m_mask[j];
			j++;
		}
		val.masked(val1.m_fill_value);
		val1.resize(val.size());
		val1 = val;
	}
	inline void Concatenate(std::vector<CMaskedArray<float>>& val1, const std::vector<CMaskedArray<float>>& val2)
	{
		std::vector<CMaskedArray<float>>  val;
		for (size_t i = 0; i < val1.size(); i++)
		{
			val.push_back(val1[i]);
			
		}
		
		for (size_t i = 0; i < val2.size(); i++)
		{
			val.push_back(val2[i]);
		}
		val1.clear();
		val1 = val;
	}
	template<typename T>
	void np_roll(std::valarray<T> &val, int valeur)
	{
		std::valarray<T> val1(val.size());
		int  value;
		if (valeur > 0) value = valeur;
		else value = val.size() + valeur;
		
		for (size_t i = 0; i < value; i++)
			val1[i] = val[val.size() - value + i];

		for (size_t i = value; i < val.size(); i++)
			val1[i] = val[i - value];
		val = val1;
	}


	// recuperer l'indice du minimum de l'array (seule la première occurence est retournée)
	inline int np_argmin(const varrayfloat& val)
	{
		float min = val.min();
		int j = 0;
		for (int i = 0; i < val.size();i++)
		{
			if (val[i] == min)
			{
				j = i;
				break;
			}
		}
		return j;
	}
	template<typename T>
	inline T ma_min(CMaskedArray<T> indata)
	{
		return indata.compressed().min();
	}
	template<typename T>
	inline T ma_max(CMaskedArray<T> indata)
	{
		return indata.compressed().max();
	}
	// recuperer l'indice du maximun de l'array (seule la première occurence est retournée)
	inline int np_argmax(const varrayfloat& val)
	{
		float max= val.max();
		int j = 0;
		for (int i = 0; i < val.size(); i++)
		{
			if (val[i] == max)
			{
				j = i;
				break;
			}
		}
		return j;
	}
	
	template<typename T>
	inline std::vector<T> np_abs(std::vector<T> data)
	{
		std::vector<T> abs_value;
		for (T value : data)
			abs_value.push_back(std::abs(value));

		return abs_value;
	}
	//Calculer l'ecart-type
	inline double stdev(varrayfloat&  indata,double mean)
	{

		double temp = 0.;
		for (int i = 0; i < indata.size(); i++)
		{
			temp += (indata[i] - mean) * (indata[i] - mean);
		}
		return std::sqrt(temp / indata.size());
	}

	inline float Round(float value, int digits)
	{
		return floor(value * pow(10, digits) + value) / pow(10, digits);
	}
	template<typename T>
	inline std::vector<T> np_ceil(std::vector<T> data)
	{
		std::vector<T> ceil_value;
		for (T value : data)
			ceil_value.push_back(T(std::ceilf(value)));

		return ceil_value;
	}
	template<typename T>
	inline std::valarray<T> np_ceil(std::valarray<T> data)
	{
		std::valarray<T> ceil_value(data.size());
		for ( size_t i = 0; i < data.size();++i)
			ceil_value[i]= T(std::ceilf(data[i]));

		return ceil_value;
	}
	template<typename T>
	inline std::valarray<T> Shape(std::vector<std::valarray<T>> data)
	{
		size_t line = data.size();
		size_t col = data[0].size();
		std::valarray<T> m_data(line*col);
		
		for (size_t i = 0; i < m_data.size();)
		{
			for (size_t j = 0; j < line; j++)
			{
				for (size_t k = 0; k < col; k++)
					m_data[i++] = data[j][k];
			}
		}
		return m_data;
	}
	
	template<typename T>
	inline CMaskedArray<T> Shape(std::valarray<std::valarray<T>> data)
	{
		size_t line = data.size();
		size_t col = data[0].size();
		CMaskedArray<T> masked_data(line*col);
		for (size_t i = 0; i < masked_data.size();)
		{
			for (size_t j = 0; j < line; j++)
			{
				for (size_t k = 0; k < col; k++)
					masked_data.m_data[i++] = data[j][k];
			}
		}
		return masked_data;
	}
	template<typename T>
	inline std::valarray<T> Shape(std::valarray<CMaskedArray<T>> data)
	{
		size_t line = data.size();
		size_t col = data[0].size();
		std::valarray<T> m_data(line*col);

		for (size_t i = 0; i < m_data.size();)
		{
			for (size_t j = 0; j < line; j++)
			{
				for (size_t k = 0; k < col; k++)
					m_data[i++] = data[j].m_data[k];
			}
		}
		return m_data;
	}

	template<typename T>
	inline std::vector<std::valarray<T>> C_reshape(std::vector<T> data,size_t col)
	{
		std::vector<valarray<T>> m_data;
		int iteration = 1;
		std::valarray<T> month(col);
		for (int i : data)
		{
			if (iteration <= col)
			{
				month[iteration - 1] = i;
				iteration++;
			}
			else
			{
				m_data.push_back(month);
				month.resize(col);
				month[0] = i;
				iteration = 2;
			}
		}
		m_data.push_back(month);

		return m_data;
	}
	
	template<typename T>
	inline std::vector<std::valarray<T>> C_reshape(std::valarray<T> data, size_t col)
	{
		std::vector<std::valarray<T>> m_data;
		int iteration = 1;
		std::valarray<T> month(col);
		for (int i : data)
		{
			if (iteration <= col)
			{
				month[iteration - 1] = i;
				iteration++;
			}
			else
			{
				m_data.push_back(month);
				month.resize(col);
				month[0] = i;
				iteration = 2;
			}
		}
		m_data.push_back(month);

		return m_data;
	}
	
	inline std::vector<CMaskedArray<float>> C_reshape( CMaskedArray<float>& a_data, size_t col)
	{
		std::vector<CMaskedArray<float>> m_data;
		int iteration = 1;
		CMaskedArray<float> dummy(a_data.m_fill_value, col,true);
		for (size_t i = 0; i < a_data.size();++i)
		{
			if (iteration <= col)
			{
				dummy.m_data[iteration - 1] = a_data.m_data[i];
				dummy.m_mask[iteration - 1] = a_data.m_mask[i];
				iteration++;
			}
			else
			{
				dummy.masked(a_data.m_fill_value);
				m_data.push_back(dummy);
				dummy.resize(col);
				dummy.m_data[0] = a_data.m_data[i];
				dummy.m_mask[0] = a_data.m_mask[i];
				iteration = 2;
			}
		}
		dummy.masked(a_data.m_fill_value);
		m_data.push_back(dummy);
		
		return m_data;
	}

	template<typename T>
	inline std::vector<std::valarray<T>> L_reshape(const std::vector<T>& data, const size_t line)
	{

		std::vector<std::valarray<T>> m_data;
		int col = int(data.size() / line);
		std::valarray<T> month(col);
		int index = 0;
		int compteur = 0;
		int iteration = 1;

		for (int i = 0; i < data.size(); i += line)
		{
			if (iteration <= col && i < data.size() && compteur<line)
			{
				month[index++] = data[i];
				iteration++;
			}
			else
			{
				if (compteur == line) break;
				m_data.push_back(month);
				compteur++;
				month.resize(col);
				index = 0;
				i = m_data.size();
				month[index++] = data[i];
				iteration = 2;
			}
			if (i + line >= data.size() && compteur<line)
			{
				i = m_data.size();
			}

		}

		return m_data;
	}
	template<typename T>
	inline std::vector<std::valarray<T>> L_reshape(const std::valarray<T>& data, const size_t line)
	{
		std::vector<std::valarray<T>> m_data;
		int col = int(data.size() / line);
		std::valarray<T> month(col);
		int index = 0;
		int compteur = 0;
		int iteration = 1;

		for (int i = 0; i < data.size(); i += line)
		{
			if (iteration <= col && i < data.size() && compteur<line)
			{
				month[index++] = data[i];
				iteration++;
			}
			else
			{
				if (compteur == line) break;
				m_data.push_back(month);
				compteur++;
				month.resize(col);
				index = 0;
				i = m_data.size();
				month[index++] = data[i];
				iteration = 2;
			}
			if (i + line >= data.size() && compteur<line)
			{
				i = m_data.size();
			}

		}

		return m_data;
	}

	inline std::vector<CMaskedArray<float>> L_reshape(CMaskedArray<float>& a_data, const size_t line)
	{
		std::vector<CMaskedArray<float>> m_data;
		int col = int(a_data.size() / line);
		CMaskedArray<float> month(0, col);
		month.masked(a_data.m_fill_value);
		size_t index = 0;
		int compteur = 0;
		int iteration = 1;

		for (size_t i = 0; i < a_data.size(); i += line)
		{
			if (iteration <= col && i < a_data.size() && compteur<line)
			{
				month.m_data[index] = a_data[i];
				month.m_mask[index++] = a_data.m_mask[i];
				iteration++;
			}
			else
			{
				if (compteur == line) break;
				m_data.push_back(month);
				compteur++;
				month.resize(col);
				index = 0;
				i = m_data.size();
				month.m_data[index] = a_data[i];
				month.m_mask[index++] = a_data.m_mask[i];
				iteration = 2;
			}
			if (i + line >= a_data.size() && compteur<line)
			{
				i = m_data.size();
			}

		}

		return m_data;
	}
	template<typename T>
	inline std::vector<std::valarray<std::pair<T, T>>> L_reshape3(std::vector<std::pair<T, T>>& data, const size_t dimline)
	{
		std::vector<std::valarray<std::pair<T, T>>> m_data;
		int index = 0;
		int compteur = 0;
		int iteration = 1;
		size_t dimCol = int(data.size() / (dimline));
		std::valarray<pair<T, T>> month(dimCol);
		std::vector<pair<int, int>>::iterator month_it = data.begin();
		for (int i = 0; i < data.size(); i ++)
		{
			if (iteration <= dimCol && i < data.size() && compteur<dimline)
			{
				month[index++] = make_pair(month_it->first, month_it->second);
				if (i + dimline < data.size()&& iteration< dimCol) std::advance(month_it, dimline);
				iteration++;
			}
			else
			{
				if (compteur == dimline) break;
				m_data.push_back(month);
				compteur++;
				month.resize(dimCol);
				index = 0;
				i = m_data.size();
				month_it = data.begin();
				std::advance(month_it, i);
				month[index++] = make_pair(month_it->first, month_it->second);
				if (i + dimline < data.size()) std::advance(month_it, dimline);
				iteration = 2;
			}
			if (i + dimline >= data.size() && compteur<dimline)
			{
				i = m_data.size();
				month_it = data.begin();
				std::advance(month_it, i + 1);
			}
			if (i == data.size() - 1) m_data.push_back(month);

		}
		
		return m_data;
	}

	inline bool IsmaskedValue(varrayfloat data, float missing_value)
	{
		varrayfloat m_data = data[data == missing_value];
		if (m_data.size() == 0) return false;
		else return true;
	}

	inline size_t CompressedSize(varrayfloat data, float missing_value)
	{
		varrayfloat m_data = data[data == missing_value];
		return m_data.size();
	}

	template<typename T>
	std::valarray<T> npDiff(std::valarray<T> data)
	{
		std::valarray<T> diff(data.size() - 1);

		for (size_t i = 0; i < diff.size(); ++i)
		{
			diff[i] = data[i + 1] - data[i];
		}
		return diff;
	}
	/*  Extraire les données correpondant à une colonne d'une matrice*/
	template<typename T>
	CMaskedArray<T> ExtractColumn(std::vector<CMaskedArray<T>> data,size_t col)
	{
		CMaskedArray<T> Exdata(data.size());
		for (size_t i = 0; i < data.size(); i++)
		{
			Exdata.m_data[i] = data[i].m_data[col];
			Exdata.m_mask[i] = data[i].m_mask[col];
		}
		Exdata.masked(data[0].m_fill_value);
		return  Exdata;
	}
	template<typename T>
	std::valarray<T> CompressedMatrice(std::valarray<CMaskedArray<T>> Matrice)
	{
		std::vector<T> data;
		
		for (size_t i = 0; i < Matrice.size(); i++)
		{
			if (Matrice[i].compressed().size()>0)
			{
				for (size_t j = 0; j < Matrice[i].size(); j++)
				{
					if (!Matrice[i].m_mask[j])
						data.push_back(Matrice[i][j]);
				}
			}
			else continue;
		}
		std::valarray<float> vec(data.size());
		std::copy(data.begin(), data.end(), std::begin(vec));
		return  vec;
	}

	template<typename T>
	std::valarray<T> CompressedMatrice(std::vector<CMaskedArray<T>> Matrice)
	{
		std::vector<T> data;

		for (size_t i = 0; i < Matrice.size(); i++)
		{
			if (Matrice[i].compressed().size()>0)
			{
				for (size_t j = 0; j < Matrice[i].size(); j++)
				{
					if (!Matrice[i].m_mask[j])
						data.push_back(Matrice[i][j]);
				}
			}
			else continue;
		}
		std::valarray<float> vec(data.size());
		std::copy(data.begin(), data.end(), std::begin(vec));
		return  vec;
	}
}
