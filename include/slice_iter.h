#pragma once
#include<iostream>
#include<valarray>
#include<algorithm>
#include<numeric>	// for inner_product
using namespace std;

// forward declarations to allow friend declarations:
template<class T> class slice_iter;
template<class T> bool operator==(const slice_iter<T>&, const slice_iter<T>&);
template<class T> bool operator!=(const slice_iter<T>&, const slice_iter<T>&);
template<class T> bool operator< (const slice_iter<T>&, const slice_iter<T>&);

template<class T> 
class slice_iter 
{
	valarray<T>* v;
	slice s;
	size_t curr;	// index of current element

	T& ref(size_t i) const { return (*v)[s.start()+i*s.stride()]; }
//        valarray<double> & const { return (*v)[]};
public:
	slice_iter(valarray<T>* vv, slice ss) :v(vv), s(ss), curr(0) { }

	slice_iter end() const
	{
		slice_iter t = *this;
		t.curr = s.size();	// index of last-plus-one element
		return t;
	}

	slice_iter& operator++() { curr++; return *this; }
	slice_iter operator++(int) { slice_iter t = *this; curr++; return t; }

	T& operator[](size_t i) { return ref(i); }		// C style subscript
	T& operator()(size_t i) { return ref(i); }		// Fortran-style subscript
	T& operator*() { return ref(curr); }			// current element

	friend bool operator==<>(const slice_iter& p, const slice_iter& q);
	friend bool operator!=<>(const slice_iter& p, const slice_iter& q);
	friend bool operator< <>(const slice_iter& p, const slice_iter& q);

};


template<class T>
bool operator==(const slice_iter<T>& p, const slice_iter<T>& q)
{
	return p.curr==q.curr
		&& p.s.stride()==q.s.stride()
		&& p.s.start()==q.s.start();
}

template<class T>
bool operator!=(const slice_iter<T>& p, const slice_iter<T>& q)
{
	return !(p==q);
}

template<class T>
bool operator<(const slice_iter<T>& p, const slice_iter<T>& q)
{
	return p.curr<q.curr
		&& p.s.stride()==q.s.stride()
		&& p.s.start()==q.s.start();
}


//-------------------------------------------------------------



// forward declarations to allow friend declarations:
template<class T> class Cslice_iter;
template<class T> bool operator==(const Cslice_iter<T>&, const Cslice_iter<T>&);
template<class T> bool operator!=(const Cslice_iter<T>&, const Cslice_iter<T>&);
template<class T> bool operator< (const Cslice_iter<T>&, const Cslice_iter<T>&);


template<class T> class Cslice_iter 
{
	valarray<T>* v;
	slice s;
	size_t curr; // index of current element
	const T& ref(size_t i) const { return (*v)[s.start()+i*s.stride()]; }
public:
	Cslice_iter(valarray<T>* vv, slice ss): v(vv), s(ss), curr(0){}
	Cslice_iter end() const
	{
		Cslice_iter t = *this;
		t.curr = s.size(); // index of one plus last element
		return t;
	}
	Cslice_iter& operator++() { curr++; return *this; }
	Cslice_iter operator++(int) { Cslice_iter t = *this; curr++; return t; }
	
	const T& operator[](size_t i) const { return ref(i); }
	const T& operator()(size_t i) const { return ref(i); }
	const T& operator*() const { return ref(curr); }

	friend bool operator==<>(const Cslice_iter& p, const Cslice_iter& q);
	friend bool operator!=<>(const Cslice_iter& p, const Cslice_iter& q);
	friend bool operator< <>(const Cslice_iter& p, const Cslice_iter& q);

};

template<class T>
bool operator==(const Cslice_iter<T>& p, const Cslice_iter<T>& q)
{
	return p.curr==q.curr
		&& p.s.stride()==q.s.stride()
		&& p.s.start()==q.s.start();
}

template<class T>
bool operator!=(const Cslice_iter<T>& p, const Cslice_iter<T>& q)
{
	return !(p==q);
}

template<class T>
bool operator<(const Cslice_iter<T>& p, const Cslice_iter<T>& q)
{
	return p.curr<q.curr
		&& p.s.stride()==q.s.stride()
		&& p.s.start()==q.s.start();
}

