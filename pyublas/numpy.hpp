//
// Copyright (c) 2008 Andreas Kloeckner
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee,
// provided that the above copyright notice appear in all copies and
// that both that copyright notice and this permission notice appear
// in supporting documentation.  The authors make no representations
// about the suitability of this software for any purpose.
// It is provided "as is" without express or implied warranty.
//




#ifndef _AFAYFYDASDFAH_PYUBLAS_HEADER_SEEN_NUMPY_HPP
#define _AFAYFYDASDFAH_PYUBLAS_HEADER_SEEN_NUMPY_HPP




#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4355 ) // 'this' : used in initializer list
#endif

#include <cstdlib>
#include <numeric>
#include <complex>
#include <pyublas/python_helpers.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_signed.hpp>
#include <boost/python.hpp>
#include <boost/python.hpp>
#include <boost/foreach.hpp>
#include <numpy/arrayobject.h>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_categories.hpp>
#include <boost/range.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_class.hpp>

namespace
{
  static struct pyublas_array_importer
  {
    static bool do_import_array()
    {
      import_array1(false);
      return true;
    }

    pyublas_array_importer()
    {
      if (!do_import_array())
        throw std::runtime_error("numpy failed to initialize");
    }
  } _array_importer;
}





namespace pyublas
{
  inline NPY_TYPES get_typenum(bool) { return NPY_BOOL; }
  // inline NPY_TYPES get_typenum(npy_bool) { return NPY_BOOL; }
  inline NPY_TYPES get_typenum(npy_byte) { return NPY_BYTE; }
  inline NPY_TYPES get_typenum(npy_ubyte) { return NPY_UBYTE; }
  inline NPY_TYPES get_typenum(npy_short) { return NPY_SHORT; }
  inline NPY_TYPES get_typenum(npy_ushort) { return NPY_USHORT; }
  inline NPY_TYPES get_typenum(npy_int) { return NPY_INT; }
  inline NPY_TYPES get_typenum(npy_uint) { return NPY_UINT; }
  inline NPY_TYPES get_typenum(npy_long) { return NPY_LONG; }
  inline NPY_TYPES get_typenum(npy_ulong) { return NPY_ULONG; }
  inline NPY_TYPES get_typenum(npy_longlong) { return NPY_LONGLONG; }
  inline NPY_TYPES get_typenum(npy_ulonglong) { return NPY_ULONGLONG; }
  inline NPY_TYPES get_typenum(npy_float) { return NPY_FLOAT; }
  inline NPY_TYPES get_typenum(npy_double) { return NPY_DOUBLE; }
  inline NPY_TYPES get_typenum(npy_cfloat) { return NPY_CFLOAT; }
  inline NPY_TYPES get_typenum(npy_cdouble) { return NPY_CDOUBLE; }
  inline NPY_TYPES get_typenum(std::complex<float>) { return NPY_CFLOAT; }
  inline NPY_TYPES get_typenum(std::complex<double>) { return NPY_CDOUBLE; }
#if HAVE_LONG_DOUBLE && (NPY_SIZEOF_LONGDOUBLE > NPY_SIZEOF_DOUBLE)
  inline NPY_TYPES get_typenum(npy_longdouble) { return NPY_LONGDOUBLE; }
  inline NPY_TYPES get_typenum(npy_clongdouble) { return NPY_CLONGDOUBLE; }
  inline NPY_TYPES get_typenum(std::complex<long double>) { return NPY_CLONGDOUBLE; }
#endif
  inline NPY_TYPES get_typenum(boost::python::object) { return NPY_OBJECT; }
  inline NPY_TYPES get_typenum(boost::python::handle<>) { return NPY_OBJECT; }
  /* NPY_STRING, NPY_UNICODE unsupported for now */

  template <class T>
  inline
  bool is_storage_compatible(PyObject *ary)
  {
    /* This piece of code works around the fact that 'int' and
     * 'long int' are the same on 32-bit machines, which can lead
     * to typenum mismatches. Therefore, for integers, we only
     * compare size and signedness.
     *
     * Also, bool and the chars are storage-compatible usually.
     */

    NPY_TYPES typenum = NPY_TYPES(PyArray_TYPE(ary));

    if (boost::is_integral<T>::value && PyArray_ISINTEGER(ary))
    {
      return (sizeof(T) == PyArray_ITEMSIZE(ary)
          && bool(boost::is_signed<T>::value)
          == bool(PyArray_ISSIGNED(ary)));
    }
    else if (typenum == NPY_BOOL && (
          boost::is_same<T, signed char>::value ||
          boost::is_same<T, unsigned char>::value))
    {
      return (sizeof(T) == PyArray_ITEMSIZE(ary)
          && bool(boost::is_signed<T>::value)
          == bool(PyArray_ISSIGNED(ary)));
    }
    else
      return typenum == get_typenum(T());
  }




  // tool functions -----------------------------------------------------------
  inline
  npy_intp size_from_dims(int ndim, const npy_intp *dims)
  {
    if (ndim != 0)
      return std::accumulate(dims, dims+ndim, 1, std::multiplies<npy_intp>());
    else
      return 1;
  }




  // ublas storage array ------------------------------------------------------
  template <class T>
  class numpy_array
  {
    private:
      // Life support for the numpy array.
      boost::python::handle<>         m_numpy_array;

    public:
      typedef std::size_t size_type;
      typedef std::ptrdiff_t difference_type;
      typedef T value_type;
      typedef const T &const_reference;
      typedef T &reference;
      typedef const T *const_pointer;
      typedef T *pointer;

      // Construction and destruction
      numpy_array()
      { }

      numpy_array(size_type n)
      {
        npy_intp dims[] = { n };
        m_numpy_array = boost::python::handle<>(
            PyArray_SimpleNew(1, dims, get_typenum(T())));
      }

      numpy_array(int ndim_, const npy_intp *dims_)
      {
        m_numpy_array = boost::python::handle<>(
            PyArray_SimpleNew(
              ndim_,
              const_cast<npy_intp *>(dims_),
              get_typenum(T())));
      }

      numpy_array(int ndim_, const npy_intp *dims_, T* data)
      {
        m_numpy_array = boost::python::handle<>(
            PyArray_SimpleNewFromData(
              ndim_,
              const_cast<npy_intp *>(dims_),
              get_typenum(T()), data));
      }

      numpy_array(size_type n, const value_type &v)
      {
        if (n)
        {
          npy_intp dims[] = { n };
          m_numpy_array = boost::python::handle<>(
              PyArray_SimpleNew(1, dims, get_typenum(T())));
          std::fill(begin(), end(), v);
        }
      }

      template<typename in_t>
      numpy_array(in_t const& in,
          typename boost::enable_if<boost::is_class<in_t> >::type *dummy=0)
      {
        if (boost::size(in))
        {
          npy_intp dims[] = { boost::size (in) };
          m_numpy_array = boost::python::handle<>(
              PyArray_SimpleNew(1, dims, get_typenum(T())));
          std::copy (boost::begin (in), boost::end (in), begin());
        }
      }

      numpy_array(const boost::python::handle<> &obj)
        : m_numpy_array(obj)
      {
        if (!obj.get())
          return;
        if (obj.get() == Py_None)
        {
          m_numpy_array = boost::python::handle<>();
          return;
        }

        if (!PyArray_Check(obj.get()))
          PYUBLAS_PYERROR(TypeError, "argument is not a numpy array");
        if (!is_storage_compatible<T>(obj.get()))
          PYUBLAS_PYERROR(TypeError, "argument is numpy array of wrong type");
        if (!PyArray_CHKFLAGS(obj.get(), NPY_ALIGNED))
            PYUBLAS_PYERROR(ValueError, "argument array is not aligned");
        if (PyArray_CHKFLAGS(obj.get(), NPY_NOTSWAPPED))
            PYUBLAS_PYERROR(ValueError, "argument array does not have native endianness");
        if (PyArray_ITEMSIZE(obj.get()) != sizeof(T))
            PYUBLAS_PYERROR(ValueError, "itemsize does not match size of target type");
      }

      numpy_array copy() const
      {
        boost::python::handle<> cp(PyArray_NewCopy(
              reinterpret_cast<PyArrayObject *>(m_numpy_array.get()), NPY_ANYORDER));
        return numpy_array(cp);
      }

    private:
      void resize_internal (size_type new_size, value_type init, bool preserve = true)
      {
        size_type old_size;
        if (m_numpy_array.get())
          old_size = size();
        else
        {
          preserve = false;
          old_size = 0;
        }

        if (new_size != old_size)
        {
          npy_intp dims[] = { new_size };
          boost::python::handle<> new_array = boost::python::handle<>(
              PyArray_SimpleNew(1, dims, get_typenum(T())));
          pointer new_data = reinterpret_cast<T *>(
              PyArray_DATA(new_array.get()));

          if (preserve)
          {
            std::copy(data(), data() + std::min(new_size, old_size), new_data);
            std::fill(new_data + std::min(new_size, old_size), new_data + new_size, init);
          }

          m_numpy_array = new_array;
        }
      }

    public:
      void resize (size_type size)
      {
        resize_internal (size, value_type(), false);
      }
      void resize (size_type size, value_type init)
      {
        resize_internal (size, init, true);
      }

      size_type size() const
      {
        if (!is_valid())
          return 0;

        if (ndim() != 0)
        {
          return end()-begin();
        }
        else
          return 1;
      }

      // metadata
      bool is_valid() const
      { return m_numpy_array.get(); }
      size_type ndim() const
      { return PyArray_NDIM(m_numpy_array.get()); }
      const npy_intp *dims() const
      { return PyArray_DIMS(m_numpy_array.get()); }
      const npy_intp dim(npy_intp i) const
      { return PyArray_DIM(m_numpy_array.get(), i); }
      const npy_intp *strides() const
      { return PyArray_STRIDES(m_numpy_array.get()); }
      const npy_intp stride(npy_intp i) const
      { return PyArray_STRIDE(m_numpy_array.get(), i); }

      npy_intp itemsize() const
      { return sizeof(T); }
      bool writable() const
      { return PyArray_ISWRITEABLE(m_numpy_array.get()); }

      // shape manipulation
      void reshape(int ndim_, const npy_intp *dims_,
          NPY_ORDER order=NPY_CORDER)
      {
        PyArray_Dims d = { const_cast<npy_intp *>(dims_), ndim_ };
        m_numpy_array = boost::python::handle<>(
            PyArray_Newshape(
              (PyArrayObject *) m_numpy_array.get(), &d, order));
      }

      // Raw data access
      T *data()
      {
        return reinterpret_cast<T *>(
            PyArray_DATA(m_numpy_array.get()));
      }

      const T *data() const
      {
        return reinterpret_cast<const T *>(
            PyArray_DATA(m_numpy_array.get()));
      }

      // Element access
      const_reference operator [] (size_type i) const
      {
        BOOST_UBLAS_CHECK(i < size(), boost::numeric::ublas::bad_index());
        return begin()[i];
      }

      reference operator [] (size_type i)
      {
        BOOST_UBLAS_CHECK(i < size(), boost::numeric::ublas::bad_index());
        return begin()[i];
      }

      // Assignment
      numpy_array &operator=(const numpy_array &a)
      {
        m_numpy_array = a.m_numpy_array;
        return *this;
      }

      numpy_array &assign_temporary(numpy_array &a)
      {
        m_numpy_array = a.m_numpy_array;
        return *this;
      }

        // Swapping
      void swap (numpy_array &a)
      {
        if (this != &a)
          std::swap(m_numpy_array, a.m_numpy_array);
      }

      friend void swap(numpy_array &a1, numpy_array &a2)
      {
        a1.swap (a2);
      }

      // Iterators simply are pointers.

      typedef const_pointer const_iterator;

    protected:

      npy_intp max_pos_stride_index() const
      {
        npy_intp current_idx = -1;
        npy_intp current_max = 0;
        for (unsigned i = 0; i < ndim(); ++i)
        {
          npy_intp si = stride(i);
          if (si > current_max)
          {
            current_max = si;
            current_idx = i;
          }
        }

        return current_idx;
      }

    public:
      const_iterator begin() const
      {
        const_iterator result = data();
        for (unsigned i = 0; i < ndim(); ++i)
        {
          const npy_intp si = stride(i)/npy_intp(sizeof(T));
          const npy_intp di = dim(i);
          if (si < 0 && di)
            result += si*(di-1);
        }

        return result;
      }

      const_iterator end() const
      {
        const npy_intp mpsi = max_pos_stride_index();

        if (mpsi != -1)
        {
          const npy_intp mps = stride(mpsi)/npy_intp(sizeof(T));
          return data() + mps*dim(mpsi);
        }
        else
          return data()+1;
      }

      typedef pointer iterator;

      iterator begin()
      {
        return const_cast<iterator>(
            const_cast<numpy_array const *>(this)->begin());
      }

      iterator end()
      {
        return const_cast<iterator>(
            const_cast<numpy_array const *>(this)->end());
      }

      // Reverse iterators
      typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
      typedef std::reverse_iterator<iterator> reverse_iterator;

      const_reverse_iterator rbegin() const
      { return const_reverse_iterator(end()); }

      const_reverse_iterator rend() const
      { return const_reverse_iterator(begin ()); }

      reverse_iterator rbegin()
      { return reverse_iterator(end()); }

      reverse_iterator rend ()
      { return reverse_iterator(begin()); }

      // Data accessor

      const boost::python::handle<> handle() const
      {
        if (is_valid())
          return m_numpy_array;
        else
          return boost::python::handle<>(
              boost::python::borrowed(Py_None));
      }

      boost::python::handle<> to_python() const
      {
        return handle();
      }
  };




  // matrix helper functions --------------------------------------------------
  inline bool is_row_major(boost::numeric::ublas::row_major_tag)
  {
    return true;
  }
  inline bool is_row_major(boost::numeric::ublas::column_major_tag)
  {
    return false;
  }

  template <class OCat, class T>
  typename numpy_array<T>::size_type get_array_size1(numpy_array<T> const &ary)
  {
    typedef numpy_array<T> mat_type;

    if (PyArray_NDIM(ary.handle().get()) != 2)
      throw std::runtime_error("ndarray->matrix converteee has dimension != 2");

    if (PyArray_STRIDE(ary.handle().get(), 1)
        == PyArray_ITEMSIZE(ary.handle().get()))
    {
      // row-major
      if (!is_row_major(OCat()))
        throw std::runtime_error("input array is not row-major (like the target type)");
      if (!PyArray_CHKFLAGS(ary.handle().get(), NPY_C_CONTIGUOUS))
        throw std::runtime_error("ndarray->matrix converteee is not C-contiguous");
    }
    else if (PyArray_STRIDE(ary.handle().get(), 0)
        == PyArray_ITEMSIZE(ary.handle().get()))
    {
      // column-major
      if (is_row_major(OCat()))
        throw std::runtime_error("input array is not column-major (like the target type)");
      if (!PyArray_CHKFLAGS(ary.handle().get(), NPY_F_CONTIGUOUS))
        throw std::runtime_error("ndarray->matrix converteee is not F-contiguous");
    }
    else
        throw std::runtime_error("input array is does not have dimension with stride==1");

    return PyArray_DIM(ary.handle().get(), 0);
  }

  template <class T>
  typename numpy_array<T>::size_type get_array_size2(numpy_array<T> const &ary)
  {
    // checking is done in size1()
    return PyArray_DIM(ary.handle().get(), 1);
  }

  template<class T, class L = boost::numeric::ublas::row_major>
  class numpy_matrix;

  template <class T, class L>
  boost::python::handle<> matrix_to_python(numpy_matrix<T, L> const &mat)
  {
    typedef numpy_matrix<T, L> mat_type;
    boost::python::handle<> orig_handle = mat.data().handle();

    npy_intp dims[] = { mat.size1(), mat.size2() };
    boost::python::handle<> result;

    if (is_row_major(typename mat_type::orientation_category()))
    {
      result = boost::python::handle<>(PyArray_New(
          &PyArray_Type, 2, dims,
          get_typenum(typename mat_type::value_type()),
          /*strides*/0,
          PyArray_DATA(orig_handle.get()),
          /* ? */ 0,
          NPY_CARRAY, NULL));
    }
    else
    {
      result = boost::python::handle<>(PyArray_New(
          &PyArray_Type, 2, dims,
          get_typenum(typename mat_type::value_type()),
          /*strides*/0,
          PyArray_DATA(orig_handle.get()),
          /* ? */ 0,
          NPY_FARRAY, NULL));
    }

    PyArray_BASE(result.get()) = boost::python::handle<>(orig_handle).release();
    return result;
  }



  // derived vector types -----------------------------------------------------
  template <class T>
  class numpy_strided_vector;

  template <class T>
  class numpy_vector;

  namespace detail
  {
    template <class Derived, class Super>
    struct vector_functionality
    {
      Derived copy() const
      {
        return Derived(static_cast<const Derived *>(this)->array());
      }

      // numpy array metadata
      bool is_valid() const
      { return static_cast<const Derived *>(this)->array().is_valid(); }
      typename Super::size_type ndim() const
      { return static_cast<const Derived *>(this)->array().ndim(); }
      const npy_intp *dims() const
      { return static_cast<const Derived *>(this)->array().dims(); }
      const npy_intp dim(npy_intp i) const
      { return static_cast<const Derived *>(this)->array().dim(i); }
      const npy_intp *strides() const
      { return static_cast<const Derived *>(this)->array().strides(); }
      const npy_intp stride(npy_intp i) const
      { return static_cast<const Derived *>(this)->array().stride(i); }
      npy_intp itemsize() const
      { return sizeof(typename Derived::value_type); }
      bool writable() const
      { return static_cast<const Derived *>(this)->array().writable(); }

      // several-d subscripts
      typename Super::value_type &sub(npy_intp i)
      { return *reinterpret_cast<typename Super::value_type *>(
          PyArray_GETPTR1(
            static_cast<const Derived *>(this)->array().handle().get(),
            i));
      }
      const typename Super::value_type &sub(npy_intp i) const
      {
        return *reinterpret_cast<const typename Super::value_type *>(
            PyArray_GETPTR1(
              static_cast<const Derived *>(this)->data().handle().get(),
              i));
      }
      typename Super::value_type &sub(npy_intp i, npy_intp j)
      {
        return *reinterpret_cast<typename Super::value_type *>(
            PyArray_GETPTR2(
              static_cast<const Derived *>(this)->array().handle().get(),
              i, j));
      }
      const typename Super::value_type &sub(npy_intp i, npy_intp j) const
      {
        return *reinterpret_cast<const typename Super::value_type *>(
            PyArray_GETPTR2(
              static_cast<const Derived *>(this)->array().handle().get(),
              i, j));
      }
      typename Super::value_type &sub(npy_intp i, npy_intp j, npy_intp k)
      {
        return *reinterpret_cast<typename Super::value_type *>(
            PyArray_GETPTR3(
              static_cast<const Derived *>(this)->array().handle().get(),
              i, j, k));
      }
      const typename Super::value_type &sub(npy_intp i, npy_intp j, npy_intp k) const
      {
        return *reinterpret_cast<const typename Super::value_type *>(
            PyArray_GETPTR3(
              static_cast<const Derived *>(this)->array().handle().get(),
              i, j, k));
      }
      typename Super::value_type &sub(npy_intp i, npy_intp j, npy_intp k, npy_intp l)
      {
        return *reinterpret_cast<typename Super::value_type *>(
            PyArray_GETPTR4(
              static_cast<const Derived *>(this)->array().handle().get(),
              i, j, k, l));
      }
      const typename Super::value_type &sub(npy_intp i, npy_intp j, npy_intp k, npy_intp l) const
      {
        return *reinterpret_cast<const typename Super::value_type *>(
          PyArray_GETPTR4(
            static_cast<const Derived *>(this)->array().handle().get(),
            i, j, k, l));
      }

      // shape manipulation
      void reshape(int ndim_, const npy_intp *dims_, NPY_ORDER order=NPY_CORDER)
      {
        static_cast<Derived *>(this)->data().reshape(ndim_, dims_, order);
      }

      boost::python::handle<> to_python() const
      {
        return static_cast<const Derived *>(this)->array().handle();
      }
    };



    template<class T>
    class numpy_vec_iterator : public boost::iterator_facade<numpy_vec_iterator<T>, T,
                                      boost::numeric::ublas::dense_random_access_iterator_tag>
    {
      private:
        typedef  boost::iterator_facade<numpy_vec_iterator, T,
                 boost::numeric::ublas::dense_random_access_iterator_tag>
                   self_t;

      public:
        typedef typename self_t::difference_type difference_type;

        numpy_vec_iterator()
        : it(0)
        { }

        numpy_vec_iterator(T *it)
        : it(it)
        { }

        // private:
        friend class boost::iterator_core_access;

        bool equal(numpy_vec_iterator const &other) const
        {
          return other.it == this->it;
        }

        void increment()
        { ++it; }

        T& dereference() const
        { return *it; }

        void decrement()
        { --it; }

        difference_type distance_to (numpy_vec_iterator const &other) const
        { return (other.it - this->it);}

        void advance (difference_type n)
        { it += n; }

        T* it;
    };
  } // end namespace detail




  template <class T>
  class numpy_vector
  : public boost::numeric::ublas::vector<T, numpy_array<T> >,
  public detail::vector_functionality<numpy_vector<T>,
  boost::numeric::ublas::vector<T, numpy_array<T> > >
  {
    private:
      typedef
        boost::numeric::ublas::vector<T, numpy_array<T> >
        super;
      typedef
        detail::vector_functionality<numpy_vector<T>,
          boost::numeric::ublas::vector<T, numpy_array<T> >
        >
        func;

    public:
      numpy_vector ()
        : super(0)
      { }

      // observe that PyObject handles are implicitly convertible
      // to numpy_array
      numpy_vector(const numpy_array<T> &s)
        : super(s.size(), s)
      { }

      numpy_vector(int ndim_, const npy_intp *dims_)
        : super(size_from_dims(ndim_, dims_),
            numpy_array<T>(ndim_, dims_))
      { }

      template<class AE>
      numpy_vector(int ndim_, const npy_intp *dims_,
          const boost::numeric::ublas::vector_expression<AE> &ae)
        : super(ae)
      {
        assert(this->size() == size_from_dims(ndim_, dims_));
        array().reshape(ndim_, dims_);
      }

      explicit
      numpy_vector(typename super::size_type size)
      : super(size)
      { }

      numpy_vector(
          typename super::size_type size,
          const typename super::value_type &init)
        : super(size, init)
      { }

      numpy_vector (const numpy_vector &v)
        : super(v)
      { }

      template<class AE>
      numpy_vector(const boost::numeric::ublas::vector_expression<AE> &ae)
      : super(ae)
      { }

      // as-ublas accessor
      super &as_ublas()
      { return *this; }

      const super &as_ublas() const
      { return *this; }

      boost::numeric::ublas::slice stride_slice() const
      {
        if (this->ndim() == 1)
        {
          npy_intp stride = this->stride(0)/npy_intp(sizeof(T));
          typename super::size_type sz = this->size();

          if (stride < 0)
          {
            return boost::numeric::ublas::slice(
                sz-1, stride, this->dim(0));
          }
          else
            return boost::numeric::ublas::slice(0, stride, this->dim(0));
        }
        else if (this->ndim() == 0)
          return boost::numeric::ublas::slice(0, 1, 1);
        else
          PYUBLAS_PYERROR(ValueError,
              "cannot generate a stride-respecting 1D slice "
              "for 2D or higher arrays");
      }

      numpy_strided_vector<T> as_strided()
      {
        return numpy_strided_vector<T>(*this, stride_slice());
      }

      boost::numeric::ublas::vector_slice<const numpy_vector>
        as_strided() const
      {
        return boost::numeric::ublas::vector_slice<const numpy_vector>
          (*this, stride_slice());
      }

      numpy_array<T> &array()
      { return super::data(); }

      numpy_array<T> const &array() const
      { return super::data(); }

      typedef detail::numpy_vec_iterator<T> iterator;
      typedef detail::numpy_vec_iterator<const T> const_iterator;

      iterator begin()
      { return iterator(array().begin()); }

      iterator end()
      { return iterator(array().end()); }

      const_iterator begin() const
      { return const_iterator(array().begin()); }

      const_iterator end() const
      { return const_iterator(array().end()); }

  };


  namespace detail
  {
    /* This is dumb, but necessary: In numpy_strided_vector, the
     * vector needs to be initialized before the slice referring
     * to it. But if the vector is a member, and the slice is a
     * base class, that won't happen. Therefore, move the vector
     * to an artificial base class that can be specified first in
     * the base-specifier-list. (cf. 12.6.2.5 in '96 working paper)
     */
    template <class V>
    class vector_holder
    {
      public:
               V m_vector;

        vector_holder(const V &v)
          : m_vector(v)
        { }
    };




    template<typename T>
    class numpy_strided_vec_iterator
    : public boost::iterator_facade<
      numpy_strided_vec_iterator<T>, T,
      boost::numeric::ublas::dense_random_access_iterator_tag>
    {
      private:
        typedef  boost::iterator_facade<
          numpy_strided_vec_iterator, T,
          boost::numeric::ublas::dense_random_access_iterator_tag>
        self_t;

      public:
        typedef typename self_t::difference_type difference_type;
        typedef typename numpy_strided_vector<T>::stride_type stride_type;

        numpy_strided_vec_iterator(T *it, stride_type s)
        : stride(s), it(it)
        { }

        // private:

        friend class boost::iterator_core_access;

        bool equal (numpy_strided_vec_iterator const& other) const
        { return other.it == this->it; }

        void increment()
        { it += stride; }

        T& dereference() const
        { return *it; }

        void decrement()
        { it -= stride; }

        difference_type
        distance_to (numpy_strided_vec_iterator const& other) const
        { return (other.it - this->it)/difference_type (stride);}

        void advance (difference_type n)
        { it += n*stride; }

      private:
        stride_type stride;
        T* it;
    };
  } // end namespace detail




  template <class T>
  class numpy_strided_vector
  : public detail::vector_holder<numpy_vector<T> >,
  public boost::numeric::ublas::vector_slice< numpy_vector<T> >,
  public detail::vector_functionality<numpy_strided_vector<T>,
    boost::numeric::ublas::vector_slice< numpy_vector<T> > >
  {
    private:
      typedef
        detail::vector_holder<numpy_vector<T> >
        vector_holder;
      typedef
        boost::numeric::ublas::vector_slice< numpy_vector<T> >
        super;

      // Make a fake 0-length array just for default constructor
      static numpy_array<T> make_fake_array() { return numpy_array<T>(0); }

    public:
      numpy_strided_vector()
        : vector_holder (make_fake_array()),
        super(this->m_vector, boost::numeric::ublas::slice(0, 1, this->m_vector.size()))
      { }

      // observe that PyObject handles are implicitly convertible
      // to numpy_array
      numpy_strided_vector(const numpy_array<T> &s)
        : vector_holder(s),
        super(this->m_vector, boost::numeric::ublas::slice(0, 1, s.size()))
      { }

      numpy_strided_vector(const numpy_strided_vector &v)
        : vector_holder(v.m_vector),
        super(this->m_vector, boost::numeric::ublas::slice(
              v.start(), v.stride(), v.size()))
      { }

      numpy_strided_vector(numpy_vector<T> &v, boost::numeric::ublas::slice const &s)
        : vector_holder(v), super(this->m_vector, s)
      { }

      template<class AE>
      numpy_strided_vector(const boost::numeric::ublas::vector_expression<AE> &ae)
      : vector_holder(ae),
      super(this->m_vector, boost::numeric::ublas::slice(0, 1, ae().size()))
      { }

      // as-ublas accessor
      super &as_ublas()
      { return *this; }

      const super &as_ublas() const
      { return *this; }

      numpy_array<T> &array()
      { return this->m_vector.data(); }

      numpy_array<T> const &array() const
      { return this->m_vector.data(); }

      typename super::difference_type stride() const
      { return super::stride(); }

      typename super::difference_type stride(npy_intp i) const
      { return this->m_vector.stride(i); }

      typedef npy_intp stride_type;

      typedef detail::numpy_strided_vec_iterator<T> iterator;
      typedef detail::numpy_strided_vec_iterator<const T> const_iterator;

      iterator begin()
      {
        if (stride() >= 0)
          return iterator(this->m_vector.array().begin(), stride());
        else
          return iterator(this->m_vector.array().end()-1, stride());
      }

      iterator end()
      {
        if (stride() >= 0)
          return iterator(this->m_vector.array().end(), stride());
        else
          return iterator(this->m_vector.array().end() - 1 - this->size(), stride());
      }

      const_iterator begin() const
      {
        if (stride() >= 0)
          return const_iterator(this->m_vector.array().begin(), stride());
        else
          return const_iterator(this->m_vector.array().end()-1, stride());
      }

      const_iterator end() const
      {
        if (stride() >= 0)
          return const_iterator(this->m_vector.array().end(), stride());
        else
          return const_iterator(this->m_vector.array().end() - 1 - this->size(), stride());
      }
  };




  // derived matrix types -----------------------------------------------------
  template<class T, class L/* = boost::numeric::ublas::row_major */>
    /* default arg declared in forward decl */
  class numpy_matrix
  : public boost::numeric::ublas::matrix<T, L, numpy_array<T> >
  {
    private:
      typedef
        boost::numeric::ublas::matrix<T, L, numpy_array<T> >
        super;

    public:
      numpy_matrix ()
      { }

      numpy_matrix(
          typename super::size_type size1,
          typename super::size_type size2)
      : super(size1, size2)
      { }

      numpy_matrix(
          typename super::size_type size1,
          typename super::size_type size2,
          const typename super::value_type &init)
      : super(size1, size2, init)
      { }

      numpy_matrix(
          typename super::size_type size1,
          typename super::size_type size2,
          const typename super::array_type &data)
      : super(size1, size2, data)
      { }

      // observe that PyObject handles are implicitly convertible
      // to numpy_array
      numpy_matrix(const typename super::array_type &data)
      : super(
          get_array_size1<typename super::orientation_category>(data),
          get_array_size2(data),
          data)
      { }

      numpy_matrix(const numpy_matrix &m)
      : super(m)
      { }

      template<class AE>
      numpy_matrix (const boost::numeric::ublas::matrix_expression<AE> &ae)
      : super(ae)
      {
      }

      numpy_matrix copy() const
      {
        return numpy_matrix(super::data().copy());
      }

      super &as_ublas()
      { return *this; }

      const super &as_ublas() const
      { return *this; }

      bool is_valid() const
      { return this->data().is_valid(); }
      boost::python::handle<> to_python() const
      { return matrix_to_python(*this); }
  };




  // conversion tags ----------------------------------------------------------
  template <class Contained>
  class invalid_ok
  {
    private:
      Contained m_contained;

    public:
      typedef Contained contained_type;

      invalid_ok(Contained c)
        : m_contained(c)
      { }

      const Contained &operator*() const
      {
        return m_contained;
      }

      Contained &operator*()
      {
        return m_contained;
      }

      const Contained *operator->() const
      {
        return &m_contained;
      }

      Contained *operator->()
      {
        return &m_contained;
      }
  };




  // data member treatment ----------------------------------------------------
  template <class T, class C>
  class by_value_rw_member_visitor
  : public boost::python::def_visitor<by_value_rw_member_visitor<T, C> >
  {
    private:
      const char *m_name;
      T C::*m_member;
      const char *m_doc;

    public:
      by_value_rw_member_visitor(const char *name, T C::*member, const char *doc = 0)
        : m_name(name), m_member(member), m_doc(doc)
      { }

      template <class Class>
      void visit(Class& cl) const
      {
        cl.add_property(m_name,
            boost::python::make_getter(m_member,
              boost::python::return_value_policy<boost::python::return_by_value>()),
            boost::python::make_setter(m_member),
            m_doc);
      }
  };

  template <class T, class C>
  by_value_rw_member_visitor<T, C> by_value_rw_member(
      const char *name, T C::*member, const char *doc = 0)
  {
    return by_value_rw_member_visitor<T, C>(name, member, doc);
  }

  template <class T, class C>
  class by_value_ro_member_visitor
  : public boost::python::def_visitor<by_value_ro_member_visitor<T, C> >
  {
    private:
      const char *m_name;
      T C::*m_member;
      const char *m_doc;

    public:
      by_value_ro_member_visitor(const char *name, T C::*member, const char *doc = 0)
        : m_name(name), m_member(member), m_doc(doc)
      { }

      template <class Class>
      void visit(Class& cl) const
      {
        cl.add_property(m_name,
            make_getter(m_member,
              boost::python::return_value_policy<boost::python::return_by_value>()),
            m_doc);
      }
  };

  template <class T, class C>
  by_value_ro_member_visitor<T, C> by_value_ro_member(
      const char *name, T C::*member, const char *doc = 0)
  {
    return by_value_ro_member_visitor<T, C>(name, member, doc);
  }
}




// interaction with boost bindings --------------------------------------------
#ifdef PYUBLAS_HAVE_BOOST_BINDINGS

#include <boost/numeric/bindings/traits/ublas_vector.hpp>




namespace boost { namespace numeric { namespace bindings { namespace traits {
  template <typename T, typename V>
  struct vector_detail_traits< pyublas::numpy_array<T>, V >
  : default_vector_traits< V, T >
  {
#ifndef BOOST_NUMERIC_BINDINGS_NO_SANITY_CHECK
    BOOST_STATIC_ASSERT(
        (boost::is_same< pyublas::numpy_array<T>,
         typename boost::remove_const<V>::type >::value) );
#endif

    typedef pyublas::numpy_array<T>                      identifier_type;
    typedef V                                            vector_type;
    typedef typename default_vector_traits<V,T>::pointer pointer;

    static pointer storage (vector_type& v) { return v.data(); }
  };
}}}}

#endif



#ifdef _MSC_VER
#pragma warning( pop )
#endif

#endif
// EMACS-FORMAT-TAG
//
// Local Variables:
// mode: C++
// eval: (c-set-style "bsd")
// eval: (c-set-offset 'access-label -2)
// eval: (c-set-offset 'inclass '++)
// c-basic-offset: 2
// tab-width: 8
// End:
