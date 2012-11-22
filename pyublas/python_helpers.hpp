//
// Copyright (c) 2008
// Andreas Kloeckner
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee,
// provided that the above copyright notice appear in all copies and
// that both that copyright notice and this permission notice appear
// in supporting documentation.  The authors make no representations
// about the suitability of this software for any purpose.
// It is provided "as is" without express or implied warranty.
//




#ifndef _AFAYFYDASDFAH_PYUBLAS_HEADER_SEEN_PYTHON_HELPERS_HPP
#define _AFAYFYDASDFAH_PYUBLAS_HEADER_SEEN_PYTHON_HELPERS_HPP




#include <boost/python.hpp>




#define PYUBLAS_PYERROR(TYPE, REASON) \
{ \
  PyErr_SetString(PyExc_##TYPE, REASON); \
  throw boost::python::error_already_set(); \
}




namespace pyublas
{
  template <typename T>
  inline PyObject *pyobject_from_new_ptr(T *ptr)
  {
    return typename boost::python::manage_new_object::apply<T *>::type()(ptr);
  }




  template <typename T>
  inline boost::python::handle<> handle_from_new_ptr(T *ptr)
  {
    return boost::python::handle<>(
        typename boost::python::manage_new_object::apply<T *>::type()(ptr));
  }




  template <typename T>
  inline boost::python::handle<> handle_from_existing_ptr(T *ptr)
  {
    return boost::python::handle<>(
        typename boost::python::reference_existing_object::apply<T *>::type()(ptr)
        );
  }




  template <typename T>
  inline boost::python::handle<> handle_from_existing_ref(T &ptr)
  {
    return boost::python::handle<>(
        typename boost::python::reference_existing_object::apply<T &>::type()(ptr)
        );
  }




  inline boost::python::handle<> handle_from_object(const boost::python::object &obj)
  {
    return boost::python::handle<>(boost::python::borrowed(obj.ptr()));
  }




  template <typename T>
  inline boost::python::handle<> handle_from_rvalue(const T &val)
  {
    boost::python::object obj(val);
    return handle_from_object(obj);
  }
}




#endif

