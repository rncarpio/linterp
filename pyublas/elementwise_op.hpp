//
// Copyright (c) 2004-2006
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




#ifndef HEADER_SEEN_PYUBLAS_ELEMENTWISE_OP_HPP
#define HEADER_SEEN_PYUBLAS_ELEMENTWISE_OP_HPP




#include <cmath>
#include <boost/numeric/ublas/functional.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>




namespace pyublas {
  template <class Functor>
  struct unary_op
  {
    template<class E> 
    static
    typename boost::numeric::ublas::vector_unary_traits<E, Functor>::result_type
    apply(const boost::numeric::ublas::vector_expression<E> &e) 
    {
      typedef typename boost::numeric::ublas::vector_unary_traits<E, Functor>
        ::expression_type expression_type;
      return expression_type(e());
    }

    template<class E>
    static
    typename boost::numeric::ublas::matrix_unary1_traits<E, Functor>::result_type
    apply(const boost::numeric::ublas::matrix_expression<E> &e) 
    {
      typedef typename boost::numeric::ublas::matrix_unary1_traits<E, Functor>
        ::expression_type expression_type;
      return expression_type(e());
    }
  };




  template <class Functor>
  struct binary_op
  {
    template<class E1, class E2> 
    static
    typename boost::numeric::ublas::vector_binary_traits<E1, E2, Functor>::result_type
    apply(
        const boost::numeric::ublas::vector_expression<E1> &e1,
        const boost::numeric::ublas::vector_expression<E2> &e2
        ) 
    {
      typedef typename boost::numeric::ublas::vector_binary_traits<E1, E2, Functor> 
        ::expression_type expression_type;
      return expression_type(e1(), e2());
    }

    template<class E1, class E2> 
    static
    typename boost::numeric::ublas::matrix_binary_traits<E1, E2, Functor>::result_type
    apply(
        const boost::numeric::ublas::matrix_expression<E1> &e1,
        const boost::numeric::ublas::matrix_expression<E2> &e2
        ) 
    {
      typedef typename boost::numeric::ublas::matrix_binary_traits<E1, E2, Functor> 
        ::expression_type expression_type;
      return expression_type(e1(), e2());
    }
  };




  namespace unary_ops
  {
    class fabs : public boost::numeric::ublas::scalar_real_unary_functor<double>
    {
      private:
        typedef boost::numeric::ublas::scalar_unary_functor<double> super;

      public:
        static super::result_type apply(super::argument_type x)
        {
          return ::std::fabs(x);
        }
    };
  }
  


  
  namespace binary_ops
  {
    template<class T1, class T2=T1>
    class max : public boost::numeric::ublas::scalar_binary_functor<T1, T2> 
    {
      private:
        typedef boost::numeric::ublas::scalar_binary_functor<T1, T2> super;

      public:
        static typename super::result_type apply(
            typename super::argument1_type t1, 
            typename super::argument1_type t2) 
        {
          if (t1 >= t2)
            return t1;
          else
            return t2;
        }
    };




    template<class T1, class T2=T1>
    class min : public boost::numeric::ublas::scalar_binary_functor<T1, T2> 
    {
      private:
        typedef boost::numeric::ublas::scalar_binary_functor<T1, T2> super;

      public:
        static typename super::result_type apply(
            typename super::argument1_type t1, 
            typename super::argument1_type t2) 
        {
          if (t1 <= t2)
            return t1;
          else
            return t2;
        }
    };
  }



  
  // square sum ---------------------------------------------------------------
  template<class V>
  struct vector_square_sum : 
    public boost::numeric::ublas::vector_scalar_real_unary_functor<V> 
  {
    private:
      typedef boost::numeric::ublas::vector_scalar_real_unary_functor<V> super;

    public:
      typedef typename super::real_type real_type;
      typedef typename super::value_type value_type;
      typedef typename super::result_type result_type;

      template<class E>
      static result_type 
      apply(const boost::numeric::ublas::vector_expression<E> &e) 
      {
        using namespace boost::numeric::ublas;

        real_type t = real_type ();
        typedef typename E::size_type vector_size_type;
        vector_size_type size (e ().size ());
        for (vector_size_type i = 0; i < size; ++ i) {
          typename super::real_type u (
              type_traits<value_type>::norm_2 (e () (i)));
          t +=  u * u;
        }
        return t;
      }

      // Dense case
      template<class D, class I>
      static result_type apply(D size, I it) 
      {
        using namespace boost::numeric::ublas;

        real_type t = real_type ();
        while (-- size >= 0) {
          real_type u (type_traits<value_type>::norm_2 (*it));
          t +=  u * u;
          ++ it;
        }
        return t;
      }

      // Sparse case
      template<class I>
      static result_type apply(I it, const I &it_end) 
      {
        using namespace boost::numeric::ublas;

        real_type t = real_type ();
        while (it != it_end) {
          real_type u (type_traits<value_type>::norm_2 (*it));
          t +=  u * u;
          ++ it;
        }
        return t;
      }
  };

  template<class E>
  inline 
  typename boost::numeric::ublas::vector_scalar_unary_traits<E, vector_square_sum<E> >::result_type
  square_sum(const boost::numeric::ublas::vector_expression<E> &e) 
  {
    typedef typename boost::numeric::ublas::
      vector_scalar_unary_traits<E, vector_square_sum<E> >::expression_type 
      expression_type;
    return expression_type(e());
  }
}




#endif
