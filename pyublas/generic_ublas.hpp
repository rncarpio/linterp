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




#ifndef HEADER_SEEN_PYUBLAS_GENERIC_UBLAS_HPP
#define HEADER_SEEN_PYUBLAS_GENERIC_UBLAS_HPP




#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>




namespace pyublas {
  namespace detail
  {
    class begin_tag { };
    class end_tag { };
  }




  template <typename ValueType, unsigned MaxElements = 2>
  class minilist
  {
  public:
    typedef ValueType value_type;
    typedef unsigned size_type;

  private:
    ValueType m_list[MaxElements];
    unsigned m_size;

  public:
    minilist()
    : m_size(0)
    { }

    minilist(const value_type &v0)
    : m_size(1)
    { 
      m_list[0] = v0;
    }

    minilist(const value_type &v0, const value_type &v1)
    : m_size(2)
    { 
      m_list[0] = v0;
      m_list[1] = v1;
    }

    size_type size() const
    {
      return m_size;
    }

    void push_back(const value_type &v)
    {
#ifndef NDEBUG
      if (m_size == MaxElements)
        throw std::runtime_error("minilist has reached max size");
#endif
      m_list[m_size++] = v;
    }

    ValueType &operator[](size_type index)
    {
      return m_list[index];
    }

    const ValueType &operator[](size_type index) const
    {
      return m_list[index];
    }
  };



  // is_vector ----------------------------------------------------------------
  template <typename UblasType>
    struct is_vector { typedef boost::mpl::false_ type; };

  template <typename ValueType, typename L>
    struct is_vector<boost::numeric::ublas::vector<ValueType, L> > { typedef boost::mpl::true_ type; };
  template <typename WrappedVector>
    struct is_vector<boost::numeric::ublas::vector_slice<WrappedVector> > { typedef boost::mpl::true_ type; };
  template <typename WrappedVector>
    struct is_vector<boost::numeric::ublas::matrix_row<WrappedVector> > { typedef boost::mpl::true_ type; };
  template <typename WrappedVector>
    struct is_vector<boost::numeric::ublas::matrix_column<WrappedVector> > { typedef boost::mpl::true_ type; };




  // matrix_iterator ----------------------------------------------------------
  template <typename MatrixType, typename _is_vector = typename is_vector<MatrixType>::type>
  class matrix_iterator :  public boost::iterator_facade<

    matrix_iterator<MatrixType>,  // Base

    typename MatrixType::value_type, // Value

    boost::forward_traversal_tag, // CategoryOrTraversal

    typename boost::mpl::if_<boost::is_const<MatrixType>,
      typename MatrixType::const_iterator2::reference,
      typename MatrixType::iterator2::reference
      >::type // Reference
   > 
  {
    typedef typename boost::mpl::if_<boost::is_const<MatrixType>,
      typename MatrixType::const_iterator1,
      typename MatrixType::iterator1
      >::type it1_t;
    typedef typename boost::mpl::if_<boost::is_const<MatrixType>,
      typename MatrixType::const_iterator2,
      typename MatrixType::iterator2
      >::type it2_t;

    it1_t       m_it1;
    it2_t       m_it2;

  public:
    matrix_iterator() { }

    matrix_iterator(MatrixType &mat, detail::begin_tag)
    : m_it1(mat.begin1()), m_it2(m_it1.begin())
    { 
      validate();
    }

    matrix_iterator(MatrixType &mat, detail::end_tag)
    : m_it1(mat.end1()), m_it2(m_it1.begin())
    { }

    minilist<typename MatrixType::size_type> index() const
    {
      return minilist<typename MatrixType::size_type>(m_it2.index1(), m_it2.index2());
    }

  private:
    friend class boost::iterator_core_access;

    void validate()
    {
      // this makes sure that the iterator points to an existing element
      while (m_it1 != m_it1().end1() && m_it2 == m_it1.end())
      {
        m_it1++;
        m_it2 = m_it1.begin();
      }
    }

    void increment() 
    {
      m_it2++;
      validate();
    }

    bool equal(matrix_iterator const& other) const
    {
      return m_it1 == other.m_it1 && m_it2 == other.m_it2;
    }

    typename it2_t::reference dereference() const 
    {
      return *m_it2; 
    }
  };




  template <typename MatrixType>
  class matrix_iterator<MatrixType, boost::mpl::true_> :  public boost::iterator_adaptor<
    matrix_iterator<MatrixType>, 
    typename MatrixType::iterator,
    typename MatrixType::value_type,
    boost::forward_traversal_tag,
    typename MatrixType::iterator::reference>
  {
    typedef
      boost::iterator_adaptor<
      matrix_iterator<MatrixType>, 
      typename MatrixType::iterator,
      typename MatrixType::value_type,
      boost::forward_traversal_tag,
      typename MatrixType::iterator::reference>
        super;

  public:
    matrix_iterator() { }

    matrix_iterator(MatrixType &mat, detail::begin_tag)
    : super(mat.begin())
    { }

    matrix_iterator(MatrixType &mat, detail::end_tag)
    : super(mat.end())
    { }

    matrix_iterator(const typename MatrixType::iterator &it)
    : super(it)
    { }

    minilist<unsigned> index() const
    {
      return minilist<unsigned>(this->base_reference().index());
    }

  private:
    friend class boost::iterator_core_access;
  };




  template <typename MatrixType>
  matrix_iterator<MatrixType> begin(MatrixType &mat)
  {
    return matrix_iterator<MatrixType>(mat, detail::begin_tag());
  }

  template <typename MatrixType>
  matrix_iterator<MatrixType> end(MatrixType &mat)
  {
    return matrix_iterator<MatrixType>(mat, detail::end_tag());
  }




  // shapes and subscripting --------------------------------------------------
  namespace detail
  {
    template <typename MatrixType>
    inline minilist<unsigned> getShape(const MatrixType &mat, boost::mpl::false_)
    {
      return minilist<unsigned>(mat.size1(), mat.size2());
    }

    template <typename MatrixType>
    inline minilist<unsigned> getShape(const MatrixType &mat, boost::mpl::true_)
    {
      return minilist<unsigned>(mat.size());
    }
  }

  template <typename MatrixType>
  inline minilist<unsigned> getShape(const MatrixType &mat)
  {
    return detail::getShape(mat, typename is_vector<MatrixType>::type());
  }




  namespace detail
  {
    template <typename MatrixType, typename IdxType>
    inline void setShape(MatrixType &mat, const minilist<IdxType> &shape, boost::mpl::false_)
    {
      mat.resize(shape[0], shape[1]);
    }

    template <typename MatrixType, typename IdxType>
    inline void setShape(MatrixType &mat, const minilist<IdxType> &shape, boost::mpl::true_)
    {
      mat.resize(shape[0]);
    }
  }

  template <typename MatrixType, typename IdxType>
  inline void setShape(MatrixType &mat, const minilist<IdxType> &shape)
  {
    detail::setShape(mat, shape, typename is_vector<MatrixType>::type());
  }






  namespace detail
  {
    template <typename MatrixType, typename IdxType>
    inline MatrixType *newWithShape(const minilist<IdxType> &shape, boost::mpl::false_)
    {
      return new MatrixType(shape[0], shape[1]);
    }

    template <typename MatrixType, typename IdxType>
    inline MatrixType *newWithShape(const minilist<IdxType> &shape, boost::mpl::true_)
    {
      return new MatrixType(shape[0]);
    }
  }



  template <typename MatrixType, typename IdxType>
  MatrixType *newWithShape(const minilist<IdxType> &shape)
  {
    return detail::newWithShape<MatrixType>(shape, typename is_vector<MatrixType>::type());
  }




  namespace detail
  {
    template <typename MatrixType, typename IdxType>
    inline void insert_element(
        MatrixType &mat,
        const minilist<IdxType> &index, 
        const typename MatrixType::value_type &value,
        boost::mpl::false_)
    {
      mat.insert_element(index[0], index[1], value);
    }

    template <typename MatrixType, typename IdxType>
    inline void insert_element(
        MatrixType &mat,
        const minilist<IdxType> &index, 
        const typename MatrixType::value_type &value,
        boost::mpl::true_)
    {
      mat.insert_element(index[0], value);
    }
  }

  template <typename MatrixType, typename IdxType>
  inline void insert_element(
      MatrixType &mat,
      const minilist<IdxType> &index, 
      const typename MatrixType::value_type &value)
  {
    detail::insert_element(mat,index, value, typename is_vector<MatrixType>::type());
  }




  namespace detail
  {
    template <typename MatrixType, typename IdxType>
    inline void set(
        MatrixType &mat,
        const minilist<IdxType> &index, 
        const typename MatrixType::value_type &value,
        boost::mpl::false_)
    {
      mat(index[0], index[1]) = value;
    }

    template <typename MatrixType, typename IdxType>
    inline void set(
        MatrixType &mat,
        const minilist<IdxType> &index, 
        const typename MatrixType::value_type &value,
        boost::mpl::true_)
    {
      mat[index[0]] = value;
    }
  }

  template <typename MatrixType, typename IdxType>
  inline void set(
      MatrixType &mat,
      const minilist<IdxType> &index, 
      const typename MatrixType::value_type &value)
  {
    detail::set(mat,index, value, typename is_vector<MatrixType>::type());
  }
}




#endif

