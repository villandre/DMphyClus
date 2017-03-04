#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

#ifndef HELPER_H
#define HELPER_H

inline uint littleCycle(uint myInt, uint cycleLength) 
{
  return myInt % cycleLength ;
}

template <typename... T>
auto zip(const T&... containers) -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
{
  auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
  auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
  return boost::make_iterator_range(zip_begin, zip_end);
}

template <typename T>
void deallocate_container(T& c)
{
  for (typename T::iterator i = c.begin(); i != c.end(); ++i)
    delete *i; 
}

#endif /* HELPER_H */