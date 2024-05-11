#include <iostream>
#include <ostream>
#include <vector>
#include <array>
template <class T, std::size_t N>
std::ostream& operator<<(std::ostream& o, const std::array<T, N>& arr)
{
  //std::copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(o, " "));
  for(auto const& i: arr) {
    o << i << " ";
  }
  o << std::endl;
  return o;
}
template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& vec)
{
  //copy(vec.cbegin(), vec.cend(), std::ostream_iterator<T>(o, " "));
  return o;
}
int main()
{
  std::array<int,4> a{0};
  std::cout << a;
  return 0;
}
