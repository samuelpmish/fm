#if 0
namespace fm {

template < Kind kind, u32 dim, typename T >
constexpr auto some_function(const matrix<kind, dim, dim, T> & A) {

  if constexpr (kind == Kind::General) {
    static_assert(always_false<T>{}, "unimplemented");
  }

  if constexpr (kind == Kind::Symmetric) {
    static_assert(always_false<T>{}, "unimplemented");
  }

  if constexpr (kind == Kind::Rotation) {
    static_assert(always_false<T>{}, "unimplemented");
  }

  if constexpr (kind == Kind::Skew) {
    static_assert(always_false<T>{}, "unimplemented");
  }

  if constexpr (kind == Kind::Diagonal) {
    static_assert(always_false<T>{}, "unimplemented");
  }

  if constexpr (kind == Kind::Isotropic) {
    static_assert(always_false<T>{}, "unimplemented");
  }

}

}
#endif