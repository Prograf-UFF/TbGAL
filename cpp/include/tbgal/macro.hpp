#ifndef __TBGAL_MACRO_HPP__
#define __TBGAL_MACRO_HPP__

#define _TBGAL_OVERLOAD_UTILS(METRIC_SPACE) \
    \
    decltype(auto) e(DefaultIndexType index) noexcept { \
        return tbgal::e<tbgal::DefaultScalarType>(METRIC_SPACE, index); \
    } \
    \
    template<typename ScalarType> \
    decltype(auto) scalar(ScalarType &&scalar) noexcept { \
        return tbgal::scalar(METRIC_SPACE, std::move(scalar)); \
    } \
    \
    template<typename... ScalarTypes, typename = std::enable_if_t<std::disjunction_v<std::bool_constant<!detail::is_iterator_v<ScalarTypes> >...> > > \
    decltype(auto) vector(ScalarTypes &&... coords) noexcept { \
        return tbgal::vector(METRIC_SPACE, std::move(coords)...); \
    } \
    \
    template<typename IteratorType, typename = std::enable_if_t<detail::is_iterator_v<IteratorType> > > \
    decltype(auto) vector(IteratorType begin, IteratorType end) noexcept { \
        return tbgal::vector(METRIC_SPACE, begin, end); \
    }

#endif // __TBGAL_MACRO_HPP__
