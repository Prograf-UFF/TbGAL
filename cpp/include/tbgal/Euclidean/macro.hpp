#ifndef __TBGAL_EUCLIDEAN_MACRO_HPP__
#define __TBGAL_EUCLIDEAN_MACRO_HPP__

#define _TBGAL_OVERLOAD_EUCLIDEAN_UTILS(METRIC_SPACE) \
    \
    template<typename... ScalarTypes, typename = std::enable_if_t<std::disjunction_v<std::bool_constant<!detail::is_iterator_v<ScalarTypes> >...> > > \
    constexpr decltype(auto) euclidean_vector(ScalarTypes &&... coords) noexcept { \
        return tbgal::vector(METRIC_SPACE, std::move(coords)...); \
    } \
    \
    template<typename IteratorType, typename = std::enable_if_t<detail::is_iterator_v<IteratorType> > > \
    constexpr decltype(auto) euclidean_vector(IteratorType begin, IteratorType end) noexcept { \
        return tbgal::vector(METRIC_SPACE, begin, end); \
    }

#endif // __TBGAL_EUCLIDEAN_MACRO_HPP__
