#ifndef __TBGAL_CONFORMAL_MACRO_HPP__
#define __TBGAL_CONFORMAL_MACRO_HPP__

#define _TBGAL_OVERLOAD_CONFORMAL_UTILS(METRIC_SPACE) \
    template<typename... ScalarTypes, typename = std::enable_if_t<std::disjunction_v<std::bool_constant<!detail::is_iterator_v<ScalarTypes> >...> > > \
    constexpr decltype(auto) euclidean_vector(ScalarTypes &&... coords) noexcept { \
        return tbgal::vector(METRIC_SPACE, std::move(coords)..., 0, 0); \
    }

//TODO [NEXT] Inclue other helpers.
//TODO [NEXT] Inclue euclidean_vector with iterators.

#endif // __TBGAL_CONFORMAL_MACRO_HPP__
