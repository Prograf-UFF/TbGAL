#ifndef __TBGAL_HOMOGENEOUS_MACRO_HPP__
#define __TBGAL_HOMOGENEOUS_MACRO_HPP__

#define _TBGAL_OVERLOAD_HOMOGENEOUS_UTILS(METRIC_SPACE) \
    \
    template<typename... ScalarTypes> \
    constexpr decltype(auto) direction(ScalarTypes &&... coords) noexcept { \
        return tbgal::vector(METRIC_SPACE, std::move(coords)..., 0); \
    } \
    \
    template<typename... ScalarTypes> \
    constexpr decltype(auto) euclidean_vector(ScalarTypes &&... coords) noexcept { \
        return tbgal::vector(METRIC_SPACE, std::move(coords)..., 0); \
    } \
    \
    template<typename... ScalarTypes> \
    constexpr decltype(auto) point(ScalarTypes &&... coords) noexcept { \
        return tbgal::vector(METRIC_SPACE, std::move(coords)..., 1); \
    }

//TODO [NEXT] Inclue other helpers.
//TODO [NEXT] Inclue euclidean_vector with iterators.

#endif // __TBGAL_HOMOGENEOUS_MACRO_HPP__
