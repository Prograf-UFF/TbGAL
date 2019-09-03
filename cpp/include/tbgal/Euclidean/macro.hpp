#ifndef __TBGAL_EUCLIDEAN_MACRO_HPP__
#define __TBGAL_EUCLIDEAN_MACRO_HPP__

#define _TBGAL_OVERLOAD_EUCLIDEAN_UTILS(METRIC_SPACE) \
    \
    template<typename... ScalarTypes> \
    constexpr decltype(auto) euclidean_vector(ScalarTypes &&... coords) noexcept { \
        return tbgal::vector(METRIC_SPACE, std::move(coords)...); \
    }

#endif // __TBGAL_EUCLIDEAN_MACRO_HPP__
