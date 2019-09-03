#ifndef __TBGAL_MACRO_HPP__
#define __TBGAL_MACRO_HPP__

#define _TBGAL_OVERLOAD_UTILS(METRIC_SPACE) \
    \
    constexpr decltype(auto) e(std::int32_t index) noexcept { \
        return tbgal::e<tbgal::DefaultScalarType>(METRIC_SPACE, index); \
    } \
    \
    template<typename ScalarType> \
    constexpr decltype(auto) scalar(ScalarType &&scalar) noexcept { \
        return tbgal::scalar(METRIC_SPACE, std::move(scalar)); \
    } \
    \
    template<typename... ScalarTypes> \
    constexpr decltype(auto) vector(ScalarTypes &&... coords) noexcept { \
        return tbgal::vector(METRIC_SPACE, std::move(coords)...); \
    }


#endif // __TBGAL_MACRO_HPP__
