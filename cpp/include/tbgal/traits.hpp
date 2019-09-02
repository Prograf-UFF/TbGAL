#ifndef __TBGAL_TRAITS_HPP__
#define __TBGAL_TRAITS_HPP__

namespace tbgal {
    
    namespace detail {

        template<typename T, typename... Rest>
        constexpr bool is_any_v = std::disjunction_v<std::bool_constant<std::is_same_v<T, Rest> >...>;

        template<typename Type>
        struct is_factoring_product;

        template<typename Type>
        using is_factoring_product_t = typename is_factoring_product<Type>::type;

        template<typename Type>
        constexpr bool is_factoring_product_v = is_factoring_product<Type>::value;

        template<typename Type>
        struct is_metric_space;

        template<typename Type>
        using is_metric_space_t = typename is_metric_space<Type>::type;

        template<typename Type>
        constexpr bool is_metric_space_v = is_metric_space<Type>::value;

        template<typename Type>
        struct is_multivector;

        template<typename Type>
        using is_multivector_t = typename is_multivector<Type>::type;

        template<typename Type>
        constexpr bool is_multivector_v = is_multivector<Type>::value;

    }

}

#endif // __TBGAL_TRAITS_HPP__
