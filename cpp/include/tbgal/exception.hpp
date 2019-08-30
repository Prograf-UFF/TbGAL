#ifndef __TBGAL_EXCEPTION_HPP__
#define __TBGAL_EXCEPTION_HPP__

namespace tbgal {

    class NotSupportedError : public std::logic_error {
    public:

        explicit NotSupportedError(const std::string &what_arg) :
            std::logic_error(what_arg) {
        }

        explicit NotSupportedError(const char *what_arg) :
            std::logic_error(what_arg) {
        }
    };

}

#endif // __TBGAL_EXCEPTION_HPP__
