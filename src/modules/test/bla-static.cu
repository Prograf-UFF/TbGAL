#include <boost/python.hpp>


class teste {
public:
    static void set_var(int var) {
        var_ = var;
    }
    int v2;
    teste() {
        if (var_ == 0) {
            throw 0;
        }
        this->v2 = var_ * 2;
    }
    static int var_;
};

int teste::var_ = 0;


/*
static void set_var(int var) {
    teste::var = var;
}
*/
BOOST_PYTHON_MODULE(teste)
{

    namespace python = boost::python;

    python::class_<teste, boost::noncopyable>("Teste")
        .def("set_var", &teste::set_var)
            .staticmethod("set_var");
//    python::def("set_var", set_var<int>);


}
