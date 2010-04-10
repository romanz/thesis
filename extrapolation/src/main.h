#ifdef BUILD_DLL
    #define DLL_EXPORT __declspec(dllexport)
#else
    #define DLL_EXPORT __declspec(dllimport)
#endif

typedef double real_t;
typedef unsigned long int uint_t;

enum EXTR_METHOD { MPE, RRE };
int DLL_EXPORT extrapolate(real_t *vectors,
    uint_t dim, uint_t num, enum EXTR_METHOD method);
