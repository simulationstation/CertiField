// Evaluate Boost.Math functions with SciPy's two Boost policies, one call per input line.
//   input:  <name> <arg1> [<arg2> ...]          (decimal or hex-float arguments)
//   output: <name> <args...> -> <value %.17g>   or   -> EXC <message>
// Built against two Boost.Math commits (the one bundled by SciPy, and develop) by
// boostcheck.py, to tell whether a SciPy defect is already fixed upstream.
#define BOOST_MATH_STANDALONE
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>

using namespace boost::math::policies;

// scipy/special/boost_special_functions.h (SciPy 1.18.1)
typedef policy<promote_float<false>, promote_double<false>, max_root_iterations<400>,
               discrete_quantile<real>> SpecialPolicy;
typedef policy<domain_error<ignore_error>, overflow_error<user_error>, evaluation_error<user_error>,
               promote_float<false>, promote_double<false>,
               discrete_quantile<integer_round_up>> StatsPolicy;

static std::string g_warn;
namespace boost { namespace math { namespace policies {
template <class T> T user_evaluation_error(const char*, const char* msg, const T& val)
{ g_warn += std::string(" WARN-eval:") + (msg ? msg : ""); return val; }
template <class T> T user_overflow_error(const char*, const char* msg, const T&)
{ g_warn += std::string(" WARN-overflow:") + (msg ? msg : ""); return 0; }
}}}

static double eval(const std::string& f, const std::vector<double>& a)
{
    namespace bm = boost::math;
    SpecialPolicy sp;
    if (f == "ibeta") return bm::ibeta(a[0], a[1], a[2], sp);
    if (f == "ibetac") return bm::ibetac(a[0], a[1], a[2], sp);
    if (f == "ibeta_inv") return bm::ibeta_inv(a[0], a[1], a[2], sp);
    if (f == "ibetac_inv") return bm::ibetac_inv(a[0], a[1], a[2], sp);
    if (f == "gamma_p") return bm::gamma_p(a[0], a[1], sp);
    if (f == "gamma_q") return bm::gamma_q(a[0], a[1], sp);
    if (f == "gamma_p_inv") return bm::gamma_p_inv(a[0], a[1], sp);
    if (f == "gamma_q_inv") return bm::gamma_q_inv(a[0], a[1], sp);
    if (f == "lgamma_q") return bm::lgamma_q(a[0], a[1], sp);
    if (f == "lgamma_p") return bm::lgamma_p(a[0], a[1], sp);
    if (f == "erf_inv") return bm::erf_inv(a[0], sp);
    if (f == "erfc_inv") return bm::erfc_inv(a[0], sp);
    if (f == "hyp1f1") return bm::hypergeometric_1F1(a[0], a[1], a[2], sp);
    if (f == "t_cdf") return bm::cdf(bm::students_t_distribution<double, SpecialPolicy>(a[0]), a[1]);
    if (f == "t_ppf") return bm::quantile(bm::students_t_distribution<double, SpecialPolicy>(a[0]), a[1]);
    if (f == "f_cdf") return bm::cdf(bm::fisher_f_distribution<double, SpecialPolicy>(a[0], a[1]), a[2]);
    if (f == "f_sf") return bm::cdf(bm::complement(bm::fisher_f_distribution<double, SpecialPolicy>(a[0], a[1]), a[2]));
    if (f == "beta_ppf") return bm::quantile(bm::beta_distribution<double, StatsPolicy>(a[1], a[2]), a[0]);
    if (f == "beta_isf") return bm::quantile(bm::complement(bm::beta_distribution<double, StatsPolicy>(a[1], a[2]), a[0]));
    if (f == "beta_sf") return bm::cdf(bm::complement(bm::beta_distribution<double, StatsPolicy>(a[1], a[2]), a[0]));
    if (f == "binom_sf") return bm::cdf(bm::complement(bm::binomial_distribution<double, StatsPolicy>(a[1], a[2]), a[0]));
    if (f == "binom_cdf") return bm::cdf(bm::binomial_distribution<double, StatsPolicy>(a[1], a[2]), a[0]);
    if (f == "binom_isf") return bm::quantile(bm::complement(bm::binomial_distribution<double, StatsPolicy>(a[1], a[2]), a[0]));
    if (f == "binom_ppf") return bm::quantile(bm::binomial_distribution<double, StatsPolicy>(a[1], a[2]), a[0]);
    throw std::invalid_argument("unknown function " + f);
}

int main()
{
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream in(line);
        std::string f, tok;
        if (!(in >> f)) continue;
        std::vector<double> a;
        while (in >> tok) a.push_back(std::strtod(tok.c_str(), nullptr));
        g_warn.clear();
        std::cout << f;
        for (double v : a) { char b[40]; std::snprintf(b, sizeof b, " %.17g", v); std::cout << b; }
        try {
            double r = eval(f, a);
            char b[40]; std::snprintf(b, sizeof b, "%.17g", r);
            std::cout << " -> " << b << g_warn << "\n";
        } catch (const std::exception& e) {
            std::cout << " -> EXC " << typeid(e).name() << ": " << e.what() << "\n";
        }
        std::cout.flush();
    }
    return 0;
}
