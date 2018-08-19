/*
    conStruct is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    conStruct is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with conStruct.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.17.0

#include <stan/model/model_header.hpp>

namespace model_oneK_namespace {

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;

typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

static int current_statement_begin__;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_oneK");
    reader.add_event(34, 34, "end", "model_oneK");
    return reader;
}

template <typename T1__, typename T2__>
Eigen::Matrix<typename boost::math::tools::promote_args<T1__, T2__>::type, Eigen::Dynamic,Eigen::Dynamic>
Cov(const int& N,
        const Eigen::Matrix<T1__, Eigen::Dynamic,1>& nugget,
        const T2__& gamma, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T1__, T2__>::type fun_scalar_t__;
    typedef fun_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        fun_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 3;
        validate_non_negative_index("parCov", "N", N);
        validate_non_negative_index("parCov", "N", N);
        Eigen::Matrix<fun_scalar_t__,Eigen::Dynamic,Eigen::Dynamic>  parCov(static_cast<Eigen::VectorXd::Index>(N),static_cast<Eigen::VectorXd::Index>(N));
        (void) parCov;  // dummy to suppress unused var warning

        stan::math::initialize(parCov, std::numeric_limits<double>::quiet_NaN());
        stan::math::fill(parCov,DUMMY_VAR__);
        current_statement_begin__ = 4;
        validate_non_negative_index("Nug_mat", "N", N);
        validate_non_negative_index("Nug_mat", "N", N);
        Eigen::Matrix<fun_scalar_t__,Eigen::Dynamic,Eigen::Dynamic>  Nug_mat(static_cast<Eigen::VectorXd::Index>(N),static_cast<Eigen::VectorXd::Index>(N));
        (void) Nug_mat;  // dummy to suppress unused var warning

        stan::math::initialize(Nug_mat, std::numeric_limits<double>::quiet_NaN());
        stan::math::fill(Nug_mat,DUMMY_VAR__);


        current_statement_begin__ = 5;
        stan::math::assign(parCov, rep_matrix(gamma,N,N));
        current_statement_begin__ = 6;
        stan::math::assign(Nug_mat, diag_matrix(nugget));
        current_statement_begin__ = 7;
        stan::math::assign(parCov, add(parCov,Nug_mat));
        current_statement_begin__ = 8;
        return stan::math::promote_scalar<fun_return_scalar_t__>(parCov);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct Cov_functor__ {
    template <typename T1__, typename T2__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T1__, T2__>::type, Eigen::Dynamic,Eigen::Dynamic>
    operator()(const int& N,
        const Eigen::Matrix<T1__, Eigen::Dynamic,1>& nugget,
        const T2__& gamma, std::ostream* pstream__) const {
        return Cov(N, nugget, gamma, pstream__);
    }
};

#include <meta_header.hpp>
 class model_oneK : public prob_grad {
private:
    int K;
    int N;
    int L;
    matrix_d obsCov;
    double varMeanFreqs;
    matrix_d LobsCov;
public:
    model_oneK(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }

    model_oneK(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }

    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning

        current_statement_begin__ = -1;

        static const char* function__ = "model_oneK_namespace::model_oneK";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        double DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        // initialize member variables
        try {
            current_statement_begin__ = 12;
            context__.validate_dims("data initialization", "K", "int", context__.to_vec());
            K = int(0);
            vals_i__ = context__.vals_i("K");
            pos__ = 0;
            K = vals_i__[pos__++];
            current_statement_begin__ = 13;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 14;
            context__.validate_dims("data initialization", "L", "int", context__.to_vec());
            L = int(0);
            vals_i__ = context__.vals_i("L");
            pos__ = 0;
            L = vals_i__[pos__++];
            current_statement_begin__ = 15;
            validate_non_negative_index("obsCov", "N", N);
            validate_non_negative_index("obsCov", "N", N);
            context__.validate_dims("data initialization", "obsCov", "matrix_d", context__.to_vec(N,N));
            validate_non_negative_index("obsCov", "N", N);
            validate_non_negative_index("obsCov", "N", N);
            obsCov = matrix_d(static_cast<Eigen::VectorXd::Index>(N),static_cast<Eigen::VectorXd::Index>(N));
            vals_r__ = context__.vals_r("obsCov");
            pos__ = 0;
            size_t obsCov_m_mat_lim__ = N;
            size_t obsCov_n_mat_lim__ = N;
            for (size_t n_mat__ = 0; n_mat__ < obsCov_n_mat_lim__; ++n_mat__) {
                for (size_t m_mat__ = 0; m_mat__ < obsCov_m_mat_lim__; ++m_mat__) {
                    obsCov(m_mat__,n_mat__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 16;
            context__.validate_dims("data initialization", "varMeanFreqs", "double", context__.to_vec());
            varMeanFreqs = double(0);
            vals_r__ = context__.vals_r("varMeanFreqs");
            pos__ = 0;
            varMeanFreqs = vals_r__[pos__++];

            // validate, data variables
            current_statement_begin__ = 12;
            check_greater_or_equal(function__,"K",K,1);
            current_statement_begin__ = 13;
            check_greater_or_equal(function__,"N",N,2);
            current_statement_begin__ = 14;
            check_greater_or_equal(function__,"L",L,(N + 1));
            current_statement_begin__ = 15;
            current_statement_begin__ = 16;
            // initialize data variables
            current_statement_begin__ = 19;
            validate_non_negative_index("LobsCov", "N", N);
            validate_non_negative_index("LobsCov", "N", N);
            LobsCov = matrix_d(static_cast<Eigen::VectorXd::Index>(N),static_cast<Eigen::VectorXd::Index>(N));
            stan::math::fill(LobsCov,DUMMY_VAR__);

            current_statement_begin__ = 20;
            stan::math::assign(LobsCov, multiply(L,obsCov));

            // validate transformed data
            current_statement_begin__ = 19;

            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 23;
            ++num_params_r__;
            current_statement_begin__ = 24;
            validate_non_negative_index("nugget", "N", N);
            num_params_r__ += N;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    ~model_oneK() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;

        if (!(context__.contains_r("gamma")))
            throw std::runtime_error("variable gamma missing");
        vals_r__ = context__.vals_r("gamma");
        pos__ = 0U;
        context__.validate_dims("initialization", "gamma", "double", context__.to_vec());
        double gamma(0);
        gamma = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,gamma);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable gamma: ") + e.what());
        }

        if (!(context__.contains_r("nugget")))
            throw std::runtime_error("variable nugget missing");
        vals_r__ = context__.vals_r("nugget");
        pos__ = 0U;
        validate_non_negative_index("nugget", "N", N);
        context__.validate_dims("initialization", "nugget", "vector_d", context__.to_vec(N));
        vector_d nugget(static_cast<Eigen::VectorXd::Index>(N));
        for (int j1__ = 0U; j1__ < N; ++j1__)
            nugget(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_lb_unconstrain(0,nugget);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable nugget: ") + e.what());
        }

        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

        T__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        try {
            // model parameters
            stan::io::reader<T__> in__(params_r__,params_i__);

            T__ gamma;
            (void) gamma;  // dummy to suppress unused var warning
            if (jacobian__)
                gamma = in__.scalar_lb_constrain(0,lp__);
            else
                gamma = in__.scalar_lb_constrain(0);

            Eigen::Matrix<T__,Eigen::Dynamic,1>  nugget;
            (void) nugget;  // dummy to suppress unused var warning
            if (jacobian__)
                nugget = in__.vector_lb_constrain(0,N,lp__);
            else
                nugget = in__.vector_lb_constrain(0,N);


            // transformed parameters
            current_statement_begin__ = 27;
            validate_non_negative_index("parCov", "N", N);
            validate_non_negative_index("parCov", "N", N);
            Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic>  parCov(static_cast<Eigen::VectorXd::Index>(N),static_cast<Eigen::VectorXd::Index>(N));
            (void) parCov;  // dummy to suppress unused var warning

            stan::math::initialize(parCov, DUMMY_VAR__);
            stan::math::fill(parCov,DUMMY_VAR__);


            current_statement_begin__ = 28;
            stan::math::assign(parCov, Cov(N,nugget,gamma, pstream__));

            // validate transformed parameters
            for (int i0__ = 0; i0__ < N; ++i0__) {
                for (int i1__ = 0; i1__ < N; ++i1__) {
                    if (stan::math::is_uninitialized(parCov(i0__,i1__))) {
                        std::stringstream msg__;
                        msg__ << "Undefined transformed parameter: parCov" << '[' << i0__ << ']' << '[' << i1__ << ']';
                        throw std::runtime_error(msg__.str());
                    }
                }
            }

            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 27;

            // model body

            current_statement_begin__ = 31;
            lp_accum__.add(normal_log<propto__>(nugget, 0, 1));
            current_statement_begin__ = 32;
            lp_accum__.add(normal_log<propto__>(gamma, varMeanFreqs, 0.5));
            current_statement_begin__ = 33;
            lp_accum__.add(wishart_log<propto__>(LobsCov, L, parCov));

        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }

        lp_accum__.add(lp__);
        return lp_accum__.sum();

    } // log_prob()

    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }


    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("gamma");
        names__.push_back("nugget");
        names__.push_back("parCov");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
        dims__.push_back(N);
        dimss__.push_back(dims__);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        vars__.resize(0);
        stan::io::reader<double> in__(params_r__,params_i__);
        static const char* function__ = "model_oneK_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double gamma = in__.scalar_lb_constrain(0);
        vector_d nugget = in__.vector_lb_constrain(0,N);
        vars__.push_back(gamma);
            for (int k_0__ = 0; k_0__ < N; ++k_0__) {
            vars__.push_back(nugget[k_0__]);
            }

        if (!include_tparams__) return;
        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        double DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {
            current_statement_begin__ = 27;
            validate_non_negative_index("parCov", "N", N);
            validate_non_negative_index("parCov", "N", N);
            matrix_d parCov(static_cast<Eigen::VectorXd::Index>(N),static_cast<Eigen::VectorXd::Index>(N));
            (void) parCov;  // dummy to suppress unused var warning

            stan::math::initialize(parCov, std::numeric_limits<double>::quiet_NaN());
            stan::math::fill(parCov,DUMMY_VAR__);


            current_statement_begin__ = 28;
            stan::math::assign(parCov, Cov(N,nugget,gamma, pstream__));

            // validate transformed parameters
            current_statement_begin__ = 27;

            // write transformed parameters
            for (int k_1__ = 0; k_1__ < N; ++k_1__) {
                for (int k_0__ = 0; k_0__ < N; ++k_0__) {
                vars__.push_back(parCov(k_0__, k_1__));
                }
            }

            if (!include_gqs__) return;
            // declare and define generated quantities



            // validate generated quantities

            // write generated quantities
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    static std::string model_name() {
        return "model_oneK";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "gamma";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= N; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "nugget" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;
        for (int k_1__ = 1; k_1__ <= N; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= N; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "parCov" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }

        if (!include_gqs__) return;
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "gamma";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= N; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "nugget" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;
        for (int k_1__ = 1; k_1__ <= N; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= N; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "parCov" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }

        if (!include_gqs__) return;
    }

}; // model

}

typedef model_oneK_namespace::model_oneK stan_model;


#endif
