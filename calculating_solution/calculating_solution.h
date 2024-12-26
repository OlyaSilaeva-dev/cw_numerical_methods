#ifndef CW_NUMERICAL_METHODS_CALCULATING_SOLUTION_H
#define CW_NUMERICAL_METHODS_CALCULATING_SOLUTION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <filesystem>

using namespace std;

double const epsilon = 1e-9;
double const iteration_limit = 1000;

struct GasProperties {
    double u_L, p_L, rho_L, gamma_L;
    double u_R, p_R, rho_R, gamma_R;
};

struct Results {
    double P;
    double U;
    double R_L, R_R;
    double D_L, D_R;
    double D_L_star, D_R_star;
};

void input_data(GasProperties &gas_props) {
    filesystem::current_path("C:/src/GitHub/cw_numerical_methods");
//    filesystem::current_path("/mnt/c/src/GitHub/cw_numerical_methods");
    string s = "input.txt";
    cout << "read file: " << s << endl;
    ifstream input(s);
    if (!input) {
        throw runtime_error("File not found");
    }

    input >> gas_props.u_L >> gas_props.p_L >> gas_props.rho_L >> gas_props.gamma_L;
    input >> gas_props.u_R >> gas_props.p_R >> gas_props.rho_R >> gas_props.gamma_R;
}

void write_output(const Results results, const string &filename) {
    ofstream out(filename);
    out << results.D_L << " " << results.D_L_star << " " << results.U << " " << results.D_R_star << " " << results.D_R
        << endl;
}

double calculate_c_k(double gamma_k, double p_k, double rho_k) {
    double tmp = gamma_k * p_k / rho_k;
    if (tmp < epsilon) {
        throw runtime_error("Speed of sound is too small or zero.");
    }
    double c_k = sqrt(tmp);
    return c_k;
}

double calculate_p_0_k(double gamma_k) {
    return (1.0 / gamma_k);
}

double calculate_pi_k(double P, double p_k) {
    return P / p_k;
}

double f(double P, double p_k, double rho_k, double gamma_k) {
    double c_k = calculate_c_k(gamma_k, p_k, rho_k);
    double pi_k = calculate_pi_k(P, p_k);

    double f;
    if (P >= p_k) {
        f = (P - p_k) /
               (rho_k * c_k * sqrt(((gamma_k + 1.0) / (2.0 * gamma_k) * pi_k + (gamma_k - 1.0) / (2.0 * gamma_k))));
    } else {
        f = (2.0 / (gamma_k - 1.0) * c_k * (pow(pi_k, ((gamma_k - 1.0) / (2.0 * gamma_k))) - 1.0));
    }
    return f;
}

double calculate_P_0(GasProperties g) {
    double c_L = calculate_c_k(g.gamma_L, g.p_L, g.rho_L);
    double c_R = calculate_c_k(g.gamma_R, g.p_R, g.rho_R);

//    cout << "\np_0_L: " << p_0_L << ", c_L: " << c_L << endl;
//    cout << "p_0_R: " << p_0_R << ", c_R: " << c_R << endl;

    return (g.p_L * g.rho_R * c_R + g.p_R * g.rho_L * c_L + (g.u_L - g.u_R) * g.rho_L * c_L * g.rho_R * c_R) /
           (g.rho_L * c_L + g.rho_R * c_R);
}

double d_f(double P, double p_k, double rho_k, double gamma_k) {
    double p_0_k = calculate_p_0_k(gamma_k);
    double c_k = calculate_c_k(gamma_k, p_k, rho_k);
    double pi_k = calculate_pi_k(P, p_k);

    if (P >= p_k) {
        return ((gamma_k + 1.0) * pi_k + (3.0 * gamma_k - 1.0)) /
               (4.0 * gamma_k * rho_k * c_k *
                pow(((gamma_k + 1.0) / (2.0 * gamma_k) * pi_k + (gamma_k - 1.0) / (2.0 * gamma_k)), 1.5));
    } else {
        return 1.0 / (gamma_k * (P + p_0_k)) * c_k * pow(pi_k, ((gamma_k - 1.0) / (2.0 * gamma_k)));
    }
}

double calculate_P(GasProperties g) {
    double P_0 = calculate_P_0(g);

    double P_prev = calculate_P_0(g);
//    if (P_prev < 0) {
//        P_prev =
//    }

    double P_cur = P_prev;

    double const_condition = epsilon * min(abs(g.p_L), abs(g.p_R));\

    int it = 0;
    do {
        P_prev = P_cur;
        P_cur = P_prev - (f(P_prev, g.p_L, g.rho_L, g.gamma_L) +
                          f(P_prev, g.p_R, g.rho_R, g.gamma_R) - (g.u_L - g.u_R)) /
                         (d_f(P_prev, g.p_L, g.rho_L, g.gamma_L) +
                          d_f(P_prev, g.p_R, g.rho_R, g.gamma_R));

        if (it > iteration_limit) {
            throw runtime_error("Iteration limit reached without convergence.");
        }
        it++;
    } while (abs(P_cur - P_prev) >= const_condition);

    return P_cur;
}

double calculate_U(double P, GasProperties g) {
    return 0.5 * (g.u_L - f(P, g.p_L, g.rho_L, g.gamma_L) + g.u_R +
                  f(P, g.p_R, g.rho_R, g.gamma_R));
}

Results calculate_other_results(double P, double U, GasProperties g) {
    Results results;
    results.P = P;
    results.U = U;

    double p_0_L = calculate_p_0_k(g.gamma_L);
    double p_0_R = calculate_p_0_k(g.gamma_R);

    double c_L = calculate_c_k(g.gamma_L, g.p_L, g.rho_L);
    double c_R = calculate_c_k(g.gamma_R, g.p_R, g.rho_R);

    if (-p_0_L < P && P < g.p_L) {
        double c_L_star = c_L + (g.gamma_L - 1.0) * (g.u_L - U) / 2.0;
        results.R_L = g.gamma_L * (P) / (c_L_star * c_L_star);
        results.D_L = g.u_L - c_L;
        results.D_L_star = g.u_L - c_L_star;

    } else if (P >= g.p_L) {
        results.R_L = g.rho_L * ((g.gamma_L + 1.0) * (P) + (g.gamma_L - 1.0) * (g.p_L)) /
                      ((g.gamma_L - 1.0) * P + (g.gamma_L + 1.0) * g.p_L);
        results.D_L = (g.p_L * g.u_L - results.R_L * U) / (g.rho_L - results.R_L);
        results.D_L_star = results.D_L;
    }

    if (-p_0_R < P && P < g.p_R) {
        double c_R_star = c_R - (g.gamma_R - 1.0) * (g.u_R - U) / 2.0;
        results.R_R = g.gamma_R * (P) / (c_R_star * c_R_star);
        results.D_R_star = U + c_R_star;
        results.D_R = g.u_R + c_R;

    } else if (P >= g.p_R) {
        results.R_R = g.rho_R * ((g.gamma_R + 1.0) * P + (g.gamma_R - 1.0) * (g.p_R)) /
                      ((g.gamma_R - 1.0) * (P) + (g.gamma_R + 1.0) * (g.p_R));
        results.D_R = (g.p_R * g.u_R - results.R_R * U) / (g.rho_R - results.R_R);
        results.D_R_star = results.D_R;
    }

    std::cout << "D_L: " << results.D_L << " D_L_star: " << results.D_L_star <<
    "\nD_R: " << results.D_R << " D_R_star: " << results.D_R_star <<
    "\nR_L: " << results.R_L << " R_R: " << results.R_R << endl;
    return results;
}

void calculate_parameters_less_D_l(double &rho_wave,double &u_wave,double &p_wave,double &gamma_wave, GasProperties g) {
    rho_wave = g.rho_L;
    u_wave = g.u_L;
    p_wave = g.p_L;
    gamma_wave = g.gamma_L;
}

///Функции расчета
void calculate_parameters_between_D_l_and_D_l_star(double &rho_wave,double &u_wave,double &p_wave,double &gamma_wave, GasProperties g,double i) {
//    double p_0_l = calculate_p_0_k(g.gamma_L);
    double c_L = calculate_c_k(g.gamma_L, g.p_L, g.rho_L);

    double c_wave = (2 * c_L + (g.gamma_L - 1.0) * (g.u_L - i)) / (g.gamma_L + 1.0);
    rho_wave = g.rho_L * pow((c_wave / c_L), 2.0 / (g.gamma_L - 1.0));
    u_wave = i + c_wave;
    p_wave = rho_wave * c_wave * c_wave / g.gamma_L;
    gamma_wave = g.gamma_L;
}

void calculate_parameters_between_D_l_star_and_U(double &rho_wave,double &u_wave,double &p_wave,double &gamma_wave, GasProperties g, Results results, double i) {
    rho_wave = results.R_L;
    u_wave = results.U;
    p_wave = results.P;
    gamma_wave = g.gamma_L;
}

void calculate_parameters_between_U_and_D_r_star(double &rho_wave,double &u_wave,double &p_wave,double &gamma_wave,double &rho_0_wave, double &c_0_wave, GasProperties g, Results results, double i) {
    rho_wave = results.R_R;
    u_wave = results.U;
    p_wave = results.P;
    gamma_wave = g.gamma_R;
}

void calculate_parameters_between_D_r_star_and_D_r(double &rho_wave,double &u_wave,double &p_wave,double &gamma_wave, GasProperties g, double i) {
    double c_R = calculate_c_k(g.gamma_R, g.p_R, g.rho_R);
    double c_wave = ((2.0 * c_R) + (g.gamma_R - 1.0) * (i - g.u_R)) / (g.gamma_R + 1);
    rho_wave = g.rho_R * pow((c_wave / c_R), (2.0 / (g.gamma_R - 1.0)));
    u_wave = i - c_wave;
    p_wave = rho_wave * c_wave * c_wave / g.gamma_R;
    gamma_wave = g.gamma_R;
}

void calculate_parameters_over_D_r(double &rho_wave,double &u_wave,double &p_wave,double &gamma_wave, GasProperties g) {
    rho_wave = g.rho_R;
    u_wave = g.u_R;
    p_wave = g.p_R;
    gamma_wave = g.gamma_R;
}

void func(GasProperties g, Results results, string filename) {
    ofstream out(filename);

    double rho_wave = 0.0, u_wave = 0.0, p_wave = 0.0, gamma_wave = 0.0;
    double rho_0_wave = 0.0, c_0_wave = 0.0; // Инициализируем переменные
    for (double i = -results.D_L - 2.0; i <= results.D_R + 2.0; i += 0.1) {
        if (i <= results.D_L) {
            calculate_parameters_less_D_l(rho_wave, u_wave, p_wave, gamma_wave, g);
        } else if (results.D_L < i && i <= results.D_L_star) {
            calculate_parameters_between_D_l_and_D_l_star(rho_wave, u_wave, p_wave, gamma_wave, g, i);
        } else if (results.D_L_star < i && i <= results.U) {
            calculate_parameters_between_D_l_star_and_U(rho_wave, u_wave, p_wave, gamma_wave, g, results, i);
        } else if (results.U < i && i <= results.D_R_star) {
            calculate_parameters_between_U_and_D_r_star(rho_wave, u_wave, p_wave, gamma_wave, rho_0_wave, c_0_wave, g, results, i);
        } else if (results.D_R_star < i && i <= results.D_R) { // Исправлено условие
            calculate_parameters_between_D_r_star_and_D_r(rho_wave, u_wave, p_wave, gamma_wave, g, i);
        } else {
            calculate_parameters_over_D_r(rho_wave, u_wave, p_wave, gamma_wave, g);
        }
        out << std::fixed << std::setprecision(3) // Округляем до 3 знаков
            << rho_wave << " " << u_wave << " " << p_wave << " " << gamma_wave
            << " " << rho_0_wave << " " << c_0_wave << endl;
    }
}

#endif //CW_NUMERICAL_METHODS_CALCULATING_SOLUTION_H
