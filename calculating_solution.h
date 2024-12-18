#ifndef CW_NUMERICAL_METHODS_CALCULATING_SOLUTION_H
#define CW_NUMERICAL_METHODS_CALCULATING_SOLUTION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

double const epsilon = 1e-9;

struct GasProperties {
    double rho_L, u_L, p_L, c_L, gamma_L, rho_0_L, c_0_L, p_0_L;
    double rho_R, u_R, p_R, c_R, gamma_R, rho_0_R, c_0_R, p_0_R;
};

struct Results {
    double P;
    double U;
    double R_L, R_R;
    double D_L, D_R;
    double D_L_star = 0, D_R_star = 0;
};

struct State {
    double rho; //плотность
    double u; //скорость
    double p; //давление
};

void input_data(GasProperties &gas_props) {
    ifstream input("input.txt");
    input >> gas_props.rho_L >> gas_props.u_L >> gas_props.p_L >> gas_props.c_L >> gas_props.gamma_L >> gas_props.rho_0_L >> gas_props.c_0_L >> gas_props.p_0_L;
    input >> gas_props.rho_R >> gas_props.u_R >> gas_props.p_R >> gas_props.c_R >> gas_props.gamma_R >> gas_props.rho_0_R >> gas_props.c_0_R >> gas_props.p_0_R;
}

void write_output(const vector<State> &results, const string &filename) {
    ofstream out(filename);
    for (const auto &state: results) {
        out << state.rho << " " << state.u << " " << state.p << endl;
    }
}

double calculate_c_k(double gamma_k, double p_k, double p_0_k, double rho_k) {
    return sqrt(gamma_k * (p_k - p_0_k) / rho_k);
}

double calculate_p_0_k(double gamma_k, double rho_0_k, double c_0_k) {
    return 1 / gamma_k * rho_0_k * c_0_k * c_0_k;
}

double calculate_pi_k(double P, double p_0_k, double p_k) {
    return (P + p_0_k) / (p_k + p_0_k);
}

double f(double P, double p_k, double rho_k, double gamma_k, double rho_0_k, double c_0_k) {
    double p_0_k = calculate_p_0_k(gamma_k, rho_0_k, c_0_k);
    double c_k = calculate_c_k(gamma_k, p_k, p_0_k, rho_k);
    double pi_k = calculate_pi_k(P, p_0_k, p_k);

    if (P >= p_k) {
        return (P - p_k) /
               (rho_k * c_k * sqrt(((gamma_k + 1.0) / (2.0 * gamma_k) * pi_k + (gamma_k - 1.0) / (2.0 * gamma_k))));
    } else {
        return (2.0 / (gamma_k - 1) * c_k * (pow((pi_k), (gamma_k - 1) / (2 * gamma_k)) - 1));
    }
}

double calculate_P_0(GasProperties g) {
    return (g.p_L * g.rho_R * g.c_R + g.p_R * g.rho_L * g.c_L + (g.u_L - g.u_R) * g.rho_L * g.c_L * g.rho_R * g.c_R) /
           (g.rho_L * g.c_L + g.rho_R * g.c_R);
}

double d_f(double P, double p_k, double rho_k, double gamma_k, double rho_0_k, double c_0_k) {
    double p_0_k = calculate_p_0_k(gamma_k, rho_0_k, c_0_k);
    double c_k = calculate_c_k(gamma_k, p_k, p_0_k, rho_k);
    double pi_k = calculate_pi_k(P, p_0_k, p_k);

    if (P >= p_k) {
        return ((gamma_k + 1) * pi_k + (3.0 * gamma_k - 1.0)) /
               (4.0 * gamma_k * rho_k * c_k *
                pow(((gamma_k + 1) / (2.0 * gamma_k) * pi_k + (gamma_k - 1) / (2.0 * gamma_k)), 1.5));
    } else {
        return 1.0 / (gamma_k * (P + p_k)) * c_k * pow(pi_k, (gamma_k - 1) / (2.0 * gamma_k));
    }
}

double calculate_P(GasProperties g) {

    vector<double> P_vec;
    double P_0 = calculate_P_0(g);
    P_vec.push_back(P_0);

    double P_prev = P_vec.back();
    double P_cur = P_vec.back();
    double const_condition = epsilon * min(abs(g.p_L), abs(g.p_R));

    do {
        P_cur = P_prev - (f(P_prev, g.p_L, g.rho_L, g.gamma_L, g.rho_0_L, g.c_0_L) +
                          f(P_prev, g.p_R, g.rho_R, g.gamma_R, g.rho_0_R, g.c_0_R) - (g.u_L = g.u_R)) /
                         d_f(P_prev, g.p_L, g.rho_L, g.gamma_L, g.rho_0_L, g.c_0_L) +
                d_f(P_prev, g.p_R, g.rho_R, g.gamma_R, g.rho_0_R, g.c_0_R);
        P_prev = P_vec.back();
        P_vec.push_back(P_cur);
    } while (abs(P_cur - P_prev) >= const_condition);

    return P_vec.back();
}

double calculate_U(double P, GasProperties g) {
    return 0.5 * (g.u_L - f(P, g.p_L, g.rho_L, g.gamma_L, g.rho_0_L, g.c_0_L) + g.u_R +
                  f(P, g.p_R, g.rho_R, g.gamma_R, g.rho_0_R, g.c_0_R));
}

Results calculate_other_results(double P, double U, GasProperties g) {
    Results results;
    results.P = P;
    results.U = U;

    if (g.p_0_L < P < g.p_L) {
        double c_L_star = g.c_L + (g.gamma_L - 1)*(g.u_L - U) / 2.0;
        results.R_L = g.gamma_L * (P + g.p_0_L) / (c_L_star * c_L_star);
        results.D_L = g.u_L - g.c_L;
        results.D_L_star = g.u_L - c_L_star;

    } else if (P >= g.p_L) {
        
    }

    if (-g.p_0_R < P < g.p_R) {

    } else if (P >= g.p_R) {

    }
}

#endif //CW_NUMERICAL_METHODS_CALCULATING_SOLUTION_H
