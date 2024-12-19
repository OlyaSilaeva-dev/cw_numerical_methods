#include "calculating_solution.h"
#include "iomanip"

int main() {
    GasProperties gasProperties;
    input_data(gasProperties);

    cout << "gamma_L: " << setprecision(9) << fixed << gasProperties.gamma_L << " gamma_R: " << gasProperties.gamma_R << std::endl;
    cout << "u_L:     " << setprecision(9) << fixed << gasProperties.u_L << " u_R:     " << gasProperties.u_R << std::endl;
    cout << "p_L:     " << setprecision(9) << fixed << gasProperties.p_L << " p_R:   " << gasProperties.p_R << std::endl;
    cout << "rho_L:   " << setprecision(9) << fixed <<gasProperties.rho_L << " rho_R:   " << gasProperties.rho_R << std::endl;


    double P = calculate_P(gasProperties);
    cout << "\nP: " << P << endl;
    double U = calculate_U(P, gasProperties);
    cout << "U: " << U << endl;

    Results results = calculate_other_results(P, U, gasProperties);
    write_output(results, "output.txt");
    func(gasProperties, results, "output_values.txt");

    return 0;
}
