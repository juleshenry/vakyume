#include <cmath>
#include <vector>
#include <complex>

#include <iostream>

std::vector<std::complex<double>> calculate_T_e(
    double P, double T_i, double S_0, double S_p,
    double p_0, double p_c, double p_s) {
    std::vector<std::complex<double>> result = {};
    std::complex<double> T_e = 5.5;
    result.push_back(T_e);
    T_e = 6.5;
    result.push_back(T_e);
    return result;
}



std::vector<double> a() {
    std::vector<double> myVector = {1.23, 4.56, 7.89};
    return myVector;
}

int main() {
    double dummy_val = 1.11;
    std::vector<std::complex<double>> myVector = calculate_T_e(
        dummy_val,dummy_val,dummy_val,dummy_val,dummy_val,dummy_val,dummy_val);
    for (int i = 0; i < myVector.size(); i++) {
        std::cout << myVector[i] << " ";
    }

    std::cout << std::endl;

    return 0;
}
