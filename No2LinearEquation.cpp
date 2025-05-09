#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

using namespace std;

// Fungsi membaca file CSV
vector<vector<string>> readCSV(const string& filename) {
    ifstream file(filename);
    vector<vector<string>> data;
    string line;

    while (getline(file, line)) {
        stringstream ss(line);
        string cell;
        vector<string> row;
        while (getline(ss, cell, ',')) {
            row.push_back(cell);
        }
        data.push_back(row);
    }

    return data;
}

// Regresi Polinomial Orde 2 (Kuadratik)
vector<double> polynomial_regression_order2(const vector<double>& x, const vector<double>& y) {
    int degree = 2;
    vector<vector<double>> A(degree + 1, vector<double>(degree + 2, 0.0));

    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < x.size(); ++k) {
                sum += pow(x[k], i + j);
            }
            A[i][j] = sum;
        }

        double sumY = 0.0;
        for (size_t k = 0; k < x.size(); ++k) {
            sumY += y[k] * pow(x[k], i);
        }
        A[i][degree + 1] = sumY;
    }

    // Eliminasi Gauss
    for (int i = 0; i <= degree; ++i) {
        double factor = A[i][i];
        for (int j = i; j <= degree + 1; ++j)
            A[i][j] /= factor;

        for (int k = 0; k <= degree; ++k) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = i; j <= degree + 1; ++j)
                    A[k][j] -= factor * A[i][j];
            }
        }
    }

    vector<double> coeffs(degree + 1);
    for (int i = 0; i <= degree; ++i)
        coeffs[i] = A[i][degree + 1];

    return coeffs;
}

// Evaluasi polinomial
double evaluatePoly(const vector<double>& coeffs, double x) {
    double result = 0.0;
    double xp = 1.0;
    for (double c : coeffs) {
        result += c * xp;
        xp *= x;
    }
    return result;
}

// Hitung R²
double calculateR2(const vector<double>& y_true, const vector<double>& y_pred) {
    double mean_true = 0.0;
    for (double y : y_true) mean_true += y;
    mean_true /= y_true.size();

    double ss_res = 0.0, ss_tot = 0.0;
    for (size_t i = 0; i < y_true.size(); ++i) {
        double err = y_true[i] - y_pred[i];
        ss_res += err * err;
        double dev = y_true[i] - mean_true;
        ss_tot += dev * dev;
    }

    return 1.0 - ss_res / ss_tot;
}

// Entry Point
int main() {
    string filename = "regression_output.csv"; // Pastikan file ini tersedia
    auto data = readCSV(filename);

    vector<double> years, population, internet;

    for (size_t i = 1; i < data.size(); ++i) {
        years.push_back(stod(data[i][0]));
        internet.push_back(stod(data[i][1]));
        population.push_back(stod(data[i][2]));
    }

    cout << fixed << setprecision(6);

    // Normalisasi input
    double base_year = 2000;
    double scaling_factor = 10;         // Setiap unit x = 10 tahun
    double population_scale = 1e6;      // Populasi dalam satuan juta jiwa

    vector<double> years_normalized;
    vector<double> population_scaled;
    for (double year : years) {
        years_normalized.push_back((year - base_year) / scaling_factor);
    }

    for (double pop : population) {
        population_scaled.push_back(pop / population_scale);
    }

    // Regresi Orde 2 untuk % pengguna internet
    vector<double> coeffs_internet = polynomial_regression_order2(years_normalized, internet);
    cout << "Persamaan Persentase Pengguna Internet:\n";
    cout << "y = ";
    cout << coeffs_internet[0] << "*x^2 + ";
    cout << coeffs_internet[1] << "*x + ";
    cout << coeffs_internet[2] << "\n";
    cout << "(x = (tahun - " << base_year << ") / " << scaling_factor << ")\n\n";

    // Regresi Orde 2 untuk populasi
    vector<double> coeffs_population = polynomial_regression_order2(years_normalized, population_scaled);
    cout << "Persamaan Populasi Indonesia (juta jiwa):\n";
    cout << "y = ";
    cout << coeffs_population[0] << "*x^2 + ";
    cout << coeffs_population[1] << "*x + ";
    cout << coeffs_population[2] << "\n";
    cout << "(x = (tahun - " << base_year << ") / " << scaling_factor << ")\n\n";

    // Prediksi ke depan
    double x_2030 = (2030 - base_year) / scaling_factor;
    double pop_2030_million = evaluatePoly(coeffs_population, x_2030);
    double pop_2030_total = pop_2030_million * population_scale;

    cout << "Prediksi populasi Indonesia di tahun 2030: " << pop_2030_total << endl;

    double x_2035 = (2035 - base_year) / scaling_factor;
    double internet_2035 = evaluatePoly(coeffs_internet, x_2035);
    
    // Batasi persentase maksimal 100%
    if (internet_2035 > 100) internet_2035 = 100;
    if (internet_2035 < 0) internet_2035 = 0;

    cout << "Prediksi % pengguna internet di tahun 2035: " << internet_2035 << "%\n\n";

    // Hitung R²
    vector<double> pred_internet, pred_population_million;
    for (size_t i = 0; i < years.size(); ++i) {
        pred_internet.push_back(evaluatePoly(coeffs_internet, years_normalized[i]));
        pred_population_million.push_back(evaluatePoly(coeffs_population, years_normalized[i]));
    }

    double r2_internet = calculateR2(internet, pred_internet);
    double r2_population = calculateR2(population_scaled, pred_population_million);

    cout << "R² untuk % pengguna internet: " << r2_internet << endl;
    cout << "R² untuk populasi: " << r2_population << endl;

    return 0;
}