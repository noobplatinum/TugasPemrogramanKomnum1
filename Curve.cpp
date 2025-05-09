#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>

using namespace std;

const string MISSING = "NaN";

// ---------- Utility ----------
bool isMissing(const string& val) {
    return val.empty() || val == MISSING;
}

// ---------- Read CSV ----------
vector<vector<string>> readCSV(const string& filename) {
    ifstream file(filename);
    string line;
    vector<vector<string>> data;

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

// ---------- Polynomial Regression ----------
vector<double> polyFit(const vector<double>& x, const vector<double>& y, int degree) {
    int n = degree + 1;
    vector<vector<double>> X(n, vector<double>(n, 0.0));
    vector<double> Y(n, 0.0);

    for (size_t i = 0; i < x.size(); ++i) {
        double xi = 1.0;
        for (int j = 0; j < n; ++j) {
            double xij = xi;
            for (int k = 0; k < n; ++k) {
                X[j][k] += xij;
                xij *= x[i];
            }
            Y[j] += xi * y[i];
            xi *= x[i];
        }
    }

    // Gaussian Elimination
    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k)
            if (fabs(X[k][i]) > fabs(X[maxRow][i]))
                maxRow = k;

        swap(X[i], X[maxRow]);
        swap(Y[i], Y[maxRow]);

        for (int k = i + 1; k < n; ++k) {
            double factor = X[k][i] / X[i][i];
            for (int j = i; j < n; ++j)
                X[k][j] -= factor * X[i][j];
            Y[k] -= factor * Y[i];
        }
    }

    vector<double> coeffs(n);
    for (int i = n - 1; i >= 0; --i) {
        coeffs[i] = Y[i];
        for (int j = i + 1; j < n; ++j)
            coeffs[i] -= X[i][j] * coeffs[j];
        coeffs[i] /= X[i][i];
    }

    return coeffs;
}

// ---------- Evaluasi Polynomial ----------
double evaluatePoly(const vector<double>& coeffs, double x) {
    double result = 0.0;
    double xp = 1.0;
    for (double c : coeffs) {
        result += c * xp;
        xp *= x;
    }
    return result;
}

// ---------- Main Fit Function ----------
void fitCSV(const string& inputFile, const string& outputFile, int degree = 2) {
    auto data = readCSV(inputFile);
    if (data.size() < 2) {
        cerr << "Insufficient data.\n";
        return;
    }

    int rows = data.size();
    int cols = data[0].size();
    vector<string> headers = data[0];

    vector<double> x;
    for (int i = 1; i < rows; ++i) {
        if (!isMissing(data[i][0]))
            x.push_back(stod(data[i][0]));
        else
            cerr << "Warning: Missing X value at row " << i << "\n";
    }

    vector<vector<string>> output = data;

    for (int col = 1; col < cols; ++col) {
        vector<double> y, validX;
        vector<int> rowIndex;

        for (int i = 1; i < rows; ++i) {
            if (!isMissing(data[i][0]) && !isMissing(data[i][col])) {
                validX.push_back(stod(data[i][0]));
                y.push_back(stod(data[i][col]));
                rowIndex.push_back(i);
            }
        }

        if (validX.size() < (size_t)(degree + 1)) {
            cerr << "Not enough data to fit column: " << headers[col] << "\n";
            continue;
        }

        vector<double> coeffs = polyFit(validX, y, degree);

        for (int i = 1; i < rows; ++i) {
            if (!isMissing(data[i][0])) {
                double xi = stod(data[i][0]);
                double yi = evaluatePoly(coeffs, xi);
                ostringstream oss;
                oss << fixed << setprecision(4) << yi;
                output[i][col] = oss.str();
            }
        }

        cout << "Fitted [" << headers[col] << "] with polynomial degree " << degree << endl;
    }

    ofstream outFile(outputFile);
    for (const auto& row : output) {
        for (size_t i = 0; i < row.size(); ++i) {
            outFile << row[i];
            if (i < row.size() - 1) outFile << ",";
        }
        outFile << "\n";
    }

    cout << "Regression fitting complete. Output written to: " << outputFile << endl;
}

// ---------- Entry ----------
int main() {
    string inputFile = "data.csv";
    string outputFile = "regression_output.csv";
    int degree = 2;

    fitCSV(inputFile, outputFile, degree);
    return 0;
}
