#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <set>
#include <map>
#include <algorithm>

using namespace std;

const string MISSING = "NaN";

bool isMissing(const string& val) {
    return val.empty() || val == MISSING;
}

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

void writeCSV(const string& filename, const vector<vector<string>>& data) {
    ofstream outFile(filename);
    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            outFile << row[i];
            if (i < row.size() - 1) outFile << ",";
        }
        outFile << "\n";
    }
}

vector<double> polyFit(const vector<double>& x, const vector<double>& y, int degree) {
    int n = degree + 1;
    vector<vector<double>> X(n, vector<double>(n, 0.0));
    vector<double> Y(n, 0.0);

    for (size_t i = 0; i < x.size(); ++i) {
        double xi = 1.0;
        for (int j = 0; j < n; ++j) {
            double xij = xi;
            for (int k = 0; k < n; ++k)
                X[j][k] += xij, xij *= x[i];
            Y[j] += xi * y[i];
            xi *= x[i];
        }
    }

    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k)
            if (fabs(X[k][i]) > fabs(X[maxRow][i])) maxRow = k;
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

double evaluatePoly(const vector<double>& coeffs, double x) {
    double result = 0.0, xp = 1.0;
    for (double c : coeffs) result += c * xp, xp *= x;
    return result;
}

// ---------- Logistic Fit for Internet Usage ----------
struct LogisticParams {
    double L = 100.0; // Upper limit
    double k = 0.2;   // Growth rate
    double x0 = 2005; // Midpoint year
};

double logistic(double x, const LogisticParams& p) {
    return p.L / (1.0 + exp(-p.k * (x - p.x0)));
}

// Very basic logistic parameter fitting using brute-force search
LogisticParams fitLogistic(const vector<double>& x, const vector<double>& y) {
    LogisticParams best;
    double bestError = 1e12;
    for (double L = 80; L <= 120; L += 2) {
        for (double k = 0.01; k <= 1.0; k += 0.01) {
            for (double x0 = 1990; x0 <= 2020; x0 += 1.0) {
                double error = 0.0;
                for (size_t i = 0; i < x.size(); ++i) {
                    double pred = L / (1.0 + exp(-k * (x[i] - x0)));
                    error += pow(pred - y[i], 2);
                }
                if (error < bestError) {
                    bestError = error;
                    best.L = L; best.k = k; best.x0 = x0;
                }
            }
        }
    }
    return best;
}

void fitCSV(const string& inputFile, const string& outputFile, int degree = 2) {
    auto data = readCSV(inputFile);
    if (data.size() < 2) return;

    vector<string> headers = data[0];
    map<int, vector<string>> yearRowMap;
    set<int> existingYears;

    for (size_t i = 1; i < data.size(); ++i) {
        if (!isMissing(data[i][0])) {
            int year = stoi(data[i][0]);
            yearRowMap[year] = data[i];
            existingYears.insert(year);
        }
    }

    for (int y = 1960; y <= 2023; ++y) {
        if (existingYears.find(y) == existingYears.end()) {
            vector<string> newRow(headers.size(), MISSING);
            newRow[0] = to_string(y);
            yearRowMap[y] = newRow;
        }
    }

    vector<vector<string>> fullData;
    fullData.push_back(headers);
    for (const auto& [year, row] : yearRowMap) {
        fullData.push_back(row);
    }

    int rows = fullData.size();
    int cols = headers.size();
    vector<vector<string>> output = fullData;

    for (int col = 1; col < cols; ++col) {
        vector<double> x, y;
        for (int i = 1; i < rows; ++i) {
            if (!isMissing(fullData[i][0]) && !isMissing(fullData[i][col])) {
                x.push_back(stod(fullData[i][0]));
                y.push_back(stod(fullData[i][col]));
            }
        }

        if (x.size() < (size_t)(degree + 1)) continue;

        if (headers[col] == "Percentage_Internet_User") {
            LogisticParams params = fitLogistic(x, y);
            for (int i = 1; i < rows; ++i) {
                if (!isMissing(fullData[i][0]) && isMissing(fullData[i][col])) {
                    double xi = stod(fullData[i][0]);
                    double yi = logistic(xi, params);
                    ostringstream oss;
                    oss << fixed << setprecision(4) << yi;
                    output[i][col] = oss.str();
                }
            }
            cout << "Fitted [" << headers[col] << "] using Logistic Regression\n";
        } else {
            auto coeffs = polyFit(x, y, degree);
            for (int i = 1; i < rows; ++i) {
                if (!isMissing(fullData[i][0]) && isMissing(fullData[i][col])) {
                    double xi = stod(fullData[i][0]);
                    double yi = evaluatePoly(coeffs, xi);
                    ostringstream oss;
                    oss << fixed << setprecision(4) << yi;
                    output[i][col] = oss.str();
                }
            }
            cout << "Fitted [" << headers[col] << "] using Polynomial Regression\n";
        }
    }

    writeCSV(outputFile, output);
    cout << "Complete. Output written to: " << outputFile << endl;
}

int main() {
    string inputFile = "data.csv";
    string outputFile = "regression_output.csv";
    int degree = 2;
    fitCSV(inputFile, outputFile, degree);
    return 0;
}

coba pake ini