#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <set>
#include <map>
#include <limits>

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

// ---------- Write CSV ----------
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

    for (int y = 1960; y <= 2022; ++y) {
        if (existingYears.find(y) == existingYears.end()) {
            vector<string> newRow(headers.size(), MISSING);
            newRow[0] = to_string(y);
            yearRowMap[y] = newRow;
        }
    }

    vector<vector<string>> fullData;
    fullData.push_back(headers);
    
    for (map<int, vector<string>>::const_iterator it = yearRowMap.begin(); it != yearRowMap.end(); ++it) {
        int year = it->first;
        const vector<string>& row = it->second;
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

        if (x.size() < (size_t)(degree + 1)) {
            cerr << "Not enough data to fit column: " << headers[col] << "\n";
            continue;
        }

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

        cout << "Fitted [" << headers[col] << "] with degree " << degree << endl;
    }

    writeCSV(outputFile, output);
    cout << "Regression fitting complete. " << outputFile << endl;
}

vector<double> fitLogistic(const vector<double>& x, const vector<double>& y) {
    double L = 100.0;  // Upper limit (100%)
    double k = 0.2;    // Growth
    double x0 = 2005;  // Tahun tengah

    // LSE
    double bestError = numeric_limits<double>::max();
    vector<double> bestParams = {L, k, x0};
    
    // Grid search 
    for (double kTry = 0.1; kTry <= 0.5; kTry += 0.05) {
        for (double x0Try = 2000; x0Try <= 2015; x0Try += 1.0) {
            double error = 0.0;
            for (size_t i = 0; i < x.size(); ++i) {
                double pred = L / (1.0 + exp(-kTry * (x[i] - x0Try)));
                error += pow(pred - y[i], 2);
            }
            
            if (error < bestError) {
                bestError = error;
                bestParams = {L, kTry, x0Try};
            }
        }
    }
    
    return bestParams;
}

double evaluateLogistic(const vector<double>& params, double x) {
    double L = params[0];  // Upper limit
    double k = params[1];  // Growth rate
    double x0 = params[2]; // Tahun tengah
    
    return L / (1.0 + exp(-k * (x - x0)));
}

void prediksiNilai(const string& inputFile, int degree = 2) 
{
    auto data = readCSV(inputFile);
    if (data.size() < 2) {
        cerr << "Insufficient data.\n";
        return;
    }

    vector<string> headers = data[0];
    
    vector<double> years, population, internetUsers;
    
    for (size_t i = 1; i < data.size(); ++i) {
        if (!isMissing(data[i][0])) {
            int year = stoi(data[i][0]);
            years.push_back(year);
            
            if (!isMissing(data[i][2])) {
                population.push_back(stod(data[i][2]));
            }
            
            if (!isMissing(data[i][1]) && stod(data[i][1]) > 0) {
                double internetVal = stod(data[i][1]);
                // H=Jika valid (internet tidak 0), masukkan ke array
                if (internetVal > 0.1) {
                    internetUsers.push_back(internetVal);
                }
            }
        }
    }
    
    // Filtered array
    vector<double> internetYears;
    for (size_t i = 1; i < data.size(); ++i) {
        if (!isMissing(data[i][0]) && !isMissing(data[i][1])) {
            double internetVal = stod(data[i][1]);
            if (internetVal > 0.1) {
                internetYears.push_back(stod(data[i][0]));
            }
        }
    }
    
    // Polynomial
    vector<double> popCoeffs = polyFit(years, population, degree);
    
    // fitLogistic untuk internet
    vector<double> internetParams = fitLogistic(internetYears, internetUsers);
    
    // Prediksi 
    double pop2030 = evaluatePoly(popCoeffs, 2030);
    double internet2035 = evaluateLogistic(internetParams, 2035);
    
    // Cek jika ada di range valid
    internet2035 = min(100.0, max(0.0, internet2035));
    
    // Output
    cout << "\n--- Preds ---\n";
    cout << "a. Tebakan Populasi Indonesia di 2030: " 
         << fixed << setprecision(0) << pop2030 << " people\n";
    cout << "b. Tebakan Persentase Pengguna Internet Indonesia di 2035: " 
         << fixed << setprecision(2) << internet2035 << "%\n";

    cout << "\n";
}

// ---------- Entry ----------
int main() {
    string inputFile = "data.csv";
    string outputFile = "regression_output.csv";
    int degree = 3; // Agar lebih akurat
    // Jika rendah, bobot tahun dengan internet 0 akan jadi terlalu besar, sehingga tidak akurat
    // Jika terlalu tinggi, overfitting
    fitCSV(inputFile, outputFile, degree);
    prediksiNilai(outputFile, degree);
    return 0;
}