#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <limits>

using namespace std;

const string MISSING_MARKER = "NaN";

// ---------- Utility Functions ----------

bool isMissing(const string& value) {
    return value.empty() || value == MISSING_MARKER;
}

// ---------- CSV Reader ----------

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

// ---------- Cubic Spline Interpolation ----------

vector<double> cubicSplineInterpolate(
    const vector<double>& xs,
    const vector<double>& ys,
    const vector<double>& queryXs
) {
    int n = xs.size() - 1;
    vector<double> a(ys);
    vector<double> b(n), d(n), h(n), alpha(n), c(n + 1), l(n + 1), mu(n + 1), z(n + 1);

    for (int i = 0; i < n; ++i)
        h[i] = xs[i + 1] - xs[i];

    for (int i = 1; i < n; ++i)
        alpha[i] = (3.0 / h[i]) * (a[i + 1] - a[i]) - (3.0 / h[i - 1]) * (a[i] - a[i - 1]);

    l[0] = 1.0; mu[0] = 0.0; z[0] = 0.0;

    for (int i = 1; i < n; ++i) {
        l[i] = 2.0 * (xs[i + 1] - xs[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n] = 1.0; z[n] = 0.0; c[n] = 0.0;

    for (int j = n - 1; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
    }

    vector<double> interpolated;
    for (double xq : queryXs) {
        int i = 0;
        while (i < n && xq > xs[i + 1]) ++i;
        double dx = xq - xs[i];
        double result = a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
        interpolated.push_back(result);
    }

    return interpolated;
}

// ---------- Main Interpolation ----------

void interpolateCSV(const string& inputFile, const string& outputFile) {
    auto data = readCSV(inputFile);
    if (data.empty()) {
        cerr << "Empty CSV file.\n";
        return;
    }

    int numCols = data[0].size();
    int numRows = data.size();
    vector<string> headers = data[0];

    vector<vector<string>> outputData = data; 
    vector<double> xValues;

    for (int row = 1; row < numRows; ++row) {
        xValues.push_back(stod(data[row][0]));
    }

    // Process internet percentage column (column 1) specially
    if (numCols > 1) {
        int col = 1; // Internet percentage column
        vector<double> knownX, knownY, queryX;
        vector<int> missingRows;

        // For internet percentage, only consider non-zero values for interpolation
        for (int row = 1; row < numRows; ++row) {
            string cell = data[row][col];
            if (isMissing(cell)) {
                queryX.push_back(xValues[row - 1]);
                missingRows.push_back(row);
            } else {
                double value = stod(cell);
                // Only use non-zero values as reference points
                if (value > 0) {
                    knownX.push_back(xValues[row - 1]);
                    knownY.push_back(value);
                } else if (row > 1) {
                    // For years with 0% internet, keep them as 0 if before first non-zero year
                    if (outputData[row][col] == "0" && knownX.empty()) {
                        continue; // Keep as 0
                    } else {
                        // Otherwise treat as missing data to interpolate
                        queryX.push_back(xValues[row - 1]);
                        missingRows.push_back(row);
                    }
                }
            }
        }

        if (knownX.size() < 2) {
            cerr << "Not enough non-zero data to interpolate internet percentage column\n";
        } else {
            vector<double> interpolatedY = cubicSplineInterpolate(knownX, knownY, queryX);

            for (size_t i = 0; i < missingRows.size(); ++i) {
                int row = missingRows[i];
                double year = xValues[row - 1];
                
                // If the year is before the first non-zero internet usage year,
                // set it to 0 rather than interpolating
                if (year < knownX[0]) {
                    outputData[row][col] = "0";
                } else {
                    // Otherwise use the interpolated value
                    double value = interpolatedY[i];
                    // Ensure interpolated values are non-negative
                    value = max(0.0, value);
                    
                    ostringstream oss;
                    oss << fixed << setprecision(4) << value;
                    outputData[row][col] = oss.str();
                }
            }
        }
    }

    // Process all other columns (like population) with normal interpolation
    for (int col = 2; col < numCols; ++col) {
        vector<double> knownX, knownY, queryX;
        vector<int> missingRows;

        for (int row = 1; row < numRows; ++row) {
            string cell = data[row][col];
            if (isMissing(cell)) {
                queryX.push_back(xValues[row - 1]);
                missingRows.push_back(row);
            } else {
                knownX.push_back(xValues[row - 1]);
                knownY.push_back(stod(cell));
            }
        }

        if (knownX.size() < 2) {
            cerr << "Not enough data to interpolate column: " << headers[col] << endl;
            continue;
        }

        vector<double> interpolatedY = cubicSplineInterpolate(knownX, knownY, queryX);

        for (size_t i = 0; i < missingRows.size(); ++i) {
            int row = missingRows[i];
            ostringstream oss;
            oss << fixed << setprecision(4) << interpolatedY[i];
            outputData[row][col] = oss.str();
        }
    }

    // Write output
    ofstream outFile(outputFile);
    for (const auto& row : outputData) {
        for (size_t i = 0; i < row.size(); ++i) {
            outFile << row[i];
            if (i < row.size() - 1) outFile << ",";
        }
        outFile << "\n";
    }

    cout << "Interpolation complete. Output written to: " << outputFile << endl;
}

// ---------- Entry ----------

int main() {
    string inputFile = "data.csv";
    string outputFile = "interpolated_output.csv";

    interpolateCSV(inputFile, outputFile);

    return 0;
}