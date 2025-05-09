#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>

using namespace std;

// Dengan gnuplot, lebih mudah dijalankan di Linux
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

void generatePlotScript(const string& csvFile, const string& scriptFile, const string& outputFile, int numCols) {
    ofstream script(scriptFile);

    script << "set datafile separator ','\n";
    script << "set terminal png size 1000,600\n";
    script << "set output '" << outputFile << "'\n";
    script << "set title 'CSV Plot'\n";
    script << "set key outside\n";
    script << "set xlabel 'X'\n";
    script << "set ylabel 'Y'\n";
    script << "plot ";

    for (int col = 2; col <= numCols; ++col) {
        script << "'" << csvFile << "' using 1:" << col << " with lines title columnheader(" << col << ")";
        if (col < numCols) script << ", ";
    }
    script << "\n";
    script.close();
}

void plotCSV(const string& csvFile, const string& imageOutput = "plot.png") {
    auto data = readCSV(csvFile);
    if (data.empty() || data[0].size() < 2) {
        cerr << "CSV must have at least two columns.\n";
        return;
    }

    int numCols = data[0].size();

    string scriptFile = "plot_script.gp";
    generatePlotScript(csvFile, scriptFile, imageOutput, numCols);

    string command = "gnuplot " + scriptFile;
    system(command.c_str());

    cout << "Plot generated: " << imageOutput << endl;
}

int main() {
    string csvFile = "regression_output.csv";  
    plotCSV(csvFile);
    return 0;
}
