#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include<algorithm>

using namespace std;

#define N_SECONDS 10


vector<double>* fill_sin(double F, double quantization_min, double quantization_max, double phi, double T) {
    vector<double>* res = new vector<double>;
    uint8_t quantization_levels_num = 64;
    int N_SAMPLES = N_SECONDS * F + 1;
    uint8_t* digital_signal = new uint8_t[N_SAMPLES];
    double k = M_PI * 2 / T; //(M_PI*2/T) - коэффициент k в синусе
    for (int i = 0; i < N_SAMPLES; i++) digital_signal[i] = (round((sin(i * k / F + phi) + 1) * (quantization_levels_num - 1))) /
        (quantization_max - quantization_min);
    for (size_t i = 0; i < N_SAMPLES; ++i)
    {
        double signal_val = double(digital_signal[i]) * (quantization_max - quantization_min) /
            (quantization_levels_num - 1) + quantization_min;
        // time = (i / F)
        res->push_back(signal_val);
    }
    delete[] digital_signal;
    return res;
}

void interpolation0(vector<double>* signal0, double new_F, double F) {
    int N_SAMPLES_NEW = N_SECONDS * new_F + 1;
    const char csv_file_name[64] = "data.csv";
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal interpolated\n";
    for (int i = 0; i < N_SAMPLES_NEW; i++) {
        csv_file << (i / new_F) << "," << signal0->at(floor(i * F / new_F)) << "\n";
    }
    csv_file.close();
    return;
}

void interpolation1(vector<double>* signal0, double new_F, double F, double quantization_min, double quantization_max) {
    uint8_t quantization_levels_num = 64;
    int N_SAMPLES_NEW = N_SECONDS * new_F + 1;
    const char csv_file_name[64] = "data.csv";
    std::ofstream csv_file;
    csv_file.open(csv_file_name);
    csv_file << "time,signal interpolated,signal interpolated + quantized\n";
    for (int i = 0; i < N_SAMPLES_NEW - 1; i++) {
        int i1 = floor(i * F / new_F);
        int i2 = i1 + 1;
        double x1 = i1 / F;
        double x2 = i2 / F;
        double y1 = signal0->at(i1);
        double y2 = signal0->at(i2);
        double k = (y2 - y1) / (x2 - x1);
        double b = y1 - k * x1;
        double signal_val = k * (i / new_F) + b; //и надо проквантовать, т.к. получились значения между уровнями квантования
        uint8_t signal_quant = (round((signal_val - quantization_min) * (quantization_levels_num - 1))) /
            (quantization_max - quantization_min);
        double res = double(signal_quant) * (quantization_max - quantization_min) /
            (quantization_levels_num - 1) + quantization_min;
        csv_file << (i / new_F) << "," << signal_val << "," << res << "\n";
    }
    csv_file.close();
    return;
}

int main() {

    double F = 10.; // частота дискретизации
    double quantization_min = -1., quantization_max = 1.;
    double phi = 2 * M_PI;
    double T = 2 * M_PI;
    vector<double>* signal = fill_sin(F, quantization_min, quantization_max, phi, T);
    interpolation0(signal, 20., F);
    delete signal;
    return 0;
}
