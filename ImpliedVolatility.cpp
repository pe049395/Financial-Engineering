#include <iostream>
#include <cmath>

class Option {
private:
    double S; // 주식 가격
    double K; // 행사 가격
    double r; // 무위험 이자율
    double T; // 잔여 만기
    double optionValue; // 옵션 가격
    char optionType; // 옵션 종류

public:
    Option(double stockPrice, double strikePrice, double riskFreeRate, double timeToMaturity, double optionValue, char optionType)
        : S(stockPrice), K(strikePrice), r(riskFreeRate), T(timeToMaturity), optionValue(optionValue), optionType(optionType) {}

    double calculateImpliedVolatility() {
        double sigma = 0.5; // 초기 추정값

        for (int i = 0; i < 100; ++i) {
            double d1 = calculateD1(sigma);
            double d2 = calculateD2(sigma);
            
            double estimatedOptionValue;
            if (optionType == 'C') {
                estimatedOptionValue = calculateCall(d1, d2);
            }
            else if (optionType == 'P') {
                estimatedOptionValue = calculatePut(d1, d2);
            }

            double vega = S * sqrt(T) * n(d1);

            double error = estimatedOptionValue - optionValue;

            if (fabs(error) < 0.001)
                break;

            sigma = sigma - error / vega;
        }

        return sigma;
    }

private:
    double calculateCall(double d1, double d2) {
        double callValue = S * N(d1) - K * exp(-r * T) * N(d2);
        return callValue;
    }

    double calculatePut(double d1, double d2) {
        double putValue = K * exp(-r * T) * N(-d2) - S * N(-d1);
        return putValue;
    }

    double N(double x) {
        return 0.5 * (1 + erf(x / sqrt(2)));
    }

    double n(double x) {
        return exp(-0.5 * pow(x, 2)) / sqrt(2 * M_PI);
    }

    double calculateD1(double sigma) {
        return (log(S / K) + (r + 0.5 * pow(sigma, 2)) * T) / (sigma * sqrt(T));
    }

    double calculateD2(double sigma) {
        return calculateD1(sigma) - sigma * sqrt(T);
    }
};

int main() {
    double S = 26874.0;
    double K = 26750.0;
    double r = 0.03;
    double T = 0.01;
    double optionValue = 188.0;
    char optionType = 'C';

    Option option(S, K, r, T, optionValue, optionType);

    double impliedVolatility = option.calculateImpliedVolatility();

    std::cout << "Implied Volatility: " << impliedVolatility << std::endl;

    return 0;
}
