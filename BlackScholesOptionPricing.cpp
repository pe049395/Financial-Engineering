#include <iostream>
#include <cmath>

class Option {
private:
    double S; // 주식 가격
    double K; // 행사 가격
    double r; // 무위험 이자율
    double T; // 잔여 만기
    double sigma; // 변동성

public:
    Option(double stockPrice, double strikePrice, double riskFreeRate, double timeToMaturity, double volatility)
        : S(stockPrice), K(strikePrice), r(riskFreeRate), T(timeToMaturity), sigma(volatility) {}

    double calculateCall() {
        double d1 = calculateD1();
        double d2 = calculateD2();
        double callValue = S * N(d1) - K * exp(-r * T) * N(d2);
        return callValue;
    }

    double calculatePut() {
        double d1 = calculateD1();
        double d2 = calculateD2();
        double putValue = K * exp(-r * T) * N(-d2) - S * N(-d1);
        return putValue;
    }

    double calculateDeltaCall() {
        double d1 = calculateD1();
        double delta = N(d1);
        return delta;
    }

    double calculateDeltaPut() {
        double d1 = calculateD1();
        double delta = N(d1) - 1;
        return delta;
    }

    double calculateGammaCall() {
        double d1 = calculateD1();
        double gamma = n(d1) / (S * sigma * sqrt(T));
        return gamma;
    }

    double calculateGammaPut() {
        double d1 = calculateD1();
        double gamma = n(d1) / (S * sigma * sqrt(T));
        return gamma;
    }

    double calculateThetaCall() {
        double d1 = calculateD1();
        double d2 = calculateD2();
        double theta = (-S * n(d1) * sigma) / (2 * sqrt(T))
                        - r * K * exp(-r * T) * N(d2);
        return theta / 365.0;
    }

    double calculateThetaPut() {
        double d1 = calculateD1();
        double d2 = calculateD2();
        double theta = (-S * n(d1) * sigma) / (2 * sqrt(T))
                        + r * K * exp(-r * T) * N(-d2);
        return theta / 365.0;
    }

    double calculateVegaCall() {
        double d1 = calculateD1();
        double vega = S * n(d1) * sqrt(T);
        return vega / 100.0;
    }

    double calculateVegaPut() {
        double d1 = calculateD1();
        double vega = S * n(d1) * sqrt(T);
        return vega / 100.0;
    }

    double calculateRhoCall() {
        double d2 = calculateD2();
        double rho = K * T * exp(-r * T) * N(d2);
        return rho / 100.0;
    }

    double calculateRhoPut() {
        double d2 = calculateD2();
        double rho = -K * T * exp(-r * T) * N(-d2);
        return rho / 100.0;
    }

private:
    double N(double x) {
        return 0.5 * (1 + erf(x / sqrt(2)));
    }

    double n(double x) {
        return exp(-0.5 * pow(x, 2)) / sqrt(2 * M_PI);
    }

    double calculateD1() {
        return (log(S / K) + (r + 0.5 * pow(sigma, 2)) * T) / (sigma * sqrt(T));
    }

    double calculateD2() {
        return calculateD1() - sigma * sqrt(T);
    }
};

int main() {
    double S = 26874.0;
    double K = 26750.0;
    double r = 0.03;
    double T = 0.01;
    double sigma = 0.2;

    Option option(S, K, r, T, sigma);

    double callValue = option.calculateCall();
    double callDelta = option.calculateDeltaCall();
    double callGamma = option.calculateGammaCall();
    double callTheta = option.calculateThetaCall();
    double callVega = option.calculateVegaCall();
    double callRho = option.calculateRhoCall();

    double putValue = option.calculatePut();
    double putDelta = option.calculateDeltaPut();
    double putGamma = option.calculateGammaPut();
    double putTheta = option.calculateThetaPut();
    double putVega = option.calculateVegaPut();
    double putRho = option.calculateRhoPut();

    std::cout << "Call Value: " << callValue << std::endl;
    std::cout << "delta: " << callDelta << std::endl;
    std::cout << "gamma: " << callGamma << std::endl;
    std::cout << "theta: " << callTheta << std::endl;
    std::cout << "vega: " << callVega << std::endl;
    std::cout << "rho: " << callRho << std::endl;

    std::cout << "Put Value: " << putValue << std::endl;
    std::cout << "delta: " << putDelta << std::endl;
    std::cout << "gamma: " << putGamma << std::endl;
    std::cout << "theta: " << putTheta << std::endl;
    std::cout << "vega: " << putVega << std::endl;
    std::cout << "rho: " << putRho << std::endl;

    return 0;
}