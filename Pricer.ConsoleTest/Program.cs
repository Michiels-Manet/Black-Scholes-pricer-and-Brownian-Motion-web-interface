using Pricer.Numerics;
using System;
using static Pricer.Numerics.BlackScholes;

Console.WriteLine();
Console.WriteLine("----- BlackScholes pricer Test -----");

// Testprogramma om de Black-Scholes formule te gebruiken en de resultaten te tonen
double S = 100.0;
double K = 100.0;
double r = 0.05;
double q = 0.0;
double T = 1.0;
double sigma = 0.20;

var call = new BlackScholes(OptionType.Call, r, T, sigma, K, S, q);
var put = new BlackScholes(OptionType.Put, r, T, sigma, K, S, q);

double callPrice = call.Price();
double putPrice = put.Price();

double callvega = call.Vega();
double putvega = put.Vega();

double leftSide = callPrice - putPrice;
double rightSide = S * Math.Exp(-q * T) - K * Math.Exp(-r * T);

Console.WriteLine($"Call price   = {callPrice}");
Console.WriteLine($"Put price    = {putPrice}");
Console.WriteLine($"Call Vega : {callvega}");
Console.WriteLine($"Put Vega : {putvega}");
Console.WriteLine($"Call - Put   = {leftSide}");
Console.WriteLine($"Parity RHS   = {rightSide}");
Console.WriteLine($"Difference   = {Math.Abs(leftSide - rightSide)}");

// Testen van de implied volatility calculator
double marketPrice = 10.450583572185565;

double impliedVol = ImpliedVolatilityCalculator.Compute(
    OptionType.Call,
    marketPrice,
    0.05,
    1.0,
    100.0,
    100.0);

Console.WriteLine($"Implied vol = {impliedVol}");

double optionPrice = 100.0 * 0.20 / Math.Sqrt(2.0 * Math.PI);

double sigmaB = ImpliedVolatilityCalculator.BachelierImpliedVolATM(
    optionPrice,
    100.0,
    1.0,
    0.05);

Console.WriteLine($"Bachelier ATM implied vol = {sigmaB}");

Console.WriteLine();
Console.WriteLine("----- Brownian Motion Test -----");

double maturityBM = 1.0;
int numberOfStepsBM = 10;

BrownianPath path = BrownianMotion.GeneratePath(maturityBM, numberOfStepsBM, seed: 42);

Console.WriteLine($"Maturity T      = {maturityBM}");
Console.WriteLine($"Number of steps = {numberOfStepsBM}");
Console.WriteLine($"dt              = {maturityBM / numberOfStepsBM:F4}");
Console.WriteLine();

Console.WriteLine("Time\tValue");
for (int i = 0; i < path.Times.Length; i++)
{
    Console.WriteLine($"{path.Times[i]:F4}\t{path.Values[i]:F6}");
}

double qv = BrownianMotion.ComputeQuadraticVariation(path);

Console.WriteLine();
Console.WriteLine($"Quadratic variation = {qv:F6}");
Console.WriteLine($"Expected around     = {maturityBM:F6}");


Console.WriteLine();
Console.WriteLine("----- Diffusive Scaling Test -----");

double maturityDS = 1.0;
int numberOfStepsDS = 100;
double sigmaDS = 1.0;

// Test 3 regimes
double[] alphas = { 1.2, 1.0, 0.8 };

foreach (double alphaDS in alphas)
{
    var dsPath = DiffusiveScaling.GeneratePath(
    maturityDS,
    numberOfStepsDS,
    sigmaDS,
    alphaDS,
    seed: 42);

    double dtDS = maturityDS / numberOfStepsDS;
    double qvDS = DiffusiveScaling.EmpiricalQuadraticVariation(dsPath);
    double theoVarDS = DiffusiveScaling.TheoreticalVarianceAtTime(maturityDS, sigmaDS, dtDS, alphaDS);
    string regimeDS = DiffusiveScaling.GetRegimeDescription(alphaDS);

    Console.WriteLine();
    Console.WriteLine($"alpha = {alphaDS}");
    Console.WriteLine(regimeDS);
    Console.WriteLine($"T              = {maturityDS}");
    Console.WriteLine($"N              = {numberOfStepsDS}");
    Console.WriteLine($"dt             = {dtDS:F4}");
    Console.WriteLine($"sigma          = {sigmaDS:F4}");
    Console.WriteLine($"Theoretical Var= {theoVarDS:F6}");
    Console.WriteLine($"Quadratic Var  = {qvDS:F6}");

    Console.WriteLine("First few points:");
    Console.WriteLine("Time\tValue");

    int maxRows = Math.Min(10, dsPath.Times.Length);
    for (int i = 0; i < maxRows; i++)
    {
        Console.WriteLine($"{dsPath.Times[i]:F4}\t{dsPath.Values[i]:F6}");
    }
}