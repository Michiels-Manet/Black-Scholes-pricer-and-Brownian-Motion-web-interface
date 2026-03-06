using System;
using MathNet.Numerics.Distributions;

namespace Pricer.Numerics
{
    public enum DiffusiveScalingRegime
    {
        Freezing,
        Brownian,
        Explosion
    }

    public class DiffusiveScalingPath
    {
        public double[] Times { get; }
        public double[] Values { get; }

        public DiffusiveScalingPath(double[] times, double[] values)
        {
            Times = times ?? throw new ArgumentNullException(nameof(times));
            Values = values ?? throw new ArgumentNullException(nameof(values));

            if (times.Length != values.Length)
                throw new ArgumentException("Times and Values must have the same length.");
        }
    }

    public static class DiffusiveScaling
    {
        public static DiffusiveScalingPath GeneratePath(
            double maturity,
            int numberOfSteps,
            double sigma, // volatiliteitsparameter
            double alpha,
            int? seed = null)
        {
            if (maturity <= 0.0)
                throw new ArgumentException("Maturity must be strictly positive.", nameof(maturity));

            if (numberOfSteps <= 0)
                throw new ArgumentException("Number of steps must be strictly positive.", nameof(numberOfSteps));

            if (sigma < 0.0)
                throw new ArgumentException("Sigma must be non-negative.", nameof(sigma));

            Random rng = seed.HasValue ? new Random(seed.Value) : new Random();

            double dt = maturity / numberOfSteps;

            // Uit de theorie:
            // Var(ΔX) = sigma^2 * (dt)^alpha
            // dus nemen we ΔX = sigma * (dt)^(alpha/2) * Z, met Z ~ N(0,1)
            double incrementScale = sigma * Math.Pow(dt, 0.5 * alpha);

            double[] times = new double[numberOfSteps + 1];
            double[] values = new double[numberOfSteps + 1];

            times[0] = 0.0;
            values[0] = 0.0;

            for (int i = 1; i <= numberOfSteps; i++)
            {
                double z = Normal.Sample(rng, 0.0, 1.0);
                double increment = incrementScale * z;

                times[i] = i * dt;
                values[i] = values[i - 1] + increment;
            }

            return new DiffusiveScalingPath(times, values);
        }

        public static double TheoreticalVarianceAtTime(double maturity, double sigma, double dt, double alpha)
        {
            if (maturity <= 0.0)
                throw new ArgumentException("Maturity must be strictly positive.", nameof(maturity));

            if (sigma < 0.0)
                throw new ArgumentException("Sigma must be non-negative.", nameof(sigma));

            if (dt <= 0.0)
                throw new ArgumentException("dt must be strictly positive.", nameof(dt));

            // Pdf:
            // Var(X_t) = sigma^2 * t * (dt)^(alpha - 1)
            return sigma * sigma * maturity * Math.Pow(dt, alpha - 1.0);
        }

        public static double EmpiricalQuadraticVariation(DiffusiveScalingPath path)
        {
            if (path is null)
                throw new ArgumentNullException(nameof(path));

            double qv = 0.0;

            for (int i = 1; i < path.Values.Length; i++)
            {
                double increment = path.Values[i] - path.Values[i - 1];
                qv += increment * increment;
            }

            return qv;
        }

        public static DiffusiveScalingRegime GetRegime(double alpha, double tolerance = 1e-12)
        {
            if (alpha > 1.0 + tolerance)
                return DiffusiveScalingRegime.Freezing;

            if (alpha < 1.0 - tolerance)
                return DiffusiveScalingRegime.Explosion;

            return DiffusiveScalingRegime.Brownian;
        }

        public static string GetRegimeDescription(double alpha, double tolerance = 1e-12)
        {
            DiffusiveScalingRegime regime = GetRegime(alpha, tolerance);

            return regime switch
            {
                DiffusiveScalingRegime.Freezing =>
                    "alpha > 1: freezing regime, variance goes to zero in the limit.",

                DiffusiveScalingRegime.Explosion =>
                    "alpha < 1: explosion regime, variance blows up in the limit.",

                _ =>
                    "alpha = 1: diffusive scaling regime, converging to Brownian motion."
            };
        }
    }
}