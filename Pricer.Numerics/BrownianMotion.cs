using System;
using System.Collections.Generic;
using MathNet.Numerics.Distributions;

namespace Pricer.Numerics
{
    public sealed class BrownianPath
    {
        public double[] Times { get; } // tijdstippen waarop de Brownian motion is geëvalueerd
        public double[] Values { get; } // waarden van de Brownian motion op die tijdstippen

        public BrownianPath(double[] times, double[] values)
        {
            Times = times ?? throw new ArgumentNullException(nameof(times)); // als times niet nul is,zet het in Times, anders gooi een foutmelding
            Values = values ?? throw new ArgumentNullException(nameof(values));

            // controleren of je evenveel tijdstippen als waarden hebt
            if (times.Length != values.Length)
                throw new ArgumentException("Times and Values must have the same length.");
        }
    }

    public static class BrownianMotion
    {
        // 1 pad genereren
        public static BrownianPath GeneratePath(double maturity, int numberOfSteps, int? seed = null)
        {
            if (maturity <= 0.0)
                throw new ArgumentException("Maturity must be strictly positive.", nameof(maturity));

            if (numberOfSteps <= 0)
                throw new ArgumentException("Number of steps must be strictly positive.", nameof(numberOfSteps));

            Random rng = seed.HasValue ? new Random(seed.Value) : new Random(); // als er een seed is, gebruik die, anders maak een nieuwe Random aan
            // met een vaste seed krijg je telkens exact hetzelfde pad terug

            double dt = maturity / numberOfSteps; // tijdstapgrootte (T/N = delta t)
            double sqrtDt = Math.Sqrt(dt); 

            double[] times = new double[numberOfSteps + 1]; // +1 omdat we ook tijdstip 0 willen opnemen
            double[] values = new double[numberOfSteps + 1];

            // startwaarde van de Brownian motion is altijd 0 op tijdstip 0 (W(0) = 0)
            times[0] = 0.0;
            values[0] = 0.0;

            for (int i = 1; i <= numberOfSteps; i++) // pad stap voor stap opbouwen (begin bij i=1 omdat i=0 al is ingevuld)
            {
                double z = Normal.Sample(rng, 0.0, 1.0); // trek een willekeurige waarde uit de standaardnormale verdeling

                times[i] = i * dt; // t_i = i*delta t
                values[i] = values[i - 1] + sqrtDt * z; // W(t_i) = W(t_{i-1}) + sqrt(delta t)*Z_i, waarbij Z_i ~ N(0,1)
            }

            return new BrownianPath(times, values);
        }

        // meerdere paden genereren door meerdere keren GeneratePath aan te roepen, met verschillende seeds
        public static List<BrownianPath> GeneratePaths(int numberOfPaths, double maturity, int numberOfSteps, int? seed = null)
        {
            if (numberOfPaths <= 0)
                throw new ArgumentException("Number of paths must be strictly positive.", nameof(numberOfPaths));

            var paths = new List<BrownianPath>(numberOfPaths); // lijst waar alle paden in opgeslagen worden
            Random masterRng = seed.HasValue ? new Random(seed.Value) : new Random(); // voor elk pad een andere seed (anders kan je identieke paden krijgen)

            for (int i = 0; i < numberOfPaths; i++) // loop over elk pad dat we willen genereren
            {
                paths.Add(GeneratePath(maturity, numberOfSteps, masterRng.Next())); // genereer een nieuw pad met een nieuwe seed en voeg het toe aan de lijst
            }

            return paths;
        }

        public static double ComputeQuadraticVariation(BrownianPath path)
        {
            if (path is null)
                throw new ArgumentNullException(nameof(path));

            double qv = 0.0;

            for (int i = 1; i < path.Values.Length; i++) // loop over alle incrementen
            {
                double increment = path.Values[i] - path.Values[i - 1]; // dW_i = W(t_i) - W(t_{i-1})
                qv += increment * increment;
            }

            return qv; // qv = sum_i (dW_i)^2
        }
        
    }
}