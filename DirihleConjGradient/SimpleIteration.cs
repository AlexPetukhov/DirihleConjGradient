using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DirihleConjGradient
{
    class SimpleIteration : Method
    {
        double tau;

        public SimpleIteration() { }

        public SimpleIteration( double Xo,
                                double Xn,
                                double Yo,
                                double Yn,
                                uint N,
                                uint M,
                                ApproximationType approximationType)
                : base(Xo, Xn, Yo, Yn, N, M, approximationType) { }

        public override void FirstStep()
        {
            for (uint i = 1; i < N; ++i)
            {
                for (uint j = 1; j < M; ++j)
                {
                    residual[i, j] = a2 * data[i, j] +
                                     h2 * (data[i - 1, j] + data[i + 1, j]) +
                                     k2 * (data[i, j - 1] + data[i, j + 1]) +
                                     function[i, j];
                }
            }
        }

        public override void InitMethod()
        {
            double lambdaMin = 4.0 / (h * h) * Math.Sin(Math.PI / (2u * N)) * Math.Sin(Math.PI / (2u * N)) +
                               4.0 / (k * k) * Math.Sin(Math.PI / (2u * M)) * Math.Sin(Math.PI / (2u * M));

            double lambdaMax = 4.0 / (h * h) * Math.Sin(Math.PI * (N - 1u) / (2.0 * N)) * Math.Sin(Math.PI * (N - 1u) / (2u * N)) +
                               4.0 / (k * k) * Math.Sin(Math.PI * (M - 1u) / (2.0 * M)) * Math.Sin(Math.PI * (M - 1u) / (2u * M));

            tau = 2.0 / (lambdaMin + lambdaMax);
        }

        public override void Run(ref uint maxIter, ref double maxAccuracy)
        {
            double accuracy;
            uint counter = 0;

            Approximate();

            for (uint i = 0u; i < N + 1; ++i)
            {
                for (uint j = 0u; j < M + 1; ++j)
                {
                    data[i, j] = V(i, j);
                    function[i, j] = Function(X(i), Y(j));
                }
            }

            do
            {
                accuracy = 0;

                FirstStep();

                for (int j = 1; j < M; ++j)
                {
                    for (int i = 1; i < N; ++i)
                    {
                        double cur = data[i, j];
                        double next = GetNextIterValue(i, j);
                        double curAcc = Math.Abs(cur - next);
                        accuracy = Math.Max(curAcc, accuracy);

                        data[i, j] = next;
                    }

                }
            } while ((maxIter > ++counter) && (accuracy >= maxAccuracy));

            maxIter = counter;
            maxAccuracy = accuracy;

        }

        public override void Step()
        {
            for (int i = 0; i <= N; i++)
            {
                for (int j = 0; j <= M; j++)
                {
                    dataPrev[i, j] = data[i, j];
                    data[i, j] = GetNextIterValue(i, j);
                }
            }
        }

        private double GetNextIterValue(int i, int j)
        {
            double nextValue = data[i, j] - tau * residual[i, j];

            return nextValue;
        }

        public override string GetMethodParameter()
        {
            return tau.ToString();
        }
    }
}
