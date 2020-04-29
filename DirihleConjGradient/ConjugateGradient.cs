using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DirihleConjGradient
{
    class ConjGradient : Method
    {
        public double[,] hvector;
        public ConjGradient() { }
        public ConjGradient(
            double Xo,
            double Xn,
            double Yo,
            double Yn,
            uint N,
            uint M,
            ApproximationType approximationType)
            : base(Xo, Xn, Yo, Yn, N, M, approximationType)
        {
            this.hvector = new double[N + 1u, M + 1u];
        }

        public double[,] Amult(ref double[,] vec)
        {
            double[,] ans = new double[N + 1, M + 1];
            for (uint i = 1; i < N; ++i)
            {
                for (uint j = 1; j < M; ++j)
                {
                    ans[i, j] = a2 * vec[i, j] +
                                     h2 * (vec[i - 1, j] + vec[i + 1, j]) +
                                     k2 * (vec[i, j - 1] + vec[i, j + 1]);
                }
            }
            return ans;
        }

        public double scalarnoe(ref double[,] v1, ref double[,] v2)
        {
            double ans = new double();

            for (uint i = 0; i < N + 1; ++i)
            {
                for (uint j = 0; j < M + 1; ++j)
                {
                    ans += v1[i, j] * v2[i, j];
                }
            }

            return ans;
        }

        public double[,] AmultB(ref double[,] vec)
        {
            double[,] ans = new double[N + 1, M + 1];
            for (uint i = 1; i < N; ++i)
            {
                for (uint j = 1; j < M; ++j)
                {
                    ans[i, j] = a2 * vec[i, j] +
                                     h2 * (vec[i - 1, j] + vec[i + 1, j]) +
                                     k2 * (vec[i, j - 1] + vec[i, j + 1]) +
                                     function[i, j];
                }
            }
            return ans;
        }

        public override void Run(ref uint maxIter, ref double maxAccuracy)
        {
            double accuracy;
            uint counter = 0u;

            Approximate();

            for (uint i = 0u; i < N + 1; ++i)
            {
                for (uint j = 0u; j < M + 1; ++j)
                {
                    data[i, j] = V(i, j);
                    function[i, j] = Function(X(i), Y(j));
                }
            }

            FirstStep();
            accuracy = 0.0;
            for (uint j = 1u; j < M; ++j)
            {
                for (uint i = 1u; i < N; ++i)
                {
                    double prev = dataPrev[i, j];
                    double cur = data[i, j];
                    double curAcc = Math.Abs(prev - cur);
                    accuracy = Math.Max(accuracy, curAcc);
                }
            }
            if ((maxIter > ++counter) && (accuracy >= maxAccuracy))
            {
                do
                {
                    accuracy = 0.0;

                    Step();

                    for (uint j = 1u; j < M; ++j)
                    {
                        for (uint i = 1u; i < N; ++i)
                        {
                            double prev = dataPrev[i, j];
                            double cur = data[i, j];
                            double curAcc = Math.Abs(prev - cur);
                            accuracy = Math.Max(accuracy, curAcc);
                        }
                    }

                } while ((maxIter > ++counter) && (accuracy >= maxAccuracy));
            }

            maxIter = counter;
            maxAccuracy = accuracy;
        }

        public override void FirstStep()
        {
            residual = AmultB(ref data);
            for (int i = 0; i <= N; i++) for (int j = 0; j <= M; j++) hvector[i, j] = -residual[i, j];

            double[,] AH = Amult(ref hvector);
            double alpha = -scalarnoe(ref hvector, ref hvector) / scalarnoe(ref AH, ref hvector);

            for (int i = 0; i <= N; i++) for (int j = 0; j <= M; j++) dataPrev[i, j] = data[i, j];

            for (int i = 0; i <= N; i++) for (int j = 0; j <= M; j++) data[i, j] = data[i, j] - alpha * hvector[i, j];
        }

        public override void InitMethod()
        {
            this.hvector = new double[N + 1u, M + 1u];
            residual = AmultB(ref data);

            for (uint i = 1; i < N; ++i)
            {
                for (uint j = 1; j < M; ++j)
                {
                    hvector[i, j] = residual[i, j];
                }
            }
        }

        public override void Step()
        {
            residual = AmultB(ref data);
            double[,] AH = Amult(ref hvector);
            double betta = scalarnoe(ref AH, ref residual) / scalarnoe(ref AH, ref hvector);
            for (int i = 0; i <= N; i++) for (int j = 0; j <= M; j++) hvector[i, j] = -residual[i, j] + betta * hvector[i, j];

            AH = Amult(ref hvector);
            double[,] _res = residual;
            for (int i = 0; i <= N; i++) for (int j = 0; j <= M; j++) _res[i, j] = -residual[i, j];
            double alpha = -scalarnoe(ref _res, ref hvector) / scalarnoe(ref AH, ref hvector);

            for (int i = 0; i <= N; i++) for (int j = 0; j <= M; j++) dataPrev[i, j] = data[i, j];
            for (int i = 0; i <= N; i++) for (int j = 0; j <= M; j++) data[i, j] = data[i, j] - alpha * hvector[i, j];

        }

        public override string GetMethodParameter()
        {
            return "-";
        }
    }
}
