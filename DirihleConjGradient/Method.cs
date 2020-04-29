using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace DirihleConjGradient
{
    enum ApproximationType
    {
        ZERO_APPROXIMATION,
        X_INTERPOLATION,
        Y_INTERPOLATION
    }

    abstract class Method
    {

        public uint N;
        public uint M;

        public double Xo;
        public double Xn;
        public double Yo;
        public double Yn;
        public double h2;
        public double k2;
        public double a2;
        public double h;
        public double k;

        public double[,] data;
        public double[,] dataPrev;
        public double[,] function;
        public double[,] residual;
        

        ApproximationType approximationType;

        public Func<double, double> mu1;
        public Func<double, double> mu2;
        public Func<double, double> mu3;
        public Func<double, double> mu4;
        public Func<double, double, double> Function;
        public Func<double, double, double> ExactFunction;

        public Method() { }

        public Method(
            double Xo, double Xn, double Yo, double Yn,
            uint N, uint M, ApproximationType approximationType)
        {
            Init(Xo, Xn, Yo, Yn, N, M, approximationType);
        }

        public abstract string GetMethodParameter();

        public void Init(
            double Xo, double Xn, double Yo, double Yn,
            uint N, uint M, ApproximationType approximationType)
        {
            this.data = new double[N + 1u, M + 1u];
            this.dataPrev = new double[N + 1u, M + 1u];
            this.function = new double[N + 1u, M + 1u];
            this.residual = new double[N + 1, M + 1];
            this.Xo = Xo;
            this.Yo = Yo;
            this.Xn = Xn;
            this.Yn = Yn;
            this.N = N;
            this.M = M;
            this.h = (this.Xn - this.Xo) / N;
            this.k = (this.Yn - this.Yo) / M;
            this.h2 = -Math.Pow(N / (Xn - Xo), 2);
            this.k2 = -Math.Pow(M / (Yn - Yo), 2);
            this.a2 = -2.0 * (h2 + k2);
            this.approximationType = approximationType;
            InitMethod();
        }

        public abstract void Run(ref uint maxIter, ref double maxAccuracy);

        public abstract void Step();

        public abstract void FirstStep();

        public abstract void InitMethod();

        public double[,] GetExactTable()
        {
            double[,] exact = new double[N + 1u, M + 1u];

            for (uint i = 0u; i < N + 1u; ++i)
            {
                for (uint j = 0u; j < M + 1u; ++j)
                {
                    exact[i, j] = ExactFunction(X(i), Y(j));
                }
            }

            return exact;
        }

        public double CalculateResidual()
        {
            double R = 0.0;

            for (uint j = 1u; j < M; ++j)
            {
                for (uint i = 1u; i < N; ++i)
                {
                    R += Math.Pow(
                         a2 * data[i, j] +
                         h2 * (data[i - 1, j] + data[i + 1, j]) +
                         k2 * (data[i, j - 1] + data[i, j + 1]) +
                         Function(X(i), Y(j)), 2.0);
                }
            }

            return Math.Sqrt(R);
        }

        protected void Approximate()
        {
            switch (approximationType)
            {
                case ApproximationType.ZERO_APPROXIMATION:
                    ZeroApproximation();
                    break;

                case ApproximationType.X_INTERPOLATION:
                    XInterpolation();
                    break;

                case ApproximationType.Y_INTERPOLATION:
                    YInterpolation();
                    break;
            }
        }

        public void ZeroApproximation()
        {
            for (uint i = 0u; i < N + 1; ++i)
            {
                for (uint j = 0u; j < M + 1; ++j)
                {
                    data[i, j] = 0.0;
                }
            }
        }

        public void XInterpolation()
        {
            for (uint i = 1u; i < N; ++i)
            {
                for (uint j = 1u; j < M; ++j)
                {
                    data[i, j] = ((Xo + i * h) - Xo) / (Xn - Xo) * mu2(Yo + j * k) +
                                 ((Xo + i * h) - Xn) / (Xo - Xn) * mu1(Yo + j * k);
                }
            }
        }

        public void YInterpolation()
        {
            for (uint i = 1u; i < N; ++i)
            {
                for (uint j = 1u; j < M; ++j)
                {
                    data[i, j] = ((Yo + j * k) - Yo) / (Yn - Yo) * mu4(Xo + i * h) +
                              ((Yo + j * k) - Yn) / (Yo - Yn) * mu3(Xo + i * h);
                }
            }
        }

        public void SetFunctions(
            Func<double, double> mu1, Func<double, double> mu2,
            Func<double, double> mu3, Func<double, double> mu4,
            Func<double, double, double> Function,
            Func<double, double, double> ExactFunction = null)
        {
            this.mu1 = mu1;
            this.mu2 = mu2;
            this.mu3 = mu3;
            this.mu4 = mu4;
            this.Function = Function;
            this.ExactFunction = ExactFunction;
        }
        public double V(uint i, uint j)
        {
            if (0u == i)
            {
                return mu1(Y(j));
            }
            if (N == i)
            {
                return mu2(Y(j));
            }
            if (0u == j)
            {
                return mu3(X(i));
            }
            if (M == j)
            {
                return mu4(X(i));
            }

            return data[i, j];
        }

        public double X(uint i) => Xo + i * h;

        public double Y(uint j) => Yo + j * k;

       
    }
}
