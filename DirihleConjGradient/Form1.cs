using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace DirihleConjGradient
{
    public partial class Form1 : Form
    {
        delegate void SetChangeButtonVisibilityCallBack(Button button, bool status);
        delegate void SetChangeTextBoxValueCallBack(TextBox textBox, string value);
        delegate void SetChangeLabelValueCallBack(Label label, string value);
        delegate void SetChangeTableValuesCallBack(DataGridView table, double[,] values);

        private double Xo;
        private double Xn;
        private double Yo;
        private double Yn;
        private double acc_max;
        private ApproximationType approximationType;
        private TableCreator tableCreator;

        private uint Nmax;
        private uint N;
        private uint M;

        private double mu1Test(double y) => Math.Exp(1.0 - Math.Pow(Xo, 2) - Math.Pow(y, 2));
        private double mu2Test(double y) => Math.Exp(1.0 - Math.Pow(Xn, 2) - Math.Pow(y, 2));
        private double mu3Test(double x) => Math.Exp(1.0 - Math.Pow(x, 2) - Math.Pow(Yo, 2));
        private double mu4Test(double x) => Math.Exp(1.0 - Math.Pow(x, 2) - Math.Pow(Yn, 2));
        private double FunctionTest(double x, double y) =>
            Math.Exp(-(x * x) - (y * y) + 1.0) * (4.0 * x * x - 2.0 + 4.0 * y * y - 2.0);
        private double ExactFunction(double x, double y) =>
            Math.Exp(1.0 - (x * x) - (y * y));

        private double mu1Main(double y) => 1.0 - Math.Pow(y, 2.0);
        private double mu2Main(double y) => (1.0 - Math.Pow(y, 2.0)) * Math.Exp(y);
        private double mu3Main(double x) => 1.0 - Math.Pow(x, 2.0);
        private double mu4Main(double x) => 1.0 - Math.Pow(x, 2.0);
        private double FunctionMain(double x, double y) => -Math.Abs(x * x - y * y);

        public Form1()
        {
            InitializeComponent();
            Nbox.Minimum = decimal.MinValue;
            Nbox.Maximum = decimal.MaxValue;
            Mbox.Minimum = decimal.MinValue;
            Mbox.Maximum = decimal.MaxValue;

            NmaxBox.Maximum = decimal.MaxValue;
            AccuracyBox.Text = "0,0001";
            NmaxBox.Value = 100;
            Nbox.Value = 10;
            Mbox.Value = 10;

            ZeroApprocsimationCheckBox.Checked = true;
            Table.RowHeadersVisible = false;
            TableMain.RowHeadersVisible = false;
            TableExact.RowHeadersVisible = false;
            TableDiffTest.RowHeadersVisible = false;
            TableHalf.RowHeadersVisible = false;
            TableDiffMain.RowHeadersVisible = false;

            label10.Text = "При решении основной задачи \n" +
                          "        с половинным шагом";

            tableCreator = new TableCreator();
        }

        private void SolveTestTask_Click(object sender, EventArgs e)
        {
            if (!ParseArguments())
            {
                MessageBox.Show("ERROR!", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                return;
            }
            
            double maxDif = 0.0;
            double maxAcc = acc_max;
            uint iterCount = Nmax;
            uint maxI = 0u;
            uint maxJ = 0u;

            Method task = new Method();
            

            task.Init(Xo, Xn, Yo, Yn, N, M, approximationType);
            task.SetFunctions(mu1Test, mu2Test, mu3Test, mu4Test, FunctionTest, ExactFunction);

            task.Run(ref iterCount, ref maxAcc);


            FindMax(
                CalculateDifferenceTableForTestTask(
                    task.data, task.GetExactTable(), task.N, task.M),
                out maxDif, out maxI, out maxJ);
            
            ChangeLabelValue(ResidualTextBox, task.CalculateResidual().ToString());
            ChangeLabelValue(IterLabel, iterCount.ToString());
            ChangeLabelValue(AccMaxLabel, maxAcc.ToString());
            ChangeLabelValue(maxDifLabel, maxDif.ToString());
            ChangeLabelValue(DotLabelTest, "Соответствует узлу x = " + Math.Abs(Math.Round(task.X(maxI), 3)).ToString() +
                                                               "  y = " + Math.Abs(Math.Round(task.Y(maxJ), 3)).ToString());

            tableCreator.Init(N + 1u, M + 1u);

            ChangeTableValues(Table, task.data);
            ChangeTableValues(TableExact, task.GetExactTable());
            ChangeTableValues(TableDiffTest,
                CalculateDifferenceTableForTestTask(task.data, task.GetExactTable(), N, M));
        }
        private void FindMax(double[,] a, out double maxDif, out uint maxX, out uint maxY)
        {
            maxDif = 0.0;
            maxX = 0u;
            maxY = 0u;

            for (uint i = 1u; i < N; ++i)
            {
                for (uint j = 1u; j < M; ++j)
                {
                    if (a[i, j] > maxDif)
                    {
                        maxDif = a[i, j];
                        maxX = i;
                        maxY = j;
                    }
                }
            }
        }
        private double[,] CalculateDifferenceTableForTestTask(double[,] solutionTable, double[,] exactTable, uint _N, uint _M)
        {
            double[,] difference = new double[_N + 1u, _M + 1u];

            for (uint i = 0; i < N + 1; ++i)
            {
                for (uint j = 0; j < M + 1; ++j)
                {
                    difference[i, j] = Math.Abs(exactTable[i, j] - solutionTable[i, j]);
                }
            }

            return difference;
        }
        private bool ParseArguments()
        {
            Xo = (double)(-1);
            Xn = (double)(1);
            Yo = (double)(-1);
            Yn = (double)(1);
            N = (uint)Nbox.Value;
            M = (uint)Mbox.Value;
            Nmax = (uint)NmaxBox.Value;
            return double.TryParse(AccuracyBox.Text, out acc_max);
        }

        private void SolveMainTask_Click(object sender, EventArgs e)
        {
            if (!ParseArguments())
            {
                MessageBox.Show("ERROR!", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                return;
            }

            double maxDif = 0.0;
            double maxAccMain = acc_max;
            double maxAccHalf = acc_max;
            uint IterMain = Nmax;
            uint IterHalf = Nmax;
            uint maxI = 0u;
            uint maxJ = 0u;


            Method main = new Method();
            Method half = new Method();

            main.Init(Xo, Xn, Yo, Yn, N, M, approximationType);
            main.SetFunctions(mu1Main, mu2Main, mu3Main, mu4Main, FunctionMain);
            main.Run(ref IterMain, ref maxAccMain);

            half.Init(Xo, Xn, Yo, Yn, N * 2u, M * 2u, approximationType);
            half.SetFunctions(mu1Main, mu2Main, mu3Main, mu4Main, FunctionMain);
            half.Run(ref IterHalf, ref maxAccHalf);

            FindMax(CalculateDifferenceTableForMainTask(main.data, half.data, N, M), out maxDif, out maxI, out maxJ);

            
            ChangeLabelValue(MaxDifLabelMain, maxDif.ToString());
            ChangeLabelValue(ResidualMainTextBox, main.CalculateResidual().ToString());
            ChangeLabelValue(IterLabelMain, IterMain.ToString());
            ChangeLabelValue(AccMaxLabelMain, maxAccMain.ToString());
            ChangeLabelValue(ResidualHalfTextBox, half.CalculateResidual().ToString());
            ChangeLabelValue(IterLabelMainHalf, IterHalf.ToString());
            ChangeLabelValue(AccMaxLabelMainHalf, maxAccHalf.ToString());
            ChangeLabelValue(DotLabel, "Соответствует узлу x = " + Math.Abs(Math.Round(main.X(maxI), 3)).ToString() +
                                                        "  y = " + Math.Abs(Math.Round(main.Y(maxJ), 3)).ToString());

            tableCreator.Init(N * 2 + 1u, M * 2 + 1u);
            ChangeTableValues(TableHalf, half.data);

            tableCreator.Init(N + 1u, M + 1u);
            ChangeTableValues(TableMain, main.data);
            ChangeTableValues(TableDiffMain, CalculateDifferenceTableForMainTask(main.data, half.data, N, M));
        }
        
        private void ChangeTextBoxValue(TextBox textBox, string value)
        {
            if (textBox.InvokeRequired)
            {
                this.Invoke(new SetChangeTextBoxValueCallBack(ChangeTextBoxValue),
                    new object[] { textBox, value });
            }
            else
            {
                
            }
        }
        
        private void ChangeLabelValue(Label label, string value)
        {
            if (label.InvokeRequired)
            {
                this.Invoke(new SetChangeLabelValueCallBack(ChangeLabelValue),
                    new object[] { label, value });
            }
            else
            {
                label.Text = value;
            }
        }
        
        private void ChangeTableValues(DataGridView table, double[,] values)
        {
            if (table.InvokeRequired)
            {
                this.Invoke(new SetChangeTableValuesCallBack(ChangeTableValues),
                    new object[] { table, values });
            }
            else
            {
                tableCreator.Fill<double>(table, values);
            }
        }

        private void AccuracyBox_TextChanged(object sender, EventArgs e)
        {

        }
        private double[,] CalculateDifferenceTableForMainTask(
            double[,] wholeStep, double[,] halfStep, uint _N, uint _M)
        {
            uint width = _M * 2;
            uint height = _N * 2;
            double[,] difference = new double[N + 1u, M + 1u];

            for (int iStep = 1, iHalf = 2; (iStep < height / 2) && (iHalf < height); ++iStep, iHalf += 2)
            {
                for (int jStep = 1, jHalf = 2; (jStep < width / 2) && (jHalf < width); ++jStep, jHalf += 2)
                {
                    difference[iStep, jStep] =
                        Math.Abs(wholeStep[iStep, jStep] - halfStep[iHalf, jHalf]);
                }
            }

            return difference;
        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }
    }
}
