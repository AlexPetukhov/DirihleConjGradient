using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace DirihleConjGradient
{
    class TableCreator
    {
        private uint width;
        private uint height;
        private uint stepH;
        private uint stepW;


        public TableCreator(uint width = 0u, uint height = 0u)
        {
            Init(width, height);
        }

        public void Init(uint width, uint height)
        {
            this.width = width;
            this.height = height;
            this.stepW = (0u == (width / 100u)) ? (1u) : (2u * (width / 100u));
            this.stepH = (0u == (height / 100u)) ? (1u) : (2u * (height / 100u));
        }

        public void Fill<T>(DataGridView table, T[,] values)
        {
            Init(table);

            for (uint i = 0u; i < table.RowCount; ++i)
            {
                table[0, (int)i].Value = height - i * stepH - 1u;

                for (uint j = 1u; j < table.ColumnCount; ++j)
                {
                    uint ii = Convert.ToUInt32(table[0, (int)i].Value);
                    uint jj = Convert.ToUInt32(table.Columns[(int)j].Name);

                    table[(int)j, (int)i].Value = values[jj, ii];
                }
            }
        }

        private void Init(DataGridView table)
        {
            table.Rows.Clear();
            table.Columns.Clear();

            table.Columns.Add("", "");
            table.Rows.Add((int)((height / stepH) - 1u));

            if (1u != stepH)
            {
                table.Rows.Add();
            }

            for (uint j = 0u; j < width; j += stepW)
            {
                table.Columns.Add(j.ToString(), j.ToString());
            }
        }
    }
}
