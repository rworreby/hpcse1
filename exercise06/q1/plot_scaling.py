import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('darkgrid')


def main():
    df = pd.read_csv('parallel_io_scaling.csv', index_col='Gridsize')
    df.plot()

    plt.title('Scaling of MPI Parallel IO on 16 Ranks')
    plt.ylabel('Time in [s]')
    plt.savefig('parallel_io_scaling.pdf', dpi=500)


if __name__ == '__main__':
    main()
