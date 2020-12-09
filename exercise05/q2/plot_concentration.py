import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()
sns.set_context('talk')


def main():
    df = pd.read_csv('diagnostics.dat', index_col=['t'])

    plt.figure(figsize=(10, 6), dpi=100)

    df.plot(xlabel='Time', ylabel='Concentration')

    plt.title('Concentration over Time')
    plt.legend(['Amount'])
    plt.tight_layout()
    plt.savefig('concentration.pdf', dpi=500)


if __name__ == '__main__':
    main()
