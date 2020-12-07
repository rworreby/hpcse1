import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set()
sns.set_context('talk')


def plot_weak_scaling(data):

    data.iloc[:4].plot(xlabel='Nodes', ylabel='Efficiency', marker='o')


def plot_perfect_scaling():
    values = [1, 1]
    locations = [1, 16]

    plt.plot(locations, values, marker='', linestyle='--')


def main():
    df = pd.read_csv('scaling_measurements.csv', index_col='index')
    ser = pd.Series(np.diag(df), index=[*df.index[:4]])
    ser = ser.iloc[0] / ser

    fig = plt.figure(figsize=(10, 8), dpi=100)
    plt.title('Weak Scaling')

    plot_perfect_scaling()
    plot_weak_scaling(ser)

    plt.legend(['Perfect', 'Implementation'])
    plt.savefig('weak_scaling.pdf', dpi=500)


if __name__ == '__main__':
    main()
