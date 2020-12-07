import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()
sns.set_context('talk')


def plot_strong_scaling(data):

    data.iloc[:4].plot(xlabel='Nodes', ylabel='Speedup', marker='o')


def plot_perfect_scaling():
    values = [1, 16]
    locations = [1, 16]

    plt.plot(values, locations, marker='', linestyle='--')


def main():
    df = pd.read_csv('scaling_measurements.csv', index_col='index')
    ser = df['1000']
    ser = ser.iloc[0] / ser

    fig = plt.figure(figsize=(10, 8), dpi=100)
    plt.title('Strong Scaling')

    plot_perfect_scaling()
    plot_strong_scaling(ser)

    plt.legend(['Perfect', 'Implementation'])
    plt.savefig('strong_scaling.pdf', dpi=500)


if __name__ == '__main__':
    main()
