import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()
sns.set_context('talk')


def main():
    df = pd.read_csv('histogram.csv', index_col=['bin'])

    plt.figure(figsize=(10, 6), dpi=100)

    df.plot(kind='bar', xlabel='Bin', ylabel='Count')
    plt.title('Concentration Distribution at Time t=0.5')
    plt.legend(['Amount'])
    plt.tight_layout()
    plt.savefig('histogram.pdf', dpi=500)


if __name__ == '__main__':
    main()
