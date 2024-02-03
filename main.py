import argparse
from AdvancedTensor.code import Code


def main():
    parser = argparse.ArgumentParser(description="This is ")
    parser.add_argument('--segment_size', type=int, help='Description for argument 1')
    parser.add_argument('--number_of_segments', type=str, help='Description for argument 2')
    args = parser.parse_args()


if __name__ == "__main__":
    main()