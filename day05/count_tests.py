import sys
import mycount

def main():

    if len(sys.argv) != 2:
        exit("Please enter in the command line: count.py filename")

    filename = sys.argv[1]
    with open(filename) as fh:
        text = fh.read()

    mycount.display(text)


main()
