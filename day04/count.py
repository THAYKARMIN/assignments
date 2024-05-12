import sys


def character_count(text):
    character_counter = {}
    for ch in text:
        if ch == "\n":
            continue
        ch = ch.lower()
        if ch not in character_counter:
            character_counter[ch] = 1
        else:
            character_counter[ch] += 1

    return character_counter


def line_count(text):
    line_counter = 1
    for line in text:
        if line == "\n":
            line_counter += 1

    return line_counter


def word_count(text):
    word_counter = 1
    for ch in text:
        if ch == " " or ch == "\n":
            word_counter += 1

    return word_counter


def display(text):

    characters = character_count(text)
    lines = line_count(text)
    words = word_count(text)

    text_without_newlines = text.replace("\n", "")
    print("Characters number:", len(text_without_newlines))
    for char, count in sorted(characters.items()):
        print(f"'{char}': {count}")
    print(f"Lines number: {lines}")
    print(f"Words number: {words}")


def main():

    if len(sys.argv) != 2:
        exit("Please enter in the command line: count.py filename")

    filename = sys.argv[1]
    with open(filename) as fh:
        text = fh.read()

    display(text)


main()
