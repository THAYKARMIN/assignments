def character_count(text):
    "counts the character number in the file"
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
    "counts how many lines in the file"
    line_counter = 1
    for line in text:
        if line == "\n":
            line_counter += 1

    return line_counter


def word_count(text):
    "counts the word number in the file"
    word_counter = 1
    for ch in text:
        if ch == " " or ch == "\n":
            word_counter += 1

    return word_counter


def display(text):
    "displays analysis of character, line, and word count of the file"

    characters = character_count(text)
    lines = line_count(text)
    words = word_count(text)

    text_without_newlines = text.replace("\n", "")
    print("Characters number:", len(text_without_newlines))
    for char, count in sorted(characters.items()):
        print(f"'{char}': {count}")
    print(f"Lines number: {lines}")
    print(f"Words number: {words}")
