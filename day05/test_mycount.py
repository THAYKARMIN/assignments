import mycount
import pytest


@pytest.fixture
def file_content():
    filename = "file.txt"
    with open(filename) as fh:
        return fh.read()


@pytest.mark.parametrize("expected", [{
    'c': 71, 'r': 137, 'a': 141, 'z': 36, 'y': 141, ' ': 731, 'i': 181, 'n': 212, 
    'l': 123, 'o': 432, 'v': 37, 'e': 245, 'b': 38, 's': 122, '(': 30, 'w': 60, 
    ',': 123, ')': 30, 't': 186, "'": 54, 'g': 147, 'h': 213, 'm': 87, 'd': 38, 
    'u': 120, '?': 4, 'k': 59, 'p': 26, 'f': 20, 'x': 4, 'j': 3, 'q': 1, '-': 2, 
    '"': 2
}])
def test_character(file_content, expected):
    result = mycount.character_count(file_content)
    assert result == expected


@pytest.mark.parametrize("expected", [95])
def test_line(file_content, expected):
    result = mycount.line_count(file_content)
    assert result == expected


@pytest.mark.parametrize("expected", [826])
def test_word(file_content, expected):
    result = mycount.word_count(file_content)
    assert result == expected
