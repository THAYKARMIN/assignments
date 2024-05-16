import random


def guess_input():
    number = random.randrange(1, 21)
    guess = input("Guess a number between 1-20: ")
    return number, guess


def handle_letters(guess, number):
    if guess == "x":
        return None
    elif guess == "n":
        return "new_game"
    elif guess == "s":
        print(f"The correct number is: {number}")
        return "cheating"
    else:
        return guess


def handle_numbers(guess, number):
    count = 0
    guess = int(guess)
    while guess != number:
        if guess > number:
            print("Your guess is too big, try again")
        elif guess < number:
            print("Your guess is too small, try again")
        count += 1
        guess = input("Guess again: ")

        guess = handle_letters(guess, number)
        if guess is None or guess == "new_game" or guess == "cheating":
            return guess

        guess = int(guess)

    return count


def guess_output(count):
    print("Congratulations, your guess is correct!")
    print(f"It took you only {count} guesses to get there!")


def main():
    while True:
        number, guess = guess_input()

        guess = handle_letters(guess, number)
        if guess is None or guess == "new_game" or guess == "cheating":
            if guess == "new_game":
                print("let's start a new game")
                continue
            else:
                break

        count = handle_numbers(guess, number)

        if type(count) == int:
            guess_output(count)
        elif count == "new_game":
            print("let's start a new game")
            continue
        else:
            break

        answer = input("Do you want to start a new game? (yes/no): ")
        if answer.lower() != "yes":
            break

    print("Game over. Thank you for playing.")


main()
