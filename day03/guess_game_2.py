import random


def guess_input():
    number = random.randrange(1, 21)
    guess = input("guess a number between 1-20: ")
    return number, guess


def guess_game():
    number, guess = guess_input()

    if guess == "x":
        return None

    if guess == "n":
        answer = input("Do you want to start a new game? (yes/no): ")
        if answer.lower() == "yes":
            return "new_game"
        else:
            guess = input("guess again: ")

    if guess == "s":
        print(f"The correct number is: {number}")
        return "cheating"

    count = 0
    guess = int(guess)
    while guess != number:

        if guess > number:
            print("your guess is too big, try again")
        elif guess < number:
            print("your guess is too small, try again")
        count += 1
        guess = input("guess again: ")

        if guess == "x":
            return None

        if guess == "n":
            answer = input("Do you want to start a new game? (yes/no): ")
            if answer.lower() == "yes":
                print("Lets start over")
                return "new_game"
            else:
                guess = input("guess again: ")

        if guess == "s":
            print(f"The correct number is: {number}")
            return

        else:
            guess = int(guess)

    else:
        return count


def guess_output(count):
    print("congratulation, your guess is correct!")
    print(f"It took you only {count} guesse to get there!")


def main():
    while True:
        guess_result = guess_game()
        if guess_result is None:
            break
        elif guess_result == "new_game":
            guess_result = guess_game()
        elif guess_result == "cheating":
            break
        guess_output(guess_result)
        answer = input("Do you want to play again? (yes/no): ")
        if answer.lower() != "yes":
            break
    print("Game over, Thank you for playing")


main()
