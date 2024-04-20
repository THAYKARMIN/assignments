import random

number = random.randrange(1, 21)
guess = int(input("guess a number between 1-20: "))
count = 0

while guess != number:
    if guess > number:
        print("your guess is too big, try again")
        guess = int(input("guess again: "))
        count = count + 1
    elif guess < number:
        print("your guess is too small, try again")
        guess = int(input("guess again: "))
        count = count + 1

if guess == number:
    print("congratulation, your guess is correct!")
    print(f"It took you only {count} guesse to get there!")
