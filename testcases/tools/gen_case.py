import random as r

# kind of bad but works for generating test cases
for i in range(36):
    print(f"{r.randint(5,10):>2}e+24 " + " ".join([f"{r.randint(-10, 10)+r.random():>3}" for x in range(4)]))