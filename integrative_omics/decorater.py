
def f1(func):
    def wrapper():
        print("Started!")
        func()
        print("Ended!")

    return wrapper

@f1
def f():
    print("Hello!")

#f1(f)
#x = f1(f)

f()

