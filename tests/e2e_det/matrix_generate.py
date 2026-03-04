import random 

tests_number = 10

for test_number in range(tests_number):
    name_of_file = f"test_{test_number + 1:02}.in"

    with open(name_of_file, 'w') as file:
        matrix_degree = random.randint(10, 10)

        test_text = str(matrix_degree) + "\n"

        elems_str = ""

        for elem_number in range(matrix_degree * matrix_degree):
            elem = random.randint(-5, 5)
            elems_str += str(elem) + " "

        elems_str += "\n"

        test_text += elems_str

        file.write(test_text)