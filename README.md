# Matrix
Проект создан в рамках курса по 
C++ Константина Владимирова. Задание состояло в том, чтобы реализовать класс, представляющий матрицу и предоставляющий метод для вычисления определителя.

## Требование
- CMake версии 3.10 (или выше)

## Как утсановить
```bash
git clone https://github.com/MaxGud10/MatrixCPP
cd Matrix
```

## Сборка
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## Использование 
Для запуска
```bash
./build/matrix
```

## Unit тесты
Для запуcка unit тестирования:
```sh
./build/tests/unit_tests
```

## End to end тестирование
Перейти в папку с тестами:
```sh
cd tests/e2e
```

Запустить end to end тестирование:
```sh
chmod +x ./end_to_end_test.sh
./end_to_end_test.sh
```