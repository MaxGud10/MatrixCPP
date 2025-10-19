# Matrix
Проектирование матрицы и поиск определителя

## Сборка

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## Использование 
Перейти в папку build:
```sh
cd build 
```

Для запуска программы:
```sh
./matrix
```
Для запуcка unit тестирования:
```sh
cd tests
./unit_tests
```

## End to end тестирование
Перейти в папку с тестами:
```sh
cd tests
```

Запустить end to end тестирование:
```sh
cd e2e
chmod +x ./end_to_end_test.sh
./end_to_end_test.sh
```