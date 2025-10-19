#!/bin/bash

files=(*.in)

if [ ${#files[@]} -eq 0 ]; then
    echo "Нет файлов с расширением .in."
    exit 1
fi

for file in "${files[@]}"; do 
    echo "Обработка файла: $file"
    echo "Определитель матрицы:"
    ../../build/matrix < "$file"
    echo "-----------------------------------"
done