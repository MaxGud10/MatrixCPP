#!/usr/bin/env bash
set -euo pipefail

BINARY="${MATRIX_CHAIN_BIN:-../../build/matrix_chain}"

if [[ ! -x "$BINARY" ]]; then
  echo "Ошибка: не найден исполняемый файл $BINARY"
  echo "Собери проект: cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build"
  exit 1
fi

shopt -s nullglob
inputs=(test_*.in)

if (( ${#inputs[@]} == 0 )); then
  echo "Нет файлов test_*.in в текущей папке"
  exit 1
fi

all_ok=1

for input_file in "${inputs[@]}"; do
  base="${input_file%.in}"
  expected_file="${base}.out"

  echo "=== ${input_file} ==="

  if [[ ! -f "$expected_file" ]]; then
    echo "Нет ожидаемого файла: $expected_file"
    all_ok=0
    continue
  fi

  output="$("$BINARY" < "$input_file")"

  python3 verify_chain_output.py "$expected_file" <<< "$output" || all_ok=0
done

echo
if (( all_ok == 1 )); then
  echo "Все e2e тесты прошли"
  exit 0
else
  echo "Есть упавшие e2e тесты"
  exit 1
fi