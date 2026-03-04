#!/usr/bin/env python3
import sys
import math

TOL = 1e-6

def parse_order(line: str):
    line = line.strip()
    if line == "":
        return []
    return [int(x) for x in line.split()]

def main():
    if len(sys.argv) != 2:
        print("Usage: verify_chain_output.py expected.out", file=sys.stderr)
        return 2

    expected_path = sys.argv[1]
    with open(expected_path, "r", encoding="utf-8") as f:
        expected_lines = [ln.rstrip("\n") for ln in f.readlines()]

    actual_lines = [ln.rstrip("\n") for ln in sys.stdin.read().splitlines()]

    if len(expected_lines) != 2:
        print(f"Bad expected file format: {expected_path} (need 2 lines)", file=sys.stderr)
        return 2

    if len(actual_lines) != 2:
        print("FAIL: program output must be exactly 2 lines:", file=sys.stderr)
        print("  line1: order", file=sys.stderr)
        print("  line2: factor", file=sys.stderr)
        print("Actual output was:", file=sys.stderr)
        for ln in actual_lines:
            print(ln, file=sys.stderr)
        return 1

    expected_order = parse_order(expected_lines[0])
    actual_order = parse_order(actual_lines[0])

    if expected_order != actual_order:
        print("FAIL: order mismatch", file=sys.stderr)
        print(f"Expected: {expected_order}", file=sys.stderr)
        print(f"Actual:   {actual_order}", file=sys.stderr)
        return 1

    try:
        expected_factor = float(expected_lines[1].strip())
        actual_factor = float(actual_lines[1].strip())
    except ValueError:
        print("FAIL: factor line is not a number", file=sys.stderr)
        print(f"Expected factor line: {expected_lines[1]!r}", file=sys.stderr)
        print(f"Actual factor line:   {actual_lines[1]!r}", file=sys.stderr)
        return 1

    if not math.isfinite(actual_factor):
        print("FAIL: factor is not finite", file=sys.stderr)
        return 1

    if abs(expected_factor - actual_factor) > TOL:
        print("FAIL: factor mismatch", file=sys.stderr)
        print(f"Expected: {expected_factor:.6f}", file=sys.stderr)
        print(f"Actual:   {actual_factor:.6f}", file=sys.stderr)
        print(f"Tolerance: {TOL}", file=sys.stderr)
        return 1

    print("OK")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())