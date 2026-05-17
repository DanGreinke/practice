#!/usr/bin/env python3
"""
change_calculator.py

Calculates the optimal number of bills and coins to return as change
when a customer pays amount P for a cost C.
"""

import sys


DENOMINATIONS = [
    ("$20 bills", 2000),
    ("$10 bills", 1000),
    ("$5 bills", 500),
    ("$1 bills", 100),
    ("quarters", 25),
    ("dimes", 10),
    ("nickels", 5),
    ("pennies", 1),
]


def calculate_change(cost: float, payment: float) -> dict:
    if payment < cost:
        raise ValueError("Payment is less than cost.")

    # Convert to cents to avoid floating-point issues
    change_cents = int(round((payment - cost) * 100))

    result = {}
    remaining = change_cents

    for name, value in DENOMINATIONS:
        count = remaining // value
        remaining = remaining % value
        if count > 0:
            result[name] = count

    return result, change_cents


def main():
    if len(sys.argv) >= 3:
        try:
            cost = float(sys.argv[1])
            payment = float(sys.argv[2])
        except ValueError:
            print("Usage: python change_calculator.py <cost> <payment>")
            print("   or: python change_calculator.py")
            sys.exit(1)
    else:
        try:
            cost = float(input("Enter cost (C): $"))
            payment = float(input("Enter payment (P): $"))
        except ValueError:
            print("Invalid input. Please enter numeric values.")
            sys.exit(1)

    try:
        breakdown, total_cents = calculate_change(cost, payment)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    total_dollars = total_cents / 100
    print(f"\nChange owed: ${total_dollars:.2f}")

    if not breakdown:
        print("Exact change — no bills or coins needed.")
        return

    print("Breakdown:")
    for name, count in breakdown.items():
        label = name if count != 1 else name.rstrip("s")
        print(f"  {count} {label}")


if __name__ == "__main__":
    main()
