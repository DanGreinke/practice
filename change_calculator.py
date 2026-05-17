#!/usr/bin/env python3
"""
Change Calculator

Given a cash payment P and a register cost C, computes the optimal change
using standard US denominations:
$20, $10, $5, $1, quarters, dimes, nickels, pennies.
"""

import sys

DENOMINATIONS = [
    (2000, "$20 bill"),
    (1000, "$10 bill"),
    (500, "$5 bill"),
    (100, "$1 bill"),
    (25, "quarter"),
    (10, "dime"),
    (5, "nickel"),
    (1, "penny"),
]


def to_cents(amount_str: str) -> int | None:
    """Parse a money string into an integer number of cents."""
    s = amount_str.strip().lstrip("$").replace(",", "")
    if not s:
        return None

    if "." in s:
        dollars, cents = s.split(".", 1)
        cents = (cents + "00")[:2]  # pad or truncate to 2 digits
    else:
        dollars = s
        cents = "00"

    try:
        return int(dollars) * 100 + int(cents)
    except ValueError:
        return None


def calculate_change(payment_cents: int, cost_cents: int):
    """Return a dict of denomination counts and the total change in cents."""
    if payment_cents < cost_cents:
        shortfall = cost_cents - payment_cents
        return None, shortfall

    change = payment_cents - cost_cents
    result = {}
    remaining = change

    for value, name in DENOMINATIONS:
        count = remaining // value
        remaining = remaining % value
        result[name] = count

    return result, change


def _pluralize(word: str, count: int) -> str:
    if count == 1:
        return word
    if word == "penny":
        return "pennies"
    return word + "s"


def display(result, total_change_cents: int):
    if result is None:
        print(
            f"Insufficient payment. Need an additional ${total_change_cents / 100:.2f}."
        )
        return

    print(f"\nTotal change: ${total_change_cents / 100:.2f}")
    print("-" * 30)
    for _, name in DENOMINATIONS:
        count = result[name]
        if count > 0:
            label = _pluralize(name, count)
            print(f"{count:>3} {label}")


def main():
    if len(sys.argv) >= 3:
        p_str, c_str = sys.argv[1], sys.argv[2]
    else:
        p_str = input("Enter payment amount (P): ")
        c_str = input("Enter cost amount   (C): ")

    p_cents = to_cents(p_str)
    c_cents = to_cents(c_str)

    if p_cents is None or c_cents is None:
        print("Error: Please enter valid numeric amounts.")
        sys.exit(1)

    result, change = calculate_change(p_cents, c_cents)
    display(result, change)


if __name__ == "__main__":
    main()
