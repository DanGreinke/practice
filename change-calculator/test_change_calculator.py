#!/usr/bin/env python3
"""Unit tests for change_calculator.py"""

import unittest
from change_calculator import calculate_change


class TestChangeCalculator(unittest.TestCase):

    def test_exact_change(self):
        result, cents = calculate_change(10.00, 10.00)
        self.assertEqual(cents, 0)
        self.assertEqual(result, {})

    def test_simple_bills(self):
        result, cents = calculate_change(0, 32.00)
        self.assertEqual(cents, 3200)
        self.assertEqual(result["$20 bills"], 1)
        self.assertEqual(result["$10 bills"], 1)
        self.assertEqual(result["$1 bills"], 2)

    def test_all_denominations(self):
        # $41.41 — uses every denomination at least once
        result, cents = calculate_change(0, 41.41)
        self.assertEqual(cents, 4141)
        self.assertEqual(result["$20 bills"], 2)
        self.assertNotIn("$10 bills", result)
        self.assertNotIn("$5 bills", result)
        self.assertEqual(result["$1 bills"], 1)
        self.assertEqual(result["quarters"], 1)
        self.assertEqual(result["dimes"], 1)
        self.assertEqual(result["nickels"], 1)
        self.assertEqual(result["pennies"], 1)

    def test_payment_less_than_cost(self):
        with self.assertRaises(ValueError):
            calculate_change(10.00, 5.00)

    def test_floating_point_safety(self):
        # Classic floating-point trap: 2.20 - 1.10 != 1.10
        result, cents = calculate_change(1.10, 2.20)
        self.assertEqual(cents, 110)
        self.assertEqual(result["$1 bills"], 1)
        self.assertEqual(result["dimes"], 1)

    def test_large_amount(self):
        result, cents = calculate_change(0, 1000.00)
        self.assertEqual(cents, 100000)
        self.assertEqual(result["$20 bills"], 50)


if __name__ == "__main__":
    unittest.main()
