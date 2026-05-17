# Change Calculator

A Python script that computes the optimal number of bills and coins to return as change when a customer pays amount **P** for a cost **C**.

## Supported Denominations

| Denomination | Value |
|---|---|
| $20 bills | $20.00 |
| $10 bills | $10.00 |
| $5 bills | $5.00 |
| $1 bills | $1.00 |
| quarters | $0.25 |
| dimes | $0.10 |
| nickels | $0.05 |
| pennies | $0.01 |

## Usage

### Interactive mode
```bash
python3 change_calculator.py
```

### Command-line mode
```bash
python3 change_calculator.py <cost> <payment>
```

### Example
```bash
$ python3 change_calculator.py 17.83 50
Change owed: $32.17
Breakdown:
  1 $20 bill
  1 $10 bill
  2 $1 bills
  1 dime
  1 nickel
  2 pennies
```

## Testing

```bash
python3 -m unittest test_change_calculator.py
```

## Features

- Converts to cents internally to avoid floating-point rounding errors
- Handles singular/plural labels automatically
- Validates inputs (no negative amounts, payment must cover cost)
