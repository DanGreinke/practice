#!/usr/bin/env python3
"""
Password Generator

Generates cryptographically secure passwords with configurable
character sets and length.
"""

import argparse
import math
import secrets
import string
import sys


AMBIGUOUS_CHARS = "0O1lI"


def build_pool(
    use_lower: bool = True,
    use_upper: bool = True,
    use_digits: bool = True,
    use_special: bool = True,
    extra_special: str = "",
    exclude_ambiguous: bool = False,
) -> str:
    """Build the character pool based on user options."""
    pool = ""
    if use_lower:
        pool += string.ascii_lowercase
    if use_upper:
        pool += string.ascii_uppercase
    if use_digits:
        pool += string.digits
    if use_special:
        pool += string.punctuation
    if extra_special:
        pool += extra_special

    if not pool:
        raise ValueError("At least one character set must be enabled.")

    if exclude_ambiguous:
        pool = "".join(ch for ch in pool if ch not in AMBIGUOUS_CHARS)

    if not pool:
        raise ValueError(
            "Character pool is empty after excluding ambiguous characters."
        )

    return pool


def generate_password(pool: str, length: int) -> str:
    """Generate a single password using secrets for cryptographic security."""
    return "".join(secrets.choice(pool) for _ in range(length))


def calculate_entropy(pool_size: int, length: int) -> float:
    """Estimate password entropy in bits."""
    return length * math.log2(pool_size)


def strength_label(entropy: float) -> str:
    """Map entropy to a human-readable strength label."""
    if entropy < 28:
        return "Very Weak"
    if entropy < 36:
        return "Weak"
    if entropy < 60:
        return "Reasonable"
    if entropy < 120:
        return "Strong"
    return "Very Strong"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate secure random passwords.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-l", "--length", type=int, default=16, help="Password length"
    )
    parser.add_argument(
        "-c", "--count", type=int, default=1, help="Number of passwords to generate"
    )
    parser.add_argument(
        "--no-lower", action="store_true", help="Exclude lowercase letters"
    )
    parser.add_argument(
        "--no-upper", action="store_true", help="Exclude uppercase letters"
    )
    parser.add_argument(
        "--no-digits", action="store_true", help="Exclude digits"
    )
    parser.add_argument(
        "--no-special", action="store_true", help="Exclude special characters"
    )
    parser.add_argument(
        "--extra-special",
        default="",
        help="Additional special characters to include",
    )
    parser.add_argument(
        "--exclude-ambiguous",
        action="store_true",
        help="Exclude visually ambiguous characters (0, O, 1, l, I)",
    )
    parser.add_argument(
        "--show-entropy",
        action="store_true",
        help="Display estimated entropy and strength for each password",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.length < 1:
        print("Error: Length must be at least 1.", file=sys.stderr)
        return 1
    if args.count < 1:
        print("Error: Count must be at least 1.", file=sys.stderr)
        return 1

    try:
        pool = build_pool(
            use_lower=not args.no_lower,
            use_upper=not args.no_upper,
            use_digits=not args.no_digits,
            use_special=not args.no_special,
            extra_special=args.extra_special,
            exclude_ambiguous=args.exclude_ambiguous,
        )
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    pool_size = len(pool)
    entropy = calculate_entropy(pool_size, args.length)
    label = strength_label(entropy)

    for i in range(args.count):
        password = generate_password(pool, args.length)
        output = password
        if args.show_entropy:
            output += f"  |  {entropy:.1f} bits  |  {label}"
        print(output)

    return 0


if __name__ == "__main__":
    sys.exit(main())
