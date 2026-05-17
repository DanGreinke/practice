# Password Generator

A small, secure command-line password generator written in Python. It uses Python's `secrets` module for cryptographically strong randomness.

## Usage

```bash
python3 password_generator.py [options]
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-l`, `--length` | Password length | `16` |
| `-c`, `--count` | Number of passwords to generate | `1` |
| `--no-lower` | Exclude lowercase letters | disabled |
| `--no-upper` | Exclude uppercase letters | disabled |
| `--no-digits` | Exclude digits | disabled |
| `--no-special` | Exclude special characters | disabled |
| `--extra-special` | Add extra special characters | none |
| `--exclude-ambiguous` | Exclude `0`, `O`, `1`, `l`, `I` | disabled |
| `--show-entropy` | Show estimated entropy & strength | disabled |

### Examples

Generate one 16-character password:
```bash
python3 password_generator.py
```

Generate three 20-character passwords:
```bash
python3 password_generator.py -l 20 -c 3
```

Generate a password without special characters and show entropy:
```bash
python3 password_generator.py --no-special --show-entropy
```

Generate a password with only letters and digits, excluding ambiguous characters:
```bash
python3 password_generator.py --exclude-ambiguous --no-special
```

## Requirements

- Python 3.6+
