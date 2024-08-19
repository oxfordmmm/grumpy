[![Tests](https://github.com/oxfordmmm/grumpy/actions/workflows/test.yaml/badge.svg)](https://github.com/oxfordmmm/grumpy/actions/workflows/test.yaml)

# grumpy
Re-implementation of [gumpy](https://github.com/oxfordmmm/gumpy) in Rust for speed

## Installation

### Rust crate
```
cargo add grumpy
```

### Python package
```
pip install bio-grumpy
```

## Tests
Running Rust unit tests
```
cargo test
```

### Coverage
Test coverage can be found with the use of `tarpaulin`.
```
# Install tarpaulin on first run
cargo install cargo-tarpaulin

# Generate an HTML coverage report
cargo tarpaulin --no-dead-code --engine llvm --out html
```