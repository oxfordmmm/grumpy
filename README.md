[![Tests](https://github.com/oxfordmmm/grumpy/actions/workflows/test.yaml/badge.svg)](https://github.com/oxfordmmm/grumpy/actions/workflows/test.yaml)

# grumpy
Re-implementation of [gumpy](https://github.com/oxfordmmm/gumpy) in Rust for speed

## Installation

### CLI
```
cargo install grumpy
```

### As a dependency for development
#### Rust crate
```
cargo add grumpy
```

#### Python package
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

## Citation
If you use `grumpy` in your work, please cite:
```
Westhead J, Baker CS, Brouard M, Colpus M, Constantinides B, Hall A, Knaggs J, Lopes Alves M, Spies R, Thai H, Surrall S, Govender K, Peto TEA, Crook DW, Omar SV, Turner R, Fowler PW
Characterising the performance of an antibiotic resistance prediction tool, gnomonicus, using a diverse testset of 2,663 Mycobacterium tuberculosis samples
Microbial Genomics 11:001592 doi:10.1099/mgen.0.001592
```

A BibTeX citation is included [here](CITATION.bib)