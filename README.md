# SA-attacher

Attach SA-tag in bam file from the output of LAST-split.

(Under development; no guarantee to work)

## USAGE

```bash
cargo build --release
cargo run --release -- /* input.bam */ > /* output.bam */
```

Tips on brew environment

```bash
export PATH=$(brew --prefix)/opt/python3/libexec/bin:$PATH 
```
