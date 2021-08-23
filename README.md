# Bamhandler

* Calculate "Primary alignment rate" and other BAM file stats.
* Attach SA-tag in bam file from the output of LAST-split.

(Under development; no guarantee to work)

## USAGE

```bash
cargo build --release
cargo run --release -- help
cargo run --release -- stats /* input.bam */ > /* output.stats */
cargo run --release -- attachsa /* input.bam */ > /* output.bam */
```

Tips on brew:

```bash
export PATH=$(brew --prefix)/opt/python3/libexec/bin:$PATH 
```
