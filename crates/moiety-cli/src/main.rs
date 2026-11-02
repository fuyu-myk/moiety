//! The `moiety` command: read condensed formulas, write them back, draw them, and
//! run the conformance suite.

mod args;
mod commands;

use std::process::ExitCode;

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    match args::parse(&argv).and_then(commands::run) {
        Ok(code) => code,
        Err(error) => {
            eprintln!("error: {error:#}");
            ExitCode::from(2)
        }
    }
}
