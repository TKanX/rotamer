#[path = "build/codegen.rs"]
mod codegen;
#[path = "build/residues/mod.rs"]
mod residues;

fn main() {
    println!("cargo:rerun-if-changed=build/");
    codegen::generate();
}
